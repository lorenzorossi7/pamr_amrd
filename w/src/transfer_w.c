//=============================================================================
// 
//  transfer.c --- functions to communicating pieces of grids between nodes
// 
//=============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "transfer_w.h"
#include "misc_w.h"
#include "pamr_w.h"
#include <mpi.h>
#include "m_util_r8_w.h"
#include <time.h>
#include <math.h>

//=============================================================================
// a few utility functions
//=============================================================================

//=============================================================================
// for debugging ...
//=============================================================================
void check_integrity(gsl *p, char *place, int first_only)
{
   int i;

   while(p)
   {
      for (i=0; i<p->g->dim; i++)
      {
         if (p->g->shape[i]<=0)
         {
            printf("shape check_integrity failed in %s: shape %i = %i\n",place,i,p->g->shape[i]);
         }
         if (p->bbox[2*i]>p->bbox[2*i+1])
         {
            printf("bbox check_integrity failed in %s: %i = [%lf,%lf]\n",place,i,p->bbox[2*i],p->bbox[2*i+1]);
         }
      }
      if (first_only) p=0; else p=p->next;
   }
}

void free_gsl(gsl *p)
{
   gsl *q;

   while(p)
   {
      q=p;
      p=p->next;
      free(q);
   }
}

int size_transfer_block(int size,int size_c)
{
   int n,i;

   for (i=0,n=0; i<curr_context->gfns; i++) 
      if (curr_context->tf_bits[i]) 
         if (curr_context->gfn_var_type[i]==PAMR_VERTEX_CENTERED) n+=size;
         else n+=size_c;

   return n;
}

gsl *clone_gsl(gsl *p, int first_only)
{
   gsl *np=0,*q=0,*pq=0;
   int i;

   while(p)
   {
      if (!(q=(gsl *)malloc(sizeof(struct gsl))))
      {
         printf("clone_gsl: error ... out of memory\n");
         free_gsl(np);
         return 0;
      }
      if (pq) pq->next=q; else np=q;
      np->last=q;
      q->next=0;
      q->g=p->g;
      for (i=0; i<p->g->dim; i++) 
      {
         q->bbox[2*i]=p->bbox[2*i];
         q->bbox[2*i+1]=p->bbox[2*i+1];
      }

      if (first_only) return np;
      pq=q;
      p=p->next;
   }
   return np;
}

//=============================================================================
// calculates the index bounding box for the current grid segment within
// the given grid ... non-aligned boundaries are rounded "up" to the
// nearest grid line
//
// NOTE: in general ibbox_c[2*i]==ibbox[2*i] (left edges line up)
//       but ibbox_c[2*i+1]=ibbox[2*i] (one EXTRA cell)
//       unless the segment spans the entire dimension, 
//       whence ibbox_c[2*i+1]=ibbox[2*i]-1
// 
//       some routines use calc_gs_ibbox and allocate a grid-structure
//       from this particular segment; the bbox there will then only
//       in general correctly refer the the vertex centered portion!
//=============================================================================
void calc_gs_ibbox(gsl *gs, int *ibbox, int *ibbox_c, real *dx)
{
   int i,d=gs->g->dim;
   grid *g=gs->g;
   int ltrace=0;
   int uflush;

   for (i=0; i<2*PAMR_MAX_DIM; i++) {
     ibbox_c[i] = 0;
     ibbox[i] = 0;
   }
   
   for (i=0; i<d; i++)
   {
      if (g->shape[i]>1)
      {
         dx[i]=(g->bbox[2*i+1]-g->bbox[2*i])/(g->shape[i]-1);
         // * plus 1 since the index is 1-based
         // * plus 0.51 for the min. boundary since in 
         //   extend_gsl (as of V2) we're always subtracting 0.5,
         //   and we want to include that point (and C truncates)
         // * plus 0.01 for the max. boundary to make sure
         //   the truncation includes that point (don't need
         //   to worry about the +0.5 in gsl_extend, as the truncation
         //   gets rid of it.
         ibbox[2*i]=(gs->bbox[2*i]-g->bbox[2*i])/dx[i]+1.51;
         ibbox[2*i+1]=(gs->bbox[2*i+1]-g->bbox[2*i])/dx[i]+1.01;

         // FP: replaced the commented-out stuff below with the
         // following identical logical, expect the test for
         // uflush now doesn't assume the coordinate "0" lines up with
         // a vertex point on all levels
         ibbox_c[2*i]=ibbox[2*i];
         // extend cell domain by 1 in the right (+ve) direction (above), *except*
         // if the right boundary segment is already flush to a vertex
         // (which at present will only happen if the right boundary
         //  is coincident with the parent grid boundary).
         uflush = 0; 
         ibbox_c[2*i+1]=ibbox[2*i+1];
         if ( ( (int) ((gs->bbox[2*i+1]-g->bbox[2*i])*2/dx[i] + 0.01) ) % 2 == 0) uflush = 1;
         if (uflush) ibbox_c[2*i+1]--;

      } else {
        ibbox[2*i]=1;
        ibbox[2*i+1]=1;
        ibbox_c[2*i]=1;
        ibbox_c[2*i+1]=1;
        dx[i]=PAMR_SMALL_DX;
      }

      if (g->shape_c[i] < 1) {
        ibbox_c[2*i]=1;
        ibbox_c[2*i+1]=0;
      }

   }

   if (ltrace)
   {
      printf("calc_gs_ibbox: g->shape   =%i,%i,%i\n",g->shape[0],g->shape[1],g->shape[2]);
      printf("               g->shape_c =%i,%i,%i\n",g->shape_c[0],g->shape_c[1],g->shape_c[2]);
      printf("               g->bbox  =%lf,%lf,%lf,%lf,%lf,%lf\n",g->bbox[0],g->bbox[1],g->bbox[2],g->bbox[3],g->bbox[4],g->bbox[5]);
      printf("               gs->bbox =%lf,%lf,%lf,%lf,%lf,%lf\n",gs->bbox[0],gs->bbox[1],gs->bbox[2],gs->bbox[3],gs->bbox[4],gs->bbox[5]);
      printf("               ibbox   =%i,%i,%i,%i,%i,%i\n",ibbox[0],ibbox[1],ibbox[2],ibbox[3],ibbox[4],ibbox[5]);
      printf("               ibbox_c =%i,%i,%i,%i,%i,%i\n",ibbox_c[0],ibbox_c[1],ibbox_c[2],ibbox_c[3],ibbox_c[4],ibbox_c[5]);
   }
}

//=============================================================================
// aligns (by first trying to expand, then trimming if not possible) 
// src bbox with dst grid lines (for coarsening)
//=============================================================================
void align_src_bbox(gsl *src, gsl *dst)
{
   int i,d=src->g->dim;
   grid *gs=src->g,*gd=dst->g;
   real dxs,dxd;

   IFG(3) printf("   << align_src_bbox >>\n");
   IFG(4) printf("before alignment, src->bbox[0,1] = [ %f, %f] \n", src->bbox[0], src->bbox[1]);

   for (i=0; i<d; i++)
   {
      if (gs->shape[i]>1 && gd->shape[i]>1)
      {
         dxs=(gs->bbox[2*i+1]-gs->bbox[2*i])/(gs->shape[i]-1);
         dxd=(gd->bbox[2*i+1]-gd->bbox[2*i])/(gd->shape[i]-1);
         src->bbox[2*i]=((int)((src->bbox[2*i]-gd->bbox[2*i]+dxs/2)/dxd))*dxd+gd->bbox[2*i];
         if (fuzz_lt(src->bbox[2*i],gs->bbox[2*i],dxs/4)) src->bbox[2*i]+=dxd;
         dst->bbox[2*i]=src->bbox[2*i];
         src->bbox[2*i+1]=((int)((src->bbox[2*i+1]-gd->bbox[2*i]+dxd-dxs/2)/dxd))*dxd+gd->bbox[2*i];
         if (fuzz_gt(src->bbox[2*i+1],gs->bbox[2*i+1],dxs/4)) src->bbox[2*i+1]-=dxd;
         dst->bbox[2*i+1]=src->bbox[2*i+1];
      }
   }

   IFG(4) printf("after alignment, src->bbox[0,1] = [ %f, %f] \n", src->bbox[0], src->bbox[1]);
}

//=============================================================================
// aligns (by first trying to expand, then trimming if not possible) 
// src bbox with src grid lines (for interpolation)
//=============================================================================
void align_src_bbox_interp(gsl *src)
{
   int i,d=src->g->dim;
   grid *gs=src->g;
   real dx, temp1, temp2;
   int i1, i2;

   IFG(3) printf("   << align_src_bbox_interp >>\n");
   IFG(4) printf("before alignment, src->bbox[0,1] = [ %f, %f] \n", src->bbox[0], src->bbox[1]);

   for (i=0; i<d; i++)
   {
      if (gs->shape[i]>1)
      {
         dx=(gs->bbox[2*i+1]-gs->bbox[2*i])/(gs->shape[i]-1);
         temp1 = (src->bbox[2*i]-gs->bbox[2*i])/dx;
         i1 = (int) temp1;
         src->bbox[2*i] = gs->bbox[2*i] + i1*dx;
         if (fuzz_lt(src->bbox[2*i],gs->bbox[2*i],dx/4)) src->bbox[2*i]+=dx;
         temp2 = (src->bbox[2*i+1]-gs->bbox[2*i])/dx;
         i2 = (int) temp2;
         if (temp2 > ((real) i2) ) src->bbox[2*i+1] = gs->bbox[2*i] + (i2+1)*dx;
         if (fuzz_gt(src->bbox[2*i+1],gs->bbox[2*i+1],dx/4)) src->bbox[2*i+1]-=dx;
      }
   }

   IFG(4) printf("after alignment, src->bbox[0,1] = [ %f, %f] \n", src->bbox[0], src->bbox[1]);
}


//=============================================================================
// alloc_gsl(g) creates a grid-segment copy of g. If interior>=1, then
// the bbox is set to the interior of g (i.e. excluding g->ghost_width);
// if AMR_reduce>0, then the grid segment is reduced by AMR_reduce points
// in size along any AMR boundary (*except* along ghost_width boundaries if
// interior is 1), up to a maximum amount of ghost_width if ghost_width is >0
//
// The next field of the returned gsl is set to zero.
//=============================================================================
gsl *alloc_gsl(grid *g, int lev, int interior, int AMR_reduce, int AMR_extra_half)
{
   gsl *p;
   int d=g->dim,i,ibbox[2*PAMR_MAX_DIM],ibbox_c[2*PAMR_MAX_DIM];
   real dx,dxt[2*PAMR_MAX_DIM];
   context *c=curr_context;

   if (!(p=(gsl *)malloc(sizeof(struct gsl))))
   {
      printf("alloc_gs: error ... out of memory\n");
      return 0;
   }

   p->g=g;
   p->next=0;
   p->last=p;

   // Initialize to make valgrind happy when it calls calc_gs_ibbox
   for (i=0; i<d; i++) {
     p->bbox[2*i]=g->bbox[2*i];
     p->bbox[2*i+1]=g->bbox[2*i+1];
   }

   for (i=0; i<d; i++) 
   {
      if (g->shape[i]>1) dx=(g->bbox[2*i+1]-g->bbox[2*i])/(g->shape[i]-1);
      else if (g->shape[i]==1) 
      { 
         dx=PAMR_SMALL_DX; printf("alloc_gsl: WARNING ... g->shape[i]=1\n"); 
      }
      else { printf("alloc_gsl: ERROR ... g->shape[0]<1\n"); return 0; }

      if (interior && g->ghost_width[2*i])
      {
         p->bbox[2*i]+=g->ghost_width[2*i]*dx;
      }
      else if (AMR_reduce>0 && !g->ghost_width[2*i])  
      {
        if (AMR_extra_half) {
          if (fuzz_gt(p->bbox[2*i],c->bbox[2*i],dx/2) || c->periodic[i]) p->bbox[2*i]+=(AMR_reduce+0.5)*dx;
        } else {
          if (fuzz_gt(p->bbox[2*i],c->bbox[2*i],dx/2) || c->periodic[i]) p->bbox[2*i]+=(AMR_reduce)*dx;
        }
      }
      if (interior && g->ghost_width[2*i+1]) 
      {
         p->bbox[2*i+1]+=-g->ghost_width[2*i+1]*dx; 
      }
      else if (AMR_reduce>0 && !g->ghost_width[2*i+1]) 
      {
        if (AMR_extra_half) {
          if (fuzz_lt(p->bbox[2*i+1],c->bbox[2*i+1],dx/2) || c->periodic[i]) p->bbox[2*i+1]+=-(AMR_reduce+0.5)*dx; 
        } else {
          if (fuzz_lt(p->bbox[2*i+1],c->bbox[2*i+1],dx/2) || c->periodic[i]) p->bbox[2*i+1]+=-(AMR_reduce)*dx; 
        }

      }

      // the following just trims virtual grids for efficiency

      if (c->periodic[i]) p->bbox[2*i]=max(p->bbox[2*i],c->bbox[2*i]); 
      if (c->periodic[i]) p->bbox[2*i+1]=min(p->bbox[2*i+1],c->bbox[2*i+1]+(2*c->ghost_width[i])*dx); 

      if (c->periodic[i] || interior || AMR_reduce) 
      {

         calc_gs_ibbox(p,ibbox,ibbox_c,dxt);
         if ((ibbox[2*i+1]-ibbox[2*i]+1)<1) 
         {
            printf("alloc_gsl ... freeing node. ghost_width=[%i,%i], shape=%i\n",g->ghost_width[2*i],g->ghost_width[2*i+1],g->shape[i]);
            free(p); 
            return 0;
         }
      }
   }

   return p;
}

//=============================================================================
// cat_gsls() concatenates two gsl's
//=============================================================================
inline gsl *cat_gsls(gsl *A, gsl *B)
{
   if (!A) return B;
   else if (!B) return A;

   A->last->next=B;
   A->last=B->last;

   return A;
}

//=============================================================================
// gsl_subtract()/gs_subtract returns a list of grid segments denoting the 
// CSG operation C=<A-B> (i.e. A and -union(B)). The splicing operation
// is implemented dimension by dimension, from 1 to d.
// 
// NOTE: grids of C inherit the 'data' structures of corresponding grids of A;
// only the bbox's of B are used here.
//
// gs_subtract(A,B) assumes A and B to be single grid-segments (i.e. it ignores
// the ->next field), while in gsl_subtract(A,B) they can be lists of 
// grid-segments.
//
// c_bdy_flag = 1 when using this routine to determine the gsl for 
// doing interpolation on an AMR boundary for a cell centered variable.
//=============================================================================
gsl *gs_subtract(gsl *A, gsl *B, int c_bdy_flag)
{
   real cut_plane[2*PAMR_MAX_DIM],dx[PAMR_MAX_DIM];
   int d=A->g->dim,i,j;
   gsl *C=0,*q;

   IFG(3) printf("   << gs_subtract >>\n");

   if (!A) return 0;
   if (!B) return clone_gsl(A,1);

   for(i=0; i<d; i++)
   {
      cut_plane[2*i]=A->bbox[2*i];
      cut_plane[2*i+1]=A->bbox[2*i+1];
      if (A->g->shape[i]>1) dx[i]=(A->g->bbox[2*i+1]-A->g->bbox[2*i])/(A->g->shape[i]-1);
      else dx[i]=PAMR_SMALL_DX;
      if (fuzz_gt(B->bbox[2*i],A->bbox[2*i+1],dx[i]/2) ||
          fuzz_gt(A->bbox[2*i],B->bbox[2*i+1],dx[i]/2)) return clone_gsl(A,1);
   }

   for(i=0; i<d; i++)
   {
      cut_plane[2*i]=max(A->bbox[2*i],B->bbox[2*i]);
      if (fuzz_gt(cut_plane[2*i],A->bbox[2*i],dx[i]/2))
      {
         q=clone_gsl(A,1); 
         if (C) { q->next=C; q->last=C->last; }
         C=q;
         for (j=0; j<d; j++)
         {
            if (i==j) 
            {
               C->bbox[2*i]=A->bbox[2*i];
               if (c_bdy_flag) {
                 C->bbox[2*i+1]=max(C->bbox[2*i],cut_plane[2*i]);
               } else {
                 C->bbox[2*i+1]=max(C->bbox[2*i],cut_plane[2*i]-dx[i]);
               }
            }
            else
            {
               C->bbox[2*j]=cut_plane[2*j];
               C->bbox[2*j+1]=cut_plane[2*j+1];
            }
         }
      }

      cut_plane[2*i+1]=min(A->bbox[2*i+1],B->bbox[2*i+1]);
      if (fuzz_lt(cut_plane[2*i+1],A->bbox[2*i+1],dx[i]/2))
      {
         q=clone_gsl(A,1); 
         if (C) { q->next=C; q->last=C->last; }
         C=q;
         for (j=0; j<d; j++)
         {
            if (i==j) 
            {
               C->bbox[2*i+1]=A->bbox[2*i+1];
               if (c_bdy_flag) {
                 C->bbox[2*i]=min(C->bbox[2*i+1],cut_plane[2*i+1]);
               } else {
                 C->bbox[2*i]=min(C->bbox[2*i+1],cut_plane[2*i+1]+dx[i]);
               }
            }
            else
            {
               C->bbox[2*j]=cut_plane[2*j];
               C->bbox[2*j+1]=cut_plane[2*j+1];
            }
         }
      }

   }
   
   return C;
}

gsl *gsl_subtract(gsl *A, gsl *B, int c_bdy_flag)
{
   gsl *C=0,*C0,*C1;

   C=clone_gsl(A,0);

   while(B)
   {
      C0=0;
      C1=C;
      while(C1)
      {
         C0=cat_gsls(C0,gs_subtract(C1,B,c_bdy_flag));
         C1=C1->next;
      }
      free_gsl(C);
      C=C0;
      B=B->next;
   }
   
   //--------------------------------------------------------------------------
   // The above operation could create a myriad of tiny grids, some of
   // which could be merged together. If we notice efficiency problems
   // due to too many grids, we may want to implement a merge operation.
   //--------------------------------------------------------------------------
   // merge_gsl(C); 
   return(C);
}

//=============================================================================
// gs_and(A,B,A0,B0) computes <A and B>, and the result is returned in 
// A0 and B0. A0 and B0 will share the same bbox (or be zero if the grids
// do not overlap), but will inherit the geometry from A and B respectively. 
//
// NOTE: no guarantees that the new shape array will make sense,
//       in particular if the overlap does not occur on cell boundaries
//       of both grids. Use the bbox for the correct geometry.
//=============================================================================
int gs_and(gsl *A, gsl *B, gsl **A0, gsl **B0)
{
   real new_bbox[2*PAMR_MAX_DIM],dxA[PAMR_MAX_DIM],dxB[PAMR_MAX_DIM],dxm[PAMR_MAX_DIM];
   int d=A->g->dim,i;

   *A0=*B0=0;

   for(i=0; i<d; i++)
   {
      if (A->g->shape[i]>1) dxA[i]=(A->g->bbox[2*i+1]-A->g->bbox[2*i])/(A->g->shape[i]-1);
      else dxA[i]=PAMR_SMALL_DX;
      if (B->g->shape[i]>1) dxB[i]=(B->g->bbox[2*i+1]-B->g->bbox[2*i])/(B->g->shape[i]-1);
      else dxB[i]=PAMR_SMALL_DX;
      dxm[i]=max(dxA[i],dxB[i]);
      new_bbox[2*i]=max(A->bbox[2*i],B->bbox[2*i]);
      new_bbox[2*i+1]=min(A->bbox[2*i+1],B->bbox[2*i+1]);
      if (fuzz_lt(new_bbox[2*i+1],new_bbox[2*i],dxm[i]/4)) return 0;
   }

   if (!(*A0=clone_gsl(A,1))) return 0;
   if (!(*B0=clone_gsl(B,1))) {free_gsl(*A0); return 0;}

   for(i=0; i<d; i++)
   {
      (*A0)->bbox[2*i]=new_bbox[2*i];
      (*A0)->bbox[2*i+1]=new_bbox[2*i+1];
      (*B0)->bbox[2*i]=new_bbox[2*i];
      (*B0)->bbox[2*i+1]=new_bbox[2*i+1];
   }

   return 1;
}

//=============================================================================
// utility routine used by build_owned_gsl() below, for interpolation
//
// extend in both +ve and -ve directions, however, in most cases 
// it will only be possible to do this in 1 direction
//
// PAMR V2 UPDATE: A similar problem arose for CC operations that extend_gsl()
// was designed to solve for VC interpolation. Namely, the CSG operations
// are based on VC bounding boxes, and so 'gaps' in ownership can 
// arise from the CC grids perspective. e.g. 
//
// * --> VC
// . --> CC
//
// ghostwidth=1
//
// grid A, ----.----*----.----*----.----*----.----*
// grid B,                              *----.----*----.----*----.----*----
//
// grid A, owned: --*----.----*----.----*
// grid B, owned:                                 *----.----*----.----*----
//
// ghostwidth=2
//
// grid A, ----.----*----.----*----.----*----.----*----.----*
// grid B,                    *----.----*----.----*----.----*----.----*----
//
// grid A, owned: --*----.----*----.----*
// grid B, owned:                                 *----.----*----.----*----
//
// So now extending by 1/2 cell in both directions to include cells
// NOTE: number of "proper" ghost cells is ghostwidth-1, and possible
//       ambiguity in who owns intermediate point (single-valuedness should
//       still be OK though). Could imagine a situation where we allow
//       for 2 ghost cells, though then we'd have a similar problem
//       with VC structures.
//=============================================================================
void extend_gsl(gsl *gs, int extend)
{
   real dx,extend_frac_p,extend_frac_m;
   grid *g;
   int d,i;

   extend_frac_p=0.5e0; extend_frac_m=0.5e0; 
   if (extend) { extend_frac_p=extend_frac_p+1; extend_frac_m=extend_frac_m+1; }

   while(gs)
   {
      g=gs->g;
      d=gs->g->dim;

      for(i=0; i<d; i++)
      {
         if (g->shape[i]>1)
         {
            dx=(g->bbox[2*i+1]-g->bbox[2*i])/(g->shape[i]-1);
            gs->bbox[2*i]=max(gs->bbox[2*i]-dx*extend_frac_m,g->bbox[2*i]);
            gs->bbox[2*i+1]=min(gs->bbox[2*i+1]+dx*extend_frac_p,g->bbox[2*i+1]);
         }
      }

      gs=gs->next;
   }

}

//=============================================================================
//  build_owned_gsl() returns the list of grid-segments owned by node 'rank'.
//  The ownership of a piece of a grid by a node is defined as follows:
//  
//  On a given level, grids are numbered (arbitrary but *unique* order,
//  as in DV) from 1. given grid number g, located on node n, 
//  the piece(s) of g owned by n is defined as
//
//      owned_bbox(g) = interior_bbox(g) - owned_bbox(g-1) - ... - owned_bbox(1)
//                    = interior_bbox(g) - interior_bbox(g-1) - ... - interior_bbox(1)
//
//  If the flag interior>=1, then the interior of a grid is the piece of the grid 
//  excluding parallel ghost zones (i.e. the g->ghost_width flag), 
//  else the 'interior' is simply the entire grid.
//  NOTE: if interior>1, then alloc_gsl is called with interior=1 and 
//        AMR_reduce=interior-1
//
//  If extend=1, then after the owned list is generated, it is extended by
//  1 grid width (if possible) along the + directions. Thus the owned
//  list will not be unique there, but this provides a simple solution
//  to the problem of missing points on the fine level in this region
//  when interpolating.
//
//  NOTE: the '-' operation above is a set-subtraction operation (i.e. 
//  as in constructive solid geometry), and cannot be simplified via elementary
//  algebraic operations.
//=============================================================================
gsl *build_owned_gsl(int rank, cgh *gh, int l, int interior, int extend)
{
   gsl *owned=0,*int_gs,*new_gsl,*int_gsl=0,*gs;
   grid *g;
   level *lev;
   context *c=curr_context;
   int cinterior,cAMR_reduce;
   int AMR_extra_half = 1;

   IFG(2) printf("   >> build_owned_gsl(l=%i,interior=%i,extend=%i)\n",l,interior,extend);

   if (interior) cinterior=1; else cinterior=0;
   if (interior>1) cAMR_reduce=interior-1; else cAMR_reduce=0;

   if (!(lev=gh->levels[l-1]))
   {
      printf("build_owned_gsl: error ... level %i does not exist\n",l);
      return 0;
   }

   g=lev->grids;

   while(g)
   {
      if (g->comm)
      {
         int_gs=alloc_gsl(g,l,cinterior,cAMR_reduce,AMR_extra_half);
         if (g->rank==rank)
         {
           new_gsl=gsl_subtract(int_gs,int_gsl,0);
           owned=cat_gsls(owned,new_gsl);
         }
         int_gsl=cat_gsls(int_gsl,int_gs);
      }

      g=g->next;
   }

   if (owned) extend_gsl(owned,extend); 

   free_gsl(int_gsl);

   IFG(2) printf("      build_owned_gsl >> \n");
   return owned;
}

//=============================================================================
// the following builds a gsl of all grids (or all grid-AMR-boundaries of
// size AMR_bdy_only, of the *union* of grids, if AMR_bdy_only>1) in gh at 
// level l
//=============================================================================
gsl *build_complete_gsl(cgh *gh, int l, int AMR_bdy_only)
{
   gsl *cgs,*cgsl=0,*igs=0,*bgsl=0,*igsl=0;
   grid *g;
   level *lev;
   context *c=curr_context;
   int i, c_bdy_flag;
   int AMR_extra_half = 0;

   IFG(3) printf("   << build_complete_gsl on l=%d >>\n",l);

   if (!(lev=gh->levels[l-1]))
   {
      printf("build_complete_gsl: error ... level %i does not exist\n",l);
      return 0;
   }

   g=lev->grids;
   
   while(g)
   {
      if (g->comm)
      {
         cgs=alloc_gsl(g,l,0,0,AMR_extra_half);

         cgsl=cat_gsls(cgsl,cgs);
         if (AMR_bdy_only>0) 
         {
            igs=alloc_gsl(g,l,0,AMR_bdy_only,AMR_extra_half);
            igsl=cat_gsls(igsl,igs);
         }
      }
      g=g->next;
   }

   if (AMR_bdy_only>0) 
   {

     // Determine c_bdy_flag by looking at whose tf_bits are set.
     c_bdy_flag=0;
     for (i=0; i<c->gfns;i++) {
       if (c->tf_bits[i] && c->gfn_var_type[i]==PAMR_CELL_CENTERED) c_bdy_flag=1;
     }

      bgsl=gsl_subtract(cgsl,igsl,c_bdy_flag);
      free_gsl(igsl); 
      free_gsl(cgsl); 
      cgsl=bgsl;
      // extend_gsl(cgsl,0); // MARTIN
   }

   return cgsl;
}

//=============================================================================
// build_gstl() builds a grid transfer list, which is the *outer* product
// of in_src with in_dst, using the AND operation. The result
// is stored in the paired list [out_src,out_dst] ... i.e., each list has
// identical structure, and corresponding elements have the same bbox.
//
// transfer's between identical grids in the current node are discarded
//=============================================================================
int build_gstl(gsl *in_src,gsl *in_dst, gsl **out_src, gsl **out_dst)
{
   gsl *d,*p;
   gsl *s,*d2;
   int i, j;

   IFG(3) printf("   << build_gstl >>\n");

   *out_src=*out_dst=0;
   i = 0;

   while(in_src)
   {
      d=in_dst;
      j = 0;
      while(d)
      {
         if (!(in_src->g->gfs && (in_src->g==d->g)))
         {
            if (gs_and(in_src,d,&s,&d2))
            {
               *out_src=cat_gsls(*out_src,s);
               *out_dst=cat_gsls(*out_dst,d2);
            }
         }
         d=d->next;
         j++;
      }
      in_src=in_src->next;
      i++;
   }
   return 1;
}

void clear_ex(real *f, real *chr, real ex, int n)
{
   int i;
   
   for (i=0; i<n; i++) if (chr[i]==ex) f[i]=0;
}

//=============================================================================
// utility routine used by data_packer and local_transfer
//=============================================================================
int calc_delta(gsl *src, gsl *dst, int *irho)
{
   real dxsm=0,dxdm=0;
   int i,dim=src->g->dim;

   *irho=0;

   for (i=0; i<dim; i++) 
   {
      if (src->g->shape[i]>1 && dst->g->shape[i]>1) 
      {
         dxsm=max(dxsm,(src->g->bbox[2*i+1]-src->g->bbox[2*i])/(src->g->shape[i]-1));
         dxdm=max(dxdm,(dst->g->bbox[2*i+1]-dst->g->bbox[2*i])/(dst->g->shape[i]-1));
      }
   }
   if (dxsm==0 || dxdm==0) return 0;
   else if (fuzz_eq(dxsm,dxdm,dxsm/100))
   {
      return 0;
   }
   else if (fuzz_gt(dxsm,dxdm,dxsm/100))
   {
      *irho=dxsm/dxdm+0.5;
      return -1;
   }
   else
   {
      *irho=dxdm/dxsm+0.5;
      return 1;
   }
}

//=============================================================================
// The following routine packs(dir=PACK) or unpacks(dir=UNPACK) data.
//
// src -> dst is the usual transfer gsl pair
// going from my_rank() to rank (PACK), or vice-versa (UNPACK). 
//
// packing reads info from src, and fills (preallocated) data array.
//
// unpacking copies data from data array to src, allocating a *new* 
// local grid src->g (i.e. DO NOT CALL WITH rank=my_rank()!), such that
// a local_transfer() will correctly execute the src -> dst after
// communication.
// NOTE: the caller is responsible for freeing this new grid!
// (use free_src_data())
//
// if 'data' is zero, then a virtual pack/unpack is performed,
// returning the size of the expected data argument.
//=============================================================================
int data_packer(real *data, gsl *src, gsl *dst, int rank, int dir)
{
   int dsize=0,size,size_c,dim=curr_context->dim,i;
   context *c=curr_context;
   int *tbits=curr_context->tf_bits;
   grid *g;
   int nx1,ny1,nz1,nx2,ny2,nz2,i1,j1,k1,i2,j2,k2,ni,nj,nk;
   int nx1_c,ny1_c,nz1_c,nx2_c,ny2_c,nz2_c,i1_c,j1_c,k1_c,i2_c,j2_c,k2_c,ni_c,nj_c,nk_c;
   int irho,ibbox[2*PAMR_MAX_DIM],ibbox_d[2*PAMR_MAX_DIM],ibbox_c[2*PAMR_MAX_DIM];
   int ibbox_d_c[2*PAMR_MAX_DIM];
   real dx[PAMR_MAX_DIM],s_bbox[2*PAMR_MAX_DIM];
   int one=1,delta;
   real *chr,*chr_c,ex;
   int do_ex;

   IFG(3) printf("    >> data_packer(data=%p, dir=%i,rank=%i, my_rank=%i) \n",data,dir,rank,c->rank);

   do_ex=0;
   if (c->excision_on)
   {
      ex=c->excised;
      do_ex=1;
   }

   while(src && dst)
   {
      if ((dir==PACK && dst->g->rank==rank && src->g->rank == c->rank) ||
          (dir==UNPACK && src->g->rank==rank && dst->g->rank == c->rank))
      {
         //--------------------------------------------------------------------
         // 1. level to level transfer and 2. interpolation (at this stage
         //    interpolation is basically a copy, so can join functions here)
         //
         //    delta=0  ... copy
         //    delta=-1 ... interpolate
         //    delta=1  ... inject
         //
         // NOTE: index bbox is based on vertex centered ... manually adjusting
         // for cell centered
         //--------------------------------------------------------------------
         
         delta=calc_delta(src,dst,&irho);

         if ((delta==0)||(delta==-1))
         {
            //--------------------------------------------------------------------
            // For interpolation:
            // if we are far from source boundaries, then we can always use
            // interior interpolation operators, hence add such a buffer here, by
            // (temporarily) adjusting grid-segment parameters
            //
            // First save copy of bbox
            //--------------------------------------------------------------------
            for (i=0; i<dim; i++) 
            {
               s_bbox[2*i]=src->bbox[2*i];
               s_bbox[2*i+1]=src->bbox[2*i+1];
            }

            // Branson:  If we first align the source bbox
            // with it's own underlying grid lines then
            // interpolation should work even without if 
            // iterp_buffer = 0.
            if (delta==-1) align_src_bbox_interp(src); 

            for (i=0; i<dim; i++) 
            {
               if (src->g->shape[i]>1)
               {
                  dx[i]=(src->g->bbox[2*i+1]-src->g->bbox[2*i])/(src->g->shape[i]-1);
                  if ((delta==-1)&&(int)((src->bbox[2*i]-src->g->bbox[2*i])/dx[i]+0.5)>=c->interp_buffer) 
                  {
                    src->bbox[2*i]-=c->interp_buffer*dx[i];
                    src->bbox[2*i]=max(src->bbox[2*i],src->g->bbox[2*i]);
                  }  // else try a smaller interp buffer
                  else
                  {
                    if ((delta==-1)&&(int)((src->bbox[2*i]-src->g->bbox[2*i])/dx[i]+0.5)>=
                        (c->interp_buffer-1) && (c->interp_buffer-1)>0)
                      {
                        src->bbox[2*i]-=(c->interp_buffer-1)*dx[i];
                        src->bbox[2*i]=max(src->bbox[2*i],src->g->bbox[2*i]);
                      }
                  }

                  if ((delta==-1)&&(int)((src->g->bbox[2*i+1]-src->bbox[2*i+1])/dx[i]+0.5)>=c->interp_buffer) 
                  {
                    src->bbox[2*i+1]+=c->interp_buffer*dx[i];
                    src->bbox[2*i+1]=min(src->bbox[2*i+1],src->g->bbox[2*i+1]);
                  }  
                  else
                  {
                    if ((delta==-1)&&(int)((src->g->bbox[2*i+1]-src->bbox[2*i+1])/dx[i]+0.5)>=
                        (c->interp_buffer-1) && (c->interp_buffer-1)>0)
                      {
                        src->bbox[2*i+1]+=(c->interp_buffer-1)*dx[i];
                        src->bbox[2*i+1]=min(src->bbox[2*i+1],src->g->bbox[2*i+1]);
                      }
                  }

               } else dx[i]=PAMR_SMALL_DX;

            }

            //---------------------------------------------------------
            // NOTE: the following size calculation is unaware of
            // cell-centered, and so will be a little wasteful, however
            // I don't think this will have a significant impact
            // on communications
            // ALSO: ibbox_c will in general have 1 extra cell volume
            // in each direction compared to ibbox, so the subsequent
            // bbox in the allocated grid structure only correctly
            // refers to the vertex centered functions!
            //---------------------------------------------------------

            calc_gs_ibbox(src,ibbox,ibbox_c,dx);
            size=size_c=1;
            for (i=0; i<dim; i++) { size*=(ibbox[2*i+1]-ibbox[2*i]+1); size_c*=(ibbox_c[2*i+1]-ibbox_c[2*i]+1); }

            dsize+=size_transfer_block(size,size_c);

            if (data) 
            {
               ni=nj=nk=i1=i2=j1=j2=k1=k2=nx1=ny1=nz1=nx2=ny2=nz2=1;
               ni_c=nj_c=nk_c=i1_c=i2_c=j1_c=j2_c=k1_c=k2_c=nx1_c=ny1_c=nz1_c=nx2_c=ny2_c=nz2_c=1;
               g=src->g;

               switch(dim)
               {
                  case 3: 
                    k1=ibbox[4];
                    k1_c=ibbox_c[4];
                    nk=nz2=ibbox[5]-ibbox[4]+1; 
                    nk_c=nz2_c=ibbox_c[5]-ibbox_c[4]+1;
                    nz1=g->shape[2]; nz1_c=g->shape_c[2];
                  case 2: 
                    j1=ibbox[2];
                    j1_c=ibbox_c[2];
                    nj=ny2=ibbox[3]-ibbox[2]+1; 
                    nj_c=ny2_c=ibbox_c[3]-ibbox_c[2]+1;
                    ny1=g->shape[1]; ny1_c=g->shape_c[1];
                  case 1: 
                    i1=ibbox[0];
                    i1_c=ibbox_c[0];
                    ni=nx2=ibbox[1]-ibbox[0]+1; 
                    ni_c=nx2_c=ibbox_c[1]-ibbox_c[0]+1;
                    nx1=g->shape[0]; nx1_c=g->shape_c[0];
               }
               if (dir==UNPACK)
               {
                  if (g->gfs) 
                  {
                     printf("data_packer: error ... unpacking to a local grid\n");
                     exit(1);
                  }
                  if (!(g=(grid *)malloc(sizeof(grid))))
                  {
                     printf("data_packer: error ... out of memory\n");
                     exit(1);
                  }
                  g->next=g->prev=0;
                  g->dim=dim;
                  g->shape[0]=nx1=ni; g->shape[1]=ny1=nj; g->shape[2]=nz1=nk;
                  g->shape_c[0]=nx1_c=ni_c; g->shape_c[1]=ny1_c=nj_c; g->shape_c[2]=nz1_c=nk_c;
                  g->coarsest=0;
                  i1=j1=k1=1;
                  i1_c=j1_c=k1_c=1;

                  if (nk == 0) {
                    printf("warning, nk = 0: ibbox[4,5] = %d, %d \n", ibbox[4], ibbox[5]);
                  }

                  for (i=0; i<dim; i++) 
                  { 
                     g->bbox[2*i]=src->g->bbox[2*i]+dx[i]*(ibbox[2*i]-1); 
                     g->bbox[2*i+1]=src->g->bbox[2*i]+dx[i]*(ibbox[2*i+1]-1); 
                     g->ghost_width[2*i]=0; g->ghost_width[2*i+1]=0;
                     g->x[i]=0; 
                     g->x_c[i]=0; 
                  }

                  g->t=0; 
                  g->ngfs=src->g->ngfs;
                  if (!(g->gfs=(real **)malloc(sizeof(real *)*g->ngfs)))
                  {
                     printf("data_packer: error ... out of memory\n");
                     exit(1);
                  }
                  for (i=0; i<g->ngfs; i++) g->gfs[i]=0;

                  if (do_ex) 
                  {
                     chr=chr_c=0;
                     if (c->amr_mask_gfn)
                     {
                        if (!(g->gfs[c->amr_mask_gfn-1]=(real *)malloc(sizeof(real)*size)))
                        {
                           printf("data_packer: error ... out of memory \n");
                           exit(1);
                        }
                        chr=g->gfs[c->amr_mask_gfn-1];
                     }
                     if (c->amr_mask_c_gfn)
                     {
                        if (!(g->gfs[c->amr_mask_c_gfn-1]=(real *)malloc(sizeof(real)*size_c)))
                        {
                           printf("data_packer: error ... out of memory \n");
                           exit(1);
                        }
                        chr_c=g->gfs[c->amr_mask_c_gfn-1];
                     }
                     c->app_fill_ex_mask(chr,chr_c,g->dim,g->shape,g->shape_c,g->bbox,c->excised);
                  }
                  g->rank=src->g->rank;
                  src->g=g;

               }


               for (i=0; i<g->ngfs; i++)
               {
                  if (tbits[i])
                  {
                     if (dir==PACK)
                     {

                        if (!(g->gfs[i])) 
                        {
                           printf("data_packer: error in tbits[]\n");
                           exit(1);
                        }
                        if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED && size)
                        {
                          dmcopy3d_(g->gfs[i],data,&nx1,&ny1,&nz1,
                                    &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&ni,&nj,&nk,
                                    &one,&one,&one,&one,&one,&one);
                          data+=size;
                        }
                        else if (c->gfn_var_type[i]==PAMR_CELL_CENTERED && size_c)
                        {
                          dmcopy3d_(g->gfs[i],data,&nx1_c,&ny1_c,&nz1_c,
                                    &nx2_c,&ny2_c,&nz2_c,&i1_c,&j1_c,&k1_c,&i2_c,&j2_c,&k2_c,&ni_c,&nj_c,&nk_c,
                                    &one,&one,&one,&one,&one,&one);
                          data+=size_c;
                        }
                     }
                     else 
                     {
                        if (g->gfs[i])
                        {
                           printf("data_packer: bug ... trying to allocate memory twice\n");
                           exit(1);
                        }

                        if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED && size)
                        {
                           if (!(g->gfs[i]=(real *)malloc(sizeof(real)*size)))
                           {
                              printf("data_packer: error ... out of memory \n");
                              exit(1);
                           }
                           dmcopy3d_(data,g->gfs[i],&nx2,&ny2,&nz2,
                              &nx1,&ny1,&nz1,&i2,&j2,&k2,&i1,&j1,&k1,&ni,&nj,&nk,
                              &one,&one,&one,&one,&one,&one);
                           data+=size;
                        }
                        else if (c->gfn_var_type[i]==PAMR_CELL_CENTERED && size_c)
                        {
                           if (!(g->gfs[i]=(real *)malloc(sizeof(real)*size_c)))
                           {
                              printf("data_packer: error ... out of memory \n");
                              exit(1);
                           }
                           dmcopy3d_(data,g->gfs[i],&nx2_c,&ny2_c,&nz2_c,
                              &nx1_c,&ny1_c,&nz1_c,&i2_c,&j2_c,&k2_c,&i1_c,&j1_c,&k1_c,&ni_c,&nj_c,&nk_c,
                              &one,&one,&one,&one,&one,&one);
                           data+=size_c;
                        }
                     }
                  }
               }
            }

            for (i=0; i<dim; i++) 
            {
               src->bbox[2*i]=s_bbox[2*i];
               src->bbox[2*i+1]=s_bbox[2*i+1];
            }
         }
         //--------------------------------------------------------------------
         // 3. coarsen . We coarsen when packing, therefore when
         // unpacking, the new 'virtual' grid gets demoted a level, and
         // local_transfer executes it as a copy.
         //
         // align_src_bbox() forces src bbox to align to dst grid lines.
         // this should not affect the restriction, as the restriction
         // routines will apply interior stencils where possible.
         // 
         // NOTE: the ibbox_c's calculated by calc_gs_ibbox will in general have 
         // 1 extra cell volume in each direction compared to ibbox, so the 
         // subsequent bbox in the allocated grid structure only correctly
         // refers to the vertex centered functions!
         //--------------------------------------------------------------------
         else if (delta==1)
         {
            align_src_bbox(src,dst); 
            calc_gs_ibbox(src,ibbox,ibbox_c,dx);
            calc_gs_ibbox(dst,ibbox_d,ibbox_d_c,dx);
            size=size_c=1;
            for (i=0; i<dim; i++) { size*=(ibbox_d[2*i+1]-ibbox_d[2*i]+1); size_c*=(ibbox_d_c[2*i+1]-ibbox_d_c[2*i]+1); }
            g=src->g;

            dsize+=size_transfer_block(size,size_c);
            if (data) 
            {
               ni=nj=nk=i1=i2=j1=j2=k1=k2=nx1=ny1=nz1=nx2=ny2=nz2=1;
               ni_c=nj_c=nk_c=i1_c=i2_c=j1_c=j2_c=k1_c=k2_c=nx1_c=ny1_c=nz1_c=nx2_c=ny2_c=nz2_c=1; 
               switch(dim)
               {
                  case 3: 
                    k1=ibbox[4];
                    k1_c=ibbox_c[4];
                    nk=nz2=(ibbox_d[5]-ibbox_d[4]+1); 
                    nk_c=nz2_c=(ibbox_d_c[5]-ibbox_d_c[4]+1); 
                    nz1=g->shape[2]; nz1_c=g->shape_c[2];
                  case 2: 
                    j1=ibbox[2];
                    j1_c=ibbox_c[2];
                    nj=ny2=(ibbox_d[3]-ibbox_d[2]+1); 
                    nj_c=ny2_c=(ibbox_d_c[3]-ibbox_d_c[2]+1); 
                    ny1=g->shape[1]; ny1_c=g->shape_c[1];
                  case 1: 
                    i1=ibbox[0];
                    i1_c=ibbox_c[0];
                    ni=nx2=(ibbox_d[1]-ibbox_d[0]+1); 
                    ni_c=nx2_c=(ibbox_d_c[1]-ibbox_d_c[0]+1); 
                    nx1=g->shape[0]; nx1_c=g->shape_c[0];
               }

               IFG(4)
               {
                 printf("inside data packer, dir = %d, rank = %d, c->rank = %d \n", dir, rank, c->rank);
                 printf("vertex centered:    nx1=%i,ny1=%i,nz1=%i,nx2=%i,ny2=%i,nz2=%i\n",nx1,ny1,nz1,nx2,ny2,nz2);
                 printf("                    i1 =%i, j1 =%i, k1 =%i, i2=%i, j2=%i, k2=%i\n",i1,j1,k1,i2,j2,k2);
                 printf("                    ni =%i, nj =%i, nk =%i\n",ni,nj,nk);
                 printf("                    g->shape=[%i,%i,%i]\n",g->shape[0],g->shape[1],g->shape[2]);
                 printf("                    gd->shape=[%i,%i,%i]\n",dst->g->shape[0],dst->g->shape[1],dst->g->shape[2]);
                 printf("cell centered:      nx1=%i,ny1=%i,nz1=%i,nx2=%i,ny2=%i,nz2=%i\n",nx1_c,ny1_c,nz1_c,nx2_c,ny2_c,nz2_c);
                 printf("                    i1 =%i, j1 =%i, k1 =%i, i2=%i, j2=%i, k2=%i\n",i1_c,j1_c,k1_c,i2_c,j2_c,k2_c);
                 printf("                    ni =%i, nj =%i, nk =%i\n",ni_c,nj_c,nk_c);
                 printf("                    g->shape=[%i,%i,%i]\n",g->shape_c[0],g->shape_c[1],g->shape_c[2]);
                 printf("                    gd->shape=[%i,%i,%i]\n",dst->g->shape_c[0],dst->g->shape_c[1],dst->g->shape_c[2]);
               }

               if (dir==UNPACK)
               {
                  if (g->gfs) 
                  {
                     printf("data_packer: error ... unpacking to a local grid\n");
                     exit(1);
                  }
                  if (!(g=(grid *)malloc(sizeof(grid))))
                  {
                     printf("data_packer: error ... out of memory\n");
                     exit(1);
                  }
                  g->next=g->prev=0;
                  g->dim=dim;
                  g->shape[0]=nx1=ni; g->shape[1]=ny1=nj; g->shape[2]=nz1=nk;
                  g->shape_c[0]=nx1_c=ni_c; g->shape_c[1]=ny1_c=nj_c; g->shape_c[2]=nz1_c=nk_c;
                  g->coarsest=src->g->coarsest;
                  i1=j1=k1=1;
                  i1_c=j1_c=k1_c=1;

                  for (i=0; i<2*dim; i++) { g->bbox[i]=src->bbox[i]; g->ghost_width[i]=0; }
                  for (i=0; i<dim; i++) { g->x[i]=0; g->x_c[i]=0; }
                  g->t=0; 
                  g->ngfs=src->g->ngfs;
                  g->comm=src->g->comm;
                  if (!(g->gfs=(real **)malloc(sizeof(real *)*g->ngfs)))
                  {
                     printf("data_packer: error ... out of memory\n");
                     exit(1);
                  }
                  for (i=0; i<g->ngfs; i++) g->gfs[i]=0;
                  g->rank=src->g->rank;
                  src->g=g;
               }

               if (do_ex) 
               {
                  chr=chr_c=0;
                  if (g->gfs[c->amr_mask_gfn-1]) chr=g->gfs[c->amr_mask_gfn-1];
                  else if (g->gfs[c->mg_mask_gfn-1]) chr=g->gfs[c->mg_mask_gfn-1];

                  if (g->gfs[c->amr_mask_c_gfn-1]) chr_c=g->gfs[c->amr_mask_c_gfn-1];
                  else if (g->gfs[c->mg_mask_c_gfn-1]) chr_c=g->gfs[c->mg_mask_c_gfn-1];

                  if (!(chr || chr_c) && dir==PACK)
                  {
                     printf("data_packer: bug ... chr and chr_c=0\n"); exit(1);
                  }
               }

               for (i=0; i<g->ngfs; i++)
               {
                  if (tbits[i])
                  {
                     if (dir==PACK)
                     {
                        if (!(g->gfs[i])) 
                        {
                           printf("data_packer: error in tbits[]\n");
                           exit(1);
                        }

                        if (!do_ex) chr=chr_c=g->gfs[i]; // dummy grid

                        if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED)
                        {

                           if (tbits[i]==PAMR_STRAIGHT_INJECT)
                             dmcopy3d_(g->gfs[i],data,&nx1,&ny1,&nz1,
                                       &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&ni,&nj,&nk,
                                       &irho,&irho,&irho,&one,&one,&one);
                           else if (tbits[i]==PAMR_HW_RESTR)
                              dmhwr3d_(g->gfs[i],data,&nx1,&ny1,&nz1, 
                              &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&ni,&nj,&nk,
                              &irho,&irho,&irho,&one,&one,&one,
                              chr,&ex,&do_ex);
                           else if (tbits[i]==PAMR_FW_RESTR)
                              dmfwr3d_(g->gfs[i],data,&nx1,&ny1,&nz1,
                              &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&ni,&nj,&nk,
                              &irho,&irho,&irho,&one,&one,&one,
                              chr,&ex,&do_ex);
                           else if (tbits[i]==PAMR_HW_RESTR_BDY_SI) // hwr in the interior, straight inject on boundary
                                dmhwrbdy3d_(g->gfs[i],data,&nx1,&ny1,&nz1,
                                &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&ni,&nj,&nk,
                                &irho,&irho,&irho,&one,&one,&one,
                                chr,&ex,&do_ex);
                           else
                           {
                              printf("data_packer: error ... unknown injection operator\n");
                              exit(1);
                           }
                           data += size;
                        }
                        else
                          {
			    // Before packing, zero in excised zones.
			   clear_ex(g->gfs[i],chr_c,ex,nx1_c*ny1_c*nz1_c);
                           if (tbits[i]==PAMR_NN_AVERAGE)
		 	      dmavg3d_c_(g->gfs[i],data,chr_c,&nx1_c,&ny1_c,&nz1_c,
                              &nx2_c,&ny2_c,&nz2_c,&i1_c,&j1_c,&k1_c,&i2_c,&j2_c,&k2_c,&ni_c,&nj_c,&nk_c,
                              &irho,&dim,&ex,&do_ex);
                           else if (tbits[i]==PAMR_NN_ADD)
                           {
                              dmadd3d_c_(g->gfs[i],data,&nx1_c,&ny1_c,&nz1_c,
                              &nx2_c,&ny2_c,&nz2_c,&i1_c,&j1_c,&k1_c,&i2_c,&j2_c,&k2_c,&ni_c,&nj_c,&nk_c,
                              &irho,&dim);
                           }
			   else
                           {
                              printf("data_packer: error ... unknown cell centered injection operator\n");
                              exit(1);
                           }
                           data += size_c;
                        }
                     }
                     else 
                     {
                        if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED)
                        {
                           if (!(g->gfs[i]=(real *)malloc(sizeof(real)*size)))
                           {
                              printf("data_packer: error ... out of memory \n");
                              exit(1);
                           }

                           dmcopy3d_(data,g->gfs[i],&nx2,&ny2,&nz2,
                                     &nx2,&ny2,&nz2,&i2,&j2,&k2,&i1,&j1,&k1,&ni,&nj,&nk,
                                     &one,&one,&one,&one,&one,&one);
                           data+=size;
                        }
                        else
                        {
                           if (!(g->gfs[i]=(real *)malloc(sizeof(real)*size_c)))
                           {
                              printf("data_packer: error ... out of memory \n");
                              exit(1);
                           }

                           dmcopy3d_(data,g->gfs[i],&nx2_c,&ny2_c,&nz2_c,
                              &nx2_c,&ny2_c,&nz2_c,&i2_c,&j2_c,&k2_c,&i1_c,&j1_c,&k1_c,&ni_c,&nj_c,&nk_c,
                              &one,&one,&one,&one,&one,&one);
                           data+=size_c;
                        }
                     }
                  }
               }
            }
         }
         else 
         {
            printf("data_packer: restrict/interpolate only implemented between adjacted levels\n");
            exit(1);
         }
      }
      src=src->next;
      dst=dst->next;
   }

   IFG(3) printf("       data_packer(dsize=%i, data=%p, dir=%i,rank=%i, my_rank=%i) << \n",dsize,data,dir,rank,c->rank);
   return dsize;
}

void free_src_data(gsl *src)
{
   int i;

   IFG(3) printf("    >> free_src_data \n");

   while(src)
   {
      if (src->g->rank != curr_context->rank && src->g->gfs) 
      {
         for (i=0; i<src->g->ngfs; i++) if (src->g->gfs[i]) free(src->g->gfs[i]);
         free(src->g->gfs);
         for (i=0; i<curr_context->dim; i++) if (src->g->x[i]) free(src->g->x[i]);
         free(src->g);
         src->g=0;
      }
      src=src->next;
   }

   IFG(3) printf("       free_src_data << \n");
}

//=============================================================================
// local_transfer executes all the copy commands implied by the (src,dst) pair,
// for all dst grid-segments with rank=my_rank()
//
// The internals of this function are very similar to pack/unpack
//
// NOTE: interpolation is currently *not* efficient for myrank->myrank
//       transfers ... should be easy to optimize if necessary (interpolates
//       entire grid then copies desired segment)
//=============================================================================
void local_transfer(gsl *src, gsl *dst)
{
   int dim=curr_context->dim,i,rank=curr_context->rank,j;
   context *c=curr_context;
   int *tbits=curr_context->tf_bits;
   grid *g1,*g2;
   int one=1,delta,ibbox[2*PAMR_MAX_DIM],ibbox_d[2*PAMR_MAX_DIM];
   int ibbox_c[2*PAMR_MAX_DIM],ibbox_d_c[2*PAMR_MAX_DIM];
   int i1i,i2i,nii,j1i,j2i,nji,k1i,k2i,nki;
   int nx1,ny1,nz1,nx2,ny2,nz2,i1,j1,k1,i2,j2,k2,ni,nj,nk,irho;
   int nib,njb,nkb,nx1b,ny1b,nz1b,i2b,j2b,k2b;
   int i1i_c,i2i_c,nii_c,j1i_c,j2i_c,nji_c,k1i_c,k2i_c,nki_c;
   int nx1_c,ny1_c,nz1_c,nx2_c,ny2_c,nz2_c,i1_c,j1_c,k1_c,i2_c,j2_c,k2_c,ni_c,nj_c,nk_c;
   int nib_c,njb_c,nkb_c,nx1b_c,ny1b_c,nz1b_c,i2b_c,j2b_c,k2b_c;
   real *data,dx[2*PAMR_MAX_DIM],dxd[2*PAMR_MAX_DIM];
   real *chr1,*chr2,ex;
   real *chr1_c,*chr2_c;
   int do_ex;

   IFG(3) printf("  >> local_transfer() rank=%i \n",rank);

   do_ex=0;
   if (c->excision_on)
   {
      ex=c->excised;
      do_ex=1;
   }

   while(src && dst)
   {
      if (dst->g->rank==rank)
      {
         delta=calc_delta(src,dst,&irho);
         if (delta==1) align_src_bbox(src,dst); 

         // Branson:  If we first align the source bbox with
         // it's own underlying grid lines then interpolation
         // should work even without if iterp_buffer = 0.
         if (delta==-1) align_src_bbox_interp(src); 

         calc_gs_ibbox(src,ibbox,ibbox_c,dx);
         calc_gs_ibbox(dst,ibbox_d,ibbox_d_c,dxd);

         g1=src->g; 
         g2=dst->g; 
         if (do_ex)
         {
            if (!(chr1=g1->gfs[c->amr_mask_gfn-1])) chr1=g1->gfs[c->mg_mask_gfn-1];
            if (!(chr2=g2->gfs[c->amr_mask_gfn-1])) chr2=g2->gfs[c->mg_mask_gfn-1];
            if (!(chr1_c=g1->gfs[c->amr_mask_c_gfn-1])) chr1_c=g1->gfs[c->mg_mask_c_gfn-1];
            if (!(chr2_c=g2->gfs[c->amr_mask_c_gfn-1])) chr2_c=g2->gfs[c->mg_mask_c_gfn-1];
            if ((!(chr1 || chr1_c) || !(chr2 || chr2_c)) && delta!=0)
            {
               printf("chr1=%p, chr2=%p, amr_mask_gfn=%i, delta=%i\n",chr1,chr2,c->amr_mask_gfn,delta);
               printf("chr1_c=%p, chr2_c=%p, amr_mask_c_gfn=%i, delta=%i\n",chr1_c,chr2_c,c->amr_mask_c_gfn,delta);
               printf("local_transfer: bug ... chr[_c|1|2]=0\n"); exit(1);
            }
         }
         ni=nj=nk=i1=i2=j1=j2=k1=k2=nx1=ny1=nz1=nx2=ny2=nz2=1;
         nib=njb=nkb=nx1b=ny1b=nz1b=i2b=j2b=k2b=1;
         nii=nji=nki=i1i=j1i=k1i=i2=j2i=k2i=1;
         ni_c=nj_c=nk_c=i1_c=i2_c=j1_c=j2_c=k1_c=k2_c=nx1_c=ny1_c=nz1_c=nx2_c=ny2_c=nz2_c=1;
         nib_c=njb_c=nkb_c=nx1b_c=ny1b_c=nz1b_c=i2b_c=j2b_c=k2b_c=1;
         nii_c=nji_c=nki_c=i1i_c=j1i_c=k1i_c=i2_c=j2i_c=k2i_c=1;
         switch(dim)
         {
           case 3: 
             k1=ibbox[4]; k1_c=ibbox_c[4];
             k2=ibbox_d[4]; k2_c=ibbox_d_c[4];
             nk=(ibbox[5]-ibbox[4]+1); 
             nk_c=(ibbox_c[5]-ibbox_c[4]+1); 
             nz1=g1->shape[2]; nz1_c=g1->shape_c[2];
             nz2=g2->shape[2]; nz2_c=g2->shape_c[2];
             nkb=(ibbox_d[5]-ibbox_d[4]+1);
             nkb_c=(ibbox_d_c[5]-ibbox_d_c[4]+1);
             if (delta==-1)
               { 
                 k2b=k2b_c=(src->bbox[4]-g1->bbox[4])/dxd[2]+1.5;  //dxd for interpolation
                 k2i=k2i_c=max(1,(((k2b-1)/irho)-c->interp_buffer)*irho+1);
                 k1i=k1i_c=(k2i-1)/irho+1;
                 nki=min(nz1,(k2b+nkb-1)/irho+1+c->interp_buffer+1)-k1i+1; nki_c=nki-1;
                 nz1b=(nki-1)*irho+1; nz1b_c=nz1b-1;
                 
                 // Since we allow different source and destination
                 // bboxes for interpolation, we need to recalculate k2b.
                 k2b=k2b_c=(dst->bbox[4]-g1->bbox[4])/dxd[2]+1.5;  
                 k2b=k2b_c=k2b-k2i+1;
               }
            case 2: 
              j1=ibbox[2]; j1_c=ibbox_c[2];
              j2=ibbox_d[2]; j2_c=ibbox_d_c[2];
              nj=(ibbox[3]-ibbox[2]+1); 
              nj_c=(ibbox_c[3]-ibbox_c[2]+1); 
              ny1=g1->shape[1]; ny1_c=g1->shape_c[1];
              ny2=g2->shape[1]; ny2_c=g2->shape_c[1];
              njb=(ibbox_d[3]-ibbox_d[2]+1); 
              njb_c=(ibbox_d_c[3]-ibbox_d_c[2]+1); 
              if (delta==-1)
                {
                  j2b=j2b_c=(src->bbox[2]-g1->bbox[2])/dxd[1]+1.5;  //dxd for interpolation
                  j2i=j2i_c=max(1,(((j2b-1)/irho)-c->interp_buffer)*irho+1);
                  j1i=j1i_c=(j2i-1)/irho+1;
                  nji=min(ny1,(j2b+njb-1)/irho+1+c->interp_buffer+1)-j1i+1; nji_c=nji-1;
                  ny1b=(nji-1)*irho+1; ny1b_c=ny1b-1;
                  
                  j2b=j2b_c=(dst->bbox[2]-g1->bbox[2])/dxd[1]+1.5;  
                  j2b=j2b_c=j2b-j2i+1;
                }
            case 1: 
              i1=ibbox[0]; i1_c=ibbox_c[0];
              i2=ibbox_d[0]; i2_c=ibbox_d_c[0];
              ni=(ibbox[1]-ibbox[0]+1); 
              ni_c=(ibbox_c[1]-ibbox_c[0]+1); 
              nx1=g1->shape[0]; nx1_c=g1->shape_c[0];
              nx2=g2->shape[0]; nx2_c=g2->shape_c[0];
              nib=(ibbox_d[1]-ibbox_d[0]+1); 
              nib_c=(ibbox_d_c[1]-ibbox_d_c[0]+1); 
              if (delta==-1)
                {
                  i2b=i2b_c=(src->bbox[0]-g1->bbox[0])/dxd[0]+1.5;  //dxd for interpolation
                  i2i=i2i_c=max(1,(((i2b-1)/irho)-c->interp_buffer)*irho+1);
                  i1i=i1i_c=(i2i-1)/irho+1;
                  nii=min(nx1,(i2b+nib-1)/irho+1+c->interp_buffer+1)-i1i+1; nii_c=nii-1;
                  nx1b=(nii-1)*irho+1; nx1b_c=nx1b-1;

                  i2b=i2b_c=(dst->bbox[0]-g1->bbox[0])/dxd[0]+1.5;  
                  i2b=i2b_c=i2b-i2i+1;
                }
         }
         for (i=0; i<g1->ngfs; i++)
         {
            if (tbits[i])
            {

              // You could have the situation where you only want to transfer vertices--no cells.
              if ((c->gfn_var_type[i]==PAMR_VERTEX_CENTERED && (ni>0) && (nj>0) && (nk>0)) ||
                  (c->gfn_var_type[i]==PAMR_CELL_CENTERED && (ni_c>0) && (nj_c>0) && (nk_c>0)) ) {

               if (!(g1->gfs[i] && g2->gfs[i]))
               {
                  printf("transfer: error in tbits[]\n");
                  printf("\nlocal tranfer: nx1=%i,ny1=%i,nz1=%i,nx2=%i,ny2=%i,nz2=%i\n",nx1,ny1,nz1,nx2,ny2,nz2);
                  printf("               nx1b=%i,ny1b=%i,nz1b=%i\n",nx1b,ny1b,nz1b);
                  printf("               i1 =%i, j1 =%i, k1 =%i, i2=%i, j2=%i, k2=%i\n",i1,j1,k1,i2,j2,k2);
                  printf("               i1i=%i, j1i=%i, k1i=%i, i2i=%i,j2i=%i,k2i=%i\n",i1i,j1i,k1i,i2i,j2i,k2i);
                  printf("               ni =%i, nj =%i, nk =%i\n",ni,nj,nk);
                  printf("               ni_c =%i, nj_c =%i, nk_c =%i\n",ni_c,nj_c,nk_c);
                  printf("               nii =%i,nji =%i,nki =%i\n",nii,nji,nki);
                  printf("               nii =%i,nji =%i,nki =%i\n",nii,nji,nki);
                  printf("               nib=%i, njb=%i, nkb=%i\n",nib,njb,nkb);
                  printf("               nib_c=%i, njb_c=%i, nkb_c=%i\n",nib_c,njb_c,nkb_c);
                  printf("               g1->shape=[%i,%i,%i]\n",g1->shape[0],g1->shape[1],g1->shape[2]);
                  printf("               g2->shape=[%i,%i,%i]\n",g2->shape[0],g2->shape[1],g2->shape[2]);
                  printf("               g1->shape_c=[%i,%i,%i]\n",g1->shape_c[0],g1->shape_c[1],g1->shape_c[2]);
                  printf("               g2->shape_c=[%i,%i,%i]\n",g2->shape_c[0],g2->shape_c[1],g2->shape_c[2]);
                  printf("               g1->bbox=[%lf,%lf,%lf,%lf,%lf,%lf]\n",g1->bbox[0],g1->bbox[1],g1->bbox[2],g1->bbox[3],g1->bbox[4],g1->bbox[5]);
                  printf("               g2->bbox=[%lf,%lf,%lf,%lf,%lf,%lf]\n",g2->bbox[0],g2->bbox[1],g2->bbox[2],g2->bbox[3],g2->bbox[4],g2->bbox[5]);
                  printf("               g1->rank=%i, gfn=%i, irho=%i\n",g1->rank,i+1,irho);
                  printf("               g2->rank=%i, gfn=%i, irho=%i\n",g2->rank,i+1,irho);
                  exit(1);
               }
               IFG(4)
               {
                  if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED)
                  {
                     printf("\nVC local tranfer: nx1=%i,ny1=%i,nz1=%i,nx2=%i,ny2=%i,nz2=%i\n",nx1,ny1,nz1,nx2,ny2,nz2);
                     printf("               nx1b=%i,ny1b=%i,nz1b=%i\n",nx1b,ny1b,nz1b);
                     printf("               i1 =%i, j1 =%i, k1 =%i, i2=%i, j2=%i, k2=%i\n",i1,j1,k1,i2,j2,k2);
                     printf("               i1i=%i, j1i=%i, k1i=%i, i2i=%i,j2i=%i,k2i=%i\n",i1i,j1i,k1i,i2i,j2i,k2i);
                     printf("               ni =%i, nj =%i, nk =%i\n",ni,nj,nk);
                     printf("               nii =%i,nji =%i,nki =%i\n",nii,nji,nki);
                     printf("               nib=%i, njb=%i, nkb=%i\n",nib,njb,nkb);
                     printf("               g1->shape=[%i,%i,%i]\n",g1->shape[0],g1->shape[1],g1->shape[2]);
                     printf("               g2->shape=[%i,%i,%i]\n",g2->shape[0],g2->shape[1],g2->shape[2]);
                     printf("               g1->rank=%i, gfn=%i, irho=%i\n",g1->rank,i+1,irho);
                  }
                  else
                  {
                     printf("\nCC local tranfer: nx1=%i,ny1=%i,nz1=%i,nx2=%i,ny2=%i,nz2=%i\n",nx1_c,ny1_c,nz1_c,nx2_c,ny2_c,nz2_c);
                     printf("               nx1b=%i,ny1b=%i,nz1b=%i\n",nx1b_c,ny1b_c,nz1b_c);
                     printf("               i1 =%i, j1 =%i, k1 =%i, i2=%i, j2=%i, k2=%i\n",i1_c,j1_c,k1_c,i2_c,j2_c,k2_c);
                     printf("               i1i=%i, j1i=%i, k1i=%i, i2i=%i,j2i=%i,k2i=%i\n",i1i_c,j1i_c,k1i_c,i2i_c,j2i_c,k2i_c);
                     printf("               ni =%i, nj =%i, nk =%i\n",ni_c,nj_c,nk_c);
                     printf("               nii =%i,nji =%i,nki =%i\n",nii_c,nji_c,nki_c);
                     printf("               nib=%i, njb=%i, nkb=%i\n",nib_c,njb_c,nkb_c);
                     printf("               g1->shape=[%i,%i,%i]\n",g1->shape_c[0],g1->shape_c[1],g1->shape_c[2]);
                     printf("               g2->shape=[%i,%i,%i]\n",g2->shape_c[0],g2->shape_c[1],g2->shape_c[2]);
                     printf("               g1->rank=%i, gfn=%i, irho=%i\n",g1->rank,i+1,irho);
                     printf("               g1->gfs[i]=%p, g2->gfs[i]=%p\n",g1->gfs[i],g2->gfs[i]);
                  }
               }

               //--------------------------------------------------------------------
               // 1. copy 
               //--------------------------------------------------------------------
               if (delta==0)
               {
                 if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED) {
                   dmcopy3d_(g1->gfs[i],g2->gfs[i],&nx1,&ny1,&nz1,
                             &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&ni,&nj,&nk,
                             &one,&one,&one,&one,&one,&one);
                 } else {
                   dmcopy3d_(g1->gfs[i],g2->gfs[i],&nx1_c,&ny1_c,&nz1_c,
                             &nx2_c,&ny2_c,&nz2_c,&i1_c,&j1_c,&k1_c,&i2_c,&j2_c,&k2_c,&ni_c,&nj_c,&nk_c,
                             &one,&one,&one,&one,&one,&one);
                 }
               }
               //--------------------------------------------------------------------
               // 2. interpolate
               //
               // To keep the interpolation routine simple, we require that
               // coarse/fine boundaries be aligned. Also, we want to use the
               // interior stencils where-ever possible. So first interpolate a 
               // larger region of temporary storage, then copy the relevant
               // part to the actual grid.
               //--------------------------------------------------------------------
               else if (delta==-1)
               {
                  if (tbits[i]>0)
                  {
                     if (!(data=(real *)malloc(sizeof(real)*nx1b*ny1b*nz1b)))
                     {
                        printf("local_transfer: out of memory ... cannot interolate\n");
                     }
                     else
                     {
                        if (!do_ex) chr1=chr1_c=g1->gfs[i]; //dummy

                        if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED)
                        {

                          dminterp3d_(g1->gfs[i],data,&nx1,&ny1,&nz1,
                                      &nx1b,&ny1b,&nz1b,&i1i,&j1i,&k1i,&one,&one,&one,&nii,&nji,&nki,
                                      &one,&one,&one,&irho,&irho,&irho,&tbits[i],
                                      chr1,&ex,&do_ex);
                          dmcopy3d_(data,g2->gfs[i],&nx1b,&ny1b,&nz1b,
                                    &nx2,&ny2,&nz2,&i2b,&j2b,&k2b,&i2,&j2,&k2,&nib,&njb,&nkb,
                                    &one,&one,&one,&one,&one,&one);
                          
                          if (do_ex)
                            {
                              clear_ex(g1->gfs[i],chr1,ex,nx1*ny1*nz1);
                              clear_ex(g2->gfs[i],chr2,ex,nx2*ny2*nz2);
                            }
                        }
                        else
                        { 
			  dminterp3d_c_(g1->gfs[i],data,&nx1_c,&ny1_c,&nz1_c,
					&nx1b_c,&ny1b_c,&nz1b_c,&i1i_c,&j1i_c,&k1i_c,&one,&one,&one,&nii_c,&nji_c,&nki_c,
					&one,&one,&one,&irho,&irho,&irho,&tbits[i],
					chr1_c,&ex,&do_ex,&dim);
			  
			  dmcopy3d_(data,g2->gfs[i],&nx1b_c,&ny1b_c,&nz1b_c,
				    &nx2_c,&ny2_c,&nz2_c,&i2b_c,&j2b_c,&k2b_c,&i2_c,&j2_c,&k2_c,&nib_c,&njb_c,&nkb_c,
				    &one,&one,&one,&one,&one,&one);
                          if (do_ex && c->amr_mask_c_gfn)
                            {
                              clear_ex(g1->gfs[i],chr1_c,ex,nx1_c*ny1_c*nz1_c);
                              clear_ex(g2->gfs[i],chr2_c,ex,nx2_c*ny2_c*nz2_c);
                            }
                        }
                        
                        free(data);
                     }
                  }
                  else
                  {
                     printf("local_transfer: unsupported interpolation order=%i\n",tbits[i]);
                     exit(1);
                  }
               }
               //--------------------------------------------------------------------
               // 3. coarsen
               //--------------------------------------------------------------------
               else if (delta==1)
               {
                  if (!do_ex) chr1=chr1_c=g1->gfs[i]; //dummy

                  if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED)
                  {
                     if (tbits[i]==PAMR_STRAIGHT_INJECT)
                        dmcopy3d_(g1->gfs[i],g2->gfs[i],&nx1,&ny1,&nz1,
                                  &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&nib,&njb,&nkb,
                                  &irho,&irho,&irho,&one,&one,&one);
                     else if (tbits[i]==PAMR_HW_RESTR)
                        dmhwr3d_(g1->gfs[i],g2->gfs[i],&nx1,&ny1,&nz1, 
                                 &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&nib,&njb,&nkb,
                                 &irho,&irho,&irho,&one,&one,&one,
                                 chr1,&ex,&do_ex);
                     else if (tbits[i]==PAMR_FW_RESTR)
                        dmfwr3d_(g1->gfs[i],g2->gfs[i],&nx1,&ny1,&nz1,
                                 &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&nib,&njb,&nkb,
                                 &irho,&irho,&irho,&one,&one,&one,
                                 chr1,&ex,&do_ex);
                     else if (tbits[i]==PAMR_HW_RESTR_BDY_SI)
                         dmhwrbdy3d_(g1->gfs[i],g2->gfs[i],&nx1,&ny1,&nz1,
                                  &nx2,&ny2,&nz2,&i1,&j1,&k1,&i2,&j2,&k2,&nib,&njb,&nkb,
                                  &irho,&irho,&irho,&one,&one,&one,
                                  chr1,&ex,&do_ex);
                     else
                     {
                        printf("local_transfer: error ... unknown injection operator\n");
                        exit(1);
                     }
                  }
                  else
                  {
		    // First, zero the gf at any excised zones.
		     if (do_ex && c->amr_mask_c_gfn)
		       {
			 clear_ex(g1->gfs[i],chr1_c,ex,nx1_c*ny1_c*nz1_c);
		       }
                     if (tbits[i]==PAMR_NN_AVERAGE)
                        dmavg3d_c_(g1->gfs[i],g2->gfs[i],chr1_c,&nx1_c,&ny1_c,&nz1_c,
                                  &nx2_c,&ny2_c,&nz2_c,&i1_c,&j1_c,&k1_c,&i2_c,&j2_c,&k2_c,&nib_c,&njb_c,&nkb_c,
                                  &irho,&dim,&ex,&do_ex);
                     else if (tbits[i]==PAMR_NN_ADD)
                     {
                        // printf("TEST: not doing NN_ADD B\n");
                        dmadd3d_c_(g1->gfs[i],g2->gfs[i],&nx1_c,&ny1_c,&nz1_c,
                                  &nx2_c,&ny2_c,&nz2_c,&i1_c,&j1_c,&k1_c,&i2_c,&j2_c,&k2_c,&nib_c,&njb_c,&nkb_c,
                                  &irho,&dim);
                     }
		     else
                     {
                        printf("local_transfer: error ... unknown cell-based injection operator\n");
                        exit(1);
                     }
		     // Now zero inside the excised zones on the coarse grid.
		     if (do_ex && c->amr_mask_c_gfn)
		       {
			 clear_ex(g2->gfs[i],chr2_c,ex,nx2_c*ny2_c*nz2_c);
		       }
                  }

               }
               else 
               {
                  printf("local_transfer: cannot coarsen/interpolate between non-adjacent levels\n");
                  exit(1);
               }
              }
           }
         }
      }
      
      src=src->next;
      dst=dst->next;
   }

   IFG(3) printf("     local_transfer << \n");
}

//=============================================================================
//  transfer() communicates source grid segments to destination grid segments.
//  
//  This function expects the transfer lists from *all* nodes (hence the 
//  array arguments) 
//
//  NOTE: this function will alter elements of the src lists (for (rank!=my_rank())
//
//  To avoid possible latency problems, we also send all of the segments
//  from the current node to a given destination in one go.
//=============================================================================
#define ONE_GIG 1000000000
#define ONE_MEG 1000000
void transfer(gsl **src, gsl **dst)
{
   context *c=curr_context;
   int i,dsize,ret;
   real *send_data[PAMR_MAX_NODES],*rec_data[PAMR_MAX_NODES];
   gsl *gs;
   MPI_Request *reqs;
   MPI_Status *stats;
   int req_no=0,num_rec=0,num_snd=0;
   static int call=0;
   int ltrace=0;
   int prev_clock,curr_clock,GB;

   call++;

   IFGL(2) printf("  >> transfer() rank=%i,size=%i \n",c->rank,c->size);

   if (!(reqs=(MPI_Request *)malloc(sizeof(MPI_Request)*c->size*2)))
   {
      printf("transfer: out of memory\n");
      exit(1);
   }
   if (!(stats=(MPI_Status *)malloc(sizeof(MPI_Status)*c->size*2)))
   {
      printf("transfer: out of memory\n");
      exit(1);
   }

   for (i=0; i<c->size; i++)
   {
      IFGL(2) { printf("syncing i=%i ... my_rank=%i\n",i,c->rank); MPI_Barrier(MPI_COMM_WORLD); }
         
      send_data[i]=0;
      rec_data[i]=0;

      //IFGL(2) { if ((c->rank==83) && i>=17)  { printf(" turing on more trace for i=83\n"); gtrace=4; }  }
      //IFGL(2) { if ((c->rank==83 || c->rank==17) && i>=17)  { printf(" turing on more trace for i=17,83\n"); gtrace=4; }  }

      // send from this node to node i
      if (i==c->rank) local_transfer(src[i],dst[i]);
      else if (dsize=data_packer(0,src[c->rank],dst[c->rank],i,PACK))
      {
         if (!(send_data[i]=(real *)malloc(sizeof(real)*dsize)))
         {
            printf("transfer: error ... out of memory\n");
            exit(1);
         }

         IFGL(4) printf("  my_rank=%i, allocated block of size %i, pointer=%p\n",c->rank,dsize,send_data[i]);

         data_packer(send_data[i],src[c->rank],dst[c->rank],i,PACK);
         // transmit it 
         IFGL(4) printf("\t\t MPI_Isend from %i to %i, dsize=%i ... \n",c->rank,i,dsize);
         if (PAMR_stats_on)
         {
            prev_clock=clock();
            GB=((sizeof(double)*dsize)/ONE_GIG);
            PAMR_comm_s_GB+=GB;
            PAMR_comm_s_B+=((sizeof(double)*dsize)-GB*ONE_GIG);
            if (PAMR_comm_s_B>=ONE_GIG) {PAMR_comm_s_GB++; PAMR_comm_s_B-=ONE_GIG;}
            PAMR_num_s++;
         }
         ret=MPI_Isend((void *)send_data[i],dsize,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&reqs[req_no++]);
         if (PAMR_stats_on)
         {
            curr_clock=clock();
            if (curr_clock<prev_clock) // skip the occasional wrap-around
            {
               PAMR_comm_s_microsecs+=(((curr_clock-prev_clock)/((real)CLOCKS_PER_SEC))*ONE_MEG);
               if (PAMR_comm_s_microsecs>=ONE_MEG) {PAMR_comm_s_secs++; PAMR_comm_s_microsecs-=ONE_MEG;}
            }
         }

         num_snd++;
      }

      IFGL(2) { printf("syncing B i=%i ... my_rank=%i\n",i,c->rank); MPI_Barrier(MPI_COMM_WORLD); }

      // receive data from node i to this node
      if ((i!=c->rank) && (dsize=data_packer(0,src[i],dst[i],i,UNPACK)))
      {
         if (!(rec_data[i]=(real *)malloc(sizeof(real)*dsize)))
         {
            printf("transfer: error ... out of memory\n");
            exit(1);
         }
         // get it 
         IFGL(4) printf("\t\t MPI_Irecv at %i from %i, dsize=%i ... \n",c->rank,i,dsize);
         if (PAMR_stats_on)
         {
            prev_clock=clock();
            GB=((sizeof(MPI_DOUBLE)*dsize)/ONE_GIG);
            PAMR_comm_r_GB+=GB;
            PAMR_comm_r_B+=((sizeof(MPI_DOUBLE)*dsize)-GB*ONE_GIG);
            if (PAMR_comm_r_B>=ONE_GIG) {PAMR_comm_r_GB++; PAMR_comm_r_B-=ONE_GIG;}
            PAMR_num_r++;
         }
         ret=MPI_Irecv((void *)rec_data[i],dsize,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&reqs[req_no++]);
         if (PAMR_stats_on)
         {
            curr_clock=clock();
            if (curr_clock<prev_clock) // skip the occasional wrap-around
            {
               PAMR_comm_r_microsecs+=(((curr_clock-prev_clock)/((real)CLOCKS_PER_SEC))*ONE_MEG);
               if (PAMR_comm_r_microsecs>=ONE_MEG) {PAMR_comm_r_secs++; PAMR_comm_r_microsecs-=ONE_MEG;}
            }
         }
         num_rec++;
      }
   }

   // wait for all requests to complete
   IFGL(4) printf("\t\t before Waitall, rank=%i, num_snd=%i, num_rec=%i, call=%i\n",c->rank,num_snd,num_rec,call);
   ret=MPI_Waitall(req_no,reqs,stats);
   if (ret)
   {
      printf("MPI_Waitall ... status array, node %i, ret=%i:\n",c->rank,ret);
      for (i=0; i<req_no; i++) printf("req %i: MPI_SOURCE=%i MPI_TAG=%i MPI_ERROR=%i\n",
                                      i,stats[i].MPI_SOURCE,stats[i].MPI_TAG,stats[i].MPI_ERROR);
      exit(1);
   }
   IFGL(4) printf("\t\t after Waitall, rank=%i, num_snd=%i, num_rec=%i, call=%i\n",c->rank,num_snd,num_rec,call);

   //now transfer all local data
   for (i=0; i<c->size; i++)
   {
      if (rec_data[i]) 
      {
 	 data_packer(rec_data[i],src[i],dst[i],i,UNPACK);
         local_transfer(src[i],dst[i]);
         free(rec_data[i]);
      }
      if (send_data[i]) free(send_data[i]);
      if (i!=c->rank) free_src_data(src[i]);
   }

   free(reqs);
   free(stats);

   IFGL(2) printf("  << transfer() rank=%i\n",c->rank);
}

//=============================================================================
// copy_cv_to_vc copies cell(vertex) centered grid functions to their
// vertex(cell) centered counterparts. The grid functions involved are
// those with their tf bits set.
//
// if curr_g_only, then assumes a valid iterator has already been set
// up, and only does the operation for the current grid
// In this case, lev, hierarchy and AMR_bdy_width are not ignored
//=============================================================================
int copy_cv_to_vc(int lev, int hierarchy, int which_tf, int AMR_bdy_width, int curr_g_only)
{
   context *c=curr_context;
   grid *g;
   int valid,i,do_ex;
   int conj_gfns[PAMR_MAX_GFNS],ret=1;
   real *chr1,*chr1_c,ex;

   fill_conjugates(conj_gfns);

   do_ex=0;
   if (c->excision_on)
   {
      ex=c->excised;
      do_ex=1;
   }

   if (curr_g_only)
   {
      if (c->iter_g && c->iter_g->rank==c->rank) valid=1; else valid=0;
   }
   else
   {
      PAMR_push_iter();
      valid=PAMR_init_s_iter(lev,hierarchy,0);
   }

   while(valid)
   {
      g=c->iter_g;
      // step to next grid already, so we know if we're at the last one,
      // to set the PAMR_SYNC flag
      if (curr_g_only) valid=0; else valid=PAMR_next_g(); 

      if (do_ex)
      {
         if (!(chr1=g->gfs[c->amr_mask_gfn-1])) chr1=g->gfs[c->mg_mask_gfn-1];
         if (!(chr1_c=g->gfs[c->amr_mask_c_gfn-1])) chr1_c=g->gfs[c->mg_mask_c_gfn-1];
      }
      for (i=0; i<c->gfns; i++)
      {
         if (g->gfs[i])
         {
            if (!do_ex) chr1=chr1_c=g->gfs[i]; //dummy
            if (c->tf_bits[i] && which_tf==PAMR_TF_C_TO_V && c->gfn_var_type[i]==PAMR_CELL_CENTERED)
            {
               if (!conj_gfns[i]) { printf("copy_vc_to_cv : error ... c->v conjugate variable undefined"); PAMR_pop_iter(); return 0; }

               dm_c_to_v_(g->gfs[i],g->gfs[conj_gfns[i]-1],chr1_c,chr1,
                          &g->shape_c[0],&g->shape_c[1],&g->shape_c[2],&g->shape[0],&g->shape[1],&g->shape[2],
                          &c->tf_bits[i],&ex,&do_ex,&g->dim);

               if (!valid && !curr_g_only)
               {
                  c->tf_bits[i]=PAMR_NO_SYNC; 
                  c->tf_bits[conj_gfns[i]-1]=PAMR_SYNC; 
               }
            }
            if (c->tf_bits[i] && which_tf==PAMR_TF_V_TO_C && c->gfn_var_type[i]==PAMR_VERTEX_CENTERED)
            {
	       if (!conj_gfns[i]) {printf("copy_vc_to_cv : error ... v->c conjugate variable undefined"); PAMR_pop_iter(); return 0; }
       
               dm_v_to_c_(g->gfs[i],g->gfs[conj_gfns[i]-1],chr1,chr1_c,
                          &g->shape[0],&g->shape[1],&g->shape[2],&g->shape_c[0],&g->shape_c[1],&g->shape_c[2],
                          &c->tf_bits[i],&ex,&do_ex,&g->dim);

               if (!valid && !curr_g_only)
               {
                  c->tf_bits[i]=PAMR_NO_SYNC; 
                  c->tf_bits[conj_gfns[i]-1]=PAMR_SYNC; 
               }
            }
         }
      }
   }

   if (!curr_g_only)
   {
      PAMR_pop_iter();

      ret=communicate(lev,lev,hierarchy,1+AMR_bdy_width,0,0,"sync");
   }
  
   return ret;
}
