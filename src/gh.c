//=============================================================================
//
//  gh.c --- grid hierarhcy functions
//
//  note: all of these routines are 'internal', in that they assume
//  a valid context exists
//
//  grid function id's, level numbers, time-levels and similar numbers 
//  communicated to users are numbered from 1
//
//=============================================================================

#include "gh.h"
#include "pamr.h"
#include "misc.h"
#include "transfer.h"
#include "io.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

context *contexts[PAMR_MAX_CONTEXTS];
context *curr_context=0;

int PAMR_parallel_io=1;

int PAMR_stats_on=0; 
int PAMR_comm_r_secs=0,PAMR_comm_r_microsecs=0; 
int PAMR_comm_r_GB=0,PAMR_comm_r_B=0;
int PAMR_comm_s_secs=0,PAMR_comm_s_microsecs=0; 
int PAMR_comm_s_GB=0,PAMR_comm_s_B=0;
int PAMR_num_s=0,PAMR_num_r=0;

/* the following globals are used by set_var_attribs() and define_var_brief() */
 
int c_in_amrh;
int c_in_mgh;
int c_num_tl;
int c_amr_inject;
int c_amr_interp;
int c_amr_bdy_interp;
int c_amr_sync;
int c_mg_inject;
int c_mg_interp;
int c_mg_sync;
int c_mg_noinj_to_amr;
int c_regrid_transfer;
int c_phys_bdy_type[2*PAMR_MAX_DIM];

//=============================================================================
// calculates dx,dt level info for the current context
// (internal routine)
//=============================================================================
void recalc_lev_info(void)
{
   int i,j,ltrace=0;
   context *cc=curr_context;

   IFG(4) printf("  >> recalc_lev_info \n");

   for (i=0; i<cc->dim; i++)
   {
      IFG(4) printf("dim=%i\nlevel\tdx\tdt\n",i);
      j=0;
      if (cc->shape[i]>1) cc->dx[j][i]=(cc->bbox[2*i+1]-cc->bbox[2*i])/(cc->shape[i]-1);
      else { cc->dx[j][i]=PAMR_SMALL_DX; printf("recalc_lev_info: WARNING ... cc->dx[j][i]~0\n"); }
      if (i==0 || (cc->dt[j]>cc->lambda*cc->dx[j][i])) cc->dt[j]=cc->lambda*cc->dx[j][i];
      while(j<(PAMR_MAX_LEVS-1))
      {
         IFG(4) printf("%i\t%lf\t%lf\n",j,cc->dx[j][i],cc->dt[j]);
         j++;
         cc->dx[j][i]=cc->dx[j-1][i]/cc->rho_sp[j-1];
         cc->dt[j]=cc->dt[j-1]/cc->rho_tm[j-1];
      }
   }

   IFG(4) printf("     recalc_lev_info << \n");
}

//=============================================================================
// allocates a new variable and updates the context accordingly 
//=============================================================================
var *new_var(char *name)
{
   var *f;
   char *new_name=0;

   IFG(4) printf("  >> new_var \n");

   f=find_var(name); if (f) return 0;

   f=(struct var*)malloc(sizeof(struct var));
   new_name=strdup(name);
   if (!f || !new_name)
   { 
      printf("new_var: error ... out of memory f=%i, new_name=%i\n",f,new_name); 
      if (f) free(f);
      if (new_name) free(new_name);
      return 0;
   }

   f->name=new_name;
   f->next=curr_context->vars;
   curr_context->vars=f;
   curr_context->num_vars++;
   IFG(4) printf("     new_var << \n");
   return f;
}

//=============================================================================
// returns the var structure corresponding to name, or 0 if not found
//=============================================================================
var *find_var(char *name)
{
   var *f;

   IFG(4) printf("  >> find_var    \n");

   f=curr_context->vars;
   while(f)
   {
      if (!(strcmp(f->name,name))) goto fin;
      f=f->next;
   }

fin:
   IFG(4) printf("     find_var << \n");
   return f;
}

//=============================================================================
// deletes all memory associated with a context (context fields are NOT reset)
//=============================================================================
void free_context_mem()
{
   var *f,*q;
   int l;

   IFG(4) printf("   >> free_context_mem\n");

   if (curr_context->MG_cgh) destroy_mgh();
   if (curr_context->tf_bits) free(curr_context->tf_bits);
   if (curr_context->curr_cgh) free_cgh(curr_context->curr_cgh,1,PAMR_MAX_LEVS);

   f=curr_context->vars;
   while(f)
   {
      q=f;
      f=f->next;
      free(q->name);
      free(q);
   }

   for (l=0; l<PAMR_MAX_LEVS; l++) 
      if (curr_context->sgh_bboxes[l]) free(curr_context->sgh_bboxes[l]);

   IFG(4) printf("      free_context_mem << \n");
}  

//=============================================================================
// utility routines 
//=============================================================================

//=============================================================================
// returns 1 if given grid function is in the AMR hierarchy
//=============================================================================
int is_in_amrh(int gfn)
{
   var *f=curr_context->vars;

   while(f)
   {
      if (f->in_amrh && (gfn >= f->sgfn) && (gfn < (f->sgfn+f->num_tl))) return 1;
      f=f->next;
   }
   return 0;
}

//=============================================================================
// returns 1 if given grid function can be allocated, as determined by status
//=============================================================================
int is_enabled(int gfn)
{
   var *f=curr_context->vars;
   int n;

   while(f)
   {
      if (f->in_amrh) n=f->num_tl; else n=0;
      if (f->in_mgh) n++;
      if ( (gfn >= f->sgfn) && (gfn < (f->sgfn+n)) ) 
      {
         if (f->status==VAR_STATUS_ON || f->status==VAR_STATUS_TURN_ON) return 1;
         else return 0;
      
      }
      f=f->next;
   }

   printf("gh.c, is_enabled(): bug ... should never get here\n");
   exit(1);
   return 0;
}

//=============================================================================
// returns sgfn if the grid function is the MG hierarchy, *AND* the variable
// has a grid function at some time within the AMR hierarhcy.
//
// Also, in num_tl, sets the number of AMR tl's of the corresponding var
//=============================================================================
int is_in_mg_and_amrh(int gfn, int *num_tl)
{
   var *f=curr_context->vars;

   while(f)
   {
      *num_tl=f->num_tl;
      if (f->in_mgh && f->in_amrh && gfn==(f->sgfn+f->num_tl)) return f->sgfn;
      f=f->next;
   }
   return 0;
}

//=============================================================================
// the following inserts a grid into a level structure, in unique order,
// using the bbox [x1,x2,y1,y2,z1,z2] for the comparison keys:
//    x1
//    y1
//    z1
//    x2
//    y2
//    z2
//=============================================================================
void insert_grid(level *lev, grid *g)
{
   grid *cg=lev->grids,*pg=0;
   int dim=g->dim;

   while(cg && (cg->bbox[0] < g->bbox[0])
            && ((dim < 2) || cg->bbox[2] < g->bbox[2])
            && ((dim < 3) || cg->bbox[4] < g->bbox[4])
            && (cg->bbox[1] < g->bbox[1])
            && ((dim < 2) || cg->bbox[3] < g->bbox[3])
            && ((dim < 3) || cg->bbox[5] < g->bbox[5])) 
   {
      pg=cg;
      cg=cg->next;
   }

   g->next=cg;
   g->prev=pg;
   if (cg) cg->prev=g;
   if (pg) pg->next=g; else lev->grids=g;
   
   return;
}

int sizeof_data(grid *g)
{
    int n;

    n=g->shape[0];
    if (g->dim>1) n*=g->shape[1];
    if (g->dim>2) n*=g->shape[2];

    return n;
}

//=============================================================================
// Saves a copy of the sgh-bboxes between levels min_lev and max_lev.
//=============================================================================
int save_sgh_bboxes(int min_lev, int max_lev, int num, int *lev, real *bbox)
{
   level *lv;
   grid *g;
   int i,l,j,n;
   context *c=curr_context;
   real *p;

   for (l=min_lev-1; l<max_lev; l++)
   {
      if (c->sgh_bboxes[l]) free(c->sgh_bboxes[l]); c->sgh_bboxes[l]=0;
      for (i=0,n=0; i<num; i++) if (lev[i]==(l+1)) n++;
      c->num_sgh_bboxes[l]=n;
      if (n)
      {
         if (!(p=c->sgh_bboxes[l]=(real *)malloc(sizeof(real)*n*2*c->dim)))
            { printf("save_bboxes: error ... out of mem\n"); return 0; }
         for (i=0; i<num; i++) 
         { 
            if (lev[i]==(l+1)) 
            {
               for (j=0; j<2*c->dim; j++) *p++=bbox[2*c->dim*i+j]; 
            }
         }
      }
   }

   return 1;
}

//=============================================================================
// builds a grid hierarchy from the given bbox list. no data is allocated,
// and the grids are not split over the nodes (i.e. build a 'sequential' grid
// hierarchy) 
//=============================================================================
cgh *build_sgh(int min_lev, int max_lev, int num, int *lev, real *bbox, real t)
{
   cgh *sgh;
   level *lv;
   grid *g;
   int i,l,j,n,d;
   context *c=curr_context;
   real *p;
   int gtrace_save;
   
   gtrace_save=gtrace;
   // gtrace=4;

   if (!(sgh=(struct cgh *)malloc(sizeof(struct cgh)))) 
      {printf("build_sgh: error ... out of mem\n"); return 0;}

   IFG(2) printf("   >> build_sgh(min_lev=%i, max_lev=%i, num=%i)\n",min_lev,max_lev,num);

   sgh->min_lev=min_lev;
   sgh->max_lev=max_lev;

   for (i=0; i<PAMR_MAX_LEVS; i++) sgh->levels[i]=0;

   for (i=0; i<num; i++)
   {
      if (lev[i]>max_lev || lev[i]<min_lev)
         {printf("build_sgh: error ... lev[i]=%i out of range\n",lev[i]); free_cgh(sgh,1,PAMR_MAX_LEVS); return 0;}
      l=lev[i]-1;
      IFG0(4) printf("adding grid %i at level %i\n",i,l+1);
      if (!(lv=sgh->levels[l]))
      {
         if (!(lv=sgh->levels[l]=(struct level *)malloc(sizeof(struct level)))) 
            {printf("build_sgh: error ... out of mem\n"); free_cgh(sgh,1,PAMR_MAX_LEVS); return 0;}
         lv->t=t;
         lv->dt=c->dt[l];
         for (j=0;j<c->dim;j++) lv->dx[j]=c->dx[l][j]; 
         lv->grids=0;
      }
      if (!(g=(struct grid *)malloc(sizeof(struct grid)))) 
         {printf("build_sgh: error ... out of mem\n"); free_cgh(sgh,1,PAMR_MAX_LEVS); return 0;}
      
      g->rank=-1;
      g->ngfs=c->gfns;
      g->comm=1;
      g->dim=c->dim;
      g->t=lv->t;
      g->gfs=0;
      g->coarsest=0;
      g->virtual_g=0;
      for (j=0; j<c->dim; j++)
      {
         g->bbox[2*j]=bbox[2*c->dim*i+2*j];
         g->bbox[2*j+1]=bbox[2*c->dim*i+2*j+1];
         g->shape[j]=(g->bbox[2*j+1]-g->bbox[2*j])/lv->dx[j]+1.5;
         g->x[j]=0;
         g->ghost_width[2*j]=0;
         g->ghost_width[2*j+1]=0;
         g->wrap[2*j]=g->wrap[2*j+1]=0;
      }
      IFG0(4) switch(g->dim)
      {
         case 1: printf("\t shape=%i ... bbox=[%16.10lf,%16.10lf]\n",g->shape[0],g->bbox[0],g->bbox[1]); break;
         case 2: printf("\t shape=%i x %i ... bbox=[%16.10lf,%16.10lf][%16.10lf,%16.10lf]\n",
                 g->shape[0],g->shape[1],g->bbox[0],g->bbox[1],g->bbox[2],g->bbox[3]); break;
         case 3: printf("\t shape=%i x %i x %i ... bbox=[%16.10lf,%16.10lf][%16.10lf,%16.10lf][%16.10lf,%16.10lf]\n",
                 g->shape[0],g->shape[1],g->shape[2],g->bbox[0],g->bbox[1],
                 g->bbox[2],g->bbox[3],g->bbox[4],g->bbox[5]); break;
      }
      insert_grid(lv,g);
   }

   //------------------------------------------------------------------
   // set the wrap flag if there are periodic boundaries
   //------------------------------------------------------------------
   for (l=min_lev; l<=max_lev; l++) if (lv=sgh->levels[l-1]) set_wrap(lv);

   IFG(3) debug_save_cgh("debug_build_sgh",sgh);

   gtrace=gtrace_save;

   IFG(2) printf("      build_sgh <<\n");
   return sgh;
}

void set_wrap(level *lv)
{
   int l,d,pm;
   context *c=curr_context;
   grid *g;
   int ltrace=0;

   IFL printf("set_wrap: periodic=[%i,%i,%i]\n",c->periodic[0],c->periodic[1],c->periodic[2]);

   for (d=1; d<=c->dim; d++)
   {
      for (pm=-1; pm<=1; pm+=2)
      {
         g=lv->grids;
         while(g)
         {
            IFL printf("  g->bbox[d=%i]=[%lf,%lf]\n",d,g->bbox[2*(d-1)],g->bbox[2*(d-1)+1]);
            if (c->periodic[d-1] && abuts(g,d,pm) && matched_abuts(lv->grids,d,-1*pm,g)) 
               g->wrap[2*(d-1)+(pm+1)/2]=1;
            else 
               g->wrap[2*(d-1)+(pm+1)/2]=0;
            IFL printf("  g->wrap[d=%i,pm=%i]=%i\n",d,pm,g->wrap[2*(d-1)+(pm+1)/2]);
            g=g->next;
         }
      }
   }
}
   
//=============================================================================
// the following two routines are used by set_wrap to determine wether
// a grid wraps around a periodic dimension
//
// abuts(g,d,pm) returns 1 if g abuts the x(d)_max (pm=1) or
//                                        x(d)_min (pm=-1) boundaries
//
// matched_abuts(gl,d,pm,gm) returns 1 if *any* of the grids in the list gl
//                                     abuts the x(d)_max (pm=1) or
//                                     x(d)_min (pm=-1) boundaries,
//                                     *and* that grid shares the same
//                                     width parameters along the non-abutting
//                                     dimensions as gm
//
// d is 1-indexed
//=============================================================================
int abuts(grid *g, int d, int pm)
{
   real dx;
   context *c=curr_context;

   if (g->shape[d-1]==0) {printf("abuts: error ... shape = 0 \n"); return 0;}
   dx=(g->bbox[2*(d-1)+1]-g->bbox[2*(d-1)])/(g->shape[d-1]);

   if (pm==-1)
   {
      if (fuzz_eq(c->bbox[2*(d-1)],g->bbox[2*(d-1)],dx/2)) return 1;
   }
   else
   {
      if (fuzz_eq(c->bbox[2*(d-1)+1],g->bbox[2*(d-1)+1],dx/2)) return 1;
   }

   return 0;
}

int matched_abuts(grid *gl, int d, int pm, grid *gm)
{
   real dx[PAMR_MAX_DIM];
   context *c=curr_context;
   grid *g;
   int d1,match;

   g=gl;

   for (d1=1; d1<=c->dim; d1++)
   {
      if (gm->shape[d1-1]==0) {printf("matched_abuts: error ... shape = 0 \n"); return 0;}
      dx[d1-1]=(gm->bbox[2*(d1-1)+1]-gm->bbox[2*(d1-1)])/(gm->shape[d1-1]);
   }

   while(g)
   {
      if (abuts(g,d,pm))
      {
         match=1;
         for (d1=1; d1<=c->dim; d1++)
         {
            if (d1!=d)
            {
               if (!(fuzz_eq(g->bbox[2*(d1-1)],gm->bbox[2*(d1-1)],dx[d1-1]/2) &&
                     fuzz_eq(g->bbox[2*(d1-1)+1],gm->bbox[2*(d1-1)+1],dx[d1-1]/2))) match=0;
            }
         }
         if (match) return 1;
      }
      g=g->next;
   }

   return 0;
}

//=============================================================================
// the following takes a sgh and eliminates all grid overlap by slicing
// away overlapping pieces of grids, and then adding appropriate 
// ghostzones to touching edges.
//
// I.e. Effectively shrinks overlap regions to the the size of the ghostzone
//
// NOTE: may need to modify the definition of ownership to exclude ghostzone
// width sized regions away from AMR boundaries, depending upon how we
// eventually implement the following (for correct behavior during 
// synchronization)
//
//=============================================================================
cgh *sgh_cull_overlap(cgh *sgh)
{
   static int first=1;

   if (first) printf("\nsgh_cull_overlap() not yet implemented\n"); 
   first=0;

   return sgh;
}

//=============================================================================
// frees all structures associate with a cgh between min_lev and max_lev 
// ... if that frees all of the cgh's levels, then the cgh structure is 
// freed too.
//=============================================================================
void free_cgh(cgh *c, int min_lev, int max_lev)
{
   level *lev;
   grid *q,*g;
   int l,i;

   IFG(2) printf("   >> free_cgh(min_lev=%i, max_lev=%i)\n",min_lev,max_lev);

   for (l=min_lev-1; l<=(max_lev-1); l++)
   {
      if (lev=c->levels[l])
      {
         g=lev->grids;
         while(g)
         {
            q=g;
            g=g->next;
            if (!q->virtual_g)
            {
               if (q->gfs)
               {
                  for (i=0; i<q->ngfs; i++) if (q->gfs[i]) pfree(q->gfs[i]);
                  free(q->gfs);
               }
               for (i=0; i<curr_context->dim; i++) if (q->x[i]) pfree(q->x[i]);
            }
            free(q);
         }
         free(lev);
         c->levels[l]=0;
      }
   }

   if (min_lev <= c->min_lev && max_lev >= c->max_lev) free(c);

   IFG(2) printf("      free_cgh <<\n");
   return;
}

//=============================================================================
// a few utility routines for compose_cgh() below
//=============================================================================
int block_size(grid *g)
{
   int size=0,i;

   size=g->shape[0];
   for (i=1; i<g->dim; i++)  size*=g->shape[i];
   return size;
}

int total_block_size(level *lev)
{
   grid *g;
   int size=0;

   g=lev->grids;
   while(g)
   {
      size+=block_size(g);
      g=g->next;
   }
   return size;
}

//=============================================================================
// This function takes an sgh (as constructed by build_sgh()) and
// creates the actually computational hiearchy, allocating all local memory
// for the grid functions; i.e. this function is responsible for distributing
// a hierarchy across the network.
//
// topology ratios ignored for now.
//
// if (alloc_mg), then *only* MG var's are allocated. 
//
// to deal with period boundaries, the compose_cgh() does the following,
// after initial parallel distribution but prior to memory allocation:
//   1)  all grids with the wrap flag set along a *max* boundary are extended
//       in that direction by 2*ghost_width-1 points (or zero if ghost_width=0).
//   2)  all boundaries that have the wrap attribute are given a ghostzone.
//=============================================================================
#define SEARCH_SIZE 5
cgh *compose_cgh(cgh *sgh, int alloc_mg)
{
   cgh *new_cgh;
   level *lev;
   grid *g,*ng;
   int i,j,k,l,dim=curr_context->dim,split_size,min_size,n,h,m;
   int nx,ny,nz,cnx,cny,cnz,mnx,mny,mnz,dn,maxnx,maxny,maxnz;
   int cmin_width,min_width;
   int sp1,sp2;
   int align,gw,rho;
   context *c=curr_context;
   real x,dx,fx,fy,fz,max_fxfy,max_fyfz,max_fxfz;
   int gtrace_save=gtrace;

   // gtrace=4;

   IFG(2) printf("   >> compose_cgh(min_lev=%i, max_lev=%i)\n",sgh->min_lev,sgh->max_lev);

   if (!(new_cgh=(struct cgh *)malloc(sizeof(struct cgh)))) 
      {printf("compose_cgh: error ... out of mem\n"); return 0;}

   if (c->gdm & PAMR_GDM_ALIGN)
      align=1;
   else
      align=0;

   new_cgh->min_lev=sgh->min_lev;
   new_cgh->max_lev=sgh->max_lev;

   for (i=0; i<PAMR_MAX_LEVS; i++) new_cgh->levels[i]=0;

   min_size=c->min_width[0];
   for (i=1; i<dim; i++)  min_size*=c->min_width[i];

   for (l=sgh->min_lev-1; l<=(sgh->max_lev-1); l++)
   {
      if (sgh->levels[l])
      {
         IFG0(4) printf("level %i:\n\n",l+1);
         if (!(lev=new_cgh->levels[l]=(struct level *)malloc(sizeof(struct level)))) 
            {printf("compose_cgh: error ... out of mem\n"); goto error;}
         lev->t=sgh->levels[l]->t;
         lev->dt=sgh->levels[l]->dt;
         for (j=0;j<c->dim;j++) lev->dx[j]=sgh->levels[l]->dx[j]; 
         lev->grids=0;
     
         //-----------------------------------------------------------------------
         // size of blocks to split grids into, if GDM = level by level
         //-----------------------------------------------------------------------
         split_size=max(min_size,total_block_size(sgh->levels[l])/c->size); 

         g=sgh->levels[l]->grids;
         while(g)
         {
            if (c->gdm & PAMR_GDM_GRID_BY_GRID) split_size=max(min_size,block_size(g)/c->size); 
            split_size=max(1,split_size);
            //--------------------------------------------------------------------
            // we want to split g into n blocks, if possible. For 2d and 3d,
            // choose a splitting nx by ny [by nz] which is 'close' to n,
            // and try to make the blocks as square/cubic as possible.
            // We also round up to the nearest n
            //--------------------------------------------------------------------
            n=max(1,(block_size(g)+split_size/2)/split_size);
            IFG0(4) printf("block_size=%i\t split_size=%i\tdesired n=%i\n",block_size(g),split_size,n);
            switch(dim)
            {
               case 1: 
                  maxnx=max(1,g->shape[0]/c->min_width[0]); maxnx=min(c->size,maxnx);
                  nx=max(1,min(n,maxnx)); ny=nz=1; 
                  IFG0(4) printf("\tsplitting [%i] into %i segments\n",g->shape[0],n);
                  IFG0(4) printf("source bbox=[%16.10lf,%16.10lf]\n",g->bbox[0],g->bbox[1]);
                  break;
               case 2: 
                  maxnx=max(1,g->shape[0]/c->min_width[0]); maxnx=min(c->size,maxnx);
                  maxny=max(1,g->shape[1]/c->min_width[1]); maxny=min(c->size,maxny);
                  fx=(double)g->shape[0]/(g->shape[0]+g->shape[1]);
                  fy=(double)g->shape[1]/(g->shape[0]+g->shape[1]);
                  nx=mnx=max(1,min(maxnx,sqrt(n)*fx/fy));
                  ny=mny=max(1,min(maxny,sqrt(n)*fy/fx));
                  nz=1;
                  dn=abs(n-nx*ny);
                  min_width=min(g->shape[0]/nx,g->shape[1]/ny);
                  for (cnx=max(1,mnx-SEARCH_SIZE); cnx<=(min(mnx+SEARCH_SIZE,maxnx)); cnx++)
                     for (cny=max(1,mny-SEARCH_SIZE); cny<=(min(mny+SEARCH_SIZE,maxny)); cny++)
                     {
                        cmin_width=min(g->shape[0]/cnx,g->shape[1]/cny);
                        if (dn>abs(n-cnx*cny) || (dn==abs(n-cnx*cny) && cmin_width>min_width))
                        {
                           dn=abs(n-cnx*cny);
                           nx=cnx;
                           ny=cny;
                           min_width=cmin_width;
                        }
                     }
                  n=nx*ny;
                  IFG0(4) printf("\tsplitting [%i x %i] into %i x %i blocks\n",
                              g->shape[0],g->shape[1],nx,ny);
                  IFG0(4) printf("source bbox=[%16.10lf,%16.10lf][%16.10lf,%16.10lf]\n",
                              g->bbox[0],g->bbox[1],g->bbox[2],g->bbox[3]);
                  break;
               case 3:
                  maxnx=max(1,g->shape[0]/c->min_width[0]); maxnx=min(c->size,maxnx);
                  maxny=max(1,g->shape[1]/c->min_width[1]); maxny=min(c->size,maxny);
                  maxnz=max(1,g->shape[2]/c->min_width[2]); maxnz=min(c->size,maxnz);
                  fx=(double)g->shape[0]/(g->shape[0]+g->shape[1]+g->shape[2]);
                  fy=(double)g->shape[1]/(g->shape[0]+g->shape[1]+g->shape[2]);
                  fz=(double)g->shape[2]/(g->shape[0]+g->shape[1]+g->shape[2]);
                  max_fxfy=max(fx,fy);
                  max_fxfz=max(fx,fz);
                  max_fyfz=max(fy,fz);
                  nx=mnx=max(1,min(maxnx,pow(n,1.0/3.0)*fx/max_fyfz));
                  ny=mny=max(1,min(maxny,pow(n,1.0/3.0)*fy/max_fxfz));
                  nz=mnz=max(1,min(maxnz,pow(n,1.0/3.0)*fz/max_fxfy));
                  IFG0(4) printf("nx,nx,nz=%i,%i,%i   fx,fy,fz=%16.10lf,%16.10lf,%16.10lf\n",nx,nx,nz,fx,fy,fz);
                  IFG0(4) printf("max_fxfy,max_fxfz,max_fyfz=%16.10lf,%16.10lf,%16.10lf\n",max_fxfy,max_fxfz,max_fyfz);
                  IFG0(4) printf("n, pow(n) = %i,%16.10lf\n",n,pow(n,1.0/3.0));
                  dn=abs(n-nx*ny*nz);
                  min_width=min(g->shape[0]/nx,g->shape[1]/ny);
                  min_width=min(min_width,g->shape[2]/nz);
                  for (cnx=max(1,mnx-SEARCH_SIZE); cnx<=(min(mnx+SEARCH_SIZE,maxnx)); cnx++)
                     for (cny=max(1,mny-SEARCH_SIZE); cny<=(min(mny+SEARCH_SIZE,maxny)); cny++)
                        for (cnz=max(1,mnz-SEARCH_SIZE); cnz<=(min(mnz+SEARCH_SIZE,maxnz)); cnz++)
                        {
                           cmin_width=min(g->shape[0]/cnx,g->shape[1]/cny);
                           cmin_width=min(cmin_width,g->shape[2]/cnz);
                           if (dn>abs(n-cnx*cny*cnz) || (dn==abs(n-cnx*cny*cnz) && cmin_width>min_width))
                           {
                              dn=abs(n-cnx*cny*cnz);
                              nx=cnx;
                              ny=cny;
                              nz=cnz;
                              min_width=cmin_width;
                           }
                        }
                  n=nx*ny*nz;
                  IFG0(4) printf("\tsplitting [%i x %i x %i] into %i x %i x %i blocks\n",
                              g->shape[0],g->shape[1],g->shape[2],nx,ny,nz);
                  IFG0(4) printf("source bbox=[%16.10lf,%16.10lf][%16.10lf,%16.10lf][%16.10lf,%16.10lf]\n",
                              g->bbox[0],g->bbox[1],g->bbox[2],g->bbox[3],g->bbox[4],g->bbox[5]);
                  break;
               default:
                  printf("compose_cgh: error ... currently only dim=1,2 or 3 supported\n");
                  goto error;
            }
            for (i=0; i<nx; i++)
            {
               for (j=0; j<ny; j++)
               {
                  for (k=0; k<nz; k++)
                  {
                     if (!(ng=(struct grid *)malloc(sizeof(struct grid)))) 
                     {
                        printf("compose_cgh: error ... out of mem\n"); 
                     }
                     ng->rank=c->n_rank++;  
                     if (c->n_rank==c->size) c->n_rank=0;
                     ng->ngfs=g->ngfs;
                     ng->comm=1;
                     ng->dim=dim;
                     ng->t=g->t;
                     ng->coarsest=g->coarsest;
                     ng->virtual_g=0;
                     if (i==0 && g->wrap[0]) ng->wrap[0]=1;
                     else ng->wrap[0]=0;

                     sp1=(g->shape[0]*i)/nx;
                     sp2=(g->shape[0]*(i+1))/nx;
                     rho=c->rho_sp[l-1];
                     gw=c->ghost_width[0];
                     if (align && l>0)
                     {
                        sp1-=(sp1%rho);
                        sp2-=(sp1%rho);
                        if (gw%rho) gw+=(rho-gw%rho);
                     }
                     sp1=max(0,sp1-gw);
                     sp2=min(g->shape[0]-1,sp2+gw);
                     // -1 here, so that a *unique* node updates each point (as
                     // much as possible):
                     if (!(align && l>0) && sp2!=(g->shape[0]-1)) sp2--;
                     if (i==(nx-1) && g->wrap[1]) 
                     {
                        ng->wrap[1]=1;
                        if (c->ghost_width[0]>0 && align) sp2+=2*c->ghost_width[0];
                        else if (c->ghost_width[0]>0) sp2+=2*c->ghost_width[0]-1;
                     }
                     else ng->wrap[1]=0;
                     ng->shape[0]=sp2-sp1+1;
                     ng->bbox[0]=g->bbox[0]+lev->dx[0]*sp1;
                     ng->bbox[1]=g->bbox[0]+lev->dx[0]*sp2;
                     if (i!=0) ng->ghost_width[0]=gw; 
                     else if (g->wrap[0]) ng->ghost_width[0]=c->ghost_width[0];
                     else ng->ghost_width[0]=0;
                     if (i!=(nx-1)) { ng->ghost_width[1]=gw; if (align && l>0) ng->ghost_width[1]++; } 
                     else if (ng->wrap[1]) ng->ghost_width[1]=c->ghost_width[0];
                     else ng->ghost_width[1]=0;
                     if (dim>1)
                     {
                        if (j==0 && g->wrap[2]) ng->wrap[2]=1;
                        else ng->wrap[2]=0;

                        sp1=(g->shape[1]*j)/ny;
                        sp2=(g->shape[1]*(j+1))/ny;
                        gw=c->ghost_width[1];
                        if (align && l>0)
                        {
                           sp1-=(sp1%rho);
                           sp2-=(sp1%rho);
                           if (gw%rho) gw+=(rho-gw%rho);
                        }
                        sp1=max(0,sp1-gw);
                        sp2=min(g->shape[1]-1,sp2+gw);
                        if (!(align && l>0) && sp2!=(g->shape[1]-1)) sp2--;
                        if (j==(ny-1) && g->wrap[3]) 
                        {
                           ng->wrap[3]=1;
                           if (c->ghost_width[1]>0 && align) sp2+=2*c->ghost_width[1];
                           else if (c->ghost_width[1]>0) sp2+=2*c->ghost_width[1]-1;
                        }
                        else ng->wrap[3]=0;
                        ng->shape[1]=sp2-sp1+1;
                        ng->bbox[2]=g->bbox[2]+lev->dx[1]*sp1;
                        ng->bbox[3]=g->bbox[2]+lev->dx[1]*sp2;
                        if (j!=0) ng->ghost_width[2]=gw; 
                        else if (g->wrap[2]) ng->ghost_width[2]=c->ghost_width[1];
                        else ng->ghost_width[2]=0;
                        if (j!=(ny-1)) { ng->ghost_width[3]=gw; if (align && l>0) ng->ghost_width[3]++; }
                        else if (ng->wrap[3]) ng->ghost_width[3]=c->ghost_width[1];
                        else ng->ghost_width[3]=0;
                     }
                     if (dim>2)
                     {
                        if (k==0 && g->wrap[4]) ng->wrap[4]=1;
                        else ng->wrap[4]=0;

                        sp1=(g->shape[2]*k)/nz;
                        sp2=(g->shape[2]*(k+1))/nz;
                        gw=c->ghost_width[2];
                        if (align && l>0)
                        {
                           sp1-=(sp1%rho);
                           sp2-=(sp1%rho);
                           if (gw%rho) gw+=(rho-gw%rho);
                        }
                        sp1=max(0,sp1-gw);
                        sp2=min(g->shape[2]-1,sp2+gw);
                        if (!(align && l>0) && sp2!=(g->shape[2]-1)) sp2--;
                        if (k==(nz-1) && g->wrap[5]) 
                        {
                           ng->wrap[5]=1;
                           if (c->ghost_width[2]>0 && align) sp2+=2*c->ghost_width[2];
                           else if (c->ghost_width[2]>0) sp2+=2*c->ghost_width[2]-1;
                        }
                        else ng->wrap[5]=0;
                        ng->shape[2]=sp2-sp1+1;
                        ng->bbox[4]=g->bbox[4]+lev->dx[2]*sp1;
                        ng->bbox[5]=g->bbox[4]+lev->dx[2]*sp2;
                        if (k!=0) ng->ghost_width[4]=gw; 
                        else if (g->wrap[4]) ng->ghost_width[4]=c->ghost_width[2];
                        else ng->ghost_width[4]=0;
                        if (k!=(nz-1)) { ng->ghost_width[5]=gw; if (align && l>0) ng->ghost_width[5]++; }
                        else if (ng->wrap[5]) ng->ghost_width[5]=c->ghost_width[2];
                        else ng->ghost_width[5]=0;
                     }
                     IFG0(4) switch(dim)
                     {
                        case 1: printf("segment %i bbox=[%16.10lf,%16.10lf], rank=%i\n",
                                       i+j*nx+k*nx*ny,ng->bbox[0],ng->bbox[1],ng->rank);
                                break;
                        case 2: printf("segment %i bbox=[%16.10lf,%16.10lf][%16.10lf,%16.10lf], rank=%i\n",i+j*nx+k*nx*ny,
                                ng->bbox[0],ng->bbox[1],ng->bbox[2],ng->bbox[3],ng->rank);
                                break;
                        case 3: printf("segment %i bbox=[%16.10lf,%16.10lf][%16.10lf,%16.10lf][%16.10lf,%16.10lf], rank=%i\n",i+j*nx+k*nx*ny,
                                ng->bbox[0],ng->bbox[1],ng->bbox[2],ng->bbox[3],ng->bbox[4],ng->bbox[5],ng->rank);
                                break;
                     }
                     for (h=0;h<c->dim; h++) ng->x[h]=0;
                     ng->gfs=0;
                     //-----------------------------------------------------------
                     // if on the local node, allocate grid memory
                     // and coordinate arrays ... also fill the coordinate arrays
                     //-----------------------------------------------------------
                     if (ng->rank==c->rank)
                     {
                        if(!(ng->gfs=(real **)malloc(sizeof(real *)*ng->ngfs)))
                           {printf("compose_cgh: error ... out of mem\n"); goto error;}
                        for (h=0;h<ng->ngfs;h++) ng->gfs[h]=0;
                        for (h=0;h<ng->ngfs;h++)
                        {
                           if (is_enabled(h+1) && 
                               ((!alloc_mg && is_in_amrh(h+1)) || (alloc_mg && !(is_in_amrh(h+1)))) )
                           {
                              // IFG0(4) printf("allocating memory for gfn %i\n",h+1);
                              if (!(ng->gfs[h]=(real *)pmalloc(sizeof(real)*block_size(ng))))
                                 {printf("compose_cgh: error ... out of mem\n"); goto error;}
                              // zero memory, for convenience
                              for (m=0; m<block_size(ng); m++) (ng->gfs[h])[m]=0;
                           }
                        }
                        for (h=0;h<c->dim; h++)
                        {
                           if (!(ng->x[h]=(real *)pmalloc(sizeof(real)*ng->shape[h])))
                              {printf("compose_cgh: error ... out of mem\n"); goto error;}
                           dx=(ng->bbox[2*h+1]-ng->bbox[2*h])/(ng->shape[h]-1);
                           for (m=0,x=ng->bbox[2*h]; m<ng->shape[h]; m++,x+=dx) (ng->x[h])[m]=x;
                        }
                     }
                     insert_grid(lev,ng);
                  }
               }
            }
            g=g->next;
         }
      }
   }

   IFG(3) debug_save_cgh("debug_compose_cgh",new_cgh);

   IFG(2) printf("      compose_cgh <<\n");
 
   gtrace=gtrace_save;
   
   return new_cgh;
error:
   free_cgh(new_cgh,1,PAMR_MAX_LEVS); 
   return 0;
}

//=============================================================================
// the following function initializes a new cgh, after regridding, by
// interpolating (up to 1-level)/transferring data from an old cgh.
//=============================================================================
void init_new_cgh(cgh *new_cgh,cgh *old_cgh)
{
   int l,ret,i;
   gsl *src_int[PAMR_MAX_NODES],*dst,*src_trans[PAMR_MAX_NODES],*gs;
   gsl *transfer_src[PAMR_MAX_NODES],*transfer_dst[PAMR_MAX_NODES];
   char db_name[256];
   context *c=curr_context;
   real t;

   IFG(2) printf("  >> init_new_cgh\n");

   // support for periodic boundaries
   for (l=max(0,new_cgh->min_lev-2); l<=(new_cgh->max_lev-1); l++)
      if (old_cgh->levels[l]) alloc_virtual_grids(old_cgh->levels[l]);

   for (l=new_cgh->min_lev-1; l<=(new_cgh->max_lev-1); l++)
   {
      if (new_cgh->levels[l])
      {
         t=new_cgh->levels[l]->t;
         dst=0;
         for (i=0; i<c->size; i++) src_int[i]=src_trans[i]=transfer_src[i]=transfer_dst[i]=0;
         PAMR_set_tf_bits(-1,PAMR_AMRH,PAMR_TF_COMPOSE);
         for (i=0; i<c->size; i++)
         {
            if (l>=1 && old_cgh->levels[l-1]) src_int[i]=build_owned_gsl(i,old_cgh,l,0,1);
            if (old_cgh->levels[l]) src_trans[i]=build_owned_gsl(i,old_cgh,l+1,0,0);
         }
         IFG(3) 
         {
            if (src_int[c->rank]) 
            {
               sprintf(db_name,"regrid_src_int_%i",c->rank); 
               debug_save_gsl(db_name,t,src_int[c->rank],-1);
            }
            if (src_trans[c->rank]) 
            {
               sprintf(db_name,"regrid_src_trans_%i",c->rank); 
               debug_save_gsl(db_name,t,src_trans[c->rank],-1);
            }
         }

         //-----------------------------------------------------------------------
         // NOTE: at this point, we could save a bit in communcation bandwidth
         //       by only interpolating (src_int '-' src_trans). However, to avoid
         //       missing the interpolation of strips of size (rho_sp[l]-1) on the
         //       fine grid, we'd need to (temporarily) 'shrink' src_trans by 
         //       rho_sp[l] points before subtracting the gsl's. Ignore this
         //       for now, and just overwrite the overlapping parts with
         //       transfered data. (NOTE: may be less of a problem now with
         //       extend flag added to build_owned_gsl ... investigate if we
         //       need to speed up this part).
         //-----------------------------------------------------------------------
   
         dst=build_complete_gsl(new_cgh,l+1,0);
         IFG0(3) {sprintf(db_name,"regrid_dst"); debug_save_gsl(db_name,t,dst,-1);}
   
         for (i=0; i<c->size; i++) if (src_int[i] && dst) build_gstl(src_int[i],dst,&transfer_src[i],&transfer_dst[i]);
   
         IFG(3) 
         {
            if (transfer_src[c->rank] && transfer_dst[c->rank]) 
            {
               sprintf(db_name,"regrid_interp_src_gsl_%i",c->rank); debug_save_gsl(db_name,t,transfer_src[c->rank],-1);
               sprintf(db_name,"regrid_interp_dst_gsl_%i",c->rank); debug_save_gsl(db_name,t,transfer_dst[c->rank],-1);
            }
         }

         transfer(transfer_src,transfer_dst);

         for (i=0; i<c->size; i++)
         {
            if (transfer_src[i]) free_gsl(transfer_src[i]);
            if (transfer_dst[i]) free_gsl(transfer_dst[i]);
            if (src_int[i]) free_gsl(src_int[i]);
            transfer_src[i]=transfer_dst[i]=0;
         }

         for (i=0; i<c->size; i++) if (src_trans[i] && dst) build_gstl(src_trans[i],dst,&transfer_src[i],&transfer_dst[i]);
         IFG(3) 
         {
            if (transfer_src[c->rank] && transfer_dst[c->rank]) 
            {
               sprintf(db_name,"regrid_trans_src_gsl_%i",c->rank); debug_save_gsl(db_name,t,transfer_src[c->rank],-1);
               sprintf(db_name,"regrid_trans_dst_gsl_%i",c->rank); debug_save_gsl(db_name,t,transfer_dst[c->rank],-1);
            }
         }
        
         transfer(transfer_src,transfer_dst);

         for (i=0; i<c->size; i++)
         {
            if (transfer_src[i]) free_gsl(transfer_src[i]);
            if (transfer_dst[i]) free_gsl(transfer_dst[i]);
            if (src_trans[i]) free_gsl(src_trans[i]);
         }
         if (dst) free_gsl(dst);
      }
   }

   for (l=max(0,new_cgh->min_lev-2); l<=(new_cgh->max_lev-1); l++)
      if (old_cgh->levels[l]) free_virtual_grids(old_cgh->levels[l]);

   IFG(2) printf("     init_new_cgh << \n");
}

//=============================================================================
// replace the set of levels in curr_cgh with those in new_cgh,
// AND deletes the replaced levels and new_cgh structure
//
// returns 0 if afterwards curr_cgh is empty (curr_cgh is free'd then!)
//=============================================================================
int incorporate_new_cgh(cgh *new_cgh, cgh *curr_cgh)
{
   int l,new_min_lev,new_max_lev;
   level *lev;

   IFG(2) printf("  >> incorporate_new_cgh\n");

   IFG(2)
   {
      printf("curr min/max lev:%i,%i . new min/max lev: %i,%i\n",curr_cgh->min_lev,curr_cgh->max_lev,new_cgh->min_lev,new_cgh->max_lev);
      for (l=0; l<PAMR_MAX_LEVS; l++) printf("l=%i. new[l]=%i, old[l]=%i\n",l,new_cgh->levels[l],curr_cgh->levels[l]);
   }

   for (l=new_cgh->min_lev-1; l<=(new_cgh->max_lev-1); l++)
   {
      lev=curr_cgh->levels[l];
      curr_cgh->levels[l]=new_cgh->levels[l];
      new_cgh->levels[l]=lev;
   }

   free_cgh(new_cgh,1,PAMR_MAX_LEVS);

   curr_cgh->min_lev=PAMR_MAX_LEVS+1;
   curr_cgh->max_lev=0;
   for (l=0; l<PAMR_MAX_LEVS; l++)
   {
      if (curr_cgh->levels[l])
      {
         curr_cgh->max_lev=l+1;
         curr_cgh->min_lev=min(curr_cgh->min_lev,l+1);
      }
   }
   if (curr_cgh->max_lev==0)
   {
      printf("PAMR: incorporate_new_cgh WARNING ...cgh deleted\n");
      free_cgh(curr_cgh,1,PAMR_MAX_LEVS);
      return 0;
   }

   for (l=curr_cgh->min_lev; l<(curr_cgh->max_lev-1); l++)
      if (!curr_cgh->levels[l]) printf("PAMR: incorporate_new_cgh WARNING ... new cgh deletes an intermediate level\n");

   IFG(2) printf("     incorporate_new_cgh <<\n");

   return 1;
}

//=============================================================================
// The following function allocates a MG level, with the same structure,
// and sharing like functions, as the AMR level
//=============================================================================
int amr_to_mgh(level *amr_lev, level **mg_lev, int tl)
{
   level *nlev;
   int i,h,m,sgfn,ctl;
   grid *g,*pg,*ng;
   context *c=curr_context;

   IFG(2) printf("  >> amr_to_mgh \n");

   if (!(*mg_lev=nlev=(level *)malloc(sizeof(level)))) { printf("amr_to_mgh: error ... out of memory\n"); return 0; }

   nlev->t=amr_lev->t;
   for (i=0; i<PAMR_MAX_DIM; i++) nlev->dx[i]=amr_lev->dx[i];
   nlev->dt=amr_lev->dt;

   g=amr_lev->grids;
   nlev->grids=0; pg=0;
   while(g)
   {
      if (!(ng=(grid *)malloc(sizeof(grid)))) { printf("amr_to_mgh: error ... out of memory\n"); return 0; }
      ng->prev=pg;
      ng->next=0;
      if (pg) pg->next=ng; else nlev->grids=ng;
      ng->dim=g->dim;
      for (i=0; i<PAMR_MAX_DIM; i++) 
      {
         ng->shape[i]=g->shape[i];
         ng->bbox[2*i]=g->bbox[2*i];
         ng->bbox[2*i+1]=g->bbox[2*i+1];
         ng->ghost_width[2*i]=g->ghost_width[2*i];
         ng->ghost_width[2*i+1]=g->ghost_width[2*i+1];
         ng->x[i]=g->x[i];
         ng->wrap[2*i]=g->wrap[2*i];
         ng->wrap[2*i+1]=g->wrap[2*i+1];
      }
      ng->t=g->t;
      ng->ngfs=g->ngfs;
      ng->comm=g->comm;
      ng->rank=g->rank;
      ng->gfs=g->gfs;
      ng->coarsest=0;
      ng->virtual_g=0;
      if (g->rank==curr_context->rank)
      {
         for (h=0; h<ng->ngfs; h++)
         {
            if (sgfn=is_in_mg_and_amrh(h+1,&ctl))
            {
               ctl=min(tl,ctl);
               ng->gfs[h]=g->gfs[sgfn-1+ctl-1];
            }
            else if (!(is_in_amrh(h+1)) && is_enabled(h+1))
            {
               if (ng->gfs[h]) { printf("amr_to_mgh: internal error ... ng->gfs[h]!=0\n"); return 0; }
   
               if (!(ng->gfs[h]=(real *)pmalloc(sizeof(real)*block_size(ng))))
                  {printf("amr_to_mgh: error ... out of mem\n"); return 0;}
               for (m=0; m<block_size(ng); m++) (ng->gfs[h])[m]=0;
            }
         }
      }
      pg=ng;
      g=g->next;
   }

   IFG(2) printf("  << amr_to_mgh\n");

   return 1;
}

//=============================================================================
// The following fills the given mg level with 'empty' grids ... i.e.
// the equivalent of build_sgh
//=============================================================================
int fill_mgh_level(level *lev, int num, real *bbox, real *island_no)
{
   grid *g;
   int i,l,j,n;
   context *c=curr_context;
   real *p;

   lev->grids=0;

   for (i=0; i<num; i++)
   {
      if (!(g=(struct grid *)malloc(sizeof(struct grid)))) 
         {printf("fill_mgh_level: error ... out of mem\n"); return 0;}

      g->rank=-1;
      g->ngfs=c->gfns;
      g->dim=c->dim;
      g->t=lev->t;
      g->gfs=0;
      g->comm=1;
      g->virtual_g=0;
      if (island_no && island_no[i]<0) g->coarsest=1; else g->coarsest=0;
      for (j=0; j<c->dim; j++)
      {
         g->bbox[2*j]=bbox[2*c->dim*i+2*j];
         g->bbox[2*j+1]=bbox[2*c->dim*i+2*j+1];
         g->shape[j]=(g->bbox[2*j+1]-g->bbox[2*j])/lev->dx[j]+1.5;
         g->x[j]=0;
         g->ghost_width[2*j]=0;
         g->ghost_width[2*j+1]=0;
         g->wrap[2*j]=0;
         g->wrap[2*j+1]=0;
      }
      insert_grid(lev,g);
   }

   set_wrap(lev);

   return 1;
}

//=============================================================================
// The following builds the mgh between AMR levels min_lev and max_lev,
// using AMR data from time level tl
//
// NOTES: 1) here we assume that the hierachy is B&O style --- i.e. finer
//        grids are entirely contained in coarser ones.
//
//        2) we do *NOT* check for alignment problems ... that is up to the
//           user when composing the AMR hierarchy
//=============================================================================
int build_mgh(int min_lev, int max_lev, int tl)
{
   int l,cl,i,n,num,j,ok,cf;
   context *c=curr_context;
   cgh *mgh,*sgh;
   cgh *new_cgh;
   real *mg_bboxes,*p,*island_no,*sgh_island_no;
   int shape[PAMR_MAX_DIM],mem=1;
   grid *g;

   IFG(2) printf("   >> build_mgh(min_lev=%i,max_lev=%i,tl=%i)\n",min_lev,max_lev,tl);

   if (!(mgh=c->MG_cgh=(cgh *)malloc(sizeof(cgh)))) { printf("build_mgh: out of memory ... \n"); return 0; }
   for (l=0; l<PAMR_MAX_LEVS; l++) { mgh->levels[l]=0; mgh->in_amrh[l]=0; }

   PAMR_thaw_tf_bits();
   PAMR_set_tf_bits(1,PAMR_MGH,PAMR_TF_MGH_INIT);
   PAMR_freeze_tf_bits(); // we don't want PAMR_inject below to reset them

   //--------------------------------------------------------------------------
   // build the mgh from the finest level down, labeling it PAMR_MAX_LEV;
   // later shift so that level 1 is the first level
   //--------------------------------------------------------------------------
   cl=PAMR_MAX_LEVS;
   mgh->max_lev=cl;

   for (l=max_lev-1; l>=(min_lev-1); l--)
   {
      cl--;
      //-----------------------------------------------------------------------
      // replicate this AMR level
      //-----------------------------------------------------------------------
      if (!(amr_to_mgh(c->curr_cgh->levels[l],&(mgh->levels[cl]),tl))) goto error;
      mgh->in_amrh[cl]=1;
      mgh->min_lev=cl+1;

      IFG(2) printf("      replicated AMR level %i->%i\n",l+1,cl+1);

      //-----------------------------------------------------------------------
      // Now add as many MG levels as needed to bridge the gap to the next 
      // AMR level, or as far beyond min_lev as we can go. 
      // It will be an error for a grid to ''vanish'' (i.e. having become 
      // coarsest at a prior level) in-between AMR levels, so we don't check 
      // for minimum sizes then.
      //-----------------------------------------------------------------------
      if (l>(min_lev-1))
      {
         n=log(c->rho_sp[l-1])/log(2)+0.5;
         if (((int)(pow(2,n)+0.5))!=c->rho_sp[l-1])
            { printf("build_mgh: error ... rho_sp(=%i) must be a power of 2\n",c->rho_sp[l-1]); mem=0; goto error; }
         n--;
      }
      else n=-1;
      cf=1;
      while(n)
      {
         cl--; cf*=2;
         if (!(sgh=(cgh *) malloc(sizeof(cgh)))) goto error;
         for (i=0; i<PAMR_MAX_LEVS; i++) sgh->levels[i]=0;
         sgh->max_lev=cl+1;
         sgh->min_lev=cl+1; 
         if (n<0)
         {
            IFG(2) printf("      n<0 ... extending; num_sgh_bboxes[%i]=%i\n",l+1,c->num_sgh_bboxes[l]);
            num=0;
            if (!(p=mg_bboxes=(real *)malloc(sizeof(real)*(2*PAMR_MAX_DIM+2)*c->num_sgh_bboxes[l])))
               { printf("build_mgh: out of memory ... \n"); return 0; }
            island_no=&p[2*PAMR_MAX_DIM*c->num_sgh_bboxes[l]];
            sgh_island_no=&p[(2*PAMR_MAX_DIM+1)*c->num_sgh_bboxes[l]];
            calc_island_no(c->sgh_bboxes[l],c->num_sgh_bboxes[l],island_no,c->curr_cgh->levels[l]->dx,cf);
            for (i=0; i<c->num_sgh_bboxes[l]; i++)
            {
               if (island_no[i]!=0)
               {
                  sgh_island_no[num]=island_no[i];
                  num++;
                  for (j=0; j<2*c->dim; j++) *p++=(c->sgh_bboxes[l])[2*i*c->dim+j];
               }
            }
         }
         else
         {
            num=c->num_sgh_bboxes[l];
            mg_bboxes=c->sgh_bboxes[l];
            sgh_island_no=0;
         }
         if (num)
         {
            if (!(sgh->levels[cl]=(struct level *)malloc(sizeof(struct level)))) 
               {free_cgh(sgh,1,PAMR_MAX_LEVS); goto error;}
            sgh->levels[cl]->t=c->curr_cgh->levels[l]->t;
            sgh->levels[cl]->dt=c->curr_cgh->levels[l]->dt; 
            for (i=0; i<PAMR_MAX_DIM; i++) sgh->levels[cl]->dx[i]=c->curr_cgh->levels[l]->dx[i]*cf;

            if (!(fill_mgh_level(sgh->levels[cl],num,mg_bboxes,sgh_island_no))) { mem=0; goto error; }
            n--;
            //--------------------------------------------------------------------
            // compose single new level
            //--------------------------------------------------------------------
            new_cgh=compose_cgh(sgh,1);
            if (c->excision_on) init_ex_mask(sgh->min_lev,sgh->max_lev,new_cgh,c->mg_mask_gfn);
            if (!(incorporate_new_cgh(new_cgh,mgh)))
            {
               printf("build_mgh: error ... incorporate_new_cgh returned 0\n"); 
               mem=0;
               c->MG_cgh=0;
               goto error;
            }
            //--------------------------------------------------------------------
            // initialize new data
            //--------------------------------------------------------------------
            mgh->min_lev=cl+1; 
            PAMR_inject(cl+2,1,PAMR_MGH);
            IFG(2) printf("      created new MG level, l=%i\n",cl+1);
         }
         if (n<0) free(mg_bboxes);
         if (!num) n=0;
         free_cgh(sgh,1,PAMR_MAX_LEVS);
      }
   }

   for (i=0; i<=(mgh->max_lev-mgh->min_lev); i++) 
   {
      mgh->levels[i]=mgh->levels[i+mgh->min_lev-1];
      mgh->in_amrh[i]=mgh->in_amrh[i+mgh->min_lev-1];
   }
   for (i=mgh->max_lev-mgh->min_lev+1; i<PAMR_MAX_LEVS; i++) { mgh->levels[i]=0; mgh->in_amrh[i]=0; }

   mgh->max_lev=mgh->max_lev-mgh->min_lev+1;
   mgh->min_lev=1;
   
   g=mgh->levels[0]->grids;
   while(g) { g->coarsest=1; g=g->next; }

   PAMR_thaw_tf_bits();

   IFG(2) printf("   << build_mgh\n");
   return 1;

error:
   if (mem) printf("build_mgh: error ... out of mem\n");
   destroy_mgh();
   return 0;
}

//=============================================================================
// frees all structures associate with a MGH, unlinks pointers to the AMRH,
// and frees the MGH
//=============================================================================
void destroy_mgh(void)
{
   level *lev;
   context *c=curr_context;
   cgh *mgh=curr_context->MG_cgh;
   grid *q,*g;
   int l,i,ctl;

   if (!mgh) return;
   IFG(4) printf("   >> destroy_mgh()\n");
   
   for (l=mgh->min_lev-1; l<=(mgh->max_lev-1); l++)
   {
      if (lev=mgh->levels[l])
      {
         g=lev->grids;
         while(g)
         {
            q=g;
            g=g->next;
            if (q->gfs)
            {
               for (i=0; i<q->ngfs; i++) 
               {  
                  if (q->gfs[i] && (!mgh->in_amrh[l] || (!is_in_amrh(i+1) && !is_in_mg_and_amrh(i+1,&ctl))))
                  {
                     pfree(q->gfs[i]);
                     q->gfs[i]=0;
                  }
                  else if (is_in_mg_and_amrh(i+1,&ctl)) q->gfs[i]=0;
               }
               if (!mgh->in_amrh[l]) free(q->gfs);
            }
            for (i=0; i<curr_context->dim; i++) if (q->x[i] && !mgh->in_amrh[l]) pfree(q->x[i]);
            free(q);
         }
         free(lev);
      }
   }

   free(mgh);
   c->MG_cgh=0;

   IFG(4) printf("      destroy_mgh <<\n");
   return;
}

//=============================================================================
// computes the "island number" for each grid --- i.e., each grid within
// a set of connected union of grids is given a number to identify which 
// connected union it is in. If the number is negative, then at least one 
// member of the union is a coarsest grid. If the number is 0, then at least 
// one member is 'smaller' than coarsest. 
//
// actual dx of level = dx*cf
//=============================================================================
void calc_island_no(real *bboxes, int num, real *island_no, real *dx, int cf)
{
   int i,j,k,delta=1,shape[PAMR_MAX_DIM],overlap;
   context *c=curr_context;

   for (i=0; i<num; i++) 
   {
      island_no[i]=i+1;
      for (j=0; j<c->dim; j++)
      {
         shape[j]=(bboxes[2*i*c->dim+2*j+1]-bboxes[2*i*c->dim+2*j])/dx[j]+1.5;
         if ( (((shape[j]-1) % cf) && ((shape[j]-1)/cf+1)>=c->MG_min_cwidth[j]) ||
              (((shape[j]-1)/cf+1)<c->MG_min_cwidth[j]) ) island_no[i]=0;
         else if (((shape[j]-1)/cf/2+1)<c->MG_min_cwidth[j]) island_no[i]=-fabs(island_no[i]);
      }
   }

   while (delta)
   {
      delta=0;
      for (i=0; i<num; i++)
      {
         for (j=i+1; j<num; j++)
         {
            if (island_no[i]!=island_no[j])
            {
               overlap=1;
               for (k=0; k<c->dim; k++)
                  if (fuzz_gt(bboxes[2*i*c->dim+2*k],bboxes[2*j*c->dim+2*k+1],dx[k]/2) ||
                      fuzz_gt(bboxes[2*i*c->dim+2*k+1],bboxes[2*j*c->dim+2*k],dx[k]/2)) overlap=0;
               if (overlap)
               {
                  delta++;
                  if (island_no[i]==0 || island_no[j]==0) { island_no[i]=island_no[j]=0; }
                  else island_no[i]=island_no[j]=min(island_no[i],island_no[j]);
               }
            }
         }
      }
   }
}

//=============================================================================
// a couple of utility routines used by merge_bboxes below
//=============================================================================
real bbox_vol(real *bbox, int dim)
{
   int i;
   real vol=1;

   for (i=0; i<dim; i++) vol*=(bbox[2*i+1]-bbox[2*i]);
   return vol;
}

int bbox_and(real *A, real *B, real *A_and_B, int dim)
{
   int i;

   for (i=0; i<dim; i++) 
   {
      A_and_B[2*i]=max(A[2*i],B[2*i]);
      A_and_B[2*i+1]=min(A[2*i+1],B[2*i+1]);
      if (A_and_B[2*i]>A_and_B[2*i+1]) return 0;
   }
   return 1;
}

void bbox_or(real *A, real *B, real *A_or_B, int dim)
{
   int i;

   for (i=0; i<dim; i++) 
   {
      A_or_B[2*i]=min(A[2*i],B[2*i]);
      A_or_B[2*i+1]=max(A[2*i+1],B[2*i+1]);
   }
}

//=============================================================================
// merges a list of bbox's, as explained in PAMR_merge_bbox_list() header.
//=============================================================================
int merge_bboxes(real *bboxes, int num, real min_eff)
{
   int dim=curr_context->dim;
   int i,j,k;
   int delta=1;
   real A_and_B[2*PAMR_MAX_DIM],A_vol,B_vol,A_and_B_vol;
   real C_vol,A_or_B[2*PAMR_MAX_DIM],A_or_B_vol,eff;

   while(delta)
   {
      delta=0;
      for (i=0; i<(num-1); i++)
      {
         for (j=i+1; j<num; j++)
         {
            A_vol=bbox_vol(&bboxes[2*i*dim],dim);
            B_vol=bbox_vol(&bboxes[2*j*dim],dim);
            if (bbox_and(&bboxes[2*i*dim],&bboxes[2*j*dim],A_and_B,dim))
            {
               A_and_B_vol=bbox_vol(A_and_B,dim);
            }
            else
            {
               A_and_B_vol=0;
            }
            bbox_or(&bboxes[2*i*dim],&bboxes[2*j*dim],A_or_B,dim);
            A_or_B_vol=bbox_vol(A_or_B,dim);
            eff=(A_vol+B_vol-A_and_B_vol)/A_or_B_vol;
            if (fuzz_gte(eff,min_eff,1.0e-5))
            {
               delta++;
               for (k=0; k<2*dim; k++) bboxes[2*i*dim+k]=A_or_B[k];
               for (k=2*j*dim; k<(2*dim*(num-1)); k++) bboxes[k]=bboxes[k+2*dim];
               j--;
               num--;
            }
         }
      }
   }

   return num;
}

//=============================================================================
// the following function calls the user defined mask definition function
// for all local grids from Lmin to Lmax ... no sync'ing is done
//=============================================================================
void init_ex_mask(int Lmin, int Lmax, cgh *gh, int gfn)
{
   context *c=curr_context;
   int L;
   grid *g;

   for (L=Lmin; L<=Lmax; L++)
   {
      if (gh->levels[L-1]) g=gh->levels[L-1]->grids; else g=0;
      while(g)
      {
         if (g->rank==c->rank)
         {
            if (!(g->gfs[gfn-1])) 
            {
               printf("init_ex_mask: error ... mask grid not allocated\n"); exit(1);
            }
            c->app_fill_ex_mask(g->gfs[gfn-1],g->dim,g->shape,g->bbox,c->excised);
         }
         g=g->next;
      }
   }
}

//=============================================================================
// functions for periodic boundary support
//
// make_virtual_copy() makes a virtual copy of a grid, sharing the
// same memory if the grid is local, but *no* perimiter coordinate arrays
// 
// alloc_virtual_grids()replicates all grids that wrap around, and those
// near the "min" boundaries that are within 2*ghostwidth of the boundary, shifting the 
// coordinates as needed, but using the same memory as the central grid.
// 
// in 1D, a grid could get duplicated 2 times
// in 2D,             "               8 times
// in 3D,             "               26 times
//
// free_virtual_grids() deletes all virtual grids in a level
//=============================================================================
grid *make_virtual_copy(grid *g)
{
   grid *ng;
   int j;

   if (!(ng=(struct grid *)malloc(sizeof(struct grid)))) 
      { printf("make_virtual_copy: error ... out of mem\n"); exit(1); }
   
   ng->prev=ng->next=0;
   ng->rank=g->rank;
   ng->ngfs=g->ngfs;
   ng->comm=g->comm;
   ng->dim=g->dim;
   ng->t=g->t;
   ng->gfs=g->gfs;
   ng->coarsest=g->coarsest;
   ng->virtual_g=1;
   for (j=0; j<g->dim; j++)
   {
      ng->bbox[2*j]=g->bbox[2*j];
      ng->bbox[2*j+1]=g->bbox[2*j+1];
      ng->shape[j]=g->shape[j];
      ng->x[j]=0;
      ng->ghost_width[2*j]=g->ghost_width[2*j];
      ng->ghost_width[2*j+1]=g->ghost_width[2*j+1];
      ng->wrap[2*j]=g->wrap[2*j];
      ng->wrap[2*j+1]=g->wrap[2*j+1];
   }

   return ng;
}

void alloc_virtual_grids(level *l)
{
   grid *g,*ng;
   context *c=curr_context;
   int d,d1,d2,d3,s1,s2,s3,g0_wrap[2*PAMR_MAX_DIM];
   real width[PAMR_MAX_DIM];

   for (d=1; d<=c->dim; d++) width[d-1]=(c->bbox[2*(d-1)+1]-c->bbox[2*(d-1)]);

   g=l->grids;
   while(g)
   {
      // we're inserting grids during the loop, so check that we don't use those
      if (!g->virtual_g)
      {
         for (d=1; d<=c->dim; d++)
         {
            if (c->periodic[d-1] &&
                (g->wrap[2*(d-1)] || (int)((g->bbox[2*(d-1)]-c->bbox[2*(d-1)]+0.5)/l->dx[d-1]) < 2*c->ghost_width[d-1])) 
            g0_wrap[2*(d-1)]=1;
            else g0_wrap[2*(d-1)]=0;
            g0_wrap[2*(d-1)+1]=g->wrap[2*(d-1)+1];
         }
         // first all single shifts
         for (d=1; d<=c->dim; d++)
         {
            if (g0_wrap[2*(d-1)])
            {
               ng=make_virtual_copy(g);
               ng->bbox[2*(d-1)]+=width[d-1];
               ng->bbox[2*(d-1)+1]+=width[d-1];
               insert_grid(l,ng);
            }
            if (g0_wrap[2*(d-1)+1])
            {
               ng=make_virtual_copy(g);
               ng->bbox[2*(d-1)]-=width[d-1];
               ng->bbox[2*(d-1)+1]-=width[d-1];
               insert_grid(l,ng);
            }
         }

         // double shifts
         for (d1=1; d1<=c->dim; d1++)
          for (d2=1+d1; d2<=c->dim; d2++)
           for (s1=0; s1<=1; s1++)
            for (s2=0; s2<=1; s2++)
             if (g0_wrap[2*(d1-1)+s1] && g0_wrap[2*(d2-1)+s2])
             {
                ng=make_virtual_copy(g);
                ng->bbox[2*(d1-1)]+=((1-2*s1)*width[d1-1]);
                ng->bbox[2*(d1-1)+1]+=((1-2*s1)*width[d1-1]);
                ng->bbox[2*(d2-1)]+=((1-2*s2)*width[d2-1]);
                ng->bbox[2*(d2-1)+1]+=((1-2*s2)*width[d2-1]);
                insert_grid(l,ng);
             }

         // triple shifts
         if (c->dim==3)
         {
            for (s1=0; s1<=1; s1++)
             for (s2=0; s2<=1; s2++)
              for (s3=0; s3<=1; s3++)
               if (g0_wrap[s1] && g0_wrap[2+s2] && g0_wrap[4+s3])
               {
                  ng=make_virtual_copy(g);
                  ng->bbox[0]+=((1-2*s1)*width[0]);
                  ng->bbox[1]+=((1-2*s1)*width[0]);
                  ng->bbox[2]+=((1-2*s2)*width[1]);
                  ng->bbox[3]+=((1-2*s2)*width[1]);
                  ng->bbox[4]+=((1-2*s3)*width[2]);
                  ng->bbox[5]+=((1-2*s3)*width[2]);
                  insert_grid(l,ng);
               }
         }
      }
      g=g->next;
   }
}

void free_virtual_grids(level *l)
{
   grid *g,*ng;

   g=l->grids;

   while(g)
   {
      ng=g->next;
      if (g->virtual_g)
      {
         if (ng) ng->prev=g->prev;
         if (g->prev) g->prev->next=ng; else l->grids=ng;
         free(g);
      }
      g=ng;
   }
}
