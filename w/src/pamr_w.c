//=============================================================================
//  pamr.c --- parallel AMR library ... API functions
//
//  V1. Copyright 2002-2007 F.Pretorius
//  V2. Copyright 2008-     F.Pretorius & B.Stephens
//=============================================================================

#include <stdio.h>
#include <bbhutil.h>
#include "pamr_w.h"
#include "gh_w.h"
#include "misc_w.h"
#include "transfer_w.h"
#include "io_w.h"
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

//=============================================================================
// The following routine sets the trace level:
//
// 0: no debugging info
// 1: prints API program flow info
// 2: 1 + some internal program flow
// 3: 2 + outputs debug sdfs
// 4: 3 + prints some local info + more internal flow
//=============================================================================
void PAMR_set_trace_lev(int lev)
{
   gtrace=lev;

   IFG(1) printf(">> PAMR_set_trace_lev(%i) <<\n",lev);
}

//=============================================================================
// initialize a context ... returns the context (>0), or 0 if an error 
// occured. All arguments passed here are required, and cannot be altered
// ... other options are set to defaults, and can be changed later.
//
// [q,q_size] ---  define memory block to use for grid function allocation 
//                 (not implemented yet!)
// dim --- spatial dimension
// shape[dim] --- shape of base level (1)
// bbox[2*dim] --- bounding box of base level (1)
//
// cp_file is a checkpoint file to load a context from disk.
//
// The c_tag/v_tag variables specify how cell-centered/vertex centered variables  
// should be defined.                                                             
//    if c_tag & v_tag are null, then *all* variables are vertex centered         
//    -- this is the default                                                      
//                                                                                
//    if c_tag (v_tag) is null but v_tag (c_tag) is not, then                     
//    then all variables ending in v_tag (c_tag) are vertex (cell)                
//    centered, and all the rest are cell (vertex) centered.                      
//                                                                                 
//    if both c_tag and v_tag are not null, then variables ending in c_tag        
//    (v_tag) are cell (vertex) centered, and variables without either            
//    appendage are vertex centered.                                              
//=============================================================================
int PAMR_init_context(real *q, int q_size, int dim, int *shape, real *bbox, 
                      char *cp_file, char *v_tag, char *c_tag)
{
   static int first=1;
   int i,cnt,ret=-10;
   int my_rank,mpi_size;
   int *levs,k,cnode,cnode_max,cp_err,g_cp_err,ltrace=0;
   real *p;

   IFG(1) printf(">> PAMR_init_context(q=%p,q_size=%i,dim=%i,...) \n",q,q_size,dim);

   if (first) 
   { 
      if (v_tag) PAMR_v_tag=strdup((const char *)v_tag); else PAMR_v_tag=0;
      if (c_tag) PAMR_c_tag=strdup((const char *)c_tag); else PAMR_c_tag=0;

      for (i=0;i<PAMR_MAX_CONTEXTS;i++) contexts[i]=0; 
      first=0; 
   }

   if (dim<1 || dim>PAMR_MAX_DIM) 
   {
      printf("PAMR_init_context: invalid dimension ... must be between 1 and %i\n",PAMR_MAX_DIM);
      return 0;
   }

   if (q || q_size)
   {
      printf("PAMR_init_context: external-memory-block option currently not supported\n");
      return 0;
   }

   i=0;
   while (contexts[i] && i<PAMR_MAX_CONTEXTS) i++;
   if (i==PAMR_MAX_CONTEXTS)
   {
      printf("PAMR_init_context: error ... PAMR_MAX_CONTEXTS=%i reached\n",PAMR_MAX_CONTEXTS);
      return 0;
   }

   MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   cnt=i+1;
   if (cp_file)
   {
      curr_context=0;
      cp_err=g_cp_err=0;
      if (PAMR_parallel_io) cnode_max=1; else cnode_max=mpi_size;
      for (cnode=0; cnode<cnode_max; cnode++)
      {
         if (PAMR_parallel_io || my_rank==cnode)
         {
            if (!(PAMR_do_cp(cp_file,my_rank,mpi_size,PAMR_CP_DIR_RESTORE))) 
            {
               printf("PAMR_init_context: error ... unable to restore state from file %s\n",cp_file);
               g_cp_err=cp_err=1;
            }
         }
         MPI_Allreduce(&cp_err,&g_cp_err,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
         if (!PAMR_parallel_io && ltrace && my_rank==0) printf(" ... PAMR_cp finished on node %i\n",cnode);
      }
      if (g_cp_err) return 0;
      contexts[cnt-1]=curr_context;
      // ======================================================================
      // recompose level by level, to keep the same t
      // ======================================================================
      for (i=0; i<PAMR_MAX_LEVS; i++) if (curr_context->num_sgh_bboxes[i])
      {
         IFL printf("my_rank=%i, recomposing level %i\n",my_rank,i+1);
         if (!(levs=(int *)malloc(sizeof(int)*curr_context->num_sgh_bboxes[i])))
         {
            printf("PAMR_init_context: error ... out of memory%s\n",cp_file);
            PAMR_free_context(cnt);
            return 0;
         }
         for (k=0; k<curr_context->num_sgh_bboxes[i]; k++) levs[k]=i+1;

         // PAMR_compose_hierarchy() allocates new memory for the sgh bbox list, 
         // hence the following:
         p=curr_context->sgh_bboxes[i];
         curr_context->sgh_bboxes[i]=0;
         if (!(PAMR_compose_hierarchy(i+1,i+1,curr_context->num_sgh_bboxes[i],
               levs,p,curr_context->curr_cgh->levels[i]->t)))
         {
            printf("PAMR_init_context: error recomposing level %i of restored hierarchy\n",i);
            PAMR_free_context(cnt);
            return 0;
         }
         free(levs); free(p);
      }
   }
   else
   {
      if (!(contexts[cnt-1]=(struct context *)malloc(sizeof(context))))
      {
         printf("PAMR_init_context: error ... out of memory\n");
         return 0;
      }
   
      curr_context=contexts[cnt-1];

      curr_context->q=q;
      curr_context->q_size=q_size;
      curr_context->dim=dim;
      curr_context->num_vars=0;
      curr_context->vars=0;
      curr_context->gfns=0;
      curr_context->tf_bits=0;
      curr_context->frozen=0;
      curr_context->curr_cgh=0;
      curr_context->MG_cgh=0;
      curr_context->n_rank=0;
      curr_context->iter_g=0;
      curr_context->iter_g_stack_top=0;
      curr_context->interp_buffer=2; // OK for 4th order and lower interpolation
      curr_context->excision_on=0;
      curr_context->rank=my_rank;
      curr_context->size=mpi_size;

      for (i=0; i<dim; i++) 
      {
         curr_context->shape[i]=shape[i];
         if (shape[i]<=1)
         {
            printf("PAMR_init_context: error ... invalid shape[%i]=%i (must be > 1)\n",i,shape[i]);
            ret=0; goto clean_up;
         }
         curr_context->bbox[i*2]=bbox[2*i];
         curr_context->bbox[i*2+1]=bbox[2*i+1];
         curr_context->ghost_width[i]=2;
         curr_context->periodic[i]=0;
         curr_context->min_width[i]=5;
         curr_context->MG_min_cwidth[i]=5;
      }

      for (i=0; i<PAMR_MAX_LEVS; i++)
      {
         curr_context->rho_sp[i]=2;
         curr_context->rho_tm[i]=2;
         curr_context->sgh_bboxes[i]=0;
         curr_context-> num_sgh_bboxes[i]=0;
      }
   
      curr_context->lambda=1;
      curr_context->top_xy=curr_context->top_xz=curr_context->top_yz=1;

      curr_context->gfn_var_type=0;
   
      recalc_lev_info();
   }

   IFG(1) printf("  PAMR_init_context, ret=%i <<\n",cnt);

   return cnt;

clean_up:
   free(curr_context);
   curr_context=contexts[cnt-1]=0;
   return ret;
}

//=============================================================================
// deletes a context and any associated memory that has been allocated
//=============================================================================
void PAMR_free_context(int num)
{
   IFG(1) printf(">> PAMR_free_context(num=%i) \n",num);

   if (num<1 || num>PAMR_MAX_CONTEXTS || (contexts[num-1]==0))
   {
      printf("PAMR_free_context: error ... invalid context number %i\n",num);
      return;
   }

   free_context_mem(); 
   free(contexts[num-1]);
   contexts[num]=0;
   curr_context=0;

   IFG(1) printf("   PAMR_free_context <<\n");
}

//=============================================================================
// changes the active context --- returns sucess or failure
// 
// disabled for now ... grid function tags (cell/vertex) are currently
// stored as globals, and will need to be remembered to make the
// switch work (did not add the tags to the context structure
// to avoid making a new cp version). 
//
// But since multiple contexts have never been
// used/tested in any case, this is not much of an issue.
// 
//=============================================================================
int switch_context(int num)
{
   IFG(1) printf(">> switch_context(num=%i) \n",num);

   printf("switch_context: multiple contexts not yet supported\n"); 
   
   IFG(1) printf("   switch_context <<\n");
   return 0;

   if (num<1 || num>PAMR_MAX_CONTEXTS || (contexts[num-1]==0))
   {
      printf("switch_context: error ... invalid context number %i\n",num);
      return 0;
   }

   curr_context=contexts[num-1];

   IFG(1) printf("   switch_context <<\n");
   return 1;
}

//=============================================================================
// The follow function is called at the initial time and during regridding
// to construct the hierarchy.
//
// levels from min_lev to max_lev inclusive are to be replaced 
// via the set of num grids defined via bbox[2*dim*num] and lev[num], all at
// time t.
//
// returns 1 for success, 0 for failure
//
// NOTE: time-levels are *NOT* considered during interpolation ... i.e.,
//       variables that are initialized from the old hierarchy are
//       only interpolated in space, and hence it is the users task
//       to do the temporal interpolation (if needed) afterwards.
//=============================================================================
int PAMR_compose_hierarchy(int min_lev, int max_lev, int num, int *lev, real *bbox, real t)
{
   cgh *new_sgh=0,*new_cgh=0;
   context *c=curr_context;
   int clmin,clmax,i;
   var *f;

   IFG(1) printf(">> PAMR_compose_hierarchy(min_lev=%i,max_lev=%i,num=%i,...,t=%lf) \n",min_lev,max_lev,num,t);
   
   if (!(c)) { printf("PAMR_compose_hierarchy: error ... no current context\n"); return 0; }
   if (c->MG_cgh) { printf("PAMR_compose_hierarchy: error ... cannot re-compose with an existing MG hierarhcy\n"); return 0; }

   c->n_rank=0;

   //--------------------------------------------------------------------------
   // If first call for this context, allocate the transfer bits
   // also, calculate which variables are 
   //--------------------------------------------------------------------------
   if (!(c->tf_bits))
   {
      if (!(c->tf_bits=(int *)malloc(sizeof(int)*c->gfns)))
      {
         printf("PAMR_compose_hierarchy: error ... out of memory\n"); return 0;
      }
      PAMR_clear_tf_bits();
   }
   if (!(c->gfn_var_type))
   {
      if (!(c->gfn_var_type=(int *)malloc(sizeof(int)*c->gfns)))
      {
         printf("PAMR_compose_hierarchy: error ... out of memory\n"); return 0;
      }
      for (i=0; i<c->gfns; i++) c->gfn_var_type[i]=PAMR_var_type(find_name(i+1));
   }
   //--------------------------------------------------------------------------
   // build_sgh converts the bbox list into a 'sequential' grid hierarchy
   // ... i.e., the hierarchy that would be used in a non-parallel code.
   // no data is allocated.
   // Also, save bboxes for future reference (when building a MGH) when
   //--------------------------------------------------------------------------
   if (!(save_sgh_bboxes(min_lev,max_lev,num,lev,bbox))) return 0;

   if (!(new_sgh=build_sgh(min_lev,max_lev,num,lev,bbox,t))) return 0;

   //--------------------------------------------------------------------------
   // if desired, eliminated grid-overlap
   //--------------------------------------------------------------------------
   if (c->gdm & PAMR_GDM_NO_OVERLAP) new_sgh=sgh_cull_overlap(new_sgh);

   //--------------------------------------------------------------------------
   // now split the sgh into a 'computational' grid hierarchy, i.e. split
   // it across the nodes of the network.
   // compose_cgh also allocates local gf memory.
   //--------------------------------------------------------------------------
   if (!(new_cgh=compose_cgh(new_sgh,0))) {free_cgh(new_sgh,1,PAMR_MAX_LEVS); return 0;}
   free_cgh(new_sgh,1,PAMR_MAX_LEVS);
   
   //--------------------------------------------------------------------------
   // transfer/interpolate data from an existing cgh, but first define excision
   // mask for new cgh if desired
   //--------------------------------------------------------------------------
   if (c->excision_on) init_ex_mask(min_lev,max_lev,new_cgh,c->amr_mask_gfn,c->amr_mask_c_gfn);

   if (c->curr_cgh) init_new_cgh(new_cgh,c->curr_cgh);

   //--------------------------------------------------------------------------
   // replace levels of curr_cgh with those of new_cgh, freeing old resources.
   //--------------------------------------------------------------------------
   if (c->curr_cgh) 
   {
      clmin=c->curr_cgh->min_lev;
      clmax=c->curr_cgh->max_lev;
      if (!(incorporate_new_cgh(new_cgh,c->curr_cgh))) c->curr_cgh=0;
   }
   else c->curr_cgh=new_cgh;

   c->iter_g=0; // reset any iterator

   //--------------------------------------------------------------------------
   // alter flags of all variables that changed status to reflect new status
   //--------------------------------------------------------------------------
   f=c->vars;
   while(f)
   {
      if (f->status==VAR_STATUS_TURN_ON) 
      {
         f->status=VAR_STATUS_ON;
         if (min_lev>clmin || max_lev<clmax) 
         {
            printf("PAMR_compose_hierarchy: error ... variable status changed to ON "
                   "but recompose not performed over entire hierarchy\n"); 
            exit(1); 
         }
      }
      if (f->status==VAR_STATUS_TURN_OFF && min_lev==clmin && max_lev==clmax) f->status=VAR_STATUS_OFF;
      f=f->next;
   }

   IFG(1) printf("   PAMR_compose_hierarchy <<\n");
   return 1;
}

//=============================================================================
// The following builds the multigrid hierarchy from time-level tl AMR grids, 
// spanning the given set of levels. If successfull return the maximum
// level in the MGH (corresponding to max_lev in the AMR hierarchy).
//=============================================================================
int PAMR_build_mgh(int min_lev, int max_lev, int tl)
{
   context *c=curr_context;
   cgh *amrh=curr_context->curr_cgh;
   real *mgh_bbox;
   int *mgh_lev;
   int *mgh_coarsest;
   int l,ret;

   IFG(1) printf(">> PAMR_build_mgh(min_lev=%i,max_lev=%i) \n",min_lev,max_lev);

   if (!(c)) { printf("PAMR_build_mgh: error ... no current context\n"); return 0; }
   if (c->MG_cgh) { printf("PAMR_build_mgh: error ... MG hierarhcy already exists\n"); return 0; }
   if (!(amrh)) { printf("PAMR_build_mgh: error ... AMR hierarchy does not exist\n"); return 0; }

   for (l=min_lev-1; l<max_lev; l++) 
   {
      if (!(amrh->levels[l])) { printf("PAMR_build_mgh: error ... level %i missing in AMR hierarchy\n",l+1); return 0; }
   }

   if (build_mgh(min_lev,max_lev,tl)) ret=c->MG_cgh->max_lev; else ret=0;

   IFG(1) printf("   PAMR_build_mgh <<\n");
   return ret;
}

//=============================================================================
// deallocate a prior MGH
//=============================================================================
int PAMR_destroy_mgh()
{
   context *c=curr_context;
   cgh *amrh=curr_context->curr_cgh;
   real *mgh_bbox;
   int *mgh_lev;
   int *mgh_coarsest;
   int l,ret;

   IFG(1) printf(">> PAMR_destroy_mgh\n");

   if (!(c)) { printf("PAMR_destroy_mgh: error ... no current context\n"); return 0; }

   destroy_mgh();

   IFG(1) printf("   PAMR_destroy_mgh <<\n");
   return 0;
}

//=============================================================================
// This is a utility routine to help users construct new grid-hierarchies
// from local truncation error estimates (for instance)
//
// global_bbox is preallocated memory of size global_num*sizeof(bbox)'s 
// ---> upon completion, global_bbox will be filled with the global list of 
// bounding boxes, and global_num will be set to the number of them.
//
// the "efficiency" of joining two bounding boxes A and B into a single,
// encompassing box C is defined as
//
// eff = vol(A union B)/vol(C)
//
// min_eff is the minimum required efficiency ( 0 <= min_eff <= 1)
//
// returns zero if fails (e.g., not enough memory given)
//
// NOTE: merge does *not* take periodic boundaries into account.
//       The reason for this is that, for the purposes of grid distribution,
//       grids can never be 'merged' across periodic boundaries. 
//       However note that PAMR won't communicate periodic information
//       for wrapping grids that don't share the same dimensions along
//       the abutting surface, so some post-processing of the list
//       may need to be made in that case.
//=============================================================================
int PAMR_merge_bboxes(real *local_bbox, int local_num, real *global_bbox, int *global_num, real min_eff)
{
   context *c=curr_context;
   int max_num=*global_num,num=0,ret=1,i;
   int local_nums[PAMR_MAX_NODES],displs[PAMR_MAX_NODES],tot_size;
   real *recvbuf;

   IFG(1) printf(">> PAMR_merge_bboxes(min_eff=%lf) \n",min_eff);

   if (!(c)) { printf("PAMR_merge_bboxes: error ... no current context\n"); return 0; }

   //--------------------------------------------------------------------------
   // First, broadcast to everyone how large each local list is
   //--------------------------------------------------------------------------
   MPI_Allgather(&local_num,1,MPI_INT,local_nums,1,MPI_INT,MPI_COMM_WORLD);

   num=0;
   for (i=0; i<c->size; i++) 
   {
      displs[i]=num*2*c->dim;
      num+=local_nums[i];
      local_nums[i]*=2*c->dim; 
   }

   if (num>(*global_num))
   {
      printf("PAMR_merge_bboxes: error ... global_bbox needs to hold at least %i bboxes\n",num);
      *global_num=0;
      return 0;
   }

   //--------------------------------------------------------------------------
   // Now, broadcast our bbox list to everyone.
   //--------------------------------------------------------------------------
   if (num)
   {
      MPI_Allgatherv(local_bbox,local_num*2*c->dim,MPI_DOUBLE,global_bbox,local_nums,displs,MPI_DOUBLE,MPI_COMM_WORLD);
      *global_num=merge_bboxes(global_bbox,num,min_eff);
   }
   else *global_num=0;

   IFG(1) printf("   PAMR_merge_bboxes <<\n");
   return ret;
}

//=============================================================================
// A utility routine that fills the bbox array (up to *num on input) with
// the sgh of level L ... returns 0 if level L is empty. Upon return *n is
// set to the actual number of bbox elements filled in.
//=============================================================================
int PAMR_get_sgh(real *bbox, int *num, int L)
{
   context *c=curr_context;
   int ret=0,i,n;

   IFG(1) printf(">> PAMR_get_sgh(lev=%i) \n",L);

   if (!(c)) { printf("PAMR_get_sgh: error ... no current context\n"); return 0; }

   if (L>PAMR_MAX_LEVS || L<0) { printf("PAMR_get_sgh: error ... level out of range\n"); return 0; }

   if (c->sgh_bboxes[L-1])
   {
      ret=1;

      n=c->num_sgh_bboxes[L-1];
      if (n>*num) n=*num;
      *num=n;

      for (i=0; i<c->dim*2*n; i++) bbox[i]=(c->sgh_bboxes[L-1])[i];
   }

   IFG(1) printf("   PAMR_get_sgh <<\n");
   return ret;
}

//=============================================================================
// called by sync/inter/inject below ... interior is passed
// to build_owned_gsl, and controls whether the ghost zones are included in 
// the definition of 'interior' or not.
// --> NOTE: see alloc_gsl() and build_owned_gsl() for latest usage of interior flag
//
// 'Owned' Grid boundaries do not overlap ... this could cause problems
// for interpolation, as the seperation between owned grids on the
// coarse level causes gaps of rho_sp grids on the fine level. To
// overcome this (with extend_owned=1), we extend the owned grid segments
// by one point along +ve and -ve directions on the coarse grids. This will cause
// a slighlty larger amount of data to be communicated, but I think the
// cost will be minimal. The alternative would be to try to reconstruct
// the coarse grids after communication ... however this may require
// additional communication in any case if the fine level is split along
// the same lines as the coarse level (for then both pieces of the fine
// grid will need to know the same info from a coarse level, i.e. a
// one-to-many mapping).
// Also, this method will require that the user specify a sufficient
// ghostwidth (AND interp_buffer) for the interpolation order if interior 
// stencils should be used at all sequential-grid interior points. This, 
// I don't think is a serious restriction, as standard FD stencils will 
// give the correct ghost width automatically.
//
// PAMR V2 UPDATE: 'extend_owned' modified for CC structures ... see discussion
// in transfer.c
//
// if (AMR_bdy_only>0), then copy to a region of size 'AMR_bdy_only' along
// the AMR boundaries only
//
// The other API-level communication function is init_new_cgh(), called by 
// compose_hierarchy()
//=============================================================================
int communicate(int l1, int l2, int hierarchy, int interior, 
                int extend_owned, int AMR_bdy_only, char *debug_tag)
{
   gsl *src[PAMR_MAX_NODES],*dst=0;
   gsl *transfer_src[PAMR_MAX_NODES],*transfer_dst[PAMR_MAX_NODES];
   context *c=curr_context;
   cgh *gh;
   int ret=0,i;
   char db_name[256];
   real t;

   IFG(1) printf(">> communicate, l1 = %d, l2 = %d \n", l1,l2);

   switch(hierarchy)
   {
      case PAMR_AMRH:
         if (!(gh=c->curr_cgh))
         { printf("communicate: error ... AMR hierarhcy doesn't exist\n"); return 0; }
         break;
      case PAMR_MGH:
         if (!(gh=c->MG_cgh))
         { printf("communicate: error ... MG hierarhcy doesn't exist\n"); return 0; }
         break;
      default:
         printf("communicate: error ... invalid hierarchy\n"); return 0;
   }

   if (!gh->levels[l2-1] || !gh->levels[l1-1])
      { printf("communicate: error ... level does not exist\n"); return 0; }

   t=gh->levels[l2-1]->t;

   for (i=0; i<c->size; i++) src[i]=transfer_src[i]=transfer_dst[i]=0;

   // support for periodic boundaries
   alloc_virtual_grids(gh->levels[l1-1]);

   for (i=0; i<c->size; i++) src[i]=build_owned_gsl(i,gh,l1,interior,extend_owned); 
   IFG(3) {sprintf(db_name,"%s_owned_gsl_%i",debug_tag,c->rank); debug_save_gsl(db_name,t,src[c->rank],-1);}

   if (!(dst=build_complete_gsl(gh,l2,AMR_bdy_only))) goto fin;
   IFG0(3) {sprintf(db_name,"%s_complete_gsl",debug_tag); debug_save_gsl(db_name,t,dst,-1);}

   for (i=0; i<c->size; i++) build_gstl(src[i],dst,&transfer_src[i],&transfer_dst[i]);
   IFG(3) {sprintf(db_name,"%s_transfer_src_%i",debug_tag,c->rank); debug_save_gsl(db_name,t,transfer_src[c->rank],-1);}
   IFG(3) {sprintf(db_name,"%s_transfer_dst_%i",debug_tag,c->rank); debug_save_gsl(db_name,t,transfer_dst[c->rank],-1);}

   transfer(transfer_src,transfer_dst);

   ret=1;

fin:
   if (dst) free_gsl(dst);
   for (i=0; i<c->size; i++)
   {
      if (src[i]) free_gsl(src[i]);
      if (transfer_src[i]) free_gsl(transfer_src[i]);
      if (transfer_dst[i]) free_gsl(transfer_dst[i]);
   }
   free_virtual_grids(gh->levels[l1-1]);

   IFG(1) printf("   communicate << \n");
   return ret;
}

//=============================================================================
// PAMR_sync(level,time levels (-1 for all), hierarchy) synchronizes
// grids at level l, communicating data in the grip-overlap regions
//
// AMR_bdy_width, if > 0, excludes a zone of size 'AMR_bdy_width' about 
// AMR boundaries from ownership
//=============================================================================
int PAMR_sync(int l, int tl, int hierarchy, int AMR_bdy_width)
{
   context *c=curr_context;
   int ret=0;

   IFG(1) printf(">> PAMR_sync(l=%i,tl=%i,hierarchy=%i)\n",l,tl,hierarchy);

   if (!(c)) { printf("PAMR_sync: error ... no current context\n"); return 0; }
   if (l<1 || l>PAMR_MAX_LEVS) { printf("PAMR_sync: error ... level %i out of range\n",l); return 0; }
   if (AMR_bdy_width<0) { printf("PAMR_sync: error ... AMR_bdy_width must be >=0\n"); return 0; }

   PAMR_set_tf_bits(tl,hierarchy,PAMR_TF_SYNC);
   ret=communicate(l,l,hierarchy,1+AMR_bdy_width,0,0,"sync");

   IFG(1) printf("   PAMR_sync <<\n");
   return ret;
}
   
//=============================================================================
// PAMR_interp(level,time levels (-1 for all), hierarchy) interpolates
// from coarse level lc to fine level lc+1, in the grid overlap region.
//=============================================================================
int PAMR_interp(int lc, int tl, int hierarchy)
{
   context *c=curr_context;
   int ret=0;

   IFG(1) printf(">> PAMR_interp(lc=%i,tl=%i,hierarchy=%i)\n",lc,tl,hierarchy);

   if (!(c)) { printf("PAMR_interp: error ... no current context\n"); return 0; }
   if (lc<1 || lc>=PAMR_MAX_LEVS) { printf("PAMR_interp: error ... level %i out of range\n",lc); return 0; }

   PAMR_set_tf_bits(tl,hierarchy,PAMR_TF_INTERP);
   ret=communicate(lc,lc+1,hierarchy,1,1,0,"interp");

   IFG(1) printf("   PAMR_interp <<\n");
   return ret;
}

//=============================================================================
// PAMR_bdy_interp(level,time levels (-1 for all)) interpolates
// from coarse level AMR boundaries at level lc to fine level lc+1, 
// in the grid overlap region (of size bdy_width on the fine grid).
// This one handles vertex centered variables only.
// AMR_bdy_width means the number of vertices across the boundary.
//=============================================================================
int PAMR_AMR_bdy_interp(int lc, int tl, int AMR_bdy_width)
{
   context *c=curr_context;
   int ret=0;

   IFG(1) printf(">> PAMR_AMR_bdy_interp(lc=%i,tl=%i)\n",lc,tl);

   if (!(c)) { printf("PAMR_AMR_bdy_interp: error ... no current context\n"); return 0; }
   if (lc<1 || lc>=PAMR_MAX_LEVS) { printf("PAMR_AMR_bdy_interp: error ... level %i out of range\n",lc); return 0; }
   if (AMR_bdy_width<1) { printf("PAMR_AMR_bdy_interp: error ... AMR_bdy_width must be > 0\n"); return 0; }

   PAMR_set_tf_bits(tl,PAMR_AMRH,PAMR_TF_BDY_INTERP);
   ret=communicate(lc,lc+1,PAMR_AMRH,1,1,AMR_bdy_width,"interp");

   IFG(1) printf("   PAMR_AMR_bdy_interp <<\n");
   return ret;
}

//=============================================================================
// PAMR_bdy_interp(level,time levels (-1 for all)) interpolates
// from coarse level AMR boundaries at level lc to fine level lc+1, 
// in the grid overlap region (of size bdy_width on the fine grid).
// This one does cell centered variables only.
// AMR_bdy_width means the number of cell widths across the boundary.
//=============================================================================
int PAMR_AMR_bdy_interp_c(int lc, int tl, int AMR_bdy_width)
{
   context *c=curr_context;
   int ret=0;

   IFG(1) printf(">> PAMR_AMR_bdy_interp(lc=%i,tl=%i)\n",lc,tl);

   if (!(c)) { printf("PAMR_AMR_bdy_interp: error ... no current context\n"); return 0; }
   if (lc<1 || lc>=PAMR_MAX_LEVS) { printf("PAMR_AMR_bdy_interp: error ... level %i out of range\n",lc); return 0; }
   if (AMR_bdy_width<1) { printf("PAMR_AMR_bdy_interp: error ... AMR_bdy_width must be > 0\n"); return 0; }

   PAMR_set_tf_bits(tl,PAMR_AMRH,PAMR_TF_BDY_INTERP_C);
   ret=communicate(lc,lc+1,PAMR_AMRH,1,1,AMR_bdy_width,"interp");

   IFG(1) printf("   PAMR_AMR_bdy_interp <<\n");
   return ret;
}
   
//=============================================================================
// the following functions transfer cell centered variables to
// their vertex centered counterparts (if a pair exists), at a given
// time level, level and in the specified hierarchy (tl=0 for MGH)
//=============================================================================
int PAMR_c_to_v(int l, int tl, int hierarchy, int AMR_bdy_width)
{
   context *c=curr_context;
   int ret=0;

   IFG(1) printf(">> PAMR_c_to_v(l=%i,tl=%i,hier=%i)\n",l,tl,hierarchy);

   if (PAMR_c_tag || PAMR_v_tag) 
   { 
      PAMR_set_tf_bits(tl,hierarchy,PAMR_TF_C_TO_V);
      ret=copy_cv_to_vc(l,hierarchy,PAMR_TF_C_TO_V,AMR_bdy_width,0);
   }
   else IFG(1) printf("   (null operation ... no PAMR_c_tag or PAMR_v_tag defined)\n");

   IFG(1) printf("   PAMR_c_to_v <<\n");
   return ret;
}

int PAMR_v_to_c(int l, int tl, int hierarchy, int AMR_bdy_width)
{
   context *c=curr_context;
   int ret=0;

   IFG(1) printf(">> PAMR_v_to_c(l=%i,tl=%i,hier=%i)\n",l,tl,hierarchy);

   if (PAMR_c_tag || PAMR_v_tag) 
   { 
      PAMR_set_tf_bits(tl,hierarchy,PAMR_TF_V_TO_C);
      ret=copy_cv_to_vc(l,hierarchy,PAMR_TF_V_TO_C,AMR_bdy_width,0);
   }
   else IFG(1) printf("   (null operation ... no PAMR_c_tag or PAMR_v_tag defined)\n");

   IFG(1) printf("   PAMR_v_to_c <<\n");
   return ret;
}

//=============================================================================
// local version of above ... specifically, assumes an iterator has
// already been set up and points to a valid grid, and will only
// do the operation for time level tl grid functions on this grid. Also, *no* 
// syncing is done afterwards. "ret" will always be one if copy_cv_to_vc
// is called ... i.e., if no iterator is set-up, nothing will happend,
// but it won't be flagged as an "error".
//=============================================================================
int PAMR_c_to_v_local(int tl, int hierarchy)
{
   context *c=curr_context;
   int ret=0;
 
   IFG(1) printf(">> PAMR_c_to_v_local(tl=%i)\n",tl);

   if (PAMR_c_tag || PAMR_v_tag) 
   { 
      PAMR_set_tf_bits(tl,hierarchy,PAMR_TF_C_TO_V);
      ret=copy_cv_to_vc(0,hierarchy,PAMR_TF_C_TO_V,0,1);
   }
   else IFG(1) printf("   (null operation ... no PAMR_c_tag or PAMR_v_tag defined)\n");

   IFG(1) printf("   PAMR_c_to_v_local <<\n");
   return ret;
}

int PAMR_v_to_c_local(int tl, int hierarchy)
{
   context *c=curr_context;
   int ret=0;

   IFG(1) printf(">> PAMR_v_to_c_local(tl=%i)\n",tl);

   if (PAMR_c_tag || PAMR_v_tag) 
   { 
      PAMR_set_tf_bits(tl,hierarchy,PAMR_TF_V_TO_C);
      ret=copy_cv_to_vc(0,0,PAMR_TF_V_TO_C,0,1);
   }
   else IFG(1) printf("   (null operation ... no PAMR_c_tag or PAMR_v_tag defined)\n");

   IFG(1) printf("   PAMR_v_to_c_local <<\n");
   return ret;
}

//=============================================================================
// inject(level,time levels (-1 for all), hierarchy) injects grid functions
// from fine level lf to level lf-1
//============================================================================
int PAMR_inject(int lf, int tl, int hierarchy)
{
   context *c=curr_context;
   int ret=0;

   IFG(1) printf(">> PAMR_inject(lf=%i,tl=%i,hierarchy=%i)\n",lf,tl,hierarchy);

   if (!(c)) { printf("PAMR_inject: error ... no current context\n"); return 0; }
   if (lf<2 || lf>PAMR_MAX_LEVS) { printf("PAMR_inject: error ... level %i out of range\n",lf); return 0; }
   
   if (hierarchy==PAMR_MGH && !c->MG_cgh->in_amrh[lf-2])
      PAMR_set_tf_bits(tl,hierarchy,PAMR_TF_INJECT_TO_MG_LEV);
   else
      PAMR_set_tf_bits(tl,hierarchy,PAMR_TF_INJECT);
   ret=communicate(lf,lf-1,hierarchy,1,0,0,"inject");

   IFG(1) printf("   PAMR_inject <<\n");
   return ret;
}

//=============================================================================
// The following set of functions define routines to access elements of
// the data structures
//=============================================================================
void PAMR_set_lambda(real lambda)
{
   IFG(1) printf(">> PAMR_set_lambda\n");

   if (!(curr_context)) { printf("PAMR_set_lambda: error ... no current context\n"); return; }
   curr_context->lambda=lambda;

   recalc_lev_info();

   IFG(1) printf("   PAMR_set_lambda <<\n");
   return;
}

void PAMR_get_lambda(real *lambda)
{
   IFG(1) printf(">> PAMR_get_lambda \n");

   if (!(curr_context)) { printf("PAMR_get_lambda: error ... no current context\n"); return; }
   *lambda=curr_context->lambda;

   IFG(1) printf("   PAMR_get_lambda <<\n");
   return;
}

void PAMR_set_rho(int *rho_sp, int *rho_tm, int num)
{
   int i;

   IFG(1) printf(">> PAMR_set_rho\n");

   if (!(curr_context)) { printf("PAMR_set_rho: error ... no current context\n"); return; }
   if (curr_context->curr_cgh) { printf("PAMR_set_rho: warning ... cannot change rho_sp after "
                                        " the initial hierarchy has been composed\n"); }
   if (num>PAMR_MAX_LEVS) { printf("PAMR_set_rho: error ... number of levels specified (%i) "
                                "is greater than the maximum allowed (%i)\n",num,PAMR_MAX_LEVS); return; }

   for (i=0; i<num; i++)
   {
      if (curr_context->curr_cgh && curr_context->rho_sp[i]!=rho_sp[i]) 
         printf("PAMR_set_rho: warning ... cannot change rho_sp after "
                " the initial hierarchy has been composed\n"); 
      else curr_context->rho_sp[i]=rho_sp[i];
      curr_context->rho_tm[i]=rho_tm[i];
   }

   recalc_lev_info();
   
   IFG(1) printf("   PAMR_set_rho <<\n");
   return;
}

void PAMR_get_rho(int *rho_sp, int *rho_tm, int num)
{
   int i;

   IFG(1) printf(">> PAMR_get_rho\n");

   if (!(curr_context)) { printf("PAMR_get_rho: error ... no current context\n"); return; }
   if (num>PAMR_MAX_LEVS) { printf("PAMR_get_rho: error ... number of levels specified (%i) "
                                "is greater than the maximum allowed (%i)\n",num,PAMR_MAX_LEVS); return; }

   for (i=0; i<num; i++)
   {
      rho_sp[i]=curr_context->rho_sp[i];
      rho_tm[i]=curr_context->rho_tm[i];
   }

   IFG(1) printf("   PAMR_get_rho <<\n");
   return;
}

void PAMR_get_dxdt(int lev, real *dx, real *dt)
{
   int i;

   IFG(1) printf(">> PAMR_get_dxdt\n");

   if (!(curr_context)) { printf("PAMR_get_dxdt: error ... no current context\n"); return; }
   if (lev>PAMR_MAX_LEVS) { printf("PAMR_get_dxdt: error ... level specified (%i) "
                                "is greater than the maximum allowed (%i)\n",lev,PAMR_MAX_LEVS); return; }

   *dt=curr_context->dt[lev-1];
   for (i=0; i<curr_context->dim; i++) dx[i]=curr_context->dx[lev-1][i];

   IFG(1) printf("   PAMR_get_dxdt <<\n");
   return;
}

void PAMR_set_top_ratios(real top_xy, real top_xz, real top_yz)
{
   IFG(1) printf(">> PAMR_set_top_ratios \n");

   if (!(curr_context)) { printf("PAMR_set_top_ratios: error ... no current context\n"); return; }

   curr_context->top_xy=top_xy;
   curr_context->top_xz=top_xz;
   curr_context->top_yz=top_yz;

   IFG(1) printf("   PAMR_set_top_ratios <<\n");
   return;
}

void PAMR_get_top_ratios(real *top_xy, real *top_xz, real *top_yz)
{
   IFG(1) printf(">> PAMR_get_top_ratios \n");

   if (!(curr_context)) { printf("PAMR_get_top_ratios: error ... no current context\n"); return; }

   *top_xy=curr_context->top_xy;
   *top_xz=curr_context->top_xz;
   *top_yz=curr_context->top_yz;

   IFG(1) printf("   PAMR_get_top_ratios <<\n");
   return;
}

void PAMR_set_ghost_width(int *ghost_width)
{
   int i;

   IFG(1) printf(">> PAMR_set_ghost_width \n");

   if (!(curr_context)) { printf("PAMR_set_ghost_width: error ... no current context\n"); return; }
   if (curr_context->curr_cgh) { printf("PAMR_set_ghost_width: error ... cannot change after "
                                        " the initial hierarchy has been composed\n"); return; }

   for (i=0; i<curr_context->dim; i++)
      curr_context->ghost_width[i]=ghost_width[i];

   IFG(1) printf("   PAMR_set_ghost_width << \n");
   return;
}

void PAMR_get_ghost_width(int *ghost_width)
{
   int i;

   IFG(1) printf(">> PAMR_get_ghost_width \n");

   if (!(curr_context)) { printf("PAMR_get_ghost_width: error ... no current context\n"); return; }

   for (i=0; i<curr_context->dim; i++)
      ghost_width[i]=curr_context->ghost_width[i];

   IFG(1) printf("   PAMR_get_ghost_width << \n");
   return;
}

void PAMR_set_periodic_bdy(int *periodic)
{
   int i;

   IFG(1) printf(">> PAMR_set_periodic_bdy \n");

   if (!(curr_context)) { printf("PAMR_set_periodic_bdy: error ... no current context\n"); return; }
   if (curr_context->curr_cgh) { printf("PAMR_set_periodic_bdy: error ... cannot change after "
                                        " the initial hierarchy has been composed\n"); return; }

   for (i=0; i<curr_context->dim; i++)
      curr_context->periodic[i]=periodic[i];

   IFG(1) printf("   PAMR_set_periodic_bdy << \n");
   return;
}

void PAMR_get_periodic_bdy(int *periodic)
{
   int i;

   IFG(1) printf(">> PAMR_get_periodic_bdy \n");

   if (!(curr_context)) { printf("PAMR_get_periodic_bdy: error ... no current context\n"); return; }

   for (i=0; i<curr_context->dim; i++)
      periodic[i]=curr_context->periodic[i];

   IFG(1) printf("   PAMR_get_periodic_bdy << \n");
   return;
}

void PAMR_set_min_width(int *min_width)
{
   int i;

   IFG(1) printf(">> PAMR_set_min_width \n");

   if (!(curr_context)) { printf("PAMR_set_min_width: error ... no current context\n"); return; }
   if (curr_context->curr_cgh) { printf("PAMR_set_min_width: error ... cannot change after "
                                        " the initial hierarchy has been composed\n"); return; }

   for (i=0; i<curr_context->dim; i++)
      curr_context->min_width[i]=min_width[i];

   IFG(1) printf("   PAMR_set_min_width <<\n");
   return;
}

void PAMR_get_min_width(int *min_width)
{
   int i;

   IFG(1) printf(">> PAMR_get_min_width \n");

   if (!(curr_context)) { printf("PAMR_get_min_width: error ... no current context\n"); return; }

   for (i=0; i<curr_context->dim; i++)
      min_width[i]=curr_context->min_width[i];

   IFG(1) printf("   PAMR_get_min_width <<\n");
   return;
}

void PAMR_set_MG_coarse_width(int *min_width)
{
   int i;

   IFG(1) printf(">> PAMR_set_MG_coarse_width \n");

   if (!(curr_context)) { printf("PAMR_set_MG_coarse_width: error ... no current context\n"); return; }
   if (curr_context->curr_cgh) { printf("PAMR_set_MG_coarse_width: error ... cannot change after "
                                        " the initial hierarchy has been composed\n"); return; }

   for (i=0; i<curr_context->dim; i++)
   {
      if (min_width[i]<2) printf("PAMR_set_MG_coarse_width: error ... width[%i]<2\n",i);
      curr_context->MG_min_cwidth[i]=max(2,min_width[i]);
   }

   IFG(1) printf("   PAMR_set_MG_coarse_width <<\n");
   return;
}

void PAMR_get_MG_coarse_width(int *min_width)
{
   int i;

   IFG(1) printf(">> PAMR_get_MG_coarse_width \n");

   if (!(curr_context)) { printf("PAMR_get_MG_coarse_width: error ... no current context\n"); return; }

   for (i=0; i<curr_context->dim; i++)
   {
      min_width[i]=curr_context->MG_min_cwidth[i];
   }

   IFG(1) printf("   PAMR_get_MG_coarse_width <<\n");
   return;
}

void PAMR_set_gdm(int method)
{
   int i;

   IFG(1) printf(">> PAMR_set_gdm\n");

   if (!(curr_context)) { printf("PAMR_set_gdm: error ... no current context\n"); return; }
   curr_context->gdm=method;

   IFG(1) printf("   PAMR_set_gdm <<\n");
   return;
}

void PAMR_get_gdm(int *method)
{
   int i;

   IFG(1) printf(">> PAMR_get_gdm\n");

   if (!(curr_context)) { printf("PAMR_set_gdm: error ... no current context\n"); return; }
   *method=curr_context->gdm;

   IFG(1) printf("   PAMR_get_gdm <<\n");
   return;
}

void PAMR_set_time(int lev, real t)
{
   context *c=curr_context;

   IFG(1) printf(">> PAMR_set_time(lev=%i,t=%lf)\n",lev,t);

   if (!(c)) { printf("PAMR_set_time: error ... no current context\n"); return; }
   if (lev>=1 && lev<=PAMR_MAX_LEVS && c->curr_cgh->levels[lev-1]) c->curr_cgh->levels[lev-1]->t=t;
   else printf("PAMR_set_time: error ... no level %i\n",lev);

   IFG(1) printf("   PAMR_set_time <<\n");
   return;
}

real PAMR_get_time(int lev)
{
   context *c=curr_context;
   real t=0;

   IFG(1) printf(">> PAMR_get_time(lev=%i)\n",lev);

   if (!(c)) { printf("PAMR_get_time: error ... no current context\n"); return 0; }
   if (lev>=1 && lev<=PAMR_MAX_LEVS && c->curr_cgh->levels[lev-1]) t=c->curr_cgh->levels[lev-1]->t;
   else printf("PAMR_get_time: error ... no level %i\n",lev);

   IFG(1) printf("   PAMR_get_time <<\n");
   return t;
}

real PAMR_tick(int lev)
{
   context *c=curr_context;
   grid *g;
   real t=0;

   IFG(1) printf(">> PAMR_tick(lev=%i)\n",lev);

   if (!(c)) { printf("PAMR_tick: error ... no current context\n"); return 0; }
   if (lev>=1 && lev<=PAMR_MAX_LEVS && c->curr_cgh->levels[lev-1]) 
   {
      c->curr_cgh->levels[lev-1]->t+=c->dt[lev-1];
      t=c->curr_cgh->levels[lev-1]->t;
      g=c->curr_cgh->levels[lev-1]->grids;
      while(g) {g->t=t; g=g->next;}
   }
   else printf("tick: error ... no level %i\n",lev);

   IFG(1) printf("   PAMR_tick <<\n");
   return t;
}

//-----------------------------------------------------------------------------
// the following swaps time levels tl1 and tl2 at level lev in the AMR
// hierarchy
//-----------------------------------------------------------------------------
int PAMR_swap_tl(int lev, int tl1, int tl2)
{
   int valid;
   static int first=1;
   var *f;
   real *p;

   context *c=curr_context;

   IFG(1) printf(">> PAMR_swap_tl(lev=%i, tl1=%i, tl2=%i)\n",lev,tl1,tl2);

   if (!(c)) { printf("PAMR_swap_tl: error ... no current context\n"); return 0; }
   if (c->MG_cgh && first) 
   { 
      first=0;
      printf("PAMR_swap_tl: WARNING (first only) ... MG hierarchy exists, and is ignored by the swap operation\n");
   }
   if (tl1<1 || tl2<1) { printf("PAMR_swap_tl: error ... time level numbers start at 1: tl1=%i,tl2=%i\n",tl1,tl2); return 0; }

   PAMR_push_iter();
   valid=PAMR_init_s_iter(lev,PAMR_AMRH,0);
   while(valid)
   {
      f=c->vars;
      while(f)
      {
         if (f->num_tl>=tl1 && f->num_tl>=tl2)
         {
            p=c->iter_g->gfs[f->sgfn-1 + tl1-1];
            c->iter_g->gfs[f->sgfn-1 + tl1-1]=c->iter_g->gfs[f->sgfn-1 + tl2-1];
            c->iter_g->gfs[f->sgfn-1 + tl2-1]=p;
         }
         f=f->next;
      }
      valid=PAMR_next_g();
   }
   PAMR_pop_iter();

   IFG(1) printf("   PAMR_swap_tl <<\n");
   return 1;
}

//-----------------------------------------------------------------------------
// PAMR_get_max/min_lev(hier) returns the current maximum/minimum level within the 
// hierarchy "hier" 
//-----------------------------------------------------------------------------
int PAMR_get_max_lev(int hier)
{
   int l=0;
   cgh *h;

   IFG(1) printf(">> PAMR_get_max_lev \n");

   if (!(curr_context)) { printf("PAMR_get_max_lev: error ... no current context\n"); return 0; }

   if (hier==PAMR_AMRH)
   {
      if (!(h=curr_context->curr_cgh)) { printf("PAMR_get_max_lev: error ... no AMR hierarchy\n"); return 0; }
   }
   else if (hier==PAMR_MGH)
   {
      if (!(h=curr_context->MG_cgh)) { printf("PAMR_get_max_lev: error ... no MG hierarchy\n"); return 0; }
   }
   else
   {
      printf("PAMR_get_max_lev: error ... invalid hierarchy\n"); return 0; 
   }

   l=PAMR_MAX_LEVS-1;
   while ((l>=0) && !(h->levels[l])) l--;

   l++;

   IFG(1) printf("   PAMR_get_max_lev <<\n");
   return l;
}

int PAMR_get_min_lev(int hier)
{
   int l=0;
   cgh *h;

   IFG(1) printf(">> PAMR_get_min_lev \n");

   if (!(curr_context)) { printf("PAMR_get_min_lev: error ... no current context\n"); return 0; }

   if (hier==PAMR_AMRH)
   {
      if (!(h=curr_context->curr_cgh)) { printf("PAMR_get_min_lev: error ... no AMR hierarchy\n"); return 0; }
   }
   else if (hier==PAMR_MGH)
   {
      if (!(h=curr_context->MG_cgh)) { printf("PAMR_get_min_lev: error ... no MG hierarchy\n"); return 0; }
   }
   else
   {
      printf("PAMR_get_min_lev: error ... invalid hierarchy\n"); return 0; 
   }

   l=0;
   while ((l<PAMR_MAX_LEVS) && !(h->levels[l])) l++;

   l++;

   IFG(1) printf("   PAMR_get_min_lev <<\n");
   return l;
}

//-----------------------------------------------------------------------------
// returns starting grid function number (from 1 ... #gfns), or 0 if an error occurs
//-----------------------------------------------------------------------------
int PAMR_def_var_brief(char *name)
{
   var *f;
   int i,num;

   IFG(1) printf(">> PAMR_def_var_brief(name=<%s>) \n",name);

   if (!(curr_context)) { printf("PAMR_def_var_[]: error ... no current context\n"); return 0; }
   if (c_in_amrh && c_num_tl<0) { printf("PAMR_def_var_[]: error ... for AMR variable,  num_tl=%i "
                                         "must be greater than 0\n",c_num_tl); return 0; }
   if (!c_in_amrh && c_num_tl>0) { printf("PAMR_def_var_[]: error ... for MG variable,  num_tl=%i "
                                         "must be 0\n",c_num_tl); return 0; }
   if (curr_context->curr_cgh) { printf("PAMR_def_var_[]: error ... cannot define a new variable after "
                                        " the initial hierarchy has been composed\n"); return 0; }
   num=0;
   if (c_in_amrh) num+=c_num_tl;
   if (c_in_mgh) num++;
   if (curr_context->gfns+num>PAMR_MAX_GFNS) { printf("PAMR_def_var_[]: maximum number of grid functions exceeded\n"); 
                                               return 0; }

   f=new_var(name); 
   if (!f) { printf("PAMR_def_var_brief: error ... variable %s already exists\n",name); return 0; }

   f->in_amrh=c_in_amrh;
   f->in_mgh=c_in_mgh;
   f->num_tl=c_num_tl;
   f->amr_inject=c_amr_inject;
   f->amr_interp=c_amr_interp;
   f->amr_bdy_interp=c_amr_bdy_interp;
   f->amr_sync=c_amr_sync;
   f->mg_inject=c_mg_inject;
   f->mg_interp=c_mg_interp;
   f->mg_sync=c_mg_sync;
   f->mg_noinj_to_amr=c_mg_noinj_to_amr;
   f->regrid_transfer=c_regrid_transfer;
   f->c_to_v=c_c_to_v;
   f->v_to_c=c_v_to_c;
   for (i=0; i<(2*curr_context->dim); i++) 
      f->phys_bdy_type[i]=c_phys_bdy_type[i];

   f->sgfn=curr_context->gfns+1;
   f->status=VAR_STATUS_ON;
   curr_context->gfns+=num;

   IFG(1) printf("   PAMR_def_var_brief <<\n");
   return f->sgfn;
}

//-----------------------------------------------------------------------------
// returns starting grid function number, or -1 if an error occurs
//-----------------------------------------------------------------------------
int PAMR_def_var_full(char *name, int in_amrh, int in_mgh, int num_tl, int amr_inject,
                      int amr_interp, int amr_bdy_interp, int amr_sync, int mg_inject, 
                      int mg_interp, int mg_sync, int mg_noinj_to_amr, int regrid_transfer, 
		      int c_to_v, int v_to_c, int *phys_bdy_type)
{
   int i,ret;

   IFG(1) printf(">> PAMR_def_var_full(name=<%s>,...) \n",name);

   c_in_amrh=in_amrh;
   c_in_mgh=in_mgh;
   c_num_tl=num_tl;
   c_amr_inject=amr_inject;
   c_amr_interp=amr_interp;
   c_amr_bdy_interp=amr_bdy_interp;
   c_amr_sync=amr_sync;
   c_mg_inject=mg_inject;
   c_mg_interp=mg_interp;
   c_mg_sync=mg_sync;
   c_regrid_transfer=regrid_transfer;
   c_mg_noinj_to_amr=mg_noinj_to_amr;
   c_c_to_v=c_to_v;
   c_v_to_c=v_to_c;
   for (i=0; i<(2*curr_context->dim); i++) c_phys_bdy_type[i]=phys_bdy_type[i];

   ret=PAMR_def_var_brief(name);

   IFG(1) printf("   PAMR_def_var_full <<\n");
   return ret;
}

//-----------------------------------------------------------------------------
// returns starting grid function number, or 0 if variable doesn't exist
//-----------------------------------------------------------------------------
int PAMR_set_var_attribs(char *name, int in_amrh, int in_mgh, int num_tl, int amr_inject,
                   int amr_interp, int amr_bdy_interp, int amr_sync, int mg_inject, 
                   int mg_interp, int mg_sync, int mg_noinj_to_amr, int regrid_transfer, 
		   int c_to_v, int v_to_c, int *phys_bdy_type)
{
   int i;
   var *f;

   IFG(1) printf(">> PAMR_set_var_attribs(name=<%s>,...) \n",name);

   if (!(curr_context)) { printf("PAMR_set_var_attribs[]: error ... no current context\n"); return 0; }

   f=find_var(name);
   if (!f) { printf("PAMR_set_var_attribs[]: error ... variable <%s> does not exist\n",name); return 0; }

   f->in_amrh=in_amrh;
   f->in_mgh=in_mgh;
   f->num_tl=num_tl;
   f->amr_inject=amr_inject;
   f->amr_interp=amr_interp;
   f->amr_bdy_interp=amr_bdy_interp;
   f->amr_sync=amr_sync;
   f->mg_inject=mg_inject;
   f->mg_interp=mg_interp;
   f->mg_sync=mg_sync;
   f->mg_noinj_to_amr=mg_noinj_to_amr;
   f->regrid_transfer=regrid_transfer;
   f->c_to_v=c_to_v;
   f->v_to_c=v_to_c;
   for (i=0; i<(2*curr_context->dim); i++) f->phys_bdy_type[i]=phys_bdy_type[i];

   IFG(1) printf("   PAMR_set_var_attribs <<\n");
   return f->sgfn;
}
//-----------------------------------------------------------------------------
// returns starting grid function number, or 0 if variable doesn't exist
//-----------------------------------------------------------------------------
int PAMR_get_var_attribs(char *name, int *in_amrh, int *in_mgh, int *num_tl, int *amr_inject,
                   int *amr_interp, int *amr_bdy_interp, int *amr_sync, int *mg_inject, 
                   int *mg_interp, int *mg_sync, int *mg_noinj_to_amr, int *regrid_transfer, 
                   int *c_to_v, int *v_to_c, int *phys_bdy_type)
{
   int i;
   var *f;

   IFG(1) printf(">> PAMR_get_var_attribs(name=<%s>,...) \n",name);

   if (!(curr_context)) { printf("PAMR_get_var_attribs[]: error ... no current context\n"); return 0; }

   f=find_var(name);
   if (!f) { printf("PAMR_get_var_attribs[]: error ... variable <%s> does not exist\n",name); return 0; }

   *in_amrh=f->in_amrh;
   *in_mgh=f->in_mgh;
   *num_tl=f->num_tl;
   *amr_inject=f->amr_inject;
   *amr_interp=f->amr_interp;
   *amr_bdy_interp=f->amr_bdy_interp;
   *amr_sync=f->amr_sync;
   *mg_inject=f->mg_inject;
   *mg_interp=f->mg_interp;
   *mg_sync=f->mg_sync;
   *mg_noinj_to_amr=f->mg_noinj_to_amr;
   *regrid_transfer=f->regrid_transfer;
   *c_to_v=f->c_to_v;
   *v_to_c=f->v_to_c;
   for (i=0; i<(2*curr_context->dim); i++) phys_bdy_type[i]=f->phys_bdy_type[i];

   IFG(1) printf("   PAMR_get_var_attribs <<\n");
   return f->sgfn;
}

//-----------------------------------------------------------------------------
// prevents variable from being allocated upon subsequent regrids, and
// and turns off transfer bits.
//-----------------------------------------------------------------------------
void PAMR_disable_var(char *name)
{
   var *f;
   int i;

   IFG(1) printf(">> PAMR_disable_var(name=<%s>,...) \n",name);

   if (!(curr_context)) { printf("PAMR_disable_var[]: error ... no current context\n"); return; }

   f=find_var(name);
   if (!f) { printf("PAMR_disable_var[]: error ... variable <%s> does not exist\n",name); return; }

   if (f->status != VAR_STATUS_OFF) f->status=VAR_STATUS_TURN_OFF;

   if (curr_context->tf_bits) {
      if (f->in_amrh) for (i=0; i<f->num_tl; i++) curr_context->tf_bits[f->sgfn-1+i]=0;
      if (f->in_mgh) curr_context->tf_bits[f->sgfn-1+f->num_tl]=0;
   }

   IFG(1) printf("   PAMR_disable_var <<\n");
}

//-----------------------------------------------------------------------------
// allows a variable to be allocated upon subsequent regrids
//-----------------------------------------------------------------------------
void PAMR_enable_var(char *name)
{
   var *f;
   int i;

   IFG(1) printf(">> PAMR_enable_var(name=<%s>,...) \n",name);

   if (!(curr_context)) { printf("PAMR_enable_var[]: error ... no current context\n"); return; }

   f=find_var(name);
   if (!f) { printf("PAMR_enable_var[]: error ... variable <%s> does not exist\n",name); return; }

   if (f->status != VAR_STATUS_ON) f->status=VAR_STATUS_TURN_ON;

   IFG(1) printf("   PAMR_enable_var <<\n");
}

//-----------------------------------------------------------------------------
// returns the type of a variable 
//-----------------------------------------------------------------------------
int PAMR_var_type(char *name)
{
   int clen=0,vlen=0,slen;

   IFG(1) printf(">> PAMR_var_type(name=<%s>,...) \n",name);

   if (!PAMR_c_tag && !PAMR_v_tag) return PAMR_VERTEX_CENTERED;

   if (PAMR_c_tag) clen=strlen(PAMR_c_tag);
   if (PAMR_v_tag) vlen=strlen(PAMR_v_tag);
   slen=strlen(name);
   if (PAMR_c_tag && slen > clen && !strcmp(PAMR_c_tag,&name[slen-clen])) return PAMR_CELL_CENTERED;
   if (PAMR_v_tag && slen > vlen && !strcmp(PAMR_v_tag,&name[slen-vlen])) return PAMR_VERTEX_CENTERED;
   if (PAMR_v_tag && !PAMR_c_tag) return PAMR_CELL_CENTERED;
   if (PAMR_c_tag && !PAMR_v_tag) return PAMR_VERTEX_CENTERED;
   return PAMR_VERTEX_CENTERED;

   IFG(1) printf("   PAMR_var_type <<\n");
}

//-----------------------------------------------------------------------------
// returns the variable type (CC or VC) of a gfn
//-----------------------------------------------------------------------------
int PAMR_gfn_var_type(int gfn)
{
   IFG(1) printf(">> PAMR_gfn_var_type(gfn=%i)\n",gfn);

   if (!(curr_context)) { printf("PAMR_gfn_var_type: error ... no current context\n"); return -1; }
   if (!(curr_context->gfn_var_type)) { printf("PAMR_gfn_var_type: error ... no hierarchy yet exists\n"); return -1; }
   if (gfn>=curr_context->gfns) { printf("PAMR_gfn_var_type: error ... invalid gfn=%i\n",gfn); return -1; }

   IFG(1) printf("   PAMR_gfn_var_type <<\n");

   return curr_context->gfn_var_type[gfn-1];
}

//-----------------------------------------------------------------------------
// returns 1 if a variable is enabled, 0 otherwise
//-----------------------------------------------------------------------------
int PAMR_is_var_enabled(char *name)
{
   var *f;
   int ret=0;

   IFG(1) printf(">> PAMR_is_var_enabled(name=<%s>,...) \n",name);

   if (!(curr_context)) { printf("PAMR_is_var_enabled[]: error ... no current context\n"); return 0; }

   f=find_var(name);
   if (!f) { printf("PAMR_is_var_enabled[]: error ... variable <%s> does not exist\n",name); return 0; }

   if (f->status == VAR_STATUS_ON) ret=1;

   IFG(1) printf("   PAMR_is_var_enabled <<\n");
   return ret;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void PAMR_set_interp_buffer(int interp_buffer)
{
   IFG(1) printf(">> PAMR_set_interp_buffer\n");

   if (!(curr_context)) { printf("PAMR_set_interp_buffer: error ... no current context\n"); return; }
   curr_context->interp_buffer=interp_buffer;

   IFG(1) printf("   PAMR_set_interp_buffer <<\n");
   return;
}

void PAMR_get_interp_buffer(int *interp_buffer)
{
   IFG(1) printf(">> PAMR_get_interp_buffer \n");

   if (!(curr_context)) { printf("PAMR_get_interp_buffer: error ... no current context\n"); return; }
   *interp_buffer=curr_context->interp_buffer;

   IFG(1) printf("   PAMR_get_interp_buffer <<\n");
   return;
}

//=============================================================================
// Accessing grids/grid-functions in the hierarchy
// All functions, except PAMR_get_g_x() and PAMR_get_g_gfs() expect the user to
// supply the storage for the returned quantity. For PAMR_get_g_x() and
// PAMR_get_g_gfs(), the user supplies an array of pointers to the relevant
// objects, and the routines fill in the pointers.
//=============================================================================
int PAMR_init_s_iter(int l, int hier, int all)
{
   context *c=curr_context;
   int ret=1;

   IFG(1) printf(">> PAMR_init_s_iter \n");
   if (!(c)) { printf("PAMR_init_s_iter: error ... no current context\n"); return 0; }
   if (l<1 || l>PAMR_MAX_LEVS) { printf("PAMR_init_s_iter: error ... invalid level %i\n",l); return 0; }

   c->iter_g=0;
   c->iter_all=all;
   c->iter_L=l;

   switch(hier)
   {
      case PAMR_MGH: if (!(c->MG_cgh)) 
                { printf("PAMR_init_s_iter: error ... no MG hierarchy\n"); return 0; }
                if (!(c->MG_cgh->levels[l-1])) 
                { printf("PAMR_init_s_iter: error ... no level %i in MG hierarchy\n",l); return 0; }
                c->iter_g=c->MG_cgh->levels[l-1]->grids;
                break;
      case PAMR_AMRH: if (!(c->curr_cgh)) 
                { printf("PAMR_init_s_iter: error ... no AMR hierarchy\n"); return 0; }
                if (!(c->curr_cgh->levels[l-1])) 
                { printf("PAMR_init_s_iter: error ... no level %i in AMR hierarchy\n",l); return 0; }
                c->iter_g=c->curr_cgh->levels[l-1]->grids;
                break;
      default: printf("PAMR_init_s_iter: error ... unknown hierarchy %i\n",hier); return 0;
   }

   if (!all) while(c->iter_g && c->iter_g->rank != c->rank) c->iter_g=c->iter_g->next;

   if (c->iter_g) ret=1; else ret=0;

   IFG(1) printf("   PAMR_init_s_iter << \n");
   return ret;
}

int PAMR_next_g()
{
   context *c=curr_context;
   int ret;

   IFG(1) printf(">> PAMR_next_g \n");
   if (!(c)) { printf("PAMR_next_g: error ... no current context\n"); return 0; }

   if (!c->iter_g) {ret=0; goto fin;}

   c->iter_g=c->iter_g->next;
   if (!c->iter_all) while(c->iter_g && c->iter_g->rank != c->rank) c->iter_g=c->iter_g->next;

   if (c->iter_g) ret=1; else ret=0;

fin: 
   IFG(1) printf("   PAMR_next_g << \n");
   return ret;
}

int PAMR_push_iter()
{
   context *c=curr_context;
   int ret;

   IFG(1) printf(">> PAMR_push_iter \n");
   if (!(c)) { printf("PAMR_push_iter: error ... no current context\n"); return 0; }

   IFG(1) printf("   PAMR_push_iter before: g,L,all,stack_top: %p,%i,%i,%i\n",c->iter_g,c->iter_L,c->iter_all,c->iter_g_stack_top);

   if (c->iter_g_stack_top==ITER_G_STACKSIZE) { printf("PAMR_push_iter: error ... stack is full\n"); return 0; }
   c->iter_g_stack[c->iter_g_stack_top]=c->iter_g;
   c->iter_L_stack[c->iter_g_stack_top]=c->iter_L;
   c->iter_all_stack[c->iter_g_stack_top++]=c->iter_all;

   IFG(1) printf("   PAMR_push_iter after: g,L,all,stack_top: %p,%i,%i,%i\n",c->iter_g,c->iter_L,c->iter_all,c->iter_g_stack_top);

   IFG(1) printf("   PAMR_push_iter << \n");
   return 1;
}

int PAMR_pop_iter()
{
   context *c=curr_context;
   int ret;

   IFG(1) printf(">> PAMR_pop_iter \n");
   if (!(c)) { printf("PAMR_pop_iter: error ... no current context\n"); return 0; }

   if (c->iter_g_stack_top==0) { printf("PAMR_pop_iter: error ... stack is empty\n"); return 0; }

   IFG(1) printf("   PAMR_pop_iter before: g,L,all,stack_top: %p,%i,%i,%i\n",c->iter_g,c->iter_L,c->iter_all,c->iter_g_stack_top);

   c->iter_g=c->iter_g_stack[--c->iter_g_stack_top];
   c->iter_L=c->iter_L_stack[c->iter_g_stack_top]; 
   c->iter_all=c->iter_all_stack[c->iter_g_stack_top];

   IFG(1) printf("   PAMR_pop_iter after: g,L,all,stack_top: %p,%i,%i,%i\n",c->iter_g,c->iter_L,c->iter_all,c->iter_g_stack_top);

   if (c->iter_g) ret=1; else ret=0;

   IFG(1) printf("   PAMR_push_iter << \n");
   return ret;
}

int PAMR_get_g_rank(int *rank)
{
   context *c=curr_context;
   grid *g;

   IFG(1) printf(">> PAMR_get_g_rank \n");

   if (!(c)) { printf("PAMR_get_g_rank: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_rank: error ... iterator not initialized\n"); return 0; }

   *rank=g->rank;

   IFG(1) printf("   PAMR_get_g_rank << \n");
   return 1;
}

int PAMR_get_g_level(int *L)
{
   context *c=curr_context;
   grid *g;

   IFG(1) printf(">> PAMR_get_g_level \n");

   if (!(c)) { printf("PAMR_get_g_level: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_level: error ... iterator not initialized\n"); return 0; }

   *L=c->iter_L;

   IFG(1) printf("   PAMR_get_g_level << \n");
   return 1;
}

int PAMR_get_g_dim(int *dim)
{
   context *c=curr_context;
   grid *g;

   IFG(1) printf(">> PAMR_get_g_dim \n");

   if (!(c)) { printf("PAMR_get_g_dim: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_dim: error ... iterator not initialized\n"); return 0; }

   *dim=g->dim;

   IFG(1) printf("   PAMR_get_g_dim << \n");
   return 1;
}

int PAMR_get_g_shape(int *shape)
{
   context *c=curr_context;
   grid *g;
   int i;

   IFG(1) printf(">> PAMR_get_g_shape \n");

   if (!(c)) { printf("PAMR_get_g_shape: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_shape: error ... iterator not initialized\n"); return 0; }

   for (i=0; i<g->dim; i++) shape[i]=g->shape[i];

   IFG(1) printf("   PAMR_get_g_shape << \n");
   return 1;
}

int PAMR_get_g_shape_c(int *shape_c)
{
   context *c=curr_context;
   grid *g;
   int i;

   IFG(1) printf(">> PAMR_get_g_shape_c \n");

   if (!(c)) { printf("PAMR_get_g_shape_c: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_shape_c: error ... iterator not initialized\n"); return 0; }

   for (i=0; i<g->dim; i++) shape_c[i]=g->shape_c[i];

   IFG(1) printf("   PAMR_get_g_shape_c << \n");
   return 1;
}

int PAMR_get_g_bbox(real *bbox)
{
   context *c=curr_context;
   grid *g;
   int i;

   IFG(1) printf(">> PAMR_get_g_bbox \n");

   if (!(c)) { printf("PAMR_get_g_bbox: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_bbox: error ... iterator not initialized\n"); return 0; }

   for (i=0; i<2*g->dim; i++) bbox[i]=g->bbox[i];

   IFG(1) printf("   PAMR_get_g_bbox << \n");
   return 1;
}

void PAMR_get_global_bbox(real *bbox)
{
   context *c=curr_context;
   int i;

   IFG(1) printf(">> PAMR_global_bbox \n");

   if (!(c)) { printf("PAMR_global_bbox: error ... no current context\n"); return; }

   for (i=0; i<2*c->dim; i++) bbox[i]=c->bbox[i];

   IFG(1) printf("   PAMR_get_global_bbox << \n");
   return;
}


int PAMR_get_g_ghost_width(int *ghost_width)
{
   context *c=curr_context;
   grid *g;
   int i;

   IFG(1) printf(">> PAMR_get_g_ghost_width \n");

   if (!(c)) { printf("PAMR_get_g_ghost_width: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_ghost_width: error ... iterator not initialized\n"); return 0; }

   for (i=0; i<2*g->dim; i++) ghost_width[i]=g->ghost_width[i];

   IFG(1) printf("   PAMR_get_g_ghost_width << \n");
   return 1;
}

int PAMR_get_g_t(real *t)
{
   context *c=curr_context;
   grid *g;

   IFG(1) printf(">> PAMR_get_g_t \n");

   if (!(c)) { printf("PAMR_get_g_t: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_t: error ... iterator not initialized\n"); return 0; }

   *t=g->t;

   IFG(1) printf("   PAMR_get_g_t << \n");
   return 1;
}

int PAMR_get_g_coarsest(int *coarsest)
{
   context *c=curr_context;
   grid *g;

   IFG(1) printf(">> PAMR_get_g_coarsest \n");

   if (!(c)) { printf("PAMR_get_g_coarsest: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_coarsest: error ... iterator not initialized\n"); return 0; }

   *coarsest=g->coarsest;

   IFG(1) printf("   PAMR_get_g_coarsest << \n");
   return 1;
}

int PAMR_set_g_comm(int comm)
{
   context *c=curr_context;
   grid *g;

   IFG(1) printf(">> PAMR_set_g_comm \n");

   if (!(c)) { printf("PAMR_set_g_comm: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_set_g_comm: error ... iterator not initialized\n"); return 0; }

   g->comm=comm;

   IFG(1) printf("   PAMR_set_g_comm << \n");
   return 1;
}

int PAMR_get_g_comm(int *comm)
{
   context *c=curr_context;
   grid *g;

   IFG(1) printf(">> PAMR_get_g_comm \n");

   if (!(c)) { printf("PAMR_get_g_comm: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_comm: error ... iterator not initialized\n"); return 0; }

   *comm=g->comm;

   IFG(1) printf("   PAMR_set_g_comm << \n");
   return 1;
}

int PAMR_get_g_ngfs(int *ngfs)
{
   context *c=curr_context;
   grid *g;

   IFG(1) printf(">> PAMR_get_g_ngfs \n");

   if (!(c)) { printf("PAMR_get_g_ngfs: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_ngfs: error ... iterator not initialized\n"); return 0; }

   *ngfs=g->ngfs;

   IFG(1) printf("   PAMR_get_g_ngfs << \n");
   return 1;
}

int PAMR_get_g_x(real **x)
{
   context *c=curr_context;
   grid *g;
   int i;

   IFG(1) printf(">> PAMR_get_g_x \n");

   if (!(c)) { printf("PAMR_get_g_x: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_x: error ... iterator not initialized\n"); return 0; }

   for (i=0; i<g->dim; i++) x[i]=g->x[i];

   IFG(1) printf("   PAMR_get_g_x << \n");
   return 1;
}

int PAMR_get_g_x_c(real **x_c)
{
   context *c=curr_context;
   grid *g;
   int i;

   IFG(1) printf(">> PAMR_get_g_x_c \n");

   if (!(c)) { printf("PAMR_get_g_x_c: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_x_c: error ... iterator not initialized\n"); return 0; }

   for (i=0; i<g->dim; i++) x_c[i]=g->x_c[i];

   IFG(1) printf("   PAMR_get_g_x_c << \n");
   return 1;
}

int PAMR_get_g_gfs(real **gfs)
{
   context *c=curr_context;
   grid *g;
   int i;

   IFG(1) printf(">> PAMR_get_g_gfs \n");

   if (!(c)) { printf("PAMR_get_g_gfs: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_gfs: error ... iterator not initialized\n"); return 0; }

   if (g->gfs) for (i=0; i<g->ngfs; i++) gfs[i]=g->gfs[i];
   else for (i=0; i<g->ngfs; i++) gfs[i]=0; // not a local grid

   IFG(1) printf("   PAMR_get_g_gfs << \n");
   return 1;
}


int PAMR_get_g_attribs(int *rank, int *dim, int *shape, int *shape_c, real *bbox,
                   int *ghost_width, real *t, int *ngfs, real **x, real **x_c, real **gfs)
{
   context *c=curr_context;
   grid *g;
   int ret=1;

   IFG(1) printf(">> PAMR_get_g_attribs \n");

   if (!(c)) { printf("PAMR_get_g_attribs: error ... no current context\n"); return 0; }
   if (!(g=c->iter_g)) { printf("PAMR_get_g_attribs: error ... iterator not initialized\n"); return 0; }

   ret*=PAMR_get_g_rank(rank);
   ret*=PAMR_get_g_dim(dim);
   ret*=PAMR_get_g_shape(shape);
   ret*=PAMR_get_g_shape_c(shape_c);
   ret*=PAMR_get_g_bbox(bbox);
   ret*=PAMR_get_g_ghost_width(ghost_width);
   ret*=PAMR_get_g_t(t);
   ret*=PAMR_get_g_ngfs(ngfs);
   ret*=PAMR_get_g_x(x);
   ret*=PAMR_get_g_x_c(x_c);
   ret*=PAMR_get_g_gfs(gfs);

   IFG(1) printf("   PAMR_get_g_attribs << \n");
   return ret;
}

//=============================================================================
// returns the requested grid-function number, or -1 if failure
//=============================================================================
int PAMR_get_gfn(char *name, int hier, int tl)
{
   context *c=curr_context;
   int offset;
   var *f;

   IFG(1) printf(">> PAMR_get_gfn \n");

   if (!(c)) { printf("PAMR_get_gfn: error ... no current context\n"); return 0; }

   f=find_var(name);
   if (!f) { printf("PAMR_get_gfn: error ... variable <%s> does not exist\n",name); return 0; }

   if (hier==PAMR_MGH && (!f->in_mgh)) 
      { printf("PAMR_get_gfn: error ... variable <%s> is not in the MG hierarchy\n",name); return 0; }
   else if (hier==PAMR_MGH) offset=f->num_tl;
   else if (hier==PAMR_AMRH && (!f->in_amrh)) 
      { printf("PAMR_get_gfn: error ... variable <%s> is not in the AMR hierarchy\n",name); return 0; }
   else if (hier!=PAMR_AMRH)
      { printf("PAMR_get_gfn: error ... unknown hierarhcy %i\n",hier); return 0; }
   else if (tl<1 || tl>f->num_tl)
      { printf("PAMR_get_gfn: error ... invalid time level %i\n",tl); return 0; }
   else offset=tl-1;

   IFG(1) printf("   PAMR_get_gfn <<\n");

   return f->sgfn+offset;
}

//=============================================================================
// transfer bit manipulation functions. 
//=============================================================================
void PAMR_freeze_tf_bits(void)
{
   IFG(1) printf(">> PAMR_freeze_tf_bits \n");

   if (!(curr_context)) { printf("PAMR_freeze_tf_bits: error ... no current context\n"); return; }
  
   curr_context->frozen=1;

   IFG(1) printf("   PAMR_freeze_tf_bits <<\n");
}

void PAMR_thaw_tf_bits(void)
{
   IFG(1) printf(">> PAMR_thaw_tf_bits \n");

   if (!(curr_context)) { printf("PAMR_thaw_tf_bits: error ... no current context\n"); return; }
   
   curr_context->frozen=0;

   IFG(1) printf("   PAMR_thaw_tf_bits <<\n");
}

void PAMR_clear_tf_bits(void)
{
   int i;

   IFG(1) printf(">> PAMR_clear_tf_bits \n");

   if (!(curr_context)) { printf("PAMR_clear_tf_bits: error ... no current context\n"); return; }
   if (!(curr_context->tf_bits)) { printf("PAMR_clear_tf_bits: error ... no hierarchy\n"); return; }

   for (i=0; i<curr_context->gfns; i++) curr_context->tf_bits[i]=0;

   IFG(1) printf("   PAMR_clear_tf_bits <<\n");
}

//=============================================================================
// NOTE: val should be set to one of the transfer types, e.g.
// STRAIGHT_INJECTION, SECOND_ORDER, NO_SYNC, etc.
// Therefore, val is operation dependent! 
//
// This operation also freezes tf bits (only 'mode' in which individual tf
// functions make sense)
//
// WARNING ... this function does *NOT* check wether corresponding gf's are
// allocated!
//=============================================================================
void PAMR_set_tf_bit(int gf, int val)
{
   IFG(1) printf(">> PAMR_set_tf_bit \n");

   if (!(curr_context)) { printf("PAMR_set_tf_bit: error ... no current context\n"); return; }
   if (!(curr_context->tf_bits)) { printf("PAMR_set_tf_bits: error ... no hierarchy\n"); return; }

   if (gf<1 || gf>curr_context->gfns) { printf("PAMR_set_tf_bits: error ... gf=%i out of range\n",gf); return; }

   PAMR_freeze_tf_bits();

   curr_context->tf_bits[gf-1]=val;

   IFG(1) printf("   PAMR_set_tf_bit <<\n");
}

//=============================================================================
// tl is the time-level ... PAMR_set to -1 for all (ingored for PAMR_MGH)
// hierarchy is one of PAMR_AMRH or PAMR_MGH
// stage is one of PAMR_TF_SYNC,PAMR_TF_COMPOSE(AMRH ONLY),PAMR_TF_INJECT,
//                 PAMR_TF_INJECT_TO_MG_LEV(PAMR_MGH ONLY),PAMR_TF_INTERPOLATE
//                 PAMR_TF_V_TO_C, or PAMR_TF_C_TO_V
//
//
// use PAMR_TF_INJECT_TO_MG_LEV instead of PAMR_TF_INJECT when injecting to a 
// 'pure' MG level (i.e., it is either coarser than the coarsest
// AMR level, or sits 'between' AMR levels).
// use PAMR_TF_MGH_INIT to only inject functions that exist in both the
// MG and AMR hieararchies.
//
// This function only works when transfer bits are thawed
//=============================================================================
void PAMR_set_tf_bits(int tl, int hierarchy, int stage)
{
   int gfs,gfe,i,op;
   context *c=curr_context;
   var *f;

   IFG(1) printf(">> PAMR_set_tf_bits\n");

   f=c->vars;

   if (!(c)) { printf("PAMR_set_tf_bits: error ... no current context\n"); return; }
   if (!(c->tf_bits)) { printf("PAMR_set_tf_bits: error ... no hierarchy\n"); return; }

   if (c->frozen) return;

   PAMR_clear_tf_bits();

   while(f)
   {

      gfs=1;gfe=0;
      op=0;
      if (f->status!=VAR_STATUS_ON) op=0;
     
      else if (hierarchy==PAMR_AMRH && f->in_amrh)
      {
         if ((stage==PAMR_TF_SYNC && (op=f->amr_sync)) ||        
             (stage==PAMR_TF_COMPOSE && (op=f->regrid_transfer)) ||        
             (stage==PAMR_TF_INJECT && (op=f->amr_inject)) ||        
             (stage==PAMR_TF_INTERP && (op=f->amr_interp)) ||
             (stage==PAMR_TF_BDY_INTERP && (op=f->amr_bdy_interp)) ||
             (stage==PAMR_TF_BDY_INTERP_C && (op=f->amr_bdy_interp)) ||
             (stage==PAMR_TF_C_TO_V && (op=f->c_to_v)) ||
             (stage==PAMR_TF_V_TO_C && (op=f->v_to_c)))
         {
            if (tl>0 && tl<=f->num_tl)
            {
               gfs=f->sgfn+(tl-1);
               gfe=gfs+1;
            }
            else if (tl<0)
            {
               gfs=f->sgfn;
               gfe=gfs+f->num_tl;
            }
         }
                         
      }
      else if (hierarchy==PAMR_MGH && f->in_mgh)
      {
         if ((stage==PAMR_TF_SYNC && (op=f->mg_sync)) ||        
             (stage==PAMR_TF_INJECT && (op=f->mg_inject) && !f->mg_noinj_to_amr) ||        
             (stage==PAMR_TF_INJECT_TO_MG_LEV && (op=f->mg_inject)) ||        
             (stage==PAMR_TF_INTERP && (op=f->mg_interp)) ||
             (stage==PAMR_TF_C_TO_V && (op=f->c_to_v)) ||
             (stage==PAMR_TF_V_TO_C && (op=f->v_to_c)) ||
             (stage==PAMR_TF_MGH_INIT && f->in_amrh && (op=f->mg_inject)))
         {
            gfs=f->sgfn+f->num_tl;
            gfe=gfs+1;
         }
      }

      for (i=gfs; i<gfe; i++) {
	if (hierarchy==PAMR_AMRH && stage==PAMR_TF_BDY_INTERP) {
	  if (c->gfn_var_type[i]==PAMR_VERTEX_CENTERED) {
	    c->tf_bits[i-1]=op;
	  } else {
	    c->tf_bits[i-1]=0;
	  }
	} else if (hierarchy==PAMR_AMRH && stage==PAMR_TF_BDY_INTERP_C) {
	  if (c->gfn_var_type[i]==PAMR_CELL_CENTERED) {
	    c->tf_bits[i-1]=op;
	  } else {
	    c->tf_bits[i-1]=0;
	  }
	} else {
	  c->tf_bits[i-1]=op;
	}
        IFG(3) printf("PAMR_set_tf_bits (stage=%i): '%s', sgfn=%i, gfn(%i) op=%i\n",stage,f->name,f->sgfn,i,op);
      }
      f=f->next;
   }

   IFG(1) printf("   PAMR_set_tf_bits<<\n");
}

int PAMR_do_save_gfn(char *name, int hier, int tl, int L, real t, char *pre_tag, char *post_tag,int close)
{
   context *c=curr_context;
   int gfn,l,ls,le;
   level *lev;
   grid *g;
   cgh *gh;
   char out_name[256];
   real t0,bbox[2*PAMR_MAX_DIM];

   IFG(1) printf(">> PAMR_save_gfn(%s,...)\n",name);

   if (!(c)) { printf("PAMR_save_gfn: error ... no current context\n"); return 0; }

   if (!(gfn=PAMR_get_gfn(name,hier,tl))) return 0;

   if (hier==PAMR_AMRH) 
   {
      sprintf(out_name,"%s%s_tl%i%s_%i",pre_tag,name,tl,post_tag,c->rank);
      gh=c->curr_cgh;
      if (!gh) { printf("PAMR_save_gfn: AMR hierarchy does not exist\n"); return 0; }
   }
   else 
   {
      sprintf(out_name,"%s%s_MG%s_%i",pre_tag,name,post_tag,c->rank);
      gh=c->MG_cgh;
      if (!gh) { printf("PAMR_save_gfn: MG hierarchy does not exist\n"); return 0; }
   }

   if (L<0) { ls=gh->min_lev-1; le=gh->max_lev-1; } else { ls=le=L-1; }

   for (l=ls; l<=le; l++)
   {
      g=gh->levels[l]->grids;
      while(g)
      {
         if (t==-1.0) t0=gh->levels[l]->t; else t0=t;
         if (g->rank==c->rank && g->gfs && g->gfs[gfn-1])
         {
            if (c->gfn_var_type[gfn-1]==PAMR_VERTEX_CENTERED)
               gft_out_bbox(out_name,t0,g->shape,g->dim,g->bbox,g->gfs[gfn-1]);
            else
            {
               bbox[0]=(g->x_c[0])[0]; bbox[1]=(g->x_c[0])[g->shape_c[0]-1];
               if (g->dim>1) { bbox[2]=(g->x_c[1])[0]; bbox[3]=(g->x_c[1])[g->shape_c[1]-1]; }
               if (g->dim>2) { bbox[4]=(g->x_c[2])[0]; bbox[5]=(g->x_c[2])[g->shape_c[2]-1]; }
               gft_out_bbox(out_name,t0,g->shape_c,g->dim,bbox,g->gfs[gfn-1]);
            }
         }
         g=g->next;
      }
   }

   if (close) gft_close(out_name);

   IFG(1) printf("   PAMR_save_gfn(%s,gfn=%i) <<\n",out_name,gfn);

   return 1;
}

//=============================================================================
// saves specified grid function to
//
// <pre_tag><name>[_tl#|_MG]<post_tag>_<rank>.sdf
//
// with time 't' ... if t==-1.0 then uses current time
//
// L is the level, -1 for all
//
// the _close version below is identical to PAMR_save_gfn(), except
// the sdf is closed afterwards (so will get over-written, instead
// of appended with subsequent calls)
//=============================================================================
int PAMR_save_gfn(char *name, int hier, int tl, int L, real t, char *pre_tag, char *post_tag)
{
   return PAMR_do_save_gfn(name,hier,tl,L,t,pre_tag,post_tag,0);
}

int PAMR_save_gfn_close(char *name, int hier, int tl, int L, real t, char *pre_tag, char *post_tag)
{
   return PAMR_do_save_gfn(name,hier,tl,L,t,pre_tag,post_tag,1);
}

//=============================================================================
// excision support functions.
//=============================================================================
int PAMR_excision_on(char *ex_mask_var, char *ex_mask_c_var,
                     void (*app_fill_ex_mask_fnc)(real *mask, real *mask_c, int dim,
                                                  int *shape, int *shape_c, real *bbox, real excised),
                     real excised, int initialize_now) 
{
   context *c=curr_context;
   var *f;

   IFG(1) printf(">> PAMR_excision_on(...)\n");

   if (!(c)) { printf("PAMR_excision_on: error ... no current context\n"); return 0; }

   c->amr_mask_gfn=c->mg_mask_gfn=c->amr_mask_c_gfn=c->mg_mask_c_gfn=0;

   if (ex_mask_var)
   {
      if (PAMR_var_type(ex_mask_var)!=PAMR_VERTEX_CENTERED) 
         { printf("PAMR_excision_on: error ... variable %s is not vertex centered\n",ex_mask_var); return 0; }

      if (!(f=find_var(ex_mask_var))) 
         { printf("PAMR_excision_on: error ... no variable %s\n",ex_mask_var); return 0; }

      if (f->status!=VAR_STATUS_ON)
         { printf("PAMR_excision_on: error ... variable %s is not active \n",ex_mask_var); return 0; }

      if (!f->in_amrh || !f->in_mgh)
         { printf("PAMR_excision_on: error ... variable %s must be defined in both the MG"
                  " and AMR hierarchies\n",ex_mask_var); return 0; }

      if (f->num_tl !=1)
         { printf("PAMR_excision_on: error ... variable %s can only have 1 time level\n",ex_mask_var); return 0; }

      c->amr_mask_gfn=f->sgfn;
      c->mg_mask_gfn=f->sgfn+1;
   }

   if (ex_mask_c_var)
   {
      if (PAMR_var_type(ex_mask_c_var)!=PAMR_CELL_CENTERED) 
         { printf("PAMR_excision_on: warning ... variable %s is not cell centered\n",ex_mask_c_var); }

      if (!(f=find_var(ex_mask_c_var))) 
         { printf("PAMR_excision_on: error ... no variable %s\n",ex_mask_c_var); return 0; }

      if (f->status!=VAR_STATUS_ON)
         { printf("PAMR_excision_on: error ... variable %s is not active \n",ex_mask_c_var); return 0; }

      if (!f->in_amrh || !f->in_mgh)
         { printf("PAMR_excision_on: error ... variable %s must be defined in both the MG"
                  " and AMR hierarchies\n",ex_mask_c_var); return 0; }

      if (f->num_tl !=1)
         { printf("PAMR_excision_on: error ... variable %s can only have 1 time level\n",ex_mask_c_var); return 0; }

      c->amr_mask_c_gfn=f->sgfn;
      c->mg_mask_c_gfn=f->sgfn+1;
   }

   c->excision_on=1;
   c->excised=excised;
   c->app_fill_ex_mask=app_fill_ex_mask_fnc;

   if (initialize_now && c->curr_cgh) init_ex_mask(c->curr_cgh->min_lev,c->curr_cgh->max_lev,
                                                   c->curr_cgh,c->amr_mask_gfn,c->amr_mask_c_gfn);
   if (initialize_now && c->MG_cgh) init_ex_mask(c->MG_cgh->min_lev,c->MG_cgh->max_lev,
                                                 c->MG_cgh,c->mg_mask_gfn,c->mg_mask_c_gfn);

   IFG(1) printf("   PAMR_excision_on() <<\n");

   return 1;
}


void PAMR_excision_off()
{
   context *c=curr_context;

   IFG(1) printf(">> PAMR_excision_off()\n");

   if (!(c)) { printf("PAMR_excision_off: error ... no current context\n"); return; }

   c->excision_on=0;

   IFG(1) printf("   PAMR_excision_off() <<\n");

   return;
}

//=================================================================================
// check-point routine
//=================================================================================
int PAMR_cp(char *cp_file)
{
   context *c=curr_context;
   int ret,cnode,cnode_max,cp_err,g_cp_err;
   int ltrace=1;

   IFG(1) printf(">> PAMR_cp(%s)\n",cp_file);

   if (!cp_file || (int)strlen(cp_file)==0) { printf("PAMR_cp: error ... no file specified\n"); return 0; }
   if (!(c)) { printf("PAMR_cp: error ... no current context\n"); return 0; }
   if (c->MG_cgh) { printf("PAMR_cp: error ... cannot currently save a MG hierarchy\n"); return 0; }

   cp_err=g_cp_err=0;
   if (PAMR_parallel_io) cnode_max=1; else cnode_max=c->size;
   for (cnode=0; cnode<cnode_max; cnode++)
   {
      if (PAMR_parallel_io || c->rank==cnode)
      {
         ret=PAMR_do_cp(cp_file,c->rank,c->size,PAMR_CP_DIR_SAVE);
         if (!ret)
         {
            printf("PAMR_init_context: error ... unable to restore state from file %s\n",cp_file);
            g_cp_err=cp_err=1;
         }
      }
      MPI_Allreduce(&cp_err,&g_cp_err,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
      if (!PAMR_parallel_io && ltrace && c->rank==0) printf(" ... PAMR_cp finished on node %i\n",cnode);
   }

   if (g_cp_err) ret=0; else ret=1;

   IFG(1) printf("   PAMR_cp <<\n");
   return ret;
}

//=================================================================================
// stats ... if turned on, all stats are cleared
//=================================================================================
void PAMR_collect_stats(int on_off)
{
   IFG(1) printf(">> PAMR_collect_stats(%i)\n",on_off);

   if (on_off)
   {
      PAMR_comm_r_secs=PAMR_comm_r_microsecs=PAMR_comm_r_GB=PAMR_comm_r_B=0;
      PAMR_comm_s_secs=PAMR_comm_s_microsecs=PAMR_comm_s_GB=PAMR_comm_s_B=0;
      PAMR_num_s=PAMR_num_r=0;
      PAMR_stats_on=1;
   }
   else PAMR_stats_on=0;

   IFG(1) printf("<< PAMR_collect_stats(%i)\n",on_off);
   return;
}

void PAMR_get_stats(real *comm_r_microsecs, real *comm_s_microsecs,
                    real *comm_r_bytes, real *comm_s_bytes,
		    int *num_r, int *num_s)
{
   IFG(1) printf(">> PAMR_get_stats\n");

   *comm_r_microsecs=PAMR_comm_r_secs/1.0e6+PAMR_comm_r_microsecs;
   *comm_s_microsecs=PAMR_comm_s_secs/1.0e6+PAMR_comm_s_microsecs;
   *comm_r_bytes=PAMR_comm_r_GB+PAMR_comm_r_B/1.0e9;
   *comm_s_bytes=PAMR_comm_s_GB+PAMR_comm_s_B/1.0e9;
   *num_r=PAMR_num_r;
   *num_s=PAMR_num_s;

   IFG(1) printf("<< PAMR_get_stats\n");
   return;
}

