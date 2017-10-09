//========================================================================================
// io.c 
//
// sdf input/output routines
//========================================================================================
#include "gh_w.h"
#include "misc_w.h"
#include "pamr_w.h"
#include "transfer_w.h"
#include "io_w.h"
#include <bbhutil.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//=================================================================================
// utility routines used by PAMR_cp
//
// the copy_block() routine : 
//
// if (*q && *p && dir==PAMR_CP_DIR_SAVE) 
//   copy n bytes from *p to *q,
//   increments *q, and returns the number of byte added to q
// else if (*q && dir==PAMR_CP_DIR_RESTORE)
//   same, but direction of copy is reversed. Also, if *p is null,
//   it is allocated before the copy
// else just returns the number of bytes.
//
// copy_type_block() below are a set of wrappers for different data types
// (returned size is always in bytes though)
//=================================================================================
int copy_block(char **p, char **q, int n, int dir)
{
   char *p0;
   int n0,n0_s,ltrace=0;

   if (n==0) return 0;

   IFL printf(" copy_block ... n=%i\n",n);

   p0=(*p);
   n0_s=n0=n;

   if (p0 && *q && dir==PAMR_CP_DIR_SAVE) while(n0--) *(*q)++=*p0++;
   else if (*q && dir==PAMR_CP_DIR_RESTORE) 
   {
      if (!p0)
      {
         if (!(p0=malloc(n0_s))) { printf("copy_int_block ... error ... out of memory\n"); exit(1); }
         *p=p0;
      }
      while(n0--) *p0++=*(*q)++;
   }

   return n0_s;
}

int copy_int_block(int **p, char **q, int n, int dir)
{
   return copy_block((char **)p,q,n*sizeof(int),dir);
}

int copy_real_block(real **p, char **q, int n, int dir)
{
   return copy_block((char **)p,q,n*sizeof(real),dir);
}

//=================================================================================
// the following function searches for the shape,bbox pair within the current
// context's set of grids ... as a sanity check, makes sure rank's match
// 
// for V2 ... matches shape to either g->shape or g->shape_c
//
// min_dx/10 is used in fuzz compares
//=================================================================================
#define IFL00 if (ltrace && my_rank==0)
grid *cp_find_g(int *shape, real *bbox, int rank, real min_dx)
{
   int i,j,jj,found;
   level *lev;
   grid *g;
   int ltrace=0;
   static int first=1;

   if (first)
   {
      first=0; IFL0 printf("(cp_find_g: using min_dx=%e)\n",min_dx);
   }

   IFL0
   {
       printf("cp_find_g ... rank=%i ... looking for grid [",rank);
       for (j=0; j<curr_context->dim; j++) printf("%i,",shape[j]);
       printf("] [");
       for (j=0; j<curr_context->dim; j++) printf("%16.10lf,%16.10lf,",bbox[2*j],bbox[2*j+1]);
       printf("]\n");
   }

   for (i=0; i<PAMR_MAX_LEVS; i++)
   {
      if (lev=curr_context->curr_cgh->levels[i])
      {
         g=lev->grids;
         while(g)
         {
            found=1;
            for (j=0; j<g->dim; j++)
            {
               if (g->shape[j]!=shape[j] && g->shape_c[j]!=shape[j]) found=0;
               if (!(fuzz_eq(g->bbox[2*j],bbox[2*j],min_dx/10))) found=0;
               if (!(fuzz_eq(g->bbox[2*j+1],bbox[2*j+1],min_dx/10))) found=0;
            }
            IFL0
            {
                printf("   cp_find_g ... rank=%i ... comparing to grid [",rank);
                for (jj=0; jj<g->dim; jj++) printf("%i,",g->shape[jj]);
                printf("] [");
                for (jj=0; jj<g->dim; jj++) printf("%16.10lf,%16.10lf,",g->bbox[2*jj],g->bbox[2*jj+1]);
                printf("]\n");
            }
            if (found) { if (g->rank!=rank) printf("cp_find_g: WARNING ... rank's don't match!\n"); return g; }
            g=g->next;
         }
      }
   }

   printf("cp_find_g: error ... grid not found in hierarchy\n");
   return 0;
}

//=================================================================================
// check-point routine
//
// Saves most internal data structures and variables to disk. Variables that
// aren't saved, such as pointers, etc., are those that are run-time dependent
// and will be reallocated upon a restore. NOTE however in some cases
// pointers are saved, but interpreted as booleans upon restore to indicate
// whether a structure should be allocated or not.
//
// NOTES: 1. currently only works with a single context
//        2. must do a recompose afterwards, as grids are *not* added
//           in the proper "insert_grid()" order.
//=================================================================================
int PAMR_do_cp(char *cp_file, int my_rank, int mpi_size, int dir)   
{
   int size,size0,i,ii,pass,pass0,*ip,level,j,jj,k,restore_rank,ok,m,err;
   int ind;
   real version=PAMR_CP_VERSION,*rp,bbox[2*PAMR_MAX_DIM];
   double *perim_coords;
   grid *g,*pg;
   int max_levs=PAMR_MAX_LEVS,ret=0,shape[PAMR_MAX_DIM];
   var *vp,*pvp;
   char *data,*name,*data0,file_name[256];
   context *c=curr_context;
   char cnames[256];
   real dtime,dx,x,min_dx;
   int ltrace=0,sizec;
   char sys_command[256];
   static int first=1;
   int sys_trace=0;
   int *itp;
   real *rtp;
   grid **gtp;
   real **rptp,***rpptp;

   if (first && sys_trace) { sprintf(sys_command,"/bin/hostname > .node.%i",my_rank); system(sys_command); }
   first=0;

   IFL00 printf("PAMR_do_cp: rank=%i, mpi_size=%i, dir=%i, file=%s\n",my_rank,mpi_size,dir,cp_file);
   if (sys_trace) printf("PAMR_do_cp: rank=%i, mpi_size=%i, dir=%i, file=%s\n",my_rank,mpi_size,dir,cp_file);

   size=0;
   data=0;
   data0=0;
   level=2;
   if (dir==PAMR_CP_DIR_RESTORE && c) { printf("PAMR_do_cp: error ... context already exists upon attempted restore\n"); goto fin; }
   if (dir==PAMR_CP_DIR_SAVE && (c->iter_g || c->iter_g_stack_top))
      printf("PAMR_do_cp: WARNING ... cannot checkpoint within an iterator (continuing though)\n"
             "node %i may be in such a state: c->iter_g=%p, c->iter_g_stack_top=%i\n",
             my_rank,c->iter_g,c->iter_g_stack_top); 

   if (dir==PAMR_CP_DIR_RESTORE) 
   {
      pass0=2; 
      sprintf(file_name,"%s_%i",cp_file,my_rank);
      if (!(gft_read_shape(file_name,1,&size0)))
      {
         sprintf(file_name,"%s_0",cp_file);
         if (!(gft_read_shape(file_name,1,&size0))) { printf("PAMR_do_cp: error opening file %s\n",file_name); goto fin;}
      }
      if (!(data0=data=malloc(size0*sizeof(double)))) { printf("PAMR_do_cp: error ... out of memory\n"); goto fin;}
      if (!(gft_read_brief(file_name,1,(double *)data))) { printf("PAMR_do_cp: error reading file %s\n",file_name); goto fin;}

      IFL00 printf("my_rank=%i, read %i bytes from %s\n",my_rank,size0*8,file_name);
   }
   else pass0=0;

   //==============================================================================================
   // for saving, a 3 pass system:
   //
   // pass 1 : calculate size of grid hierarchy information
   // pass 2 : fill "data" with grid hierarchy, and on all nodes save it to level 1 of cp_n.sdf
   // pass 3 : save all local grid (but no perimeter arrays) information to levels 2,... of cp_n.sdf
   //
   // restore in 3 stages, but only use a single pass of above:
   //
   // stage 1 (above): all nodes read level 1 of either cp_n.sdf, or cp_0.sdf. The first option
   // allows for restoration on nodes which only have access to a local file system,
   // but the restriction here is that restoration must occur to the same MPI size computation.
   // (you can get around this by manually copying files to 'simulate' a global fs)
   // The second option is for a global file system, and the restore can go to any number
   // of nodes.
   //
   // stage 2 (pass 3 below): restore the hierarchy ... rank -> rank % N, where
   // there are N nodes in the current run.
   //
   //
   // stage 3 : node n reads the data, in
   // levels 2,... from *all* of the following: cp_n.sdf, cp_[n+N-1].sdf, cp_[n+2*N-1].sdf, ...
   //==============================================================================================

   for (pass=pass0; pass<3; pass++)
   {
      IFL00 printf("   rank=%i: pass=%i\n",my_rank,pass);
      rp=&version; size+=copy_real_block(&rp,&data,1,dir);
      if (version != PAMR_CP_VERSION && version != PAMR_CP_VERSION_P1)   
      {
         printf("PAMR_do_cp: error ... unknown version number %lf\n",version);
         goto fin;
      }
      IFL00 printf("   version=%lf\n:",version);

      //-----------------------------------------------------------------------
      // PAMR_CP_VERSION_P1 did not have the last two mask gfn's, nor pointer
      // to var types
      //-----------------------------------------------------------------------
      if (version > PAMR_CP_VERSION_P1) 
      {
         size+=copy_block((char **)&c,&data,sizeof(struct context),dir);
      }
      else
      {
         size+=copy_block((char **)&c,&data,sizeof(struct context)-3*sizeof(int),dir);
         if (dir==PAMR_CP_DIR_RESTORE) { c->mg_mask_c_gfn=c->amr_mask_c_gfn=0; }
      }
      if (dir==PAMR_CP_DIR_RESTORE)
      {
         curr_context=c;
         c->vars=0;
         c->tf_bits=0;
         c->curr_cgh=0;
         c->MG_cgh=0;
         c->iter_g=0;
         c->iter_g_stack_top=0;
         if (c->excision_on && my_rank==0) printf("\nPAMR_cp: WARNING ... excision turned off on restore ... calling program must re-declare the fill_ex_mask function\n\n");
         c->excision_on=0;
         c->app_fill_ex_mask=0;
         c->gfn_var_type=0;
      }
      if (my_rank==0) printf("ghost_width=[%i,%i,%i]\n",c->ghost_width[0],c->ghost_width[1],c->ghost_width[2]);


      ip=&max_levs; size+=copy_int_block(&ip,&data,1,dir);
      if (max_levs!=PAMR_MAX_LEVS) 
      {
         printf("PAMR_do_cp: error ... PAMR_MAX_LEVS=%i disagrees with value(%i) in cp file\n",PAMR_MAX_LEVS,max_levs);
         goto fin;
      }

      for (i=0; i<max_levs; i++)
      {
         if (dir==PAMR_CP_DIR_RESTORE) c->sgh_bboxes[i]=0;
         size+=copy_real_block(&c->sgh_bboxes[i],&data,c->num_sgh_bboxes[i]*2*c->dim,dir);
      }

      for (i=0, vp=c->vars, pvp=0; i<c->num_vars; i++)
      {
         //-----------------------------------------------------------------------
         // PAMR_CP_VERSION_P1 did not have the c_to_v,v_to_c flags
         //-----------------------------------------------------------------------
         if (version > PAMR_CP_VERSION_P1) 
         {
            size+=copy_block((char **)&vp,&data,sizeof(struct var),dir);
         }
         else
         {
            size+=copy_block((char **)&vp,&data,sizeof(struct var)-2*sizeof(int),dir);
            if (dir==PAMR_CP_DIR_RESTORE) { vp->c_to_v=vp->v_to_c=0; }
         }
         if (dir==PAMR_CP_DIR_RESTORE) { name=data; vp->name=0;}  else name=vp->name;
         size+=copy_block(&(vp->name),&data,(int)strlen(name)+1,dir);
         if (dir==PAMR_CP_DIR_RESTORE)
         {
            if (pvp) pvp->next=vp; else c->vars=vp;
            vp->next=0;
            pvp=vp;
         }
         vp=vp->next;
      }

      size+=copy_int_block(&c->tf_bits,&data,c->gfns,dir);
      if (version > PAMR_CP_VERSION_P1) size+=copy_int_block(&c->gfn_var_type,&data,c->gfns,dir);
      else if (dir==PAMR_CP_DIR_RESTORE) c->gfn_var_type=0;
      size+=copy_block((char **)&c->curr_cgh,&data,sizeof(struct cgh),dir);

      for (i=0; i<PAMR_MAX_LEVS; i++)
      {
         if (c->curr_cgh->levels[i])
         {
            min_dx=c->dx[i][0]; // for fuzzy comparison in find_g
            if (dir==PAMR_CP_DIR_RESTORE) c->curr_cgh->levels[i]=0;
            size+=copy_block((char **)&c->curr_cgh->levels[i],&data,sizeof(struct level),dir);
            g=c->curr_cgh->levels[i]->grids;
            IFL00 printf("   my_rank=%i ... level %i, first g=%p\n",my_rank,i+1,g);
            pg=0;
            while(g)
            {
               IFL00 printf("     my_rank=%i ... level %i, curr g=%p\n",my_rank,i+1,g);
               if (dir==PAMR_CP_DIR_RESTORE)
               {
                  if (!(g=(struct grid *)malloc(sizeof(struct grid))))
                     { printf("PAMR_do_cp: error ... out of memory \n"); goto fin; }
               }
               //--------------------------------------------------------------
               // need to split up cp of grid structures:
               //
               // size+=copy_block((char **)&g,&data,sizeof(struct grid),dir);
               //
               // due to added info in V2
               //--------------------------------------------------------------
               gtp=&(g->next); size+=copy_block((char **)&(gtp),&data,sizeof(struct grid*),dir);
               gtp=&(g->prev); size+=copy_block((char **)&(gtp),&data,sizeof(struct grid*),dir);
               itp=&(g->dim); size+=copy_block((char **)&(itp),&data,sizeof(int),dir);
               itp=&(g->shape[0]); size+=copy_block((char **)&(itp),&data,sizeof(int)*PAMR_MAX_DIM,dir);
               if (version > PAMR_CP_VERSION_P1) 
               {
                  itp=&(g->shape_c[0]); size+=copy_block((char **)&(itp),&data,sizeof(int)*PAMR_MAX_DIM,dir);
               }
               else
               {
                  for (j=0; j<PAMR_MAX_DIM; j++) g->shape_c[j]=0;
               }
               rtp=&(g->bbox[0]); size+=copy_block((char **)&(rtp),&data,sizeof(real)*2*PAMR_MAX_DIM,dir);
               itp=&(g->ghost_width[0]); size+=copy_block((char **)&(itp),&data,sizeof(int)*2*PAMR_MAX_DIM,dir);
               itp=&(g->wrap[0]); size+=copy_block((char **)&(itp),&data,sizeof(int)*2*PAMR_MAX_DIM,dir);
               rtp=&(g->t); size+=copy_block((char **)&(rtp),&data,sizeof(real),dir);
               rptp=&(g->x[0]); size+=copy_block((char **)&(rptp),&data,sizeof(real *)*PAMR_MAX_DIM,dir);
               if (version > PAMR_CP_VERSION_P1) 
               {
                  rptp=&(g->x_c[0]); size+=copy_block((char **)&(rptp),&data,sizeof(real *)*PAMR_MAX_DIM,dir);
               }
               else
               {
                  for (j=0; j<PAMR_MAX_DIM; j++) g->x_c[j]=0;
               }
               itp=&(g->ngfs); size+=copy_block((char **)&(itp),&data,sizeof(int),dir);
               rpptp=&(g->gfs); size+=copy_block((char **)&(rpptp),&data,sizeof(real *),dir);
               itp=&(g->rank); size+=copy_block((char **)&(itp),&data,sizeof(int),dir);
               itp=&(g->coarsest); size+=copy_block((char **)&(itp),&data,sizeof(int),dir);
               itp=&(g->comm); size+=copy_block((char **)&(itp),&data,sizeof(int),dir);
               itp=&(g->virtual_g); size+=copy_block((char **)&(itp),&data,sizeof(int),dir);

               if (dir==PAMR_CP_DIR_RESTORE) 
               {
                  if (pg) pg->next=g; else c->curr_cgh->levels[i]->grids=g;
                  g->prev=pg;
                  g->rank=g->rank%mpi_size;
                  for (j=0; j<PAMR_MAX_DIM; j++) g->x[j]=g->x_c[j]=0;
                  if (g->rank==my_rank)
                  {
                     if (!(g->gfs=(real **)malloc(sizeof(real *)*g->ngfs))) 
                        { printf("PAMR_do_cp: error ... out of memory \n"); goto fin; }
                     for (j=0; j<g->ngfs; j++) g->gfs[j]=0;
                  }
                  else g->gfs=0;
               }
               else if (pass==2 && dir==PAMR_CP_DIR_SAVE && g->gfs && !g->virtual_g)
               {
		  for (j=0; j<g->ngfs; j++)
                  {
                     IFL00 printf("   my_rank=%i, g=%p, j=%i, f=%p, shape=[%i,%i]\n",
                                my_rank,g,j,g->gfs[j],g->shape[0],g->shape[1]);
                     if (g->gfs[j]) 
                     {
                        if (c->gfn_var_type[j]==PAMR_VERTEX_CENTERED)
                        {
                          if (!(gft_out_bbox(file_name,0.0e0,g->shape,g->dim,g->bbox,g->gfs[j]))) 
                           { printf("PAMR_do_cp: error saving grid information to cp_file=%s\n",file_name); goto fin; }
                        }
                        else
                        {
                          if (!(gft_out_bbox(file_name,0.0e0,g->shape_c,g->dim,g->bbox,g->gfs[j])))  // important! ... same bbox as VC
                           { printf("PAMR_do_cp: error saving grid information to cp_file=%s\n",file_name); goto fin; }
                        }
                     }
                     level++;
                  }
                  IFL00 printf("   my_rank=%i ... done with g=%p ... g->next=%p\n",my_rank,g,g->next);
               }
               pg=g; g=g->next;
            }
            IFL00 printf("   my_rank=%i ... done with level %i\n",my_rank,i+1);
         }
      }

      // round up to double
      size+=(sizeof(double)-size%(sizeof(double )));

      if (pass==0 && dir==PAMR_CP_DIR_SAVE)
         if (!(data0=malloc(size))) { printf("PAMR_do_cp: error ... out of memory\n"); goto fin; }
      data=data0;
      if (pass==1 && dir==PAMR_CP_DIR_SAVE)
      {
         size/=sizeof(double);
         sprintf(file_name,"%s_%i",cp_file,my_rank);
         if (!(gft_out_brief(file_name,0.0,&size,1,(double *) data))) { printf("PAMR_do_cp: error saving file %s\n",cp_file); goto fin;}
      }
      size=0;
   }

   free(data0);
   data0=0;

   // when reading in grids, we assume that, starting at sdf level 2,
   // all the grid functions are stored in the sdf the same order
   // that they are required in the grid structure

   if (dir==PAMR_CP_DIR_RESTORE)
   {
      restore_rank=my_rank;
      sprintf(file_name,"%s_%i",cp_file,restore_rank);
      level=2;
      ok=gft_read_shape(file_name,level,shape);
      while(ok)
      {
         IFL00 printf("my_rank=%i, restore_rank=%i, attempting to restore from %s\n",my_rank,restore_rank,file_name);
         while(ok)
         {
            for (i=0; i<c->gfns; i++)
            {
               if (is_in_amrh(i+1))
               {
		 // Note that size and sizec will now depend on the level inside the checkpoint file.
		 // That's because cell- and vertex- centered gfs can be stored in successive levels.
		 gft_read_shape(file_name,level,shape);
		 size=shape[0];
		 for (ii=1; ii<c->dim; ii++) size*=shape[ii];
		 sizec=shape[0]; for (ii=1; ii<c->dim; ii++) sizec+=shape[ii];

		 if (!(perim_coords=(double *)malloc(sizeof(double)*sizec))) { printf("PAMR_do_cp: error ... out of memory\n"); goto fin; }

                  IFL00 printf("   my_rank = %i ... reading gf %i\n",my_rank,i+1);
                  if (!(data=pmalloc(sizeof(double)*size))) { printf("PAMR_do_cp: error ... out of memory\n"); goto fin; }

                  if (!(err=gft_read_full(file_name,level++,shape,cnames,c->dim,&dtime,perim_coords,(double *)data))) 
                     { printf("PAMR_do_cp: error reading level %i data from %s ... err=%i\n",level-1,file_name,err); goto fin;}

                  bbox[0]=perim_coords[0]; bbox[1]=perim_coords[shape[0]-1];

                  if (c->dim>1) bbox[2]=perim_coords[shape[0]]; bbox[3]=perim_coords[shape[0]+shape[1]-1];

                  if (c->dim>2) { bbox[4]=perim_coords[shape[0]+shape[1]]; bbox[5]=perim_coords[shape[0]+shape[1]+shape[2]-1]; }
                  if (c->dim>3) { printf("PAMR_do_cp: error ... dim>3 not yet supported\n"); goto fin;}

                  if (!(g=cp_find_g(shape,bbox,my_rank,min_dx))) { printf("PAMR_do_cp: error in check point file\n"); goto fin; }
                  if (!(g->gfs)) { printf("PAMR_do_cp: error in check point file (grid is non-local)\n"); goto fin; }
                  g->gfs[i]=(real *)data;

                  for (j=0;j<2*c->dim;j++) g->bbox[j]=bbox[j];
                  if (!g->x[0])
                  {
                     for (j=0;j<c->dim; j++)
                     {
                        if (!(g->x[j]=(real *)pmalloc(sizeof(real)*g->shape[j])))
                           {printf("PAMR_do_cp: error ... out of memory\n"); goto fin; }
                        dx=(g->bbox[2*j+1]-g->bbox[2*j])/(g->shape[j]-1);
                        for (m=0,x=g->bbox[2*j]; m<g->shape[j]; m++,x+=dx) (g->x[j])[m]=x;
                     }
                  }
                  if (!g->x_c[0])
                  {
                     for (j=0;j<c->dim; j++)
                     {
                        if (!(g->x_c[j]=(real *)pmalloc(sizeof(real)*g->shape_c[j])))
                           {printf("PAMR_do_cp: error ... out of memory\n"); goto fin; }
                        dx=(g->bbox[2*j+1]-g->bbox[2*j])/(g->shape[j]-1);
                        for (m=0,x=g->bbox[2*j]+dx/2; m<g->shape_c[j]; m++,x+=dx) (g->x_c[j])[m]=x;
                     }
                  }
                  free(perim_coords);
               }
            }
            ok=gft_read_shape(file_name,level,shape);
         }
         restore_rank+=mpi_size;
         sprintf(file_name,"%s_%i",cp_file,restore_rank);
         level=2;
         ok=gft_read_shape(file_name,level,shape);
      }
   }

   // this is important for the integrity of a saved file if the job terminates abnormally
   if (dir==PAMR_CP_DIR_SAVE) gft_close(file_name); 

   IFL00 printf("PAMR_do_cp finished ... my_rank=%i\n",my_rank);

   c->rank=my_rank;
   c->size=mpi_size;

   ret=1;

fin:
   if (data0) free(data0);
   return ret;
}

//========================================================================================
// the following routine produces a 'rank' grid function of the cgh. 
// (only rank 0 node)
//========================================================================================
#define DEBUG_CLOSE 0
void debug_save_cgh(char *name, cgh *gh)
{
   context *c=curr_context;
   level *lev;
   grid *g;
   int l,i;

   real *data;

   if (c->rank!=0 || gh==0) return;

   IFG(2) printf("   >> debug_save_cgh(name=%s,...)\n",name);

   for (l=gh->min_lev-1; l<=(gh->max_lev-1); l++)
   {
      if (lev=gh->levels[l])
      {
         g=lev->grids;
         while(g)
         {
            if (!(data=(real *)malloc(sizeof(real)*sizeof_data(g))))
            {
               printf("debug_save_cgh: error ... out of memory\n");
               return;
            }
            for (i=0;i<sizeof_data(g);i++) data[i]=g->rank;
            gft_out_bbox(name,lev->t,g->shape,g->dim,g->bbox,data);
            free(data);
            g=g->next;
         }
      }
   }

   if (DEBUG_CLOSE) gft_close(name);

   IFG(2) printf("      debug_save_cgh <<\n");
}

//========================================================================================
// the following routine produces a 'rank' grid function of the given gsl. 
// (<0 for all)
//========================================================================================
void debug_save_gsl(char *name, real t, gsl *g, int rank)
{
   int l,i,n;

   real *data,dx[2*PAMR_MAX_DIM];
   int ibbox[2*PAMR_MAX_DIM],ibbox_c[2*PAMR_MAX_DIM],shape[PAMR_MAX_DIM];

   IFG(2) printf("   >> debug_save_gsl(name=%s,gsl=%p...)\n",name,g);

   while(g)
   {
      if (rank<0 || (g->g->rank==rank))
      {
         calc_gs_ibbox(g,ibbox,ibbox_c,dx);
         n=1;
         for (i=0; i<g->g->dim; i++)
         {
            shape[i]=(ibbox[2*i+1]-ibbox[2*i]+1);
            n*=shape[i];
         }
         if (!(data=(real *)malloc(sizeof(real)*n)))
         {
            printf("debug_save_gsl: error ... out of memory\n");
            return;
         }
         for (i=0;i<n;i++) data[i]=g->g->rank; 
         gft_out_bbox(name,t,shape,g->g->dim,g->bbox,data);
         free(data);
      }
      g=g->next;
   }

   if (DEBUG_CLOSE) gft_close(name);

   IFG(2) printf("      debug_save_gsl <<\n");
}
