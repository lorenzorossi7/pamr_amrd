//          on a level) overwhich to zero the TRE. May not work properly
// cls.c --- clustering routines
// (two stage --- first find clusters, then adujst for MG)
//=============================================================================

#include <stdio.h>
#include "cls_w.h"
#include "globals_w.h"
#include "mg_w.h"
#include "evolve_w.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>

//=============================================================================
// calculates the region of high TRE (l2 norm of the set of f_tre variables
// ... this routine does *NOT* compute f_tre) then sets the tre 
// variable to 1 where TRE>TRE_max, 0 elsewhere.
//
// options: TRE_norm ... normalizes each tre by the variables' norm 
//          (over the level), if non-zero
//        : TRE_buffer ... number of points to expand tre by (calculated
//          point-by-point, so this can be used in leau of a smoothing option
//        : TRE_ibc_buffer & TRE_ibc_a_buffer ... interior-boundary-clear buffer, i.e. the region 
//          adjacent to interior boundaries (computed over the union of grids
//          on a level) overwhich to zero the TRE. May not work properly
//          if TRE_ibc_buffer > minimum size of a grid. TRE_ibc_buffer is applied
//          *BEFORE* TRE_buffer, and TRE_ibc_a_buffer is applied *AFTER*
//        : TRE_exc_buffer ... region adjacent to excised cells overwhich to clear 
//          (NOTE: only clears if size of adjacent high TRE region is less then 
//                 TRE_buffer)
//          TRE (only for levels > TRE_exc_buffer_lmin)
//        : TRE_sgpbh ... single grid per black hole 
//        : TRE_ibcp_buffer ... same as TRE_ibc_buffer, but computed
//           grid-by-grid (and independent of periodic or not)
//          NOTE: TRE_ibc_buffer/TRE_ibcp_buffer/ are applied *BEFORE* TRE_buffer.
//
// NOTE: we also enforce the min_width here (by setting the local TRE_buffer),
//       and require that TRE_buffer>=rho_sp. This is for two reasons ... first,
//       doing it here we can more easily ensure that clusters don't outgrow 
//       there bounds, and second, by requiring it to be >=rho_sp we can make
//       sure that the the action of injecting the TRE from a child to the parent
//       should always work in keeping a proper hierarchy.
//=============================================================================
void calc_tre(int L1, int L2)
{
   int L,valid,i,j,k,n,ltrace=0,exc_buf0;
   real tre0;
   static real lcount=-1;
   real tre_norm[AMRD_MAX_VARS];
   int rho_sp[PAMR_MAX_LEVS],rho_tm[PAMR_MAX_LEVS],l_TRE_buffer;

   // lcount++; use current time

   IFL printf(">> calc_tre(%i,%i)\n",L1,L2);

   PAMR_get_rho(rho_sp,rho_tm,PAMR_MAX_LEVS);

   for (L=L1; L<=L2; L++)
   {
      IFLR for (i=0; i<AMRD_num_f_tre_vars; i++) PAMR_save_gfn(AMRD_f_tre_vars[i],PAMR_AMRH,1,L,lcount,"debug_","_bare");

      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      while(valid)
      {
         ev_ldptr();
         if (AMRD_tre_gfn) for (i=0; i<AMRD_g_size; i++) AMRD_tre[i]=0;
	 if (AMRD_tre_c_gfn) for (i=0; i<AMRD_g_size_c; i++) AMRD_tre_c[i]=0;
         valid=PAMR_next_g();
      }

      //-----------------------------------------------------------------------
      // zero tre variables within TRE_ibc_buffer of interior boundaries
      //-----------------------------------------------------------------------
      if (AMRD_TRE_ibc_buffer>0)
      {
         set_cmask_bdy(L,PAMR_AMRH);
         if (AMRD_using_cc_tre) set_gfn_sync(AMRD_cmask_c_gfn); else set_gfn_sync(AMRD_cmask_gfn);
         PAMR_sync(L,0,PAMR_AMRH,0); // thawing at end
         //--------------------------------------------------------------------
         // in case there are overlapping (sequential) grids
         //--------------------------------------------------------------------
         IFLR PAMR_save_gfn("cmask",PAMR_AMRH,1,L,lcount,"debug_","_cmask_before_extend");
         if (AMRD_TRE_ibc_buffer>0 && !AMRD_using_cc_tre) extend_cmask_bdy(L,PAMR_AMRH,AMRD_TRE_ibc_buffer); 

         IFLR for (i=0; i<AMRD_num_f_tre_vars; i++) PAMR_save_gfn(AMRD_f_tre_vars[i],PAMR_AMRH,1,L,lcount,"debug_","_after_cmask");
         valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
         while(valid)
         {
            ev_ldptr();
            if (!AMRD_using_cc_tre) zero_f_tre_ibc(); 
	    if (AMRD_using_cc_tre) zero_f_tre_ibc_c();
            valid=PAMR_next_g();
         }
      }

      IFLR for (i=0; i<AMRD_num_f_tre_vars; i++) PAMR_save_gfn(AMRD_f_tre_vars[i],PAMR_AMRH,1,L,lcount,"debug_","_after_ibc");
      IFLR PAMR_save_gfn("cmask",PAMR_AMRH,1,L,lcount,"debug_","_cmask_after_ibc");
   }

   if (!(AMRD_num_f_tre_vars))
   {
      printf("calc_tre: warning ... no variables set to calculate tre\n");
      return;
   }

   //--------------------------------------------------------------------------
   // compute TRE, and propagate from finest to coarsest, to 
   // account for the (unusual) situation where a TRE might be small on
   // a part of a coarse grid, but not on a fine grid
   //--------------------------------------------------------------------------

   for (L=L2; L>=L1; L--)
   {
      if (L>AMRD_TRE_exc_buffer_lmin && AMRD_do_ex) exc_buf0=AMRD_TRE_exc_buffer; else exc_buf0=0;

      if (L>1) l_TRE_buffer=max(AMRD_TRE_buffer,rho_sp[L-2]); else l_TRE_buffer=AMRD_TRE_buffer;
      if (L<L2) 
      {
         // NOTE: this is the flagged TRE we're injecting
         if (!AMRD_using_cc_tre) set_gfn_in(AMRD_tre_gfn,PAMR_STRAIGHT_INJECT);
	 if (AMRD_using_cc_tre) set_gfn_in(AMRD_tre_c_gfn,PAMR_NN_AVERAGE);
         PAMR_inject(L+1,0,PAMR_AMRH);
	 PAMR_thaw_tf_bits();
      }

      if (AMRD_TRE_norm)
      {
         for (j=0; j<AMRD_num_f_tre_vars; j++) 
         {
            tre_norm[j]=approx_l2norm(AMRD_AMR_f_gfn[0][AMRD_f_tre_fn[j]],L,PAMR_AMRH);
            if (tre_norm[j]==0) tre_norm[j]=1;
         }
      }
      else
      {
         for (j=0; j<AMRD_num_f_tre_vars; j++) tre_norm[j]=1;
      }
      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      //-------------------------------------------------------------------------
      // to distribute l_TRE_buffer and TRE_exc_buffer accross processor
      // boundaries, we "bleed" the tre to 0 in two stages, repeating the
      // algorithm twice, first for the excision buffer then the tre buffer
      //-------------------------------------------------------------------------
      while(valid)
      {
         ev_ldptr();

	 if (!AMRD_using_cc_tre) 
	 {
	   for (i=0; i<AMRD_g_size; i++) 
	   {
	     tre0=0;
	     for (j=0; j<AMRD_num_f_tre_vars; j++) tre0+=(AMRD_f_tre[j])[i]*(AMRD_f_tre[j])[i]/tre_norm[j]/tre_norm[j];
	     tre0=sqrt(tre0/AMRD_num_f_tre_vars);
	     
	     if (tre0>AMRD_TRE_max || AMRD_tre[i]>0.5) AMRD_tre[i]=1+l_TRE_buffer; else AMRD_tre[i]=0;
	     
	     if (exc_buf0 && AMRD_chr[i]==AMRD_ex) AMRD_tre[i]=-exc_buf0;
	   }

	   if (exc_buf0) bleed_tre(exc_buf0,-1); // exc only
	   else bleed_tre(l_TRE_buffer,1);
	 }

	 if (AMRD_using_cc_tre) 
	 {

	   for (i=0; i<AMRD_g_size_c; i++) 
	   {
	     tre0=0;
	     for (j=0; j<AMRD_num_f_tre_vars; j++) tre0+=(AMRD_f_tre[j])[i]*(AMRD_f_tre[j])[i]/tre_norm[j]/tre_norm[j];
	     tre0=sqrt(tre0/AMRD_num_f_tre_vars);
	     
	     if (tre0>AMRD_TRE_max || AMRD_tre_c[i]>0.0625) AMRD_tre_c[i]=1+l_TRE_buffer; else AMRD_tre_c[i]=0;
	     
	     if (exc_buf0 && AMRD_chr_c[i]==AMRD_ex) AMRD_tre_c[i]=-exc_buf0;
	   }

	   if (exc_buf0) bleed_tre_c(exc_buf0,-1); // exc only
	   else bleed_tre_c(l_TRE_buffer,1);
	 }

         valid=PAMR_next_g();
      }
      if (!AMRD_using_cc_tre) set_gfn_sync(AMRD_tre_gfn); else set_gfn_sync(AMRD_tre_c_gfn);
      PAMR_sync(L,0,PAMR_AMRH,0); // thawing at end

      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      while(valid)
      {
         ev_ldptr();
         if (exc_buf0 && !AMRD_using_cc_tre) bleed_tre(exc_buf0,-1); 
         if (exc_buf0 && AMRD_using_cc_tre) bleed_tre_c(exc_buf0,-1); 
         if (!AMRD_using_cc_tre) bleed_tre(l_TRE_buffer,1);
         if (AMRD_using_cc_tre) bleed_tre_c(l_TRE_buffer,1);
         if (!exc_buf0 && !AMRD_using_cc_tre) {
	   for (i=0; i<AMRD_g_size; i++) if (AMRD_tre[i]>0.5) AMRD_tre[i]=1; else AMRD_tre[i]=0;
	 }
         if (!exc_buf0 && AMRD_using_cc_tre) {
	   for (i=0; i<AMRD_g_size_c; i++) if (AMRD_tre_c[i]>0.0625) AMRD_tre_c[i]=1; else AMRD_tre_c[i]=0;
	 }
         valid=PAMR_next_g();
      }
      if (!AMRD_using_cc_tre) set_gfn_sync(AMRD_tre_gfn); else set_gfn_sync(AMRD_tre_c_gfn);
      PAMR_sync(L,0,PAMR_AMRH,0); // thawing at end

      if (exc_buf0)
      {
         valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
         while(valid)
         {
            ev_ldptr();
            if (!AMRD_using_cc_tre) {
	      bleed_tre(l_TRE_buffer,1);
	      for (i=0; i<AMRD_g_size; i++) if (AMRD_tre[i]>0.5) AMRD_tre[i]=1; else AMRD_tre[i]=0;
	    }
            if (AMRD_using_cc_tre) {
	      bleed_tre_c(l_TRE_buffer,1);
	      for (i=0; i<AMRD_g_size_c; i++) if (AMRD_tre_c[i]>0.5) AMRD_tre_c[i]=1; else AMRD_tre_c[i]=0;
	    }
            valid=PAMR_next_g();
         }
         if (!AMRD_using_cc_tre) set_gfn_sync(AMRD_tre_gfn); else set_gfn_sync(AMRD_tre_c_gfn);
         PAMR_sync(L,0,PAMR_AMRH,0); // thawing at end
      }

      //-----------------------------------------------------------------------
      // zero tre variables within TRE_ibc_a_buffer of interior boundaries
      //-----------------------------------------------------------------------
      if (AMRD_TRE_ibc_a_buffer>0)
      {
         set_cmask_bdy(L,PAMR_AMRH);
         if (AMRD_using_cc_tre) set_gfn_sync(AMRD_cmask_c_gfn); else set_gfn_sync(AMRD_cmask_gfn);
         PAMR_sync(L,0,PAMR_AMRH,0); // thawing at end
         //--------------------------------------------------------------------
         // in case there are overlapping (sequential) grids
         //--------------------------------------------------------------------
         if (AMRD_TRE_ibc_a_buffer>0 && !AMRD_using_cc_tre) extend_cmask_bdy(L,PAMR_AMRH,AMRD_TRE_ibc_a_buffer); 

         valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
         while(valid)
         {
            ev_ldptr();
            if (!AMRD_using_cc_tre) zero_tre_ibc_a(); 
	    if (AMRD_using_cc_tre) zero_tre_ibc_a_c();
            valid=PAMR_next_g();
         }
         if (!AMRD_using_cc_tre) set_gfn_sync(AMRD_tre_gfn); else set_gfn_sync(AMRD_tre_c_gfn);
         PAMR_sync(L,0,PAMR_AMRH,0); // thawing at end
      }

      if (!AMRD_using_cc_tre) IFLR PAMR_save_gfn("AMRD_tre",PAMR_AMRH,1,L,lcount,"debug_","_flagged");
      if (AMRD_using_cc_tre) IFLR PAMR_save_gfn("AMRD_tre_c",PAMR_AMRH,1,L,lcount,"debug_","_flagged");
   }

   PAMR_thaw_tf_bits();

   IFL printf("<< calc_tre\n");
}

void bleed_tre(int max_buf, int min_max)
{
   int i,j,k,n,i0,j0,k0,ind;
   real adj_tre_max,adj_tre_min;

   for (n=0; n<max_buf; n++)
   {
      for (i=0; i<AMRD_g_Nx; i++)
      {
         for (j=0; j<AMRD_g_Ny; j++)
         {
            for (k=0; k<AMRD_g_Nz; k++)
            {
               adj_tre_max=0;
               adj_tre_min=0;
               for (i0=max(0,i-1); i0<=min(AMRD_g_Nx-1,i+1); i0++)
               {
                  for (j0=max(0,j-1); j0<=min(AMRD_g_Ny-1,j+1); j0++)
                  {
                     for (k0=max(0,k-1); k0<=min(AMRD_g_Nz-1,k+1); k0++)
                     {
                        if (i!=i0 || j!=j0 || k!=k0) 
                        {
                           adj_tre_max=max(adj_tre_max,AMRD_tre[i0+j0*AMRD_g_Nx+k0*AMRD_g_Nx*AMRD_g_Ny]);
                           adj_tre_min=min(adj_tre_min,AMRD_tre[i0+j0*AMRD_g_Nx+k0*AMRD_g_Nx*AMRD_g_Ny]);
                        }
                     }
                  }
               }
               ind=i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny;
               if ((min_max<0) && (adj_tre_min<0 || AMRD_tre[ind]<0)) AMRD_tre[ind]=min(AMRD_tre[ind],adj_tre_min+1);
               if (min_max>0) AMRD_tre[ind]=max(AMRD_tre[ind],adj_tre_max-1);
            }
         }
      }
   }
}

// A cell-centered version.
void bleed_tre_c(int max_buf, int min_max)
{
   int i,j,k,n,i0,j0,k0,ind;
   real adj_tre_max,adj_tre_min;

   for (n=0; n<max_buf; n++)
   {
      for (i=0; i<AMRD_g_Nx_c; i++)
      {
         for (j=0; j<AMRD_g_Ny_c; j++)
         {
            for (k=0; k<AMRD_g_Nz_c; k++)
            {
               adj_tre_max=0;
               adj_tre_min=0;
               for (i0=max(0,i-1); i0<=min(AMRD_g_Nx_c-1,i+1); i0++)
               {
                  for (j0=max(0,j-1); j0<=min(AMRD_g_Ny_c-1,j+1); j0++)
                  {
                     for (k0=max(0,k-1); k0<=min(AMRD_g_Nz_c-1,k+1); k0++)
                     {
                        if (i!=i0 || j!=j0 || k!=k0) 
                        {
                           adj_tre_max=max(adj_tre_max,AMRD_tre_c[i0+j0*AMRD_g_Nx_c+k0*AMRD_g_Nx_c*AMRD_g_Ny_c]);
                           adj_tre_min=min(adj_tre_min,AMRD_tre_c[i0+j0*AMRD_g_Nx_c+k0*AMRD_g_Nx_c*AMRD_g_Ny_c]);
                        }
                     }
                  }
               }
               ind=i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c;
               if ((min_max<0) && (adj_tre_min<0 || AMRD_tre_c[ind]<0)) AMRD_tre_c[ind]=min(AMRD_tre_c[ind],adj_tre_min+1);
               if (min_max>0) AMRD_tre_c[ind]=max(AMRD_tre_c[ind],adj_tre_max-1);
            }
         }
      }
   }
}

void zero_tre_ibc_a(void)
{
   int i,j,k,ip,jp,kp,im,jm,km,n;

   for (i=0; i<AMRD_g_Nx; i++)
   {
      ip=min(i+AMRD_TRE_ibc_a_buffer,AMRD_g_Nx-1);
      im=max(i-AMRD_TRE_ibc_a_buffer,0);
      for (j=0; j<AMRD_g_Ny; j++)
      {
         jp=min(j+AMRD_TRE_ibc_a_buffer,AMRD_g_Ny-1);
         jm=max(j-AMRD_TRE_ibc_a_buffer,0);
         for (k=0; k<AMRD_g_Nz; k++)
         {
            kp=min(k+AMRD_TRE_ibc_a_buffer,AMRD_g_Nz-1);
            km=max(k-AMRD_TRE_ibc_a_buffer,0);
            if ( (ip==(AMRD_g_Nx-1) && AMRD_cmask[ip+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                 (im==0 && AMRD_cmask[im+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                 ( AMRD_g_dim>1 &&
                  ( (jp==(AMRD_g_Ny-1) && AMRD_cmask[i+jp*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                    (jm==0 && AMRD_cmask[i+jm*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ) ) ||
                 ( AMRD_g_dim>2 && 
                  ( (kp==(AMRD_g_Nz-1) && AMRD_cmask[i+j*AMRD_g_Nx+kp*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                    (km==0 && AMRD_cmask[i+j*AMRD_g_Nx+km*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ) ) )
                 AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0;
         }
      }
   }
}

// The cell centered version.
void zero_tre_ibc_a_c(void)
{
   int i,j,k,ip,jp,kp,im,jm,km,n;

   // first TRE_ibc_buffer
   for (i=0; i<AMRD_g_Nx_c; i++)
   {
      ip=min(i+AMRD_TRE_ibc_a_buffer,AMRD_g_Nx_c-1);
      im=max(i-AMRD_TRE_ibc_a_buffer,0);
      for (j=0; j<AMRD_g_Ny_c; j++)
      {
         jp=min(j+AMRD_TRE_ibc_a_buffer,AMRD_g_Ny_c-1);
         jm=max(j-AMRD_TRE_ibc_a_buffer,0);
         for (k=0; k<AMRD_g_Nz_c; k++)
         {
            kp=min(k+AMRD_TRE_ibc_a_buffer,AMRD_g_Nz_c-1);
            km=max(k-AMRD_TRE_ibc_a_buffer,0);
            if ( (ip==(AMRD_g_Nx_c-1) && AMRD_cmask_c[ip+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                 (im==0 && AMRD_cmask_c[im+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                 ( AMRD_g_dim>1 &&
                  ( (jp==(AMRD_g_Ny_c-1) && AMRD_cmask_c[i+jp*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                    (jm==0 && AMRD_cmask_c[i+jm*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ) ) ||
                 ( AMRD_g_dim>2 && 
                  ( (kp==(AMRD_g_Nz_c-1) && AMRD_cmask_c[i+j*AMRD_g_Nx+kp*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                    (km==0 && AMRD_cmask_c[i+j*AMRD_g_Nx+km*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ) ) )
              AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0;
         }
      }
   }
}

void zero_f_tre_ibc(void)
{
   int i,j,k,ip,jp,kp,im,jm,km,n;

   // first TRE_ibc_buffer
   for (i=0; i<AMRD_g_Nx; i++)
   {
      ip=min(i+AMRD_TRE_ibc_buffer,AMRD_g_Nx-1);
      im=max(i-AMRD_TRE_ibc_buffer,0);
      for (j=0; j<AMRD_g_Ny; j++)
      {
         jp=min(j+AMRD_TRE_ibc_buffer,AMRD_g_Ny-1);
         jm=max(j-AMRD_TRE_ibc_buffer,0);
         for (k=0; k<AMRD_g_Nz; k++)
         {
            kp=min(k+AMRD_TRE_ibc_buffer,AMRD_g_Nz-1);
            km=max(k-AMRD_TRE_ibc_buffer,0);
            if ( (ip==(AMRD_g_Nx-1) && AMRD_cmask[ip+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                 (im==0 && AMRD_cmask[im+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                 ( AMRD_g_dim>1 &&
                  ( (jp==(AMRD_g_Ny-1) && AMRD_cmask[i+jp*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                    (jm==0 && AMRD_cmask[i+jm*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ) ) ||
                 ( AMRD_g_dim>2 && 
                  ( (kp==(AMRD_g_Nz-1) && AMRD_cmask[i+j*AMRD_g_Nx+kp*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                    (km==0 && AMRD_cmask[i+j*AMRD_g_Nx+km*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ) ) )
                 for (n=0; n<AMRD_num_f_tre_vars; n++) (AMRD_f_tre[n])[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0;
         }
      }
   }

   // now TRE_ibcp_buffer
   for (i=0; i<AMRD_g_Nx; i++)
   {
      ip=min(i+AMRD_TRE_ibcp_buffer,AMRD_g_Nx-1);
      im=max(i-AMRD_TRE_ibcp_buffer,0);
      for (j=0; j<AMRD_g_Ny; j++)
      {
         jp=min(j+AMRD_TRE_ibcp_buffer,AMRD_g_Ny-1);
         jm=max(j-AMRD_TRE_ibcp_buffer,0);
         for (k=0; k<AMRD_g_Nz; k++)
         {
            kp=min(k+AMRD_TRE_ibcp_buffer,AMRD_g_Nz-1);
            km=max(k-AMRD_TRE_ibcp_buffer,0);
            if ( (ip==(AMRD_g_Nx-1) || im==0) ||
                 ( AMRD_g_dim>1 &&
                  ( jp==(AMRD_g_Ny-1) || jm==0 ) ) ||
                 ( AMRD_g_dim>2 && 
                  ( kp==(AMRD_g_Nz-1) || km==0 ) ) )
                 for (n=0; n<AMRD_num_f_tre_vars; n++) (AMRD_f_tre[n])[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0;
         }
      }
   }

   // now TRE_exc_buffer ... i.e., zero tre where cmask is off
   for (i=0; i<AMRD_g_size; i++)
   {
      if (AMRD_cmask[i]==AMRD_CMASK_OFF) 
         for (n=0; n<AMRD_num_f_tre_vars; n++) (AMRD_f_tre[n])[i]=0;
   }
}

// The cell centered version.
void zero_f_tre_ibc_c(void)
{
   int i,j,k,ip,jp,kp,im,jm,km,n;

   // first TRE_ibc_buffer
   for (i=0; i<AMRD_g_Nx_c; i++)
   {
      ip=min(i+AMRD_TRE_ibc_buffer,AMRD_g_Nx_c-1);
      im=max(i-AMRD_TRE_ibc_buffer,0);
      for (j=0; j<AMRD_g_Ny_c; j++)
      {
         jp=min(j+AMRD_TRE_ibc_buffer,AMRD_g_Ny_c-1);
         jm=max(j-AMRD_TRE_ibc_buffer,0);
         for (k=0; k<AMRD_g_Nz_c; k++)
         {
            kp=min(k+AMRD_TRE_ibc_buffer,AMRD_g_Nz_c-1);
            km=max(k-AMRD_TRE_ibc_buffer,0);
            if ( (ip==(AMRD_g_Nx_c-1) && AMRD_cmask_c[ip+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                 (im==0 && AMRD_cmask_c[im+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                 ( AMRD_g_dim>1 &&
                  ( (jp==(AMRD_g_Ny_c-1) && AMRD_cmask_c[i+jp*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                    (jm==0 && AMRD_cmask_c[i+jm*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ) ) ||
                 ( AMRD_g_dim>2 && 
                  ( (kp==(AMRD_g_Nz_c-1) && AMRD_cmask_c[i+j*AMRD_g_Nx+kp*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ||
                    (km==0 && AMRD_cmask_c[i+j*AMRD_g_Nx+km*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF) ) ) )
              for (n=0; n<AMRD_num_f_tre_vars; n++) (AMRD_f_tre[n])[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0;
         }
      }
   }

   // now TRE_ibcp_buffer
   for (i=0; i<AMRD_g_Nx_c; i++)
   {
      ip=min(i+AMRD_TRE_ibcp_buffer,AMRD_g_Nx_c-1);
      im=max(i-AMRD_TRE_ibcp_buffer,0);
      for (j=0; j<AMRD_g_Ny_c; j++)
      {
         jp=min(j+AMRD_TRE_ibcp_buffer,AMRD_g_Ny_c-1);
         jm=max(j-AMRD_TRE_ibcp_buffer,0);
         for (k=0; k<AMRD_g_Nz_c; k++)
         {
            kp=min(k+AMRD_TRE_ibcp_buffer,AMRD_g_Nz_c-1);
            km=max(k-AMRD_TRE_ibcp_buffer,0);
            if ( (ip==(AMRD_g_Nx_c-1) || im==0) ||
                 ( AMRD_g_dim>1 &&
                  ( jp==(AMRD_g_Ny_c-1) || jm==0 ) ) ||
                 ( AMRD_g_dim>2 && 
                  ( kp==(AMRD_g_Nz_c-1) || km==0 ) ) )
              for (n=0; n<AMRD_num_f_tre_vars; n++) (AMRD_f_tre[n])[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0;
         }
      }
   }

   // now TRE_exc_buffer ... i.e., zero tre where cmask is off
   for (i=0; i<AMRD_g_size_c; i++)
   {
      if (AMRD_cmask_c[i]==AMRD_CMASK_OFF) 
         for (n=0; n<AMRD_num_f_tre_vars; n++) (AMRD_f_tre[n])[i]=0;
   }
}

//=============================================================================
// simple_cls surrounds tre by minimal bounding boxes, and then merges
// all clusters within cls_merge_dist grid-points in any direction. 
//
// for periodic dimensions, 'wrapping' clusters are merged.
//
// returns a (global) list of bboxes and levels ... caller must free after use,
//
// NOTE: this routine zeroes tre in the calculation process
//=============================================================================
#define MAX_BHS 6
void simple_cls(int L1, int L2, real **gbbox, int **glev, int *gnum)
{
   int i0,j0,k0,ltrace=0,i;
   real lbbox[2*PAMR_MAX_DIM*MAX_CLUSTERS];
   int lnum,lgnum=0,*lev,L,valid,plgnum;
   real *lc,*gc,*gc0,*pgc;
   static real lcount=-1;
   real bh_bbox[2*MAX_BHS*PAMR_MAX_DIM];
   int num_bhs=0;

   lcount++;

   *gnum=0;

   if (!(*gbbox=(real *)malloc(sizeof(real)*2*PAMR_MAX_DIM*MAX_CLUSTERS)))
      AMRD_stop("simple_cls: out of memory","");
   if (!(*glev=(int *)malloc(sizeof(int)*MAX_CLUSTERS)))
      AMRD_stop("simple_cls: out of memory","");

   lev=*glev;

   gc=*gbbox;
   pgc=0;
   plgnum=0;
   for (L=L1; L<=L2; L++)
   {
      lc=lbbox;
      lnum=0;
      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      while(valid)
      {
         ev_ldptr();
         if (!AMRD_using_cc_tre) while(find_seed(&i0,&j0,&k0)) { if (add_cluster(i0,j0,k0,lc)) {lc+=2*AMRD_dim; lnum++;} }
         if (AMRD_using_cc_tre) while(find_seed_c(&i0,&j0,&k0)) { if (add_cluster_c(i0,j0,k0,lc)) {lc+=2*AMRD_dim; lnum++;} }
         valid=PAMR_next_g();
      }

      IFLR debug_save_bbox_list(lbbox,0,L,lnum,lcount,"simple_cls_bf_merge");

      lgnum=MAX_CLUSTERS-*gnum; // how much cluster memory we have left
      if (!(PAMR_merge_bboxes(lbbox,lnum,gc,&lgnum,1))) AMRD_stop("PAMR merge error","");
      IFLR printf("rank=%i, lgnum after PAMR_merge_bboxes=%i\n",my_rank,lgnum);

      IFLR debug_save_bbox_list(gc,0,L,lgnum,lcount,"simple_cls_af_PAMR_merge");

      // if desired, make sure black holes are covered by single clusters

      merge_cls(gc,&lgnum,L);

      if (AMRD_TRE_sgpbh && AMRD_do_ex)
      {
         app_fill_bh_bboxes(bh_bbox,&num_bhs,MAX_BHS);
         if (num_bhs>0 && num_bhs<=MAX_BHS) extend_ex_cls(gc,lgnum,L,bh_bbox,num_bhs);
         merge_cls(gc,&lgnum,L);
      }

      clean_cls(gc,&lgnum,L,pgc,plgnum);

      IFLR printf("rank=%i, lgnum after merge_cls=%i\n",my_rank,lgnum);

      IFLR debug_save_bbox_list(gc,0,L,lgnum,lcount,"simple_cls_af_merge");

      pgc=gc; plgnum=lgnum;

      gc+=2*AMRD_dim*lgnum;
      *gnum+=lgnum;
      while(lgnum--) *(lev++)=L;
   }

   return;
}

int find_seed(int *i0, int *j0, int *k0)
{
   int i,j,k;

   for (i=0; i<AMRD_g_Nx; i++)
   {
      for (j=0; j<AMRD_g_Ny; j++)
      {
         for (k=0; k<AMRD_g_Nz; k++)
         {
            if (AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]>0.5)
            {
               *i0=i; *j0=j; *k0=k; return 1;
            }
         }
      }
   }

   return 0;
}

// A cell centered version
int find_seed_c(int *i0, int *j0, int *k0)
{
   int i,j,k;

   for (i=0; i<AMRD_g_Nx_c; i++)
   {
      for (j=0; j<AMRD_g_Ny_c; j++)
      {
         for (k=0; k<AMRD_g_Nz_c; k++)
         {
            if (AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]>0.0625)
            {
               *i0=i; *j0=j; *k0=k; return 1;
            }
         }
      }
   }

   return 0;
}

int add_cluster(int i0, int j0, int k0, real *bbox)
{
   int i,j,k,n[2*PAMR_MAX_DIM],np,imin,jmin,kmin,imax,jmax,kmax,lnp;

   AMRD_tre[i0+j0*AMRD_g_Nx+k0*AMRD_g_Nx*AMRD_g_Ny]=0;
   imin=imax=i0;
   jmin=jmax=j0;
   kmin=kmax=k0;
   np=1;
   while(np)
   {
      np=0;
      if (imin>0)
      {
         i=imin-1; lnp=0;
         for (j=jmin; j<=jmax; j++)
            for (k=kmin; k<=kmax; k++)
               if (AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]>0.5) { AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0; lnp++; }
         if (lnp) { imin--; np+=lnp; }
      }
      if (imax<(AMRD_g_Nx-1))
      {
         i=imax+1; lnp=0;
         for (j=jmin; j<=jmax; j++)
            for (k=kmin; k<=kmax; k++)
               if (AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]>0.5) { AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0; lnp++; }
         if (lnp) { imax++; np+=lnp; }
      }
      if (jmin>0)
      {
         j=jmin-1; lnp=0;
         for (i=imin; i<=imax; i++)
            for (k=kmin; k<=kmax; k++)
               if (AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]>0.5) { AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0; lnp++; }
         if (lnp) { jmin--; np+=lnp; }
      }
      if (jmax<(AMRD_g_Ny-1))
      {
         j=jmax+1; lnp=0;
         for (i=imin; i<=imax; i++)
            for (k=kmin; k<=kmax; k++)
               if (AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]>0.5) { AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0; lnp++; }
         if (lnp) { jmax++; np+=lnp; }
      }
      if (kmin>0)
      {
         k=kmin-1; lnp=0;
         for (i=imin; i<=imax; i++)
            for (j=jmin; j<=jmax; j++)
               if (AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]>0.5) { AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0; lnp++; }
         if (lnp) { kmin--; np+=lnp; }
      }
      if (kmax<(AMRD_g_Nz-1))
      {
         k=kmax+1; lnp=0;
         for (i=imin; i<=imax; i++)
            for (j=jmin; j<=jmax; j++)
               if (AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]>0.5) { AMRD_tre[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=0; lnp++; }
         if (lnp) { kmax++; np+=lnp; }
      }
   }

   n[0]=imin; n[1]=imax; 
   n[2]=jmin; n[3]=jmax; 
   n[4]=kmin; n[5]=kmax; 

   // min/max in case we have periodic boundaries, and drop slivers or 'smaller'

   for (i=0; i<AMRD_g_dim; i++) 
   { 
      bbox[2*i]=max(AMRD_base_bbox[2*i],(AMRD_g_x[i])[n[2*i]]); 
      bbox[2*i+1]=min(AMRD_base_bbox[2*i+1],(AMRD_g_x[i])[n[2*i+1]]); 
      if ((bbox[2*i+1]-bbox[2*i])<(((AMRD_g_x[i])[1]-(AMRD_g_x[i])[0])/2)) return 0;
      // if periodic then any cluster within 2*ghost_width of left boundary must be extended
      // to include this zone
      if (AMRD_periodic[i] && ((bbox[2*i]-AMRD_base_bbox[2*i])/AMRD_g_dx[i] < 2*AMRD_g_ghost_width[i]))
         bbox[2*i]=AMRD_base_bbox[2*i];
      if (AMRD_periodic[i] && ((bbox[2*i+1]-AMRD_base_bbox[2*i])/AMRD_g_dx[i] < 2*AMRD_g_ghost_width[i]))
         bbox[2*i+1]=AMRD_base_bbox[2*i]+2*AMRD_g_ghost_width[i]*AMRD_g_dx[i];
   }

   return 1;
}

int add_cluster_c(int i0, int j0, int k0, real *bbox)
{
   int i,j,k,n[2*PAMR_MAX_DIM],np,imin,jmin,kmin,imax,jmax,kmax,lnp;

   AMRD_tre_c[i0+j0*AMRD_g_Nx+k0*AMRD_g_Nx*AMRD_g_Ny]=0;
   imin=imax=i0;
   jmin=jmax=j0;
   kmin=kmax=k0;
   np=1;
   while(np)
   {
      np=0;
      if (imin>0)
      {
         i=imin-1; lnp=0;
         for (j=jmin; j<=jmax; j++)
            for (k=kmin; k<=kmax; k++)
               if (AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]>0.5) { AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0; lnp++; }
         if (lnp) { imin--; np+=lnp; }
      }
      if (imax<(AMRD_g_Nx_c-1))
      {
         i=imax+1; lnp=0;
         for (j=jmin; j<=jmax; j++)
            for (k=kmin; k<=kmax; k++)
               if (AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]>0.5) { AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0; lnp++; }
         if (lnp) { imax++; np+=lnp; }
      }
      if (jmin>0)
      {
         j=jmin-1; lnp=0;
         for (i=imin; i<=imax; i++)
            for (k=kmin; k<=kmax; k++)
               if (AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]>0.5) { AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0; lnp++; }
         if (lnp) { jmin--; np+=lnp; }
      }
      if (jmax<(AMRD_g_Ny_c-1))
      {
         j=jmax+1; lnp=0;
         for (i=imin; i<=imax; i++)
            for (k=kmin; k<=kmax; k++)
               if (AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]>0.5) { AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0; lnp++; }
         if (lnp) { jmax++; np+=lnp; }
      }
      if (kmin>0)
      {
         k=kmin-1; lnp=0;
         for (i=imin; i<=imax; i++)
            for (j=jmin; j<=jmax; j++)
               if (AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]>0.5) { AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0; lnp++; }
         if (lnp) { kmin--; np+=lnp; }
      }
      if (kmax<(AMRD_g_Nz_c-1))
      {
         k=kmax+1; lnp=0;
         for (i=imin; i<=imax; i++)
            for (j=jmin; j<=jmax; j++)
               if (AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]>0.5) { AMRD_tre_c[i+j*AMRD_g_Nx_c+k*AMRD_g_Nx_c*AMRD_g_Ny_c]=0; lnp++; }
         if (lnp) { kmax++; np+=lnp; }
      }
   }

   // Now we're translating to vertex-centered boxes.
   n[0]=imin; n[1]=imax+1; 
   n[2]=jmin; n[3]=jmax+1; 
   n[4]=kmin; n[5]=kmax+1; 

   // min/max in case we have periodic boundaries, and drop slivers or 'smaller'

   for (i=0; i<AMRD_g_dim; i++) 
   { 
      bbox[2*i]=max(AMRD_base_bbox[2*i],(AMRD_g_x[i])[n[2*i]]); 
      bbox[2*i+1]=min(AMRD_base_bbox[2*i+1],(AMRD_g_x[i])[n[2*i+1]]); 
      if ((bbox[2*i+1]-bbox[2*i])<(((AMRD_g_x[i])[1]-(AMRD_g_x[i])[0])/2)) return 0;
      // if periodic then any cluster within 2*ghost_width of left boundary must be extended
      // to include this zone
      if (AMRD_periodic[i] && ((bbox[2*i]-AMRD_base_bbox[2*i])/AMRD_g_dx[i] < 2*AMRD_g_ghost_width[i]))
         bbox[2*i]=AMRD_base_bbox[2*i];
      if (AMRD_periodic[i] && ((bbox[2*i+1]-AMRD_base_bbox[2*i])/AMRD_g_dx[i] < 2*AMRD_g_ghost_width[i]))
         bbox[2*i+1]=AMRD_base_bbox[2*i]+2*AMRD_g_ghost_width[i]*AMRD_g_dx[i];
   }

   return 1;
}

//=============================================================================
// cleans up the cluster list ... currently this means
//
// 1) deletes clusters that are too small
//
// 2) if a parent bbox list is supplied (pbbox,pnum), trims/deletes
//    any clusters that don't have a parent
//=============================================================================
void clean_cls(real *bbox, int *num, int lev, real *pbbox, int pnum)
{
   int i,j,k,del;
   real dx[PAMR_MAX_DIM],dt;
   int min_width[PAMR_MAX_DIM],found;

   PAMR_get_dxdt(lev+1,dx,&dt); // using lev+1, as these grids are going to be refined
   PAMR_get_min_width(min_width);

   // trim clusters that overlap a parent (or have none)
   if (pbbox)
   {
      i=0;
      while(i<*num)
      {
         // find parent
         k=0;
         found=0;
         while (k<pnum && !found)
         {
            found=1;
            for (j=0; j<AMRD_dim; j++) 
               if ( (((bbox[2*(j+AMRD_dim*i)+1]-pbbox[2*(j+AMRD_dim*k)])/dx[j]+1) < 1) ||
                    (((pbbox[2*(j+AMRD_dim*k)+1]-bbox[2*(j+AMRD_dim*i)])/dx[j]+1) < 1) ) found=0;
            if (!found) k++;
         }
         if (!found)
         {
            if (my_rank==0) printf("clean_cls: WARNING ... child cluster has no parent ... deleting\n");
            bbox[2*(AMRD_dim*i)+1]=bbox[2*(AMRD_dim*i)]-min_width[0]; // set x-width so that it's deleted below
         }
         else
         {
            for (j=0; j<AMRD_dim; j++) 
            {
               if ( (bbox[2*(j+AMRD_dim*i)]-pbbox[2*(j+AMRD_dim*k)]) < (-dx[j]/10))
               {
                  if (my_rank==0) printf("clean_cls: WARNING ... child cluster overlaps (left) parent cluster ... trimming [%e,%e,%e]\n",
                     bbox[2*(j+AMRD_dim*i)],pbbox[2*(j+AMRD_dim*k)],dx[j]);
                  bbox[2*(j+AMRD_dim*i)]=pbbox[2*(j+AMRD_dim*k)];
               }
               if ( (bbox[2*(j+AMRD_dim*i)+1]-pbbox[2*(j+AMRD_dim*k)+1]) > dx[j]/10)
               {
                  if (my_rank==0) printf("clean_cls: WARNING ... child cluster overlaps (right) parent cluster ... trimming [%e,%e,%e]\n",
                     bbox[2*(j+AMRD_dim*i)+1],pbbox[2*(j+AMRD_dim*k)+1],dx[j]);
                     
                  bbox[2*(j+AMRD_dim*i)+1]=pbbox[2*(j+AMRD_dim*k)+1];
               }
            }
         }
         i++;
      }
   }

   i=0;
   while(i<*num)
   {
      del=0;
      for (j=0; j<AMRD_dim; j++) 
         if (((bbox[2*(j+AMRD_dim*i)+1]-bbox[2*(j+AMRD_dim*i)])/dx[j]+1) < min_width[j]) del=1;

      if (del)
      {
         (*num)--;
         for (j=2*AMRD_dim*i; j<(2*AMRD_dim*(*num)); j++) bbox[j]=bbox[j+2*AMRD_dim];
      }
      else i++;
   }

}

//=============================================================================
// merges clusters within cls_merge_dist of one another
// (set cls_merge_dist=0 to have no merging)
//=============================================================================
void merge_cls(real *bbox, int *num, int lev)
{
   real *c1,*c2;
   int i1,i2,j,k,overlap,ltrace=0,j0,delta,i3,swap;
   real dx[PAMR_MAX_DIM],dt;
   real dw[PAMR_MAX_DIM];

   if (AMRD_cls_merge_dist<0) AMRD_stop("merge_cls: error ... cls_merge_dist must be >=0","");
   if (AMRD_cls_merge_dist==0) return;

   PAMR_get_dxdt(lev,dx,&dt);

   IFLR printf("\tmerge_cls: rank=%i, AMRD_cls_merge_dist=%i, num=%i, AMRD_dim=%i, dx=[%lf,%lf,%lf]\n",
              my_rank,AMRD_cls_merge_dist,*num,AMRD_dim,dx[0],dx[1],dx[2]);

   for (j=0; j<AMRD_dim; j++) dw[j]=(AMRD_base_bbox[2*j+1]-AMRD_base_bbox[2*j]);

   i1=0;
   
   delta=1;
   while(delta)
   {
      delta=0;
      i1=0;
      while(i1<*num)
      {
         i2=i1;
         while(i2<*num)
         {
            if (i1!=i2)
            {
               overlap=1;
               for (j=0; j<AMRD_dim; j++) 
                  if ( ! ((((bbox[2*(j+AMRD_dim*i2)] - bbox[2*(j+AMRD_dim*i1)+1])/dx[j]+1) < AMRD_cls_merge_dist) && 
                          (((bbox[2*(j+AMRD_dim*i1)] - bbox[2*(j+AMRD_dim*i2)+1])/dx[j]+1) < AMRD_cls_merge_dist) )) overlap=0;
               if (overlap)
               {
                  (*num)--;
                  for (j=0; j<AMRD_dim; j++) 
                  {
                     bbox[2*(j+AMRD_dim*i1)]=min(bbox[2*(j+AMRD_dim*i1)],bbox[2*(j+AMRD_dim*i2)]);
                     bbox[2*(j+AMRD_dim*i1)+1]=max(bbox[2*(j+AMRD_dim*i1)+1],bbox[2*(j+AMRD_dim*i2)+1]);
                  }
                  for (k=2*AMRD_dim*i2; k<(2*AMRD_dim*(*num)); k++) bbox[k]=bbox[k+2*AMRD_dim];
                  i2=i1;
                  delta++;
               }
            }

            // check if grids overlap along periodic directions ... in that case extend both grids
            // to effectively merge (a grid can overlap with itself in the periodic direction)
            
            for (j0=0; j0<AMRD_dim; j0++)
            {
               if (AMRD_periodic[j0])
               {
                  // first check by shifting i1, then i2
                  swap=2;
                  while(swap--)
                  {
                     i3=i1; i1=i2; i2=i3;

                     overlap=1;
                     for (j=0; j<AMRD_dim; j++) 
                     {
                        if ((j!=j0) &&
                            ( ! ((((bbox[2*(j+AMRD_dim*i2)] - bbox[2*(j+AMRD_dim*i1)+1])/dx[j]+1) < AMRD_cls_merge_dist) && 
                                (((bbox[2*(j+AMRD_dim*i1)] - bbox[2*(j+AMRD_dim*i2)+1])/dx[j]+1) < AMRD_cls_merge_dist) ))) overlap=0;
                        if ((j==j0) &&
                            ( ! ((((bbox[2*(j+AMRD_dim*i2)] - bbox[2*(j+AMRD_dim*i1)+1]-dw[j])/dx[j]+1) < AMRD_cls_merge_dist) && 
                                (((bbox[2*(j+AMRD_dim*i1)]+dw[j] - bbox[2*(j+AMRD_dim*i2)+1])/dx[j]+1) < AMRD_cls_merge_dist) ))) overlap=0;
                     }
   
                     if (overlap)
                     {
                        for (j=0; j<AMRD_dim; j++) 
                        {
                           if (j!=j0)
                           {
                              if (bbox[2*(j+AMRD_dim*i1)]!=bbox[2*(j+AMRD_dim*i2)]) delta++;
                              if (bbox[2*(j+AMRD_dim*i1)+1]!=bbox[2*(j+AMRD_dim*i2)+1]) delta++;
                              bbox[2*(j+AMRD_dim*i1)]=bbox[2*(j+AMRD_dim*i2)]=min(bbox[2*(j+AMRD_dim*i1)],bbox[2*(j+AMRD_dim*i2)]);
                              bbox[2*(j+AMRD_dim*i1)+1]=bbox[2*(j+AMRD_dim*i2)+1]=max(bbox[2*(j+AMRD_dim*i1)+1],bbox[2*(j+AMRD_dim*i2)+1]);
                           }
                           else
                           {
                              if (bbox[2*(j+AMRD_dim*i1)]!=AMRD_base_bbox[2*j]) delta++;
                              if (bbox[2*(j+AMRD_dim*i2)+1]!=AMRD_base_bbox[2*j+1]) delta++;
                              bbox[2*(j+AMRD_dim*i1)]=AMRD_base_bbox[2*j];
                              bbox[2*(j+AMRD_dim*i2)+1]=AMRD_base_bbox[2*j+1];
                           }
                        }
                     }
                  }
               }
            }
            i2++;
         } 
         i1++;
      }
   }
}

//=============================================================================
// if any of bbox intersects any of (bh_bbox+TRE_buffer), 
// then corresponding bh_bbox is or'd with bbox (unless the result is not
// contained within any parent sgh)
//
// !! NOTE: bh_bbox is *not* expected to be aligned with grid boundaries,
//          hence we align it here
//
//          ALSO ... routine currently assumes non-overlapping sgh clusters
//=============================================================================
#define CLS_MAX_GRIDS 128
void extend_ex_cls(real *bbox, int num, int lev, real *bh_bbox, int num_bhs)
{
   real *c1,*c2;
   int i1,i2,j,k,overlap,ltrace=0,ind1,ind2,num_pc,i,contained,any_contained;
   real dx[PAMR_MAX_DIM],dt,delta,sgh[2*PAMR_MAX_DIM*CLS_MAX_GRIDS];
   static int first=1;

   if (first && AMRD_TRE_buffer<1) printf("extend_ex_cls: warning ... AMRD_TRE_buffer should be >=1 with TRE_sgpbh");
   first=0;
   if (AMRD_cls_merge_dist<1) AMRD_stop("extend_ex_cls: cls_merge_dist must be >=1 when TRE_sgpbh is on","");

   PAMR_get_dxdt(lev,dx,&dt);

   IFLR printf("extend_ex_cls: num=%i, num_bhs=%i\n",num,num_bhs);

   num_pc=CLS_MAX_GRIDS;
   if (!PAMR_get_sgh(sgh,&num_pc,lev)) AMRD_stop("extend_ex_cls ... error ... level empty!!!\n","");

   for (i1=0; i1<num; i1++)
   {
      for (i2=0; i2<num_bhs; i2++)
      {
         overlap=1;
         if (AMRD_dim>=2 && ltrace && my_rank==0)
         {
            ind1=2*AMRD_dim*i1;
            ind2=2*AMRD_dim*i2;
            printf("extend_ex_cls: bbox[%i]=[%lf,%lf],[%lf,%lf]",i1,bbox[ind1],bbox[ind1+1],bbox[ind1+2],bbox[ind1+3]);
            if (AMRD_dim>2) printf(",[%lf,%lf]",bbox[ind1+4],bbox[ind1+5]);
            printf("\n               bh_bbox[%i]=[%lf,%lf],[%lf,%lf]",i2,bh_bbox[ind2],bh_bbox[ind2+1],bh_bbox[ind2+2],bh_bbox[ind2+3]);
            if (AMRD_dim>2) printf("[%lf,%lf]",bh_bbox[ind2+4],bh_bbox[ind2+5]);
            printf("\n");
         }
         for (j=0; j<AMRD_dim; j++)
            if ( ! ((((bh_bbox[2*(j+AMRD_dim*i2)] - bbox[2*(j+AMRD_dim*i1)+1])/dx[j]+1) < AMRD_TRE_buffer) && 
                    (((bbox[2*(j+AMRD_dim*i1)] - bh_bbox[2*(j+AMRD_dim*i2)+1])/dx[j]+1) < AMRD_TRE_buffer) )) overlap=0;

         if (overlap)
         {
            any_contained=0;
            for (i=0; i<num_pc && !any_contained; i++)
            {
               contained=1;
               // allow bh's to touch physical boundaries (for axisymmetry, black strings, etc.)
               for (j=0; j<AMRD_dim; j++)
                  if ( 
                       ( ( (int)((bh_bbox[2*(j+AMRD_dim*i2)] - sgh[2*(j+AMRD_dim*i)])/dx[j]+0.5) < 0) ||
                         ( ( (int)((bh_bbox[2*(j+AMRD_dim*i2)] - sgh[2*(j+AMRD_dim*i)])/dx[j]+0.5) < AMRD_TRE_buffer) &&
                           ( (int)((bh_bbox[2*(j+AMRD_dim*i2)] - AMRD_base_bbox[2*j])/dx[j]+0.5) > 0) ) ) || 
                       ( ( (int)((sgh[2*(j+AMRD_dim*i)+1] - bh_bbox[2*(j+AMRD_dim*i2)+1])/dx[j]+0.5) < 0) ||
                         ( ( (int)((sgh[2*(j+AMRD_dim*i)+1] - bh_bbox[2*(j+AMRD_dim*i2)+1])/dx[j]+0.5) < AMRD_TRE_buffer) &&
                           ( (int)((AMRD_base_bbox[2*j+1] - bh_bbox[2*(j+AMRD_dim*i2)+1])/dx[j]+0.5) > 0) ) ) ) contained=0;
               if (contained) any_contained=1;

               if (AMRD_dim>=2 && ltrace && my_rank==0)
               {
                  ind1=2*AMRD_dim*i;
                  printf("             : sgh[%i]=[%lf,%lf],[%lf,%lf]",i,sgh[ind1],sgh[ind1+1],sgh[ind1+2],sgh[ind1+3]);
                  if (AMRD_dim>2) printf("[%lf,%lf]",sgh[ind1+4],sgh[ind1+5]);
                  printf(" (contained=%i)\n",contained);
               }
            }
            if (!any_contained) overlap=0;
         }

         if (overlap)
         {
            for (j=0; j<AMRD_dim; j++) 
            {
               delta=bbox[2*(j+AMRD_dim*i1)]-(bh_bbox[2*(j+AMRD_dim*i2)]-AMRD_TRE_buffer*dx[j]);
               if (delta>0) bbox[2*(j+AMRD_dim*i1)]-=dx[j]*((int)(delta/dx[j]+0.5));
               bbox[2*(j+AMRD_dim*i1)]=max(bbox[2*(j+AMRD_dim*i1)],AMRD_base_bbox[2*j]);

               delta=(bh_bbox[2*(j+AMRD_dim*i2)+1]+AMRD_TRE_buffer*dx[j])-bbox[2*(j+AMRD_dim*i1)+1];
               if (delta>0) bbox[2*(j+AMRD_dim*i1)+1]+=dx[j]*((int)(delta/dx[j]+0.5));
               bbox[2*(j+AMRD_dim*i1)+1]=min(bbox[2*(j+AMRD_dim*i1)+1],AMRD_base_bbox[2*j+1]);
            }
         }
         if (AMRD_dim>=2 && ltrace && my_rank==0)
         {
            printf("overlap=%i: bbox[%i]=[%lf,%lf],[%lf,%lf]",overlap,i1,bbox[ind1],bbox[ind1+1],bbox[ind1+2],bbox[ind1+3]);
            if (AMRD_dim>2) printf("[%lf,%lf]",bbox[ind1+4],bbox[ind1+5]);
            printf("\n");
         }
      }
   }
}

//=============================================================================
// a width table, defining the allowable dimension width's for grid
// dimensions that can be fully coarsened within the MG hierarhcy -->
// ([min_width | min_width+1 | ... | 2*(min_width-1)] - 1)*2^n + 1, 
// n=0,1,2,...
//=============================================================================

#define AMRD_MAX_WTAB_NUM 60

int wtab[PAMR_MAX_DIM][AMRD_MAX_WTAB_NUM];
int wtab_min_width[PAMR_MAX_DIM][AMRD_MAX_WTAB_NUM];
int wtab_n[PAMR_MAX_DIM][AMRD_MAX_WTAB_NUM];

void fill_wtab(void)
{
   int size,min_size,i,j,k,min_width[PAMR_MAX_DIM],n;
   static int first=1;
   int ltrace=0;

   if (!(first)) return; first=0;

   PAMR_get_MG_coarse_width(min_width);

   for (j=0; j<AMRD_dim; j++)
   {
      wtab[j][0]=min_width[j];
      wtab_n[j][0]=0;
      wtab_min_width[j][0]=min_width[j];
      for (i=1; i<AMRD_MAX_WTAB_NUM; i++) 
      {
         wtab[j][i]=(wtab[j][i-1]-1)*2+1;
         wtab_n[j][i]=i;
         wtab_min_width[j][i]=min_width[j];
      }
      size=min_size=min_width[j]+1;
      n=0;
      while(min_size<(min_width[j]*2-1))
      {
         i=0;
         while(i<AMRD_MAX_WTAB_NUM)
         {
            while(i<(AMRD_MAX_WTAB_NUM-1) && size>wtab[j][i]) i++;
            for (k=AMRD_MAX_WTAB_NUM-1; k>i; k--) 
            {
               wtab[j][k]=wtab[j][k-1];
               wtab_n[j][k]=wtab_n[j][k-1];
               wtab_min_width[j][k]=wtab_min_width[j][k-1];
            }
            wtab[j][i]=size;
            wtab_n[j][i]=n;
            wtab_min_width[j][i]=min_size;
            size=(size-1)*2+1; n++; i++;
         }
         min_size++;
         size=min_size;
         n=0;
      }
      
      if (ltrace)
      {
         printf("\nwidth table for dim %i:\n\nwidth\t\tmin_width\t\t n\n=====\t\t=========\t===\n",j);
         for (i=0; i<AMRD_MAX_WTAB_NUM; i++) printf("%i\t\t%i\t\t%i\n",wtab[j][i],wtab_min_width[j][i],wtab_n[j][i]);
         printf("----------------------------------------------------------\n");
      }
   }
}

//=============================================================================
// The following adjusts a cluster list for consistency with MG, by aligning
// them along grid-intervals of 
//
// a) 2^net_rhosp for overlapping clusters, where net_rhosp is the 
//    net refinement ratio between the given level and the coarsest
//    level of the given cluster-island.
//
// b) and adjusts the sizes so that the smallest dimension of a cluster
// factors into a wtab[] width
//
// cls_align_mode determines how clusters are adjusted to force alignment.
// current modes:
//
// CLS_ALIGN_SHRINK : clusters are shrunk 
// CLS_ALIGN_EXPAND : clusters are expanded 
//
// This routine assumes the clusters are saved level-by-level, in order,
// in gbbox, from L1 to L2; AND, the clusters at level i
// will become grids at level i+1.
//
// If (AMRD_num_elliptic_vars==0) then no alignment is done
//
// NOTES: the current algorithm is rather simple ... i.e. no gaurantees
//        that grid-overlap/isolation will be preserved, etc.
//        ALSO, insensitive to parent grid positions ... there may be
//        (almost certainly will be with overlapping grids) circumstances
//        where this routine will move a child outside of the domain
//        of a parent ... improve this routine when those situations arise.
//=============================================================================
void adjust_cls(int L1, int L2, real **gbbox, int **glev, int *gnum)
{
   int L,ltrace=0,plgnum;
   real *bbox=*gbbox,*p,*pgc;
   int *island_no,*net_rhosp_n,i,j,is,ie,align_mod,inside,k;
   int num=*gnum,*lev=*glev,*l,tnum=0,delta,shape,rem,*q;
   int rho_sp[PAMR_MAX_LEVS],rho_tm[PAMR_MAX_LEVS];
   real dx[PAMR_MAX_DIM],dt,ds;
   static int lcount=-1;

   lcount++;

   IFLR printf(">> adjust_cls(%i,%i,.,.,%i)\n",L1,L2,*gnum);

   if (AMRD_num_elliptic_vars==0) return;

   PAMR_get_rho(rho_sp,rho_tm,PAMR_MAX_LEVS);

   if (AMRD_cls_align_mode!=CLS_ALIGN_EXPAND && AMRD_cls_align_mode!=CLS_ALIGN_SHRINK) 
      AMRD_stop("adjust_cls: unknown cls_align_mode","");

   fill_wtab();

   if (!(island_no=(int *)malloc(sizeof(int)*(num)))) AMRD_stop("adjust_cls: out of memory\n","");
   if (!(net_rhosp_n=(int *)malloc(sizeof(int)*(num)))) AMRD_stop("adjust_cls: out of memory\n","");

   tnum=0;
   l=lev;
   pgc=0;
   plgnum=0;
   for (L=L1; L<=L2; L++)
   {
      PAMR_get_dxdt(L+1,dx,&dt);

      num=0; while(tnum<*gnum && *l==L) { l++; num++; tnum++; }
      if (num==0) AMRD_stop("adjust_cls: error ... a level is empty of clusters\n","");

      cls_island_info(bbox,num,island_no,net_rhosp_n,L);

      for (i=0; i<num; i++)
      {
         for (j=0; j<AMRD_dim; j++)
         {
            shape=(bbox[2*i*AMRD_dim+2*j+1]-bbox[2*i*AMRD_dim+2*j])/dx[j]+1.5;

            align_mod=pow(2,net_rhosp_n[i]);

            rem=(shape-1)%align_mod;

            //-----------------------------------------------------------------
            // 1. Adjust size to allowable one:
            //-----------------------------------------------------------------
            if ( rem && ( AMRD_cls_align_mode==CLS_ALIGN_EXPAND || 
                         ((shape-rem-1)/align_mod < wtab[j][0]) ) )
               delta=align_mod-rem;
            else
               delta=-rem;

            //-----------------------------------------------------------------
            // 2. Align grid : move edge closest to an alignment point 
            // (defined relative to the global bounding boxes),
            // and trim/expand other edge. For isolated clusters, delta
            // could be much larger that align_mod ... in that case trim/expand
            // both edges
            //-----------------------------------------------------------------
            if (island_no[i]<=num) align_mod=rho_sp[L-1];
            is=(bbox[2*i*AMRD_dim+2*j]-AMRD_base_bbox[2*j])/dx[j]+1.5;
            ie=(bbox[2*i*AMRD_dim+2*j+1]-AMRD_base_bbox[2*j])/dx[j]+1.5;
            if (((is-1)%align_mod) < (align_mod-((ie-1)%align_mod)))
            {
               //--------------------------------------------------------------
               // move to 'is'
               //--------------------------------------------------------------
               bbox[2*i*AMRD_dim+2*j]-=dx[j]*((is-1)%align_mod);
               bbox[2*i*AMRD_dim+2*j+1]-=dx[j]*((is-1)%align_mod);
               if (delta/align_mod>=2)
               {
                  bbox[2*i*AMRD_dim+2*j]-=(int)(delta/align_mod/2)*align_mod*dx[j];
                  bbox[2*i*AMRD_dim+2*j+1]+=(int)(delta/align_mod/2)*align_mod*dx[j];
                  delta-=(delta/abs(delta))*(int)((delta/align_mod/2))*align_mod*2;
               }
               bbox[2*i*AMRD_dim+2*j+1]+=delta*dx[j];
            }
            else
            {
               //--------------------------------------------------------------
               // move to 'ie'
               //--------------------------------------------------------------
               bbox[2*i*AMRD_dim+2*j]+=dx[j]*(align_mod-(ie-1)%align_mod);
               bbox[2*i*AMRD_dim+2*j+1]+=dx[j]*(align_mod-(ie-1)%align_mod);
               if (delta/align_mod>=2)
               {
                  bbox[2*i*AMRD_dim+2*j]-=(int)(delta/align_mod/2)*align_mod*dx[j];
                  bbox[2*i*AMRD_dim+2*j+1]+=(int)(delta/align_mod/2)*align_mod*dx[j];
                  delta-=(delta/abs(delta))*(int)((delta/align_mod/2))*align_mod*2;
               }
               bbox[2*i*AMRD_dim+2*j]-=delta*dx[j];
            }
            //-----------------------------------------------------------------
            // The only check we currently do w.r.t parent grids, is 
            // wether the grid extends beyond the computational domain. This
            // one is important, becase TRE_ibc_buffer does not apply near 
            // physical boundaries.
            //-----------------------------------------------------------------
            ds=bbox[2*i*AMRD_dim+2*j]-AMRD_base_bbox[2*j];
            if (ds<(-dx[j]/2))
            {
               bbox[2*i*AMRD_dim+2*j]+=(-ds);
               bbox[2*i*AMRD_dim+2*j+1]+=(-ds);
            }
            ds=AMRD_base_bbox[2*j+1]-bbox[2*i*AMRD_dim+2*j+1];
            if (ds<(-dx[j]/2))
            {
               bbox[2*i*AMRD_dim+2*j]-=(-ds);
               bbox[2*i*AMRD_dim+2*j+1]-=(-ds);
            }
         }
      }
      //-----------------------------------------------------------------------
      // now, eliminate any redundent boxes ... this should not occur
      // frequently if simple_cls() produced cluster list, for it
      // doesn't produce redundent clusters
      //-----------------------------------------------------------------------
      i=0;
      while(i<num)
      {
         j=0;
         while(j<num)
         {
            if (i!=j)
            {
               // check wether i is entirely in j
               inside=1;
               for (k=0; k<AMRD_dim; k++) 
                  if ((bbox[2*j*AMRD_dim+2*k]-bbox[2*i*AMRD_dim+2*k])>dx[k]/2 || 
                      (bbox[2*i*AMRD_dim+2*k+1]-bbox[2*j*AMRD_dim+2*k+1])>dx[k]/2) inside=0;
               if (inside) {lev[i]=-1; j=num;} // mark for deletion, actually delete later
            }
            j++;
         }
         i++;
      }

      // sometimes adjustment will trim a parent grid so that the
      // child grid isn't entirely contained within the parent ... the
      // following trims it back (at the risk of having a grid
      // not consistent with width-table).

      clean_cls(bbox,&num,L,pgc,plgnum);
      pgc=bbox;plgnum=num;

      bbox+=(2*AMRD_dim*num);
      lev+=num;
   }

   p=bbox=*gbbox; q=lev=*glev;
   num=0;
   for (i=0; i<*gnum; i++)
   {
      if ((*glev)[i]!=-1) 
      {
         for (j=0;j<2*AMRD_dim;j++) *bbox++=*p++; 
         *lev++=*q++;
         num++; 
      }
      else 
      {
         p+=2*AMRD_dim;
         q++;
      }
   }
   *gnum=num;

   IFLR debug_save_bbox_list(*gbbox,*glev,0,*gnum,lcount,"adjust_cls_");

   free(island_no);
   free(net_rhosp_n);

   IFLR printf("<< adjust_cls(%i,%i,.,.,%i)\n",L1,L2,*gnum);
}

//=============================================================================
// the following fills in island_no and net_rhosp_n (number of refinements
// for the corresponding island) ... island numbers <=num are isolated clusters
//=============================================================================
void cls_island_info(real *bboxes, int num, int *island_no, int *net_rhosp_n, int L)
{
   int i,j,k,delta=1,shape,overlap;
   real dx[PAMR_MAX_DIM],dt;

   PAMR_get_dxdt(L+1,dx,&dt);
   
   for (i=0; i<num; i++)
   {
      island_no[i]=i+1;
      net_rhosp_n[i]=100;
      for (j=0; j<AMRD_dim; j++) 
      {
         shape=(bboxes[2*i*AMRD_dim+2*j+1]-bboxes[2*i*AMRD_dim+2*j])/dx[j]+1.5;
         k=0;
         if (shape<wtab[j][k]) AMRD_stop("cls_island_info: error ... shape < min_width\n","");
         while(shape>wtab[j][k]) k++;
         if (shape==wtab[j][k] || AMRD_cls_align_mode==CLS_ALIGN_EXPAND)
         {
            net_rhosp_n[i]=min(net_rhosp_n[i],wtab_n[j][k]);
         }
         else 
         {
            net_rhosp_n[i]=min(net_rhosp_n[i],wtab_n[j][k-1]);
         }
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
               for (k=0; k<AMRD_dim; k++)
                  if ((bboxes[2*i*AMRD_dim+2*k]-bboxes[2*j*AMRD_dim+2*k+1]>dx[k]/2) ||
                      (bboxes[2*i*AMRD_dim+2*k+1]-bboxes[2*j*AMRD_dim+2*k]>dx[k]/2)) overlap=0;
               if (overlap)
               {
                  delta++;
                  island_no[i]=island_no[j]=min(island_no[i],island_no[j]);
                  if (island_no[i]<=num) island_no[i]=island_no[j]=island_no[i]+num;
                  net_rhosp_n[i]=net_rhosp_n[j]=min(net_rhosp_n[i],net_rhosp_n[j]);
               }
            }
         }
      }
   }
}

//=============================================================================
// clears all f_tre variables on level L
//=============================================================================
void zero_f_tre(int L)
{
   int valid,i,j;

   valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
   while(valid)
   {
      ev_ldptr();
      // Achtung:  this loop will need to be reorganized if there
      // are both cell and vertex centered tre variables.
      // Make sure to check the other occurrences of AMRD_g_size
      // in this file.
      for (j=0; j<AMRD_num_f_tre_vars; j++) {
	if (PAMR_var_type(AMRD_f_tre_vars[j])==PAMR_VERTEX_CENTERED) {
	  for (i=0; i<AMRD_g_size; i++) {
	    (AMRD_f_tre[j])[i]=0;
	  }
	} else {
	  for (i=0; i<AMRD_g_size_c; i++) {
	    (AMRD_f_tre[j])[i]=0;
	  }
	}
      }
      valid=PAMR_next_g();
   }
}
