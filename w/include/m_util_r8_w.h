#ifndef _M_UTIL_R8
#define _M_UTIL_R8 1
/*======================================================================================== */
/* m_util_r8.h                                                                             */
/*                                                                                         */
/* c-prototypes for the fortran library                                                    */
/*======================================================================================== */
#include "internal_opts_w.h"

// following are vertex-centered routines

void dmcopy3d_(real *u1, real *u2, int *nx1, int *ny1, int *nz1, int *nx2, int *ny2,
              int *nz2, int *i1, int *j1, int *k1, int *i2, int *j2, int *k2,
              int *ni, int *nj, int *nk,
              int *si1, int *sj1, int *sk1, int *si2, int *sj2, int *sk2);

void dminterp3d_(real *u1, real *u2, int *nx1, int *ny1, int *nz1, 
                int *nx2, int *ny2, int *nz2, int *i1, int *j1, int *k1, 
                int *i2, int *j2, int *k2, int *ni, int *nj, int *nk, 
                int *si1, int *sj1, int *sk1, int *irho, int *jrho, int *krho, 
                int *ord, real *chr1, real *ex, int *do_ex);

void dmhwr3d_(real *u1, real *u2, int *nx1, int *ny1, int *nz1, 
             int *nx2, int *ny2, int *nz2, int *i1, int *j1, int *k1, 
             int *i2, int *j2, int *k2, int *ni, int *nj, int *nk, 
             int *si1, int *sj1, int *sk1, int *si2, int *sj2, int *sk2,
             real *chr1, real *ex, int *do_ex);

void dmhwrbdy3d_(real *u1, real *u2, int *nx1, int *ny1, int *nz1,
             int *nx2, int *ny2, int *nz2, int *i1, int *j1, int *k1,
             int *i2, int *j2, int *k2, int *ni, int *nj, int *nk,
             int *si1, int *sj1, int *sk1, int *si2, int *sj2, int *sk2,
             real *chr1, real *ex, int *do_ex);

void dmfwr3d_(real *u1, real *u2, int *nx1, int *ny1, int *nz1, 
             int *nx2, int *ny2, int *nz2, int *i1, int *j1, int *k1, 
             int *i2, int *j2, int *k2, int *ni, int *nj, int *nk, 
             int *si1, int *sj1, int *sk1, int *si2, int *sj2, int *sk2,
             real *chr1, real *ex, int *do_ex);

void dmdiss1d_(real *f, real *work, real *eps, int *do_bdy, 
               int *phys_bdy_type, int *even, int *odd,
               real *mask, real *mask_off, int *nx, 
               real *chr, real *ex, int *do_ex, int *use_6th);

void dmdiss3d_(real *f, real *work, real *eps, int *do_bdy, 
               int *phys_bdy_type, int *even, int *odd,
               real *mask, real *mask_off, int *nx, int *ny, int *nz,
               real *chr, real *ex, int *do_ex, int *use_6th, int *stride);

void dmdiss3d_64th_(real *f, real *work, real *eps, int *do_bdy,
               int *phys_bdy_type, int *even, int *odd,
               real *mask, real *mask_off, int *nx,int *ny, int *nz,
               real *chr, real *ex, int *do_ex, int *ord, int *stride);
 
void linbnd3d_(real *f, int *phys_bdy, int *which_bnd, int *Nx, int *Ny, int*Nz,
               real *chr, real *ex, int *do_ex);

void dmrepop3d1_(real *f, real *chr, real *ex, int *io, int *nx, int *ny, int *nz);

void dmrepop3dc1_(real *chr, real *ex, int *nx, int *ny, int *nz);

// following are cell-centered routines

void dminterp3d_c_(real *u1, real *u2, int *nx1, int *ny1, int *nz1, 
                int *nx2, int *ny2, int *nz2, int *i1, int *j1, int *k1, 
                int *i2, int *j2, int *k2, int *ni, int *nj, int *nk, 
                int *si1, int *sj1, int *sk1, int *irho, int *jrho, int *krho, 
                int *ord, real *chr1, real *ex, int *do_ex, int *dim);

void dmrepop3d1_c_(real *f, real *chr, real *ex, int *io, int *nx, int *ny, int *nz);

void dm_v_to_c_(real *f_v,real *f_c,real *chr_v,real *chr_c,
		int *nxv, int *nyv, int *nzv, int *nxc, int *nyc, int *nzc,
		int *ord, real *ex, int *do_ex, int *dim);

void dm_c_to_v_(real *f_c,real *f_v,real *chr_c,real *chr_v,
		int *nxc, int *nyc, int *nzc, int *nxv, int *nyv, int *nzv,
		int *ord, real *ex, int *do_ex, int *dim);

void dmavg3d_c_( real *u_f, real *u_c, real *chr, int *nx_f, int *ny_f,
		 int *nz_f, int *nx_c, int *ny_c, int *nz_c, int *i_f, 
		 int *j_f, int *k_f, int *i_c, int *j_c, int *k_c, int *ni_c,
		 int *nj_c, int *nk_c, int *rho, int *dim, real *ex, int *do_ex);

void dmadd3d_c_( real *u_f, real *u_c, int *nx_f, int *ny_f,
		 int *nz_f, int *nx_c, int *ny_c, int *nz_c, int *i_f, 
		 int *j_f, int *k_f, int *i_c, int *j_c, int *k_c, int *ni_c,
		 int *nj_c, int *nk_c, int *rho, int *dim);

#endif /* _M_UTIL_R8 */
