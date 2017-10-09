#ifndef _M_UTIL_R8
#define _M_UTIL_R8 1
/*======================================================================================== */
/* m_util_r8.h                                                                             */
/*                                                                                         */
/* c-prototypes for the fortran library                                                    */
/*======================================================================================== */
#include "internal_opts.h"

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

void dmfwr3d_(real *u1, real *u2, int *nx1, int *ny1, int *nz1, 
             int *nx2, int *ny2, int *nz2, int *i1, int *j1, int *k1, 
             int *i2, int *j2, int *k2, int *ni, int *nj, int *nk, 
             int *si1, int *sj1, int *sk1, int *si2, int *sj2, int *sk2,
             real *chr1, real *ex, int *do_ex);

void dmdiss3d_(real *f, real *work, real *eps, int *do_bdy, 
               int *phys_bdy_type, int *even, int *odd,
               real *mask, real *mask_off, int *nx, int *ny, int *nz,
               real *chr, real *ex, int *do_ex);

void linbnd3d_(real *f, int *phys_bdy, int *which_bnd, int *Nx, int *Ny, int*Nz,
               real *chr, real *ex, int *do_ex);

void dmrepop3d1_(real *f, real *chr, real *ex, int *io, int *nx, int *ny, int *nz);

void dmrepop3dc1_(real *chr, real *ex, int *nx, int *ny, int *nz);

#endif /* _M_UTIL_R8 */
