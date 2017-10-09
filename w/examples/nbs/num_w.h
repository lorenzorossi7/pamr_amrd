#ifndef _NUM
#define _NUM 1
/*======================================================================================== */
/* num.h                                                                                   */
/*                                                                                         */
/* c-prototypes for the fortran routines                                                   */
/*======================================================================================== */
void gauss3d_(real *f, real *A, real *r0, real *delta, real *xu0, real *yu0, real *zu0,
              real *ex, real *ey, real *ez, real *x, real *y, real *z, int *Nx, int *Ny, int *Nz);

void lop_(real *LV, real *V, real *phi_r, real *phi_i, real *cmask, 
          real *x, real *y, real *z, int *Nx, int *Ny, int *Nz);

void residual_(real *res, real *rhs, real *V, real *phi_r, real *phi_i, real *cmask, 
               real *x, real *y, real *z, real *norm, int *Nx, int *Ny, int *Nz);

void relax_(real *V, real *rhs, real *phi_r, real *phi_i, real *cmask, int *phys_bdy,
            real *x, real *y, real *z, real *norm, int *Nx, int *Ny, int *Nz);

void se_ires__(real *res, real *l2norm, real *phi_rn, real *phi_rnp1, real *phi_in, real *phi_inp1, 
               real *Vn, real *Vnp1, real *cmask, real *x, real *y, real *z, real *dt, int *Nx, int *Ny, int *Nz);

void phi_1step_cnc__(real *phi_rn, real *phi_rnp1, real *phi_in, real *phi_inp1, real *Vn, real *Vnp1, 
                     real *cmask, real *x, real *y, real *z, real *dt, int *Nx, int *Ny, int *Nz);
#endif /*_NUM */
