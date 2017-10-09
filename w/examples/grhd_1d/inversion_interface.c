#include "u2p_util.h"

int Utoprim_2d(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
               FTYPE gdet, FTYPE prim[NPR], FTYPE GAMMA);

void inversion_interface_(int *retval, double *rho, double *P, 
			  double *vx, double *vy, double *vz, 
			  double *rho_s, double *s_t, 
			  double *s_x, double *s_y, double *s_z,
			  double *g_tt, double *g_tx, double *g_ty, double *g_tz,
			  double *g_xx, double *g_xy, double *g_xz, 
			  double *g_yy, double *g_yz, double *g_zz, 
			  double *gup_tt, double *gup_tx, double *gup_ty, double *gup_tz,
			  double *gup_xx, double *gup_xy, double *gup_xz, 
			  double *gup_yy, double *gup_yz, double *gup_zz, double *detg,
			  double *GAMMA_in)

{

  FTYPE U[NPR], Prim[NPR];
  FTYPE gdet, gcon[NDIM][NDIM], gcov[NDIM][NDIM];
  int   i,j,k,answer, bad_answer, ret;
  double tilde_ux, tilde_uy, tilde_uz, ut;
  double beta_sq;

  // try introducting local variables.
  double rhol, Pl, vxl, vyl, vzl;
  double rho_sl, s_tl, s_xl, s_yl, s_zl;
  double g_ttl, g_txl, g_tyl, g_tzl, g_xxl, g_xyl, g_xzl, g_yyl, g_yzl, g_zzl;
  double gup_ttl, gup_txl, gup_tyl, gup_tzl, gup_xxl, gup_xyl, gup_xzl, gup_yyl, gup_yzl, gup_zzl;
  double detgl,GAMMA;

  GAMMA = *GAMMA_in;

  rhol = *rho;
  Pl   = *P;
  vxl  = *vx;
  vyl  = *vy;
  vzl  = *vz;
  rho_sl = *rho_s;
  s_tl   = *s_t;
  s_xl   = *s_x;
  s_yl   = *s_y;
  s_zl   = *s_z;

  g_ttl = *g_tt;
  g_txl = *g_tx;
  g_tyl = *g_ty;
  g_tzl = *g_tz;
  g_xxl = *g_xx;
  g_xyl = *g_xy;
  g_xzl = *g_xz;
  g_yyl = *g_yy;
  g_yzl = *g_yz;
  g_zzl = *g_zz;

  gup_ttl = *gup_tt;
  gup_txl = *gup_tx;
  gup_tyl = *gup_ty;
  gup_tzl = *gup_tz;
  gup_xxl = *gup_xx;
  gup_xyl = *gup_xy;
  gup_xzl = *gup_xz;
  gup_yyl = *gup_yy;
  gup_yzl = *gup_yz;
  gup_zzl = *gup_zz;

  detgl = *detg;

  /*********************************************************************
    Assign metric: 
  ***********************************************************************/

/*   printf("inside inversion interface\n");  */
/*   printf("rho = %f\n",rhol);  */
/*   printf("P   = %f\n",Pl);  */
/*   printf("vx  = %f\n",vxl);  */
/*   printf("vy  = %f\n",vyl);  */
/*   printf("vz  = %f\n",vzl); */
/*   printf("rhos= %f\n",rho_sl); */
/*   printf("s_t  = %f\n",s_tl); */
/*   printf("s_x  = %f\n",s_xl); */
/*   printf("s_y  = %f\n",s_yl); */
/*   printf("s_z  = %f\n",s_zl); */

  gdet = *detg;

  gcon[0][0] = gup_ttl;
  gcon[0][1] = gup_txl;
  gcon[0][2] = gup_tyl;
  gcon[0][3] = gup_tzl;
  gcon[1][1] = gup_xxl;
  gcon[1][2] = gup_xyl;
  gcon[1][3] = gup_xzl;
  gcon[2][2] = gup_yyl;
  gcon[2][3] = gup_yzl;
  gcon[3][3] = gup_zzl;

  // symmetry
  gcon[1][0] = gcon[0][1];
  gcon[2][0] = gcon[0][2];
  gcon[3][0] = gcon[0][3];
  gcon[2][1] = gcon[1][2];
  gcon[3][1] = gcon[1][3];
  gcon[3][2] = gcon[2][3];

  gcov[0][0] = g_ttl;
  gcov[0][1] = g_txl;
  gcov[0][2] = g_tyl;
  gcov[0][3] = g_tzl;
  gcov[1][1] = g_xxl;
  gcov[1][2] = g_xyl;
  gcov[1][3] = g_xzl;
  gcov[2][2] = g_yyl;
  gcov[2][3] = g_yzl;
  gcov[3][3] = g_zzl;

  // symmetry
  gcov[1][0] = gcov[0][1];
  gcov[2][0] = gcov[0][2];
  gcov[3][0] = gcov[0][3];
  gcov[2][1] = gcov[1][2];
  gcov[3][1] = gcov[1][3];
  gcov[3][2] = gcov[2][3];
 
  /*********************************************************************
    Assign the primitive variables guess 
  ***********************************************************************/
  // Calculate the \tilde{u}^i from v^i.
  ut = g_ttl + 2.0*g_txl*vxl + 2.0*g_tyl*vyl + 2.0*g_tzl*vzl + 
    g_xxl*vxl*vxl + 2.0*g_xyl*vxl*vyl + 2.0*g_xzl*vxl*vzl + 
    g_yyl*vyl*vyl + 2.0*g_yzl*vyl*vzl + g_zzl*vzl*vzl;
  if (ut >= 0.0) {
    printf("Error superluminal!\n");
    ut = g_ttl;
  }
  ut = sqrt(-1.0/ut);
  
  tilde_ux = ut*(vxl - gup_txl/ gup_ttl);
  tilde_uy = ut*(vyl - gup_tyl/ gup_ttl);
  tilde_uz = ut*(vzl - gup_tzl/ gup_ttl);

  Prim[RHO]    = rhol;
  Prim[UU]     = Pl/(GAMMA-1.0);
  Prim[UTCON1] = tilde_ux;
  Prim[UTCON2] = tilde_uy;
  Prim[UTCON3] = tilde_uz;
  Prim[BCON1]  = 0.0;
  Prim[BCON2]  = 0.0;
  Prim[BCON3]  = 0.0;

  /*********************************************************************
    Set the conserved variables :
  ***********************************************************************/

  U[RHO] = rho_sl;
  U[UU] = s_tl + rho_sl;
  U[UTCON1] = s_xl;
  U[UTCON2] = s_yl;
  U[UTCON3] = s_zl;
  U[BCON1] = 0.0;
  U[BCON2] = 0.0;
  U[BCON3] = 0.0;

  // Call the inversion procedure.

  *retval = Utoprim_2d(    U, gcov, gcon, gdet, Prim, GAMMA); 

  // Assign output values.
  *rho = Prim[RHO];
  *P = (GAMMA-1.0)*Prim[UU];
  tilde_ux = Prim[UTCON1];
  tilde_uy = Prim[UTCON2];
  tilde_uz = Prim[UTCON3];

  // Calculate three velocities.
  beta_sq = g_xxl*gup_txl*gup_txl + 2.0*g_xyl*gup_tyl*gup_txl + 
    2.0*g_xzl*gup_txl*gup_tzl + g_yyl*gup_tyl*gup_tyl + 
    2.0*g_yzl*gup_tyl*gup_tzl + g_zzl*gup_tzl*gup_tzl;

  beta_sq = beta_sq / (gup_ttl*gup_ttl);
  ut = g_xxl*tilde_ux*tilde_ux + 2.0*g_xyl*tilde_uy*tilde_ux + 
    2.0*g_xzl*tilde_ux*tilde_uz + g_yyl*tilde_uy*tilde_uy + 
    2.0*g_yzl*tilde_uy*tilde_uz + g_zzl*tilde_uz*tilde_uz;
  ut = ut + 1.0;
  ut = sqrt(ut/(-g_ttl + beta_sq));
	    
  *vx = tilde_ux/ut + gup_txl/ gup_ttl;
  *vy = tilde_uy/ut + gup_tyl/ gup_ttl;
  *vz = tilde_uz/ut + gup_tzl/ gup_ttl;

  return;
  
}

