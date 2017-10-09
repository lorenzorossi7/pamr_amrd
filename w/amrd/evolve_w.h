#ifndef _EVOLVE_H
#define _EVOLVE_H
/*============================================================================= */
/* evolve.h --- evolution routines                                              */
/*============================================================================= */

/*============================================================================= */
/* global variables that need to be saved for check-pointing                    */
/*============================================================================= */
extern int AMRD_stack[2*PAMR_MAX_LEVS],AMRD_top;
extern real rgs_t;
extern int rgs_L1;
extern int rgs_Lf;

void evolve(int coarse_steps, int ic);
void regrid(int lev, int resolve, int trace);
void tstep(int lev, int ic);
void tstep_hydro(int lev, int ic);
void compute_ssh_tre(int lev);
void compute_ssh_tre_hydro(int lev);
void set_mg_brs_vars(int L1, int L2);
void update_mg_extrap_vars(int L1, int Lf);
void set_cmask_oldgs(real *bbox, int num);
void generate_id(void);
void init_rest(void);
void id_pl_2_save(void);
void id_pl_2_init(real tx);
void flip_dt(void);
void set_gfn_sync(int gfn);
void set_gfn_in(int gfn, int in);
void set_regrid_sync();
void ev_ldptr(void);
void evo_dump(int L1, int L2, char *tag, int iter, int which);
void interp_AMR_bdy(real* f, real *work, int rhosp);
void calc_global_norms();
void ibnd_interp(int L);
void set_amr_bdy_2nd_ord(int lev, real tnp1, real dt, real tnp1_lm1, real dt_lm1, real frac);
void set_amr_bdy_2nd_ord_hydro(int lev, real tnp1, real dt, real tnp1_lm1, real dt_lm1, real frac);
void set_amr_bdy_3rd_ord(int lev, real tnp1, real dt, real tnp1_lm1, real dt_lm1, real frac);
int *get_ifc_mask(real *fc_mask, int size);

#endif /* _EVOLVE_H */
