#ifndef _MG_H
#define _MG_H
/*============================================================================= */
/* mg.h --- multigrid routines                                                  */
/*============================================================================= */

/*============================================================================= */
/* internal routines                                                            */
/*============================================================================= */
void mg_ldptr();
void clear_MG_vars(int L);
void set_cmask_child(int L, int hier);
void set_cmask_bdy(int L, int hier);
void set_cmask_bdy_w1(int L, int hier);
void extend_cmask_bdy(int L, int hier, int num);
int is_coarsest(int L);
void set_comm_coarsest(int L);
void reset_comm(int L);
void compute_MG_rhs(int L, int save_tre);
void compute_cgc(int L);
void set_res_sync(void);
void mg_dump(int L, char *tag, int giter, int which);
void sync_cnst_data(int L);
void interp_cnst_data(int Lc);
void inject_cnst_data(int Lf);

/*============================================================================= */
/* solves the set of elliptic equations from levels L..Lf at time level tl      */
/*============================================================================= */
void solve_elliptics(int L, int tl, int save_tre);

/*============================================================================= */
/* 1 vcycle ... returns the largest residual on the finest level                */
/*============================================================================= */
real vcycle(int giter, int save_tre);

#endif /* _MG_H */
