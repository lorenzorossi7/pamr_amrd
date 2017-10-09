#ifndef _PAMR_H
#define _PAMR_H
/*============================================================================= */
/*                                                                              */
/*  pamr.h --- parallel AMR library                                             */
/*                                                                              */
/*  Copyright 2002/3/4/5 F.Pretorius                                              */
/*                                                                              */
/*============================================================================= */

#include "internal_opts.h"

/*====================================================================================================== */
/* API --- see pamr.c for some comments on arguments                                                     */
/*====================================================================================================== */
int PAMR_init_context(real *q, int q_size, int dim, int *shape, real *bbox, char *cp_file, int cp_rank);
int PAMR_compose_hierarchy(int min_lev, int max_lev, int num, int *lev, real *bbox, real t);
int PAMR_build_mgh(int min_lev, int max_lev, int tl);
int PAMR_destroy_mgh();
void PAMR_free_context(int num);
int PAMR_merge_bboxes(real *local_bbox, int local_num, real *global_bbox, int *global_num, real min_eff);
int PAMR_get_sgh(real *bbox, int *num, int L);

void PAMR_set_lambda(real lambda);
void PAMR_get_lambda(real *lambda);
void PAMR_set_rho(int *rho_sp, int *rho_tm, int num);
void PAMR_get_rho(int *rho_sp, int *rho_tm, int num);
void PAMR_get_dxdt(int lev, real *dx, real *dt);
void PAMR_set_top_ratios(real top_xy, real top_xz, real top_yz);
void PAMR_get_top_ratios(real *top_xy, real *top_xz, real *top_yz);
void PAMR_set_ghost_width(int *ghost_width);
void PAMR_get_ghost_width(int *ghost_width);
void PAMR_set_min_width(int *min_width);
void PAMR_get_min_width(int *min_width);
void PAMR_set_MG_coarse_width(int *min_width);
void PAMR_get_MG_coarse_width(int *min_width);
void PAMR_set_interp_buffer(int interp_buffer);
void PAMR_get_interp_buffer(int *interp_buffer);
void PAMR_get_global_bbox(real *bbox);
void PAMR_set_periodic_bdy(int *periodic);
void PAMR_get_periodic_bdy(int *periodic);

/*============================================================================= */
/* current grid distribution methods                                            */
/* The following #define parameters are switches, where the default is 'off'.   */
/* So set "method" to the bitwise OR (|) of all desired on-flags.               */
/*                                                                              */
/* GDM_GRID_BY_GRID : distribute across network grid-by-grid, as opposed to     */
/*                    level-by-level (default)                                  */
/* GDM_ALIGN        : forces *all* grid boundaries to be aligned with           */
/*                    a parent grid, and does so by increasing the              */
/*                    ghostwidth, if needed                                     */
/* GDM_NO_OVERLAP   : clips all overlapping grids in the sequential hierarchy   */
/*                    prior to distribution (i.e., the only overlap after       */
/*                    compose() will be the ghostzones).                        */
/*                    At this stage this is the 'safest' option for in situations*/
/*                    where there is grid overlap (in terms of guaranteeing that*/
/*                    interior, overlap boundaries are treated correctly as     */
/*                    interior points, and not AMR boundaries).                 */
/*                                                                              */
/* rule of thumb: use level-by-level if (#nodes) >> (typical # grids/level)     */
/*                else use GDM_GRID_BY_GRID                                     */
/*============================================================================= */
#define PAMR_GDM_GRID_BY_GRID 1
#define PAMR_GDM_ALIGN 2
#define PAMR_GDM_NO_OVERLAP 4

void PAMR_set_gdm(int method);
void PAMR_get_gdm(int *method);

/*============================================================================= */
/* variable definition                                                          */
/*============================================================================= */
int PAMR_def_var_brief(char *name);
int PAMR_def_var_full(char *name, int in_amrh, int in_mgh, int num_tl, int amr_inject,
    int amr_interp, int amr_bdy_interp, int amr_sync, int mg_inject, int mg_interp, 
    int mg_sync, int mg_noinj_to_amr, int regrid_transfer, int *phys_bdy_type);
int PAMR_set_var_attribs(char *name, int in_amrh, int in_mgh, int num_tl, int amr_inject,
    int amr_interp, int amr_bdy_interp, int amr_sync, int mg_inject, int mg_interp, 
    int mg_sync, int mg_noinj_to_amr, int regrid_transfer, int *phys_bdy_type);
int PAMR_get_var_attribs(char *name, int *in_amrh, int *in_mgh, int *num_tl, int *amr_inject,
    int *amr_interp, int *amr_bdy_interp, int *amr_sync, int *mg_inject, int *mg_interp, 
    int *mg_sync, int *mg_noinj_to_amr, int *regrid_transfer, int *phys_bdy_type);
void PAMR_disable_var(char *name);
void PAMR_enable_var(char *name);
int PAMR_is_var_enabled(char *name);

void PAMR_set_trace_lev(int lev);
/*============================================================================= */
/* allowed variable attributes                                                  */
/*============================================================================= */
/* injection transfer bits */
#define PAMR_NO_INJECT 0
#define PAMR_STRAIGHT_INJECT 1
#define PAMR_HW_RESTR 2
#define PAMR_FW_RESTR 3

/* interpolation transfer bits.  */
/* Also used to specify the order of interpolation to use during a  */
/* transfer operation. */
#define PAMR_NO_INTERP 0
#define PAMR_SECOND_ORDER 2
#define PAMR_FOURTH_ORDER 4

/* post-regrid transfer bits; for values > 1 uses interpolation types */
#define PAMR_NO_TRANSFER 0

/* synchronization transfer bits (i.e. within a level during an iteration) */
#define PAMR_NO_SYNC 0
#define PAMR_SYNC 1

/* physical boundary types */
#define PAMR_UNKNOWN 0
#define PAMR_EVEN 1
#define PAMR_ODD 2

/* hierarchy indentifiers */
#define PAMR_AMRH 1
#define PAMR_MGH 2

/*============================================================================= */
/* transfer bit manipulation                                                    */
/*============================================================================= */
void PAMR_freeze_tf_bits(void);
void PAMR_thaw_tf_bits(void);
void PAMR_clear_tf_bits(void);
void PAMR_set_tf_bit(int gf, int val);
void PAMR_set_tf_bits(int tl, int hierarchy, int stage);

/*'stage' variables for set_gf_transfer_bits() */
#define PAMR_TF_SYNC 0
#define PAMR_TF_COMPOSE 1
#define PAMR_TF_INJECT 2
#define PAMR_TF_INJECT_TO_MG_LEV 3
#define PAMR_TF_INTERP 4
#define PAMR_TF_BDY_INTERP 5
#define PAMR_TF_MGH_INIT 6
#define PAMR_TF_C_TO_V 7
#define PAMR_TF_V_TO_C 8

/*============================================================================= */
/* communication functions (in addition to compose_hierarchy() )                */
/*============================================================================= */
int PAMR_sync(int l, int tl, int hierarchy, int AMR_bdy_width);
int PAMR_inject(int lf, int tl, int hierarchy);
int PAMR_interp(int lc, int tl, int hierarchy);
int PAMR_AMR_bdy_interp(int lc, int tl, int AMR_bdy_width);

/*============================================================================= */
/* time functions                                                               */
/*============================================================================= */
real PAMR_tick(int lev);
void PAMR_set_time(int lev, real t);
real PAMR_get_time(int lev);
int PAMR_swap_tl(int lev, int tl1, int tl2);

/*============================================================================= */
/* excision support functions.                                                  */
/*============================================================================= */
int PAMR_excision_on(char *ex_mask_var, 
                     void (*app_fill_ex_mask_fnc)(real *mask, int dim, int *shape, real *bbox, real excised),
                     real excised, int initialize_now);
void PAMR_excision_off();

/*============================================================================= 
/* Accessing grids/grid-functions in the hierarchy                              */
/* All functions, except get_gf_x() and get_gf_gfs() expect the user to         */
/* supply the storage for the returned quantity. For get_gf_x() and             */
/* get_gf_gfs(), the user supplies an array of pointers to the relevant         */
/* objects, and the routines fill in the pointers.                              */
/*============================================================================= */
int PAMR_init_s_iter(int l, int hier, int all);
int PAMR_next_g();
int PAMR_push_iter();
int PAMR_pop_iter();
int PAMR_get_g_attribs(int *rank, int *dim, int *shape, real *bbox, 
                  int *ghost_width, real *t, int *ngfs, real **x, real **gfs);
int PAMR_get_g_rank(int *rank);
int PAMR_get_g_dim(int *dim);
int PAMR_get_g_shape(int *shape);
int PAMR_get_g_bbox(real *bbox);
int PAMR_get_g_ghost_width(int *ghost_width);
int PAMR_get_g_t(real *t);
int PAMR_get_g_ngfs(int *ngfs);
int PAMR_get_g_x(real **x);
int PAMR_get_g_gfs(real **gfs);
int PAMR_get_g_coarsest(int *coarsest);
int PAMR_get_gfn(char *name, int hier, int tl);
int PAMR_get_g_level(int *L);

int PAMR_get_max_lev(int hier);
int PAMR_get_min_lev(int hier);

int PAMR_get_g_comm(int *comm);
int PAMR_set_g_comm(int comm);

/*============================================================================= */
/* utility io functions                                                         */
/*                                                                              */
/* names are <pre_tag><name>[_tl#|_MG]<post_tag>_<rank>.sdf                     */
/*============================================================================= */
int PAMR_save_gfn(char *name, int hier, int tl, int L, real t, char *pre_tag, char *post_tag);
int PAMR_save_gfn_close(char *name, int hier, int tl, int L, real t, char *pre_tag, char *post_tag);

/*============================================================================= */
/* Check-points grid-hierarchy to file with tag cp_name.                        */
/*============================================================================= */
int PAMR_cp(char *cp_name, int cp_rank);

/*============================================================================= */
/* Stats                                                                        */
/*============================================================================= */
void PAMR_collect_stats(int on_off);
void PAMR_get_stats(real *comm_r_microsecs, real *comm_s_microsecs,
		    real *comm_r_bytes, real *comm_s_bytes,
		    int *num_r, int *num_s);

#endif /* _PAMR_H */
