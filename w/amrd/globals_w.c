//=============================================================================
// globals.c --- global variable declarations ... see globals.h for a
// description of some of these
//=============================================================================

#include "globals_w.h"

const char *line_break="===============================================================================";

int (*app_id)(void)=0;
void (*app_var_pre_init)(char *pfile)=0;
void (*app_var_post_init)(char *pfile)=0;
void (*app_AMRH_var_clear)(void)=0;
void (*app_free_data)(void)=0;
void (*app_t0_cnst_data)(void)=0;
real (*app_evo_residual)(void)=0;
real (*app_MG_residual)(void)=0;
void (*app_evolve)(int iter, int *ifc_mask)=0;
real (*app_MG_relax)(void)=0;
void (*app_L_op)(void)=0;
void (*app_pre_io_calc)(void)=0;
void (*app_scale_tre)(void)=0;
void (*app_post_regrid)(void)=0;
void (*app_post_tstep)(int L)=0;
void (*app_pre_tstep)(int L)=0;
void (*app_fcs_var_clear)(int flag, int *ifc_mask)=0;
void (*app_flux_correct)(void)=0;
void (*app_post_flux_correct)(void)=0;
void (*app_fill_ex_mask)(real *mask, real *mask_c, int dim,
			 int *shape, int *shape_c, real *bbox, real excised);
void (*app_elliptic_vars_t0_init)(void)=0;
void (*app_fill_bh_bboxes)(real *bbox, int *num, int max_num);

void (*app_user_cp)(int save_restore, char *data)=0;

// Branson messing around.
void (*app_post_evolve)(int)=0;

int AMRD_cp_data_size,AMRD_cp_restart;
real AMRD_cp_delta_t_hrs,AMRD_cp_first_delta_t_hrs;
char *AMRD_cp_save_fname;
char *AMRD_cp_restore_fname;

int AMRD_num_hyperbolic_vars;
int AMRD_num_elliptic_vars_t0;
int AMRD_num_elliptic_vars;
int AMRD_num_AMRH_work_vars;
int AMRD_num_MGH_work_vars;
int AMRD_num_AMRH_work_in_MGH_vars;
int AMRD_num_fc_vars;
int AMRD_num_past_bdy_interp_vars;
int AMRD_num_work_repop_vars;
int AMRD_num_AMRH_work_inject_vars;

char **AMRD_hyperbolic_vars;
char **AMRD_elliptic_vars_t0;
char **AMRD_elliptic_vars;
char **AMRD_AMRH_work_vars;
char **AMRD_MGH_work_vars;
char **AMRD_AMRH_work_in_MGH_vars;
char **AMRD_fc_vars;
char **AMRD_past_bdy_interp_vars;
char **AMRD_work_repop_vars;
char **AMRD_AMRH_work_inject_vars;

//Variables that have weighted average injection
int AMRD_num_inject_wavg_vars;
char **AMRD_inject_wavg_vars;
int AMRD_num_MG_inject_wavg_vars;

char *app_name;

int AMRD_dim;        // spatial dimension
int AMRD_num_evo_tl; // number of evolution time levels
int AMRD_ic_n; // which time level to initialize at t=0

int AMRD_steps;     
int AMRD_lsteps[PAMR_MAX_LEVS];     
int AMRD_evo_max_iter;
int AMRD_evo_min_iter;
int AMRD_MG_start_iter;
int AMRD_MG_max_iter;
int AMRD_MG_min_iter;
int AMRD_MG_max_citer;
int AMRD_MG_pre_swp;
int AMRD_MG_pst_swp;
int AMRD_np1_initial_guess;
int AMRD_MG_reinterp_bdy;

int AMRD_base_shape[PAMR_MAX_DIM];
real AMRD_base_bbox[PAMR_MAX_DIM*2];
int AMRD_max_lev;
real AMRD_t0;

real AMRD_evo_tol;  
real AMRD_MG_tol;
real AMRD_MG_crtol;
real AMRD_MG_w0_r,AMRD_MG_w0_i;

int AMRD_MG_extrap_method;
int AMRD_evo_ssc;
real AMRD_MG_eps_c;

int AMRD_id_method;
int AMRD_id_pl_method,AMRD_id_pl_steps,AMRD_id_user_mg_pl;
real AMRD_id_pl_lambda;

real lambda;

real AMRD_TRE_max;
int AMRD_using_cc_tre;
//int AMRD_num_TRE_max_hydro;
//real AMRD_TRE_max_hydro[PAMR_MAX_LEVS];
int AMRD_TRE_norm;
int AMRD_TRE_buffer;
int AMRD_TRE_ibc_buffer;
int AMRD_TRE_ibc_a_buffer;
int AMRD_TRE_exc_buffer;
int AMRD_TRE_exc_buffer_lmin;
int AMRD_TRE_ibcp_buffer;
int AMRD_TRE_sgpbh;
int AMRD_regrid_interval;
int AMRD_regrid_min_lev;
int AMRD_skip_frg;
int can_regrid;

int AMRD_regrid_script;
char *AMRD_regrid_script_name;

char *AMRD_restart_file;
char *AMRD_cp_tag;
real AMRD_cp_times[AMRD_MAX_CP_TIMES];
char *AMRD_save_tag;
char AMRD_curr_cp_tag='A';

int *AMRD_save_ivec[PAMR_MAX_LEVS];
int *AMRD_save_ivec0,AMRD_num_save_mg_vars,AMRD_num_save_n_vars[AMRD_MAX_TIMES];
char *AMRD_vt_times;
char **AMRD_save_n_vars[AMRD_MAX_TIMES];
char **AMRD_save_mg_vars;

int AMRD_save_iter_mod;

int AMRD_save_all_finer;

char AMRD_coord_names[PAMR_MAX_DIM][AMRD_MAX_VNAME_SIZE];

int AMRD_echo_params;
int AMRD_MG_trace;
int AMRD_ID_DV_trace;
int AMRD_MG_DV_trace;
real AMRD_MG_DV_trace_t_on;
real AMRD_MG_DV_trace_t_off;
int AMRD_evo_trace;
int AMRD_evo_DV_trace;
int AMRD_regrid_trace;
real AMRD_evo_DV_trace_t_on;
real AMRD_evo_DV_trace_t_off;

int pamr_context=0;
int my_rank;

int AMRD_cls_method;
int AMRD_cls_merge_dist;
int AMRD_cls_align_mode;

int AMRD_num_tnp1_diss_vars;
int AMRD_num_tn_diss_vars;
char **AMRD_tnp1_diss_vars;
char **AMRD_tn_diss_vars;
int AMRD_tn_diss_pbt[AMRD_MAX_VARS][PAMR_MAX_DIM*2];
int AMRD_tnp1_diss_pbt[AMRD_MAX_VARS][PAMR_MAX_DIM*2];
real AMRD_eps_diss,AMRD_tn_eps_diss,AMRD_tnp1_eps_diss;
int AMRD_repop_diss_bdy;
int AMRD_diss_bdy;
int AMRD_diss_freq;
int AMRD_diss_all_past;
int AMRD_diss_use_6th_order;
int AMRD_diss_stride;
int AMRD_num_rg_diss_vars;
int AMRD_tnum_rg_diss_vars;
char **AMRD_rg_diss_vars;
int AMRD_rg_diss_pbt[AMRD_MAX_VARS][PAMR_MAX_DIM*2];
real AMRD_rg_eps_diss;

// linearly interpolate i-1 physical/amr/both boundaries:
int AMRD_num_tnp1_liipb_vars;
int AMRD_num_tnp1_liiab_vars;
int AMRD_num_tnp1_liibb_vars;
char **AMRD_tnp1_liipb_vars;
char **AMRD_tnp1_liiab_vars;
char **AMRD_tnp1_liibb_vars;

int AMRD_num_interp_AMR_bdy_vars;
int AMRD_interp_AMR_bdy_order;
char **AMRD_interp_AMR_bdy_vars;

int AMRD_calc_global_var_norms;
int AMRD_global_var_norm_type;
real AMRD_global_var_norms[PAMR_MAX_GFNS];
real AMRD_global_var_norm_floor;

real AMRD_ex;
int AMRD_do_ex;

char *AMRD_c_tag, *AMRD_v_tag;

int AMRD_AMR_bdy_width, AMRD_AMR_bdy_width_c;

int AMRD_resid_in_evo;

//=============================================================================
// "internal" globals 
//=============================================================================
int AMRD_MG_res_gfn[AMRD_MAX_VARS],AMRD_MG_lop_gfn[AMRD_MAX_VARS],AMRD_MG_brs_gfn[AMRD_MAX_VARS],
    AMRD_MG_rhs_gfn[AMRD_MAX_VARS],AMRD_MG_f_gfn[AMRD_MAX_VARS],AMRD_fcs_gfn[AMRD_MAX_VARS];
real *AMRD_MG_res[AMRD_MAX_VARS],*AMRD_MG_lop[AMRD_MAX_VARS],*AMRD_MG_brs[AMRD_MAX_VARS],
     *AMRD_MG_rhs[AMRD_MAX_VARS],*AMRD_MG_f[AMRD_MAX_VARS];
int AMRD_cmask_gfn,AMRD_cmask_mg_gfn,AMRD_tre_gfn,AMRD_mg_tre_gfn,AMRD_tre_c_gfn;
real *AMRD_cmask,*AMRD_cmask_mg,*AMRD_tre,*AMRD_mg_tre,*AMRD_tre_c;

int AMRD_cmask_c_gfn;  // I don't think we'll need the cell centered cmask for any
real *AMRD_cmask_c;    // multigrid stuff.

// The flux correction mask
int AMRD_fc_mask_gfn;
real *AMRD_fc_mask;

int AMRD_chr_gfn,AMRD_chr_mg_gfn;
real *AMRD_chr,*AMRD_chr_mg;

int AMRD_chr_c_gfn,AMRD_chr_mg_c_gfn;
real *AMRD_chr_c,*AMRD_chr_mg_c;

//For weighted averaging (e.g. by volume)
int AMRD_wavg_c_gfn,AMRD_wavg_mg_c_gfn;

int AMRD_AMR_f_gfn[AMRD_MAX_TIMES][AMRD_MAX_VARS];
real *AMRD_AMR_f[AMRD_MAX_TIMES][AMRD_MAX_VARS];

int AMRD_AMR_mgf_gfn[AMRD_MAX_TIMES][AMRD_MAX_VARS];
real *AMRD_AMR_mgf[AMRD_MAX_TIMES][AMRD_MAX_VARS];

int AMRD_AMR_mgf_extrap_tm1_gfn[AMRD_MAX_VARS];
real *AMRD_AMR_mgf_extrap_tm1[AMRD_MAX_VARS];

int AMRD_AMR_mgf_extrap_tm2_gfn[AMRD_MAX_VARS];
real *AMRD_AMR_mgf_extrap_tm2[AMRD_MAX_VARS];

int AMRD_work_repop_gfn[AMRD_MAX_VARS];
real *AMRD_work_repop[AMRD_MAX_VARS];

int AMRD_AMRH_work_inject_gfn[AMRD_MAX_VARS];
int AMRD_AMRH_work_inject_op[AMRD_MAX_VARS];

int AMRD_inject_wavg_gfn[AMRD_MAX_VARS];
int AMRD_MG_inject_wavg_gfn[AMRD_MAX_VARS];

int AMRD_past_bdy_interp_vars_gfn[AMRD_MAX_TIMES][AMRD_MAX_VARS];
int AMRD_past_bdy_interp_op[AMRD_MAX_VARS];

int AMRD_tnp1_diss_gfn[AMRD_MAX_VARS];
real *AMRD_tnp1_diss[AMRD_MAX_VARS];
int AMRD_tn_diss_gfn[AMRD_MAX_VARS];
real *AMRD_tn_diss[AMRD_MAX_VARS];

int AMRD_ex_repop1_sgfn[AMRD_MAX_VARS];
int AMRD_ex_repop2_sgfn[AMRD_MAX_VARS];
int AMRD_ex_repop3_sgfn[AMRD_MAX_VARS];
int AMRD_ex_repop4_sgfn[AMRD_MAX_VARS];

int AMRD_tnp1_liipb_gfn[AMRD_MAX_VARS];
int AMRD_tnp1_liiab_gfn[AMRD_MAX_VARS];
int AMRD_tnp1_liibb_gfn[AMRD_MAX_VARS];
real *AMRD_tnp1_liipb[AMRD_MAX_VARS];
real *AMRD_tnp1_liiab[AMRD_MAX_VARS];
real *AMRD_tnp1_liibb[AMRD_MAX_VARS];

int AMRD_rg_diss_gfn[AMRD_MAX_VARS*AMRD_MAX_TIMES];
real *AMRD_rg_diss[AMRD_MAX_VARS*AMRD_MAX_TIMES];

int AMRD_interp_AMR_bdy_f_gfn[AMRD_MAX_VARS];
real *AMRD_interp_AMR_bdy_f[AMRD_MAX_VARS];

int AMRD_num_f_tre_vars;
char *AMRD_f_tre_vars[AMRD_MAX_VARS];
int AMRD_f_tre_gfn[AMRD_MAX_VARS],AMRD_f_tre_fn[AMRD_MAX_VARS];
real *AMRD_f_tre[AMRD_MAX_VARS];

real *AMRD_g_x[PAMR_MAX_DIM],*AMRD_g_x_c[PAMR_MAX_DIM],AMRD_g_dx[PAMR_MAX_DIM];
int AMRD_g_shape[PAMR_MAX_DIM],AMRD_g_shape_c[PAMR_MAX_DIM],AMRD_g_ghost_width[2*PAMR_MAX_DIM],AMRD_g_rank,AMRD_g_dim,AMRD_g_size,AMRD_g_size_c,AMRD_g_Nx,AMRD_g_Ny,AMRD_g_Nz;
int AMRD_g_Nx_c,AMRD_g_Ny_c,AMRD_g_Nz_c;
int AMRD_g_phys_bdy[2*PAMR_MAX_DIM];
real AMRD_g_bbox[2*PAMR_MAX_DIM];

int AMRD_periodic[PAMR_MAX_DIM];

int AMRD_w1_gfn;
real *AMRD_w1;
int AMRD_w1_c_gfn;
real *AMRD_w1_c;

int AMRD_tre_valid[PAMR_MAX_LEVS];

int AMRD_num_MG_cnst_data_vars;
int AMRD_MG_cnst_data_gfn[AMRD_MAX_VARS];
real *AMRD_MG_cnst_data[AMRD_MAX_VARS];
char **AMRD_MG_cnst_data_vars;

int AMRD_num_free_data_vars;
int AMRD_free_data_gfn[AMRD_MAX_VARS];
real *AMRD_free_data[AMRD_MAX_VARS];
char **AMRD_free_data_vars;

int AMRD_elliptic_vars_t0_gfn[AMRD_MAX_VARS];

int AMRD_is_t0=1,AMRD_state;

int AMRD_num_ex_repop1_vars;
char **AMRD_ex_repop1_vars;

int AMRD_num_ex_repop2_vars;
char **AMRD_ex_repop2_vars;

int AMRD_num_ex_repop3_vars;
char **AMRD_ex_repop3_vars;

int AMRD_num_ex_repop4_vars;
char **AMRD_ex_repop4_vars;

int AMRD_max_t_interp_order;
real *AMRD_t_interp_substeps;
int AMRD_num_t_interp_substeps;

int AMRD_re_interp_width;
int AMRD_re_interp_width_c;

int AMRD_c_to_v_in_tstep;
int AMRD_v_to_c_in_tstep;

int AMRD_curr_MG_iter=0;

int AMRD_magic_cookie=0;
