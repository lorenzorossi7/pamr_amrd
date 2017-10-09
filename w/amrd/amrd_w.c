//=============================================================================
// AMR driver program. Copyright 2002 F.Pretorius.
//
// "main" module.
//=============================================================================

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "globals_w.h"
#include "io_w.h"
#include "util_w.h"
#include "evolve_w.h"

void amrd(int argc, char **argv,
          int (*app_id_f)(void),
          void (*app_var_pre_init_f)(char *pfile),
          void (*app_var_post_init_f)(char *pfile),
          void (*app_AMRH_var_clear_f)(void),
          void (*app_free_data_f)(void),
          void (*app_t0_cnst_data_f)(void),
          real (*app_evo_residual_f)(void),
          real (*app_MG_residual_f)(void),
          void (*app_evolve_f)(int iter, int *ifc_mask),
          real (*app_MG_relax_f)(void),
          void (*app_L_op_f)(void),
          void (*app_pre_io_calc_f)(void),
          void (*app_scale_tre_f)(void),
          void (*app_post_regrid_f)(void),
          void (*app_post_tstep_f)(int L),
          void (*app_fill_ex_mask_f)(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised),
          void (*app_fill_bh_bboxes_f)(real *bbox, int *num, int max_num))
{
   int size,Lf,i,used_app_id;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&size);

   if (argc!=2)
   {
      printf("usage: %s param_file\n",argv[0]);
      goto end;
   }

   app_id=app_id_f;
   app_var_pre_init=app_var_pre_init_f;
   app_var_post_init=app_var_post_init_f;
   app_AMRH_var_clear=app_AMRH_var_clear_f;
   app_free_data=app_free_data_f;
   app_t0_cnst_data=app_t0_cnst_data_f;
   app_evo_residual=app_evo_residual_f;
   app_MG_residual=app_MG_residual_f;
   app_evolve=app_evolve_f;
   app_MG_relax=app_MG_relax_f;
   app_L_op=app_L_op_f;
   app_pre_io_calc=app_pre_io_calc_f;
   app_scale_tre=app_scale_tre_f;
   app_post_regrid=app_post_regrid_f;
   app_post_tstep=app_post_tstep_f;
   app_fill_ex_mask=app_fill_ex_mask_f;
   app_fill_bh_bboxes=app_fill_bh_bboxes_f;

   app_name=0;
   AMRD_str_param(argv[1],"app_name",&app_name,1);

   if (my_rank==0) 
      printf("\n%s\nAMR driver for %s, on %i nodes\n%s\n\n",line_break,app_name,size,line_break);

   AMRD_state=AMRD_STATE_PRE_INIT;
   app_var_pre_init(argv[1]);

   init_context(argv[1]);

   AMRD_state=AMRD_STATE_GENERATE_ID;
   app_var_post_init(argv[1]);

   AMRD_state=AMRD_STATE_GENERATE_ID;
   if (!AMRD_cp_restart && !(app_id())) generate_id();
   else if (AMRD_cp_restart)
   {
      Lf=PAMR_get_max_lev(PAMR_AMRH);
      for (i=1; i<=Lf; i++) {
          if (i>1) set_fc_mask(i);
          call_app(app_post_regrid,i,PAMR_AMRH); // since restarting involves a regrid
      }
   }

   AMRD_state=AMRD_STATE_EVOLVE;
   if (AMRD_steps>0) evolve(AMRD_steps,0);

   if (my_rank==0) printf("\nEvolution complete\n\n");
   
end:
   if (pamr_context) PAMR_free_context(pamr_context);
   clean_up();
   MPI_Finalize();
}

/*=============================================================================*/
/* The following function defines a post-regrid hook                           */
/*=============================================================================*/
void amrd_set_app_pre_tstep_hook(void (*app_pre_tstep_f)(int L))
{
   app_pre_tstep=app_pre_tstep_f;
}

/*=============================================================================*/
/* The following function defines the hook function for initializing elliptic  */
/* t0 variables                                                                */
/*=============================================================================*/
void amrd_set_elliptic_vars_t0_init(void (*app_elliptic_vars_t0_init_f)(void))
{
   app_elliptic_vars_t0_init=app_elliptic_vars_t0_init_f;
}

/*=============================================================================*/
/* The following function defines the check-point hook                         */
/*=============================================================================*/
void amrd_set_app_user_cp_hook(void (*app_pre_tstep_f)(int save_restore, char *data), int cp_data_size)
{
   app_user_cp=app_pre_tstep_f;
   if (cp_data_size<0) AMRD_stop("amrd_set_app_user_cp_hook: error ... cp_data_size must be >=0\n","");
   AMRD_cp_data_size=cp_data_size;
}

/*=============================================================================*/
/* The following hook functions are for hydro                                  */
/*=============================================================================*/
void amrd_set_app_fcs_var_clear_hook(void (*app_fcs_var_clear_f)(int, int *))
{
  app_fcs_var_clear=app_fcs_var_clear_f;
}

void amrd_set_app_flux_correct_hook(void (*app_flux_correct_f)())
{
   app_flux_correct=app_flux_correct_f;
}

void amrd_set_app_post_flux_correct_hook(void (*app_post_flux_correct_f)())
{
   app_post_flux_correct=app_post_flux_correct_f;
}

