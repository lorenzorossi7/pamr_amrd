#ifndef _IO_H
#define _IO_H 1
/*======================================================================================== */
/* io.h                                                                                    */
/*                                                                                         */
/* sdf input/output routines                                                               */
/*======================================================================================== */

void debug_save_gh(char *name, cgh *gh);
void debug_save_cgh(char *name, cgh *gh);
void debug_save_gsl(char *name, real t, gsl *g, int rank);

#define PAMR_CP_DIR_SAVE 0
#define PAMR_CP_DIR_RESTORE 1

int PAMR_do_cp(char *cp_file, int my_rank, int mpi_size, int dir);

#endif /* _IO_H */
