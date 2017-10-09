#ifndef _IO_H
#define _IO_H
/*============================================================================= */
/* io.h --- reading/writing parameters, files, etc.                             */
/*============================================================================= */

/*============================================================================= */
/* The following reads parameters associated with amrd, from pfile. At the      */
/* same time, the context is initialized.                                       */
/*============================================================================= */
void init_context(char *pfile);

/*============================================================================= */
/* various save-to-disk functions (some such debug related functions            */
/* are in evolve.c, mg.c, ...)                                                  */
/*============================================================================= */
void save0(int iter);
void saveL(int lev, int iter);
int amrd_do_cp(char *AMRD_cp_restore_fname, int dir);
#endif /* _IO_H */
