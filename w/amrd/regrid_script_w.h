#ifndef _REGRID_SCRIPT_H
#define _REGRID_SCRIPT_H
/*============================================================================= */
/* regrid_script.h --- reading/writing regridding script                        */
/*============================================================================= */

#include <stdio.h>

extern FILE *rgs_stream;

/*============================================================================= */
/* Reads the next bit of regridding info from the specified file.               */
/* Returns 1 if OK.                                                             */
/*                                                                              */
/* returns:                                                                     */
/*                                                                              */
/* time : time of regrid                                                        */
/* L1   : coarsest level overwhich hierarchy can change                         */
/* Lf   : finest level of new hierarhcy                                         */
/* gnum : number of bboxes                                                      */
/* glev : level #'s                                                             */
/* gbbox : bbox's                                                               */
/*                                                                              */
/* caller must free() gnum,glev,gbbox.                                          */
/*============================================================================= */
int rgs_read_next(real *time, int *L1, int *Lf, int *gnum, int **glev, real **gbbox);

/*============================================================================= */
/* The 'inverse' of rgs_read_next()                                             */
/*============================================================================= */
int rgs_write_next(real time, int L1, int Lf, int gnum, int *glev, real *gbbox);

#endif /* _REGRID_SCRIPT_H */
