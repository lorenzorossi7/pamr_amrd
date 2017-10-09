#ifndef _UTIL_H
#define _UTIL_H
#include "mpi.h"
#include "amrd.h"
/*============================================================================= */
/* util.h --- misc. utility routines                                            */
/*                                                                              */
/* NOTE ... some of the prototypes for functions in util.c are in               */
/* amrd.h, so that applications can use them                                    */
/*============================================================================= */

#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))

#define IFL if(ltrace)
#define IFL0 if(ltrace && my_rank==0)
/*============================================================================= */
/* clean up stuff before we stop(), or end normally                             */
/*============================================================================= */
void clean_up(void);

/*============================================================================= */
/* error handler                                                                */
/*============================================================================= */
void AMRD_stop(char *err1, char *err2);

/*============================================================================= */
/* misc. routines                                                               */
/*============================================================================= */
void call_app(void (*f)(void), int L, int hier);
void call_app_1int(void (*f)(int arg1), int arg1, int L, int hier);
real call_rr_app(real (*f)(void), int L, int hier, int comm_only, MPI_Op op);

real approx_l2norm(int gfn, int L, int hier);

void debug_save_bbox_list(real *bbox, int *lev, int glev, int num, real t, char *tag);
void gh_stats(void);
int total_mem(void);

#endif /* _UTIL_H */
