#ifndef _TRANSFER_H
#define _TRANSFER_H
/*============================================================================= */
/*                                                                              */
/*  transfer.h --- data structures and functions related to communicating       */
/*  pieces of grids between nodes                                               */
/*                                                                              */
/*============================================================================= */

#include "gh.h"

/*============================================================================= */
/* grid-segment-list structure                                                  */
/*============================================================================= */
typedef struct gsl
{
   struct gsl *next;
   real bbox[2*PAMR_MAX_DIM]; /* bounding box of segment */
   grid *g;                   /* the grid to which this segment refers */
} gsl;

void free_gsl(gsl *p);
gsl *clone_gsl(gsl *p,int first_only);
gsl *alloc_gsl(grid *g, int dst_lev, int interior, int AMR_reduce);
gsl *cat_gsls(gsl *A, gsl *B);
gsl *build_owned_gsl(int rank, cgh *gh, int lev, int interior, int extend);
gsl *build_complete_gsl(cgh *gh, int lev, int AMR_bdy_only);
gsl *gs_subtract(gsl *A, gsl *B);
gsl *gsl_subtract(gsl *A, gsl *B);
int gs_and(gsl *A, gsl *B, gsl **A0, gsl **B0);
int build_gstl(gsl *in_src,gsl *in_dst, gsl **out_src,gsl **out_dst);
void transfer(gsl **src, gsl **dst);
int num_transfer_bits(void);
void local_transfer(gsl *src, gsl *dst);
void free_src_data(gsl *src);
void calc_gs_ibbox(gsl *gs, int *ibbox, real *dx);
int data_packer(real *data, gsl *src, gsl *dst, int rank, int dir);
#define PACK 1
#define UNPACK 2
void MPI_send_gsl(gsl *src, gsl *dst);
int MPI_receive_gsl(gsl **src, gsl **dst);

#endif /* _TRANSFER_H */
