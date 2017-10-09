#ifndef _CLS_H
#define _CLS_H
/*============================================================================= */
/* cls.h --- clustering routines                                                */
/*============================================================================= */

#include "pamr_w.h"
#include "util_w.h"

#define MAX_CLUSTERS 1000
#define CLS_SIMPLE 0
#define CLS_ALIGN_SHRINK 0
#define CLS_ALIGN_EXPAND 1

void simple_cls(int L1, int L2, real **gbbox, int **glev, int *gnum);
// void simple_cls_hydro(int L1, int L2, real **gbbox, int **glev, int *gnum);
// void calc_tre_hydro(int L1, int L2);
void calc_tre(int L1, int L2);
void zero_f_tre(int L);
void adjust_cls(int L1, int L2, real **gbbox, int **glev, int *gnum);

/*============================================================================= */
/* internal routines                                                            */
/*============================================================================= */
void bleed_tre(int max_buf, int min_max);
void bleed_tre_c(int max_buf, int min_max);
void zero_f_tre_ibc(void);
void zero_f_tre_ibc_c(void);
void zero_tre_ibc_a(void);
void zero_tre_ibc_a_c(void);
int find_seed(int *i0, int *j0, int *k0);
int add_cluster(int i0, int j0, int k0, real *bbox);
int find_seed_c(int *i0, int *j0, int *k0);
int add_cluster_c(int i0, int j0, int k0, real *bbox);
void cls_island_info(real *bboxes, int num, int *island_no, int *net_rhosp_n, int L);
void merge_cls(real *bbox, int *num, int lev);
void clean_cls(real *bbox, int *num, int lev, real *pbbox, int pnum);
void extend_ex_cls(real *bbox, int num, int lev, real *bh_bboxes, int num_bhs);

#endif /* _CLS_H */
