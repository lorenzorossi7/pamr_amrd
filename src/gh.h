#ifndef _GH_H
#define _GH_H
/*============================================================================= */
/*                                                                              */
/*  gh.h --- grid hierarhcy structure definition and functions                  */
/*                                                                              */
/*============================================================================= */

#include "internal_opts.h"

#define ITER_G_STACKSIZE 10

/*============================================================================= */
/* 'var' is a variable-structure                                                */
/*                                                                              */
/* NOTE: the flags amr_inject ... regrid_transfer are ONLY used when            */
/*       the automatic transfer-bit set mechanism is used. The transfer routines*/
/*       use the transfer bits to determine what kind of                        */
/*       interpolation/restriction/etc operator to use, and hence the           */
/*       interpretation of a transfer bit is operation dependent. (This         */
/*       is because at that level the transfer operations are not aware of      */
/*       AMR or MG hierarchies, or what particular function initiation the      */
/*       copy operation).                                                       */
/*                                                                              */
/* regarding status:                                                            */
/* VAR_STATUS_ON  : variable is allocated                                       */
/* VAR_STATUS_OFF : variable is not allocated                                   */
/* VAR_STATUS_TURN_ON : variable is currently not allocated, but will be        */
/*                      during next compose                                     */
/* VAR_STATUS_TURN_OFF: variable is currently allocated, but will not be        */
/*                      during next compose                                     */
/*============================================================================= */
#define VAR_STATUS_ON 1
#define VAR_STATUS_OFF 2
#define VAR_STATUS_TURN_ON 3
#define VAR_STATUS_TURN_OFF 4
typedef struct var
{
   struct var *next;
   char *name;     
   int in_amrh;         /* whether the variable exists in the AMR hierarchy or not */
   int in_mgh;          /*                        "           MG       " */
   int status;          /* one of STATUS_... flages above */
   int num_tl;          /* number of time-levels (in AMR hierarchy), from 1 .. num_tl */
   int sgfn;            /* starting grid-function number */
   int amr_inject;      /* injection operator to use within AMR hiearchy */
   int amr_interp;      /*  interpolation " */
   int amr_bdy_interp;  /*  interpolation " when setting AMR boundary conditions */
   int amr_sync;        /*  synchronization  " */
   int mg_inject;       /* Same flags for variables in the AMR hierarchy, except */
   int mg_interp;       /*  can turn off (with mg_noinj_to_amr=1) injection into */
   int mg_sync;         /*  AMR levels (i.e. so that source functions from the */
   int mg_noinj_to_amr; /*  evolution equations are used verbatim when possible) */
   int regrid_transfer; /* wether to transfer after a regrid */
   int phys_bdy_type[2*PAMR_MAX_DIM]; /* to guide the interpolation */
} var;

/* the following globals are used by set_var_attribs() and define_var_brief() */

extern int c_in_amrh;
extern int c_in_mgh;
extern int c_num_tl;
extern int c_amr_inject;
extern int c_amr_interp;
extern int c_amr_bdy_interp;
extern int c_amr_sync;
extern int c_mg_inject;
extern int c_mg_interp;
extern int c_mg_sync;
extern int c_mg_noinj_to_amr;
extern int c_regrid_transfer;
extern int c_phys_bdy_type[2*PAMR_MAX_DIM];

/*============================================================================= */
/* grid structure ... some info is redundent, such as dim & t, but adding       */
/* here for convenience                                                         */
/*============================================================================= */
typedef struct grid
{
   struct grid *next; 
   struct grid *prev; 
   int dim;
   int shape[PAMR_MAX_DIM];
   real bbox[2*PAMR_MAX_DIM];
   int ghost_width[2*PAMR_MAX_DIM]; /* ghostzone width at each boundary */
   int wrap[2*PAMR_MAX_DIM];      /* wether the corresponding boundary wraps around a periodic dimension */
   real t;
   real *x[PAMR_MAX_DIM];         /* pointers to perimeter coordinate arrays */
   int ngfs;                      /* number of grid functions */
   real *(*gfs);                  /* pointer to an array of pointers to grid function data; */
   int rank;                      /* MPI rank of processor containing this grid */
   int coarsest;                  /* flag denoting wether this is a coarsest "island" for MG */
   int comm;                      /* flag denoting wether this grid partakes in communcations */
   int virtual_g;                 /* a 'virtual' grid ... for periodic boundaries */
} grid;

/*============================================================================= */
/* the level structure, used within a cgh                                       */
/*============================================================================= */
typedef struct level
{
   real t,dt;              /* time info */
   real dx[PAMR_MAX_DIM];  /* spatial discretization info */
   grid *grids;            /* pointer to the list of grids at this level */
} level;

/*============================================================================= */
/* computational grid hierarchy structure                                       */
/*============================================================================= */
typedef struct cgh
{
   int min_lev;                   /* levels[min_lev..max_lev] are used.   */
   int max_lev;
   level *levels[PAMR_MAX_LEVS];    
   int in_amrh[PAMR_MAX_LEVS];    /* for the mgh, to keep track of which levels are shared with the amrh */
} cgh;

/*============================================================================= */
/* A 'context' is an instance of a distributed, adaptive grid hiearhcy          */
/*                                                                              */
/* NOTE: that rho_sp is same for all dimensions is hard-coded in several        */
/* routines. Also, grid widths of 1 will cause problems with some of the        */
/* communication routines (>1 is assumed in most cases)                         */
/*============================================================================= */
typedef struct context
{
   real *q;                         /* if non-zero, then use q[q_size] block for gf storage   */
   int q_size;

   /* the following define the structure of the grid-hierarchy */
   
   int dim;                         /* spatial dimension */
   int rho_sp[PAMR_MAX_LEVS];       /* spatial refinement ratio --- same for all dimensions  */
   int rho_tm[PAMR_MAX_LEVS];       /* temporal refinement ratio */
   real dx[PAMR_MAX_LEVS][PAMR_MAX_DIM]; /* level defined by discretization scale dx of first coordinate */
   real dt[PAMR_MAX_LEVS];      
   real shape[PAMR_MAX_DIM];        /* geometry of base level   */
   real bbox[2*PAMR_MAX_DIM];      
   real *sgh_bboxes[PAMR_MAX_LEVS]; /* a copy of the sequential, AMR grid-hierarchy */
   int num_sgh_bboxes[PAMR_MAX_LEVS];
   real lambda;                     /* courant factor */
   int num_vars;                    /* total number of variables */
   var *vars;                       /* pointer to a linked list of (num_vars) variable structures */
   int gfns;                        /* total number of grid functions */
   int *tf_bits;                    /* pointer to an array of transfer bits, one for each grid function */
   int frozen;                      /* whether the transfer bits are frozen or not. */

   /* the current computational-grid hierachies */
   
   cgh *curr_cgh;
   cgh *MG_cgh;

   /* parallel options */

   int rank;                        /* MPI rank */
   int size;                        /* MPI size */
   int n_rank;                      /* the next node to which data should be allocated  */
   real top_xy;                     /* desired topology ratios for splitting grids */
   real top_xz;
   real top_yz; 
   int ghost_width[PAMR_MAX_DIM];   /* ghost-widths to add for splitting grids */
   int periodic[PAMR_MAX_DIM];      /* whether given boundary is periodic or not */
   int min_width[PAMR_MAX_DIM];     /* the minimum allowed shape of a grid */
   int gdm;                         /* grid distribution method */
   int interp_buffer;               /* during interpolation, an optional buffer zones to pad interior sources */

   /* MG options */
   
   int MG_min_cwidth[PAMR_MAX_DIM]; /* the smallest allowed size of a coarsest grid  */

   /* iterator state, and stack for nested looping */

   grid *iter_g;
   int iter_all;
   int iter_L;
   grid *iter_g_stack[ITER_G_STACKSIZE];
   int iter_all_stack[ITER_G_STACKSIZE];
   int iter_L_stack[ITER_G_STACKSIZE];
   int iter_g_stack_top;

   /* excision functions */

   int excision_on;                 /* if true, we excise */
   void (*app_fill_ex_mask)(real *mask, int dim, int *shape, real *bbox, real excised); 
                                    /* user function to call to fill excision mask */
   real excised;                    /* value denoting an excised grid point */
   int mg_mask_gfn;                 /* mask variable grid function numbers */
   int amr_mask_gfn;             

} context;

extern context *curr_context;
extern context *contexts[PAMR_MAX_CONTEXTS];

extern int PAMR_parallel_io; /* whether to do disk i/o in parallel or not ... currently only works with cp stuff */

extern int PAMR_stats_on; /* collect a variety of stats ... only the two below so far */
extern int PAMR_comm_r_secs,PAMR_comm_r_microsecs; /* time spent and memory used communicating via Isend and Irecieve in transfer.c */
extern int PAMR_comm_r_GB,PAMR_comm_r_B; 
extern int PAMR_comm_s_secs,PAMR_comm_s_microsecs; /* time spent and memory used communicating via Isend and Irecieve in transfer.c */
extern int PAMR_comm_s_GB,PAMR_comm_s_B;
extern int PAMR_num_s,PAMR_num_r;

/* see .c file for function comments */

int is_in_amrh(int gfn);
int is_in_mg_and_amrh(int gfn, int *num_tl);
var *new_var(char *name);
var *find_var(char *name);
void recalc_lev_info(void);
void free_context_mem();
int save_sgh_bboxes(int min_lev, int max_lev, int num, int *lev, real *bbox);
cgh *build_sgh(int min_lev, int max_lev, int num, int *lev, real *bbox, real t);
cgh *sgh_cull_overlap(cgh *sgh);
void free_cgh(cgh *c, int min_lev, int max_lev);
cgh *compose_cgh(cgh *sgh, int alloc_mg);
void init_new_cgh(cgh *new_cgh, cgh *curr_cgh);
int incorporate_new_cgh(cgh *new_cgh, cgh *curr_cgh);
void insert_grid(level *lev, grid *g);
int amr_to_mgh(level *amr_lev, level **mg_lev, int tl);
int fill_mgh_level(level *lev, int num, real *bbox, real *island_no);
int build_mgh(int min_lev, int max_lev, int tl);
void destroy_mgh(void);
void calc_island_no(real *bboxes, int num, real *island_no, real *dx, int cf);
int sizeof_data(grid *g);
int merge_bboxes(real *bboxes, int num, real min_eff);
void init_ex_mask(int Lmin, int Lmax, cgh *gh, int gfn);
void set_wrap(level *lv);
int abuts(grid *g, int d, int pm);
int matched_abuts(grid *gl, int d, int pm, grid *gm);
grid *make_virtual_copy(grid *g);
void alloc_virtual_grids(level *l);
void free_virtual_grids(level *l);

/* alternative memory allocation option not yet supported */
#define pmalloc(x) malloc(x)
#define pfree(x) free(x)

#endif /* _GH_H */
