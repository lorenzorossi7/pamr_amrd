#ifndef _MISC_H
#define _MISC_H 1
/*======================================================================================== */
/* misc.h                                                                                  */
/*                                                                                         */
/* miscellaneous stuff                                                                     */
/*======================================================================================== */

#define IFL if (ltrace)
#define IFL0 if (ltrace && curr_context->rank==0)

extern int gtrace;

#define IFG(x) if (gtrace>=(x))
#define IFG0(x) if (gtrace>=(x) && curr_context->rank==0)

#define IFGL(x) if (gtrace>=(x) || ltrace)
#define IFGL0(x) if ((gtrace>=(x) || ltrace) && curr_context->rank==0)

int fuzz_eq(double x1, double x2, double fuzz);
int fuzz_lt(double x1, double x2, double fuzz);
int fuzz_lte(double x1, double x2, double fuzz);
int fuzz_gt(double x1, double x2, double fuzz);
int fuzz_gte(double x1, double x2, double fuzz);

#define max(a, b)  (((a) > (b)) ? (a) : (b)) 
#define min(a, b)  (((a) < (b)) ? (a) : (b)) 

#endif /* _MISC_H */
