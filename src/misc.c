//=============================================================================
// misc.c --- miscellaneous routines
//=============================================================================

#include "misc.h"
#include <math.h>

// 'global' trace variable
int gtrace=0;

// fuzzy logic
// 'fuzz' is assumed positive

int fuzz_eq(double x1, double x2, double fuzz)
{
   if (fabs(x2-x1)<fuzz) return 1; else return 0;
}
 
int fuzz_lt(double x1, double x2, double fuzz)
{
   if (x1<-fuzz+x2) return 1; else return 0;
}
 
int fuzz_lte(double x1, double x2, double fuzz)
{
   if (x1<=fuzz+x2) return 1; else return 0;
}

int fuzz_gt(double x1, double x2, double fuzz)
{
   if (x1>fuzz+x2) return 1; else return 0;
} 

int fuzz_gte(double x1, double x2, double fuzz)
{
   if (x1>=-fuzz+x2) return 1; else return 0;
} 
