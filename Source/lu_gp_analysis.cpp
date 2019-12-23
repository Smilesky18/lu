# include <stdio.h>
# include <stdlib.h>
# include "../Include/lu.h"
# include <sys/time.h>
# include <float.h>
# define MICRO_IN_SEC 1000000.00

/* Time Stamp */
double microtime()
{
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv,&tz);

        return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}
/* a ?= b */
bool equal( double a, double b )
{
  if ( Abs(a-b) < 0.001)
  {
    return true;
  }
  else
  {
    return false;
  }
}
/* return |x| */
double Abs(double x)
{
  return x < 0 ? -x : x;
}

void* lu_gp_analysis(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U, double *l_data, double *u_data)
{
    
}

