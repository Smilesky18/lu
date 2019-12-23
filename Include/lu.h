#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
// #include <slu_ddefs.h>
#include <cs.h>



/*! \brief Driver routines */
extern double* lu_gp_sparse(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_sparse_avx2(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *);
extern void* lu_gp_sparse_row_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int * );
extern double* lu_gp_sparse_supernode_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_sparse_supernode_computing_blas(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_sparse_sn_u(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int *);
extern double Abs(double );
extern double microtime();


