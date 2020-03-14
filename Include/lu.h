#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
// #include <slu_ddefs.h>
//#include <cs.h>



/*! \brief Driver routines */
//extern int min( int , int )
extern double* lu_gp_sparse(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_sparse_parallel(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int, int *, int *);
extern double* lu_gp_sparse_v2(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_sparse_blocking(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *);
//extern double* lu_gp_sparse_row_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_sparse_avx2(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *);
extern void* lu_gp_sparse_row_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int * );
extern double* lu_gp_sparse_supernode_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_sparse_supernode_computing_parallel(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int, int *, int *);
extern double* lu_gp_sparse_supernode_computing_v2_parallel(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int, int *, int *, int *);
extern double* lu_gp_sparse_supernode_computing_blas(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int );
extern double* lu_gp_sparse_supernode_computing_nontri_blas(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int, int );
extern double* lu_gp_sparse_sn_u(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int *);
extern double Abs(double );
extern double microtime();
extern double* lu_gp_sparse_supernode_sparse_row_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int );
extern double* lu_gp_sparse_supernode_sparse_column_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int );
extern double* lu_gp_sparse_supernode_dense_row_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int );
extern double* lu_gp_sparse_supernode_dense_column_computing(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int );
extern double* lu_gp_u_sn(double *, int *, int *, int , int , int , int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_u_sn_all_dense_vec(double *, int *, int *, int , int , int , int *, int *, int *, int *, int *, int *, int *, int *, int *);
extern double* lu_gp_u_sn_new(double *, int *, int *, int, int, int, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
