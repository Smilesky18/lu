# include <stdio.h>
# include <stdlib.h>
# include "../Include/lu.h"
# include <sys/time.h>
# include <float.h>
# include <immintrin.h>
//# include <avx2intrin.h>
//# include <avxintrin.h>
# define MICRO_IN_SEC 1000000.00

/* Time Stamp */
/*double microtime()
{
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv,&tz);

        return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}*/
/* a ?= b */
/*bool equal( double a, double b )
{
  if ( Abs(a-b) < 0.001)
  {
    return true;
  }
  else
  {
    return false;
  }
}*/
/* return |x| */
/*double Abs(double x)
{
  return x < 0 ? -x : x;
}*/

  typedef __attribute__((aligned(64))) union
  {
    __m256d vec;
    double ptr_vec[4];
  }v2df_t;

  typedef union
  {
    __m128i vec;
    int ptr_vec[4];
  }v2if_t;
double* lu_gp_sparse_avx2(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U)
{
  int sum_l = 0, sum_u = 0, row, sum_for = 0;
  double *L, *U, *xx;
  int *row_index, row_column;
  double U_diag;
  double start, end;
  int j, k, current_column;
  int i;
  L = ( double * )_mm_malloc(sizeof(double) * nzl, 64);
  U = ( double * )_mm_malloc(sizeof(double) * nzu, 64);
  xx = ( double *)_mm_malloc(sizeof(double) * n, 64 );
  
  /* Array xx initialization*/
  for ( i = 0; i < n; i++ )
  {
    xx[i] = 0;
    L[xa_L[i]] = 1.0;
  }
  

  v2df_t v_l, v_row, v_mul, v_sub, v;
  v2if_t vi;
  int column_divid, m;
  double temp;

  /* column-oriented G/P algorithm without partial pivoting */
  for ( k = 0; k < n; k++ )
  {
    current_column = perm_c[k];
    
    /* xx[] = A[:,current_column] */
    for ( j = xa[current_column]; j < xa[current_column+1]; j++ )
    {
      xx[perm_r[asub[j]]] = a[j];
    }
    row_column = xa_U[k+1] - xa_U[k] - 1;
    row_index = (int *)malloc( sizeof(int) * row_column);
    for ( j = 0; j < row_column; j++ )
    {
        row_index[j] = asub_U[j+xa_U[k]];
    }
    
    /* L[:,0~k]*xx = A[:,current_column], solve for xx*/
    for ( j = 0; j < row_column; j++ )
    {
        row = row_index[j];
        column_divid = ( xa_L[row+1] - xa_L[row] - 1 ) / 4;
/*        v_row.ptr_vec[0] = xx[row];
        v_row.ptr_vec[1] = xx[row];
        v_row.ptr_vec[2] = xx[row];
        v_row.ptr_vec[3] = xx[row];*/
        v_row.vec = _mm256_set1_pd(xx[row]);
        for ( i = xa_L[row]+1; i < xa_L[row]+1+column_divid*4; i+=4 )
        {
/*            v_l.ptr_vec[0] = L[i];
            v_l.ptr_vec[1] = L[i+1];
            v_l.ptr_vec[2] = L[i+2];
            v_l.ptr_vec[3] = L[i+3]; */
            v_l.vec = _mm256_loadu_pd(&L[i]);
            
     //       v_l.vec = _mm256_set1_pd(1);
     //       v.vec = _mm256_set1_pd(i);
/*            vi.ptr_vec[0] = asub_L[i];
            vi.ptr_vec[1] = asub_L[i+1];
            vi.ptr_vec[2] = asub_L[i+2];
            vi.ptr_vec[3] = asub_L[i+3]; */
            vi.vec = _mm_loadu_si128((__m128i const *)&asub_L[i]);

/*            v.ptr_vec[0] = xx[vi.ptr_vec[0]];
            v.ptr_vec[1] = xx[vi.ptr_vec[1]];
            v.ptr_vec[2] = xx[vi.ptr_vec[2]];
            v.ptr_vec[3] = xx[vi.ptr_vec[3]]; */
            v.vec = _mm256_i32gather_pd(&xx[0], vi.vec, 8);

/*            v_mul.ptr_vec[0] = v_l.ptr_vec[0] * v_row.ptr_vec[0];
            v_mul.ptr_vec[1] = v_l.ptr_vec[1] * v_row.ptr_vec[1];
            v_mul.ptr_vec[2] = v_l.ptr_vec[2] * v_row.ptr_vec[2];
            v_mul.ptr_vec[3] = v_l.ptr_vec[3] * v_row.ptr_vec[3]; 
  //          v_mul.vec = _mm256_mul_pd(v_l.vec, v_row.vec);

            v_sub.ptr_vec[0] = v.ptr_vec[0] - v_mul.ptr_vec[0];
            v_sub.ptr_vec[1] = v.ptr_vec[1] - v_mul.ptr_vec[1];
            v_sub.ptr_vec[2] = v.ptr_vec[2] - v_mul.ptr_vec[2];
            v_sub.ptr_vec[3] = v.ptr_vec[3] - v_mul.ptr_vec[3]; */
    //        v_sub.vec = _mm256_sub_pd(v.vec, v_mul.vec); 
            v_sub.vec = _mm256_fnmadd_pd(v_l.vec, v_row.vec, v.vec);

         //   v_sub.vec = v.vec;
         //   _mm256_store_pd(&xx[asub_L[i]], v_sub.vec);
            xx[asub_L[i]] = v_sub.ptr_vec[0];
            xx[asub_L[i+1]] = v_sub.ptr_vec[1];
            xx[asub_L[i+2]] = v_sub.ptr_vec[2];
            xx[asub_L[i+3]] = v_sub.ptr_vec[3]; 

        }
        for ( m = xa_L[row]+1+column_divid*4; m < xa_L[row+1]; m++ )
        {
            xx[asub_L[m]] -=  xx[row]*L[m];
        }
    }

    /* solve for U[:,k]*/
    for ( i = xa_U[k]; i < xa_U[k+1]; i++ )
    {
      U[i] = xx[asub_U[i]];
      xx[asub_U[i]] = 0;
    }

    /* solve for L[:,k] */
    U_diag = U[i-1];
    for ( i = xa_L[k]+1; i < xa_L[k+1]; i++ )
    {
      L[i] = xx[asub_L[i]] / U_diag;
      xx[asub_L[i]] = 0;
    } 
  }
   /* solve for Ly = b and Ux = y */
   double *y, *x;
   y = ( double *)malloc( sizeof( double ) * n );
   x = ( double *)malloc( sizeof( double ) * n );


  for ( i = 0; i < n; i++ )
  {
    y[i] = 1.0;
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = xa_L[i]+1; j < xa_L[i+1]; j++ )
    {
        y[asub_L[j]] -= y[i] * L[j];
    }
  }

  //x[n-1] = y[n-1];
  for ( i = 0; i < n; i++ )
  {
    x[i] = y[i];
  }
  x[n-1] = y[n-1]/U[nzu-1];
  for ( i = n-1; i > 0; i-- )
  {
    for ( j = xa_U[i]; j < xa_U[i+1]-1; j++ )
    {
        x[asub_U[j]] -= x[i] *U[j];
    }
    x[i-1] = x[i-1]/U[xa_U[i]-1];
  } 
  
  return x;
}

