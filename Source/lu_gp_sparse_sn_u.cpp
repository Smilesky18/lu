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
// double Abs(double x)
// {
//   return x < 0 ? -x : x;
// }

double* lu_gp_sparse_sn_u(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U, int *sn_record, int *index)
{
  int sum_l = 0, sum_u = 0, row, row_else, sum_for = 0, ii, pack;
  double *L, *U, *xx, *xx_else;
  int element_num;
  int *row_index, row_column, *row_index_else, row_column_else;
  double U_diag;
  double start, end;
  int j, k, current_column, current_column_else;
  int i;
  FILE *fp_L_Error = fopen("Error_Result/L_Error.txt", "w");
  FILE *fp_U_Error = fopen("Error_Result/U_Error.txt", "w");
  L = ( double * )malloc(sizeof(double *) * nzl);
  U = ( double * )malloc(sizeof(double *) * nzu);
//   xx = ( double *)malloc(sizeof(double) * n );
  
  /* Array xx initialization*/
  for ( i = 0; i < n; i++ )
  {
//      xx[i] = 0;
     L[xa_L[i]] = 1.0;
  }

  //double sum_temp = 0;
  /* column-oriented G/P algorithm without partial pivoting */
//   start = microtime();
  for ( k = 0; k < n; k+=pack )
  {
     if ( sn_record[k] == -1 )
     {
        current_column = perm_c[k];
        xx = ( double *)malloc(sizeof(double) * n );
//         memset(xx, 0, sizeof(int) * n );
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
            for ( i = xa_L[row]+1; i < xa_L[row+1]; i++ )
            {
                xx[asub_L[i]] -=  xx[row]*L[i];
            }
        }

        /* solve for U[:,k]*/
        for ( i = xa_U[k]; i < xa_U[k+1]; i++ )
        {
            U[i] = xx[asub_U[i]];
            // if ( xx[asub_U[i]] == 0 ) printf("yes\n");
            xx[asub_U[i]] = 0;
        }

        /* solve for L[:,k] */
        U_diag = U[i-1];
        for ( i = xa_L[k]+1; i < xa_L[k+1]; i++ )
        {
            L[i] = xx[asub_L[i]] / U_diag;
            // if ( xx[asub_L[i]] == 0 ) printf("yes yes\n");
            xx[asub_L[i]] = 0;
        } 
        pack =1;
        free(xx);
     }
     else
     {
         element_num = index[sn_record[k]] * n;
//          printf("index[sn_record[%d] = %d\n", k, index[sn_record[k]]);
//          printf("k = %d element_num = %d\n", k, element_num);
         xx_else = ( double *)malloc(sizeof(double) * element_num );
         for ( i = 0; i < index[sn_record[k]]; i++ )
         {
             current_column_else = perm_c[k+i]; 
             for ( j = xa[current_column_else]; j < xa[current_column_else+1]; j++ )
             {
                 xx_else[perm_r[asub[j]]+i*n] = a[j];
             }
         }
        row_column_else = xa_U[k+1] - xa_U[k] - 1;
        row_index_else = (int *)malloc( sizeof(int) * row_column_else);
        for ( j = 0; j < row_column_else; j++ )
        {
            row_index_else[j] = asub_U[j+xa_U[k]];
        }
        for ( ii = 0; ii < index[sn_record[k]]; ii++ )
        {
            for ( j = 0; j < row_column_else; j++ )
            {
                row_else = row_index_else[j];
                for ( i = xa_L[row_else]+1; i < xa_L[row_else+1]; i++ )
                {
                    xx_else[asub_L[i]+ii*n] -=  xx_else[row_else+ii*n]*L[i];
                }
            }
        }
        for ( ii = 0; ii < index[sn_record[k]]; ii++ )
        {
            for ( i = xa_U[k+ii]; i < xa_U[k+ii+1]; i++ )
            {
                U[i] = xx_else[asub_U[i]+ii*n];
                // if ( xx[asub_U[i]] == 0 ) printf("yes\n");
                xx_else[asub_U[i]+ii*n] = 0;
            }

            /* solve for L[:,k] */
            U_diag = U[i-1];
//             printf("U_diag = %lf\n", U_diag);
            for ( i = xa_L[k+ii]+1; i < xa_L[k+ii+1]; i++ )
            {
                L[i] = xx_else[asub_L[i]+ii*n] / U_diag;
                // if ( xx[asub_L[i]] == 0 ) printf("yes yes\n");
                xx_else[asub_L[i]+ii*n] = 0;
            } 
        }
        pack = index[sn_record[k]];
        free(xx_else);
    }
  }

//    for ( i = 0; i < nzu; i++ )
//    {
//        printf("U[%d] = %lf\n", i, U[i]);
// }
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

