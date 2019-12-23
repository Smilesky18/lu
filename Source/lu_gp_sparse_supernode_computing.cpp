# include <stdio.h>
# include <stdlib.h>
# include "../Include/lu.h"
# include <sys/time.h>
# include <float.h>
# define MICRO_IN_SEC 1000000.00

/* Time Stamp */
// double microtime()
// {
//         struct timeval tv;
//         struct timezone tz;
//         gettimeofday(&tv,&tz);
// 
//         return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
// }
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
int belong ( int x, int *array, int length )
{
    int i;
    for ( i = 0; i < length; i++ )
    {
        if ( x == array[i] ) return 1;
    }
    return 0;
}

double* lu_gp_sparse_supernode_computing(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U, int *sn_record)
{
  int row;
  double *L, *U, *xx;
  int *row_index;
  double row_column;
  double U_diag;
  int j, k, current_column;
  int i;
  L = ( double * )malloc(sizeof(double *) * nzl);
  U = ( double * )malloc(sizeof(double *) * nzu);
  xx = ( double *)malloc(sizeof(double) * n );
  
  /* Array xx initialization*/
  for ( i = 0; i < n; i++ )
  {
    xx[i] = 0;
    L[xa_L [i]] = 1.0;
  }

  /* column-oriented G/P algorithm without partial pivoting */
  int *index;
  int row_number_sn, row_num;
  double row_count;
  int pack, r, t, p, w, ww, qq, ss, qqq, sss;
  int row_else1, row_else2, row_if;
  double *dense_vec, sum = 0;
  int row_sn_start, row_sn_end;
  index = ( int * )malloc((sizeof(int) * n));
  memset(index, 0, sizeof(int) * n);
  double temp;
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
    /*for ( j = 0; j < row_column; j++ )
    {
        row_index[j] = asub_U[j+xa_U[k]];
    }*/
    for ( j = 0; j < row_column; j++ )
    {
        row_index[j] = asub_U[j+xa_U[k]];
        row = asub_U[j+xa_U[k]];
        if ( sn_record[row] != -1 ) 
        {
            index[sn_record[row]]++;
        }
    }
    /* L[:,0~k]*xx = A[:,current_column], solve for xx*/
    for ( j = 0; j < row_column; j+=pack )
    {
        row_if = row_index[j];
        temp = xx[row_if];
        /*if ( index[sn_record[row_if]] == 1 || sn_record[row_if] == 0 )
        {
        //    sum2++;
            for ( i = xa_L[row_if]+1; i < xa_L[row_if+1]; i++ )
            {
                xx[asub_L[i]] -=  xx[row_if]*L[i];
            }
            index[sn_record[row_if]] = 0;
            pack = 1;
        }*/
        if ( sn_record[row_if] == -1 )
        {
            for ( i = xa_L[row_if]+1; i < xa_L[row_if+1]; i++ )
            {
                xx[asub_L[i]] -=  temp*L[i];
            }
            pack = 1;
        }
        else if ( index[sn_record[row_if]] == 1 )
        {
            for ( i = xa_L[row_if]+1; i < xa_L[row_if+1]; i++ )
            {
                xx[asub_L[i]] -=  temp*L[i];
            }
            index[sn_record[row_if]] = 0;
            pack = 1;
        }
        else
        {   
            row_number_sn = index[sn_record[row_if]];
            row_sn_start = row_if;
            row_sn_end = row_index[j+row_number_sn-1];
            index[sn_record[row_if]] = 0;
            row_num = xa_L[row_sn_end+1] - xa_L[row_sn_end] - 1;
            row_count = row_num * row_number_sn;
            dense_vec = ( double *) malloc(sizeof(double) * row_count);
            ww = xa_L[row_sn_start+1] - row_num;
            for ( qqq = 0; qqq < row_number_sn; qqq++ )
            {
                row_else1 = row_index[j+qqq];
                for ( sss = xa_L[row_else1]+1; sss < xa_L[row_else1+1]-row_num; sss++ )
                {
                    xx[asub_L[sss]] -=  xx[row_else1]*L[sss];
                }
            }
            for ( r = 0; r < row_number_sn; r++ )
            {
                row_else2 = row_index[j+r];
                w = xa_L[row_else2+1] - row_num;
                for ( t = 0; t < row_num; t++ )
                {
                    dense_vec[r+row_number_sn*t] = L[w+t];
                }
            }

            for ( qq = 0; qq < row_num; qq++ )
            {
                p = qq * row_number_sn;
                for ( ss = 0; ss < row_number_sn; ss++ )
                {
                    sum += dense_vec[p+ss]*xx[row_index[j+ss]];
                }
                xx[asub_L[ww+qq]] -= sum;
                sum = 0;
            }
            pack = row_number_sn;
            free(dense_vec);
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

