# include <stdio.h>
# include <stdlib.h>
# include "../Include/lu.h"
# include <float.h>
# include <omp.h>

double* lu_gp_sparse_parallel(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U, int level, int *order_k_xa, int *order_k_asub)
{
  int row;
  double *L, *U, *xx;
  int *row_index, row_column;
  double U_diag;
  int j, k, current_column;
  int i;
  double temp;
  int n_xx = 32 * n;
  L = ( double * )malloc(sizeof(double) * nzl);
  U = ( double * )malloc(sizeof(double) * nzu);
  xx = ( double *)malloc(sizeof(double) * n_xx);
  memset( xx, 0, sizeof(double)*n_xx );
  
  /* Array xx initialization*/
  for ( i = 0; i < n; i++ )
  {
    L[xa_L[i]] = 1.0;
  }

  int qq, rr, ss, tt, zz;
  int jj, jjj, jjjj, ii, iii, iiii;

  for ( qq = 0; qq <= level; qq++ )
  {
	  rr = order_k_xa[qq];
	  ss = order_k_xa[qq+1];
	  tt = ss - rr;
	  if ( tt > 30 )
	   {
		  #pragma omp parallel for private(k, current_column, jj, jjj, ii, iii, row_column, row, temp, i, U_diag) num_threads(32) schedule(static)
		  for ( zz = 0; zz < tt; zz++ )
			{
				int thread_num = omp_get_thread_num();
				k = order_k_asub[zz+rr];
				current_column = perm_c[k];
				for ( jj = xa[current_column]; jj < xa[current_column+1]; jj++ )
				{	
					xx[ perm_r[asub[jj]] + thread_num*n ] = a[jj];
				}
				row_column = xa_U[k+1] - xa_U[k] - 1;	
				for ( jjj = 0; jjj < row_column; jjj++ )
				{
					row = asub_U[jjj+xa_U[k]];
					temp = xx[row + thread_num*n];
					for ( i = xa_L[row]+1; i < xa_L[row+1]; i++ )
					{
						xx[asub_L[i] + thread_num*n] -=  temp*L[i];
					}
				}	
				/* solve for U[:,k]*/
				for ( ii = xa_U[k]; ii < xa_U[k+1]; ii++ )
				{
					U[ii] = xx[asub_U[ii] + thread_num*n];
					xx[asub_U[ii] + thread_num*n] = 0;
				}
				/* solve for L[:,k] */
				U_diag = U[ii-1];
				for ( iii = xa_L[k]+1; iii < xa_L[k+1]; iii++ )
				{
					L[iii] = xx[asub_L[iii] + thread_num*n] / U_diag;
					xx[asub_L[iii] + thread_num*n] = 0;
				}
			}
	    }
		else
		{
			for ( zz = 0; zz < tt; zz++ )
			{
				k = order_k_asub[zz+rr];
				current_column = perm_c[k];
				for ( jj = xa[current_column]; jj < xa[current_column+1]; jj++ )
				{	
					xx[ perm_r[asub[jj]] ] = a[jj];
				}
				row_column = xa_U[k+1] - xa_U[k] - 1;	
				for ( jjj = 0; jjj < row_column; jjj++ )
				{
					row = asub_U[jjj+xa_U[k]];
					temp = xx[row];
					for ( i = xa_L[row]+1; i < xa_L[row+1]; i++ )
					{
						xx[asub_L[i]] -=  temp*L[i];
					}
				}	
				/* solve for U[:,k]*/
				for ( ii = xa_U[k]; ii < xa_U[k+1]; ii++ )
				{
					U[ii] = xx[asub_U[ii]];
					xx[asub_U[ii]] = 0;
				}
				/* solve for L[:,k] */
				U_diag = U[ii-1];
				for ( iii = xa_L[k]+1; iii < xa_L[k+1]; iii++ )
				{
					L[iii] = xx[asub_L[iii]] / U_diag;
					xx[asub_L[iii]] = 0;
				}
			}
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

