#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Include/lu.h"
#include "../Include/nicslu.h"
#include "../Include/nicslu_util.h"
#include <sys/time.h>
#include <stdbool.h>
#define MICRO_IN_SEC 1000000.00

int belong( int column, int *xa_belong, int *asub_belong, int *length_belong)
{
	int i, temp = length_belong[column];

	for ( i = xa_belong[column]; i < xa_belong[column+1] - 1; i++ )
	{
		if ( length_belong[ asub_belong[i] ] > temp ) temp = length_belong[ asub_belong[i] ];
	}
	return temp+1;
}

bool equal( double a, double b )
{
	  if ( Abs(a-b) < 0.000001 )
	  {
		  return true;
	  }
	  else
	  {
		  return false;
	  }
}
double Abs(double x)
{
	  return x < 0 ? -x : x;
} 

void BubbleSort( int *arr, int size )
{
    int i, j, tmp;
    for ( i = 0; i < size - 1; i++ )
    {
        for ( j = 0; j < size - i - 1; j++ )
        {
            if ( arr[j] > arr[j+1] )
            {
                tmp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = tmp;
            }
        }
    }
}
/* Time Stamp */
double microtime()
{
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv,&tz);

        return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}
int main(int argc[], char *argv[])
{
	int ret;
	uint__t n, nnz, i;
	real__t *ax;
	uint__t *ai, *ap;
	SNicsLU *nicslu;
	real__t *x, *b, err;
	ax = NULL;
	ai = NULL;
	ap = NULL;
	x = NULL;
	b = NULL;

	nicslu = (SNicsLU *)malloc(sizeof(SNicsLU));
	NicsLU_Initialize(nicslu);

	ret = NicsLU_ReadTripletColumnToSparse(argv[1], &n, &nnz, &ax, &ai, &ap);
	if (ret != NICS_OK) goto EXIT;

	x = (real__t *)malloc(sizeof(real__t)*(n+n));
	b = x + n;
	for (i=0; i<n+n; ++i) x[i] = 1.;

	NicsLU_CreateMatrix(nicslu, n, nnz, ax, ai, ap);
	nicslu->cfgf[0] = 1.;

//	for ( i = 0; i <= n; i++ ) printf("%d\n", ap[i]);

	printf("****************%s***************\n", argv[1]);
	NicsLU_Analyze(nicslu);
//	printf("analysis time: %.8g\n", nicslu->stat[0]);

	NicsLU_Factorize(nicslu);
//	printf("factorization time: %.8g\n", nicslu->stat[1]);

	int *rowperm = nicslu->row_perm;
	int *colperm = (int *)malloc(sizeof(int) * n);
	for ( i = 0; i < n; i++ ) colperm[i] = nicslu->pivot_inv[nicslu->col_perm_inv[i]];
	int *asub_L, *xa_L, *asub_U, *xa_U;
	void *lu_array = nicslu->lu_array;
	int *ulen = nicslu->ulen, *llen = nicslu->llen;
	size_t *up = nicslu->up;
	int u_row_length, l_row_length, *u_row_index, *l_row_index, lnz = n, unz = n;
	int k;
    xa_L = (int *)malloc(sizeof(int) * (n+1));	
	xa_U = (int *)malloc(sizeof(int) * (n+1));
	memset(xa_L, 0, sizeof(int) * (n+1));
	memset(xa_U, 0, sizeof(int) * (n+1));
	for ( i = 1; i <= n; i++ )
	{
		xa_L[i] = xa_L[i-1] + ulen[i-1] + 1;
		xa_U[i] = xa_U[i-1] + llen[i-1] + 1;

/*		lnz += ulen[i];
		unz += llen[i];
		u_row_length = ulen[i];
		l_row_length = llen[i];
		l_row_index = (int *)(((char *)lu_array) + up[i] + u_row_length*12);
		u_row_index = (int *)(((char *)lu_array) + up[i]);

//		for ( k = 0; k < u_row_length; k++ ) printf("%d %d\n", u_row_index[k], i);
		for ( k = 0; k < l_row_length; k++ ) printf("%d %d\n", l_row_index[k], i); */
	}
	asub_L = (int *)malloc(sizeof(int) * xa_L[n]);
	asub_U = (int *)malloc(sizeof(int) * xa_U[n]);
	for ( i = 0; i < n; i++ )
	{
		asub_L[ xa_L[i] ] = i;
		u_row_length = ulen[i];
		l_row_length = llen[i];
		u_row_index = (int *)(((char *)lu_array) + up[i]);
		l_row_index = (int *)(((char *)lu_array) + up[i] + u_row_length*12);
		for ( k = 0; k < u_row_length; k++ )
		{
			asub_L[k + xa_L[i] + 1] = u_row_index[k]; 
		}
		for ( k = 0; k < l_row_length; k++ )
		{
			asub_U[k + xa_U[i]] = l_row_index[k]; 
		}
		asub_U[xa_U[i+1] - 1] = i;
	}

	int size, size_u;
	for ( i = 0; i < n; i++ )
    {
        size = xa_L[i+1] - xa_L[i]; 
        BubbleSort( &asub_L[xa_L[i]], size );
		
		size_u = xa_U[i+1] - xa_U[i]; 
        BubbleSort( &asub_U[xa_U[i]], size_u );
    }

	NicsLU_ReFactorize(nicslu, ax);
//	printf("re-factorization time: %.8g\n", nicslu->stat[2]);

	NicsLU_Solve(nicslu, x);
//	printf("substitution time: %.8g\n", nicslu->stat[3]);

	printf("Time of nicslu is: %.8g\n", nicslu->stat[2]+nicslu->stat[3]);
	printf("NNZ(L+U-I): %ld\n", nicslu->lu_nnz);

	/*parallel level*/
	int *length = ( int *)malloc( sizeof(int) * n);
	memset(length, -1, sizeof(int) * n);
	for ( i = 0; i < n; i++ )
	{
		if ( llen[i] == 0 ) length[i] = 0;
		else
		{
			length[i] = belong(i, xa_U, asub_U, length);
		}
	}

	int max_level = 0;
	for ( i = 0; i < n; i++ )
	{
		if ( length[i] > max_level ) max_level = length[i];
	}

	printf("level = %d\n", max_level);
	int *level = ( int *)malloc( sizeof(int) *( max_level+1 ) );
	memset(level, 0, sizeof(int) *(  max_level+1 ));
	for ( i = 0; i < n; i++ )
	{
		level[length[i]]++;
	}
	
	int *xa = ( int *)malloc( sizeof(int) * (max_level+2));
	memset(xa, 0, sizeof(int) * (max_level+2));
	
	for ( i = 1; i <= max_level+1; i++ )
	{
		xa[i] = xa[i-1] + level[i-1];
	}

	int *xa_trans = (int *)malloc( sizeof(int) * (max_level+2));
	for ( i = 0; i < max_level+2; i++ ) xa_trans[i] = xa[i];
	int *asub_U_level = (int *)malloc(sizeof(int) * n);
	for ( i = 0; i < n; i++ )
	{
		asub_U_level[ xa[length[i]]++ ] = i;
	}
	
     double *l_gp, *l_gp_blas, *l_gp_sn_parallel, *l_gp_sn_v2_parallel;
    int blas_error = 0, sn_error = 0, sn_error_v2 = 0;
    double time_blas_start, time_blas_end, time_sn_parallel_start, time_sn_parallel_end, t1, t2;
	double start_GP, end_GP;
  
	start_GP = microtime();	
    l_gp = lu_gp_sparse(ax, ai, ap, n, xa_L[n], xa_U[n], rowperm, colperm, asub_L, xa_L, asub_U, xa_U);
    end_GP = microtime() - start_GP;
    printf("Time of LU_GP: %lf\n", end_GP);

    time_blas_start = microtime();
    l_gp_blas = lu_gp_sparse_parallel(ax, ai, ap, n, xa_L[n], xa_U[n], rowperm, colperm, asub_L, xa_L, asub_U, xa_U, max_level,
	xa_trans, asub_U_level);
    time_blas_end = microtime() - time_blas_start;
    printf("Time of LU_GP_parallel: %lf\n", time_blas_end); 
	
/*	time_sn_parallel_start = microtime();
    l_gp_sn_parallel = lu_gp_sparse_supernode_computing_parallel(Ax, Ai, Ap, n, nzl, nzu, q, pinv, Li, Lp, Ui, Up,sn_record, max_level, xa_trans, asub_U );
    time_sn_parallel_end = microtime() - time_sn_parallel_start;
    printf("Time of LU_GP__SN_parallel: %lf\n", time_sn_parallel_end); */
   
/*	t1 = microtime();
    l_gp_sn_v2_parallel = lu_gp_sparse_supernode_computing_v2_parallel(Ax, Ai, Ap, n, nzl, nzu, q, pinv, Li, Lp, Ui, Up,sn_record, max_level, xa_trans, asub_U, sn_record_end);
    t2 = microtime() - t1;
    printf("Time of LU_GP_SN_v2_parallel: %lf\n", t2); */

	double *x_lu = (double *)malloc(sizeof(double) * n);
	for ( i = 0; i < n; i++ )
	{
		x_lu[rowperm[i]] = l_gp[i];
	}
	double *x_lu_parallel = (double *)malloc(sizeof(double) * n);
	for ( i = 0; i < n; i++ )
	{
		x_lu_parallel[rowperm[i]] = l_gp_blas[i];
	}
    for ( i = 0; i < n; i++ )
    {
         if ( !equal(x_lu[i], x[i]) ) 
         {
             //printf("nicslu[%d] = %.16f me[%d] = %.16f\n", i, x[i], i, x_lu[i]);
             blas_error++;
         }
		 if ( !equal(x_lu_parallel[i], x[i]) ) 
         {
             //printf("nicslu[%d] = %.16f me_parallel[%d] = %.16f\n", i, x[i], i, x_lu_parallel[i]);
             sn_error++;
         }
    }
	
    printf("sequential error results are: %d \n", sn_error); 
    printf("parallel error results are: %d \n", sn_error); 
	
EXIT:
	NicsLU_Destroy(nicslu);
	free(ax);
	free(ai);
	free(ap);
	free(nicslu);
	free(x);
	return 0;
}
