# include <stdio.h>
# include <stdlib.h>
# include "../Include/lu.h"
# include <sys/time.h>
# include <float.h>
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
}
/* return |x| */
/*double Abs(double x)
{
  return x < 0 ? -x : x;
}*/

void* lu_gp_sparse_row_computing(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U)
{
  int sum_l = 0, sum_u = 0, row, sum_for = 0;
  int *Li, *Lp, *ptr_l_row, *row_count_l, row_num_l;
  double *L, *U, *xx;
  int *row_index, row_column;
  double U_diag;
  double *y, *x;
  int j, k, current_column;
  int i;
  FILE *fp_L_Error = fopen("Error_Result/L_Error.txt", "w");
  FILE *fp_U_Error = fopen("Error_Result/U_Error.txt", "w");
  L = ( double * )malloc(sizeof(double *) * nzl);
  U = ( double * )malloc(sizeof(double *) * nzu);
  Li = ( int * )malloc(sizeof(int *) * nzl);
  Lp = ( int * )malloc(sizeof(int *) * n+1);
  y = ( double * )malloc(sizeof(double *) * n);
  x = ( double * )malloc(sizeof(double *) * n);
  row_count_l = ( int * )malloc(sizeof(int *) * n);
  xx = ( double *)malloc(sizeof(double) * n );
  ptr_l_row = ( int * )malloc(sizeof(int *) * n);
 
  memset(L, 0, sizeof(double)*nzl);
  memset(U, 0, sizeof(double)*nzu);
 
  /*for ( i = 0; i < 23; i++ )
  {
    printf("a[%d] = %d\n", i, a[i]);
  }*/
  /*preprocessing - column stored of L -> row stored*/
    for ( i = 0; i < n; i++ )
    {
        Lp[i] = 0;
        row_count_l[i] = 0;
        ptr_l_row[i] = 0;
    }
    Lp[n] = 0;
    for ( i = 0; i < nzl; i++ )
    {
        row_count_l[asub_L[i]]++;
    }
    for ( i = 1; i <= n; i++ )
    {
        Lp[i] = Lp[i - 1] + row_count_l[i - 1];
    }
    /*for ( i = 0; i <= n; i++ )
    {
        printf("Lp[%d] = %d\n", i, Lp[i]);
    }*/
    for ( i = 0; i < n; i++ )
    {
        for ( j = xa_L[i]; j < xa_L[i+1]; j++ )
        {
            row_num_l = asub_L[j];
      //      printf("i = %d j = %d xa_L[%d] = %d xa_L[%d] = %d\n", i, j, i, xa_L[i], i+1, xa_L[i+1]);
            Li[Lp[row_num_l]+ptr_l_row[row_num_l]] = i;
            ptr_l_row[row_num_l]++;
        }
    }
    /*for ( i = 0; i < nzl; i++ )
    {
        printf("Li[%d] = %d\n", i, Li[i]);
    }*/

  /* Array xx initialization*/
    for ( i = 0; i < n; i++ )
    {
        xx[i] = 0;
        ptr_l_row[i] = 0;
  //      L[Lp[i+1]-1] = 1.0;
    }

  int row_column_U, row_column_L, p_u, p_l, row_l, row_u;
  double accle = 0;
  double com_count = 0, if1_count = 0, if2_count = 0;
  double start, end, time1, time2, sum_time = 0;
  double if1_time_start, if2_time_start, if1_time_end, if2_tine_end, if1_sum_time = 0, if2_sum_time = 0;
  /* verify the correctness */
  /* column-oriented G/P algorithm without partial pivoting */
//   start = microtime();
  for ( k = 0; k < n; k++ )
  {
    current_column = perm_c[k];
    
    /* xx[] = A[:,current_column] */
    for ( j = xa[current_column]; j < xa[current_column+1]; j++ )
    {
      xx[perm_r[asub[j]]] = a[j];
    }
    /*printf("*************current_column is: %d**************\n", k);
    for ( i = 0; i < n; i++ )
    {
        printf("xx[%d] = %lf\n", i, xx[i]);
    }*/
    row_column_L = xa_L[k+1] - xa_L[k] - 1;
    row_column_U = xa_U[k+1] - xa_U[k];
    row_column = row_column_L + row_column_U;
    row_index = (int *)malloc( sizeof(int *) * row_column);
    for ( j = 0; j < row_column_U; j++ )
    {
        row_index[j] = asub_U[j+xa_U[k]];
    }
    for ( j = row_column_U; j < row_column; j++ )
    {
        row_index[j] = asub_L[j+xa_L[k]-row_column_U+1];
    }
    
    //time1 = microtime();
    for ( j = 0; j < row_column; j++ )
    {
        row = row_index[j];
      //  printf("k = %d BEGIN: xx[%d] = %lf\n", k, row, xx[row]);
        p_l = 0;
        p_u = 0;
        row_l = Lp[row+1] - Lp[row];
        row_u = xa_U[k+1] - xa_U[k];
        //time1 = microtime();
        while ( p_l < row_l && p_u < row_u )
        {
            if ( Li[Lp[row]+p_l] == asub_U[xa_U[k]+p_u] ) 
            {
                //com_count++;
                 accle += L[Lp[row]+p_l] * xx[asub_U[xa_U[k]+p_u]];
                 p_l++;
                 p_u++;
            }
//             if1_time_start = microtime();
            if ( Li[Lp[row]+p_l] < asub_U[xa_U[k]+p_u] ) 
            {
                //if1_count++;
                 p_l++;
            }
//             if1_time_end = microtime() - if1_time_start;
//             if1_sum_time += if1_time_end;

            //if2_time_start = microtime();
            if ( Li[Lp[row]+p_l] > asub_U[xa_U[k]+p_u] ) 
            {
                //if2_count++;
                 p_u++;
            } 
            //if2_tine_end = microtime() - if2_time_start;
            //if2_sum_time += if2_tine_end;
        }
        //time2 = microtime() - time1;
        //sum_time = sum_time + time2;
        xx[row] = xx[row] - accle;
//        printf("k = %d AFTER: xx[%d] = %lf\n", k, row, xx[row]);
        accle = 0; 
    }
    //time2 = microtime() - time1;
    //sum_time = sum_time + time2;

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
      L[Lp[asub_L[i]]+ptr_l_row[asub_L[i]]] = xx[asub_L[i]] / U_diag;
      ptr_l_row[asub_L[i]]++;
      xx[asub_L[i]] = 0;
    }
  }
//   end = microtime() - start;
//   printf("GP time is: %lf\n", end);
  //printf("while time is: %lf\n", sum_time);
//   printf("if1 time is: %lf\n", if1_sum_time);
  //printf("if2 time is: %lf\n", if2_sum_time);
  //printf("计算次数为: %lf\n", com_count);
  //printf("if1次数为: %lf\n", if1_count);
  //printf("if2次数为: %lf\n", if2_count);
  for ( i = 0; i < n; i++ )
  {
     L[Lp[i+1]-1] = 1.0;
  }

  /* solve for Ly = b and Ux = y */
//   y[0] = 1.0;
//   double y_sum = 0;
//   for ( i = 1; i < n; i++ )
//   {
//     for ( j = Lp[i]; j < Lp[i+1]-1; j++ )
//     {
//         y_sum += L[j] * y[Li[j]];
//     }
//     y[i] = 1 - y_sum;
//     y_sum = 0;
//   }
//   y[n-1] = y[n-1]/U[nzu-1];
//   for ( i = n-1; i > 0; i-- )
//   {
//     for ( j = xa_U[i]; j < xa_U[i+1]-1; j++ )
//     {
//         y[asub_U[j]] -= U[j] * y[i];
//     }
//     y[i-1] = y[i-1]/U[xa_U[i]-1];
//   }
// //   /*printf("**********************************The solution is: \n");
//   for ( i = 0; i < n; i++ )
//   {
//     //x[i] = y[i]/U[xa_U[i+1]-1];
//     printf("x[%d] = %lf\n", i, y[i]);
//   }
// 
  /* Check L value*/
//   for ( i = 0; i < nzl; i++ )
//   {
//       /*if ( equal(L[i], l_data[i]) ) continue;
//       else 
//       {
//           sum_l++; 
//           fprint f(fp_L_Error, "Correct: %lf Error: %lf\n", l_data[i], L[i]);
//           
//       }*/
//       fprintf(fp_L_Error, "L[%d] =  %lf\n", i, L[i]);
// //       printf("L[%d] = %lf\n", i, L[i]);
//   }
//   
//   
//   /* Check U value*/
//   for ( i = 0; i < nzu; i++ )
//   {
//       /*if ( equal(u_data [i], U[i]) ) continue;
//       else
//       {
//           sum_u++;
//           fprintf(fp_U_Error, "Correct: %lf Error: %lf\n", u_data[i], U[i]);
//           
//       }*/
//       fprintf(fp_U_Error, "U[%d] =  %lf\n", i, U[i]);
// //       printf("U[%d] = %lf\n", i, U[i]);
//   }


  /* Print check results information*/
  /*if ( sum_l == 0 && sum_u == 0 ) printf("Correct Results!\n");
  else 
  {
      printf("Num of Errors in L: %d\n", sum_l);
      printf("Num of Errors in U: %d\n", sum_u);
      printf("Notice: For Error Results, please check the detailed information in Error_Result directory!\n");
  }*/
  
  
  /* Print Time Information */
 // printf("\n********************MY LU Decomposition: Time Statisistics********************\n");
//   printf("Calculation Time: %lf\n", sum_cal);
//   printf("Calculation Counts: %d\n", sum_calculation);
   //   printf("Counts of OUTER-FOR LOOP: %d\n", sum_for);
  //    //printf("aaa = %d\n", aaa);
  
}

