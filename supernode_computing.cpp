#include "Include/lu.h"
#include <stdio.h>
#include <stdlib.h>
# include <cs_demo.h>

int min ( int a, int b )
{
    if ( a > b ) return b;
    return a;
}

int detect ( int *asub,int *xa, int lower_col, int higher_col )
{
    int lower_col_ptr = xa[lower_col + 1] - 1;
    int higher_col_ptr = xa[higher_col + 1] - 1;
    int lower_col_count = xa[lower_col + 1] - xa[lower_col];
    int higher_col_count = xa[higher_col + 1] - xa[higher_col];
    int count = min(lower_col_count, higher_col_count) - 1;
    int i;

    for ( i = 0; i < count; i++ )
    {
        if ( asub[lower_col_ptr] == asub[higher_col_ptr] )
        {
            lower_col_ptr--;
            higher_col_ptr--;
             continue;
        }
        else
        {
          break;  
         }
    }

    if ( i == count )
    {
        /* relaxed supernode */
        /*if ( asub[lower_col_ptr] == lower_col || asub[lower_col_ptr] == higher_col )
        {
            if ( asub[higher_col_ptr] == higher_col ) return 1;
             return 0;
        }*/
        /* strict supernode */
        if ( asub[lower_col_ptr] == higher_col && asub[higher_col_ptr] == higher_col ) return 1;
         return 0;
    }
    return 0;
}
int detect_nontri ( int *asub,int *xa, int lower_col, int higher_col, int length )
{
    int lower_col_ptr = xa[lower_col + 1] - 1;
    int higher_col_ptr = xa[higher_col + 1] - 1;
    int i;

    for ( i = 0; i < length; i++ )
    {
        if ( asub[lower_col_ptr] == asub[higher_col_ptr] )
        {
            lower_col_ptr--;
            higher_col_ptr--;
            continue;
        }
        else
        {
          break;  
        }
    }

    if ( i == length )
    {
        return 1;
    }
    return 0;
}
/*int detect_U ( int *asub,int *xa, int lower_col, int higher_col )
{
    in t lower_col_ptr = xa[lower_col];
    int higher_col_ptr = xa[higher_col];
    int lower_col_count = xa[lower_col + 1] - xa[lower_col];
    int higher_col_count = xa[higher_col + 1] - xa[higher_col];
    int count = min(lower_col_count, higher_col_count) - 1;
    int i;

    for ( i = 0; i < count; i++ )
    {
        if ( asub[lower_col_ptr] == asub[higher_col_ptr] )
        {
            lower_col_ptr++;
            higher_col_ptr++;
            continue;
        }
        else
        {
          break;  
        }
    }

    if ( i == count )
    {*/
        /* relaxed supernode */
        /*if ( asub[lower_col_ptr] == lower_col || asub[lower_col_ptr] == higher_col )
        {
            if ( asub[higher_col_ptr] == higher_col ) return 1;
            return 0;
        }*/
        /* strict supernode */
//         if ( asub[lower_col_ptr] == lower_col && asub[higher_col_ptr] == lower_col ) return 1;
/*        if ( asub[lower_col_ptr] == lower_col && asub[higher_col_ptr] == higher_col ) return 1;
    }
    return 0;
}*/

/*int detect_similar_U ( int *asub,int *xa, int lower_col, int higher_col )
{
    int lower_col_ptr = xa[lower_col];
    int higher_col_ptr = xa[higher_col];
    int lower_col_count = xa[lower_col + 1] - xa[lower_col];
    int higher_col_count = xa[higher_col + 1] - xa[higher_col];
//     int count = min(lower_col_count, higher_col_count) - 1;
    int i;

    if ( lower_col_count == higher_col_count )
    {
        if ( lower_col_count <= 4 ) return 0;
        for ( i = 0; i < lower_col_count-1; i++ )
        {
            if ( asub[lower_col_ptr] == asub[higher_col_ptr] )
            {
                lower_col_ptr++;
                higher_col_ptr++;
                continue;
            }
            else
            {
                return 0;  
            }
        }
        return 1;
    }
    else
    {
        return 0;
    }
} */

bool equal( double a, double b )
{
  if ( Abs(a-b) < 0.001  )
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

int main( int argc, char *argv[] )
{
     /* Read the file1 and store  it by CSC: <a, asub, xa> */
    FILE *file1, *file2, *file3;
    FILE *file4, *file5;
    int *pinv, *q;
    int m, n, i, j = 0;
    int l_i, l_j, sum = 0, u_i, u_j, a_i, a_j;
    int nzl, nzu;
    int *Ai, *Ap, *Li, *Ui, *Lp, *Up;
    int *p, *p_u, *p_a;
    double *Ax, a_x;
    int nnz;
    double start, end, start_GP, end_GP, start_GP_sn, end_GP_sn, start_GP_sn_u, end_GP_sn_u;
    int *sn_record, *sn_record_u;
    file1 = fopen(argv[1], "r"); 
    file2 = fopen(argv[2], "r");
    file3 = fopen(argv[3], "r");
    file4 = fopen(argv[4], "r");
    file5 = fopen(argv[5], "r");
    int row_thresold = atoi(argv[6]);
    int column_thresold = atoi(argv[7]);
 //   nzl = atoi(argv[6]);
 //   nzu = atoi(argv[7]);
//     
//     printf("nzl = %d nzu = %d\n", nzl, nzu );
//     fscanf ( file1, "%d %d %lf\n", &m, &n, &nnz);
    char ch[100];
    do
    {
       // fscanf ( file1, "%[^\n]%*c", &ch);
       // fscanf ( file1, "%c", &ch );
       // getline(&ch, &len, file1);
       fgets(ch, 100, file1);
    }while ( ch[0] == '%' );
    sscanf ( ch, "%d %d %d\n", &m, &n, &nnz);
fscanf( file2, "%d\n", &nzl);
fscanf( file3, "%d\n", &nzu);
//     fscanf( file2, "%d\n", &nzl);

    Ap = ( int * )malloc( sizeof( int ) * n+1 ) ;
    Ai = ( int * )malloc( sizeof( int ) * nnz ) ;
    Ax = ( double * )malloc( sizeof( double ) * nnz ) ;
    p_a = ( int * )malloc( sizeof( int ) * n+1 ) ;
    Lp = ( int * )malloc( sizeof( int ) * n+1 ) ;
    Up = ( int * )malloc( sizeof( int ) * n+1 ) ;
    p = ( int * )malloc( sizeof( int ) * n ) ;
    p_u = ( int * )malloc( sizeof( int ) * n ) ;
    Li = ( int * )malloc( sizeof( int ) * nzl ) ;
    Ui = ( int * )malloc( sizeof( int ) * nzu ) ;
   sn_record = ( int * )malloc( sizeof( int ) * n );
   sn_record_u = ( int * )malloc( sizeof( int ) * n );
    memset(Lp, 0, sizeof(int)*(n+1));
    memset(Up, 0, sizeof(int)*(n+1));
    memset(p, 0, sizeof(int)*n);
    memset(p_u, 0, sizeof(int)*n);
    memset(Ap, 0, sizeof(int)*(n+1));
    memset(p_a, 0, sizeof(int)*n);
   memset(sn_record, -1, sizeof(int)*n);
   memset(sn_record_u, -1, sizeof(int)*n);

    while ( fscanf ( file1, "%d %d %lf\n", &a_i, &a_j, &a_x) == 3 )
    {
         Ai[sum] = a_i-1; 
        Ax[sum] = a_x;
        p_a[a_j-1]++;
        sum++ ;
    }
    for ( i = 1; i <=n; i++ )
    {
        Ap[i] = Ap[i-1]+p_a[i-1];
    }
    fclose(file1);

    sum = 0;
    while ( fscanf ( file2, "%d %d\n", &l_i, &l_j) == 2 )
    {
        Li[sum] = l_i; 
        p[l_j]++;
          sum++;
    }
    for ( i = 1; i <=n; i++ )
    {
        Lp[i] = Lp[i-1]+p[i-1];
    }
/*    for ( i = 0; i < n; i++ )
    {
        //printf("%d\n", p[i]);
        if ( p[i] > 8)  sum_l_greater_4++;
    }
    printf("row: %d greater4: %d\n", n, sum_l_greater_4);*/
    fclose(file2);

    sum = 0;
    while ( fscanf ( file3, "%d %d\n", &u_i, &u_j) == 2 )
    {
          Ui[sum] = u_i;
        p_u[u_j]++;
         sum++ ;
    }
    for ( i = 1; i <=n; i++ )
    {
        Up[i] = Up[i-1]+p_u[i-1];
    }
    fclose(file3);

    pinv = ( int * )malloc( sizeof( int ) * n ) ;
    while ( fscanf(file4, "%d\n", &i) == 1 )
    {
         pinv[j] = i;  
         j++; 
    }

    q = ( int * )malloc( sizeof( int ) * n ) ;
    j = 0;
    while ( fscanf(file5, "%d\n", &i) == 1 )
    {
          q[j] = i;  
        j++; 
    }
    
    
    /* similar columns detected in U */
//     int *similar_u;
//     similar_u = (int *)malloc( sizeof( int ) * n );
//     memset(similar_u, 0, sizeof(int)*n);
//     for ( i = 100; i < n; i++ )
//     {
//         if ( similar_u[i] == 0)
//         {
//             for ( j = i+1; j < n; j++ )
//             {
//                 if ( detect_similar_U(Ui, Up, i, j) )
//                 {
//                     similar_u[i]++;
//                     printf(" i = %d j = %d\n ", i, j);
//                     similar_u[j] = -1;
//                 }
//             }
//         }
//         
//     }
//     
//     for ( i = 100; i < n; i++ )
//     {
//         if ( similar_u[i] > 0 ) printf("similar_u[%d] = %d\n", i, similar_u[i]);
//     }
    
    
    
    
    
    
     /* supernode detect */
    int *sn_start = (int *)malloc( sizeof( int ) * n );
    int *sn_end = (int *)malloc( sizeof( int ) * n );
    memset(sn_start, 0, sizeof(int)*n);
    memset(sn_end, 0, sizeof(int)*n);
    int sn_sum = 0;
    int sn_sum_final = 0;

    int *sn_start_u = (int *)malloc( sizeof( int ) * n );
    int *sn_end_u = (int *)malloc( sizeof( int ) * n );
    memset(sn_start_u, 0, sizeof(int)*n);
    memset(sn_end_u, 0, sizeof(int)*n);
    int sn_sum_u = 0;
    int sn_sum_final_u = 0;
    
    for ( i = 0; i < n-1; i++ )
    {
        if ( p[i] == 1 || p[i+1] == 1 ) continue; 
        else
        { 
             if ( detect_nontri(Li, Lp, i, i+1, row_thresold) )
            {
                 sn_start[sn_sum] = i;
                sn_end[sn_sum] = i + 1;
                 sn_sum++;
            }
        } 
    }
    for ( i = sn_sum-1; i > 0; i-- )
    {
        if ( sn_start[i] == sn_end[i-1] )
        {
            sn_end[i-1] = sn_end[i];
            sn_start[i] = 0;
             sn_end[i] = 0;
        } 
    } 
    for ( i = 0; i < sn_sum; i++ )
    {
        if ( sn_end[i] == 0 ) continue;
        if ( sn_end[i] - sn_start[i] < column_thresold ) continue;
        else 
         {
            for ( j = sn_start[i]; j <= sn_end[i]; j++ )
            {
                sn_record[j] = sn_sum_final;
            }
             sn_sum_final++;
         }
    } 
    
//     for ( i = 0; i < n-1; i++ )
//     {
//         if ( p_u[i] == 1 || p_u[i+1] == 1 ) continue; 
//         else
//         { 
//              if ( detect_U(Ui, Up, i, i+1) )
//             {
//                  sn_start_u[sn_sum_u] = i;
//                 sn_end_u[sn_sum_u] = i + 1;
//                  sn_sum_u++;
//             }
//         } 
//     }
//     for ( i = sn_sum_u-1; i > 0; i-- )
//     {
//         if ( sn_start_u[i] == sn_end_u[i-1] )
//         {
//             sn_end_u[i-1] = sn_end_u[i];
//             sn_start_u[i] = 0;
//             sn_end_u[i] = 0;
//         } 
//     } 
//     for ( i = 0; i < sn_sum_u; i++ )
//     {
//         if ( sn_end_u[i] == 0 ) continue;
//         //if ( sn_end_u[i] - sn_start_u[i] == 1 ) continue;
//         else 
//          {
//             for ( j = sn_start_u[i]; j <= sn_end_u[i]; j++ )
//             {
//                 sn_record_u[j] = sn_sum_final_u;
//             }
//              sn_sum_final_u++;
//          }
//     } 
   
//     int *index = (int *)malloc( sizeof( int ) * sn_sum_final_u );
//     memset(index, 0, sizeof(int) * sn_sum_final_u);
//     for ( i = 0; i < n; i++ )
//     {
//         if ( sn_record_u[i] != -1 )
//         {
//             index[sn_record_u[i]]++;
//         }
//     }
    
    printf("defined number of supernodes of L are: %d\n", sn_sum_final);

    double *l_gp, *l_gp_sn, *x, *l_gp_sn_u;
    double sum_error1 = 0, sum_error2 = 0;;
    x = ( double * )malloc( sizeof(double) * n);   
    
     start_GP_sn = microtime();
     l_gp_sn = lu_gp_sparse_supernode_computing(Ax, Ai, Ap, n, nzl, nzu, q, pinv, Li, Lp, Ui, Up, sn_record);
     end_GP_sn = microtime() - start_GP_sn;
     printf("Time of LU_GP_SN: %lf\n", end_GP_sn);
    
     start_GP_sn_u = microtime();
     l_gp_sn_u = lu_gp_sparse_supernode_computing_blas(Ax, Ai, Ap, n, nzl, nzu, q, pinv, Li, Lp, Ui, Up, sn_record);
     end_GP_sn_u = microtime() - start_GP_sn_u;
     printf("Time of LU_GP_SN_blas: %lf\n", end_GP_sn_u);
    
    
    start_GP = microtime();
    l_gp = lu_gp_sparse(Ax, Ai, Ap, n, nzl, nzu, q, pinv, Li, Lp, Ui, Up);
    end_GP = microtime() - start_GP;
    printf("Time of LU_GP: %lf\n", end_GP);
 
//     start_GP_sn = microtime();
//     l_gp_sn = lu_gp_sparse_avx2(Ax, Ai, Ap, n, nzl, nzu, q, pinv, Li, Lp, Ui, Up);
//     end_GP_sn = microtime() - start_GP_sn;
//     printf("Time of LU_GP_avx2: %lf\n", end_GP_sn);
    
    for ( i = 0; i < n; i++ )
    {
        if ( !equal ( l_gp[i], l_gp_sn[i] ) ) 
        {
//             printf("result[%d] = %lf sn_u_result[%d] = %lf\n", i, l_gp[i], i, l_gp_sn[i]);
            sum_error1++;
        }
        
        if ( !equal ( l_gp[i], l_gp_sn_u[i] ) ) 
        {
//              printf("result[%d] = %lf sn_u_blas_result[%d] = %lf\n", i, l_gp[i], i, l_gp_sn_u[i]);
            sum_error2++;
        }
    }
//     
   printf("error results are: %lf %lf\n", sum_error1, sum_error2);
//     for ( i = 0; i < n; i++ )
//     {
//         x[q[i]] = l_gp[i];
//     }
// 
//     for ( i = 0; i < n; i++ )
//     {
//         printf("x[%d] = %lf\n", i, x[i]);
//     }
//     
    return 0;
}







