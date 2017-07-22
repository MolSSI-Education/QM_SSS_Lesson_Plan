#include <lawrap/blas.h>
#include <lawrap/lapack.h>

int main()
{
   double A[10][10], B[10][10], C[10][10], w[10];
   int i, j, info;

   for (i = 0;i < 10;i++)
   {
       for (j = 0;j < 10;j++)
       {
           A[i][j] = i+j;
           B[i][j] = 1;
       }
   }

   dgemm('N', 'T', 10, 10, 10, 2.0, A, 10, B, 10, 0.0, C, 10);
   info = dsyev('V', 'U', 10, A, 10, w);
   printf("info %d\n", info);

   return info;
}

