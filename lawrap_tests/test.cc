#include <vector>
#include <cstdio>

#include <lawrap/blas.h>
#include <lawrap/lapack.h>

int main()
{
    std::vector<double> w(10);
    std::vector<double> A(100), B(100), C(100);
 
    for (int i = 0;i < 10;i++)
    {
        for (int j = 0;j < 10;j++)
        {
            A[i*10 + j] = i+j;
            B[i*10 + j] = 1;
        }
    }

    LAWrap::gemm('N', 'T', 10, 10, 10,
                 2.0, A.data(), 10,
                      B.data(), 10,
                 0.0, C.data(), 10);
   
    int info = LAWrap::heev('V', 'U', 10, A.data(), 10, w.data());
    printf("info %d\n", info);

    return info;
}
