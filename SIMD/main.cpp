#include<iostream>
#include<emmintrin.h>
#include<time.h>
#include<Windows.h>
#include<immintrin.h>
using namespace std;
const int N = 1024;
float A[N][N] = {0};

void init()
{
    for(int i=0;i<N;i++)
    {
        A[i][i]=1.0;
        for(int j=i+1;j<N;j++)
            A[i][j]=rand();
    }
    for(int k=0;k<N;k++)
        for(int i=k+1;i<N;i++)
            for(int j=0;j<N;j++)
                A[i][j]+=A[k][j];
}

//平凡算法
void simple()
{
    for (int k= 0; k < N; k++)
    {
        for (int j = k+1;j<N; j++)
            A[k][j] = A[k][j] / A[k][k];

        A[k][k] = 1.0;

        for (int i = k + 1; i < N; i++)
        {
            for (int j = k + 1; j < N; j++)
                A[i][j] = A[i][j] - A[i][k] * A[k][j];

            A[i][k] = 0;
        }
    }
    return;
}

//SSE
void sse_gauss()
{
    __m128 t0,t1,t2,t3;
    for(int k=0;k<N;k++)
    {
        float temp1[4]={A[k][k],A[k][k],A[k][k],A[k][k]};
        t0=_mm_loadu_ps(temp1);
        int j;
        for(j=k+1;j+3<N;j+=4)//j+4的值为下一次load操作的起始位置
        {
            t1=_mm_loadu_ps(A[k]+j);
            t2=_mm_div_ps(t1,t0);
            _mm_storeu_ps(A[k]+j,t2);
        }
        for(;j<N;j++)
            A[k][j]/=A[k][k];

        A[k][k]=1.0;
        for(int i=k+1;i<N;i++)
        {
            float temp2[4]={A[i][k],A[i][k],A[i][k],A[i][k]};
            t0=_mm_loadu_ps(temp2);
            int j;
            for(j=k+1;j+3<N;j+=4)
            {
                t1=_mm_loadu_ps(A[k]+j);
                t2=_mm_loadu_ps(A[i]+j);
                t3=_mm_mul_ps(t0,t1);
                t2=_mm_sub_ps(t2,t3);
                _mm_storeu_ps(A[i]+j,t2);
            }
            for(;j<N;j++)
                 A[i][j]-=A[i][k]*A[k][j];
            A[i][k]=0.0;
        }
    }
}

//AVX
void avx_gauss()
{
    __m256 t1, t2, t3;
    for (int k = 0; k < N; k++)
    {
        float temp1[8] = { A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k] };
        t1 = _mm256_loadu_ps(temp1);
        int j = k + 1;
        for (j; j < N - 7; j += 8)
        {
            t2 = _mm256_loadu_ps(A[k] + j);
            t3 = _mm256_div_ps(t2, t1);
            _mm256_storeu_ps(A[k] + j, t3);
        }
        for (j; j < N; j++)
            A[k][j] = A[k][j] / A[k][k];

        A[k][k] = 1.0;

        for (int i = k + 1; i < N; i++)
        {
            float temp2[8] = { A[i][k], A[i][k], A[i][k], A[i][k], A[i][k], A[i][k], A[i][k], A[i][k] };
            t1 = _mm256_loadu_ps(temp2);
            j = k + 1;
            for (j; j < N - 7; j += 8)
            {
                t2 = _mm256_loadu_ps(A[i] + j);
                t3 = _mm256_loadu_ps(A[k] + j);
                t3 = _mm256_mul_ps(t1, t3);
                t2 = _mm256_sub_ps(t2, t3);
                _mm256_storeu_ps(A[i] + j, t2);
            }
            for (j; j < N; j++)
                A[i][j] = A[i][j] - A[i][k] * A[k][j];

            A[i][k] = 0;
        }
    }
}

int main()
{
    long long head1, tail1 , freq1;
    init();
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq1 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head1);
    simple();
    QueryPerformanceCounter((LARGE_INTEGER *)&tail1 );
    cout<<"simple:"<<(tail1-head1)*1000.0/freq1<<"ms"<<endl;

    long long head2, tail2 , freq2;
    init();
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq2 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head2);
    sse_gauss();
    QueryPerformanceCounter((LARGE_INTEGER *)&tail2 );
    cout << "sse:" <<(tail2 - head2) * 1000.0/freq2<<"ms"<<endl;

    long long head3, tail3 , freq3;
    init();
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq3 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head3);
    avx_gauss();
    QueryPerformanceCounter((LARGE_INTEGER *)&tail3 );
    cout << "avx:" <<(tail3 - head3) * 1000.0/freq3<<"ms"<<endl;

    return 0;
}

