#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <windows.h>
#include <pmmintrin.h>
#include <omp.h>

using namespace std;

const int N = 2000;

float A[N][N];
float test[N][N];

const int thread_count = 4;
long long head, tail, freq;        // timers

void init(float test[][N])
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

void reset(float A[][N], float test[][N])
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = test[i][j];
}

void simple(float A[][N])
{
	for (int k = 0; k < N; k++)
	{
		for (int j = k + 1; j < N; j++)
			A[k][j] = A[k][j] / A[k][k];
		A[k][k] = 1.0;
		for (int i = k + 1; i < N; i++)
		{
			for (int j = k + 1; j < N; j++)
				A[i][j] = A[i][j] - A[i][k] * A[k][j];
			A[i][k] = 0;
		}
	}
}

void omp(float A[][N])
{
    #pragma omp parallel num_threads(thread_count)
	for (int k = 0; k < N; k++)
	{
	    #pragma omp for schedule(static)
		for (int j = k + 1; j < N; j++)
		{
		    A[k][j] = A[k][j] / A[k][k];
		}
		A[k][k] = 1.0;
		#pragma omp for schedule(static)
		for (int i = k + 1; i < N; i++)
		{
			for (int j = k + 1; j < N; j++)
				A[i][j] = A[i][j] - A[i][k] * A[k][j];
			A[i][k] = 0;
		}
	}
}

void omp_dynamic(float A[][N])
{
    #pragma omp parallel num_threads(thread_count)
	for (int k = 0; k < N; k++)
	{
	    #pragma omp for schedule(dynamic, 24)
		for (int j = k + 1; j < N; j++)
        {
            A[k][j] = A[k][j] / A[k][k];
        }
		A[k][k] = 1.0;
		#pragma omp for schedule(dynamic, 24)
		for (int i = k + 1; i < N; i++)
		{
			for (int j = k + 1; j < N; j++)
				A[i][j] = A[i][j] - A[i][k] * A[k][j];
			A[i][k] = 0;
		}
	}
}


void omp_sse_dynamic(float A[][N])
{
    __m128 t1, t2, t3;
    #pragma omp parallel for shared(A) private(t1, t2, t3) schedule(dynamic, 24)
    for (int k = 0; k < N; k++)
    {
        float temp1[4] = {A[k][k], A[k][k], A[k][k], A[k][k]};
        t1 = _mm_loadu_ps(temp1);

        for (int j = k + 1; j < N - 3; j += 4)
        {
            t2 = _mm_loadu_ps(A[k] + j);
            t3 = _mm_div_ps(t2, t1);
            _mm_storeu_ps(A[k] + j, t3);
        }
        for (int j = N - (N % 4); j < N; j++)
        {
            A[k][j] = A[k][j] / A[k][k];
        }

        A[k][k] = 1.0;
        #pragma omp parallel for schedule(dynamic, 24)
        for (int i = k + 1; i < N; i++)
        {
            float temp2[4] = {A[i][k], A[i][k], A[i][k], A[i][k]};
            t1 = _mm_loadu_ps(temp2);

            for (int j = k + 1; j < N - 3; j += 4)
            {
                t2 = _mm_loadu_ps(A[i] + j);
                t3 = _mm_loadu_ps(A[k] + j);
                t3 = _mm_mul_ps(t1, t3);
                t2 = _mm_sub_ps(t2, t3);
                _mm_storeu_ps(A[i] + j, t2);
            }

            for (int j = N - (N % 4); j < N; j++)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }

            A[i][k] = 0;
        }
    }
}


int main()
{
    init(test);
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq);	// similar to CLOCKS_PER_SEC

	reset(A, test);
	QueryPerformanceCounter((LARGE_INTEGER *)&head);	// start time
	simple(A);
	QueryPerformanceCounter((LARGE_INTEGER *)&tail);	// end time
	cout << "simple: " << (tail - head) * 1000.0 / freq << "ms" << endl;


	reset(A, test);
	QueryPerformanceCounter((LARGE_INTEGER *)&head);	// start time
	omp(A);
	QueryPerformanceCounter((LARGE_INTEGER *)&tail);	// end time
	cout << "blocked OpenMP: " << (tail - head) * 1000.0 / freq << "ms" << endl;

	reset(A, test);
	QueryPerformanceCounter((LARGE_INTEGER *)&head);	// start time
	omp_dynamic(A);
	QueryPerformanceCounter((LARGE_INTEGER *)&tail);	// end time
	cout << "dynamic OpenMP: " << (tail - head) * 1000.0 / freq << "ms" << endl;


	reset(A, test);
	QueryPerformanceCounter((LARGE_INTEGER *)&head);	// start time
	omp_sse_dynamic(A);
	QueryPerformanceCounter((LARGE_INTEGER *)&tail);	// end time
	cout << "sse_dynamic OpenMP: " << (tail - head) * 1000.0 / freq << "ms" << endl;

	cout << endl;
    return 0;
}
