#include <iostream>
#include <windows.h>
using namespace std;
const int n = 1000;
int B[n][n]={0};
int A[n]={0};
int sum[n]={0};
void init()
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            B[i][j]=i+j;
        }
    }
    for(int i=0;i<n;i++)
    {
        A[i]=i*10;
    }
}
void col_func()
{
    for(int i = 0; i < n; i++)
    {
        sum[i] = 0;
        for(int j = 0; j < n; j++)
        {
            sum[i] += B[j][i]*A[j];
        }
    }
}
void row_func()
{
    for(int i=0;i<n;i++)
    {
        sum[i]=0;
    }
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n;i++)
        {
            sum[i]+=B[j][i]*A[j];
        }

    }
}
int main()
{
    long long head1, tail1 , freq1;
    init();
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq1 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head1);
    for(int i=0;i<1000;i++)
    {
        col_func();
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail1 );
    cout<<"Col:"<<(tail1-head1)*1.0/freq1<<"ms"<<endl;

    long long head2, tail2 , freq2;
    init();
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq2 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head2);
    for(int i=0;i<1000;i++)
    {
        row_func();
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail2 );
    cout << "Row:" <<(tail2 - head2) * 1.0/freq2<<"ms"<<endl;
    return 0;
}
