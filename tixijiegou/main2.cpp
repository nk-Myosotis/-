#include<iostream>
#include<windows.h>

using namespace std;
const int n =100000;
int a[n];
int sum=0;
void init()
{
    sum=0;
    for(int i=0;i<n;i++)
    {
        a[i]=i;
    }
}
//平凡算法
void trivialAlgorithm()
{
    sum=0;
    for(int i=0;i<n;i++)
    {
        sum+=a[i];
    }
}
//多链路式优化算法
void optimizationAlgorithm()
{
    int sum1 = 0;
    int sum2 = 0;
    for (int i = 0;i < n; i += 2)
    {
        sum1 += a[i];
        sum2 += a[i + 1];
    }
    sum = sum1 + sum2;
}
int main()
{
    long long head1, tail1 , freq1;
    init();
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq1 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head1);
    for(int i=0;i<1000;i++)
    {
        trivialAlgorithm();
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail1 );
    cout <<"trivial:"<< (tail1-head1) *1.0/freq1<<"ms"<<endl;

    long long head2, tail2 , freq2;
    init();
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq2 );
    QueryPerformanceCounter((LARGE_INTEGER *)&head2);
    for(int i=0;i<1000;i++)
    {
        optimizationAlgorithm();
    }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail2);
    cout <<"optimization:"<< (tail2-head2) *1.0/freq2<<"ms"<<endl;
    return 0;
}
