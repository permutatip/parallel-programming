#include<stdio.h>
#include<assert.h>

typedef unsigned long long ull;

//----linux timer----
#ifdef __linux
#include<sys/time.h>
struct timeval tm_start,tm_end;
#define set_time(x) gettimeofday(&x,0)
double get_time()
{
    suseconds_t t1=tm_end.tv_usec-tm_start.tv_usec;
    time_t t2=tm_end.tv_sec-tm_start.tv_sec;
    return t1*0.001 + t2*1000;
}
#endif

//----windows timer----
#ifdef _WIN32
#include<windows.h>
LARGE_INTEGER tm_start,tm_end,tm_freq;
#define set_time(x) QueryPerformanceCounter(&x)
double get_time()
{
    QueryPerformanceFrequency(&tm_freq);
    long long t1=tm_end.QuadPart-tm_start.QuadPart;
    long long t2=tm_freq.QuadPart;
    double t=(double)t1/t2*1000;
    return t;
}
#endif
// no timer
#ifndef set_time
#error "Not linux or windows, can not excute."
#endif

ull queen_bit(int n,int row,ull col,ull d1,ull d2)
{
    if(row==n) return 1;
    ull count=0;
    ull ava=((1<<n)-1)&(~(col|d1|d2));
    while(ava)
    {
        ull p=ava&(-ava);
        ava=ava-p;
        count+=queen_bit(n,row+1,col|p,((1<<n)-1)&((d1|p)<<1),(d2|p)>>1);
    }
    return count;
}

int main()
{
    for(int size=11;size<=16;size++)
    {
        set_time(tm_start);
        ull cnt=queen_bit(size,0,0,0,0);
        set_time(tm_end);
        double t=get_time();
        printf("size:%d ans:%15llu time:%20.3f\n",size,cnt,t);
    }
    return 0;
}