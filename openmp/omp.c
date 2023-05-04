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

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
typedef unsigned long long ull;
ull queen(int n,int row,ull col,ull d1,ull d2);
ull qhalf(int n);

int task_cnt;

ull queen(int n,int row,ull col,ull d1,ull d2)
{
    if(row==n) return 1;
    ull count=0;
    ull ava=((1<<n)-1)&(~(col|d1|d2));
    while(ava)
    {
        ull p=ava&(-ava);
        ava=ava-p;
        count+=queen(n,row+1,col|p,((1<<n)-1)&((d1|p)<<1),(d2|p)>>1);
    }
    return count;
}

ull qhalf(int n)
{
    ull count=0;
    int end_id=n/2-1;
    for(int i=0;i<=end_id;i++)
    {
        count+=2*queen(n,1,1<<i,1<<(i+1),(i>0)?1<<(i-1):0);
    }
    if(n%2) count+=queen(n,1,1<<(n/2),1<<(n/2+1),1<<(n/2-1));
    return count;
}

static const ull correct_ans[]={
1,1,0,0,2,10,4,40,92,//0~8
352,724,2680,14200,73712,365596,2279184,//9~15
14772512,95815104,666090624,4968057848//16~19
};

int main(int argc,char*argv[])
{
    assert(argc==2);
    int thrd_cnt=atoi(argv[1]);
    for(int size=11;size<=16;size++)
    {
        set_time(tm_start);
        ull cnt=qhalf(size);
        set_time(tm_end);
        double t0=get_time();
        printf("common  size:%d ans:%15llu time:%20.3fms\n",size,cnt,t0);
        assert((unsigned long)size<sizeof(correct_ans)/sizeof(correct_ans[0])&&cnt==correct_ans[size]);

        //thread method
        set_time(tm_start);//thread alloc begin
        //alloc task
    }
    return 0;
}