#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#define MAX_NUM 1000
#define MAX_SIZE 10000

int n;
int arr[MAX_SIZE],result[MAX_SIZE];
int result2[MAX_SIZE];
int mat[MAX_SIZE][MAX_SIZE];
//----linux timer----
#ifdef __linux
#include<sys/time.h>
struct timeval tm_start,tm_end;
#define set_time(x) gettimeofday(&x,0)
double get_time()
{
    suseconds_t t1=tm_end.tv_usec-tm_start.tv_usec;
    time_t t2=tm_end.tv_sec-tm_start.tv_sec;
    return t1*0.001 + t2*1000.0;
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
    double t=(double)t1/t2*1000.0;
    return t;
}
#endif

void traditional_mul()
{
    set_time(tm_start);
    for(int i=0;i<n;i++)
    {
        result[i]=0;
        for(int j=0;j<n;j++)
        {
            result[i]+=mat[j][i]*arr[j];
        }
    }
    set_time(tm_end);
    printf("traditional time:%fms\n",get_time());
}

void cache_opt_mul()
{
    set_time(tm_start);
    for(int i=0;i<n;i++)
    {
        result2[i]=0;
    }
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n;i++)
        {
            result2[i]+=mat[j][i]*arr[j];
        }
    }
    set_time(tm_end);
    printf("cache opt time:%fms\n",get_time());
}



int main(int argc,char** argv)
{
    srand(time(0));
    if(argc<2)
    {
        printf("usage:%s pro_size\n",argv[0]);
        return 0;
    }
    n=atoi(argv[1]);
    if(n>MAX_SIZE||n<1)
    {
        printf("problem size between 1 and %d\n",MAX_SIZE);
        return 0;
    }
    for(int i=0;i<n;i++)
    {
        arr[i]=rand()%MAX_NUM;
        for(int j=0;j<n;j++)
        {
            mat[i][j]=rand()%MAX_NUM;
        }
    }
    
    traditional_mul();
    cache_opt_mul();

    for(int i=0;i<n;i++)
    {
        if(result[i]!=result2[i])
        {
            printf("result diff at %d\n",i);
        }
    }

    return 0;
}