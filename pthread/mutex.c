#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<time.h>
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

#define THRD_CNT 4
int mutex_sum=0;
pthread_mutex_t mu;

void* hello(void* thrd_arg)
{
    int thrd_num=*(int*)thrd_arg;
    pthread_mutex_trylock(&mu);
    printf("hello from thread%d\n",thrd_num);
    mutex_sum+=thrd_num;
    pthread_mutex_unlock(&mu);
    return NULL;
}

int main()
{
    pthread_t* handle=malloc(THRD_CNT*sizeof(pthread_t));
    pthread_mutex_init(&mu,0);
    for(int i=0;i<THRD_CNT;i++) pthread_create(handle+i,NULL,hello,&i);
    for(int i=0;i<THRD_CNT;i++) pthread_join(handle[i],NULL);
    printf("ans=%d\n",mutex_sum);
    pthread_mutex_destroy(&mu);
    free(handle);
    return 0;
}