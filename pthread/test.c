#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<unistd.h>
#include<time.h>
#include<semaphore.h>
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
sem_t s_parent;
sem_t s_child;
int mutex_sum=0;

void* hello(void* thrd_arg)
{
    int thrd_num=*(int*)thrd_arg;
    printf("hello1 from thread%d\n",thrd_num);
    sem_post(&s_parent);
    sem_wait(&s_child);
    printf("hello2 from thread%d\n",thrd_num);
    return NULL;
}

int main()
{
    pthread_t* handle=malloc(THRD_CNT*sizeof(pthread_t));
    sem_init(&s_parent,0,0);
    sem_init(&s_child,0,0);
    for(int i=0;i<THRD_CNT;i++) pthread_create(&handle[i],NULL,hello,(void*)&i);
    for(int i=0;i<THRD_CNT;i++) sem_wait(&s_parent);
    printf("hello from main1\n");
    for(int i=0;i<THRD_CNT;i++) sem_post(&s_child);
    for(int i=0;i<THRD_CNT;i++) pthread_join(handle[i],NULL);
    printf("hello from main2\n");
    sem_destroy(&s_child);
    sem_destroy(&s_parent);
    free(handle);
    return 0;
}