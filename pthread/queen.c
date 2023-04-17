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
#include<pthread.h>
typedef unsigned long long ull;
ull queen(int n,int row,ull col,ull d1,ull d2);
typedef struct task_arg
{
    int n;
    int row;
    ull col,d1,d2;
}task_arg;
typedef struct thrd_arg
{
    int thrd_id;
    int beg_task_id,end_task_id;
}thrd_arg;
pthread_mutex_t mu;
task_arg* all_task;
int task_cnt;
ull global_cnt;
void* thrd_run(void* arg)
{
    thrd_arg t_arg=*(thrd_arg*)arg;
    ull ans=0;
    for(int i=t_arg.beg_task_id;i<=t_arg.end_task_id;i++)
    {
        ans+=queen(all_task[i].n,all_task[i].row,all_task[i].col,
        all_task[i].d1,all_task[i].d2);
    }
    pthread_mutex_lock(&mu);
    printf("thrd %d ans:%llu\n",t_arg.thrd_id,ans);
    global_cnt+=ans;
    pthread_mutex_unlock(&mu);
    return NULL;
}

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

int main(int argc,char*argv[])
{
    assert(argc==2);
    int thrd_cnt=atoi(argv[1]);
    for(int size=13;size<=13;size++)
    {
        set_time(tm_start);
        ull cnt=queen(size,0,0,0,0);
        set_time(tm_end);
        double t=get_time();
        printf("  common size:%d ans:%15llu time:%20.3fms\n",size,cnt,t);

        //thread method
        pthread_mutex_init(&mu,0);
        task_cnt=size*size-3*size+2;
        all_task=malloc(sizeof(task_arg)*task_cnt);
        global_cnt=0;
        int id=0;
        //alloc task
        for(int i=0;i<size;i++)
        {
            for(int j=0;j<size;j++)
            {
                if(i==j||i-j==1||i-j==-1) continue;
                all_task[id].n=size;
                all_task[id].row=2;
                all_task[id].col=(ull)((1<<i)|(1<<j));
                all_task[id].d1=(ull)((1<<(i+2))|(1<<(j+1)));
                all_task[id].d2=(ull)((1<<(i-2))|(1<<(j-1)));
                id++;
            }
        }
        assert(id==task_cnt);
        //thread begin
        pthread_t* handle=malloc(sizeof(pthread_t)*thrd_cnt);
        int task_per_thrd=task_cnt/thrd_cnt;
        for(int i=1;i<=thrd_cnt;i++)
        {
            thrd_arg ta={.thrd_id=i,
                .beg_task_id=task_per_thrd*(i-1),
                .end_task_id=(i==thrd_cnt)?task_cnt-1:task_per_thrd*i-1};
            // printf("thrd %d: start task:%d, end task:%d\n",
            // ta.thrd_id,ta.beg_task_id,ta.end_task_id);
            pthread_create(&handle[i-1],NULL,thrd_run,&ta);
        }
        for(int i=0;i<thrd_cnt;i++)
            pthread_join(handle[i],NULL);

        printf("%dthread size:%d ans:%15llu\n",thrd_cnt,size,global_cnt);
        pthread_mutex_destroy(&mu);
        free(all_task);
        free(handle);
    }
    return 0;
}