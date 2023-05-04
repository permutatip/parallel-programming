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
pthread_t* handle;
task_arg* all_task;
thrd_arg* all_args;
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

static const ull correct_ans[]={
1,1,0,0,2,10,4,40,92,//0~8
352,724,2680,14200,73712,365596,2279184,//9~15
14772512,95815104,666090624,4968057848//16~19
};

int main(int argc,char*argv[])
{
    assert(argc==2);
    int thrd_cnt=atoi(argv[1]);
    for(int size=18;size<=18;size++)
    {
        // set_time(tm_start);
        // ull cnt=queen(size,0,0,0,0);
        // set_time(tm_end);
        // double t0=get_time();
        // printf(" common size:%d ans:%15llu time:%20.3fms\n",size,cnt,t0);
        // assert((unsigned long)size<sizeof(correct_ans)/sizeof(correct_ans[0])&&correct_ans[size]==cnt);

        //thread method
        // set_time(tm_start);//thread alloc begin
        pthread_mutex_init(&mu,0);
        task_cnt=size*size-3*size+2;
        all_task=malloc(sizeof(task_arg)*task_cnt);
        all_args=malloc(sizeof(thrd_arg)*thrd_cnt);
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
        // set_time(tm_end);//thread alloc end
        
        set_time(tm_start);//thread excute begin
        handle=malloc(sizeof(pthread_t)*thrd_cnt);
        int task_per_thrd=task_cnt/thrd_cnt;
        for(int i=0;i<thrd_cnt;i++)
        {
            all_args[i].beg_task_id=task_per_thrd*i;
            all_args[i].end_task_id=(i==thrd_cnt-1)?task_cnt-1:task_per_thrd*(i+1)-1;
            all_args[i].thrd_id=i;
            pthread_create(&handle[i],NULL,thrd_run,&all_args[i]);
        }
        for(int i=0;i<thrd_cnt;i++) pthread_join(handle[i],NULL);
        
        pthread_mutex_destroy(&mu);
        free(all_task);
        free(all_args);
        free(handle);

        set_time(tm_end);//thread excute end
        double t2=get_time();
        printf("%dthread size:%d ans:%15llu time:%20.3fms\n",thrd_cnt,size,global_cnt,t2);
        assert((unsigned long)size<sizeof(correct_ans)/sizeof(correct_ans[0])&&correct_ans[size]==global_cnt);
    }
    return 0;
}