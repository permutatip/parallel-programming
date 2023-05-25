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

#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>

typedef unsigned long long ull;
ull queen(int n,int row,ull col,ull d1,ull d2);
ull queen_cal(int n);
ull qhalf(int n);
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
ull queen_cal(int n)
{
    return queen(n,0,0,0,0);
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
14772512,95815104,666090624,4968057848,//16~19
39029188884,314666222712//20,21
};
int main(int argc,char* argv[])
{
    int myid,numprocs;
    // int namelen;
    // char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    // MPI_Get_processor_name(name,&namelen);
    int size,beg_id,end_id;
    for(size=8;size<=9;size++)
    {
        int size_per_id=size/numprocs;
        beg_id=size_per_id*myid;
        end_id=(myid==numprocs-1)?(size-1):(beg_id+size_per_id-1);
        ull task_ans=0;
        for(int j=beg_id;j<=end_id;j++)
        {
            ull col_arg=1<<j;
            ull d1_arg=1<<(j+1);
            ull d2_arg=(j==0)?0:j-1;
            task_ans+=queen(size,1,col_arg,d1_arg,d2_arg);
        }
        printf("size:%d ans:%llu from id%d begin%d end%d\n",
        size,task_ans,myid,beg_id,end_id);
        ull* ans_buffer=(ull*)malloc(sizeof(ull)*numprocs);
        if(myid==0)
        {
            ull final_ans=task_ans;
            for(int i=1;i<numprocs;i++)
            {
                MPI_Recv(&task_ans,1,MPI_UNSIGNED_LONG_LONG,
                i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                final_ans+=task_ans;
            }
            printf("final ans:%18llu of size %d\n",final_ans,size);
        }
        else
        {
            MPI_Send(&task_ans,1,MPI_UNSIGNED_LONG_LONG,
            0,0,MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}