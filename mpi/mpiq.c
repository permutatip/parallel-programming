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
#include<assert.h>
typedef unsigned long long ull;


ull shift(int t,int size)
{
    if(t<0||t>=size) return 0;
    return 1<<t;
}
ull queen(int n,int row,ull col,ull d1,ull d2);
ull serial_queen(int n);
ull serial_half_queen(int n);

ull coarse_grained_mpi_queen(int size,int myid,int numprocs)
{
    int size_per_id=(size+1)/numprocs;
    int beg_id=size_per_id*myid;
    int end_id=(myid==numprocs-1)?(size-1):(beg_id+size_per_id-1);
    ull task_ans=0;
    for(int j=beg_id;j<=end_id;j++)
    {
        ull col_arg=shift(j,size);
        ull d1_arg=shift(j+1,size);
        ull d2_arg=shift(j-1,size);
        task_ans+=queen(size,1,col_arg,d1_arg,d2_arg);
    }
    ull final_ans=task_ans;
    ull token=0;
    if(myid==0)
    {
        for(int i=1;i<numprocs;i++)
        {
            MPI_Recv(&token,1,MPI_UNSIGNED_LONG_LONG,
            i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            final_ans+=token;
        }
        // printf("final coarse ans:%18llu of size %d\n",final_ans,size);
    }
    else
    {
        //send subtask ans to task 0
        MPI_Send(&task_ans,1,MPI_UNSIGNED_LONG_LONG,
        0,0,MPI_COMM_WORLD);
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    return final_ans;
}

ull fine_grained_mpi_queen(int size,int myid,int numprocs)
{
    int num_of_task=size*size-3*size+2;
    int size_per_id=(num_of_task+1)/numprocs;
    int beg_id=size_per_id*myid;
    int end_id=(myid==numprocs-1)?(num_of_task-1):(beg_id+size_per_id-1);
    ull task_ans=0;
    for(int k=beg_id;k<=end_id;k++)
    {
        int row1,row2;
        if(k<=size-3)
        {
            row1=0,row2=k+2;
        }
        else
        {
            row1=(k-size+2)/(size-3)+1;
            row2=(k-size+2)%(size-3);
            if(row2>=row1-1) row2+=3;
            if(row1==size) row1=size-1,row2=size-3;
        }
        // printf("row1:%d row2:%d id:%d\n",row1,row2,k);
        ull col_arg=shift(row1,size)|shift(row2,size);
        ull d1_arg=shift(row1+1,size)|shift(row2+2,size);
        ull d2_arg=shift(row1-1,size)|shift(row2-2,size);
        task_ans+=queen(size,2,col_arg,d1_arg,d2_arg);
    }
    ull final_ans=task_ans;
    ull token=0;
    if(myid==0)
    {
        for(int i=1;i<numprocs;i++)
        {
            MPI_Recv(&token,1,MPI_UNSIGNED_LONG_LONG,
            i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            final_ans+=token;
        }
    }
    else
    {
        //send subtask ans to task 0
        MPI_Send(&task_ans,1,MPI_UNSIGNED_LONG_LONG,
        0,0,MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return final_ans;
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
ull serial_queen(int n)
{
    return queen(n,0,0,0,0);
}
ull serial_half_queen(int n)
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
    double serial_time[30];
    for(int i=0;i<30;i++) serial_time[i]=-1.0;
    int size;
    for(size=11;size<=15;size++)
    {
        set_time(tm_start);
        ull ans1=serial_queen(size);
        set_time(tm_end);
        assert(ans1==correct_ans[size]);
        serial_time[size]=get_time();
    }

    int myid,numprocs;
    MPI_Init(0,0);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    
    for(size=11;size<=15;size++)
    {
        double time1,time2;
        if(myid==0) set_time(tm_start);
        ull ans2=coarse_grained_mpi_queen(size,myid,numprocs);
        if(myid==0) 
        {
            set_time(tm_end);
            assert(ans2==correct_ans[size]);
            time1=get_time();
        }

        if(myid==0) set_time(tm_start);
        ull ans3=fine_grained_mpi_queen(size,myid,numprocs);
        if(myid==0) 
        {
            set_time(tm_end);
            assert(ans3==correct_ans[size]);
            time2=get_time();
        }

        if(myid==0)
        {
            printf("++++++size %d++++++\n",size);
            printf("serial        :%10.3fms\n",serial_time[size]);
            printf("coarse-grained:%10.3fms\n",time1);
            printf("fine-grained  :%10.3fms\n",time2);
        }
    }
    MPI_Finalize();
    
    return 0;
}