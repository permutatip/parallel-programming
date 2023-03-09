#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#define MAX_NUM 1000
#define MAX_SIZE 100000

int n;
int a[MAX_SIZE];
int b[MAX_SIZE];//as a helper
FILE* fp;

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

void two_ways_add()
{
    int result=0;
    int i=0;
    set_time(tm_start);
    for(;i+1<n;i+=2)
    {
        result+=a[i];
        result+=a[i+1];
    }
    for(;i<n;i++) result+=a[i];
    set_time(tm_end);
    // printf("two ways time:%fms\n",get_time());
    fprintf(fp,"%f\t",get_time());
    printf("two ways res: %d\n",result);
}

void two_ways_new_add()
{
    int result=0,res1=0,res2=0;
    int i=0;
    set_time(tm_start);
    for(;i+1<n;i+=2)
    {
        res1+=a[i];
        res2+=a[i+1];
    }
    for(;i<n;i++) res1+=a[i];
    result=res1+res2;
    set_time(tm_end);
    // printf("two ways time:%fms\n",get_time());
    fprintf(fp,"%f\t",get_time());
    printf("two ways res: %d\n",result);
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
        a[i]=rand()%MAX_NUM;
    }
    
    fp=fopen("add_stat2.txt","a");
    if(fp==NULL)
    {
        printf("file open error\n");
        return 0;
    }

    // fseek(fp,0,0);
    long pos=ftell(fp);
    if(pos==-1)
    {
        printf("pos error\n");
        return 0;
    }
    else if(pos==0)
    {
        fprintf(fp,"\t%s\t%s\n","two_ways","two_ways_new");
    }

    fprintf(fp,"%s\t",argv[0]);
    two_ways_add();
    two_ways_new_add();

    fprintf(fp,"size:%d\n",n);
    fclose(fp);
    return 0;
}