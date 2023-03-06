#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#define MAX_NUM 1000
#define MAX_SIZE 100000

int n;
int a[MAX_SIZE];
int b[MAX_SIZE];//as a helper
int offset=0;
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

void traditional_add()
{
    int result=0;
    set_time(tm_start);
    for(int i=0;i<n;i++)
    {
        result+=a[i];
    }
    set_time(tm_end);
    // printf("traditional time:%fms\n",get_time());
    fprintf(fp,"%f\t",get_time());
    printf("traditional res: %d\n",result);
}

void two_ways_add()
{
    int result=0;
    set_time(tm_start);
    for(int i=0;i<n;i+=2)
    {
        result+=a[i];
        result+=a[i+1];
    }
    set_time(tm_end);
    // printf("two ways time:%fms\n",get_time());
    fprintf(fp,"%f\t",get_time());
    printf("two ways res: %d\n",result);
}

int recursion(int size)
{
    if(size==1)
    {
        return b[0];
    }
    else
    {
        for(int i=0;i<size/2;i++)
        {
            b[i]+=b[size-i-1];
        }
        size=(size+1)/2;
        recursion(size);
    }
}

void recursion_add()
{
    int result=0;
    memset(b,0,sizeof(b));
    for(int i=0;i<n;i++) b[i]=a[i];

    set_time(tm_start);
    result=recursion(n);
    set_time(tm_end);

    // printf("recursion time:%fms\n",get_time());
    fprintf(fp,"%f\t",get_time());
    printf("recursion res: %d\n",result);
}

void loop2_add()
{
    int result=0;
    memset(b,0,sizeof(b));
    for(int i=0;i<n;i++) b[i]=a[i];


    set_time(tm_start);
    for(int m=n;m>1;m=(m+1)/2)
    {
        for(int i=0;i<m/2;i++)
        {
            b[i]+=b[m-1-i];
        }
    }
    result=b[0];
    set_time(tm_end);
    // printf("loop2 time:%fms\n",get_time());
    fprintf(fp,"%f\t",get_time());
    printf("loop2 res: %d\n",result);
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
    
    fp=fopen("add_stat_bin.csv","ab");
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
        fprintf(fp,"\t%s\t%s\t%s\t%s\n","traditional","two_ways","recursion","loop2");
    }

    fprintf(fp,"%s\t",argv[0]);
    traditional_add();
    two_ways_add();
    recursion_add();
    loop2_add();

    fprintf(fp,"size:%d\n",n);
    fclose(fp);
    return 0;
}