#include<iostream>
#include<cassert>

using namespace std;
typedef unsigned long long ull;


int board[100];
bool line[100],dia1[200],dia2[200];
int size;
//----linux timer----
#ifdef __linux
#include<sys/time.h>
struct timeval tm_start,tm_end;
#define set_time(x) gettimeofday(&x,0)
double get_time()
{
    suseconds_t t1=tm_end.tv_usec-tm_start.tv_usec;
    time_t t2=tm_end.tv_sec-tm_start.tv_sec;
    return t1*0.000001 + t2;
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
    double t=(double)t1/t2;
    return t;
}
#endif

//use global:cnt size board line dia1 dia2
ull queen_all(int pos,int siz)
{
    if(pos==siz) return 1;
    ull cnt=0;
    for(int i=0;i<siz;i++)
    {
        if(!line[i]&&!dia1[i-pos+siz]&&!dia2[i+pos])
        {
            board[pos]=i;
            line[i]=dia1[i-pos+siz]=dia2[i+pos]=1;
            cnt+=queen_all(pos+1,siz);
            line[i]=dia1[i-pos+siz]=dia2[i+pos]=0;
            board[pos]=0;
        }
    }
    return cnt;
}

ull queen_bit(int n,int row,ull col,ull d1,ull d2)
{
    if(row==n) return 1;
    ull count=0;
    ull ava=((1<<n)-1)&(~(col|d1|d2));
    while(ava)
    {
        ull p=ava&(-ava);
        ava=ava-p;
        count+=queen_bit(n,row+1,col|p,(d1|p)<<1,(d2|p)>>1);
    }
    return count;
}

int main()
{
    ios::sync_with_stdio(0);
    for(int i=0;i<100;i++) board[i]=line[i]=0;
    for(int i=0;i<200;i++) dia1[i]=dia2[i]=0;

    for(size=11;size<=16;size++)
    {
        set_time(tm_start);
        ull cnt=queen_all(0,size);
        set_time(tm_end);
        double time=get_time();

        set_time(tm_start);
        ull cnt2=queen_bit(size,0,0,0,0);
        set_time(tm_end);
        double time2=get_time();

        assert(cnt2==cnt);
        cout<<size<<" : "<<cnt<<" answers\t\t";
        cout<<"normal: "<<time<<"\t\t";
        cout<<"bit: "<<time2<<"\n";
    }
    return 0;
}