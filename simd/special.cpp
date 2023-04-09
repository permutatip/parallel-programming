#include<fstream>
#include<sstream>
#include<iostream>
#include<cstring>
#include<cassert>
#include<iomanip>
#include<algorithm>

#include<immintrin.h>

//----linux timer----
#ifdef __linux
#include <sys/time.h>
struct timeval tm_start, tm_end;
#define set_time(x) gettimeofday(&x, 0)
double get_time()
{
    suseconds_t t1 = tm_end.tv_usec - tm_start.tv_usec;
    time_t t2 = tm_end.tv_sec - tm_start.tv_sec;
    return t1 * 0.001 + t2 * 1000.0;
}
#endif
//----windows timer----
#ifdef _WIN32
#include <windows.h>
LARGE_INTEGER tm_start, tm_end, tm_freq;
#define set_time(x) QueryPerformanceCounter(&x)
double get_time()
{
    QueryPerformanceFrequency(&tm_freq);
    long long t1 = tm_end.QuadPart - tm_start.QuadPart;
    long long t2 = tm_freq.QuadPart;
    double t = (double)t1 / t2 * 1000.0;
    return t;
}
#endif
// no timer
#ifndef set_time
#error "Not linux or windows, can not excute."
#endif

#define COL 130
#define ELI_NUM 22
#define TOBE_ELI_NUM 8
#define ARR_LEN (COL/31+1)

struct mat_row
{
    /* data */
    int dat[ARR_LEN];
    int first;
};

mat_row tobe_eli[TOBE_ELI_NUM],elitor[COL+1],ans[TOBE_ELI_NUM];

char line[10000];
const char path1[]="./Groebner/1_130_22_8/1.txt";//消元子
const char path2[]="./Groebner/1_130_22_8/2.txt";//被消元行
const char path3[]="./Groebner/1_130_22_8/3.txt";//消元结果

void show_mat(mat_row mat[],int len)
{
    for(int i=0;i<len;i++)
    {
        if(mat[i].first==0)
        {
            std::cout<<"\n";
            continue;
        }
        for(int j=0;j<ARR_LEN;j++)
        {
            std::cout<<mat[i].dat[j]<<" ";
        }
        std::cout<<"\n";
    }
}

void read_data()
{
    std::fstream fs(path1,fs.in);
    std::string file1;
    std::memset(elitor,0,sizeof(elitor));
    while(std::getline(fs,file1))//消元子
    {
        int temp_i;
        int row_num=0;
        std::istringstream iss(file1);
        while(iss>>temp_i)
        {
            if(row_num==0) 
            {
                row_num=temp_i;
                elitor[row_num].first=row_num;
            }
            int a=temp_i/31,b=temp_i%31;
            elitor[row_num].dat[a]|=(1<<b);
        }
    }
    fs.close();

    fs.open(path2,fs.in);
    std::memset(tobe_eli,0,sizeof(tobe_eli));
    for(int i=0;i<TOBE_ELI_NUM;i++)//被消元行
    {
        int temp_i;
        std::getline(fs,file1);
        std::istringstream iss(file1);
        while(iss>>temp_i)
        {
            int a=temp_i/31,b=temp_i%31;
            if(tobe_eli[i].first==0) tobe_eli[i].first=temp_i;
            tobe_eli[i].dat[a]|=(1<<b);
        }
    }
    assert(fs.eofbit);
    fs.close();

    fs.open(path3,fs.in);
    std::memset(ans,0,sizeof(ans));
    for(int i=0;i<TOBE_ELI_NUM;i++)//结果
    {
        int temp_i;
        std::getline(fs,file1);
        std::istringstream iss(file1);
        while(iss>>temp_i)
        {
            int a=temp_i/31,b=temp_i%31;
            if(ans[i].first==0) ans[i].first=temp_i;
            ans[i].dat[a]|=(1<<b);
        }
    }
    assert(fs.eofbit);
    fs.close();
}

void row_xor(int eli_row,int tobe_eli_row)
{
    for(int i=0;i<ARR_LEN;i++)
    {
        tobe_eli[tobe_eli_row].dat[i]^=elitor[eli_row].dat[i];
    }
    for(int i=ARR_LEN-1;i>=0;i--)
    {
        // for(int j=30;j>=0;j--)
        // {
        //     if(tobe_eli[tobe_eli_row].dat[i]&(1<<j))
        //     {
        //         tobe_eli[tobe_eli_row].first=i*31+j;
        //         return;
        //     }
        // }
        //https://www.zhihu.com/question/35361094
        if(tobe_eli[tobe_eli_row].dat[i])
        {
            int t=std::__lg(tobe_eli[tobe_eli_row].dat[i]);
            tobe_eli[tobe_eli_row].first=i*31+t;
            return;
        }
    }
    tobe_eli[tobe_eli_row].first=0;
}

void sse128_xor(int eli_row,int tobe_eli_row)
{
    int ii=0;
    for(;ii+4<ARR_LEN;ii+=4)
    {
        __m128i tobe_pack=_mm_load_si128((__m128i*)tobe_eli[tobe_eli_row].dat+ii);
        __m128i elitor_pack=_mm_load_si128((__m128i*)elitor[eli_row].dat+ii);
        tobe_pack=_mm_xor_si128(tobe_pack,elitor_pack);
        _mm_storeu_si128((__m128i*)tobe_eli[tobe_eli_row].dat+ii,tobe_pack);
    }
    for(;ii<ARR_LEN;ii++)
    {
        tobe_eli[tobe_eli_row].dat[ii]^=elitor[eli_row].dat[ii];
    }
    // for(int i=0;i<ARR_LEN;i++)
    // {
    //     tobe_eli[tobe_eli_row].dat[i]^=elitor[eli_row].dat[i];
    // }
    for(int i=ARR_LEN-1;i>=0;i--)
    {
        if(tobe_eli[tobe_eli_row].dat[i])
        {
            int t=std::__lg(tobe_eli[tobe_eli_row].dat[i]);
            tobe_eli[tobe_eli_row].first=i*31+t;
            return;
        }
    }
    tobe_eli[tobe_eli_row].first=0;
}

double common()
{
    set_time(tm_start);
    for(int i=0;i<ELI_NUM;i++)
    {
        while(tobe_eli[i].first)
        {
            int f=tobe_eli[i].first;
            if(elitor[f].first)
                row_xor(f,i);
            else
            {
                elitor[f]=tobe_eli[i];
                break;
            }
        }
    }
    set_time(tm_end);
    return get_time();
}

void sse128()
{
    for(int i=0;i<ELI_NUM;i++)
    {
        while(tobe_eli[i].first)
        {
            int f=tobe_eli[i].first;
            if(elitor[f].first)
                sse128_xor(f,i);
            else
            {
                elitor[f]=tobe_eli[i];
                break;
            }
        }
    }
}

void check_ans()
{
    for(int i=0;i<TOBE_ELI_NUM;i++)
    {
        assert(tobe_eli[i].first==ans[i].first);
        for(int j=0;j<ARR_LEN;j++)
        {
            assert(tobe_eli[i].dat[j]==ans[i].dat[j]);
        }
    }
}

int main()
{
    std::ios::sync_with_stdio(0);
    read_data();
    double time1=common();
    check_ans();
    std::cout<<"time1="<<std::fixed<<std::setprecision(3)<<time1<<"ms"<<std::endl;

    // read_data();
    // sse128();
    // check_ans();
    return 0;
}