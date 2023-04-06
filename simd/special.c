#include<stdio.h>
#include<stdlib.h>
#include<time.h>

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


void special()
{
    //
}

int main()
{
    srand(time(0));
    return 0;
}