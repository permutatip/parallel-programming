#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<assert.h>
#include<math.h>

#define MAT_TYPE float

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

MAT_TYPE** alloc_mat(int size)
{
    MAT_TYPE** m=(MAT_TYPE**)malloc(sizeof(MAT_TYPE*)*size);
    for(int i=0;i<size;i++)
    {
        m[i]=(MAT_TYPE*)malloc(sizeof(MAT_TYPE)*size);
        for(int j=0;j<size;j++)
        {
            if(i<j) m[i][j]=0.0;
            else if(i==j) m[i][j]=1.0;
            else m[i][j]=1.0*rand();
        }
    }
    return m;
}

void shuffle_mat(int size,MAT_TYPE** mat)
{
    for(int k=0;k<size;k++)
    {
        for(int i=k+1;i<size;i++)
        {
            for(int j=0;j<size;j++)
            {
                mat[i][j]+=mat[k][j];
            }
        }
    }
}

MAT_TYPE** copy_mat(int size,MAT_TYPE** mat)
{
    MAT_TYPE** m=(MAT_TYPE**)malloc(sizeof(MAT_TYPE*)*size);
    for(int i=0;i<size;i++)
    {
        m[i]=(MAT_TYPE*)malloc(sizeof(MAT_TYPE)*size);
        for(int j=0;j<size;j++)
        {
            m[i][j]=mat[i][j];
        }
    }
    return m;
}

MAT_TYPE** common_guass(int size,MAT_TYPE** mat)
{
    MAT_TYPE** m=copy_mat(size,mat);

    set_time(tm_start);
    for(int k=0;k<size;k++)
    {
        assert(m[k][k]!=0);
        for(int j=k+1;j<size;j++)
        {
            m[k][j]=m[k][j]/m[k][k];//if m[k][k]==0?
        }
        m[k][k]=1.0;
        for(int i=k+1;i<size;i++)
        {
            for(int j=k+1;j<size;j++)
            {
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0;
        }
    }

    set_time(tm_end);
    printf("common guass: size:%d time:%fms\n",size,get_time());
    return m;
}

MAT_TYPE** simd_alianed_guass(int size,MAT_TYPE** mat)
{
    MAT_TYPE** m=copy_mat(size,mat);

    set_time(tm_start);
    for(int k=0;k<size;k++)
    {
        assert(m[k][k]!=0);
        for(int j=k+1;j<size;j++)
        {
            m[k][j]=m[k][j]/m[k][k];//if m[k][k]==0?
        }
        m[k][k]=1.0;
        for(int i=k+1;i<size;i++)
        {
            for(int j=k+1;j<size;j++)
            {
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0;
        }
    }
    set_time(tm_end);
    printf("simd guass: size:%d time:%fms\n",size,get_time());
    return m;
}

int comp_mat(int size,MAT_TYPE** m1,MAT_TYPE** m2)
{
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            if(m1[i][j]!=m2[i][j]) 
            {
                printf("mat diff at (%d,%d)\n",i,j);
                return 0;
            }
        }
    }
    return 1;
}

int comp_mat_with_iden(int size,MAT_TYPE** mat)
{
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<size;j++)
        {
            MAT_TYPE s=(i==j?1.0:0.0);
            if(fabs(mat[i][j]-s)>1e-10)
            {
                printf("mat not identifier at (%d,%d)\n",i,j);
                return 0;
            }
        }
    }
    return 1;
}

int main()
{
    srand(time(0));
    int n=100;

    while(n<10001)
    {
        MAT_TYPE** mat=alloc_mat(n);
        MAT_TYPE** m2=copy_mat(n,mat);
        assert(comp_mat(n,mat,m2));

        shuffle_mat(n,mat);
        MAT_TYPE** mat_g=common_guass(n,mat);
        assert(comp_mat_with_iden(n,mat_g));
        if(n<1000) n+=100;
        else n+=1000;
    }
    return 0;
}