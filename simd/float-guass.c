#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>

#include <immintrin.h> 
// #include<xmmintrin.h>
// #include<pmmintrin.h>

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

float **alloc_mat(int size)
{
    float **m = (float **)malloc(sizeof(float *) * size);
    for (int i = 0; i < size; i++)
    {
        m[i] = (float *)malloc(sizeof(float) * size);
        for (int j = 0; j < size; j++)
        {
            if (i < j)
                m[i][j] = 0.0;
            else if (i == j)
                m[i][j] = 1.0;
            else
                m[i][j] = 1.0 * rand();
        }
    }
    return m;
}

void shuffle_mat(int size, float **mat)
{
    for (int k = 0; k < size; k++)
    {
        for (int i = k + 1; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                mat[i][j] += mat[k][j];
            }
        }
    }
}

float **copy_mat(int size, float **mat)
{
    float **m = (float **)malloc(sizeof(float *) * size);
    for (int i = 0; i < size; i++)
    {
        m[i] = (float *)malloc(sizeof(float) * size);
        for (int j = 0; j < size; j++)
        {
            m[i][j] = mat[i][j];
        }
    }
    return m;
}

float **common_guass(int size, float **mat)
{
    float **m = copy_mat(size, mat);

    set_time(tm_start);
    for (int k = 0; k < size; k++)
    {
        assert(m[k][k] != 0);
        for (int j = k + 1; j < size; j++)
        {
            m[k][j] = m[k][j] / m[k][k]; // if m[k][k]==0?
        }
        m[k][k] = 1.0;
        for (int i = k + 1; i < size; i++)
        {
            for (int j = k + 1; j < size; j++)
            {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            }
            m[i][k] = 0;
        }
    }

    set_time(tm_end);
    printf("common guass:\tsize:%d\t\ttime:%fms\n", size, get_time());
    return m;
}

float **simd_alianed128_guass(int size, float **mat)
{
    float **m = copy_mat(size, mat);

    set_time(tm_start);

    for (int k = 0; k < size; k++)
    {
        assert(m[k][k] != 0);
        for (int j = k + 1; j < size; j++)
        {
            m[k][j] = m[k][j] / m[k][k];
        }
        m[k][k] = 1.0;

        assert(size%4==0);//to simplify the problem
        for (int i = k + 1; i < size; i++)
        {
            int pad_size=(size-k-1)/4;
            int beg_ind=size-pad_size*4;
            //not aligned part
            for(int j=k+1;j<beg_ind;j++)
            {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            }
            //aligned part
            __m128 m_ik=_mm_set_ps1(m[i][k]);//const part
            __m128* m_k=malloc(sizeof(__m128)*pad_size);//row k
            __m128* m_i=malloc(sizeof(__m128)*pad_size);//row i
            for (int j = 0; j < pad_size; j++)
            {
                m_k[j]=_mm_load_ps(m[k]+beg_ind+4*j);
                m_i[j]=_mm_load_ps(m[i]+beg_ind+4*j);

                // m[i][j] = m[i][j] - m[i][k] * m[k][j];
                m_k[j]=_mm_mul_ps(m_ik,m_k[j]);
                m_i[j]=_mm_sub_ps(m_i[j],m_k[j]);

                _mm_store_ps(m[i]+beg_ind+4*j,m_i[j]);
            }
            m[i][k] = 0;
        }
    }
    set_time(tm_end);
    printf("aligned128 guass:\tsize:%d\t\ttime:%fms\n", size, get_time());
    return m;
}

float **simd_alianed256_guass(int size, float **mat)
{
    float **m = copy_mat(size, mat);

    set_time(tm_start);

    for (int k = 0; k < size; k++)
    {
        assert(m[k][k] != 0);
        for (int j = k + 1; j < size; j++)
        {
            m[k][j] = m[k][j] / m[k][k];
        }
        m[k][k] = 1.0;

        assert(size%8==0);//to simplify the problem
        for (int i = k + 1; i < size; i++)
        {
            int pad_size=(size-k-1)/8;
            int beg_ind=size-pad_size*8;
            //not aligned part
            for(int j=k+1;j<beg_ind;j++)
            {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            }
            //aligned part
            __m256 m_ik=_mm256_set1_ps(m[i][k]);//const part
            __m256* m_k=malloc(sizeof(__m256)*pad_size);//row k
            __m256* m_i=malloc(sizeof(__m256)*pad_size);//row i
            for (int j = 0; j < pad_size; j++)
            {
                m_k[j]=_mm256_load_ps(m[k]+beg_ind+4*j);
                m_i[j]=_mm256_load_ps(m[i]+beg_ind+4*j);

                // m[i][j] = m[i][j] - m[i][k] * m[k][j];
                m_k[j]=_mm256_mul_ps(m_ik,m_k[j]);
                m_i[j]=_mm256_sub_ps(m_i[j],m_k[j]);

                _mm256_store_ps(m[i]+beg_ind+4*j,m_i[j]);
            }
            m[i][k] = 0;
        }
    }
    set_time(tm_end);
    printf("aligned256 guass: size:%d\t\t\ttime:%fms\n", size, get_time());
    return m;
}

int comp_mat(int size, float **m1, float **m2)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (m1[i][j] != m2[i][j])
            {
                printf("mat diff at (%d,%d)\n", i, j);
                return 0;
            }
        }
    }
    return 1;
}

int comp_mat_with_iden(int size, float **mat)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            float s = (i == j ? 1.0 : 0.0);
            if (fabs(mat[i][j] - s) > 1e-10)
            {
                printf("mat not identifier at (%d,%d)\n", i, j);
                return 0;
            }
        }
    }
    return 1;
}

void test_size(int n)
{
    float **mat = alloc_mat(n);
    shuffle_mat(n, mat);


    float **mat_g = common_guass(n, mat);
    float **mat2_g = simd_alianed128_guass(n, mat);
    float **mat3_g = simd_alianed256_guass(n, mat);
    assert(comp_mat_with_iden(n, mat_g));
    assert(comp_mat_with_iden(n, mat2_g));
    assert(comp_mat_with_iden(n, mat3_g));
}

int main()
{
    srand(time(0));
    int n = 120;

    while (n < 1001)
    {
        test_size(n);

        if (n < 1000)
            n += 80;
        else
            n += 1000;
    }
    return 0;
}