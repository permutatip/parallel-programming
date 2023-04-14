#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <math.h>

#include <immintrin.h>

#define EPS 1e-6

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

float *alloc_mat(int size)
{
    float *m = (float *)aligned_alloc(32,sizeof(float) * size * size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i > j)
                m[i*size+j] = 0.0;
            else if (i == j)
                m[i*size+j] = 1.0;
            else
                m[i*size+j] = rand()%10;
        }
    }
    return m;
}

void free_mat(float *mat)
{
    if(mat)
    {
        free(mat);
        mat=NULL;
    }
}

void shuffle_mat(int size, float *mat)
{
    for (int k = 0; k < size; k++)
    {
        for (int i = k + 1; i < size ; i++)
        {
            float r = (rand()%3-1)/2.0;
            for (int j = 0; j < size; j++)
            {
                mat[i*size+j] += r * mat[k*size+j];
            }
        }
    }
}

void print_mat(int size, float *mat)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%6.3f ", mat[i*size+j]);
        }
        printf("\n");
    }
}

float *copy_mat(int size, float *mat)
{
    float *m = (float *)aligned_alloc(32, sizeof(float) * size * size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            m[i*size+j] = mat[i*size+j];
        }
    }
    return m;
}

float *common_guass(int size, float *mat)
{
    float *m = copy_mat(size, mat);

    set_time(tm_start);
    for (int k = 0; k < size; k++)
    {
        for (int j = k + 1; j < size; j++)
        {
            m[k*size+j] = m[k*size+j] / m[k*size+k];
        }
        m[k*size+k] = 1.0;
        for (int i = k + 1; i < size; i++)
        {
            for (int j = k + 1; j < size; j++)
            {
                m[i*size+j] = m[i*size+j] - m[i*size+k] * m[k*size+j];
            }
            m[i*size+k] = 0;
        }
        // print_mat(size,m);
    }
    set_time(tm_end);
    printf("common:%d %.3fms\n", size, get_time());
    return m;
}

float *simd_ualianed128_guass(int size, float *mat)
{
    float *m = copy_mat(size, mat);
    set_time(tm_start);
    for (int k = 0; k < size; k++)
    {
        for (int j = k + 1; j < size; j++)
        {
            m[k*size+j] = m[k*size+j] / m[k*size+k];
        }
        m[k*size+k] = 1.0;

        for (int i = k + 1; i < size; i++)
        {
            __m128 m_ik = _mm_set_ps1(m[i*size+k]);              // const part
            int jj=k+1;
            for(;jj+4<size;jj+=4)
            {
                __m128 m_kj = _mm_loadu_ps(m+ k*size + jj);
                __m128 m_ij = _mm_loadu_ps(m+ i*size + jj);
                // m[i*size+j] = m[i*size+j] - m[i*size+k] * m[k*size+j];
                m_kj = _mm_mul_ps(m_ik, m_kj);
                m_ij = _mm_sub_ps(m_ij, m_kj);
                _mm_storeu_ps(m+ i*size + jj, m_ij);
            }
            for (;jj<size; jj++)
            {
                m[i*size+jj] = m[i*size+jj] - m[i*size+k] * m[k*size+jj];
            }
            m[i*size+k] = 0;
        }
    }
    set_time(tm_end);
    printf("sse128:%d %.3fms\n", size, get_time());
    return m;
}

float *simd_ualianed256_guass(int size, float *mat)
{
    float *m = copy_mat(size, mat);
    set_time(tm_start);
    for (int k = 0; k < size; k++)
    {
        for (int j = k + 1; j < size; j++)
        {
            m[k*size+j] = m[k*size+j] / m[k*size+k];
        }
        m[k*size+k] = 1.0;
        for (int i = k + 1; i < size; i++)
        {
            int jj=k+1;
            __m256 m_ik = _mm256_set1_ps(m[i*size+k]);           // const part
            for(;jj+8<size;jj+=8)
            {
                __m256 m_kj = _mm256_loadu_ps(m+ k*size + jj);
                __m256 m_ij=_mm256_loadu_ps(m+ i*size+jj);
                // m[i*size+j] = m[i*size+j] - m[i*size+k] * m[k*size+j];
                m_kj=_mm256_mul_ps(m_ik,m_kj);
                m_ij=_mm256_sub_ps(m_ij,m_kj);
                _mm256_storeu_ps(m+ i*size+jj,m_ij);
            }
            for(;jj<size;jj++)
            {
                m[i*size+jj] = m[i*size+jj] - m[i*size+k] * m[k*size+jj];
            }
            m[i*size+k] = 0;
        }
    }
    set_time(tm_end);
    printf("avx256:%d %.3fms\n", size, get_time());
    return m;
}

int comp_mat(int size, float *m1, float *m2)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (fabs(m1[i*size+j] - m2[i*size+j]) > EPS)
            {
                printf("size%d mat diff at (%d,%d)\n", size, i, j);
                printf("m1[i*size+j]=%f m2[i*size+j]=%f\n",m1[i*size+j],m2[i*size+j]);
                return 0;
            }
        }
    }
    return 1;
}

void test_size(int n)
{
    float *mat = alloc_mat(n);
    float *mat2 = copy_mat(n, mat);
    assert(comp_mat(n, mat, mat2));
    shuffle_mat(n, mat);

    float *mat_g = common_guass(n, mat);
    float *mat2_g = simd_ualianed128_guass(n, mat);
    assert(comp_mat(n, mat_g, mat2_g));
    float *mat3_g = simd_ualianed256_guass(n, mat);
    assert(comp_mat(n, mat_g, mat3_g));

    free_mat(mat);
    free_mat(mat2);
    free_mat(mat_g);
    free_mat(mat2_g);
    free_mat(mat3_g);
}

int main()
{
    srand(time(0));
    int n = 256;
    while(n<=1024)
    {
        test_size(n);
        n*=2;
    }
    return 0;
}