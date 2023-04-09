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

double *alloc_mat(int size)
{
    double *m = (double *)aligned_alloc(32,sizeof(double) * size * size);
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

void free_mat(double *mat)
{
    if(mat)
    {
        free(mat);
        mat=NULL;
    }
}

void shuffle_mat(int size, double *mat)
{
    for (int k = 0; k < size; k++)
    {
        for (int i = k + 1; i < size ; i++)
        {
            double r = (rand()%3-1)/2.0;
            for (int j = 0; j < size; j++)
            {
                mat[i*size+j] += r * mat[k*size+j];
            }
        }
    }
}

void print_mat(int size, double *mat)
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

double *copy_mat(int size, double *mat)
{
    double *m = (double *)aligned_alloc(32, sizeof(double) * size * size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            m[i*size+j] = mat[i*size+j];
        }
    }
    return m;
}

double *common_guass(int size, double *mat)
{
    double *m = copy_mat(size, mat);

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

double *simd_alianed128_guass(int size, double *mat)
{
    double *m = copy_mat(size, mat);
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
            int pad_size = (size - k - 1) / 4;
            int beg_ind = size - pad_size * 4;
            // not aligned part
            for (int j = k + 1; j < beg_ind; j++)
            {
                m[i*size+j] = m[i*size+j] - m[i*size+k] * m[k*size+j];
            }
            // aligned part
            __m128d m_ik = _mm_set_pd1(m[i*size+k]);              // const part
            for (int j = 0; j < pad_size; j++)
            {
                __m128d m_kj = _mm_load_pd(m+ k*size + beg_ind + 4 * j);
                __m128d m_ij = _mm_load_pd(m+ i*size + beg_ind + 4 * j);
                // m[i*size+j] = m[i*size+j] - m[i*size+k] * m[k*size+j];
                m_kj = _mm_mul_pd(m_ik, m_kj);
                m_ij = _mm_sub_pd(m_ij, m_kj);
                _mm_store_pd(m+ i*size + beg_ind + 4 * j, m_ij);
            }
            m[i*size+k] = 0;
        }
    }
    set_time(tm_end);
    printf("sse128:%d %.3fms\n", size, get_time());
    return m;
}

double *simd_alianed256_guass(int size, double *mat)
{
    double *m = copy_mat(size, mat);
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
            int pad_size = (size - k - 1) / 8;
            int beg_ind = size - pad_size * 8;
            // not aligned part
            for (int j = k + 1; j < beg_ind; j++)
            {
                m[i*size+j] = m[i*size+j] - m[i*size+k] * m[k*size+j];
            }
            // aligned part
            __m256d m_ik = _mm256_set1_pd(m[i*size+k]);           // const part
            for (int j = 0; j < pad_size; j++)
            {
                __m256d m_kj = _mm256_load_pd(m+ k*size + beg_ind + 8 * j);
                __m256d m_ij=_mm256_load_pd(m+ i*size+beg_ind+8*j);
                // m[i*size+j] = m[i*size+j] - m[i*size+k] * m[k*size+j];
                m_kj=_mm256_mul_pd(m_ik,m_kj);
                m_ij=_mm256_sub_pd(m_ij,m_kj);
                _mm256_store_pd(m+ i*size+beg_ind+8*j,m_ij);
            }
            m[i*size+k] = 0;
        }
    }
    set_time(tm_end);
    printf("avx256:%d %.3fms\n", size, get_time());
    return m;
}

int comp_mat(int size, double *m1, double *m2)
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
    double *mat = alloc_mat(n);
    double *mat2 = copy_mat(n, mat);
    assert(comp_mat(n, mat, mat2));
    shuffle_mat(n, mat);

    double *mat_g = common_guass(n, mat);
    double *mat2_g = simd_alianed128_guass(n, mat);
    assert(comp_mat(n, mat_g, mat2_g));
    double *mat3_g = simd_alianed256_guass(n, mat);
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
    int n = 40;

    while (n <= 1000)
    {
        test_size(n);
        n+=80;
    }
    return 0;
}