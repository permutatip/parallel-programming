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

float **alloc_mat(int size)
{
    float **m = malloc(sizeof(float *) * size);
    for (int i = 0; i < size; i++)
    {
        m[i] = (float *)malloc(sizeof(float) * size);
        for (int j = 0; j < size; j++)
        {
            if (i > j)
                m[i][j] = 0.0;
            else if (i == j)
                m[i][j] = 1.0;
            else
                m[i][j] = 0.5 * (rand() % 10 - 5);
        }
    }
    return m;
}

void shuffle_mat(int size, float **mat)
{
    for (int k = 0; k < size; k++)
    {
        for (int i = k + 1; i < size ; i++)
        {
            float r = rand()%3 - 1;
            for (int j = 0; j < size; j++)
            {
                mat[i][j] += r * mat[k][j];
            }
        }
    }
}

void print_mat(int size, float **mat)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%6.3f ", mat[i][j]);
        }
        printf("\n");
    }
}

float **copy_mat(int size, float **mat)
{
    float **m = malloc(sizeof(float *) * size);
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
        if (fabs(m[k][k]) <= EPS)
        {
            printf("m[%d][%d] of size%d is zero\n", k, k, size);
            exit(1);
        }
        for (int j = k + 1; j < size; j++)
        {
            m[k][j] = m[k][j] / m[k][k];
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
        // print_mat(size,m);
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
        if (fabs(m[k][k]) <= EPS)
        {
            printf("m[%d][%d] of size%d is zero\n", k, k, size);
            exit(1);
        }
        for (int j = k + 1; j < size; j++)
        {
            m[k][j] = m[k][j] / m[k][k];
        }
        m[k][k] = 1.0;

        assert(size % 4 == 0); // to simplify the problem
        for (int i = k + 1; i < size; i++)
        {
            int pad_size = (size - k - 1) / 4;
            int beg_ind = size - pad_size * 4;
            // not aligned part
            for (int j = k + 1; j < beg_ind; j++)
            {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            }
            // aligned part
            __m128 m_ik = _mm_set_ps1(m[i][k]);              // const part
            __m128 *m_k = malloc(sizeof(__m128) * pad_size); // row k
            __m128 *m_i = malloc(sizeof(__m128) * pad_size); // row i
            for (int j = 0; j < pad_size; j++)
            {
                m_k[j] = _mm_load_ps(m[k] + beg_ind + 4 * j);
                m_i[j] = _mm_load_ps(m[i] + beg_ind + 4 * j);

                // m[i][j] = m[i][j] - m[i][k] * m[k][j];
                m_k[j] = _mm_mul_ps(m_ik, m_k[j]);
                m_i[j] = _mm_sub_ps(m_i[j], m_k[j]);

                _mm_store_ps(m[i] + beg_ind + 4 * j, m_i[j]);
            }
            free(m_k);
            free(m_i);
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

        assert(size % 8 == 0); // to simplify the problem
        for (int i = k + 1; i < size; i++)
        {
            int pad_size = (size - k - 1) / 8;
            int beg_ind = size - pad_size * 8;
            // not aligned part
            for (int j = k + 1; j < beg_ind; j++)
            {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            }
            // aligned part
            __m256 m_ik = _mm256_set1_ps(m[i][k]);           // const part
            __m256 *m_k = malloc(sizeof(__m256) * pad_size); // row k
            __m256 *m_i = malloc(sizeof(__m256) * pad_size); // row i
            for (int j = 0; j < pad_size; j++)
            {
                // assert(((unsigned int))%32==0);
                printf("%p\n", m[k] + beg_ind + 8 * j);
                m_k[j] = _mm256_loadu_ps(m[k] + beg_ind + 8 * j);
                // m_i[j]=_mm256_loadu_ps(m[i]+beg_ind+8*j);

                // m[i][j] = m[i][j] - m[i][k] * m[k][j];
                // m_k[j]=_mm256_mul_ps(m_ik,m_k[j]);
                // m_i[j]=_mm256_sub_ps(m_i[j],m_k[j]);

                // _mm256_storeu_ps(m[i]+beg_ind+8*j,m_i[j]);
            }
            m[i][k] = 0;
            free(m_k);
            free(m_i);
        }
    }
    set_time(tm_end);
    printf("aligned256 guass:\tsize:%d\t\ttime:%fms\n", size, get_time());
    return m;
}

int comp_mat(int size, float **m1, float **m2)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (fabs(m1[i][j] - m2[i][j]) > EPS)
            {
                printf("size%d mat diff at (%d,%d)\n", size, i, j);
                printf("m1[i][j]=%f m2[i][j]=%f\n",m1[i][j],m2[i][j]);
                return 0;
            }
        }
    }
    return 1;
}
void print_test_size(int n)
{
    float **mat = alloc_mat(n);
    float **mat2 = copy_mat(n, mat);
    print_mat(n,mat);
    printf("--------------------\n");
    assert(comp_mat(n, mat, mat2));
    shuffle_mat(n, mat);
    print_mat(n,mat);

    float **mat_g = common_guass(n, mat);
    assert(comp_mat(n, mat2, mat_g));
    float **mat2_g = simd_alianed128_guass(n, mat);
    assert(comp_mat(n, mat2, mat2_g));
    // float **mat3_g = simd_alianed256_guass(n, mat);
    // assert(comp_mat_with_iden(n, mat3_g));
    free(mat), free(mat2), free(mat_g);
    free(mat2_g);
    // free(mat3_g);
}
void test_size(int n)
{
    float **mat = alloc_mat(n);
    float **mat2 = copy_mat(n, mat);
    // print_mat(n,mat);
    assert(comp_mat(n, mat, mat2));
    shuffle_mat(n, mat);
    // print_mat(n,mat);

    float **mat_g = common_guass(n, mat);
    assert(comp_mat(n, mat2, mat_g));
    float **mat2_g = simd_alianed128_guass(n, mat);
    assert(comp_mat(n, mat2, mat2_g));
    // float **mat3_g = simd_alianed256_guass(n, mat);
    // assert(comp_mat_with_iden(n, mat3_g));
    free(mat), free(mat2), free(mat_g);
    free(mat2_g);
    // free(mat3_g);
}

int main()
{
    srand(time(0));
    int n = 40;

    test_size(8);
    print_test_size(24);
    // print_test_size(16);

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