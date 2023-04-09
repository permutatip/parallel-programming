#include<stdio.h>
#include<stdlib.h>

#include<immintrin.h>

// void test_mm256()
// {
//     __m256 pad[2][10];
//     float p[2][80];
//     for(int i=0;i<80;i++) p[1][i]=1.3f*i,p[0][i]=0.6f*i;
//     printf("%p %p\n",p[0],p[0]+8);
//     for(int i=0;i<10;i++) pad[0][i]=_mm256_loadu_ps(p[0]+i*8);
//     for(int i=0;i<10;i++) pad[1][i]=_mm256_loadu_ps(p[1]+i*8);
//     for(int i=0;i<10;i++)
//     {
//         pad[0][i]=_mm256_add_ps(pad[0][i],pad[0][i]);
//         pad[1][i]=_mm256_add_ps(pad[1][i],pad[1][i]);
//     }
//     for(int i=0;i<10;i++) _mm256_storeu_ps(p[0]+i*8,pad[0][i]);
//     for(int i=0;i<10;i++) _mm256_storeu_ps(p[1]+i*8,pad[1][i]);
// }

/*
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
            __m128 m_ik = _mm_set_ps1(m[i*size+k]);              // const part
            __m128 *m_k = malloc(sizeof(__m128) * pad_size); // row k
            __m128 *m_i = malloc(sizeof(__m128) * pad_size); // row i
            for (int j = 0; j < pad_size; j++)
            {
                m_k[j] = _mm_load_ps(m+ k*size + beg_ind + 4 * j);
                m_i[j] = _mm_load_ps(m+ i*size + beg_ind + 4 * j);

                // m[i*size+j] = m[i*size+j] - m[i*size+k] * m[k*size+j];
                m_k[j] = _mm_mul_ps(m_ik, m_k[j]);
                m_i[j] = _mm_sub_ps(m_i[j], m_k[j]);

                _mm_store_ps(m+ i*size + beg_ind + 4 * j, m_i[j]);
            }
            free(m_k);
            free(m_i);
            m[i*size+k] = 0;
        }
*/

int main()
{
    // float* m=(float*)aligned_alloc()
    return 0;
}