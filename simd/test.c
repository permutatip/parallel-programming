#include<stdio.h>
#include<stdlib.h>

#include<immintrin.h>

void test_mm256()
{
    __m256 pad[2][10];
    float p[2][80];
    for(int i=0;i<80;i++) p[1][i]=1.3f*i,p[0][i]=0.6f*i;
    printf("%p %p\n",p[0],p[0]+8);
    for(int i=0;i<10;i++) pad[0][i]=_mm256_loadu_ps(p[0]+i*8);
    for(int i=0;i<10;i++) pad[1][i]=_mm256_loadu_ps(p[1]+i*8);
    for(int i=0;i<10;i++)
    {
        pad[0][i]=_mm256_add_ps(pad[0][i],pad[0][i]);
        pad[1][i]=_mm256_add_ps(pad[1][i],pad[1][i]);
    }
    for(int i=0;i<10;i++) _mm256_storeu_ps(p[0]+i*8,pad[0][i]);
    for(int i=0;i<10;i++) _mm256_storeu_ps(p[1]+i*8,pad[1][i]);
}

int main()
{
    double** m;

    return 0;
}