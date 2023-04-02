#include<stdio.h>
#include<stdlib.h>

#include<immintrin.h>

int main()
{
    __m128 pad[10];
    float p[2][41];
    for(int i=0;i<41;i++) p[1][i]=1.0f*i,p[0][i]=0.5f*i;
    for(int i=0;i<10;i++) pad[i]=_mm_load_ps(p[1]+i*4+1);
    for(int i=0;i<10;i++)
    {
        pad[i]=_mm_add_ss(pad[i],pad[i]);
    }
    for(int i=0;i<10;i++) _mm_store_ss(p[1]+i*4+1,pad[i]);
    for(int i=0;i<41;i++) printf("%.4f ",p[1][i]);
    printf("\n");

    return 0;
}