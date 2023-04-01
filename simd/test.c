#include<stdio.h>
#include<stdlib.h>

#include<immintrin.h>

int main()
{
    __m128 pad[10];
    float p[40];
    for(int i=0;i<40;i++) p[i]=1.0f*i;
    for(int i=0;i<9;i++) pad[i]=_mm_load_ss(p+i*4+1);
    for(int i=0;i<9;i++)
    {
        pad[i]=_mm_add_ss(pad[i],pad[i]);
    }
    for(int i=0;i<9;i++) _mm_store_ss(p+i*4+1,pad[i]);
    for(int i=0;i<40;i++) printf("%.4f ",p[i]);
    printf("\n");

    return 0;
}