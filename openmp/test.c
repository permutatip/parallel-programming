#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#ifndef _OPENMP
#error "this is a openmp file"
#endif
#include<omp.h>
typedef unsigned long long ull;

// ull queen(int n,int row,ull col,ull d1,ull d2)
// {
//     if(row==n) return 1;
//     ull count=0;
//     ull ava=((1<<n)-1)&(~(col|d1|d2));
//     while(ava)
//     {
//         ull p=ava&(-ava);
//         ava=ava-p;
//         count+=queen(n,row+1,col|p,((1<<n)-1)&((d1|p)<<1),(d2|p)>>1);
//     }
//     return count;
// }
int main(int argc,char* argv[])
{
    assert(argc==2);
    int thrd_cnt=atoi(argv[1]);
    assert(thrd_cnt>0&&thrd_cnt<10);
    int sum=0;

    #pragma omp parallel for\
        reduction(+:sum)
    for(int i=0;i<10;i++)
    {
        sum+=i;
    }
    printf("sum=%d\n",sum);

    return 0;
}