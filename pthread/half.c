
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
typedef unsigned long long ull;
ull queen(int n,int row,ull col,ull d1,ull d2)
{
    if(row==n) return 1;
    ull count=0;
    ull ava=((1<<n)-1)&(~(col|d1|d2));
    while(ava)
    {
        ull p=ava&(-ava);
        ava=ava-p;
        count+=queen(n,row+1,col|p,((1<<n)-1)&((d1|p)<<1),(d2|p)>>1);
    }
    return count;
}
ull qhalf(int n)
{
    ull count=0;
    /*
    5 0 1 2 3 4:0~1
    4 0 1 2 3 :0~1
    */
    int end_id=n/2-1;
    for(int i=0;i<=end_id;i++)
    {
        count+=2*queen(n,1,1<<i,1<<(i+1),(i>0)?1<<(i-1):0);
    }
    if(n%2) count+=queen(n,1,1<<(n/2),1<<(n/2+1),1<<(n/2-1));
    return count;
}
int main()
{
    for(int k=10;k<=15;k++)
    {
        ull ans1=queen(k,0,0,0,0);
        ull ans2=qhalf(k);
        printf("size:%d common:%llu half:%llu\n",k,ans1,ans2);
    }
    return 0;
}