#include<mpi.h>
#include<stdio.h>

int main()
{
    int myid,numprocs,namelen;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(0,0);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Get_processor_name(name,&namelen);
    printf("hello from processor%d of %d from %s\n",myid,numprocs,name);
    MPI_Finalize();
    return 0;
}