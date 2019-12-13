#include <iostream>
#include <mpi.h>

#define DIM 1000

using namespace std;

void Mat_vect_mult(
                   double local_A[],
                   double local_x[],
                   double local_y[],
                   int local_m,
                   int n,
                   int local_n,
                   MPI_Comm comm){
    double* x = (double *)malloc(n * sizeof(double));
    int local_i,j;
    MPI_Allgather(local_x,local_n,MPI_DOUBLE,x,local_n,MPI_DOUBLE,comm);
    for(local_i = 0;local_i < local_m;local_i++){
        local_y[local_i] = 0.0;
        for(j = 0;j< n;j++){
            local_y[local_i] += local_A[local_i*n + j] * x[j];
        }
    }
    free(x);
}

int main(int argc,char *argv[])
{
    int mynode, totalnodes;
    //int sum,startval,endval,accum;
    double A[DIM * DIM],x[DIM],y[DIM];
    int m = DIM;
    int n = DIM;
    for(int i = 0;i < m;i++){
        x[i] = 1.0;
        for(int j = 0;j < n;j++){
            A[i*n + j] = 1.0;
        }
    }
    double time;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
    if(mynode == 0){
        time = MPI_Wtime();
    }
    Mat_vect_mult(A,x,y,m,n,totalnodes,MPI_COMM_WORLD);
    if(mynode == 0){
        time = MPI_Wtime() - time;
        cout << "Tempo de programa efetivo nodo 0: " << time << endl;
    }
    /*sum = 0;
    startval = 1000 * mynode / totalnodes + 1;
    endval = 1000 * (mynode + 1) / totalnodes;
    for(int i = startval;i <= endval; i++){
        sum = sum + i;
    }
    MPI_Reduce(&sum,&accum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(mynode == 0){
        time = MPI_Wtime() - time;
        cout << "Tempo de programa efetivo nodo 0: " << time << endl;
    }
    MPI_Finalize();*/
    return 0;
}
