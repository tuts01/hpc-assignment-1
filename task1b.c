
#define MAXITER 1000
#define N	8000
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>

/* ----------------------------------------------------------------*/

int main(int argc, char** argv) {
    int	   i,j,k,green,blue,loop;
    int	   ierr; //Error number
    int    rank; //Process rank
    int    nump; //Num of processes
    float  *x, *a;
    FILE   *fp;
    double t, t1, t2, t3;

    float complex   z, kappa;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nump);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char* outfile_name = (char*)calloc(6, sizeof(char));
    sprintf(outfile_name, "out%d", rank);
    FILE* outfile = fopen(outfile_name, "a");

    a=(float*) malloc((N*N/10)*sizeof(float));

    if(rank == 0) x=(float *) malloc(N*N*sizeof(float));
    
    t = MPI_Wtime();

    for (loop=0 + (rank * N*N/10); loop< (rank+1)*N*N/10; loop++) {
	    i=loop%N;
	    j=loop/N;

	    z=kappa= (4.0*(i-N/2))/N + (4.0*(j-N/2))/N * I;
	
	    k=1;
	    while ((cabs(z)<=2) && (k++<MAXITER)) 
	        z= z*z + kappa;
	  
	    a[loop - (rank *N*N/10)]= log((float)k) / log((float)MAXITER);
    }
    t1 = MPI_Wtime() - t;
    
    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime() - t1 - t;

    MPI_Gather(a, (N*N/10), MPI_FLOAT, x, (N*N/10), MPI_FLOAT, 0, MPI_COMM_WORLD);
    t3 = MPI_Wtime() - t2 - t1 - t;



/* ----------------------------------------------------------------*/

    fprintf(outfile, "Computation time:\t%f\nWait time:\t\t%f\nCommunication time:\t%f\n\n", t1, t2, t3);
    fclose(outfile);

/* ----------------------------------------------------------------*/

    free(a);
    if(rank==0) free(x);
    MPI_Finalize();
}
