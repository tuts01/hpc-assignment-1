
#define MAXITER 1000
#define N	8000
#define MASTER	0
#include <stdbool.h>
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
    int    chunksize; //Chunksize
    float  *x, *a;
    FILE   *fp;

    float complex   z, kappa;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nump);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(nump < 2)
    {
        puts("Error:\tNot Enough Processes\n");
        MPI_Finalize();
        exit(1);
    }

    MPI_Status status;

    chunksize = atoi(argv[1]);
    //a=(float*) malloc((N*N/10)*sizeof(float));
    if(rank == MASTER)
    {
        a=(float *) malloc(chunksize*sizeof(float));
        x=(float *) malloc(N*N*sizeof(float));

        int startNum = 0;
        int endNum = chunksize - 1;

        int index = 0;

        for(int i = 1; i <= nump; i++)
        {
            MPI_Send(&startNum, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&endNum, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

            startNum += chunksize;
            endNum += chunksize;
        }

        while(endNum < N*N)
        {
            MPI_Recv(a, chunksize, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            MPI_Send(&startNum, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
            MPI_Send(&endNum, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);

            if( (N*N - endNum) > chunksize)
            {
                for(int i = 0; i < chunksize; i++)
                {
                    x[index] = a[i];
                    index++;
                }

                startNum += chunksize;
                endNum += chunksize;
            }
            else
            {
                startNum += chunksize;
                endNum = N*N;

                MPI_Recv(a, chunksize, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

                MPI_Send(&startNum, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
                MPI_Send(&endNum, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);

                for(int i = 0; i < (endNum - startNum); i++)
                {
                    x[index] = a[i];
                    index++;
                }
            }
        }

        for(int i = 1; i <= nump; i++)
        {

        }

        free(a);
        
    }
    else //Worker
    {
        int startNum;
        bool done = false;

        while(!done)
        {
            a=(float *) malloc(chunksize*sizeof(float));
            MPI_Recv(&startNum, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, NULL);

            if(startNum < 0) done = true;
            else
            {
                MPI_Recv(&endNum, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, NULL);

                for (loop = startNum; loop < endNum; loop++) {
        	        i=loop%N;
        	        j=loop/N;

        	        z=kappa= (4.0*(i-N/2))/N + (4.0*(j-N/2))/N * I;

        	        k=1;
        	        while ((cabs(z)<=2) && (k++<MAXITER))
        	            z= z*z + kappa;

        	        a[loop - startNum]= log((float)k) / log((float)MAXITER);
                }

                MPI_Send(a, chunksize, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD, NULL);
            }

            free(a);
        }
    }


/*    for (loop=0 + (rank * N*N/10); loop< (rank+1)*N*N/10; loop++) {
	    i=loop%N;
	    j=loop/N;

	    z=kappa= (4.0*(i-N/2))/N + (4.0*(j-N/2))/N * I;

	    k=1;
	    while ((cabs(z)<=2) && (k++<MAXITER))
	        z= z*z + kappa;

	    a[loop - (rank *N*N/10)]= log((float)k) / log((float)MAXITER);
    }
*/
//    MPI_Gather(a, (N*N/10), MPI_FLOAT, x, (N*N/10), MPI_FLOAT, 0, MPI_COMM_WORLD);

/* ----------------------------------------------------------------*/


    if(rank == MASTER){

        printf("Writing mandelbrot.ppm\n");
        fp = fopen ("mandelbrot.ppm", "w");
        fprintf (fp, "P3\n%4d %4d\n255\n", N, N);

        for (loop=0; loop<N*N; loop++)
	    if (x[loop]<0.5) {
	        green= (int) (2*x[loop]*255);
                fprintf (fp, "%3d\n%3d\n%3d\n", 255-green,green,0);
	    } else {
	        blue= (int) (2*x[loop]*255-255);
                fprintf (fp, "%3d\n%3d\n%3d\n", 0,255-blue,blue);
	    }

        fclose(fp);

    }

/* ----------------------------------------------------------------*/

    //free(a);
    if(rank==MASTER) free(x);
    MPI_Finalize();
}
