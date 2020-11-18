# include <iostream>
# include <cstdlib>
# include <mpi.h>


int main( int nargs, char* argv[] )
{
	
	MPI_Init( &nargs, &argv );
	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
	int nbp;
	MPI_Comm_size(globComm, &nbp);
	int rank;
	int value;
	MPI_Comm_rank(globComm, &rank);
	MPI_Status status;
	

		if (rank==0) {
			fprintf(stderr, "\nEntrez une nouvelle valeur=");
            scanf("%d",&value);
            fprintf(stderr, "Le processus %d recoit la valeur (%d)\n",rank,value);

            if (nbp>1) {
                MPI_Send(&value, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
            }
       }
       else {
             MPI_Recv(&value, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
             value+=1;
             fprintf(stderr, "Le processus %d recoit la valeur (%d)\n",rank, value);
             if (rank<nbp-1) {
                 MPI_Send(&value, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
             }
             if (rank==nbp-1) {
                 MPI_Send(&value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                 value+=1;
                 fprintf(stderr, "Le processus 0 recoit la valeur (%d)\n", value);
             }
         }
       MPI_Barrier(MPI_COMM_WORLD);
	
	
	MPI_Finalize();
	return EXIT_SUCCESS;
}
