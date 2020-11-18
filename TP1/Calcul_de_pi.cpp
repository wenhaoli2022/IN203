# include <chrono>
# include <random>
# include <cstdlib>
# include <sstream>
# include <string>
# include <fstream>
# include <iostream>
# include <iomanip>
# include <mpi.h>

// Attention , ne marche qu'en C++ 11 ou supérieur :
double approximate_pi( unsigned long nbSamples ) 
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = beginning.time_since_epoch();
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution <double> distribution ( -1.0 ,1.0);
    unsigned long nbDarts = 0;
    // Throw nbSamples darts in the unit square [-1 :1] x [-1 :1]
    for ( unsigned sample = 0 ; sample < nbSamples ; ++ sample ) {
        double x = distribution(generator);
        double y = distribution(generator);
        // Test if the dart is in the unit disk
        if ( x*x+y*y<=1 ) nbDarts ++;
    }
    // Number of nbDarts throwed in the unit disk
    double ratio = double(nbDarts)/double(nbSamples);
    return 4*ratio;
}

int main( int nargs, char* argv[] )
{
	MPI_Init( &nargs, &argv );
	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
	int nbp;
	MPI_Comm_size(globComm, &nbp);
	int rank;
	MPI_Comm_rank(globComm, &rank);

	// Rajout de code....
	int nb_point = 100000; // Nombre de points total (Ici on choisit 1000.)
	unsigned long nbSamples = (unsigned long)(nb_point/nbp);;
	double pi_part = approximate_pi(nbSamples); // La valeur de pi calculee par le processus present
	printf("La valeur de pi calculee par le processus %d: %f\n", rank, pi_part);
	double sum = 0; // La somme des valeurs de pi_part
	double pi = 0; // La valeur finale de pi = sum/nbp
	if(rank==0){
		if(nb_point==0){
			MPI_Finalize();
			return EXIT_SUCCESS;
		}
		sum += pi_part;
		for(int i=1;i<nbp;i++){
			MPI_Recv(&pi_part, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += pi_part;
		}
		pi = sum/nbp;
		printf("La valeur moyenne de pi est %f. (il y a %d points au total.)\n", pi, nb_point);
	}
	else{
		MPI_Send(&pi_part, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
