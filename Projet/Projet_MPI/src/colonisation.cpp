#include <cstdlib>
#include <string>
#include <iostream>
#include <SDL2/SDL.h>        
#include <SDL2/SDL_image.h>
#include <fstream>
#include <ctime>
#include <iomanip>      // std::setw
#include <chrono>
#include <mpi.h>

#include "parametres.hpp"
#include "galaxie.hpp"
 
int main(int argc, char ** argv)
{
    MPI_Init( &argc, &argv);
	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
	int nbp;
	MPI_Comm_size(globComm, &nbp);
	int rank;
	MPI_Comm_rank(globComm, &rank);
    
    char commentaire[4096];
    int width, height;
    SDL_Event event;
    SDL_Window   * window;

    parametres param;


    std::ifstream fich("parametre.txt");
    fich >> width;
    fich.getline(commentaire, 4096);
    fich >> height;
    fich.getline(commentaire, 4096);
    fich >> param.apparition_civ;
    fich.getline(commentaire, 4096);
    fich >> param.disparition;
    fich.getline(commentaire, 4096);
    fich >> param.expansion;
    fich.getline(commentaire, 4096);
    fich >> param.inhabitable;
    fich.getline(commentaire, 4096);
    fich.close();
    
    galaxie g(width, height, param.apparition_civ);
    galaxie g_next(width, height);
    
    int nbSamples = (int)(height/(nbp-1));
    int nbRest = height-(nbp-1)*nbSamples;
    
    if(rank==0){
		std::cout << "Resume des parametres (proba par pas de temps): " << std::endl;
		std::cout << "\t Chance apparition civilisation techno : " << param.apparition_civ << std::endl;
		std::cout << "\t Chance disparition civilisation techno: " << param.disparition << std::endl;
		std::cout << "\t Chance expansion : " << param.expansion << std::endl;
		std::cout << "\t Chance inhabitable : " << param.inhabitable << std::endl;
		std::cout << "Proba minimale prise en compte : " << 1./RAND_MAX << std::endl;
		std::srand(std::time(nullptr));
    
		int deltaT = (20*52840)/width;
		std::cout << "Pas de temps : " << deltaT << " années" << std::endl;

		std::cout << std::endl;
		
		std::chrono::time_point<std::chrono::system_clock> start, end1, end2;
		unsigned long long temps = 0;
		SDL_Init(SDL_INIT_TIMER | SDL_INIT_VIDEO);

		window = SDL_CreateWindow("Galaxie", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                              width, height, SDL_WINDOW_SHOWN);
		galaxie_renderer gr(window);
		
		int count=0;
		std::chrono::duration<double> sum1,sum2,avg1,avg2;

		while (1) {
			start = std::chrono::system_clock::now();

			for(int i=1;i<nbp;i++){
				if(nbp==2) MPI_Recv(&g.data()[0], width*height, MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				else{
					if(i==1) MPI_Recv(&g.data()[0], width*(nbSamples+1), MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					else if(i==nbp-1) MPI_Recv(&g.data()[width*(height-nbRest-nbSamples-1)], width*(nbSamples+nbRest+1), MPI_BYTE, nbp-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					else if(i!=1 && i!=nbp-1) MPI_Recv(&g.data()[width*((i-1)*nbSamples-1)], width*(nbSamples+2), MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}
			end1 = std::chrono::system_clock::now();
			gr.render(g);
			
			end2 = std::chrono::system_clock::now();
			
			std::chrono::duration<double> elaps1 = end1 - start;
			std::chrono::duration<double> elaps2 = end2 - end1;
        
			temps += deltaT;
			std::cout << "Temps passe : "
                  << std::setw(10) << temps << " années"
                  << std::fixed << std::setprecision(3)
                  << "  " << "|  CPU(ms) : calcul " << elaps1.count()*1000
                  << "  " << "affichage " << elaps2.count()*1000
                  << "\r" << std::flush;
            
            for(int i=1;i<nbp;i++) MPI_Send(&g.data()[0], width*height, MPI_BYTE, i, 2, MPI_COMM_WORLD);
            
            if(count<100){
				count++;
				sum1+=elaps1;
				sum2+=elaps2;
			}
            if(count==100){
				avg1 = sum1/count;
				avg2 = sum2/count;
				
				std::cout << "en moyenne : CPU(ms) : calcul " << avg1.count()*1000
                  << "  " << "affichage " << avg2.count()*1000
                  << "\r" << std::endl;
                
                count++;
			}
              
			//_sleep(1000);
			if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
				std::cout << std::endl << "The end" << std::endl;
				MPI_Finalize();
				break;
			}
		}
	}
	
	if(rank!=0){
		while (1){
			if(nbp==2){
				mise_a_jour(param, width, height, 0, height, g.data(), g_next.data());
				g_next.swap(g);
				MPI_Send(&g.data()[0], width*height, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
			}
			else{
				if(rank==1){ 
					mise_a_jour(param, width, height, 0, nbSamples+1, g.data(), g_next.data());
					g_next.swap(g);
					MPI_Send(&g.data()[0], width*(nbSamples+1), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
				}
				else if(rank==nbp-1){ 
					mise_a_jour(param, width, height, height-nbRest-nbSamples-1, height, g.data(), g_next.data());
					g_next.swap(g);
					MPI_Send(&g.data()[width*(height-nbRest-nbSamples-1)], width*(nbSamples+nbRest+1), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
				}
				else if(rank!=1 && rank!=nbp-1){
					mise_a_jour(param, width, height, (rank-1)*nbSamples-1, rank*nbSamples+1, g.data(), g_next.data());
					g_next.swap(g);
					MPI_Send(&g.data()[width*((rank-1)*nbSamples-1)], width*(nbSamples+2), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
				}
			}
			MPI_Recv(&g.data()[0], width*height, MPI_BYTE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

    
    SDL_DestroyWindow(window);
    SDL_Quit();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
