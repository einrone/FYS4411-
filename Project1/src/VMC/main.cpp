#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "Parameters/parameters.h"
#include "simulation.h"
#include "system.h"
#include "mpi.h"

using namespace std;

int main(int nargs, char *args[])
{
    //Init MPI
    int numprocs, my_rank;
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    double StartTime = MPI_Wtime();


    //Enables Eigen to do matrix operations in parallel
    Eigen::initParallel();

    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");

    cout << "MC steps, particles: " << Parameters::MC_cycles << ", " << Parameters::N << endl;

    //Initiates system and simulation
    System * system = new System();
    Simulation * simulation = new Simulation(system);
    simulation->initiate();


    //Code for using gradient decent to find the optimal alpha, and gets the energy for that alpha
    /*double optimal_alpha = simulation->gradient_descent(0.4);
    std::cout << "Correct a: " << optimal_alpha <<std::endl;
    std::cout << "Running simulation with optimal alpha." << std::endl;
    simulation->run(my_rank,optimal_alpha);*/


    //Runs the simulation for the alpha interval defined in the Parameters class
    simulation->run(my_rank);

    //Finds the onebody density for a given alpha
    //simulation->oneBodyDensity(0.5,0.05,0.,6.);




    double EndTime = MPI_Wtime();
    double TotalTime = EndTime-StartTime;

    if ( my_rank == 0 )  cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;

    MPI_Finalize ();


    return 0;
}
