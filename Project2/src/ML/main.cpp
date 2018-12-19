#include <iostream>
#include "Eigen/Dense"
#include "system.h"
#include "Parameters/parameters.h"
#include "simulation.h"
#include <ctime>

using namespace std;

void distribute_weights_and_biases(Eigen::ArrayXd &);

int main(int argc, char *argv[])
{
    //Enables Eigen to do matrix operations in parallel
    Eigen::initParallel();

    //Starts the clock
    std::clock_t begin = clock();





    //Reads the parameter file
    Parameters::read_parameters("../input/parameters.txt");



    //Makes a random initial vector holding the parameters (biases and weights)
    Eigen::ArrayXd test_parameters = Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);

    //Distributes the parameters with normal distribution.
    distribute_weights_and_biases(test_parameters);


    System * system = new System();
    Simulation * simulation = new Simulation(system);


    //Optimize the parameters with stochastic descent.
    Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);


    //Due some strange bug, the first time run() is run it may return an extremly large answer
    //But is run() is run one more time, this will not happen. We have yet to find the cause...
    simulation->run(done);
    simulation->run(done);

    //Takes the time and deletes pointers.
    std::clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << "The Simulation took " << elapsed_secs << " sec." << std::endl;

    delete simulation;
    delete system;


    //For testing different learning rates and Ns
    /*

    double Ns[5] = {1,2,3,4,5};
    double rates[5] = {0.32,0.31,0.3,0.2,0.1};


    for(int i = 0;i<5;i++){
        for(int j = 0; j<5;j++){
            std::cout << Ns[j] << " "<< rates[i] << std::endl;
            Parameters::N = Ns[j];
            Parameters::learning_rate = rates[i];
            System * system = new System();
            Simulation * simulation = new Simulation(system);
            Eigen::ArrayXd test_parameters = 0.01*Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
            distribute_weights_and_biases(test_parameters);
            Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);
            delete simulation;
            delete system;
        }

    }
    */


    //For testing different sigmas. Not used in raport
    /*
    double sigma = 0.5;
    for(int i = 0;i<50;i++){
        Parameters::sigma = sigma;
        System * system = new System();
        Simulation * simulation = new Simulation(system);
        Eigen::ArrayXd test_parameters = Eigen::ArrayXd::Zero(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
        distribute_weights_and_biases(test_parameters);
        Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);
        simulation->run(done);
        simulation->run(done);
        delete simulation;
        delete system;

        sigma += 0.02;
    }
    */


    //For testing different dx
    /*
    //double dx[7] = {1.5,1.25,1,0.75,0.5,0.25,0.1}; //For brute force
    double dx[7] = {1,0.5,0.1,0.05,0.01,0.005,0.001}; //For importance
    Eigen::ArrayXd test_parameters = Eigen::ArrayXd::Zero(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
    distribute_weights_and_biases(test_parameters);
    //Eigen::ArrayXd test_parameters = 0.01*Eigen::ArrayXd::Random(Parameters::P*Parameters::dimension + Parameters::N + Parameters::P*Parameters::dimension*Parameters::N);
    for(int i = 0;i<7;i++){
        Parameters::dx = dx[i];
        System * system = new System();
        Simulation * simulation = new Simulation(system);
        Eigen::ArrayXd done = simulation->stochastic_descent(test_parameters);
        simulation->run(done);
        simulation->run(done);
        delete simulation;
        delete system;
    }
    */

    return 0;
}

void distribute_weights_and_biases(Eigen::ArrayXd & array){
    int size = array.size();


    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::normal_distribution<double> distribution(0,0.5);
    for(int i = 0;i<size;i++){
        array[i] = distribution(gen);
    }
}
