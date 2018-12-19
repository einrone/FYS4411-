#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <cmath>
#include "system.h"
#include "Parameters/parameters.h"
#include "DataDump/datadump.h"


class Simulation
{
public:
    Simulation(System *m_system);
    void initiate();
    void run(int rank);
    void run(int rank,double alpha);

    double gradient_descent(double alpha_0);
    double compute_local_energy_derivative(double alpha);
    void oneBodyDensity(double optimal_alpha, double r_step, double r_min, double r_max);
private:
    System *system;
    double alpha_step;
    double alpha_min;
    double alpha_max;

    double total_energy;

    int MC_cycles;
    int N;

    double energy;


};

#endif // SIMULATION_H
