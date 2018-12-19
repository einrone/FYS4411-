#ifndef SYSTEM_H
#define SYSTEM_H
#include "Eigen/Dense"
#include "Parameters/parameters.h"
#include <random>
#include <time.h>
#include <math.h>
#include <iostream>

class System
{
public:
    System();

    //Matrices and vectors holding information about the system
    Eigen::MatrixXd r;
    Eigen::MatrixXd next_r;
    Eigen::MatrixXd distance;
    Eigen::MatrixXd next_distance;
    Eigen::VectorXd quantum_force_vector;
    Eigen::VectorXd quantum_force_vector_new;


    //Functions for making and updating the grid/system
    void make_grid(double m_alpha);
    void make_move_and_update(int move);
    void update();
    void update_expectation();
    void update_next_distance(int move);
    void distribute_particles_interacting();
    void distribute_particles_noninteracting();


    //Variables for dE_L/dAlpha
    double expectation_local_energy;
    double expectation_derivative;
    double expectation_derivative_energy;
    double expectation_local_energy_squared;


    //Returns various wavefunctions and energies
    double get_wavefunction();
    double get_probability_ratio(int move);
    double get_probability();
    double get_local_energy();
    double get_local_energy_interacting();
    double get_local_energy_noninteracting();
    const Eigen::MatrixXd get_position() const;
    double check_acceptance_and_return_energy(int);

    //Holds the number of accepted moves
    int number_accept;

    //Saves all the variables from the parameters to save time
    const int N = Parameters::N;
    const int dimension = Parameters::dimension;
    const double beta = Parameters::beta;
    const double dx = Parameters::dx;
    const double omega = Parameters::omega;
    const double omega_z = Parameters::omega_z;
    const double a = Parameters::a;
    const double D = Parameters::D;
    const bool is_numerical = Parameters::numerical;


    //Vector and double used for holding temp values
    Eigen::VectorXd temp_r;
    double temp_value;
    Eigen::MatrixXd *r_temp_pointer;


    //First try on using function pointers to speed up the implementation
    //of importance sampling and interacting. Because of laziness we stop halfway
    //and only used if-statements...
    void (System::*wavefunction_function_pointer)(const int);
    double (System::*compute_energy_numeric)();
    double (System::*compute_local_energy)();


    //Functions for calulations involving wavefunction
    double alpha;
    double wavefunction_value;
    double wavefunction_probability;
    double local_energy;
    double h;


    //Updates energy and wavefunctions
    double phi_exponant(const Eigen::VectorXd &r);
    void update_wavefunction(const int move);
    void update_probability_ratio();
    double udiv(int,int);
    double udivdiv(int,int);
    double update_wavefunction_interacting_f(const int);
    void update_wavefunction_interacting(const int);
    void update_wavefunction_noninteracting(const int move);
    void quantum_force(int move);
    double greens_function_ratio(int move);
    double f(double);
    double calculate_energy_numerically();



    std::random_device rd;
    std::mt19937_64 gen;
    std::normal_distribution<double> distribution;




};

#endif // SYSTEM_H
