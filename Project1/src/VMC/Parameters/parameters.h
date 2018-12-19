#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Parameters
{
public:
    static Parameters* m_instance;

    static void read_parameters(std::string location);
    static int MC_cycles;
    static bool MC_cycles_set;

    static int dimension;
    static bool dimension_set;


    static int N;
    static bool N_set;


    static double alpha_max;
    static bool alpha_max_set;
    static double alpha_min;
    static bool alpha_min_set;
    static int alpha_num;
    static bool alpha_num_set;

    static double beta;
    static bool beta_set;

    static double omega;
    static bool omega_set;
    static double omega_z;
    static bool omega_z_set;

    static double a;
    static bool a_set;

    static double dx;
    static bool dx_set;

    static double D;
    static bool D_set;

    static bool numerical;
    static bool numerical_set;


};

#endif // PARAMETERS_H
