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

    static int P;
    static bool P_set;

    static int dimension;
    static bool dimension_set;

    static int N;
    static bool N_set;

    static double omega;
    static bool omega_set;

    static double sigma;
    static bool sigma_set;

    static double D;
    static bool D_set;

    static double dx;
    static bool dx_set;

    static bool interacting;
    static bool interacting_set;

    static bool numerical;
    static bool numerical_set;  

    static bool gibbs;
    static bool gibbs_set;

    static double learning_rate;
    static bool learning_rate_set;
};

#endif // PARAMETERS_H
