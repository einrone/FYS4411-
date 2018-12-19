#include "parameters.h"


/*
Loads the parameters from a parameter txt file. This class is static, so that the parameters only need to be loaded
once. The parameters in then retrieved by using Parameters::name_of_variables.

In the txt file the parameters are given simply as

parameter_name parameter_value #Comment

'#' works as a marker for comments (a simple space between the value and comment may work
but use # to be safe.

The parameter name has to coinside with the names in this class. If unknow name is found in the
file, or a parameter if not set in the file, an error will be given.
*/

void Parameters::read_parameters(std::string location){
    std::ifstream infile;
    std::string line;

    infile.open(location);
    if(infile.fail()){
        std::cout << "Could not find parameter file at " << location << std::endl;
        exit(EXIT_FAILURE);
    }
    infile.clear();
    infile.seekg(0,std::ios::beg);
    for (std::string l; getline(infile,l);)
    {
        std::stringstream ss(l);

        std::string name;
        double variable;

        ss >> name >> variable;
        //std::cout << "Ting: " << name[0] << std::endl;

        if (name.front() == '#'){
            continue;
        }
        else if(name == "MC_cycles"){
            MC_cycles = variable;
            if(MC_cycles>0){
                MC_cycles_set=true;
            }
        }
        else if(name == "N"){
            N = variable;
            if(N>0){
                N_set=true;
            }
        }

        else if(name == "dimensions"){
            dimension = variable;
            if(dimension>0 && dimension<4){
                dimension_set=true;
            }
        }

        else if(name == "alpha_min"){
            alpha_min = variable;
            alpha_min_set=true;
        }
        else if(name == "alpha_max"){
            alpha_max = variable;
            alpha_max_set=true;
        }
        else if(name == "alpha_num"){
            alpha_num = variable;
            alpha_num_set=true;
        }

        else if(name == "beta"){
            beta = variable;
            beta_set=true;
        }
        else if(name == "omega"){
            omega = variable;
            omega_set=true;
        }

        else if(name == "omega_z"){
            omega_z=variable;
            omega_z_set=true;
        }

        else if(name == "a"){
            a=variable;
            a_set=true;
        }

        else if(name == "dx"){
            dx=variable;
            dx_set=true;
        }
        else if(name == "D"){
            D=variable;
            D_set=true;
        }
        else if(name == "numerical"){
            numerical=variable;
            numerical_set=true;
        }


        else{
            std::cout << "Unknonw Variable found: " << name << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if(!MC_cycles_set){
        std::cout << "MC cycles not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!dimension_set){
        std::cout << "Dimensions not set or invalid!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!N_set){
        std::cout << "Number of particles not set or invalid!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!alpha_max_set){
        std::cout << "Max alpha not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!alpha_min_set){
        std::cout << "Min alpha not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!alpha_num_set){
        std::cout << "Numbers of alpha not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!beta_set){
        std::cout << "Beta not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!omega_set){
        std::cout << "Omega not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!omega_z_set){
        std::cout << "Omega_z not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!a_set){
        std::cout << "a not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!dx_set){
        std::cout << "dx not set!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!D_set){
        std::cout << "D not set!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!numerical_set){
        std::cout << "numerical not set!" << std::endl;
        exit(EXIT_FAILURE);
    }


}



bool Parameters::MC_cycles_set=false;
bool Parameters::dimension_set=false;
bool Parameters::N_set=false;
bool Parameters::alpha_max_set=false;
bool Parameters::alpha_min_set=false;
bool Parameters::alpha_num_set=false;
bool Parameters::beta_set=false;
bool Parameters::omega_set=false;
bool Parameters::omega_z_set=false;
bool Parameters::a_set=false;
bool Parameters::dx_set=false;
bool Parameters::D_set=false;
bool Parameters::numerical=false;

int Parameters::MC_cycles = 0;
int Parameters::dimension=0;
int Parameters::N=0;
double Parameters::alpha_min = 0;
double Parameters::alpha_max = 0;
int Parameters::alpha_num = 0;
double Parameters::beta = 0;
double Parameters::omega = 0;
double Parameters::omega_z = 0;
double Parameters::a = 0;
double Parameters::dx=0;
double Parameters::D=0;
bool Parameters::numerical_set=false;




