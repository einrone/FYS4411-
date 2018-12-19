#include "parameters.h"


/*
Loads the parameters from a parameter txt file. This class is static, so that the parameters only need to be loaded
once. The parameters in then retrieved by using Parameters::name_of_variables.

In the txt file the parameters are given simply as

parameter_name parameter_value #Comment

'#' works as a marker for comments (a simple space between the value and comment may work
but use # to be safe.

The parameter name has to coincide with the names in this class. If unknow name is found in the
file, or a parameter is not set in the file, an error will be given.
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

        if (name.front() == '#'){
            continue;
        }
        else if(name == "MC_cycles"){
            MC_cycles = variable;

            if(MC_cycles>0){
                MC_cycles_set=true;
            }
        }

        else if(name == "P"){
            P = variable;
            if(P>0){
                P_set=true;
            }
        }

        else if(name == "dimensions"){
            dimension = variable;
            if(dimension>0 && dimension<4){
                dimension_set=true;
            }
        }

        else if(name == "N"){
            N = variable;
            if(N>0){
                N_set=true;
            }
        }

        else if(name == "D"){
            D=variable;
            D_set=true;
        }

        else if(name == "omega"){
            omega = variable;
            omega_set=true;
        }

        else if(name == "sigma"){
            sigma = variable;
            sigma_set=true;
        }


        else if(name == "dx"){
            dx=variable;
            dx_set=true;
        }

        else if(name == "interacting"){
            interacting=variable;
            interacting_set=true;
        }

        else if(name == "numerical"){
            numerical=variable;
            numerical_set=true;
        }

        else if (name == "gibbs"){
            gibbs=variable;
            gibbs_set=true;
        }

        else if (name == "learning_rate"){
            learning_rate=variable;
            learning_rate_set=true;
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

    else if(!P_set){
        std::cout << "P not set or invalid!" << std::endl;
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
    else if(!omega_set){
        std::cout << "Omega not set!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!sigma_set){
        std::cout << "Sigma not set!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!D_set){
        std::cout << "D not set!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!dx_set){
        std::cout << "dx not set!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!interacting_set){
        std::cout << "interacting not set!" << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(!numerical_set){
        std::cout << "numerical not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!gibbs_set){
        std::cout << "gibbs not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

    else if(!learning_rate_set){
        std::cout << "Learning rate not set!" << std::endl;
        exit(EXIT_FAILURE);
    }

}



bool Parameters::MC_cycles_set=false;
bool Parameters::P_set=false;
bool Parameters::N_set=false;
bool Parameters::dimension_set=false;
bool Parameters::omega_set=false;
bool Parameters::sigma_set=false;
bool Parameters::D_set=false;
bool Parameters::dx_set=false;
bool Parameters::interacting_set=false;
bool Parameters::numerical_set=false;
bool Parameters::gibbs_set=false;
bool Parameters::learning_rate_set=false;

int Parameters::MC_cycles = 0;
int Parameters::P = 0;
int Parameters::dimension=0;
int Parameters::N=0;
double Parameters::omega = 0;
double Parameters::sigma = 0;
double Parameters::D=0;
double Parameters::dx=0;
double Parameters::learning_rate = 0;
bool Parameters::interacting=false;
bool Parameters::numerical=false;
bool Parameters::gibbs=false;





