#include "system.h"

/*
This code runs the heavy-duty computational parts,
including updating the wavefunction and the energy.
*/
System::System()
{

    r = Eigen::MatrixXd(Parameters::dimension,Parameters::N);
    next_r = Eigen::MatrixXd(Parameters::dimension,Parameters::N);

    distance.resize(Parameters::N,Parameters::N);
    next_distance.resize(Parameters::N,Parameters::N);

    quantum_force_vector.resize(Parameters::dimension);
    quantum_force_vector_new.resize(Parameters::dimension);

    std::random_device rd;
    gen = std::mt19937_64(rd());
    distribution = std::normal_distribution<double>(0.0,1.0);

    h=1e-6; //Step size for numerical derivative
    number_accept = 0;

    if(Parameters::a !=0){
	   //We use some function pointers to avoid some if statements in the functions
           System::wavefunction_function_pointer=&System::update_wavefunction_interacting;
           System::compute_local_energy=&System::get_local_energy_interacting;

       }
     else{

           System::wavefunction_function_pointer=&System::update_wavefunction_noninteracting;
           System::compute_local_energy=&System::get_local_energy_noninteracting;

    }

    //Sets the seed of rand() to the current time
    srand(time(NULL));

}

void System::make_grid(double m_alpha){
    /*
    This function initializes everything
    to prepare for a simulation with variational
    parameter m_alpha.
    */
    alpha = m_alpha;
    number_accept = 0;
    //Sets all positions to a random position [-1,1]
    if(a!=0){
        distribute_particles_interacting();
    }
    else{
        distribute_particles_noninteracting();
    }

    next_r = r;
    update();
    wavefunction_value=get_wavefunction();
}

void System::distribute_particles_interacting(){
    /*
    This function distributes the particles
    in the interacting case, ensuring that none
    of them start out closer than a from one another.
    */
    bool safe_distance = false;
    for(int i = 0;i<N;i++){
        safe_distance = false;
        while(!safe_distance){
            for(int j = 0;j<dimension;j++){
                if(D!=0){
                    r(j,i) = distribution(gen)*sqrt(dx);
                }else{
                    r(j,i) = 2*(static_cast<double>(rand())/RAND_MAX - 0.5);
                }
            }
            safe_distance = true;
            for(int n = 0;n<i;n++){
                if((r.col(n) - r.col(i)).norm() <= a){
                    safe_distance = false;
                    break;
                }
            }
        }
    }
}
void System::distribute_particles_noninteracting(){
   /*
   This distribution is easier than the one above,
   as we do not need to take a minimum distance
   into account.
   */
    for(int i = 0;i<N;i++){
        for(int j = 0;j<dimension;j++){
            if(D!=0){
                r(j,i) = distribution(gen)*sqrt(dx);
            }else{
                r(j,i) = 2*(static_cast<double>(rand())/RAND_MAX - 0.5);
            }
        }
    }
}


void System::update(){
    /*
    Updates the distance
    maxtrix after a new move
    */
    double temp_value = 0;
    for(int i = 0; i<N;i++){
        for(int j = 0;j<i;j++){
            temp_value = (r.col(i)- r.col(j)).norm();
            distance(i,j) = temp_value;
            distance(j,i) = temp_value;
        }
    }

    next_distance = distance;
    if(D!=0){
        quantum_force(0); //Re-initialize quantum force
    }
    update_expectation();
}

void System::make_move_and_update(const int move){
    /*
    Takes a random move as input and 
    proposes the new step.
    */
    if(D!=0){
        quantum_force(move);
    }
    double random_nr = 0;
    for(int i = 0; i<dimension; i++){
        if(D!=0){
            random_nr = sqrt(dx)*distribution(gen) + quantum_force_vector(i)*dx*D;
        }else{
            random_nr = dx*2*(static_cast<double>(rand())/RAND_MAX - 0.5);
        }

        next_r(i,move) = r(i,move) +  random_nr;
    }

    if(D!=0 || a!=0){
        update_next_distance(move);
    }




}

void System::update_next_distance(int move){
    /*
    Updates the next_distance matrix for
    to use in importance samplinf.
    */
    double dist = 0;
    for(int i = 0;i<N;i++){
        dist = (next_r.col(i)- next_r.col(move)).norm();

        next_distance(i,move) = dist;
        next_distance(move,i) = dist;
    }

}

void System::update_expectation(){
    /*
    Used for the gradient descent method
    */
    expectation_derivative=0;
    expectation_derivative_energy=0;
    expectation_local_energy=0;
    expectation_local_energy_squared=0;
    wavefunction_probability=0;
    wavefunction_value=get_wavefunction();
}

double System::check_acceptance_and_return_energy(int move){
    /*
    Performs the Metropolis test
    and updates the energy
    */


    //Random value [0,1]
    double temp_value = (double)rand()/RAND_MAX;

    //If r is less than the acceptance prob, r is updated to the new r
    if(temp_value <= get_probability_ratio(move)){
        update_wavefunction(move);
        r.col(move) = next_r.col(move);

        distance.col(move) = next_distance.col(move);
        distance.row(move) = next_distance.row(move);

        number_accept++;

    }
    else{
        next_r.col(move) = r.col(move);
        next_distance.col(move) = distance.col(move);
        next_distance.row(move) = distance.row(move);

    }
    if(is_numerical){
        return calculate_energy_numerically();
    }
    else{
        return get_local_energy();
    }


}


double System::phi_exponant(const Eigen::VectorXd &r){
    /*
    Smart way to update wavefunction
    */
    temp_value = 0;

    for(int i = 0;i<dimension;i++){
        if(i == 2){
            //Multiplices beta to the z-componant
            temp_value += beta*r(i)*r(i);
        }
        else{
            temp_value += r(i)*r(i);
        }
    }
    return -alpha*temp_value;
}



double System::f(double dist){
    /*
    Computes Jastrow factor
    */
    double function;
    if(dist <= a){
        function = 0;
    }
    else{
        function = 1 - a/dist;
    }

    return function;
}

double System::get_probability_ratio(int move){
    /*
    Smart way to compute probability ratio, 
    which accounts for interacting and importance sampling
    */

    double temp_value = phi_exponant(r.col(move)); //Stores the probability before move
    double temp_value2 = phi_exponant(next_r.col(move)); //Stores the probability of move
    double f_part = 1;
    double green_part = 1;
    if(a!=0){ //Interacting if needed
        f_part *= update_wavefunction_interacting_f(move);
    }
    if(D!=0){ //Importance sampling if needed
        green_part = greens_function_ratio(move);

    }

    return exp(2*(temp_value2-temp_value))*f_part*f_part*green_part;
}

double System::get_wavefunction(){
    /*
    Computes the wavefunction (complete)
    */
    double temp_value = 0; //Stores the exponants of phi
    double f_part = 1;
    Eigen::VectorXd temp_r;
    for(int i = 0;i<N;i++){
        temp_r = r.col(i);
        temp_value += phi_exponant(temp_r);
        if(a != 0){
            for(int other = i+1;other<N;other++){
                f_part*= f(distance(other,i));
            }
        }

    }
    return exp(temp_value)*f_part;
}

void System::update_wavefunction(const int move){
    return (this->*wavefunction_function_pointer)(move);
}

void System::update_wavefunction_noninteracting(const int move){
    /*
    Smart way to update wavefunction given move
    */
    wavefunction_value*=exp(phi_exponant(next_r.col(move))-phi_exponant(r.col(move)));

}

void System::update_wavefunction_interacting(const int move){
    /*
    Smart way to update wavefunction given move
    */
    wavefunction_value*=exp(phi_exponant(next_r.col(move))-phi_exponant(r.col(move)))*update_wavefunction_interacting_f(move);

}


double System::update_wavefunction_interacting_f(const int move){
    /*
    Returns the ratio of the jastrow Factor
    */
    double second_factor_of_psi=1;
    for(int i = 0; i < N; i++){
            if(i != move){
                 second_factor_of_psi*= f(next_distance(i,move))/f(distance(i,move));
        }
    }
    return second_factor_of_psi;
}

double System::get_probability(){
    double temp_value = get_wavefunction();
    return temp_value*temp_value;
}


double System::get_local_energy(){
    return (this->*compute_local_energy)();
}

//Returns the local energy. Does also compute the derivative of E_L used in the gradien descent.
double System::get_local_energy_noninteracting(){
        double total_energy = 0;
        double temp_value = 0;
        double factor1_noB = -2*(dimension)*alpha*N ;
        double factor1_B = -2*alpha*(dimension - 1)*N -  2*alpha*beta*N;
        double factor2 = 4*alpha*alpha;
        double pot_factor = 0.5*omega*omega;
        double r_i_annen = 0;
        double wavefunction_derivative_value=0;
        double omega_ratio = omega_z/omega;


        for(int k = 0;k<N;k++){
            for(int i = 0; i<dimension;i++){
                if(i==2){
                    temp_value += r(i,k)*r(i,k)*beta*beta;
                    wavefunction_derivative_value+=beta*r(i,k)*r(i,k);
                    r_i_annen += omega_ratio*omega_ratio*r(i,k)*r(i,k);
                }
                else{
                    temp_value += r(i,k)*r(i,k);
                    wavefunction_derivative_value+=r(i,k)*r(i,k );
                    r_i_annen += r(i,k)*r(i,k);
                }


            }
        }

        if(dimension >= 3){
            total_energy = factor1_B + factor2*temp_value;
        }
        else{
            total_energy = factor1_noB + factor2*temp_value;
        }

        wavefunction_derivative_value*=-1;

        temp_value=-0.5*total_energy+ pot_factor*r_i_annen;
        expectation_local_energy+=temp_value;
        expectation_local_energy_squared+=temp_value*temp_value;
        expectation_derivative+=wavefunction_derivative_value;
        expectation_derivative_energy+=(wavefunction_derivative_value)*temp_value;


    return temp_value;
}
//u' . This assumes the distance between particles is more than a because of the placement.
double System::udiv(int idx1,int idx2){
    return a/(distance(idx1,idx2)*distance(idx1,idx2) - a*distance(idx1,idx2));;
    }

//u''
double System::udivdiv(int idx1,int idx2){
    return (a*a - 2*a*distance(idx1,idx2))/(
                (distance(idx1,idx2)*distance(idx1,idx2)-
                 a*distance(idx1,idx2))*(distance(idx1,idx2)*distance(idx1,idx2)-a*distance(idx1,idx2)));}

//Returns the horrible local energy for interacting... Does also compute the derivative of E_L
double System::get_local_energy_interacting(){
    double sec_fac = 0;
    double trd_fac = 0;
    double frt_fac = 0;

    double total_energy = 0;
    double temp_value = 0;
    temp_r = Eigen::VectorXd::Zero(dimension);
    Eigen::VectorXd grad_r = Eigen::VectorXd::Zero(dimension);
    Eigen::VectorXd temp_difference_vector = Eigen::VectorXd::Zero(dimension);
    double factor1_noB = -2*(dimension)*alpha*N ;
    double factor1_B = -2*alpha*(dimension -1)*N -2*alpha*beta*N;
    double factor2 = 4*alpha*alpha;
    double pot_factor = 0.5*omega*omega;
    double r_i_annen = 0;
    double wavefunction_derivative_value=0;
    double omega_ratio = omega_z/omega;



    for(int k = 0;k<N;k++){
        for(int i = 0; i<dimension;i++){
            if(i==2){
                temp_value += r(i,k)*r(i,k)*beta*beta;
                r_i_annen += omega_ratio*omega_ratio*r(i,k)*r(i,k);
                wavefunction_derivative_value+=beta*r(i,k)*r(i,k);
            }
            else{
                temp_value += r(i,k)*r(i,k);
                wavefunction_derivative_value+=r(i,k)*r(i,k);
                r_i_annen += r(i,k)*r(i,k);
            }
        }
    }


    for(int k = 0;k<N;k++){
        temp_r = r.col(k);
        for(int dim = 0;dim<dimension;dim++){
            if(dim == 2){
                grad_r(dim) = temp_r(dim)*beta;
            }
            else{
                grad_r(dim) = temp_r(dim);
            }
        }
        grad_r *= -4*alpha;
        for(int j = 0;j<N;j++){
            if(j!=k){
                temp_difference_vector += (r.col(k)-r.col(j))/(distance(k,j))*udiv(k,j);
            }
        }

        sec_fac += grad_r.dot(temp_difference_vector);
    }

    for(int k = 0; k<N;k++){
        for(int j = 0; j<N;j++){
            for(int i = 0; i<N;i++){
                if(k!=i && k!=j){
                    trd_fac += (r.col(k) - r.col(i)).dot(r.col(k) - r.col(j))/(distance(k,j)*distance(k,i))*udiv(k,i)*udiv(k,j);
                }
           }
        }
    }
    for(int k = 0;k<N;k++){
        for(int j = 0;j<N;j++){
            if(j!=k){
                frt_fac += udivdiv(k,j) + 2.0/(distance(k,j))*udiv(k,j);
            }
        }
    }

    if(dimension >= 3){
        total_energy = factor1_B + factor2*temp_value + sec_fac + trd_fac + frt_fac;
    }
    else{
        total_energy = factor1_noB + factor2*temp_value + sec_fac + trd_fac + frt_fac;
    }


    wavefunction_derivative_value*=-1;
    temp_value= -0.5*total_energy+ pot_factor*r_i_annen;
    expectation_local_energy+=temp_value;
    expectation_local_energy_squared+=temp_value*temp_value;
    expectation_derivative+=wavefunction_derivative_value;
    expectation_derivative_energy+=(wavefunction_derivative_value)*temp_value;

    return temp_value;
}

//Computes the quantum force
void System::quantum_force(int move){
    double grad_value = 0;
    double grad_value_new = 0;




    for(int j=0; j<dimension;j++){

        if(a!=0){
            for(int k=0; k<N;k++){
                if(k !=move){
                    if(distance(move,k)>a){
                        grad_value+=2*a*(r(j,move)-r(j,k))/(distance(move,k)*distance(move,k)*distance(move,k) - a*distance(move,k)*distance(move,k)) ;//distance(move,k)*udiv(move,k);
                    }
                    if(next_distance(move,k)>a){
                        grad_value_new+=2*a*(next_r(j,move)-next_r(j,k))/(next_distance(move,k)*next_distance(move,k)*next_distance(move,k) - a*next_distance(move,k)*next_distance(move,k));//next_distance(move,k)*udiv(move,k);
                    }
                 }
            }
        }
        if(j==2){
           quantum_force_vector(j)=-4.0*alpha*beta*r(j,move)+grad_value;
           quantum_force_vector_new(j)=-4.0*alpha*beta*next_r(j,move)+grad_value;

        }
        else{
            quantum_force_vector(j)=-4.0*alpha*r(j,move)+grad_value;
            quantum_force_vector_new(j)=-4.0*alpha*next_r(j,move)+grad_value;
        }
        grad_value=0;
        grad_value_new=0;
    }


}

//Returns the ratio of the green functions
double System::greens_function_ratio(int move)
{
    double value=0;

    quantum_force(move);

    value = -(r.col(move) - next_r.col(move) -D*dx*quantum_force_vector_new).squaredNorm()+(next_r.col(move)-r.col(move) - D*dx*quantum_force_vector).squaredNorm();
    value /= 4*D*dx;

    return exp(value);
}


//Calculates the energy numerically. Does NOT calculate the derivative of E_L
double System::calculate_energy_numerically(){
    double wavefunction_value_plus = 0;
    double wavefunction_value_minus = 0;

    double kinetic_energy = 0;
    double potential_energy = 0;
    double omega_ratio = omega_z/omega;
    wavefunction_value=get_wavefunction();
    for(int i = 0; i<N;i++){
        for(int j = 0; j<dimension;j++){
            if (j==2){
                potential_energy += omega_ratio*omega_ratio*r(j,i)*r(j,i);
            }
            else{
                potential_energy+=omega*omega*r(j,i)*r(j,i);
            }


            r(j,i)+=h;
            wavefunction_value_plus=get_wavefunction();
            r(j,i)-=2*h;
            wavefunction_value_minus=get_wavefunction();
            r(j,i)+=h;
            kinetic_energy -= (wavefunction_value_plus+wavefunction_value_minus - 2*wavefunction_value)/h/h/wavefunction_value;

        }
    }
    return (0.5*potential_energy+0.5*(kinetic_energy));
}





