#include "system.h"

/*
This code runs the heavy-duty computational parts,
including updating the wavefunction and the energy.
*/
System::System()
{

    r = Eigen::MatrixXd(dimension,P);
    next_r = Eigen::MatrixXd(dimension,P);

    distance.resize(P,P);
    next_distance.resize(P,P);

    quantum_force_vector.resize(Parameters::dimension);
    quantum_force_vector_new.resize(Parameters::dimension);


    a_bias.resize(M);
    b_bias.resize(N);
    weights.resize(M,N);

    X.resize(M);
    X_next.resize(M);
    H.resize(N);

    std::random_device rd;
    gen = std::mt19937_64(rd());
    distribution = std::normal_distribution<double>(0.0,1.0);

    h=1e-6; //Step size for numerical derivative
    number_accept = 0;



    //Most of the difference in the energy between gibbs and MC is due to this factor
    if(gibbs){
        gibbs_factor=0.5;
    }
    else{
        gibbs_factor=1.;
    }

    if(D!=0){
        System::greens_pointer=&System::greens_function_ratio;
    }
    else{
        System::greens_pointer=&System::greens_function_ratio_none;
    }

    //Sets the seed of rand() to the current time
    srand(time(NULL));
}


void System::make_grid(Eigen::ArrayXd &parameters)
{
    /*
    This function initializes everything
    to prepare for a simulation using the biases and
    weights given in the array given as an argument.
    */


    //The biases and weight as held in one 1d vector/array, so this will make a vector for the biases for visible/hidden nodes
    //and a matrix for weights
    Eigen::VectorXd par = (Eigen::VectorXd) parameters;
    Eigen::VectorXd w_flatten(M*N);
    for(int i = 0;i<par.size();i++){
        if(i<M){
            a_bias(i) = par(i);
        }else if(i>=M && i<(N+M)){
            b_bias(i-(M)) = par(i);
        }else{
            w_flatten(i-(M+N)) = par(i);
        }
    }

    Eigen::Map<Eigen::MatrixXd> M2(w_flatten.data(),M,N);

    weights = M2;


    number_accept = 0;
    distribute_particles();
    X_next = X;
    update();
    wavefunction_value=get_wavefunction();

}

void System::distribute_particles(){
   /*
   Distributes particles
   */

    for(int i = 0;i<P*dimension;i++){
            if(D!=0){
                X(i) = 0.5*distribution(gen)*sqrt(dx);
            }else{
                X(i) = 0.5*distribution(gen);
        }
    }

}

void System::update_X_next(int move){
    //Moves one particles

    double random_nr = 0;
    if(D!=0){
        quantum_force(move);
    }


    for(int i = 0; i<dimension; i++){
        if(D!=0){
            random_nr = sqrt(dx)*distribution(gen) + quantum_force_vector(i)*dx*D;
        }else{
            random_nr = dx*2*(static_cast<double>(rand())/RAND_MAX - 0.5);
        }

        X_next(move*dimension+i) = X(move*dimension+i) +  random_nr;
    }
 }



void System::update(){
    /*
    Updates the distance
    maxtrix.
    */

    double dist = 0;
    for(int i = 0; i<P;i++){
        for(int j = 0;j<i;j++){
            for(int dim=0;dim<dimension;dim++){
                dist+=(X(i*dimension+dim)-X(j*dimension+dim))*(X(i*dimension+dim)-X(j*dimension+dim));
            }

            dist=sqrt(dist);
            distance(i,j) = dist;
            distance(j,i) = dist;
            dist = 0;
        }
        temp_value=0;
    }


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

    update_X_next(move);

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
        wavefunction_value = get_wavefunction();
        X = X_next;

        number_accept++;

    }
    else{

        X_next = X;

    }
    if(is_interacting){
        update();
    }

    if(is_numerical){
        return calculate_energy_numerically();
    }
    else{
        return get_local_energy();
    }


}


//Function for gibbs sampling
double System::gibbs_sample_and_return_energy(){
    sample_h();
    sample_x();
    if(is_interacting){
        update();
    }
    if(is_numerical){
        return calculate_energy_numerically();
    }
    else{
        return get_local_energy();
    }

}




double System::get_probability_ratio(int move){
    /*
    Returns probability ratio
    */

    double wavefunction_old=get_wavefunction();

    double wavefunction_new=get_wavefunction_next();

    return (wavefunction_new*wavefunction_new)/(wavefunction_old*wavefunction_old)*greens_factor(move);

}





double System::get_wavefunction(){
    /*
    Computes the wavefunction (complete)
    */
    double wave_function_first_part = 0;
    double wave_function_second_part = 1;
    double exp_factor=0;

    for(int i=0;i<M;i++){
        wave_function_first_part+=-(X(i) - a_bias(i))*(X(i) - a_bias(i));
    }

    wave_function_first_part = exp(wave_function_first_part/(2*sigma_squared));

    for(int j = 0;j<N;j++){
        exp_factor=0;
        for(int i=0;i<M;i++){
            exp_factor+=X(i)*weights(i,j);
        }
        wave_function_second_part *= (1+exp(b_bias(j)+(1.0/sigma_squared)*exp_factor));
    }

    if(gibbs){
        return sqrt(wave_function_first_part*wave_function_second_part);
    }
    else{
        return wave_function_first_part*wave_function_second_part;
    }
}


double System::get_wavefunction_next(){
    /*
    Computes the next wavefunction (complete)
    */
    double wave_function_first_part = 0;
    double wave_function_second_part = 1;
    double exp_factor=0;

    for(int i=0;i<M;i++){
        wave_function_first_part+=-(X_next(i) - a_bias(i))*(X_next(i) - a_bias(i));
    }

    wave_function_first_part = exp(wave_function_first_part/(2*sigma_squared));

    for(int j = 0;j<N;j++){
        exp_factor=0;
        for(int i=0;i<M;i++){
            exp_factor+=X_next(i)*weights(i,j);
        }
        wave_function_second_part *= (1+exp(b_bias(j)+(1.0/sigma_squared)*exp_factor));
    }

    if(gibbs){
        return sqrt(wave_function_first_part*wave_function_second_part);
    }
    else{
        return wave_function_first_part*wave_function_second_part;
    }
}



double System::get_probability(){
    double temp_value = get_wavefunction();
    return temp_value*temp_value;
}




//Returns the local energy. Does also compute the derivative of E_L used in the gradien descent.
double System::get_local_energy(){
    double derivative_of_log_psi=0;
    double second_derivative_of_log_psi=0;
    double exp_factor=0;
    double denominator_factor=0;
    double potential_energy=0;
    double repulsive_interaction=0;
    double potential_factor=0.5*omega*omega;
    Eigen::VectorXd x_weight_product(N);
    double derivative_loc=0;



    for(int k=0;k<M;k++){
           derivative_loc=(a_bias(k)-X(k))/(sigma_squared);
           second_derivative_of_log_psi+=-1.0/(sigma_squared);

           x_weight_product=(1.0/sigma_squared)*(weights.transpose()*X);

           for(int j=0;j<N;j++){
                exp_factor=exp(-b_bias(j)-x_weight_product(j));
                derivative_loc+=weights(k,j)/(sigma_squared*(1+exp_factor));
                denominator_factor = sigma_squared*sigma_squared*(1+exp_factor)*(1+exp_factor);
                second_derivative_of_log_psi+=(weights(k,j)*weights(k,j))*exp_factor/denominator_factor;
           }
           derivative_of_log_psi+=derivative_loc*derivative_loc*gibbs_factor*gibbs_factor;

        }

    if(is_interacting){
        for(int i=0;i<P;i++){
            for(int j=0;j<i;j++){
                repulsive_interaction+=1.0/distance(i,j);
            }
        }
    }



    for(int k=0;k<M;k+=dimension){
        for(int dim=0;dim<dimension;dim++)
            potential_energy+=X(k+dim)*X(k+dim);
    }

   return -0.5*(derivative_of_log_psi+second_derivative_of_log_psi*gibbs_factor)+potential_factor*potential_energy+repulsive_interaction;
}



//Function for derivatives of biases and weights
double System::d_psi_da(int k){
    return (X[k]-a_bias[k])/sigma_squared;
}

double System::d_psi_db(int k){
    double exp_factor=0;
    for(int i=0;i<M;i++){
        exp_factor+=X(i)*weights(i,k);
    }
    return gibbs_factor*1.0/(1+exp(-b_bias(k)-(1.0/sigma_squared)*exp_factor));
}

double System::d_psi_dw(int k, int l){
    double exp_factor=0;
    for(int i=0;i<M;i++){
        exp_factor+=X(i)*weights(i,l);
    }
    return gibbs_factor*X(k)/(sigma_squared*(1+exp(-b_bias(l)-(1.0/sigma_squared)*exp_factor)));
}

double System::d_psi_da_log(int k){
    return gibbs_factor*(X(k)-a_bias(k))/(2.0*sigma_squared);
}

double System::d_psi_db_log(int k){
    double exp_factor=0;
    for(int i=0;i<M;i++){
        exp_factor+=X(i)*weights(i,k);
    }
    return 1.0/(2.0*(1+exp(-b_bias(k)-(1.0/sigma_squared)*exp_factor)));
}

double System::d_psi_dw_log(int k, int l){
    double exp_factor=0;
    for(int i=0;i<M;i++){
        exp_factor+=X(i)*weights(i,l);
    }
    return gibbs_factor*X(k)/(2.0*sigma_squared*(1+exp(-b_bias(l)-(1.0/sigma_squared)*exp_factor)));
}




//Computes the quantum force
void System::quantum_force(int move){
    double grad_value = 0;
    double grad_value_next = 0;
    double exp_factor = 0;
    double exp_factor_next = 0;
    Eigen::VectorXd x_weight_product(N);
    Eigen::VectorXd x_weight_product_next(N);


    x_weight_product=(1.0/sigma_squared)*(weights.transpose()*X);
    x_weight_product_next=(1.0/sigma_squared)*(weights.transpose()*X_next);


    for(int dim=0;dim<dimension;dim++){
        grad_value+=(a_bias(move*dimension+dim)-X(move*dimension+dim));
        grad_value_next+=(a_bias(move*dimension+dim)-X_next(move*dimension+dim));
        grad_value/=sigma_squared;
        grad_value_next/=sigma_squared;
        for(int j=0;j<N;j++){
            exp_factor=exp(-b_bias(j)-x_weight_product(j));
            exp_factor_next=exp(-b_bias(j)-x_weight_product_next(j));
            grad_value+=weights(move,j)/(sigma_squared*(1+exp_factor));
            grad_value_next+=weights(move,j)/(sigma_squared*(1+exp_factor_next));
       }
       quantum_force_vector(dim)=gibbs_factor*grad_value;
       quantum_force_vector_new(dim)=gibbs_factor*grad_value_next;
       grad_value=0;
       grad_value_next=0;

    }
}

//Returns the ratio of the green functions
double System::greens_function_ratio(int move)
{
    double value=0;
    Eigen::VectorXd x_move(dimension);
    Eigen::VectorXd x_move_next(dimension);

    for(int dim=0;dim<dimension;dim++){
        x_move(dim)=X(dimension*move+dim);
        x_move_next(dim)=X_next(dimension*move+dim);

    }


    quantum_force(move);

    value = -(x_move - x_move_next -D*dx*quantum_force_vector_new).squaredNorm()+(x_move_next-x_move - D*dx*quantum_force_vector).squaredNorm();
    value /= 4*D*dx;

    return exp(value);
}

double System::greens_factor(const int move){
    return (this->*greens_pointer)(move);
}

double System::greens_function_ratio_none(int move){
    return 1.0;
}



//Calculates the energy numerically. Does NOT calculate the derivative of E_L
double System::calculate_energy_numerically(){
    double wavefunction_value_plus = 0;
    double wavefunction_value_minus = 0;

    double kinetic_energy = 0;
    double potential_energy = 0;
    double omega_ratio = 1;// omega_z/omega;
    wavefunction_value=get_wavefunction();
    for(int i = 0; i<M;i++){
            potential_energy+=omega*omega*X(i)*X(i);


            X(i)+=h;
            wavefunction_value_plus=get_wavefunction();
            X(i)-=2*h;
            wavefunction_value_minus=get_wavefunction();
            X(i)+=h;
            kinetic_energy -= (wavefunction_value_plus+wavefunction_value_minus - 2*wavefunction_value)/h/h/wavefunction_value;

    }
    return (0.5*potential_energy+0.5*(kinetic_energy));
}




//Sampling of Hidden Nodes
void System::sample_h(){
    std::uniform_real_distribution<double> sampling_distribution(0,1);

    Eigen::VectorXd backwards_H = -b_bias - (X.transpose()*weights).transpose()/sigma_squared;

    double coin_toss = 0;
    double sample_prob = 0;
    for(int i = 0;i<N;i++){
        sample_prob = 1./(1+exp(backwards_H(i)));
        coin_toss = 2*(static_cast<double>(rand())/RAND_MAX - 0.5);
        if(sample_prob > coin_toss){
            H(i) = 1.;
        }
        else{
            H(i) = 0.;
        }
    }
}


//Sampling of Visible Nodes
void System::sample_x(){
    std::normal_distribution<double> sampling_distribution;
    Eigen::VectorXd backwards_X = a_bias + weights*H;

    for(int i = 0;i<M;i++){
        sampling_distribution = std::normal_distribution<double>(backwards_X(i),sigma);
        X(i) = sampling_distribution(gen);
    }
}




