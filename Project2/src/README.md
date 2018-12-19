# Project 2/src

This folder contains the source code for our project.

## Content

### 1: Analyzing Tools
This folder contains the Python scripts used to analyze the output of our data.

### 2: ML
This folder contains our main C++ code.

### 3: input
This folder contains our parameter file, that can be updated to 
change the parameters in our script.

## Use
### Setting the inputs
The most important parameters to change is
* **MC_cycles** The number of monte carlo steps. Please use *2^n* for some integer *n*
* **N** The number of hidden nodes
* **P** Number of particles.
* **dimension** The number of dimension.
* **numerical** 1 or 0. Determines if the energy will be calculated analytically or numerically.
* **dx** The step length used in both brute force and importance sampling.
* **D** The diffusion constant. When set to *0*, this will disable importance sampling, and use brute force. When using
importance sampling, set this to *0.5*
* **learning_rate** The learning rate used in the stochastic gradient descent-
* **gibbs** 1 or 0. Whether to use Gibbs sampling of Monte Carlo


One can use comments in the parameter file by using *#* in front of the comment.


### Running the Simulation
The program is compiled through QT, by opening the pro file in the QT editor and compiling. Running the code as is will do one simulation using the parameters found in *parameters.txt*. One can also test the simulations with different *learning rates* and *N*, or different *dx*. The instruction of which lines to comment/uncomment to run which type of simulation if found in **main**. One will also need to uncomment some lines in *run()* in *simulation.cpp* if one wants the simulation with different *dx* to output different files for different *dx*. Instruction for which lines to uncomment is found in *simulation.cpp*.


### Analysis
There are three different analysis programs: 

The first will plot the results of the stochastic gradient descent. After doing the simulation, one should only need to run
```
python3 rate_and_N_analyzer.py
```

The second is to make plots for different *dx*. Make sure to (un)comment some lines in this code if one want to use importance sampling (instructions in the code):
```
python3 dx_analyzer.py
```

The last program simply takes the energy from one simulation, calculates the average and error, and print them:
```
python3 energy_analyzer.py
```
Note: When using the **energy_analyzer.py** and **rate_and_N_analyzer.py**, use that *MC_cycles* times *number of process* is on the form *2^n*, for some integer *n>0*. This is because blocking needs data with such a length.
