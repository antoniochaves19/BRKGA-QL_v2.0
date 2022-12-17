
# BRKGA-QL: An Adaptive and near Parameter-free BRKGA using Q-Learning

This is an implementation of the Biased Random-key Genetic Algorithm with Q-Learning (BRKGA-QL) to solve combinatorial optmization problems.

The C++ code of this algorithm has been designed to be easy of reuse. Users can only implement specific functions (read, decoder and local search). 

Here we have the BRKGA-QL version 2.0 code.


## References

When using this algorithm in academic studies, please refer to the following works:

BRKGA-QL version 1.0
[1] Chaves, A.A. and Lorena, L.H.N. (2021)
An Adaptive and near Parameter-free BRKGA using Q-Learning, In: 2021 IEEE Congress on Evolutionary Computation (CEC), 2021, Kraków, Poland. Available in: http://dx.doi.org/10.1109/cec45853.2021.9504766.

BRKGA-QL version 2.0
[2] Chaves et al. (2022)
Hybrid Metaheuristic for an Integrated Dial-a-Ride Problem with Common Carrier and Public Transportation.
Available here in technical report form: https://github.com/antoniochaves19/BRKGA-QL_v2.0/blob/main/CIRRELT-2022-32.pdf

## Scope

This code has been designed to solve the Traveling Salesman Problem (TSP). Users need to configure only Problem.cpp file.


## Running the algorithm

* Enter the Program directory: `cd Program`
* Run the make command: `make`
* Run the BRKGA_QL: `./runTest`

* Or compile via terminal: `g++ -std=c++11 -o runTest BRKGA_QL.cpp Problem.cpp -O3 -openmp` (or -fopenmp if Linux)


## Code structure

The code structure is documented in [1] and organized in the following manner:

* **SPECIFIC_CODE:**
    * **Problem.cpp**: Contains data structure of the problem, the input function, the decoders of the BRKGA-QL, and the local search heuristics.

* **GENERAL_CODE:**
    * **BRKGA_QL.cpp**: Contains all of the BRKGA-QL algorithm's population mechanisms and the main function to start the algorithm.
    * **Data.h**: Represents the data structures of BRKGA and QL.
    * **Define.h**: Stores the global variables of the BRKGA-QL.
    * **Output.h**: Stores the outputs functions, including the best solution found and statistical analysis of the BRKGA-QL.

## File testScenario.csv is the input data problem and each line consists in:

- Instance Name
- Run mode (0 = debug, prints in the screen; 1 = run, prints in files)
- Number of implemented decoders
- Number of implemented local search heuristics (0 if local search is not available)
- Maximum rumtime (in seconds)
- Maximum number of runs
- Number of threads used in OpenMP
- Optimal or lower bound solution (if is known), 0 otherwise

Users need to create a folder named "Instances/ProblemName", where the instances must be.
