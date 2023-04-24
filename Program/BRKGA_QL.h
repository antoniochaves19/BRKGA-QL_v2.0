#ifndef _BRKGA_QL_H
#define _BRKGA_QL_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <cstring>
#include <string.h>
#include <algorithm>
#include <sys/time.h>
#include <utility>  // pair
#include <numeric>  // iota
#include <map>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip> //graph
#include <sstream> //graph
#include <fstream> //graph
#include <numeric>

#include "Data.h"
#include "Define.h"
#include "Output.h"

#include "Problem.h"


//****************************** General Functions **********************************

/************************************************************************************
 Method: BRKGA_QL()
 Description: Apply the method BRKGA_QL to solve the problem
*************************************************************************************/
void BRKGA_QL();

/************************************************************************************
 Method: updateBestSolution()
 Description: Update the best solution found during the run
*************************************************************************************/
void updateBestSolution(TSol s);

/************************************************************************************
 Method: InitiateQTable()
 Description: Initiate the Q-Table with random values
*************************************************************************************/
void InitiateQTable();

/************************************************************************************
 Method: ChooseAction()
 Description: Choose actions and update the parameters
*************************************************************************************/
void ChooseAction();

/************************************************************************************
 Method: UpdatePopulationSize()
 Description: Update the population size with new value of p
*************************************************************************************/
void UpdatePopulationSize();

/************************************************************************************
 Method: UpdateQLParameters(currentTime)
 Description: Update the parameters of the Q-Learning method
*************************************************************************************/
void SetQLParameters(float currentTime);

/************************************************************************************
 Method: UpdateQTable()
 Description: Update the values of Q-Table
*************************************************************************************/
void UpdateQTable();
 
/************************************************************************************
 Method: CREATE INITIAL SOLUTIONS
 Description: create a initial chromossom with random keys
*************************************************************************************/
TSol CreateInitialSolutions();

/************************************************************************************
 Method: ChaoticInd
 Description: create a solution between a mutant and a elite individual
*************************************************************************************/
void ChaoticInd(TSol &s);

/************************************************************************************
 Method: PARAMETRICUNIFORMCROSSOVER
 Description: create a new offspring with parametric uniform crossover
*************************************************************************************/
TSol ParametricUniformCrossover(int elitesize, int popSize);

/************************************************************************************
 Method: PEARSON CORRELATION
 Description: calculate the Pearson correlation coefficient between two chromossoms
*************************************************************************************/
double PearsonCorrelation(std::vector <TVecSol> s1, std::vector <TVecSol> s2);

/************************************************************************************
 Metodo: IC(TSol Pop)
 Description: apply clustering method to find promising solutions in the population
*************************************************************************************/
void IC();

/************************************************************************************
 Method: LP
 Description: Apply Label Propagation to find communities in the population
*************************************************************************************/
void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas);

/************************************************************************************
 Method: PROMISINGLP
 Description: Find the promising solutions to represent the communities
*************************************************************************************/
void PromisingLP();

/************************************************************************************
Method: FREE MEMORY
Description: free memory of global vector
*************************************************************************************/
void FreeMemory();

/************************************************************************************
 Method: RANDOMICO
 Description: Generate a double random number between min and max
*************************************************************************************/
double randomico(double min, double max);

/************************************************************************************
 Method: IRANDOMICO
 Description: Generate a int random number between min and max
*************************************************************************************/
int irandomico(int min, int max);

#endif