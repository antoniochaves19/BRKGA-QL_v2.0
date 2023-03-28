#ifndef _DEFINE_H
#define _DEFINE_H
#include <vector>

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF BRKGA-QL --------


std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

// Input File
int debug = 1;                              // 0 - run mode      		    1 - debug mode
int MAXTIME = 1;                            // maximum runtime
int MAXRUNS =  1;                           // maximum number of runs of the method
unsigned MAX_THREADS = 1;            		// number of threads
float OPTIMAL = 0;                          // optimal solution (if it is known)
int numDecoders;                            // number of decoders
int numLS;                                  // 0 - without local search     > k - number of local search heuristics


// Run
char instance[100];                         // name of instance
int numLP = 0;                              // number of LP calls

// computational time (unix systems)
struct timeval Tstart, Tend, Tbest;   

//BRKGA
int n;                                      // size of cromossoms
int p;          	                        // size of population
double pe;              	                // fraction of population to be the elite-set
double pm;          	                    // fraction of population to be replaced by mutants
double rhoe;             	                // probability that offspring inherit an allele from elite parent

double sigma;                               // pearson correlation factor

const double PI = 3.14159265;               // pi

std::vector <TSol> Pop;                     // current population
std::vector <TSol> PopInter;               	// intermediary population

TSol bestSolution;                          // best solution found in the A-BRKGA


// Sort TSol by objective function
bool sortByFitness(const TSol &lhs, const TSol &rhs) { return lhs.ofv < rhs.ofv; }


// Q-Learning parameters
double epsilon;                             // greed choice possibility
double lf;                                  // learning factor
double df;                                  // discount factor
double R;                                   // reward
double qTotal;                              // q*

// list of actions RL
int sizeP []    = {233, 377, 610, 987, 1597, 2584};
double Pe[]     = {0.10, 0.15, 0.20, 0.25, 0.30}; 
double Pm[]     = {0.01, 0.02, 0.03, 0.04, 0.05}; 
double Rhoe[]   = {0.55, 0.60, 0.65, 0.70, 0.75, 0.80}; 
 

// number of parameters in Q-table
const int par = 4;

// actions
int a0 = 0,                                 // p (current action)
    a1 = 0,                                 // pe (current action)
    a2 = 0,                                 // pm (current action)
    a3 = 0;                                 // rhoe (current action)

float Qmax = 0;

std::vector <std::vector <TQ> > Q;          // Q-Table

// QL parameters (auxiliar)
float epsilon_max = 1.0;                    // maximum epsilon 
float epsilon_min = 0.1;                    // minimum epsilon
int Ti = 1;                                 // number of epochs performed
int restartEpsilon = 1;                     // number of restart epsilon

#endif