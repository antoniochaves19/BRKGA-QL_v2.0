#include "BRKGA_QL.h"

/************************************************************************************
								MAIN FUNCTION AREA
*************************************************************************************/
int main()
{ 
    // file with test instances and input data
	FILE *arqProblems;
    arqProblems = fopen ("testScenario.csv", "r"); 

    if (arqProblems == NULL){
        printf("\nERROR: File testScenario.csv not found\n");
        getchar();
        exit(1);
    }

    char nameTable[100];

    //read first line of arqProblems file
    fgets(nameTable, sizeof(nameTable), arqProblems); 

    // best solution that is saved in out file
    TSol sBest;

	// run the BRKGA-QL for all test instances
	while (!feof(arqProblems))
	{
		// read the name of instances, debug mode, local search module, maximum time, maximum number of runs, maximum number of threads
		fscanf(arqProblems,"%s %d %d %d %d %d %d %f", nameTable, &debug, &numDecoders, &numLS, &MAXTIME, &MAXRUNS, &MAX_THREADS, &OPTIMAL);
        strcpy(instance,nameTable);
        
		//read data of the instance
        ReadData(nameTable, n);

        double foBest = INFINITY,
               foAverage = 0;

        float timeBest = 0,
              timeTotal = 0;

        std::vector <double> ofvs;
        ofvs.clear();

        // best solutions found in MAXRUNS
        sBest.ofv = INFINITY;

		// run BRKGA-QL MaxRuns for each instance
        printf("\n\nInstance: %s \nRun: ", instance);
        for (int j=0; j<MAXRUNS; j++)
        {
            // fixed seed
            if (debug == 1)
                srand(j+1); 
            else
                srand(time(NULL));

            printf("%d ", j+1);
            
            gettimeofday(&Tstart, NULL);
            gettimeofday(&Tend, NULL);
            gettimeofday(&Tbest, NULL);

            // best solution found in this run
            bestSolution.ofv = INFINITY;

            // execute the evolutionary method
            BRKGA_QL();

            gettimeofday(&Tend, NULL);

            // store the best solution found in MAXRUNS
            if (bestSolution.ofv < sBest.ofv)
                sBest = bestSolution;

            // calculate best and average results
            if (bestSolution.ofv < foBest)
                foBest = bestSolution.ofv;

            foAverage += bestSolution.ofv;

            // fitness of each solution found in the runs
            ofvs.push_back(bestSolution.ofv);

            timeBest += ((Tbest.tv_sec  - Tstart.tv_sec) * 1000000u + Tbest.tv_usec - Tstart.tv_usec) / 1.e6;
            timeTotal += ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 
        }

        // create a .csv file with average results
        foAverage = foAverage / MAXRUNS;
        timeBest = timeBest / MAXRUNS;
        timeTotal = timeTotal / MAXRUNS;

        if (!debug)
        {
        	WriteSolution(sBest, n, timeBest, timeTotal, instance);
        	WriteResults(foBest, foAverage, ofvs, timeBest, timeTotal, instance);
        }
        else
        {
            WriteSolutionScreen(sBest, n, timeBest, timeTotal, instance);
        }

        // free memory with problem data
        FreeMemoryProblem();

        // free memory of BRKGA-QL components
        FreeMemory();
    }

    fclose(arqProblems);
    return 0;
}


/************************************************************************************
			                  GENERAL FUNCTIONS
*************************************************************************************/

void BRKGA_QL()
{
    // free memory of BRKGA components
    FreeMemory();
 
    // initialize Q-Table
    InitiateQTable();
    
    // number of restart epsilon
    restartEpsilon = 1;  

    // maximum epsilon  
    epsilon_max = 1.0;  
    
    // initialize population
    Pop.clear(); 
    PopInter.clear(); 

    // define the population size with the higher value of P
    p = sizeP[sizeof(sizeP)/sizeof(sizeP[0]) - 1]; 

    Pop.resize(p);
    PopInter.resize(p);

    // Create the initial chromosomes with random keys
    for (int i=0; i<p; i++)
    {
        TSol ind = CreateInitialSolutions(); 
        Decoder(ind);
        Pop[i] = PopInter[i] = ind;
    }
    
    // sort population in increase order of fitness
    sort(Pop.begin(), Pop.end(), sortByFitness);

    // save the best solution found
    updateBestSolution(Pop[0]);
    
    // useful local variables
    int numGenerations = 0;             // number of generations
    int bestGeneration = 0;             // generation in which found the best solution
    double bestFitness = Pop[0].ofv;    // best fitness found in past generation
    float currentTime = 0;              // computational time of the search process
    int sumLS = 0;                      // number of local search applied in each generation
    int noImprov = 0;                   // number of generations without improvement in the best solution
    
    // run the evolutionary process until stop criterion
    while(1)
    {
    	// number of generations
        numGenerations++;

        // number of generations without improvement in the best solution
        noImprov++;

        // *************************** BRKGA-QL **************************************
        // set Q-Learning parameters                                              //**
        SetQLParameters(currentTime);                                             //**
        //                                                                        //**
        // choose a action (value) for each state (parameter)                     //**
        ChooseAction();                                                           //**
        //                                                                        //**
        // update population size                                                 //**
        UpdatePopulationSize();                                                   //**
        // ***************************************************************************


        // **************************** BRKGA ****************************************
        // define parameters for classic BRKGA                                    //**
        //p       = 987;                                                          //**
        //pe      = 0.20;                                                         //**
        //pm      = 0.05;                                                         //**
        //rhoe    = 0.70;                                                         //** 
        // ***************************************************************************


        // The 'Pe' best chromosomes are maintained, so we just copy these into PopInter:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i=0; i<(int)(p*pe); i++){
            // copy the chromosome for next generation
            PopInter[i] = Pop[i]; 
        }  

        // We'll mate 'P - Pe' pairs; initially, i = p*pe, so we need to iterate until i < p:
        #pragma omp parallel for num_threads(MAX_THREADS)
        for (int i = (int)(p*pe); i < p; i++){            
            // Parametric uniform crossover with mutation
            PopInter[i] = ParametricUniformCrossover((int)(p*pe), p);
 
            // Calculate the fitness of new chromosomes
            Decoder(PopInter[i]); 
        }
                
        // Update the current population
        Pop = PopInter;   

        // Sort population in increase order of fitness
        sort(Pop.begin(), Pop.end(), sortByFitness);
        updateBestSolution(Pop[0]);

        // Reward
        R = 0;

        // We improve the best fitness in the current population 
        if (Pop[0].ofv < bestFitness){

            // The reward function is based on improvement of the current best fitness and binary reward
             R = 1 + 100000*(((bestFitness/Pop[0].ofv)-1)/(p));

            // The immediate reward values are limited by some constant.
            if (R > 2) R = 2;

            bestFitness = Pop[0].ofv;
            bestGeneration = numGenerations;
            noImprov = 0;
        }

        // Print reward of each generation
        if (debug) {
            FILE *arquivo;
            arquivo = fopen("../Results/Reward.txt","a");
            fprintf(arquivo, "\n%d \t %.3lf",numGenerations, R);
            fclose(arquivo);
        }
        
        // Update the Q-Table values
        if (R > 0){
            UpdateQTable();
        }


        // ********************* LOCAL SEARCH IN COMMUNITIES *******************
        sumLS = 0;

        // Verify if there are local search heuristics 
        if (numLS > 0){   

            //apply local search when BRKGA found a new better solution or n*pe generations without improvements
            if (R >= 1 || noImprov > (int)n*pe){

                // restart the count of generations without improvements (local search)
                noImprov = 0;

                // Identify commuties in the Elite with Label Propagation method
	            IC();

	            std::vector <int> promisingSol; 
                promisingSol.clear();

	            for (int i=0; i < (int)(p*pe); i++) {
                    
                    // insert the individual index in the promising list
	                if (Pop[i].promising == 1){
	                	promisingSol.push_back(i);
	                }
                    
                    // generate caotic individual (crossover between one elite and one mutant)
                    else if (i > 0){
                        ChaoticInd(Pop[i]);
                        Decoder(Pop[i]);

                        // set flag as 0 to permit new local search
                        Pop[i].flag = 0;
                    }
	            }

	            #pragma omp parallel for num_threads(MAX_THREADS)
                for (unsigned int i=0; i < promisingSol.size(); i++){

                    // local search not influence the evolutionary process
                    TSol s = Pop[promisingSol[i]];
                    LocalSearch(s);
                    updateBestSolution(s);                    

                    // set flag as 1 to prevent new local search in the same solution
                    Pop[promisingSol[i]].flag = 1;
                }

                sumLS = promisingSol.size();
                promisingSol.clear();

                sort(Pop.begin(), Pop.end(), sortByFitness);
                updateBestSolution(Pop[0]);
	        }
	    }
        // *********************************************************************


        // ******************************* SHAKING *****************************
        if ((numGenerations - bestGeneration) > n) {
            
            if (debug) printf("\n\nShaking elite and reset non-elite...\n\n");

            // reset the number of generations without improvement
            bestGeneration = numGenerations;

            // Shake the elite set
            float shaking_type = 0.0;
            int intensity = n*rhoe;
            for(int e = 0; e < (int)(p*pe); e++) {
                for(int k = 0; k < intensity; k++) {
                    shaking_type = irandomico(1,4);
                    int i = irandomico(0, n - 1);
                    if(shaking_type == 1){
                        // Invert value
                        Pop[e].vec[i].rk = 1.0 - Pop[e].vec[i].rk;
                    }
                    else 
                    if (shaking_type == 2){
                        // Swap two random positions
                        int j = irandomico(0, n - 1);
                        double temp = Pop[e].vec[i].rk;
                        Pop[e].vec[i].rk = Pop[e].vec[j].rk;
                        Pop[e].vec[j].rk = temp;
                    }
                    else
                    if(shaking_type == 3){
                        // Change to random value
                        Pop[e].vec[i].rk = randomico(0,1);
                    }
                    i = irandomico(0, n - 2);
                    if(shaking_type == 4){
                        // Swap with neighbor
                        double temp = Pop[e].vec[i].rk;
                        Pop[e].vec[i].rk = Pop[e].vec[i+1].rk;
                        Pop[e].vec[i+1].rk = temp;
                    }
                }
                Decoder(Pop[e]);
            }

            // reset the non-elite chromosomes
            for (int i=(int)(p*pe); i<p; i++){
                Pop[i] = CreateInitialSolutions();
                Decoder(Pop[i]);
            }

            sort(Pop.begin(), Pop.end(), sortByFitness);
            updateBestSolution(Pop[0]);
            bestFitness = Pop[0].ofv;
        }
        // *********************************************************************

        // print screen 
        if (debug){
            printf("\nGeneration: %3d [%4d - %3d(%.2lf) (%3d) - (%.2lf) - (%.2lf)] \t %.2lf (%.2lf)  \t %.2lf [%.4lf] \t %.4lf \t %.4lf",
                        numGenerations, p, (int)(p*pe), pe, sumLS, pm, rhoe, bestSolution.ofv, bestSolution.vec[n].rk, bestFitness, R, epsilon, lf);
        }

        // terminate the evolutionary process in MAXTIME
        gettimeofday(&Tend, NULL);
        currentTime = ((Tend.tv_sec  - Tstart.tv_sec) * 1000000u + Tend.tv_usec - Tstart.tv_usec) / 1.e6; 
        
        // stop criterium
        if (currentTime >= MAXTIME || bestSolution.ofv <= OPTIMAL){
            break;
        }
    }
}

void Decoder(TSol &s)
{
    // copy the random-key sequence of current solution 
    TSol temp = s;

    // define decoder function based in the random-key of position n+1
    int dec = floor(s.vec[n].rk*numDecoders)+1;
    switch (dec)
    {
        case 1: 
            Dec1(s, n);
            break;

        case 2: 
            Dec2(s, n);
            break;
        
        case 3: 
            Dec3(s, n);
            break;

        case 4: 
            Dec4(s, n);
            break;

        case 5: 
            Dec5(s, n);
            break;

        default:
            break;
    }

    // return initial random-key sequence and maintain the solution sequence
    for (int i=0; i<n; i++){
        s.vec[i].rk = temp.vec[i].rk;
    }
}

void LocalSearch(TSol &s)
{
    // ***** we use a Random Variable Neighborhood Descent (RVND) as local search ****

    // current neighborhood
	int k = 1;

    // predefined number of neighborhood moves
    std::vector <int> NSL;
    std::vector <int> NSLAux;
    
    for (int i=1; i<=numLS; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

	while (!NSL.empty())
	{
        // current objective function
        double foCurrent = s.ofv;

        // randomly choose a neighborhood
        int pos = rand() % NSL.size();
        k = NSL[pos];

        switch (k)
        {
            case 1: 
                LS1(s, n); 
                break;

            case 2:
                LS2(s, n); 
                break;

            case 3:
                LS3(s, n); 
                break;

            case 4:
                LS4(s, n); 
                break;

            default:
                break;
        }

        // we better the current solution
        if (s.ofv < foCurrent)
        {
            // refresh NSL
            NSL.clear();
            NSL = NSLAux;
        }
        else
        {
            // Remove N(n) from NSL
            NSL.erase(NSL.begin()+pos);
        }
	} //end while
}

void updateBestSolution(TSol s)
{
    // save the best solution found 
    if (s.ofv < bestSolution.ofv)
    {
        bestSolution = s;
        gettimeofday(&Tbest, NULL);
    }
}

void InitiateQTable()
{
    // initialize the Q-Table values at 0
    Q.clear();
    Q.resize(par);

    qTotal = 0.0;

    // rho_e
    for (unsigned int j=0; j<sizeof(Rhoe)/sizeof(Rhoe[0]); j++)
    {
        TQ aux;
        aux.S = 0;
        aux.pVar = Rhoe[j];
        aux.q = 0;
        aux.k = 0;
        aux.kImp = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    // p
    for (unsigned int j=0; j<sizeof(sizeP)/sizeof(sizeP[0]); j++)
    {
        TQ aux;
        aux.S = 1;
        aux.pVar = sizeP[j];
        aux.q = 0;
        aux.k = 0;
        aux.kImp = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    // pm
    for (unsigned int j=0; j<sizeof(Pm)/sizeof(Pm[0]); j++)
    {
        TQ aux;
        aux.S = 2;
        aux.pVar = Pm[j];
        aux.q = 0;
        aux.k = 0;
        aux.kImp = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }

    // pe
    for (unsigned int j=0; j<sizeof(Pe)/sizeof(Pe[0]); j++)
    {
        TQ aux;
        aux.S = 3;
        aux.pVar = Pe[j];
        aux.q = 0;
        aux.k = 0;
        aux.kImp = 0;

        Q[aux.S].push_back(aux);
        qTotal += aux.q;
    }                                    
}

void ChooseAction()
{
    // choose actions for each state from Q-Table using epsilon-Greedy policy
    for (int i=0; i<par; i++)
    {
        int a = 0, 
            aAux = 0;

        // set variable a with the current action
        switch (i)
        {
            case 0:
                a = a0;
                break;
        
            case 1:
                a = a1;
                break;
            
            case 2:
                a = a2;
                break; 

            case 3:
                a = a3;
                break;
        }       
                
        // found actions with the highest value of Q(i,-).q 
        double bQ = -INFINITY;  
        for (unsigned int j=0; j<Q[i].size(); j++)
        {
            if (Q[i][j].q > bQ)
            {
                bQ = Q[i][j].q;
                aAux = j;
            }
            else
            if (Q[i][j].q == bQ && randomico(0,1) >= 0.5) // trie randomly
            {
                aAux = j;
            }

            // update the best future reward
            if (Q[i][j].q > Qmax)
                Qmax = Q[i][j].q;
        }

        // epsilon-greedy policy
        if (randomico(0,1) <= 1-epsilon) 
        {
            // choose the action with highest Q value
            a = aAux;
        }
        else
        {
            // choose a randonly selected action (value)
            a = irandomico(0,Q[i].size()-1);
        }

        // set new action
        switch (i)
        {
            case 0:
                a0 = a;
                break;
        
            case 1:
                a1 = a;
                break;
            
            case 2:
                a2 = a;
                break; 

            case 3:
                a3 = a;
                break;
        }

        // update number of choices state i and action a
        Q[i][a].k++;
    }
        
    // update parameters with actions 
    rhoe    = Q[0][a0].pVar;   
    p       = Q[1][a1].pVar;
    pm      = Q[2][a2].pVar;
    pe      = Q[3][a3].pVar;  
}

void SetQLParameters(float currentTime)
{
    // **** define epsilon ****

    // restart epsilon once Ti epochs are performed (Ti is 10% of the runtime)
    Ti = MAXTIME * 0.1;
    if (currentTime >= restartEpsilon * Ti){
        restartEpsilon++;

        // cosine decay with warm restart
        epsilon_max = epsilon_max - 0.1;
        if (epsilon_max < epsilon_min)
            epsilon_max = epsilon_min;
        epsilon = epsilon_max;
    }
    else {
        epsilon = epsilon_min + 0.5 * (epsilon_max - epsilon_min) * (1 + cos((((int)currentTime%Ti)/(float)(Ti))*M_PI));
    }
    
    // *** define learning rate ***

    // initialy, a higher priority is given to the newly gained information (exploration mode)
    // then, we decrement lf and have a higher priority for the existing information in Q-Table (exploitation mode)
    lf = 1 - (0.9 * currentTime / MAXTIME); 

    // *** define discount rate ***

    // we look for a higher, long-term reward
    df = 0.8;
}

void UpdateQTable()
{ 
    // set the current action for each state
    for (int s=0; s<par; s++)
    {
        int a; 

        switch (s)
        {
            case 0:
                a = a0;
                break;
            
            case 1:
                a = a1;
                break;

            case 2:
                a = a2;
                break;

            case 3:
                a = a3;
                break;
        }
        
        qTotal -= Q[s][a].q;

        // Q-Learning
        // Q(s,a) is incremented when the action leads to a state, in which there exists an action such that the best possible Q-value and
        // the reward R is greater than current value of Q(s,a).
        // i.e., the old value of Q(s,a) was too pessimistic 
        // df*Qmax is the target Q-value

        Q[s][a].q = Q[s][a].q + lf*(R + df*Qmax - Q[s][a].q); 
       
        Q[s][a].kImp++;
        qTotal += Q[s][a].q;
    }
}

void UpdatePopulationSize()
{
    // *** define the new population size

    // size of the current population
    int oldPsize = Pop.size();

    // proportional pruning 
    if (oldPsize > p){

        // copy the current population
        PopInter = Pop;

        // define new size of Pop
        Pop.resize(p);

        // select the elite chromosomes
        for (int i=0; i<(int)(p*pe); i++){
            // copy p*pe best chromosomes
            Pop[i] = PopInter[i];
        }

        // select the non-elite chromosomes
        int pos = (int)(pe*oldPsize);
        for (int i=(int)(p*pe); i<p; i++){
            // copy the chromosome
            Pop[i] = PopInter[pos];
            pos++;
        }

        // clean intermediate population
        PopInter.clear();
        PopInter.resize(p);
    }
    
    // generate new chromosomes 
    else if (oldPsize < p){

        // define new size of Pop
        Pop.resize(p);

        for (int k = oldPsize; k < p; k++)
        {
        	Pop[k] = ParametricUniformCrossover((int)(oldPsize*pe), oldPsize-1);
            Decoder(Pop[k]);
        }

        // sort new population
        sort(Pop.begin(), Pop.end(), sortByFitness);
        updateBestSolution(Pop[0]);
        
        // clean intermediate population
        PopInter.clear();
        PopInter.resize(p);
    }
}

TSol CreateInitialSolutions()
{
	TSol s;
	TVecSol aux;

    s.vec.clear();

	// create a random-key for each allelo (consider decoder type in the n-th random-key)
	for (int j = 0; j < n+1; j++)
	{
        aux.rk  = randomico(0,1);  // random value between [0,1[
        aux.sol = 0;
        s.vec.push_back(aux);
	}

    // flag to control the local search memory
    s.flag = 0;

	return s;
}

void ChaoticInd(TSol &s)
{
    // generate a caotic individual
    for (int k=0; k<n+1; k++)
    {      
        if (randomico(0,1) > rhoe)
           s.vec[k].rk = randomico(0,1);
    }

    // set the flag of local search as zero
    s.flag = 0;
}

TSol ParametricUniformCrossover(int eliteSize, int popSize)
{	
	TSol s;

    int eliteParent = irandomico(0, eliteSize - 1);                 // one chromosome from elite set
    int nonEliteParent = irandomico(0, popSize-1);                  // one chromosome from entire population

    // best fit parent is parent elite
    if (Pop[eliteParent].ofv > Pop[nonEliteParent].ofv){
        int temp = eliteParent;
        eliteParent = nonEliteParent;
        nonEliteParent = temp;
    }

	// create a new offspring
	s.vec.resize(n+1);

    // Mate: including decoder gene in the n-th rk 
    for(int j = 0; j < n+1; j++)
    {
        // mutation
        if (randomico(0,1) < pm)
        {
            s.vec[j].rk = randomico(0,1);
        }
        else
        {
            //copy alelos of top chromossom of the new generation
            if (randomico(0,1) < rhoe){
                s.vec[j].rk = Pop[eliteParent].vec[j].rk;
            }
            else{
                s.vec[j].rk = Pop[nonEliteParent].vec[j].rk;
            }
        }
    }

    // set the flag of local search as zero
    s.flag = 0;

    return s;
}

double PearsonCorrelation(std::vector <TVecSol> X, std::vector <TVecSol> Y)
{
    double correlation = 0;
    double sumXY = 0;
    double sumX2 = 0;
    double sumY2 = 0;
    double sumX = 0;
    double sumY = 0;

    for(int j=0; j<n; j++)
    {
        sumX += X[j].rk;
        sumX2 += X[j].rk * X[j].rk;
        sumXY += X[j].rk * Y[j].rk;
        sumY += Y[j].rk;
        sumY2 += Y[j].rk * Y[j].rk;
    }

    //Pearson
    correlation= ((n*sumXY) - (sumX*sumY) ) / (sqrt( (n*sumX2 - sumX*sumX) * (n*sumY2 - sumY*sumY) ));
    return correlation;
}

void IC() 
{
    int Tpe = (int)p*pe;
    std::vector<std::vector<std::pair<int, double> > > listaArestas(Tpe, std::vector<std::pair<int, double> >());

	// create weighted (pearson correlation) graph
	int entrouAresta = 0;
	double pearson = 0.0;
	for (int i = 0; i < Tpe - 1; i++) {
		for (int j = i + 1; j < Tpe; j++)
		{
			pearson = PearsonCorrelation(Pop[i].vec, Pop[j].vec);
			if (pearson > 0.6) {
				entrouAresta++;
				listaArestas[i].push_back(std::make_pair(j, pearson));
				listaArestas[j].push_back(std::make_pair(i, pearson));
			}
		}
	}

	// apply clustering method
	LP(listaArestas);

	PromisingLP();
    listaArestas.clear();
}

void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas)
{
    int nk = listaArestas.size();

	// Create vector with visit order
	std::vector<int> ordemVisita(nk);
	iota(ordemVisita.begin(), ordemVisita.end(), 0);

	// initialize each node with its own label
	for (int i = 0; i < nk; i++)
		Pop[i].label = i;

	int iteracao = 1;
	int labelVizinho, melhorLabel;
	double melhorPeso;
	std::map<int, double> totalLabels;
	std::map<int, double>::iterator it;

	int movimentos = 1;
	while (movimentos) {
		movimentos = 0;
		random_shuffle(ordemVisita.begin(), ordemVisita.end());
		for (auto idVertice : ordemVisita) {

			// Calculate the weigth of the labels
			totalLabels.clear();
			for (auto idVizinho : listaArestas[idVertice]) {
				labelVizinho = Pop[idVizinho.first].label;
				it = totalLabels.find(labelVizinho);
				if (it != totalLabels.end()) {
					it->second += idVizinho.second;
				}
				else {
					totalLabels[labelVizinho] = idVizinho.second;
				}
			}

			// Best label is itself initially
			melhorLabel = Pop[idVertice].label;
			melhorPeso = std::numeric_limits<double>::min();
			for (auto totais : totalLabels) {
				if (totais.second > melhorPeso) {
					melhorLabel = totais.first;
					melhorPeso = totais.second;
				}
			}

			if (melhorLabel != Pop[idVertice].label) {
				Pop[idVertice].label = melhorLabel;
				movimentos = 1;
			}
		}
		iteracao++;
	}

    ordemVisita.clear();
}

void PromisingLP()
{
    int Tpe = (int)p*pe;
    std::vector <int> grupos;
	int tamanhoGrupos = 0;

	// initialize promisings solutions
	for (int i = 0; i < Tpe; i++)
		Pop[i].promising = 0;

	// save labels defined by LP in groups
	int achei;

    for (int i = 0; i < Tpe; i++)
	{
		achei = 0;
		for (unsigned int j = 0; j < grupos.size(); j++)
		{
			if (Pop[i].label == grupos[j])
				achei = 1;
		}
		if (achei == 0)
		{
			tamanhoGrupos++;
			grupos.push_back(Pop[i].label);
		}
	}

	// find the best solution in the group (with flag = 0)
	for (unsigned int j = 0; j < grupos.size(); j++)
	{
		double menorFO = INFINITY;
		int localMenor = -1;
		int local = -1;
		for (int i = 0; i < Tpe; i++)
		{
			if (Pop[i].label == grupos[j])
			{
				// find the best solution of the group
				if (local == -1)
					local = i;

				// we not apply local search in this solution yet
                if (Pop[i].ofv < menorFO && Pop[i].flag == 0) 
				{
					menorFO = Pop[i].ofv;
					localMenor = i;
				}
			}
		}

		if (localMenor == -1)
			localMenor = local;

		if (Pop[localMenor].label != -1)
			Pop[localMenor].promising = 1;
	}
}

void FreeMemory()
{
    //methods
    Pop.clear();
    PopInter.clear();
    Q.clear();
}

double randomico(double min, double max)
{
	return ((double)(rand()%10000)/10000.0)*(max-min)+min;
}

int irandomico(int min, int max)
{
	return (int)randomico(0,max-min+1.0) + min;
}