/************************************************************************************
									IO Functions
*************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Data.h"
 
void WriteSolutionScreen(TSol s, int n, float timeBest, float timeTotal, char instance[])
{
	printf("\n\n\n Instance: %s \nsol: ", instance);
	for (int i=0; i<n; i++)
		printf("%d ", s.vec[i].sol);

    printf("\nDecoder: %.2lf",s.vec[n].rk);
	printf("\nfo: %.5lf",s.ofv);
	printf("\nTotal time: %.3f",timeTotal);
	printf("\nBest time: %.3f\n\n",timeBest);

	// print Q-Table
	for (unsigned int q=0; q<Q.size(); q++){
		printf("\n");
		if (q == 0)
			printf("rho_e => \t");
		if (q == 1)
			printf("    p => \t");
		if (q == 2)
			printf("  p_m => \t");
		if (q == 3)
			printf("  p_e => \t");

		for (unsigned int j=0; j<Q[q].size(); j++){
			if (q==1)
				printf("[%d %.0lf %.3lf %d %d] \t ", Q[q][j].S, Q[q][j].pVar, Q[q][j].q, Q[q][j].k, Q[q][j].kImp);
			else
				printf("[%d %.2lf %.3lf %d %d] \t ", Q[q][j].S, Q[q][j].pVar, Q[q][j].q, Q[q][j].k, Q[q][j].kImp);
		}
	}
	printf("\n\n");
}

void WriteSolution(TSol s, int n, float timeBest, float timeTotal, char instance[])
{
	FILE *arq;
    arq = fopen("../Results/Solutions.txt","a");

	if (!arq)
	{
		printf("\n\nFile not found Solutions.txt!!!");
		getchar();
		exit(1);
	}

    fprintf(arq,"\n\nInstance: %s", instance);
	fprintf(arq,"\nSol: ");
	for (int i=0; i<n; i++)
		fprintf(arq,"%d ", s.vec[i].sol);
	fprintf(arq,"\nFO: %.5lf", s.ofv);
  	fprintf(arq,"\nBest time: %.3f",timeBest);
	fprintf(arq,"\nTotal time:%.3f \n",timeTotal);

	fclose(arq);
}

void WriteResults(double ofv, double ofvAverage, std::vector <double> ofvs, float timeBest, float timeTotal, char instance[])
{
	FILE *arq;
    arq = fopen("../Results/Results.csv","a");

	if (!arq)
	{
		printf("\n\nFile not found Results.xls!!!");
		getchar();
		exit(1);
	}

	fprintf(arq,"\n%s", instance);
    fprintf(arq,"\t%d", (int)ofvs.size());
    for (unsigned int i=0; i<ofvs.size(); i++){
        fprintf(arq,"\t%lf", ofvs[i]);
	}
	fprintf(arq,"\t%lf", ofv);
	fprintf(arq,"\t%lf", ofvAverage);
	fprintf(arq,"\t%.3f", timeBest);
	fprintf(arq,"\t%.3f", timeTotal);

	fclose(arq);
}