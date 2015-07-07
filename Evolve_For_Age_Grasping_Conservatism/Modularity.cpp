// Modularity.cpp : Defines the entry point for the console application.
//

#include "stdio.h"
#include "time.h"

#include "afpo.h"

int main(int argc, char* argv[])
{

	if ( argc < 4 ) {

		printf("missing command line parameters.\n");
		exit(0);
	}

	int randomSeed = atoi(argv[1]);

	int selectForAttractors = atoi(argv[2]);

	int useEvolvedBodies = atoi(argv[3]); 


	printf("Working on ./Modularity %d %d %d\n",randomSeed,selectForAttractors,useEvolvedBodies);

	srand(randomSeed);

	AFPO *afpo = new AFPO(randomSeed,useEvolvedBodies);

	afpo->Evolve(randomSeed,selectForAttractors,useEvolvedBodies);

	delete afpo;
	afpo = NULL;

	printf("\n");

	return 0;
}

