// Modularity.cpp : Defines the entry point for the console application.
//

#include "stdio.h"
#include "time.h"

#include "afpo.h"

int main(int argc, char* argv[])
{

	int randomSeed = 0;

        if ( argc==1 )

		randomSeed = time(NULL);
        else
		randomSeed = atoi(argv[1]);

        srand(randomSeed);

	AFPO *afpo = new AFPO;

	afpo->Evolve(randomSeed);

	delete afpo;
	afpo = NULL;

	return 0;
}

