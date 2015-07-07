#ifndef _AFPO_CPP
#define _AFPO_CPP

#include "genome.h"
#include "afpo.h"

extern int POP_SIZE;
extern int NUM_GENERATIONS;

AFPO::AFPO(void) {

	currGeneration = 0;

	nextAvailableID = 0;

	numNonDominated = 0;

	genomes = new GENOME * [POP_SIZE]; 

	for (int i=0;i<POP_SIZE;i++)

		genomes[i] = new GENOME(nextAvailableID++);

        objects = new OBJECT * [4];

        for (int i=0;i<4;i++)

                objects[i] = new OBJECT(i);

	robot = new ROBOT;
}

AFPO::~AFPO(void) {

	if ( robot ) {

		delete robot;
		robot = NULL;
	}

	if ( genomes ) {

		for (int i=0;i<POP_SIZE;i++) {

			delete genomes[i];
			genomes[i] = NULL;
		}

		delete[] genomes;
		genomes = NULL;
	}

       if ( objects ) {

                for (int i=0;i<4;i++) {
                        delete objects[i];
                        objects[i] = NULL;
                }
                delete[] objects;
                objects = NULL;
        }
}

void AFPO::Evolve(int randomSeed) {

	for (currGeneration=0;currGeneration<NUM_GENERATIONS;currGeneration++) {

		Evaluate_Genomes(randomSeed);

		Find_Pareto_Front();

		Sort_By_Dominated();

		Count_NumNonDominated();

		Sort_By_Fitness();

		// if ( (currGeneration%100)==0 ) {

		if ( currGeneration==(NUM_GENERATIONS-1) ) {

			Print();
			genomes[0]->Evaluate(randomSeed,objects,true,robot);

		}
		Age_Everyone();

		Fill_Empty_Slots();

		Inject_Newcomer();
	}
}

// --------------------- Private methods ------------------

void AFPO::Age_Everyone(void) {

	for (int i=0;i<POP_SIZE;i++)

		genomes[i]->Set_Age( genomes[i]->Get_Age() + 1 );
}

void AFPO::Count_NumNonDominated(void) {

	numNonDominated = 0;

	for (int i=0;i<POP_SIZE;i++)

		numNonDominated = numNonDominated + (genomes[i]->Is_Dominated()==false);
}

void AFPO::Evaluate_Genomes(int randomSeed) {

	for (int i=numNonDominated;i<POP_SIZE;i++)

		genomes[i]->Evaluate(randomSeed, objects,false,robot);
}

void AFPO::Fill_Empty_Slots(void) {

        for (int i=numNonDominated;i<POP_SIZE;i++) {

		int genomeToCopy = rand() % numNonDominated;

		genomes[i]->CopyFrom( genomes[genomeToCopy] , nextAvailableID++ );

		genomes[i]->Mutate();
        }
}

void AFPO::Find_Pareto_Front(void) {

	for (int i=0;i<POP_SIZE;i++)

		genomes[i]->Set_Dominated(false);

	for (int i=0;i<POP_SIZE-1;i++)

		for (int j=i+1;j<POP_SIZE;j++) {

	                if ( genomes[j]->IsDominatedBy(genomes[i]) )

                                genomes[j]->Set_Dominated(true);

                        else if ( genomes[i]->IsDominatedBy(genomes[j]) )

                                genomes[i]->Set_Dominated(true);
		}
}

void AFPO::Inject_Newcomer(void) {

	delete genomes[POP_SIZE-1];

	genomes[POP_SIZE-1] = new GENOME(nextAvailableID++);
}

void AFPO::Print(void) {

	printf("[g: %d] \t",currGeneration);

	printf("[nnd: %d] \t",numNonDominated);

	//for (int i=0;i<numNonDominated;i++)
	//	genomes[i]->Print();

	genomes[0]->Print();
}

void AFPO::Sort_By_Dominated(void) {

	for (int c = 0 ; c < ( POP_SIZE - 1 ); c++) {

		for (int d = 0 ; d < POP_SIZE - c - 1; d++) {

			if ( (genomes[d]->Is_Dominated()==true) && (genomes[d+1]->Is_Dominated()==false) ) {
        
				GENOME *swap       = genomes[d];
        			genomes[d] = genomes[d+1];
        			genomes[d+1] = swap;
				swap = NULL;
      			}
    		}
  	}
}

void AFPO::Sort_By_Fitness(void) {

        for (int c = 0 ; c < ( numNonDominated - 1 ); c++) {

                for (int d = 0 ; d < numNonDominated - c - 1; d++) {

			if ( genomes[d]->Get_Fitness() < genomes[d+1]->Get_Fitness() ) {

                                GENOME *swap       = genomes[d];
                                genomes[d] = genomes[d+1];
                                genomes[d+1] = swap;
                                swap = NULL;
                        }
                }
        }
}

#endif
