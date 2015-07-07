/* ---------------------------------------------------
   FILE:     bodyPlan.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#ifndef _POPULATION_CPP
#define _POPULATION_CPP

#include "pop.h"
#include "simParams.h"

extern int NUM_GENES;

extern SIM_PARAMS *simParams;
extern int GENS_FOR_PHASE_I;
extern int GENS_FOR_PHASE_II;
extern int GENS_FOR_PHASE_III;
extern int GENS_FOR_PHASE_IV;
extern int GENS_FOR_PHASE_V;
extern int GENS_FOR_PHASE_VI;
extern int GENS_FOR_PHASE_VII;
extern int GENS_FOR_PHASE_VIII;

extern int ARM_GENES;
extern int ARM_WRIST_GENES;

extern int PRINT_EVERY;

extern int FITNESS_CASES;

POPULATION::POPULATION(int ps) {

	popSize = ps;

	GRNs = new GENETIC_REGULATORY_NETWORK * [popSize];

	for (int i=0;i<popSize;i++)
		GRNs[i] = new GENETIC_REGULATORY_NETWORK(NUM_GENES);

	CreatePattern_I();
	CreatePattern_II();
	CreatePattern_III();
	CreatePattern_IV();
        CreatePattern_V();
        CreatePattern_VI();
        CreatePattern_VII();
        CreatePattern_VIII();


	fits = new MATRIX(1,GENS_FOR_PHASE_I+GENS_FOR_PHASE_II+GENS_FOR_PHASE_III+GENS_FOR_PHASE_IV+GENS_FOR_PHASE_V+GENS_FOR_PHASE_VI+GENS_FOR_PHASE_VII+GENS_FOR_PHASE_VIII,0.0);

	success = false;
}

POPULATION::~POPULATION(void) {

//	char fileName[100];
//	sprintf(fileName,"Data/runData_%d.dat",simParams->randSeed);
//	fits->WriteToFile(fileName,false);

	delete fits;

	DestroyPattern_I();
	DestroyPattern_II();
	DestroyPattern_III();
	DestroyPattern_IV();
        DestroyPattern_V();
        DestroyPattern_VI();
        DestroyPattern_VII();
        DestroyPattern_VIII();


	for (int i=0;i<popSize;i++) {

		delete GRNs[i];
		GRNs[i] = NULL;
	}

	delete[] GRNs;
	GRNs = NULL;
}

void POPULATION::Evolve(void) {

        Evaluate(referencePattern_I,0,FITNESS_CASES);
        Evaluate(referencePattern_IV,1,FITNESS_CASES);
        Evaluate(referencePattern_II,2,FITNESS_CASES);
        Evaluate(referencePattern_III,3,FITNESS_CASES);
}

void POPULATION::Evolve(int environmentIndex, int caseIndex) {

	Evaluate(referencePattern_I,0,caseIndex);
        Evaluate(referencePattern_IV,1,caseIndex);
        Evaluate(referencePattern_II,2,caseIndex);
        Evaluate(referencePattern_III,3,caseIndex);
}

// --------------------------- Private methods ----------------------------

double POPULATION::BestConnectionsWithin(int startGroup1, int endGroup1, int startGroup2, int endGroup2) {

	return( GRNs[0]->MeanConnectionsWithin(startGroup1,endGroup1,startGroup2,endGroup2) );
}

void POPULATION::ComputeFitnesses(int firstEnvironment, int lastEnvironment, int caseIndex) {

	for (int i=0;i<popSize;i++)

		GRNs[i]->ComputeFitness(firstEnvironment,lastEnvironment,caseIndex);
}

void POPULATION::CreatePattern_I(void) {
        referencePattern_I = new MATRIX(2,NUM_GENES,0.0);
        for (int j=0;j<NUM_GENES;j=j+2) {
                referencePattern_I->Set(0,j  ,-1);
                referencePattern_I->Set(0,j+1,-1);
	}
        referencePattern_I->Print(0);
        printf("\n");
}

void POPULATION::CreatePattern_II(void) {
        referencePattern_II = new MATRIX(2,NUM_GENES,0.0);
        for (int j=0;j<NUM_GENES;j=j+2) {
                referencePattern_II->Set(0,j  ,-1);
                referencePattern_II->Set(0,j+1,+1);
        }
        referencePattern_II->Print(0);
        printf("\n");
}

void POPULATION::CreatePattern_III(void) {
        referencePattern_III = new MATRIX(2,NUM_GENES,0.0);
        for (int j=0;j<NUM_GENES;j=j+2) {
                referencePattern_III->Set(0,j  ,+1);
                referencePattern_III->Set(0,j+1,-1);

        }
        referencePattern_III->Print(0);
        printf("\n");
}

void POPULATION::CreatePattern_IV(void) {
        referencePattern_IV = new MATRIX(2,NUM_GENES,0.0);
        for (int j=0;j<NUM_GENES;j=j+2) {
                referencePattern_IV->Set(0,j  ,+1);
                referencePattern_IV->Set(0,j+1,+1);

        }
        referencePattern_IV->Print(0);
        printf("\n");
}

void POPULATION::CreatePattern_V(void) {
        referencePattern_V = new MATRIX(2,NUM_GENES,0.0);
        for (int j=0;j<NUM_GENES;j=j+3) {
                referencePattern_V->Set(0,j  ,-1);
                referencePattern_V->Set(0,j+1,-1);
                referencePattern_V->Set(0,j+2,-1);
        }
        referencePattern_V->Print(0);
        printf("\n");
}

void POPULATION::CreatePattern_VI(void) {
        referencePattern_VI = new MATRIX(2,NUM_GENES,0.0);
        for (int j=0;j<NUM_GENES;j=j+3) {
                referencePattern_VI->Set(0,j  ,+1);
		referencePattern_VI->Set(0,j+1,-1);
                referencePattern_VI->Set(0,j+2,-1);
        }
        referencePattern_VI->Print(0);
        printf("\n");
}

void POPULATION::CreatePattern_VII(void) {
        referencePattern_VII = new MATRIX(2,NUM_GENES,0.0);
        for (int j=0;j<NUM_GENES;j=j+3) {
                referencePattern_VII->Set(0,j  ,-1);
                referencePattern_VII->Set(0,j+1,-1);
                referencePattern_VII->Set(0,j+2,+1);
        }
        referencePattern_VII->Print(0);
        printf("\n");
}

void POPULATION::CreatePattern_VIII(void) {
        referencePattern_VIII = new MATRIX(2,NUM_GENES,0.0);
        for (int j=0;j<NUM_GENES;j=j+3) {
                referencePattern_VIII->Set(0,j  ,+1);
                referencePattern_VIII->Set(0,j+1,-1);
                referencePattern_VIII->Set(0,j+2,+1);
        }
        referencePattern_VIII->Print(0);
        printf("\n");
}

void POPULATION::DestroyPattern_I(void) {

	delete referencePattern_I;
	referencePattern_I = NULL;
}

void POPULATION::DestroyPattern_II(void) {

	delete referencePattern_II;
	referencePattern_II = NULL;
}

void POPULATION::DestroyPattern_III(void) {

        delete referencePattern_III;
        referencePattern_III = NULL;
}

void POPULATION::DestroyPattern_IV(void) {
        delete referencePattern_IV;
        referencePattern_IV = NULL;
}

void POPULATION::DestroyPattern_V(void) {
        delete referencePattern_V;
        referencePattern_V = NULL;
}
void POPULATION::DestroyPattern_VI(void) {
        delete referencePattern_VI;
        referencePattern_VI = NULL;
}
void POPULATION::DestroyPattern_VII(void) {
        delete referencePattern_VII;
        referencePattern_VII = NULL;
}
void POPULATION::DestroyPattern_VIII(void) {
        delete referencePattern_VIII;
        referencePattern_VIII = NULL;
}


void POPULATION::DistributeCumulativeFitness(void) {

	double cumulativeFitness = GRNs[0]->fitness;

	for (int i=1;i<popSize;i++) {

		cumulativeFitness = cumulativeFitness + GRNs[i]->fitness;
		GRNs[i]->fitness = cumulativeFitness;
	}
}

void POPULATION::Evaluate(MATRIX *refPattern, int environmentIndex, int caseIndex) {

	GRNs[0]->Evaluate(refPattern,environmentIndex,caseIndex);
}

double POPULATION::Fitness_Best(void) {

	return( GRNs[0]->fitness );
}

double POPULATION::Fitness_Mean(void) {

	double mean = 0.0;

	for (int i=0;i<popSize;i++)

		mean = mean + GRNs[i]->fitness;

	mean = mean / double(popSize);

	return( mean );
}

double POPULATION::MeanConnectionsWithin(int startGroup1, int endGroup1, int startGroup2, int endGroup2) {

        double mean = 0.0;
        double num = 0.0;

        for (int i=0;i<popSize;i++) {

                mean = mean + GRNs[i]->MeanConnectionsWithin(startGroup1,endGroup1,startGroup2,endGroup2);
                num = num + 1.0;
        }

        return( mean / num );
}

void POPULATION::Misc(int g) {

	fits->Set(0,g,Fitness_Mean());

	if ( (g%PRINT_EVERY)==0 ) {

		printf("[p: %d] ",simParams->phase);

                printf("[g: %d] [bestF: %4.4f] ",
                                g,
                                Fitness_Best());

                int startMod1    = 0;
                int endMod1      = ARM_GENES-1;
                int startMod2    = endMod1+1;
                int endMod2      = NUM_GENES-1;
                printf("[%d_%d_%d_%d: %3.3f] [%d_%d_%d_%d: %3.3f] [%d_%d_%d_%d: %3.3f] ",
                                                      startMod1,endMod1,startMod1,endMod1,
                                BestConnectionsWithin(startMod1,endMod1,startMod1,endMod1),
                                                      startMod2,endMod2,startMod2,endMod2,
                                BestConnectionsWithin(startMod2,endMod2,startMod2,endMod2),
                                                      startMod1,endMod1,startMod2,endMod2,
                                BestConnectionsWithin(startMod1,endMod1,startMod2,endMod2));

                startMod1    = 0;
                endMod1      = ARM_WRIST_GENES-1;
                startMod2    = endMod1+1;
                endMod2      = NUM_GENES-1;
                printf("[%d_%d_%d_%d: %3.3f] [%d_%d_%d_%d: %3.3f] [%d_%d_%d_%d: %3.3f] \n",
                                                      startMod1,endMod1,startMod1,endMod1,
                                BestConnectionsWithin(startMod1,endMod1,startMod1,endMod1),
                                                      startMod2,endMod2,startMod2,endMod2,
                                BestConnectionsWithin(startMod2,endMod2,startMod2,endMod2),
                                                      startMod1,endMod1,startMod2,endMod2,
                                BestConnectionsWithin(startMod1,endMod1,startMod2,endMod2));
	}

	fitBest = Fitness_Best();
	fitMean = Fitness_Mean();

	NormalizeFitnesses();

	DistributeCumulativeFitness();
}

void POPULATION::Mutate(void) {

	// Elitism
	for (int i=1;i<popSize;i++)
//	for (int i=0;i<popSize;i++)
		GRNs[i]->Mutate();
}

void POPULATION::NormalizeFitnesses(void) {


	double totalFitness = 0.0;

/*
	// fitness by rank
	for (int i=0;i<popSize;i++)
		GRNs[i]->fitness = (popSize-1)-i;
*/
	
	for (int i=0;i<popSize;i++)

		totalFitness = totalFitness + GRNs[i]->fitness;

	for (int i=0;i<popSize;i++) {

		GRNs[i]->fitness = GRNs[i]->fitness / totalFitness;
	}
}

void POPULATION::PrintGRN(int grnID) {

	for (int i=0;i<popSize;i++)

		if ( GRNs[i]->ID==grnID )
			GRNs[i]->Print();
}

void POPULATION::PrintGenomes(void) {

	for (int i=0;i<popSize;i++)

		printf("[%d] [f: %3.3f]\n",i,GRNs[i]->fitness);

	char ch = getchar();
}

void POPULATION::Reset(void) {

	for (int i=0;i<popSize;i++)

		GRNs[i]->Reset();
}

void POPULATION::SaveBestGRN(int fileIndex) {

	char fileName[100];

	sprintf(fileName,"Data/GRN_%d_%d.dat",simParams->randSeed,fileIndex);

	GRNs[0]->Save(fileName);
}

void POPULATION::SaveConnectivities(void) {

        int startMod1;
        int endMod1  ;
        int startMod2;
        int endMod2  ;

	char fileName[100];
	sprintf(fileName,"Data/runData_%d.dat",simParams->randSeed);
	ofstream *outFile = new ofstream(fileName,ios::app);

	(*outFile) << envIndex << " ";
	(*outFile) << fitCase << " ";
	(*outFile) << g       << " ";
	(*outFile) << fitBest << " ";
	(*outFile) << fitMean << "   ";

        startMod1    = 0;
        endMod1      = ARM_GENES-1;
        startMod2    = endMod1+1;
        endMod2      = NUM_GENES-1;

        (*outFile) << BestConnectionsWithin(startMod1,endMod1,startMod1,endMod1) << " ";
        (*outFile) << BestConnectionsWithin(startMod2,endMod2,startMod2,endMod2) << " ";
        (*outFile) << BestConnectionsWithin(startMod1,endMod1,startMod2,endMod2) << "   ";

        startMod1    = 0;
        endMod1      = ARM_WRIST_GENES-1;
        startMod2    = endMod1+1;
        endMod2      = NUM_GENES-1;

        (*outFile) << BestConnectionsWithin(startMod1,endMod1,startMod1,endMod1) << " ";
        (*outFile) << BestConnectionsWithin(startMod2,endMod2,startMod2,endMod2) << " ";
        (*outFile) << BestConnectionsWithin(startMod1,endMod1,startMod2,endMod2) << "   ";

        startMod1    = 0;
        endMod1      = ARM_GENES-1;
        startMod2    = endMod1+1;
        endMod2      = NUM_GENES-1;

	(*outFile) << MeanConnectionsWithin(startMod1,endMod1,startMod1,endMod1) << " ";
	(*outFile) << MeanConnectionsWithin(startMod2,endMod2,startMod2,endMod2) << " ";
	(*outFile) << MeanConnectionsWithin(startMod1,endMod1,startMod2,endMod2) << "   ";

        startMod1    = 0;
        endMod1      = ARM_WRIST_GENES-1;
        startMod2    = endMod1+1;
        endMod2      = NUM_GENES-1;

        (*outFile) << MeanConnectionsWithin(startMod1,endMod1,startMod1,endMod1) << " ";
        (*outFile) << MeanConnectionsWithin(startMod2,endMod2,startMod2,endMod2) << " ";
        (*outFile) << MeanConnectionsWithin(startMod1,endMod1,startMod2,endMod2) << "\n";


	outFile->close();
}

void POPULATION::SaveGRNs(int fileIndex) {

	char fileName[100];

	for (int i=0;i<1;i++) {

		sprintf(fileName,"Data/GRN_%d_%d_%d.dat",simParams->randSeed,fileIndex,i);
		GRNs[i]->Save(fileName);
	}
}

void POPULATION::Select(void) {

	GENETIC_REGULATORY_NETWORK **copiedGRNs;

	copiedGRNs = new GENETIC_REGULATORY_NETWORK * [popSize];

	// Elitism
	copiedGRNs[0] = new GENETIC_REGULATORY_NETWORK(GRNs[0]);
	for (int i=1;i<popSize;i++) {
//	for (int i=0;i<popSize;i++) {

		double selectThreshold = simParams->Rand(0.0,1.0);
		int currentGenomeIndex = 0;

		while (	(currentGenomeIndex<popSize) && 
				(GRNs[currentGenomeIndex]->fitness<selectThreshold) )

				currentGenomeIndex++;

		if ( currentGenomeIndex==popSize )
			currentGenomeIndex = popSize-1;

		copiedGRNs[i] = new GENETIC_REGULATORY_NETWORK(GRNs[currentGenomeIndex]);
	}

	for (int i=0;i<popSize;i++) {
		delete GRNs[i];
		GRNs[i] = copiedGRNs[i];
		copiedGRNs[i] = NULL;
	}

	delete[] copiedGRNs;
	copiedGRNs = NULL;
}

void POPULATION::Sort(void) {

	// Sort by fitness
	int n = popSize;
	while ( n>1 ) {

		int newN = 0;

		for (int i=0;i<n-1;i++) {

			if ( GRNs[i]->fitness <= GRNs[i+1]->fitness ) {

				GENETIC_REGULATORY_NETWORK *tmp = GRNs[i];
				GRNs[i] = GRNs[i+1];
				GRNs[i+1] = tmp;
				tmp = NULL;
				newN = i+1;
			}
		}
		n = newN;
	}

/*
	// For GRNs with same fitness, give more fitness to the one
	// with fewer updates.
	n = popSize;
	while ( n>1 ) {

		int newN = 0;

		for (int i=0;i<n-1;i++) {

			if ( 	(GRNs[i]->fitness == GRNs[i+1]->fitness) &&
				(GRNs[i]->meanUpdates > GRNs[i+1]->meanUpdates) ) {

				GENETIC_REGULATORY_NETWORK *tmp = GRNs[i];
				GRNs[i] = GRNs[i+1];
				GRNs[i+1] = tmp;
				tmp = NULL;
				newN = i+1;
			}
		}
		n = newN;
	}
*/
}

int POPULATION::Success(void) {

	return( success );
}

double POPULATION::Updates_Best(void) {

	return( GRNs[0]->meanUpdates );
}

double POPULATION::Updates_Mean(void) {

	double mean = 0.0;

	for (int i=0;i<popSize;i++)

		mean = mean + GRNs[i]->meanUpdates;

	mean = mean / double(popSize);

	return( mean );
}

#endif
