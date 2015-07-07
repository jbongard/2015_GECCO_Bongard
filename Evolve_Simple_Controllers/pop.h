/* ---------------------------------------------------
   FILE:     bodyPlan.h
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#ifndef _POPULATION_H
#define _POPULATION_H

#include "grn.h"

class POPULATION {

public:
	int popSize;
	GENETIC_REGULATORY_NETWORK **GRNs;
	MATRIX *referencePattern_I;
	MATRIX *referencePattern_II;
	MATRIX *referencePattern_III;
	MATRIX *referencePattern_IV;
        MATRIX *referencePattern_V;
        MATRIX *referencePattern_VI;
        MATRIX *referencePattern_VII;
        MATRIX *referencePattern_VIII;

	MATRIX *fits;
	double fitBest;
	double fitMean;
	int success;
	int fitCase;
	int envIndex;

	int g;

public:
	POPULATION(int ps);
	~POPULATION(void);
	void Evolve(void);
	void Evolve(int environmentIndex, int caseIndex);

private:
	double BestConnectionsWithin(int startGroup1, int endGroup1, int startGroup2, int endGroup2);
	void ComputeFitnesses(int firstEnvironment, int lastEnvironment, int caseIndex);
	void CreatePattern_I(void);
	void CreatePattern_II(void);
	void CreatePattern_III(void);
	void CreatePattern_IV(void);
        void CreatePattern_V(void);
        void CreatePattern_VI(void);
        void CreatePattern_VII(void);
        void CreatePattern_VIII(void);
	void DestroyPattern_I(void);
	void DestroyPattern_II(void);
	void DestroyPattern_III(void);
	void DestroyPattern_IV(void);
        void DestroyPattern_V(void);
        void DestroyPattern_VI(void);
        void DestroyPattern_VII(void);
        void DestroyPattern_VIII(void);
	void DistributeCumulativeFitness(void);
	void Evaluate(MATRIX *refPattern, int environmentIndex, int caseIndex);
	double Fitness_Best(void);
	double Fitness_Mean(void);
        double MeanConnectionsWithin(int startGroup1,int endGroup1,int startGroup2,int endGroup2);
	void Misc(int g);
	void Mutate(void);
	void NormalizeFitnesses(void);
	void PrintGenomes(void);
	void PrintGRN(int grnID);
	void Reset(void);
	void SaveBestGRN(int fileIndex);
	void SaveConnectivities(void);
	void SaveGRNs(int fileIndex);
	void Select(void);
	void Sort(void);
	int  Success(void);
	double Updates_Best(void);
	double Updates_Mean(void);
};

#endif
