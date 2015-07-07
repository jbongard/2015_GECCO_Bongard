/* ---------------------------------------------------
   FILE:     bodyPlan.h
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "matrix.h"
#include "robot.h"
#include "simParams.h"

#ifndef _GENETIC_REGULATORY_NETWORK_H
#define _GENETIC_REGULATORY_NETWORK_H

class GENETIC_REGULATORY_NETWORK {

public:
	int numGenes;
	double fitness;
//	MATRIX *adjMatrix;
	MATRIX *fitnessCases;
	double meanUpdates;
	MATRIX *phenotype;
	int numEdges;
	MATRIX *edgeList;
	int ID;
	ROBOT *robot;

private:
	double D;
	int lastTimeStep;
	int cycleLength;

public:
	GENETIC_REGULATORY_NETWORK(SIM_PARAMS *simParams, int ng);
	GENETIC_REGULATORY_NETWORK(SIM_PARAMS *simParams, GENETIC_REGULATORY_NETWORK *parent);
	~GENETIC_REGULATORY_NETWORK(void);
	void ComputeFitness(int firstEnvironment, int lastEnvironment, int caseIndex);
	void Evaluate(MATRIX *referencePattern, int environmentIndex, int caseIndex);
        double MeanConnectionsWithin(int startGroup1,int endGroup1,int startGroup2,int endGroup2);
	void Mutate(SIM_PARAMS *simParams);
	int  Perfect(void);
	void Perturb(int j);
	void Print(void);
	void Print(int decimalPlaces);
	void PrintPhenotype(void);
	void Reset(void);
	void Robot_PrintResponse(MATRIX *referencePattern, int environmentIndex, int caseIndex);
	void Robot_SaveResponse(MATRIX *referencePattern, int environmentIndex, int caseIndex);
	void Save(char *fileName);

private:
	void ComputeFitnessCase(int fitnessCase, 
							int lastTimeStep, 
							MATRIX *referencePattern,
							int environmentIndex);
	void Connection_Add(void);
	void Connection_Add(int u);
	void Connection_Remove(int u);
	void EdgeList_AddEdge(int j, int i, double val);
	int  EdgeList_In(int i, int j);
	void EdgeList_Print(void);
	void EdgeList_RegulateGenes(int timeStep);
	void EdgeList_Remove(int i, int j);
	int  InAnAttractor(int i1, int i2);
	void Initialize(MATRIX *referencePattern);
	void MutateGene(int u);
	void RegulateGene(int geneIndex, int timeStep);
	void RegulateGenes(int timeStep);
	void Robot_SetPosition(void);
	int  Update(void);
};

#endif
