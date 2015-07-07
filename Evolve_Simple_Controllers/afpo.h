#include "genome.h"
#include "robot.h"

#ifndef _AFPO_H
#define _AFPO_H

class AFPO {

private:
	GENOME **genomes;
	OBJECT **objects;
	int    numNonDominated;
	int    nextAvailableID;
	int    currGeneration;
	ROBOT  *robot;

public:
	AFPO(void);
	~AFPO(void);
	void Evolve(int randomSeed);

private:
	void Age_Everyone(void);
	void Count_NumNonDominated(void);
	void Delete_Dominated_Genomes(void);
	void Evaluate_Genomes(int randomSeed);
	void Fill_Empty_Slots(void);
	void Find_Pareto_Front(void);
	void Inject_Newcomer(void);
	void Print(void);
	void Sort_By_Dominated(void);
	void Sort_By_Fitness(void);
};

#endif
