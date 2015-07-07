#include "genome.h"
#include "robot.h"

#ifndef _AFPO_H
#define _AFPO_H

class AFPO {

private:
	GENOME **genomes;
	OBJECT **objects;
	int    numNonDominated;
	int    numWithMaximalGraspingAbility;
	int    nextAvailableID;
	int    currGeneration;
	ROBOT  *robot;

public:
	AFPO(int randomSeed, int useEvolvedBodies);
	~AFPO(void);
	void Evolve(int randomSeed, int selectForAttractors, int useEvolvedBodies);

private:
	void Age_Everyone(void);
	void Count_NumNonDominated(void);
        void Count_NumWithMaximalGraspingAbility(void);
	void Delete_Dominated_Genomes(void);
	void Evaluate_Genomes(int randomSeed, int selectForAttractors, int useEvolvedBodies);
	void Fill_Empty_Slots(void);
	void Find_Pareto_Front(void);
	void Inject_Newcomer(void);
	void Print(void);
        void Sort(void);
        void Sort_By_Pose_Quality(void);
        void Sort_By_Dominated(void);
        void Sort_By_Grasping_Ability(void);
};

#endif
