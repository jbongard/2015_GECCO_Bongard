#include "rbn.h"
#include "matrix.h"
#include "object.h"
#include "robot.h"

#include <vector>
using std::vector;

#ifndef _GENOME_H
#define _GENOME_H

extern int NUM_OBJECTS;
extern int NUM_BODY_PARTS;
extern int NUM_NODES;

class GENOME {

public:
	RANDOM_BOOLEAN_NETWORK    *rbn;
	//MATRIX *angles;
	vector<vector<vector<double> > > angles;
        vector<vector<vector<double> > > nodeValues;
	MATRIX *segmentLengths;
	double fitness;
	int    age;
	int    dominated;
	int    ID;
	double graspingAbility;
	double attractorAttainment;
	double modularity;
	double angleSimilarities;
	double nodeDifferences;
	double poseQuality;

public:
	GENOME(int newID);
	~GENOME(void);
        void    CopyFrom(GENOME *parent, int newID);
	void    Evaluate(int randomSeed, OBJECT **objects, int saveToFile, ROBOT *robot, int selectForAttractors, int useEvolvedBodies);
	int     Get_Age(void);
        double  Get_Fitness(void);
	int     Get_ID(void);
	int     Is_Dominated(void);
        int     IsDominatedBy(GENOME *otherGenome);
	void    Mutate(void);
	void    Print(void);
	void    Set_Age(int newAge);
	void	Set_Dominated(int dominatedState);

private:
	double  Angle_For_BodyPart(int angleIndex, int bodyPartIndex);
	double  Angle_Robustness(int objectIndex, int firstPerturbation, int secondPerturbation);
	double  Arm_Angle_Similarities(int firstObject, int secondObject);
        double  Arm_Angle_Similarities_Across_Perturbations(int objectIndex, int perturbationIndex);
	double  Arm_Node_Differences(int firstObject, int secondObject);
	void    Compute_Angle_Similarities(int numObjectsReached, int numPerturbationsReached);
	void    Compute_Fitness(OBJECT *object, ROBOT *robot);
	void	Compute_Node_Differences(int numObjectsReached, int numPerturbationsReached);
	void	Create_Angles(void);
	void	Create_NodeValues(void);
	void    Create_Random_RBN(void);
	void	Create_Random_SegmentLengths(void);
	void    Ensure_Valid_Segment_Lengths(void);
	void	Evaluate_Against_One_Object_With_Perturbation(int randomSeed, int objectIndex, int perturbationIndex, OBJECT *object, int saveToFile, ROBOT *robot, int useEvolvedBodies);
	double	Get_RandomSegmentLength(void);
        MATRIX *Get_SegmentLengths(void);
	int     Grasping_Failed(void);
	double  Hand_Angle_Similarities(int firstObject, int secondObject);
        double  Hand_Angle_Similarities_Across_Perturbations(int objectIndex, int perturbationIndex);
	double  Hand_Node_Differences(int firstObject, int secondObject);
	void	Mutate_SegmentLengths(void);
	void	Print_NodeValues(int objectIndex, int perturbationIndex);
	double  Rand_Gaussian(void);
	void	Save(int randomSeed);
	void    Save_Robot(int randomSeed, int objectIndex, ROBOT *robot);
};

#endif
