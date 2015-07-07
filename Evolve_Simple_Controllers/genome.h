#include "matrix.h"
#include "object.h"
#include "robot.h"

#ifndef _GENOME_H
#define _GENOME_H

class GENOME {

private:
	MATRIX *angles;
	MATRIX *segmentLengths;
	double fitness;
	int    age;
	int    dominated;
	int    ID;
	double angleSimilarities;

public:
	GENOME(int newID);
	~GENOME(void);
        void    CopyFrom(GENOME *parent, int newID);
	void    Evaluate(int randomSeed, OBJECT **objects, int saveToFile, ROBOT *robot);
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
	double  Arm_Angle_Similarities(int firstObject, int secondObject);
	void    Compute_Angle_Similarities(void);
	void    Compute_Fitness(OBJECT *object, ROBOT *robot);
	void    Create_Random_Angles(void);
	void	Create_Random_SegmentLengths(void);
	void    Ensure_Valid_Segment_Lengths(void);
	void	Evaluate_Against_One_Object(int randomSeed, int objectIndex, OBJECT *object, int saveToFile, ROBOT *robot);
        MATRIX *Get_Angles(void);
	double	Get_RandomSegmentLength(void);
        MATRIX *Get_SegmentLengths(void);
	double  Hand_Angle_Similarities(int firstObject, int secondObject);
	void	Mutate_Angles(void);
	void	Mutate_SegmentLengths(void);
	double  Rand_Gaussian(void);
        void    Save(int randomSeed);
        void    Save_Robot(int randomSeed, int objectIndex, ROBOT *robot);
};

#endif
