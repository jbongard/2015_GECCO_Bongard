#ifndef _GENOME_CPP
#define _GENOME_CPP

#include "stdio.h"
#include "stdlib.h"

#include "genome.h"
#include "matrix.h"
#include "robot.h"

#include <random>

extern int NUM_BODY_PARTS;
extern int NUM_POSSIBLE_ANGLES;

GENOME::GENOME(int newID) {

	ID = newID;

	Create_Random_Angles();

	Create_Random_SegmentLengths();

	fitness = 0.0;

	age = 0;

	dominated = false;
}

GENOME::~GENOME(void) {

	if ( angles ) {
		delete angles;
		angles = NULL;
	}

        if ( segmentLengths ) {
                delete segmentLengths;
                segmentLengths = NULL;
        }
}

void GENOME::CopyFrom(GENOME *parent, int newID) {

        ID = newID;

        angles->CopyFrom( parent->Get_Angles() );

        segmentLengths->CopyFrom( parent->Get_SegmentLengths() );

        fitness = 0.0;

        age = parent->age;

        dominated = false;
}


void GENOME::Evaluate(int randomSeed, OBJECT **objects, int saveToFile, ROBOT *robot) {

	double overallFitness = 0.0;

	for (int i=0;i<4;i++) {

		Evaluate_Against_One_Object(randomSeed,i,objects[i],saveToFile,robot);

		overallFitness = overallFitness + fitness;
	}

	// overallFitness = overallFitness * Hand_Angle_Similarities();

	// overallFitness = overallFitness * Arm_Angle_Similarities();

	fitness = overallFitness / 4.0;

	Compute_Angle_Similarities();

	if ( saveToFile )

		Save(randomSeed);
}

int GENOME::Get_Age(void) {

	return age;
}

double GENOME::Get_Fitness(void) {

	return fitness;
}

int GENOME::Get_ID(void) {

	return ID;
}

int  GENOME::Is_Dominated(void) {

	return dominated;
}

int  GENOME::IsDominatedBy(GENOME *otherGenome) {

        if ( fitness <= otherGenome->fitness ) {

                if ( age >= otherGenome->age ) {

                        if (    (fitness == otherGenome->fitness) &&

                                (age == otherGenome->age) ) {

                                        return ( ID < otherGenome->ID ); // newer genomes dominate older ones.
                        }
                        else {
                                return ( true );
                        }
                }
        }

        return ( false );
}

void GENOME::Mutate(void) {

	int numAngleMutations = rand()%3;

	for (int m=0;m<numAngleMutations;m++)

		Mutate_Angles();

	int numSegmentLengthMutations = rand()%3;

	for (int m=0;m<numSegmentLengthMutations;m++)

		Mutate_SegmentLengths();
}

void GENOME::Print(void) {

	printf("[f: %50.50f] \t",fitness);

	printf("[d: %d] \t",dominated);

        printf("[a: %d] \t",age);

        printf("[as: %5.5f] \t",angleSimilarities);

	printf("\n");
}

void GENOME::Set_Age(int newAge) {

	age = newAge;
}

void GENOME::Set_Dominated(int dominatedState) {

	dominated = dominatedState;
}

// ----------- Private methods ------------------

double GENOME::Angle_For_BodyPart(int angleIndex, int bodyPartIndex) {

	double minAngle, maxAngle;

	if ( bodyPartIndex==3 ) { // Left proximate phalange

		minAngle = -70.0;
		maxAngle = -0.0;
	}

        else if ( bodyPartIndex==4 ) { // Right proximate phalange

                minAngle = +0.0;
                maxAngle = +70.0;
        }

        else if ( bodyPartIndex==5 ) { // Left distal phalange

                minAngle = +0.0;
                maxAngle = +70.0;
        }

        else if ( bodyPartIndex==6 ) { // Right distal phalange

                minAngle = -70.0;
                maxAngle = -0.0;
        }

        else { // All other body parts

                minAngle = -70.0;
                maxAngle = +70.0;
	}

	double angleGap = (maxAngle - minAngle) / (double(NUM_POSSIBLE_ANGLES) - 1.0);

        return minAngle + double(angleIndex)*angleGap;

}

double GENOME::Arm_Angle_Similarities(int firstObject, int secondObject) {

        double equalities = 0.0;
        double counts = 0.0;

        for (int j=0;j<=2;j++) {

                if ( angles->Get(firstObject,j) == angles->Get(secondObject,j) )

                        equalities++;

                counts++;
        }

        return equalities / counts;
}

void GENOME::Compute_Angle_Similarities(void) {

        double armArm02 = Arm_Angle_Similarities(0,2);

        double armArm13 = Arm_Angle_Similarities(1,3);

        double handHand01 = Hand_Angle_Similarities(0,1);

        double handHand23 = Hand_Angle_Similarities(2,3);

        angleSimilarities = (armArm02 + armArm13 + handHand01 + handHand23) / 4.0;
}

void GENOME::Compute_Fitness(OBJECT *object, ROBOT *robot) {

	// Fitness must range within [0,1]

	if ( robot->Claws_Cross() ) {

		fitness = 0.0;
		return;
	}

	double leftFromObj  = robot->Get_LeftFingerTip_Distance_From_Circumference_Of(object);

        double rightFromObj = robot->Get_RightFingerTip_Distance_From_Circumference_Of(object);

	double term1 = 1.0 / (1.0 + leftFromObj);

	double term2 = 1.0 / (1.0 + rightFromObj);

	double diffBetFingertipDistAndCircumference = robot->Get_Diff_Between_Fingertip_Dist_And_Circumference_Of(object);
	
	double term3 = 1.0 / (1.0 + diffBetFingertipDistAndCircumference);

	fitness = term1 * term2 * term3;
}

void GENOME::Create_Random_Angles(void) {

        angles = new MATRIX(4,NUM_BODY_PARTS,0.0);

	for (int i=0;i<4;i++) {

		for (int j=0;j<NUM_BODY_PARTS;j++) {

			int angleIndex = rand() % NUM_POSSIBLE_ANGLES;

			angles->Set(i,j,Angle_For_BodyPart(angleIndex,j));
		}
	}
}

void GENOME::Create_Random_SegmentLengths(void) {

	segmentLengths = new MATRIX(1,NUM_BODY_PARTS,0.0);

	for (int j=0;j<NUM_BODY_PARTS;j++)

		segmentLengths->Set(0,j,Get_RandomSegmentLength());

	Ensure_Valid_Segment_Lengths();
}

void GENOME::Ensure_Valid_Segment_Lengths(void) {

	for (int j=0;j<NUM_BODY_PARTS;j++) {

		if ( segmentLengths->Get(0,j) < 0.5 )

			segmentLengths->Set(0,j,0.5);

		if ( segmentLengths->Get(0,j) > 3.0 )

			segmentLengths->Set(0,j,3.0);
	}
}

void GENOME::Evaluate_Against_One_Object(int randomSeed, int objectIndex, OBJECT *object, int saveToFile, ROBOT *robot) {

        for (int j=0;j<NUM_BODY_PARTS;j++) {

                robot->SetAngle(j,angles->Get(objectIndex,j));
                robot->SetLength(j,segmentLengths->Get(0,j));
        }

        robot->UpdateAngles();
        robot->ComputePositions();

        Compute_Fitness(object,robot);

	if ( saveToFile )

		Save_Robot(randomSeed,objectIndex,robot);

}

MATRIX *GENOME::Get_Angles(void) {

        return angles;
}

double GENOME::Get_RandomSegmentLength(void) {

	double ln = float(rand())/float(RAND_MAX); // ln in [0,1]

	ln = ln * 2.5; // ln in [0,2.5]

	ln = ln + 0.5; // ln in [0.5,3]

	return ln;
}

MATRIX *GENOME::Get_SegmentLengths(void) {

        return segmentLengths;
}

double GENOME::Hand_Angle_Similarities(int firstObject, int secondObject) {

        double equalities = 0.0;
        double counts = 0.0;

        for (int j=3;j<=6;j++) {

                if ( angles->Get(firstObject,j) == angles->Get(secondObject,j) )

                        equalities++;

                counts++;
        }

        return equalities / counts;
}

void GENOME::Mutate_Angles(void) {

	int i = rand() % 4;
        int j = rand() % NUM_BODY_PARTS;

        int newAngleIndex = rand() % NUM_POSSIBLE_ANGLES;

	angles->Set(i,j,Angle_For_BodyPart(newAngleIndex,j));
}

void GENOME::Mutate_SegmentLengths(void) {

	int geneToMutate = rand() % NUM_BODY_PARTS;

	segmentLengths->Set(0,geneToMutate, segmentLengths->Get(0,geneToMutate) + Rand_Gaussian() );

	Ensure_Valid_Segment_Lengths();
}

double GENOME::Rand_Gaussian(void) {

	float x1, x2, w, y1, y2;
 
         do {
                 x1 = 2.0 * float(rand())/float(RAND_MAX) - 1.0;
                 x2 = 2.0 * float(rand())/float(RAND_MAX) - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );

	 return( x1 * w );
}

void GENOME::Save(int randomSeed) {

        char fileName[100];

        sprintf(fileName,"Data/results_%d.txt",randomSeed);

        ofstream *outFile = new ofstream(fileName,ios::app);

        (*outFile) << "g: " << fitness << " ";

        (*outFile) << "as: " << angleSimilarities << " ";

        (*outFile) << "\n";

        outFile->close();

        delete outFile;
        outFile = NULL;
}

void GENOME::Save_Robot(int randomSeed, int objectIndex, ROBOT *robot) {

        char fileName[500];
        sprintf(fileName,"Data/Robot_Matrix_%d_%d.txt",randomSeed,objectIndex);

        robot->Save(fileName);
}

#endif
