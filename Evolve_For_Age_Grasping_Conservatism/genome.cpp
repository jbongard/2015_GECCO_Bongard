#ifndef _GENOME_CPP
#define _GENOME_CPP

#include "stdio.h"
#include "stdlib.h"

#include "genome.h"
#include "matrix.h"
#include "robot.h"

#include <random>

extern int NUM_OBJECTS;

extern int NUM_BODY_PARTS;
extern int NUM_POSSIBLE_ANGLES;
extern int NUM_UPDATES;
extern int NUM_NODES;

extern double MIN_ANGLE;
extern double MAX_ANGLE;

extern int NUM_PERTURBATIONS;

GENOME::GENOME(int newID) {

	ID = newID;

	Create_Angles();

	Create_NodeValues();

	Create_Random_RBN(); 

	Create_Random_SegmentLengths();

	fitness = 0.0;

	age = 0;

	dominated = false;

	graspingAbility = 0.0;

	attractorAttainment = 0.0;

	modularity = 0.0;
}

GENOME::~GENOME(void) {

	if ( rbn ) {

		delete rbn;
		rbn = NULL;
	}

        if ( segmentLengths ) {
                delete segmentLengths;
                segmentLengths = NULL;
        }
}

void GENOME::CopyFrom(GENOME *parent, int newID) {

        ID = newID;

	rbn->CopyFrom( parent->rbn );

        segmentLengths->CopyFrom( parent->Get_SegmentLengths() );

        fitness = 0.0;

        age = parent->age;

        dominated = false;

        graspingAbility = 0.0;

        attractorAttainment = 0.0;

	modularity = 0.0;
}

void GENOME::Evaluate(int randomSeed, OBJECT **objects, int saveToFile, ROBOT *robot, int selectForAttractors, int useEvolvedBodies) {

	fitness = 0.0;

	attractorAttainment = 0;

	int failed = false;

        int objectIndex = 0;

	int perturbationIndex = 0;

        while ( (perturbationIndex<NUM_PERTURBATIONS) && (failed==false) ) {

		objectIndex = 0;

        	while ( (objectIndex<NUM_OBJECTS) && (failed==false) ) {

                	Evaluate_Against_One_Object_With_Perturbation(randomSeed,objectIndex,perturbationIndex,objects[objectIndex],saveToFile,robot,useEvolvedBodies);

			fitness = fitness + graspingAbility;

			if ( 	(robot->Claws_Are_Crossed()) ||

				(rbn->In_Point_Attractor()==false) ||

				(graspingAbility < 0.9) )

				failed = true;

			else {
				attractorAttainment++;
				objectIndex++;
			}
		}

		if ( failed==false )

			perturbationIndex++;
	}

	if ( objectIndex == NUM_OBJECTS )

		objectIndex = NUM_OBJECTS - 1;

	if ( perturbationIndex == NUM_PERTURBATIONS )

		perturbationIndex = NUM_PERTURBATIONS - 1;

	Compute_Angle_Similarities(objectIndex,perturbationIndex);

	Compute_Node_Differences(objectIndex,perturbationIndex);

	poseQuality = angleSimilarities * nodeDifferences;

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

                if ( angleSimilarities <= otherGenome->angleSimilarities ) {

                        if ( age >= otherGenome->age ) {

                                if (    (fitness == otherGenome->fitness) &&

                                        (angleSimilarities == otherGenome->angleSimilarities) &&

                                        (age == otherGenome->age) ) {

                                        return ( ID < otherGenome->ID ); // newer genomes dominate older ones.
                                }
                                else {
                                        return ( true );
                                }
                        }
                }
        }

        return ( false );
}

void GENOME::Mutate(void) {

	rbn->Mutate();

	int numSegmentLengthMutations = rand()%3;

	for (int m=0;m<numSegmentLengthMutations;m++)

		Mutate_SegmentLengths();

	modularity = rbn->Get_Modularity();
}

void GENOME::Print(void) {

        printf("[f: %50.50f] \t",fitness); 

	printf("[d: %d] \t",dominated);

        printf("[a: %d] \t",age);

	printf("[m: %5.5f] \t",modularity);

        printf("[g: %5.5f] \t",graspingAbility);

	printf("[aa: %5.5f] \t",attractorAttainment);

        printf("[as: %5.5f] \t",angleSimilarities);

        printf("[nd: %5.5f] \t",nodeDifferences);

        printf("[pq: %5.5f] \t",poseQuality);


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

		minAngle = MIN_ANGLE;
		maxAngle = 0.0;
	}

        else if ( bodyPartIndex==4 ) { // Right proximate phalange

                minAngle = 0.0;
                maxAngle = MAX_ANGLE;
        }

        else if ( bodyPartIndex==5 ) { // Left distal phalange

                minAngle = 0.0;
                maxAngle = MAX_ANGLE;
        }

        else if ( bodyPartIndex==6 ) { // Right distal phalange

                minAngle = MIN_ANGLE;
                maxAngle = 0.0;
        }

        else { // All other body parts

                minAngle = MIN_ANGLE;
                maxAngle = MAX_ANGLE;
	}

	double angleGap = (maxAngle - minAngle) / (double(NUM_POSSIBLE_ANGLES) - 1.0);

	double returnVal = minAngle + double(angleIndex)*angleGap;

	return returnVal;

}

double GENOME::Angle_Robustness(int objectIndex, int firstPerturbation, int secondPerturbation) {

        double equalities = 0.0;
        double counts = 0.0;

        for (int j=0;j<NUM_BODY_PARTS;j++) {

                if ( angles[objectIndex][j][firstPerturbation] == angles[objectIndex][j][secondPerturbation] )

                        equalities++;

                counts++;
        }

        return equalities / counts;
}

double GENOME::Arm_Angle_Similarities(int firstObject, int secondObject) {

        double equalities = 0.0;
	double counts = 0.0;

        for (int j=0;j<3;j++) {

		if ( angles[firstObject][j][0] == angles[secondObject][j][0] )

			equalities++;

		counts++;
	}

	return equalities / counts;
}


double GENOME::Arm_Angle_Similarities_Across_Perturbations(int objectIndex, int perturbationIndex) {

	double distances = 0.0;

	for (int j=0;j<=2;j++)

       		distances = distances + fabs( angles[objectIndex][j][perturbationIndex-1] - angles[objectIndex][j][perturbationIndex] );

	return ( 1.0 / (1.0 + distances) );
}

double GENOME::Arm_Node_Differences(int firstObject, int secondObject) {

        double differences = 0.0;
        double counts = 0.0;

        for (int j=0;j<6;j++) {

                if ( nodeValues[firstObject][j][0] != nodeValues[secondObject][j][0] )

                        differences++;

                counts++;
        }

        return differences / counts;
}

void GENOME::Compute_Angle_Similarities(int numObjectsReached, int numPerturbationsReached) {

	// How robust is the robot to perturbations in its sensors?

	angleSimilarities = 0.0;
	
	double count = 0.0;

	if ( numPerturbationsReached > 0 ) {

		for (int p=1;p<=numPerturbationsReached;p++)
	
			for (int o=0;o<=numObjectsReached;o++) {

				angleSimilarities = angleSimilarities + Angle_Robustness(o,0,p);
				count++;

				// for each object reached, how similar are the 
				// angles in each body part, to the angle of that
				// body part for the zeroth perturbation?
			}
	}

	// How conservative is the movement of the robot in different environments?

	if ( numObjectsReached > 0 ) {

		// In the first and second environments, the object is large,
		// so make sure the hand does the same thing.

		angleSimilarities = angleSimilarities + Hand_Angle_Similarities(0,1);
		count++;
	}

        if ( numObjectsReached > 1 ) {

		// In the first and third environments, the object is on the left,
		// so make sure the arm does the same thing.

                angleSimilarities = angleSimilarities + Arm_Angle_Similarities(0,2);
                count++;
        }

       if ( numObjectsReached > 2 ) {

		// In the second and fourth environments, the object is on the right,
		// so make sure the arm does the same thing.

                angleSimilarities = angleSimilarities + Arm_Angle_Similarities(1,3); 
                count++;

		// In the third and fourth environements, the object is small,
		// so make sure the hand does the same thing.

                angleSimilarities = angleSimilarities + Hand_Angle_Similarities(2,3); 
                count++;
        }

	if ( count == 0.0 )

		angleSimilarities = 1.0;
	else
		angleSimilarities = angleSimilarities / count;
}

void GENOME::Compute_Node_Differences(int numObjectsReached, int numPerturbationsReached) {

	nodeDifferences = 0.0;

	double count = 0.0;

        // How different is the movement of the robot in different environments?

        if ( numObjectsReached > 0 ) {

                // In the first and second environments, the object is at different positions, 
                // so make sure the arm does different things. 

                nodeDifferences = nodeDifferences + Arm_Node_Differences(0,1);
                count++;
        }

        if ( numObjectsReached > 1 ) {

                // In the first and third environments, the object is different sizes, 
                // so make sure the hand does different things.

                nodeDifferences = nodeDifferences + Hand_Node_Differences(0,2);
                count++;
        }

       if ( numObjectsReached > 2 ) {

                // In the second and fourth environments, the object is different sizes, 
                // so make sure the hand does different things.

                nodeDifferences = nodeDifferences + Hand_Node_Differences(1,3);
                count++;

                // In the third and fourth environements, the object is at different positions, 
                // so make sure the arm does different things.

                nodeDifferences = nodeDifferences + Arm_Node_Differences(2,3);
                count++;
        }

        if ( count == 0.0 )

                nodeDifferences = 1.0;
        else
                nodeDifferences = nodeDifferences / count;

}

void GENOME::Compute_Fitness(OBJECT *object, ROBOT *robot) {

	// Fitness must range within [0,1]

	double leftFromObj  = robot->Get_LeftFingerTip_Distance_From_Circumference_Of(object);

        double rightFromObj = robot->Get_RightFingerTip_Distance_From_Circumference_Of(object);

	double term1 = 1.0 / (1.0 + leftFromObj);

	double term2 = 1.0 / (1.0 + rightFromObj);

	double diffBetFingertipDistAndCircumference = robot->Get_Diff_Between_Fingertip_Dist_And_Circumference_Of(object);
	
	double term3 = 1.0 / (1.0 + diffBetFingertipDistAndCircumference);

	graspingAbility = term1 * term2 * term3;
}

void GENOME::Create_Angles(void) {

	angles.resize(NUM_OBJECTS);

	for (int i=0; i < NUM_OBJECTS ; i++ ) {

		angles[i].resize(NUM_BODY_PARTS);

		for (int j=0 ; j < NUM_BODY_PARTS ; j++ )

			angles[i][j].resize(NUM_PERTURBATIONS);
	}
}

void GENOME::Create_NodeValues(void) {

        nodeValues.resize(NUM_OBJECTS);

        for (int i=0; i < NUM_OBJECTS ; i++ ) {

                nodeValues[i].resize(NUM_NODES);

                for (int j=0 ; j < NUM_NODES ; j++ )

                        nodeValues[i][j].resize(NUM_PERTURBATIONS);
        }
}

void GENOME::Create_Random_RBN(void) {

	rbn = new RANDOM_BOOLEAN_NETWORK();
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

void GENOME::Evaluate_Against_One_Object_With_Perturbation(int randomSeed, int objectIndex, int perturbationIndex, OBJECT *object, int saveToFile, ROBOT *robot, int useEvolvedBodies) {

	rbn->Set_Initial_Conditions(objectIndex);

	if ( perturbationIndex == 0 ) {

		// Don't perturb the RBN.
	}

	else if ( perturbationIndex <= NUM_NODES ) {

		// Perturb just one node.

		rbn->Perturb_Node(perturbationIndex-1);
	}

	else {
		perturbationIndex = perturbationIndex - NUM_NODES;

		int p1 = 0;
		int p2 = 0;
		int k = 0; 

		while ( (p1<(NUM_NODES-1)) && (k<perturbationIndex) ) {

			p2 = p1+1;

			while ( (p2<NUM_NODES) && (k<perturbationIndex) ) {

				k++;

				if ( k<perturbationIndex )
					p2++;
			}

			if ( k<perturbationIndex )
				p1++;
		}

		rbn->Perturb_Node(p1);
               	rbn->Perturb_Node(p2);
	}

	rbn->Update();

	for (int j=0;j<NUM_NODES;j++)

		nodeValues[objectIndex][j][perturbationIndex] = rbn->Get_NodeValue(j);

//	fitness = 1.0 / (1.0 + rbn->Distance_From_Target_Pattern(objectIndex));

        for (int j=0;j<NUM_BODY_PARTS;j++) {

		int angleIndex = rbn->Get_Angle_Index(j);

		double angle = Angle_For_BodyPart(angleIndex,j);

		angles[objectIndex][j][perturbationIndex] = angle;

                robot->SetAngle(j, angle );

		if ( useEvolvedBodies )

                	robot->SetLength(j,segmentLengths->Get(0,j));
        }

        robot->UpdateAngles();
        robot->ComputePositions();

        Compute_Fitness(object,robot);

	if ( saveToFile && (perturbationIndex==0) ) {

		rbn->Print();

		Save_Robot(randomSeed,objectIndex,robot);
	}
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

int    GENOME::Grasping_Failed(void) {

	return( fitness < 0.9 );
}

double GENOME::Hand_Angle_Similarities(int firstObject, int secondObject) {

        double equalities = 0.0;
        double counts = 0.0;

        for (int j=3;j<NUM_BODY_PARTS;j++) {

                if ( angles[firstObject][j][0] == angles[secondObject][j][0] )

                        equalities++;

                counts++;
        }

        return equalities / counts;
}

double GENOME::Hand_Angle_Similarities_Across_Perturbations(int objectIndex, int perturbationIndex) {

        double distances = 0.0;

        for (int j=3;j<=6;j++)

                distances = distances + fabs( angles[objectIndex][j][perturbationIndex-1] - angles[objectIndex][j][perturbationIndex] );

        return( 1.0 / (1.0 + distances) );
}

double GENOME::Hand_Node_Differences(int firstObject, int secondObject) {

        double differences = 0.0;
        double counts = 0.0; 

        for (int j=6;j<NUM_NODES;j++) {

                if ( nodeValues[firstObject][j][0] != nodeValues[secondObject][j][0] )

                        differences++;

                counts++;
        }
                        
        return differences / counts;
}

void GENOME::Mutate_SegmentLengths(void) {

	int geneToMutate = rand() % NUM_BODY_PARTS;

	segmentLengths->Set(0,geneToMutate, segmentLengths->Get(0,geneToMutate) + Rand_Gaussian() );

	Ensure_Valid_Segment_Lengths();
}

void   GENOME::Print_NodeValues(int objectIndex, int perturbationIndex) {

	for (int o=0;o<=objectIndex;o++)

		for (int p=0;p<=perturbationIndex;p++) {

			for (int n=0;n<NUM_NODES;n++) {

                                if ( n==(2*3) )
                                        printf(" ");

				if ( nodeValues[o][n][p] == 1.0 )

					printf("P");
				else
					printf("N");
			}
		}

	char ch = getchar();
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

        (*outFile) << "m: " << rbn->Get_Modularity() << " ";

        (*outFile) << "g: " << graspingAbility << " ";

        (*outFile) << "aa: " << attractorAttainment << " ";

	(*outFile) << "as: " << angleSimilarities << " ";

        (*outFile) << "nd: " << nodeDifferences << " ";

        (*outFile) << "pq: " << poseQuality << " ";

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
