/* ---------------------------------------------------
   FILE:     bodyPlan.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#ifndef _GENETIC_REGULATORY_NETWORK_CPP
#define _GENETIC_REGULATORY_NETWORK_CPP

#include "grn.h"
#include "math.h"
#include "simParams.h"

extern int		NUM_GENES;
extern int		START_NUM_CONNECTIONS;
extern double	PERTURBATION_PROBABILITY;
extern double	MUTATION_PROBABILITY;
extern int		MAX_CYCLE_LENGTH;
extern int		MAX_GRN_UPDATES;
extern int		FITNESS_CASES;
extern int		NUM_PHASES;

GENETIC_REGULATORY_NETWORK::GENETIC_REGULATORY_NETWORK(SIM_PARAMS *simParams, int ng) {

	ID = simParams->currentAvailableID++;

	numGenes = ng;

//	adjMatrix = new MATRIX(numGenes,numGenes,0.0);

	numEdges = 0;
	edgeList = new MATRIX(numGenes*numGenes,3,-1.0);

	for (int i=0;i<START_NUM_CONNECTIONS;i++)

		Connection_Add();
	
	fitnessCases = new MATRIX(NUM_PHASES,FITNESS_CASES,0.0);

	fitness = 0.0;

	meanUpdates = 0.0;

	robot = new ROBOT();
}

GENETIC_REGULATORY_NETWORK::GENETIC_REGULATORY_NETWORK(SIM_PARAMS *simParams, GENETIC_REGULATORY_NETWORK *parent) {

        ID = simParams->currentAvailableID++;

	numGenes = parent->numGenes;

//	adjMatrix = new MATRIX(parent->adjMatrix);

	numEdges = parent->numEdges;
	edgeList = new MATRIX(parent->edgeList);

	fitnessCases = new MATRIX(NUM_PHASES,FITNESS_CASES,0.0);

	fitness = 0.0;

	meanUpdates = 0.0;

	robot = new ROBOT();
}

GENETIC_REGULATORY_NETWORK::~GENETIC_REGULATORY_NETWORK(void) {

	delete robot;
	robot = NULL;

	delete fitnessCases;
	fitnessCases = NULL;

//	delete adjMatrix;
//	adjMatrix = NULL;

	delete edgeList;
	edgeList = NULL;
}

void GENETIC_REGULATORY_NETWORK::ComputeFitness(int firstEnvironment, int lastEnvironment, int caseIndex) {

	double g = 0.0;
	double num = 0.0;

	if ( caseIndex>1 ) {

		for (int i=firstEnvironment;i<=3;i++) {

			for (int j=1;j<caseIndex;j++) {

                        	if ( i>0 ) {
                                	if ( fitnessCases->Get(i-1,j-1)>0.999 )
                                        	g = g + fitnessCases->Get(i,j-1);
                        	}
                        	else
                                	g = g + fitnessCases->Get(i,j-1);
                        	num++;
			}
		}
	}
                for (int i=firstEnvironment;i<=lastEnvironment;i++) {

			int j = caseIndex;

                                if ( i>0 ) {
                                        if ( fitnessCases->Get(i-1,j-1)>0.999 ) {
                                                g = g + fitnessCases->Get(i,j-1);
					}
                               }
                               else {
                                        g = g + fitnessCases->Get(i,j-1);
				}
                                num++;

                }

/*	
	for (int i=firstEnvironment;i<=lastEnvironment;i++) {

		int numCases;

		if ( caseIndex==1 )
			numCases = caseIndex;
		else {
			if ( i<=lastEnvironment )
				numCases = caseIndex;
			else
				numCases = caseIndex-1;
		}

		for (int j=0;j<numCases;j++) {

			if ( i>0 ) {
				if ( fitnessCases->Get(i-1,j)>0.999 )
					g = g + fitnessCases->Get(i,j);
			}
			else
				g = g + fitnessCases->Get(i,j);

			num++;
		}
	}

*/
	g = g / num;

	meanUpdates = meanUpdates / num;

//	fitness = 1.0 - exp(-3.0*g);
	fitness = g;
}

void GENETIC_REGULATORY_NETWORK::Evaluate(MATRIX *referencePattern, int environmentIndex, int caseIndex) {

	int f=0;
		for (int j=0;j<caseIndex;j++) {

                	Initialize(referencePattern);

			if ( j>0 )
				Perturb(j-1);

	                lastTimeStep = 0;
	                Robot_SetPosition();

			printf("%d %d\n",environmentIndex,j);
			robot->Print();

        		char fileName[100];
        		sprintf(fileName,"Data/robot_%d_%d.dat",environmentIndex,j);

                	robot->Save(fileName);
		}
}
double GENETIC_REGULATORY_NETWORK::MeanConnectionsWithin(int startGroup1, int endGroup1,
                                                         int startGroup2, int endGroup2) {

        double mean = 0.0;
        double num = 0.0;

	for (int i=0;i<numEdges;i++) {

//		int source = edgeList->Get(i,0);
		int source = int(edgeList->vals[i*edgeList->width]);

//		int target = edgeList->Get(i,1);
		int target = int(edgeList->vals[i*edgeList->width+1]);

		if (	(source>=startGroup1) &&
			(source<=endGroup1) &&
			(target>=startGroup2) &&
			(target<=endGroup2) )
			mean++;

		if (	(source>=startGroup2) &&
			(source<=endGroup2) &&
			(target>=startGroup1) &&
			(target<=endGroup1) )
			mean++;
	}

	num = 	2 * 
		( endGroup1 - startGroup1 + 1 ) * 
		( endGroup2 - startGroup2 + 1 );

        return(mean/num);
}

void GENETIC_REGULATORY_NETWORK::Mutate(SIM_PARAMS *simParams) {

	int mutated = false;

	for (int u=0;u<numGenes;u++)

		if ( simParams->Rand(0.0,1.0) < MUTATION_PROBABILITY ) {

			MutateGene(u);
			mutated = true;
		}
}

void GENETIC_REGULATORY_NETWORK::Perturb(int j) {

	if ( simParams->geneValueHistory->Get(0,j) == 1.0 )
		simParams->geneValueHistory->Set(0,j,-1.0);
	else
		simParams->geneValueHistory->Set(0,j,+1.0);
}

int  GENETIC_REGULATORY_NETWORK::Perfect(void) {

	return( fitness == (1.0 - exp(-3.0*1.0)) );
}

void GENETIC_REGULATORY_NETWORK::Print(void) {

	printf("[ID: %d] \t",ID);
	printf("[F: %4.4f] \t",fitness);
	printf("\n");

	Print(0);
}

void GENETIC_REGULATORY_NETWORK::Print(int decimalPlaces) {

//	adjMatrix->Print(decimalPlaces);
	edgeList->Print(decimalPlaces);
}

void GENETIC_REGULATORY_NETWORK::PrintPhenotype(void) {

	phenotype->Print(1);
}

void GENETIC_REGULATORY_NETWORK::Reset(void) {

	fitness = 0.0;
	meanUpdates = 0.0;
        D = 0.0;
        lastTimeStep = 0;
        cycleLength = 0;
}

void GENETIC_REGULATORY_NETWORK::Robot_PrintResponse(MATRIX *referencePattern, int environmentIndex, int caseIndex) {

	Initialize(referencePattern);

        if ( caseIndex>1 )
                Perturb(caseIndex-2);

	lastTimeStep = Update();
	if ( lastTimeStep == MAX_GRN_UPDATES )
		lastTimeStep = MAX_GRN_UPDATES-1;
	Robot_SetPosition();

	robot->Print();
//	simParams->geneValueHistory->Print(0,lastTimeStep,0,NUM_GENES-1);

}

void GENETIC_REGULATORY_NETWORK::Robot_SaveResponse(MATRIX *referencePattern, int environmentIndex, int caseIndex) {

        char fileName[100];
        sprintf(fileName,"Data/robot_%d_%d.dat",simParams->randSeed,environmentIndex);

        Initialize(referencePattern);

//        lastTimeStep = Update();
        int timeStep = 1;
        int inAnAttractor = false;
	lastTimeStep=0;
        while ( (!inAnAttractor) && (timeStep<MAX_GRN_UPDATES) ) {
		Robot_SetPosition();
		robot->Save(fileName);
                EdgeList_RegulateGenes(timeStep);
                cycleLength = 1;
                while ( (!inAnAttractor) && (cycleLength<MAX_CYCLE_LENGTH) ) {
                        if ( timeStep>=cycleLength )
                                inAnAttractor = InAnAttractor(timeStep-cycleLength,timeStep);
                        if ( !inAnAttractor )
                                cycleLength++;
                }
                if ( !inAnAttractor ) {
                        timeStep++;
			lastTimeStep++;
		}
	}
}


void GENETIC_REGULATORY_NETWORK::Save(char *fileName) {

	ofstream *outFile = new ofstream(fileName);

//	adjMatrix->Write(outFile);
	edgeList->Write(outFile);

	outFile->close();
	delete outFile;
	outFile = NULL;
}

// ----------------- Private methods -------------------

void GENETIC_REGULATORY_NETWORK::ComputeFitnessCase(	int fitnessCase, 
							int lastTimeStep, 
							MATRIX *referencePattern,
							int environmentIndex) {
/*
	D = 0.0;

	for (int g=0;g<NUM_GENES;g++) {

		//double desiredGeneValue = referencePattern->Get(1,g);
		double desiredGeneValue = referencePattern->vals[1*referencePattern->width+g];

		//double actualGeneValue  = simParams->geneValueHistory->Get(lastTimeStep,g);
		double actualGeneValue = 
		simParams->geneValueHistory->vals[lastTimeStep*simParams->geneValueHistory->width+g];

		double diff = fabs( desiredGeneValue - actualGeneValue );

		D = D + diff;
	}

	double DMax = 2.0*NUM_GENES;

	if ( lastTimeStep == (MAX_GRN_UPDATES-1) )
		D = DMax;

	if ( cycleLength > 1 )
		D = DMax;

//	fitnessCases->Set(environmentIndex,fitnessCase,pow(1.0-D/DMax,5.0));
	fitnessCases->vals[environmentIndex*fitnessCases->width+fitnessCase] = pow(1.0-D/DMax,5.0);
*/

	double circleX;
	double circleY;
	double circleRadius;

	if ( environmentIndex==0 ) {
		circleX      = -3.268;
		circleY      = -0.853;
		circleRadius = 1.207;
	}
	else if ( environmentIndex==2 ) {
                circleX      = +3.268;
                circleY      = -0.853;
                circleRadius = 1.207;
	}
	else if ( environmentIndex==3 ) {
                circleX      = -2.793;
                circleY      = -0.379;
                circleRadius = 0.536;
	}
	else if ( environmentIndex==1 ) {
                circleX      = +2.793;
                circleY      = -0.379;
                circleRadius = 0.536;

	}
/*
        else if ( environmentIndex==4 ) {
                circleX = -3.5;
                circleRadius = 0.75;
        }
        else if ( environmentIndex==5 ) {
                circleX = +3.5;
                circleRadius = 0.75;
        }
        else if ( environmentIndex==6 ) {
                circleX = -3.0;
                circleRadius = 0.943; 
        }
        else if ( environmentIndex==7 ) {
                circleX = +3.0;
                circleRadius = 0.943;
        }
*/
	double fit = 0.0;

/*
        // Minimize distance from wrist to circle center 
                double fingerX = robot->robot->Get(2,3);
                double fingerY = robot->robot->Get(2,4);

                double diffX = (fingerX-circleX)*(fingerX-circleX);
                double diffY = (fingerY-circleY)*(fingerY-circleY);
                double dist = sqrt( diffX + diffY );
                double distFromCenter = fabs( 0 - dist );

                fit = fit - distFromCenter;
*/

	// Minimize distance from fingertips to circle rim
	for (int i=5;i<=6;i++) {

		double fingerX = robot->robot->Get(i,3);
		double fingerY = robot->robot->Get(i,4);

		double diffX = (fingerX-circleX)*(fingerX-circleX);
		double diffY = (fingerY-circleY)*(fingerY-circleY);
		double dist = sqrt( diffX + diffY );
		double distFromRim = fabs( circleRadius - dist );

//		if ( distFromRim > 0.01 )
//			distFromRim = +8.0;

		fit = fit + distFromRim;
	}
	fit = fit/2.0;
/*
        double minDist = 1000.0;
        // Maximize the minimum distance between any two proximal phalanges 
        for (int i=3;i<=3;i++) {

                double iX = robot->robot->Get(i,3);
                double iY = robot->robot->Get(i,4);

                for (int j=i+1;j<=4;j++) {
                        double jX = robot->robot->Get(j,3);
                        double jY = robot->robot->Get(j,4);
                        double diffX = (iX-jX)*(iX-jX);
                        double diffY = (iY-jY)*(iY-jY);
                        double dist = sqrt( diffX + diffY );
                        if ( dist < minDist )
                                minDist = dist;
                }
        }
        fit = fit + minDist;
*/

/*
	double minDist = 1000.0;
	// Maximize the minimum distance between any two finger tips.
	for (int i=5;i<=5;i++) {	

		double iX = robot->robot->Get(i,3);
		double iY = robot->robot->Get(i,4);

		for (int j=i+1;j<=6;j++) {
			double jX = robot->robot->Get(j,3);
			double jY = robot->robot->Get(j,4);
			double diffX = (iX-jX)*(iX-jX);
			double diffY = (iY-jY)*(iY-jY);
			double dist = sqrt( diffX + diffY );
			if ( dist < minDist )
				minDist = dist;
		}
	}
	fit = fit + minDist;
*/

	// Fingers crossed?
	if ( robot->robot->Get(4,2) <= robot->robot->Get(3,2) )
		fit = +8.0;
        else if ( lastTimeStep == (MAX_GRN_UPDATES-1) )
		fit = +8.0;
        else if ( cycleLength > 1 )
                fit = +8.0; 

	fit = pow(1.0-fit/8.0,1.0);

/*
	if ( fit > -0.0012 )
		fit = 1.0;
	else
		fit = 0.0;
*/

	fitnessCases->vals[environmentIndex*fitnessCases->width+fitnessCase] = fit;	
}

void GENETIC_REGULATORY_NETWORK::Connection_Add(void) {

	int i = simParams->RandInt(0,numGenes-1);
	int j = simParams->RandInt(0,numGenes-1);

//	while ( adjMatrix->Get(j,i) != 0.0 ) {
	while ( EdgeList_In(j,i) ) {

		i = simParams->RandInt(0,numGenes-1);
		j = simParams->RandInt(0,numGenes-1);
	}

	if ( simParams->FlipCoin() ) {

//		adjMatrix->Set(j,i,+1.0);
		EdgeList_AddEdge(j,i,+1.0);
	}

	else {
//		adjMatrix->Set(j,i,-1.0);
		EdgeList_AddEdge(j,i,-1.0);
	}
}

void GENETIC_REGULATORY_NETWORK::Connection_Add(int u) {

	int j = simParams->RandInt(0,numGenes-1);

//	while ( adjMatrix->Get(j,u) != 0.0 ) {
	while ( EdgeList_In(j,u) ) {

		j = simParams->RandInt(0,numGenes-1);
	}

	if ( simParams->FlipCoin() ) {

//		adjMatrix->Set(j,u,+1.0);
		EdgeList_AddEdge(j,u,+1.0);
	}
	else {
//		adjMatrix->Set(j,u,-1.0);
		EdgeList_AddEdge(j,u,-1.0);
	}
}

void GENETIC_REGULATORY_NETWORK::Connection_Remove(int u) {

	int j = simParams->RandInt(0,numGenes-1);

//	while ( adjMatrix->Get(j,u) == 0.0 ) {
	while ( !EdgeList_In(j,u) ) {

		j = simParams->RandInt(0,numGenes-1);
	}

//	adjMatrix->Set(j,u,0.0);
	EdgeList_Remove(j,u);
}

void GENETIC_REGULATORY_NETWORK::EdgeList_AddEdge(int j, int i, double val) {

	edgeList->Set(numEdges,0,j);
	edgeList->Set(numEdges,1,i);
	edgeList->Set(numEdges,2,val);

	numEdges++;
}

int  GENETIC_REGULATORY_NETWORK::EdgeList_In(int i, int j) {

        int found = false;
        int currEdge = 0;

        while ( (!found) && (currEdge<numEdges) ) {

		if (	(edgeList->vals[currEdge*edgeList->width]==i) &&
			(edgeList->vals[currEdge*edgeList->width+1]==j) )
                        found = true;
                else
                        currEdge++;
        }

        return( found );
}

void GENETIC_REGULATORY_NETWORK::EdgeList_Print(void) {

        for (int i=0;i<numEdges;i++) {

                printf("%1.0f \t",edgeList->Get(i,0));
                printf("%1.0f \t",edgeList->Get(i,1));
                printf("%1.0f \n",edgeList->Get(i,2));

        }
}

void GENETIC_REGULATORY_NETWORK::EdgeList_RegulateGenes(int timeStep) {

        for (int j=0;j<NUM_GENES;j++)
                simParams->geneValueHistory->vals[timeStep*simParams->geneValueHistory->width + j] = 0.0;

        for (int i=0;i<numEdges;i++) {

                int source    = int(edgeList->vals[i*edgeList->width + 0]);
                int target    = int(edgeList->vals[i*edgeList->width + 1]);

                double influence =  edgeList->vals[i*edgeList->width + 2];

                double geneValue = simParams->geneValueHistory->vals[
                                        (timeStep-1)*simParams->geneValueHistory->width + source];

                double amountOfRegulation = geneValue * influence;

                simParams->geneValueHistory->vals[timeStep*simParams->geneValueHistory->width + target] =
                simParams->geneValueHistory->vals[timeStep*simParams->geneValueHistory->width + target] +
                amountOfRegulation;
        }

        for (int j=0;j<NUM_GENES;j++)

                if ( simParams->geneValueHistory->vals[timeStep*simParams->geneValueHistory->width + j] > 0 )

                        simParams->geneValueHistory->vals[timeStep*simParams->geneValueHistory->width + j] = +1.0;
                else
                        simParams->geneValueHistory->vals[timeStep*simParams->geneValueHistory->width + j] = -1.0;
}

void GENETIC_REGULATORY_NETWORK::EdgeList_Remove(int i, int j) {

	int found = false;
	int edgeIndex = 0;

	while ( !found ) {
		if ( 	(edgeList->Get(edgeIndex,0)==i) &&
			(edgeList->Get(edgeIndex,1)==j) )
			found = true;
		else
			edgeIndex++;
	}

        for (int i=edgeIndex;i<(numEdges-1);i++) {

                edgeList->CopyRow(i+1,i);
        }
        numEdges--;
}

int  GENETIC_REGULATORY_NETWORK::InAnAttractor(int i1, int i2) {

	int inAnAttractor = true;
	int geneIndex = 0;

	while ( (inAnAttractor) && (geneIndex<NUM_GENES) ) {

		if (	(simParams->geneValueHistory->vals[i1*simParams->geneValueHistory->width + geneIndex]) !=
				(simParams->geneValueHistory->vals[i2*simParams->geneValueHistory->width + geneIndex]) )

			 inAnAttractor = false;
		else
			geneIndex++;
	}

	return( inAnAttractor );
}

void GENETIC_REGULATORY_NETWORK::Initialize(MATRIX *referencePattern) {

	simParams->geneValueHistory->CopyRow(0,referencePattern,0);
}

void GENETIC_REGULATORY_NETWORK::MutateGene(int u) {

	double N = numGenes;

	double ru = 0.0;

	for (int j=0;j<numGenes;j++)
//		if ( adjMatrix->Get(j,u) != 0.0 )
		if ( EdgeList_In(j,u) )
			ru++;

	double pu = (4*ru) / (4*ru + (N-ru) ); // Specialization paper, Eq. (5)

	if ( simParams->Rand(0.0,1.0) < pu )
		Connection_Remove(u);

	if ( simParams->Rand(0.0,1.0) < (1.0-pu) )
		Connection_Add(u);
}

void GENETIC_REGULATORY_NETWORK::RegulateGene(int geneIndex, int timeStep) {

/*
	simParams->geneValueHistory->Set(timeStep,geneIndex,0.0);

	for (int i=0;i<NUM_GENES;i++) {

		double geneValue = 
		simParams->geneValueHistory->vals[(timeStep-1)*simParams->geneValueHistory->width + i];

		double regulation = adjMatrix->vals[i*adjMatrix->width + geneIndex];

		double amountOfRegulation = geneValue * regulation;

		simParams->geneValueHistory->vals[    
		timeStep*simParams->geneValueHistory->width + geneIndex] = 
		simParams->geneValueHistory->vals[timeStep*simParams->geneValueHistory->width + geneIndex] + 
		amountOfRegulation;

	}

	if ( simParams->geneValueHistory->vals[timeStep*simParams->geneValueHistory->width + geneIndex] > 0 )

		simParams->geneValueHistory->vals[ timeStep*simParams->geneValueHistory->width + geneIndex] 
		= +1.0;

	else
		simParams->geneValueHistory->vals[ timeStep*simParams->geneValueHistory->width + geneIndex] 
		= -1.0;
*/
}

void GENETIC_REGULATORY_NETWORK::RegulateGenes(int timeStep) {

	for (int j=0;j<NUM_GENES;j++)

		RegulateGene(j,timeStep);
}

void   GENETIC_REGULATORY_NETWORK::Robot_SetPosition(void) {

	int currentBodyPart = 0;
	double leftGC, rightGC;

	for (int j=0;j<NUM_GENES;j=j+2) {

		int b1 = int(simParams->geneValueHistory->Get(lastTimeStep,j));
		int b2 = int(simParams->geneValueHistory->Get(lastTimeStep,j+1));

//		double gc = double(simParams->GrayCode(b1,b2,b3));
		double gc = double(simParams->BinaryToInteger(b1,b2));

/*
                if ( (gc==3) || (gc==4) )
                        gc=3.5;
*/
//		if ( (gc==1) || (gc==2) )
//			gc=1.5;
/*
		if ( currentBodyPart==3 ) { // left proximate phalange 

			if ( gc >= 7-1 )
				gc=7-2;

			leftGC = gc;
		}
		if ( currentBodyPart==4 ) { // right proximate phalange 
			if ( gc <= leftGC )
				if ( leftGC==3.5 )
					gc = 5;
				else
					gc = leftGC + 1; 
			rightGC = gc;
		}
*/

/*
if ( currentBodyPart== 0 )
	gc = 3;
if ( currentBodyPart== 1 )
        gc = 3;
if ( currentBodyPart== 2 )
        gc = 3;
if ( currentBodyPart== 3 )
        gc = 0;
if ( currentBodyPart== 4 )
        gc = 3;
if ( currentBodyPart== 5 )
        gc = 0;
if ( currentBodyPart== 6 )
        gc = 3;
*/

//		double angle = simParams->Scale(gc,0,7,-45,45);
		double angle = simParams->Scale(gc,0,3,-45,45);

		robot->SetAngle(currentBodyPart,angle);

		currentBodyPart++;
	}

	robot->UpdateAngles();
	robot->ComputePositions();
}

int  GENETIC_REGULATORY_NETWORK::Update(void) {

	int timeStep = 1;
	int inAnAttractor = false;

	while ( (!inAnAttractor) && (timeStep<MAX_GRN_UPDATES) ) {
		
//		RegulateGenes(timeStep);
		EdgeList_RegulateGenes(timeStep);

		cycleLength = 1;
		
		while ( (!inAnAttractor) && (cycleLength<MAX_CYCLE_LENGTH) ) {

			if ( timeStep>=cycleLength )
				inAnAttractor = InAnAttractor(timeStep-cycleLength,timeStep);

			if ( !inAnAttractor )
				cycleLength++;
		}

		if ( !inAnAttractor )
			timeStep++;
	}

	return( timeStep );
}

#endif
