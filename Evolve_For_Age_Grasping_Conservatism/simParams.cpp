/* ---------------------------------------------------
   FILE:     simParams.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 5, 2000
	FUNCTION: This class contains all miscellaneous
				 data and functions for this simulation.
 -------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"

#ifndef _SIM_PARAMS_CPP
#define _SIM_PARAMS_CPP

#include "simParams.h"

extern int		RANDOM_SEED;

SIM_PARAMS::SIM_PARAMS(int argc, char **argv) {

	startTime = time(NULL);

	currentAvailableID = 0;

	randSeed = RANDOM_SEED;

	ParseParameters(argc,argv);

	srand(randSeed);

	char fileName[100];
	sprintf(fileName,"Data/runData_%d.dat",randSeed);
	outFile = new ofstream(fileName);
	outFile->close();

	phase = -1;
}

SIM_PARAMS::~SIM_PARAMS(void) {

}

int SIM_PARAMS::BinaryToInteger(int b1, int b2) {

        if (            (b1==-1) && (b2==-1) )
                return( 0 );
        else if (       (b1==-1) && (b2==+1) )
                return( 1 );
        else if (       (b1==+1) && (b2==-1) )
                return( 2 );
        else // if (       (b1==+1) && (b2==+1) )
                return( 3 );
}

void  SIM_PARAMS::CloseDataFiles(void) {

}

void  SIM_PARAMS::DirectoryMake(char *dirName) {

	char command[500];

	sprintf(command,"rm -r %s",dirName);
	system(command);

	sprintf(command,"mkdir %s",dirName);
	system(command);
}

void  SIM_PARAMS::FileCreate(char *fileName) {

}

void  SIM_PARAMS::FileDelete(char *fileName) {

	char command[100];

	if ( FileExists(fileName) ) {
		sprintf(command,"rm %s",fileName);
		system(command);
	}
}

int   SIM_PARAMS::FileExists(char *fileName) {

	int exists;

	ifstream *inFile = new ifstream(fileName,ios::in);

	if ( inFile->good() )
		exists = true;

	else
		exists = false;

	inFile->close();
	delete inFile;
	inFile = NULL;

	return( exists );
}

void SIM_PARAMS::FileRename(char *src, char *dest) {

	char command[200];

	sprintf(command,"mv %s %s",src,dest);

	system(command);
}

int SIM_PARAMS::FlipCoin(void) {

	return( Rand(0.0,1.0) < 0.5 );
}

int SIM_PARAMS::GrayCode(int b1, int b2) {

	if (		(b1==-1) && (b2==-1) )
		return( 0 );
	else if ( 	(b1==-1) && (b2==+1) )
		return( 1 );
	else if ( 	(b1==+1) && (b2==+1) )
		return( 2 );
	else if ( 	(b1==+1) && (b2==-1) )
		return( 3 );

	return( 0 );
}

int SIM_PARAMS::GrayCode(int b1, int b2, int b3) {

        if (            (b1==-1) && (b2==-1) && (b3==-1) )
                return( 0 );
        else if (       (b1==-1) && (b2==-1) && (b3==+1) )
                return( 1 );
        else if (       (b1==-1) && (b2==+1) && (b3==+1) )
                return( 2 );
        else if (       (b1==-1) && (b2==+1) && (b3==-1) )
                return( 3 );
        else if (       (b1==+1) && (b2==+1) && (b3==-1) )
                return( 4 );
        else if (       (b1==+1) && (b2==+1) && (b3==+1) )
                return( 5 );
        else if (       (b1==+1) && (b2==-1) && (b3==+1) )
                return( 6 );
        else // if (       (b1==+1) && (b2==-1) && (b3==-1) )
                return( 7 );
}


double SIM_PARAMS::Min(double x, double y, double z) {

	if ( (x<=y) && (x<=z) )
		return(x);
	else
		if ( (y<=x) && (y<=z) )
			return(y);
		else
			return(z);
}

void  SIM_PARAMS::ParseParameters(int argc, char **argv) {

	int currParam;

	for(currParam=0;currParam<argc;currParam++) {

		if ( strcmp(argv[currParam],"-r") == 0 )
			randSeed = atoi(argv[currParam+1]);
	}
}

double SIM_PARAMS::Rand(double min, double max) {

	double zeroToOne = ((double)rand()) / RAND_MAX;
	double returnVal;

	returnVal = (zeroToOne * (max-min)) + min;

	if ( returnVal < min )
		returnVal = min;

	if ( returnVal > max )
		returnVal = max;

	return returnVal;
}

double SIM_PARAMS::RandGaussian(void) {

	double w = 1.01;
	double x1, x2;

	while ( w >= 1.0 ) {
		x1 = 2.0*Rand(0,1) - 1;
		x2 = 2.0*Rand(0,1) - 1;
		w = x1*x1 + x2*x2;
	}	
	w = sqrt( (-2.0*log(w))/w );
	
	return( x1*w );
}

int SIM_PARAMS::RandInt(int min, int max) {

	int val;

	if ( min == max )
		val = min;
	else
		val = (rand() % (max-min+1)) + min;

	if ( val < min )
		val = min;

	if ( val > max )
		val = max;

	return( val );
}

double SIM_PARAMS::Scale(double value, double min1, double max1,
								 double min2, double max2) {

	if ( min1 < 0 )
		value = value - min1;
	else
		value = value + min1;

	return( (value*(max2-min2)/(max1-min1)) + min2 );
}

void SIM_PARAMS::WaitForFile(char *fileName) {

	while ( !FileExists(fileName) );
}

#endif
