/* ---------------------------------------------------
   FILE:     simParams.h
	AUTHOR:   Josh Bongard
	DATE:     October 5, 2000
	FUNCTION: This class contains all miscellaneous
				 data and functions for this simulation.
 -------------------------------------------------- */

#include "iostream"
#include "fstream"

using namespace std;

#ifndef _SIM_PARAMS_H
#define _SIM_PARAMS_H

class SIM_PARAMS {

public:
	int		randSeed;

	time_t		startTime;

	ofstream *outFile;
	int currentAvailableID;
	
	int phase;

public:
	SIM_PARAMS(int argc, char **argv);
	~SIM_PARAMS(void);
	int    BinaryToInteger(int b1, int b2);
	void   CloseDataFiles(void);
	void   DirectoryMake(char *dirName);
	void   FileCreate(char *fileName);
	void   FileDelete(char *fileName);
	int    FileExists(char *fileName);
	void   FileRename(char *src, char *dest);
	int    FlipCoin(void);
	int    GrayCode(int b1, int b2);
        int    GrayCode(int b1, int b2, int b3);
	double Min(double x, double y, double z);
	void   ParseParameters(int argc, char **argv);
	double Rand(double min, double max);
	double RandGaussian(void);
	int    RandInt(int min, int max);
	double Scale(double value, double min1, double max1,
		         double min2, double max2);
	void   WaitForFile(char *fileName);
};

#endif
