#ifndef _OBJECT_CPP
#define _OBJECT_CPP

#include "math.h"
#include "stdio.h"

#include "object.h"

extern int NUM_OBJECTS;

OBJECT::OBJECT(int objectIndex) {

	int numPositions = int( pow( double(NUM_OBJECTS) , 0.5 ) );

        int numSizes = numPositions;

	int currIndex = 0;

	for (double size = 1.0 ; size >= 0.5 ; size = size - ((1.0 - 0.5) / (double(numSizes)-1.0)) ) {

		for (double position = 2.5 ; position >= -2.5 ; position = position - ((2.5 - -2.5) / (double(numPositions)-1.0)) ) {

			if ( currIndex == objectIndex ) {

				x = position;
				radius = size;
			}

			currIndex = currIndex + 1;
		}
	}

	y = radius;
}

OBJECT::~OBJECT(void) {

}

double OBJECT::Get_X(void) {

	return x;
}

double OBJECT::Get_Y(void) {

	return y;
}

double OBJECT::Get_Radius(void) {

	return radius;
}

void OBJECT::Print(int objectIndex) {

	printf("%d %3.3f %3.3f %3.3f\n",objectIndex,x,y,radius);

}

#endif
