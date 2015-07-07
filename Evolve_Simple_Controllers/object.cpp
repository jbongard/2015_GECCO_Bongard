#ifndef _OBJECT_CPP
#define _OBJECT_CPP

#include "object.h"

OBJECT::OBJECT(int objectIndex) {

	if ( objectIndex==0 ) {
	
		x = 2.5;
		radius = 1.0;
	}

        else if ( objectIndex==1 ) {

                x = -2.5;
                radius = 1.0;
        }

        else if ( objectIndex==2 ) {

                x = 2.5;
                radius = 0.5;
        }

        else {

                x = -2.5;
                radius = 0.5;
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

#endif
