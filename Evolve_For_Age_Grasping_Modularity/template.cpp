#ifndef _HILLCLIMBER_CPP
#define _HILLCLIMBER_CPP

#include "matrix.h"
#include "hillClimber.h"

HILLCLIMBER::HILLCLIMBER(void) {

	parent = NULL;

	child = NULL;
}

HILLCLIMBER::~HILLCLIMBER(void) {

	if ( parent ) {
		delete parent;
		parent = NULL;
	}

        if ( child ) {
                delete child;
                child = NULL;
        }
}

#endif
