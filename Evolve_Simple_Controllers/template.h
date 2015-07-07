#include "matrix.h"

#ifndef _HILLCLIMBER_H
#define _HILLCLIMBER_H

class HILLCLIMBER {

public:
	MATRIX *parent;
	MATRIX *child;

public:
	HILLCLIMBER(void);
	~HILLCLIMBER(void);
};

#endif
