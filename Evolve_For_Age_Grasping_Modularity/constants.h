#ifndef _CONSTANTS_H
#define _CONSTANTS_H

// AFPO constants

int			POP_SIZE			= 100;

int			MAX_NON_DOMINATED_GENOMES	= int(double(POP_SIZE)/4.0);

int			NUM_GENERATIONS			= 1000000;
// int			NUM_GENERATIONS			= 1000;

double			MUTATION_PROBABILITY		= 0.05;

// Environment constants

int			NUM_OBJECTS			= 4;

// Robot constants

int			NUM_BODY_PARTS			= 3 + 2 + 2;

double			MIN_ANGLE			= -70.0;
double			MAX_ANGLE			= +70.0;

// Random boolean network constants 

int                     NUM_NODES_PER_BODY_PART         = 2;

int			NUM_NODES			= NUM_BODY_PARTS * NUM_NODES_PER_BODY_PART;

int			START_NUM_CONNECTIONS		= 2 * NUM_NODES;

int			NUM_POSSIBLE_ANGLES		= pow(2,NUM_NODES_PER_BODY_PART);

int			NUM_UPDATES			= 4;

// int			NUM_PERTURBATIONS		= 1 + NUM_NODES + int( ( double(NUM_NODES) * ( double(NUM_NODES) - 1.0 ) ) / 2.0 );

int			NUM_PERTURBATIONS		= 1 + NUM_NODES;

#endif
