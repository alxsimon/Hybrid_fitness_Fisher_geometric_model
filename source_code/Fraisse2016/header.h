using namespace std;

#ifndef HEADER_H
#define HEADER_H

#include <vector>
#include <iostream>
#include "MersenneTwister.h"


/*------------------------------- Global variables ------------------------------------*/
extern MTRand rnd;
extern FILE * fichierE;


/* The population will be represented by a table of the following elements */
struct population
{
	int N;
	double z;
	vector<vector<double> > fixedmuts;
};


/*------------------------------ Simulation settings ----------------------------------*/
#define fichierLecture "parameters"			// name of parameter file
#define NUMDAUGHTERPOPS 2					// number of daughter populations
#define PARENTALPOP 0						// index for parental population
#define DAUGHTERPOP1 1						// index for the first daugther population	
#define DAUGHTERPOP2 2						// index for the second daugther population

#define PIE 3.14159265359

#endif
