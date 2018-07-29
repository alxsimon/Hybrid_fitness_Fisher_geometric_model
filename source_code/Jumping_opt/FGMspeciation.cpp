using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <vector>
#include "header.h"
#include "MersenneTwister.h"

/*----------------------------------------------------------------------------
 --Recursion--
 Simulates an initial population of size N evolving under stabilizing selection in a n-trait space.
 After reaching mutation-selection-drift equilibirum, the initial population is split into two allopatric populations of size N/2,
 which experience several cycles (ncycles) of environmental change (diff1 and diff2).
 During each cycle nmuts mutations appeared, and some reached fixation allowing populations to track the optima.
 The selection coefficient of mutations that reach fixation and their effect on each trait is recorded.

 Mutations are rare (we consider the fate of mutations one at a time), unique (infinite sites) and affect all traits (complete pleiotropy).
 Their magnitude (mutsize) is drawn from some distribution.

 --Model parameters--
 N: size of the initial population (halved after split)
 n: number of traits under selection
 selcoef: approximate mean selection coefficient at the optimum
 diff1: change in optimum along the zeroth axis, in the first daugther population
 diff2: change in optimum along the zeroth axis, in the second daugther population
 k: curvature of the fitness function
 nmuts: expected number of mutants per cycle
 ncycles: number of cycles
 MVN: Will we use Multivariate Normal mutation (MVN=1) or an exponential distribution (MVN=0) ?
 sigma: overall strength of selection
 replic: number of replicates
 -------------------------------------------------------------------------------*/


/*---------------------------- Global variables --------------------------------*/
MTRand rnd;
FILE * fichierE;

/*--------------------------- Utilitary functions ------------------------------*/
/* Opens input file */
void ouvrirFichierE()
{
	fichierE = fopen(fichierLecture,"r");
}
/* Reads parameter values from input file, returns 1 if end of input file, else returns 0 */
bool lireFichier(int &Nr, int &nr, double &selcoefr, double &diff1r, double &diff2r, double &kr,
				 double &nmutsr, int &ncyclesr, int &MVNr, double &sigmar, int &replicr)
{
	int x;
	bool term;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
	//each parameter set must start with *
	if (x == EOF)
		term = true;
	else
	{
		fscanf(fichierE,"%d ",&Nr);
		fscanf(fichierE,"%d ",&nr);
		fscanf(fichierE,"%lf ",&selcoefr);
		fscanf(fichierE,"%lf ",&diff1r);
		fscanf(fichierE,"%lf ",&diff2r);
		fscanf(fichierE,"%lf ",&kr);
		fscanf(fichierE,"%lf ",&nmutsr);
		fscanf(fichierE,"%d ",&ncyclesr);
		fscanf(fichierE,"%d ",&MVNr);
		fscanf(fichierE,"%lf ",&sigmar);
		fscanf(fichierE,"%d ",&replicr);

		term = false;
	}
	return term;
}
/* Samples a random number from a Normal distribution */
double gasdev()
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (iset == 0) {
		do {
			v1=2.0*rnd.rand()-1.0;
			v2=2.0*rnd.rand()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
/* Gets exponential  */
double rexp(double mymean)
{
	return(-mymean*log(rnd.rand()));
}
/* Gets factorial */
double factorial(double n)
{
	return (n == 1.0 || n == 0.0) ? 1 : factorial(n - 1.0) * n;
}

/*------------------------------- Internal functions --------------------------------------*/
/* Gets optimum which moves in an oscillating fashion */
double get_opt(double t,double maxval,double nmutsv, double kv)
{
	double tmp;

	if(maxval==0)
	{
		return(0.0);
	} else {
		tmp = maxval/2.0 + sin(2.0*PIE*t/nmutsv-PIE/2)/maxval/2.0;
		if(tmp>0.0)
			return(pow(tmp,2.0/kv));
		if(tmp<0.0)
			return(-pow(-tmp,2.0/kv));
		if(tmp==0.0)
			return(0.0);
	}
}
/* Gets mutation size */
double get_mutsize_MVN(int nv, double kv, double s_barv)
{
	double sigv;
	sigv = s_barv * 2.0;
	sigv *= exp(lgamma((double)nv/2.0) - lgamma(((double)nv+kv)/2.0));
	sigv = pow (sigv, 2.0/kv);
	sigv = sqrt(sigv/2.0);
	return(sigv);
}
double get_mutsize_EXP(double kv, double s_barv)
{
	double lamv;
	lamv = 2.0 * s_barv / factorial(kv);
	lamv = pow (lamv, 1.0/kv);
	return(lamv);
}
/* Generates a new random mutation */
void generatemutation(vector<double> & currentmutation, int nv, double mutsize, int MVNv, double kv)
{
	if(MVNv) //Multivariate normal
	{
		for(int i = 0; i<nv; i++)
			currentmutation[i] = mutsize * gasdev();
	} else { //Transformed exponential
		double r = 0.0;
		double curr_mutsize;
		for(int i = 0; i<nv; i++)
		{
			currentmutation[i] = gasdev();
			r += currentmutation[i]*currentmutation[i];
		}
		r = sqrt(r);
		curr_mutsize = pow(rexp(mutsize),2.0/kv); // This line does the transformation
		curr_mutsize /= r;

		for(int i = 0; i<nv; i++)
			currentmutation[i] *= curr_mutsize;
	}
}
/* Gets trait values */
double get_traitvals(population & pop, double optval,vector<double> & traitvals)
{
	int nmut = pop.fixedmuts.size();

	double z;
	int n = int(pop.fixedmuts[0].size());

	//Zeroth trait has a moving optimum (optval)
	z = 0.0;
	for(int j=0; j<nmut; j++) {
		z += pop.fixedmuts[j][0];
	}
	traitvals[0] = z-optval;

	//Remaining traits have a fixed optimum (0.0) - stabilizing selection
	for(int i=1; i<n; i++){
		z = 0.0;
		for(int j=0; j<nmut; j++) {
			z += pop.fixedmuts[j][i];
		}
		traitvals[i] = z;
	}
}
/* Calculates the squared distance from the optimum */
double calc_dist2(population & pop, double optval)
{
	int nmut = pop.fixedmuts.size();
	if(nmut > 0)
	{
		double z;
		double ztot = 0.0;
		int n = int(pop.fixedmuts[0].size());

		//Zeroth trait has a moving optimum (optval)
		z = 0.0;
		for(int j=0; j<nmut; j++) {
			z += pop.fixedmuts[j][0];
		}
		ztot += (z-optval) * (z-optval);

		//Remaining traits have a fixed optimum (0.0) - stabilizing selection
		for(int i=1; i<n; i++){
			z = 0.0;
			for(int j=0; j<nmut; j++) {
				z += pop.fixedmuts[j][i];
			}
			ztot += z*z;
		}
		return(ztot);
	}
	return(optval*optval);
}
/* Gets fitness */
double getfitness(double z,double kv, double sigmav)
{
	return(exp(-sigmav * pow(z, kv/2.0)));
}
/* Saves selection coefficient of mutations fixed and their effect on each trait */
void save_selcoef(ofstream * fout, int which_cycle, int which_pop, int which_time, double optval, double d2b, double d2a, double wb,double wa, double s, vector<double> & currentmutation, int nv)
{
	(*fout) << which_cycle << "\t" << which_pop << "\t" << which_time << "\t" << optval << "\t" << sqrt(d2b) << "\t" << sqrt(d2a) << "\t" << wb << "\t" << wa << "\t" << s;
	for(int i=0; i<nv;i++)
		(*fout) << "\t" << currentmutation[i];
	(*fout) << "\n";
}

/*------------------------------- The recursion --------------------------------------*/
void recursion(int Nv, int nv, double selcoefv, double diff1v, double diff2v, double kv,
				double nmutsv, int ncyclesv, int MVNv, int repv, double sigmav)
{
	//Variables
	int i, j, k;
	double mutant_z, s, mutsize;
	vector<double>currentmutation(nv,0.0);
	population * pop = new population[3];
	double optval;

	if(MVNv) { //Multivariate normal
		mutsize = get_mutsize_MVN(nv, kv, selcoefv);
	} else { //Exponential
		mutsize = get_mutsize_EXP(2.0, selcoefv); // r_bar
	}

	//Initialises populations
	pop[PARENTALPOP].N = Nv;
	pop[PARENTALPOP].z = 0.0;

	pop[DAUGHTERPOP1].N = int(Nv/2);
	pop[DAUGHTERPOP2].N = Nv - pop[DAUGHTERPOP1].N;

	double t[NUMDAUGHTERPOPS]; // Current time within each cycle
	int which_pop[NUMDAUGHTERPOPS]    = {DAUGHTERPOP1,DAUGHTERPOP2}; // Used to reference populations
	double which_opt[NUMDAUGHTERPOPS] = {diff1v,diff2v}; // Used to reference optima

	//For time length measure
	time_t debut, fin;
	struct tm *ptr;
	debut = time(0);

	char myfilestem[256];
	stringstream filestem;
	filestem << Nv << "_n" << nv << "_sigma" << sigmav << "_s" << selcoefv << "_k" << kv << "_Ncycles" << ncyclesv << "_Nmuts" << nmutsv << "_opt1_" << diff1v << "_opt2_" << diff2v << "_MVN" << MVNv << "_cycle" << ncyclesv << "_" << repv << ".txt";
	filestem >> myfilestem;

	char nomFichier3[256];
	stringstream nomF3;
	nomF3 << "SelCoef_N" << myfilestem;
	nomF3 >> nomFichier3;
	ofstream fout3;
	fout3.open(nomFichier3);
	fout3 << "cycle\tpop\ttcycle\to\tzb\tza\twb\twa\ts";
	for(j=0;j<nv;j++)
		fout3 << "\t" << "dz" << j;
	fout3 << "\n";


	//Evolves the initial population
	for(i = 0; i < nmutsv; i++)
	{
		//(1) Generates a new mutation
		generatemutation(currentmutation,nv,mutsize,MVNv,kv);
		pop[PARENTALPOP].fixedmuts.push_back(currentmutation);

		//(2) Calculates the new position in the phenotypic space
		mutant_z = calc_dist2(pop[PARENTALPOP],0.0);

		//(3) Calculates its selective coefficient
		s = getfitness(mutant_z,kv,sigmav) / getfitness(pop[PARENTALPOP].z,kv,sigmav) - 1.0;

		//(4) The new mutation fixes or not given its standard probablity of fixation
		if( rnd.rand() < (1.0 - exp(-2.0*s))/(1.0 - exp(-2.0*pop[PARENTALPOP].N*s)) ) {
			save_selcoef(&fout3, 0, PARENTALPOP, 0.0, 0.0, pop[PARENTALPOP].z, mutant_z, getfitness(pop[PARENTALPOP].z,kv,sigmav), getfitness(mutant_z,kv,sigmav), s, currentmutation, nv);
			pop[PARENTALPOP].z = mutant_z;
		} else {
			pop[PARENTALPOP].fixedmuts.pop_back();
		}
	}

	//Sets up daughter populations
	pop[DAUGHTERPOP1].fixedmuts = pop[DAUGHTERPOP2].fixedmuts = pop[PARENTALPOP].fixedmuts;
	t[0] = t[1] = 0;

	//Evolves the daughter populations
	for(j = 0; j < ncyclesv; j++)
	{
		for(k = 0; k < NUMDAUGHTERPOPS; k++)
		{
			if(j>0) //Generates the first waiting time for the first cycle
				t[k] = t[k]-nmutsv;
			else
				t[k] = rexp(1.0);

			while(t[k] < nmutsv) //We have not reached the end of the cycle
			{
				//Calculates current distance from optimum
				//optval = get_opt(t[k],which_opt[k],nmutsv,kv);
				optval = which_opt[k];
				pop[which_pop[k]].z = calc_dist2(pop[which_pop[k]],optval);

				//(1) Generates a new mutation
				generatemutation(currentmutation,nv,mutsize,MVNv,kv);
				pop[which_pop[k]].fixedmuts.push_back(currentmutation);

				//(2) Calculates the new position in the phenotypic space
				mutant_z = calc_dist2(pop[which_pop[k]],which_opt[k]);

				//(3) Calculates its selective coefficient
				s = getfitness(mutant_z,kv,sigmav) / getfitness(pop[which_pop[k]].z,kv,sigmav) - 1.0;

				//(4) The new mutation fixes or not given its standard probablity of fixation
				if( rnd.rand() < (1.0 - exp(-2.0*s))/(1.0 - exp(-2.0*pop[which_pop[k]].N*s)) ) {
					save_selcoef(&fout3, j+1, which_pop[k], t[k], optval, pop[which_pop[k]].z, mutant_z, getfitness(pop[which_pop[k]].z,kv,sigmav), getfitness(mutant_z,kv,sigmav), s, currentmutation, nv);
					pop[which_pop[k]].z = mutant_z;
				} else {
					pop[which_pop[k]].fixedmuts.pop_back();
				}
				//Adds wait before next mutation
				t[k] += rexp(1.0);
			}
		}
	}
	fin = time(0);
}
/*-------------------------------- Main function -------------------------------------*/
int main()
{

	//Variables
	int Nt, MVN, n, ncycles, replic;
	double sigma, nmuts, diff1, diff2, k, mutsize;

	//Opens files
	bool fin;
	ouvrirFichierE();
	fin = false;
	int rep;
	do
	{
		//Reads parameter values
		fin = lireFichier(Nt, n, mutsize, diff1, diff2, k, nmuts, ncycles, MVN, sigma, replic);

		if (!fin)
		{
			for(rep = 1; rep <= replic; rep++)
			{
                cout << "replicate number    " << rep << endl;
				//Runs simulation
				recursion(Nt, n, mutsize, diff1, diff2, k, nmuts, ncycles,MVN, rep, sigma);
			}
		}
	} while (!fin);

	//Closes files
	fclose(fichierE);
    return 0 ;
}
