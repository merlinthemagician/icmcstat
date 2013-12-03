#ifndef CHANGE
#define CHANGE
/*
 * bernoulliOneChangePointMCMC.h
 *
 * MCMC moves for calculating the distribution of changepoint
 * locations and the distribution of open probabilities for the
 * segments defined by the changepoints. It is assumed that open and
 * closed events within segments are distributed according to
 * Bernoulli distributions data set whose parameter p changes
 * instantly from one segment to another.
 *
 * Ivo Siekmann, 2/11/2010
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include "mp_parameter.h"
#include "likelihood.h"


/* Either moves changepoint locations or probabilities. Does NOT add
   or remove changepoints. */
void sampleMixed(const gsl_rng *r, 
		 parameters* p, int nPar,
		 intparameters *ip, int nIPar);

/* Birth/Death proposal or moving changepoints. */
void sampleBirthDeath(const gsl_rng *r,
		      parameters *p, int actP,
		      int *propP, int nPar,
		      intparameters *ip, int actIP,
		      int *propIP, int nIPar);

/* Set number of data points in trace */
void setNdata(int n);

/* Get number of data points in trace */
int getNdata();

/* Set number of data points in trace */
const int *getData();

/* Get number of data points in trace */
void setData(int *CDF, int n);

/* Number of successes observed between k0 and k1 */
int succInterval(const int *CDF, int k0, int k1);

/* Number of failures observed between k0 and k1 */
int failInterval(const int *CDF, int k0, int k1);

/* Get scaling factor for birth/death move */
double getProposalScale();
#endif
