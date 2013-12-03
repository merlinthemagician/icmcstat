#ifndef CP_LIKELIHOOD
#define CP_LIKELIHOOD
/*
 * cp_likelihood.c
 *
 * Calculates likelihood for an arbitrary number of changepoint
 * locations and open probabilities.
 *
 * Ivo Siekmann, 20/11/2013
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>

/* Calculates likelihood for given locations of changepoints and open
   probabilities. */
double dLikelihood(const parameters *p, 
		   int actP, int propP, int nP, 
		   const intparameters *ip, 
		   int actIP, int propIP, int nIP, 
		   likelihood *L);

#endif
