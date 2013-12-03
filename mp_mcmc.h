#ifndef MPMCMC
#define MPMCMC
/*
 *  mp_mcmc.h
 *
 *
 * Miss Piggy MCMC: Implementation of Metropolis-Hastings
 * 
 *
 * Ivo Siekmann, 2/11/2010
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "mp_parameter.h"
#include "likelihood.h"

/* Was the last move just accepted? */
int parameterUpdated();

/* Metropolis-Hastings sampler: output to files p_fp, ip_fp, diff
   prior and posterior, initial values for p0 and i0, iterations */
void iterateDiff(FILE *p_fp, FILE *ip_fp, FILE *like_fp,
		 double (*dPrior)(const parameters *p, int nP,
				  const intparameters *ip, int nIP,
				  likelihood *L),
		 double (*dPosterior)(const parameters *p, int nP,
				      const intparameters *ip, int nIP,
				      likelihood *L),
		 const parameters *p0, int nP0,
		 const intparameters *ip0, int nIP0,
		 int seed, int nIter);

/* Metropolis-Hastings sampler: output to files p_fp, ip_fp, diff
   prior and posterior, initial values for p0 and i0, iterations */
void iterateDoubleDiff(FILE *p_fp, FILE *like_fp,
		       double (*dPrior)(const parameters *p, int nP,
				  const intparameters *ip, int nIP),
		       double (*dPosterior)(const parameters *p, int nP,
				  const intparameters *ip, int nIP),
		       parameters *p0, int nP0,
		       int seed, int nIter);

/* SHOULD BE PROGRAMMED AS A CALL TO THE MORE GENERAL iterateDiff() */
/* Metropolis-Hastings sampler: output to files p_fp, prior and
   posterior, initial values for p0, iterations */
void iterateDouble(FILE *p_fp, FILE *like_fp,
		   double (*dPrior)(const parameters *p, int nP, 
				    likelihood *L),
		   double (*dPosterior)(const parameters *p, int nP, 
					likelihood *L),
		   parameters *p0, int nP0,
		   int seed, int nIter);

/* Metropolis-Hastings sampler: output to files p_fp, ip_fp, diff
 * prior and posterior, initial values for p0 and i0,
 * iterations. Convention for prior/posterior: Both return the
 * difference and set values in likelihood which is given as an
 * argument.
 */
void iterate(FILE *p_fp, FILE *ip_fp, FILE *like_fp,
	     double (*dPrior)(const parameters *p, int nP, 
			      const intparameters *ip, int nIP, 
			      likelihood *L),
	     double (*dPosterior)(const parameters *p, int nP, 
				  const intparameters *ip, int nIP,
				  likelihood *L),
		 const parameters *p0, int nP0,
		 const intparameters *ip0, int nIP0,
	     int seed, int nIter);

/* change output routine */
void setOutput(void (*newOutput)(FILE *fp, int iteration, 
				 const parameters *p, int nP));

/* outputs all paramters in file fp */
void generalOutput(FILE *fp, int iteration, const parameters *p, int nP);

/* outputs all parameters in file fp, starting from index k0 */
void generalOutputK0(FILE *fp, int iteration, 
		     const parameters *p, int k0, int nP);

/* change output routine */
void setPrintTitles(void (*newPT)(FILE *fp,const parameters *p, int nP));

/*Output titles and initial values*/
void generalPrintTitles(FILE *p_fp, const parameters *p, int nP);

/* Set new sample generator */
void setSample(void (*newSample)(const gsl_rng *r, 
		      parameters *p, int,
				 intparameters *ip, int));

/* Sample parameters according to constraints */
void sampleConstrained(const gsl_rng *r, 
		       parameters *p, int nP,
		       intparameters *ip, int nIP);

/* sample one randomly chosen parameter, constrained by max and min */
void sampleOneDoubleConstrained(const gsl_rng *r, 
				parameters* p, int nPar,
				intparameters *ip, int nIPar);

/* sets routine for copying double parameters */
void setCopyPar(void (*copy)(parameters *, int));

/* Metropolis-Hastings sampler: output to files p_fp, ip_fp, diff
 * prior and posterior, initial values for p0 and i0,
 * iterations. Convention for prior/posterior: Both return the
 * difference and set values in likelihood which is given as an
 * argument.
 *
 * sample: Number of parameters in proposal is given by values in pointers
 */
void iterateRJ(char *k_prefix,
	       char *p_prefix, char *ip_prefix, char *like_prefix,
	       void (*sample)(const gsl_rng *r, 
			      parameters *p, int actP, 
			      int *propP, int nP,
			      intparameters *ip, int actIP, 
			      int *propIP, int nIP),
	       double (*dPrior)(const parameters *p, 
				int actP, int propP, int nP, 
				const intparameters *ip, 
				int actIP, int propIP, int nIP, 
				likelihood *L),
	       double (*dLikelihood)(const parameters *p, 
				     int actP, int propP, int nP, 
				     const intparameters *ip, 
				     int actIP, int propIP, int nIP,
				     likelihood *L),
	       const parameters *p0, int actP0, int minP, int nP0,
	       const intparameters *ip0, int actIP0, int minIP, int nIP0,
	       int seed, int nIter);

/* SHOULD BE PROGRAMMED AS A CALL TO THE MORE GENERAL iterateDiff() */
/* Metropolis-Hastings sampler: output to files p_fp, prior and
   posterior, initial values for p0, iterations */
/* sample returns 0 if p1 was moved and 1 if p2 was moved */
void iterateDoubleTwalk(FILE *p_fp, FILE *like_fp,
			int (*sample)(const gsl_rng *r, 
				      parameters *p1, int,
				      parameters *p2, int),
			double (*prior)(const parameters *p, int nP, 
					likelihood *L),
			double (*likely)(const parameters *p, int nP, 
					 likelihood *L),
			parameters *p1, int nP1,
			parameters *p2, int nP2,
			int seed, int nIter);

/* Set temperature for simulated tempering */
void setTemperature(double Tpar, double Tprop);

double getTpar();


double getTprop();
#endif
