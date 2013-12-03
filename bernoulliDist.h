#ifndef BERNOULLI
#define BERNOULLI
/*
 * bernoulliDist.h
 *
 * Generating Bernoulli distributed CDFs, possibly with one or more
 * change points, analytical distribution for change point.
 * 
 *
 * Ivo Siekmann, 17/10/2010
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>

/* Prints n values from int-vector v to file fp starting from index
   n0 */
void fprintInts(FILE *fp, const int* v, int n0, int n);

/* Generates a sequence of length L starting from index k0 where 0,1
   are sampled with probability p from a random generator ran */
void bernoulliProcess(const gsl_rng *ran, int *B, double p,  int k0, int L);

/* computes cumulative distribution function from Bernoulli sequence
   dist of length N */
void bernoulliCDF(int *cdf, const int *dist, int n);

/* Converts data set to binomial CDF by thresholding */
void double2CDF(int *cdf, const double *data, 
		double threshold, int n);

/* Likelihood of different choices of k0, Bayesian hyperparameters
   alpha and beta for beta prior */
void distK0(const int *CDF, double * dist, int n, double alpha, double beta);

/* Generates test data for nCP change points and saves them in a trace
   CDF of length nData */
void generateTestData(const gsl_rng*r, int *CDF, 
		      const int *k, const double *p, 
		      int nCP, int nData);
#endif
