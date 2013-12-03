/*
 * bernoulliDist.c
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

#ifndef HPC
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#else
#include "gsl/gsl_rng.h"
#include "gsl/gsl_sf.h"
#endif

#define OUT stdout
#define ERR stderr

/* Prints n values from int-vector v to file fp starting from index
   n0 */
void fprintInts(FILE *fp, const int* v, int n0, int n) {
  int i;

  for(i=n0; i<n0+n-1; i++) {
    fprintf(fp, "%i\t",v[i]);
    /* fprintf(fp, "%i",v[i]); */
  }
  fprintf(fp,"%i\n",v[n0+n-1]);
}

/* Generates a sequence of length L starting from index k0 where 0,1
   are sampled with probability p from a random generator ran */
void bernoulliProcess(const gsl_rng *ran, int *B, double p,  int k0, int L) {
  int k;

  for (k=k0; k< k0+L; k++ ) {
    B[k] = (gsl_rng_uniform(ran) < p) ? 1:0;
  }
}

/* computes cumulative distribution function from Bernoulli sequence
   dist of length N */
void bernoulliCDF(int *cdf, const int *dist, int n) {
  int i;
  cdf[0]=dist[0];
  for(i=1; i<n; i++) {
    cdf[i]=cdf[i-1]+dist[i];
  }
}

/* Converts data set to binomial CDF by thresholding */
void double2CDF(int *cdf, const double *data, 
		double threshold, int n) {
  int i;
  cdf[0]=data[0]<threshold?1:0;
  for(i=1; i<n; i++) {
    cdf[i]=cdf[i-1]+(data[i]<threshold?1:0);
  }
}

/* Analytical solution for change point for given CDF and change point
   k0 as well as hyperparameters alpha and beta */
double probChangePoint(const int *CDF, int k0, int n, 
		       double alpha, double beta) { 
  int succIn1=CDF[k0], failIn1 = k0-CDF[k0];
  int succIn2=CDF[n-1]-CDF[k0], failIn2 = n-k0-succIn2;
  int status;
  double p1,p2;
  gsl_sf_result result;

  /* fprintf(OUT, "Successes until k0 = %i: %i\n", k0, succIn1); */
  /* fprintf(OUT, "Failures: %i\n", failIn1); */
  /* fprintf(OUT, "\nSuccesses after k0 = %i: %i\n", k0, succIn2); */
  /* fprintf(OUT, "Failures: %i\n", failIn2); */

  status = gsl_sf_beta_e (succIn1+alpha+1, failIn1+beta+1, &result);
  if(status != GSL_SUCCESS) {
    fprintf(ERR, "Evaluation of beta function B(%f,%f) failed.\n",
	    succIn1+alpha+1, failIn1+beta+1);
    exit(1);
  }
  p1=result.val;

  status = gsl_sf_beta_e (succIn2+alpha+1, failIn2+beta+1, &result);
  if(status != GSL_SUCCESS) {
    fprintf(ERR, "Evaluation of beta function B(%f,%f) failed.\n",
	    succIn2+alpha+1, failIn2+beta+1);
    exit(1);
  }
  p2=result.val;

  return p1*p2;
}

/* Likelihood of different choices of k0 */
void distK0(const int *CDF, double * dist, int n, double alpha, double beta) {
  int i;
  double sum=0;
  for(i=0;i<n; i++) {
    dist[i]=probChangePoint(CDF, i, n, alpha, beta);
    sum +=dist[i];
  }
  for(i=0;i<n; i++) {
    dist[i]/=sum;
  }
}

/* Generates test data for nCP change points and saves them in a trace
   CDF of length nData */
void generateTestData(const gsl_rng*r, int *CDF, 
		      const int *k, const double *p, 
		      int nCP, int nData) {
  int i;
  int *process=malloc(nData*sizeof(int));
    fprintf(OUT, "generateTestData():\n");
    for(i=0; i<nCP; i++) {
      fprintf(OUT, "\tk%i=%i", i+1, k[i]);
    }
    fprintf(OUT, "\n");
    
    for(i=0; i<nCP; i++) {
      fprintf(OUT, "p(%i)=%f\t", i, p[i]);
    }
    fprintf(OUT, "p(%i)=%f\n", nCP, p[nCP]);

  for(i=0; i<=nCP; i++) {
    int k0=(i-1>=0)?k[i-1]:0;
    int k1=(i<nCP)?k[i]:nData;
    bernoulliProcess(r, process, p[i], k0, k1-k0);
  }
  bernoulliCDF(CDF, process, nData);
  free(process);
}
