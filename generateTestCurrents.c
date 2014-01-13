/*
 *  generateTestsCurrents.c
 *
 *  Simulates ion channel currents as a Bernoulli process. 
 * 
 *
 * Ivo Siekmann, 19/12/2013
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Maximum number of change points */
#define MAX 10

/* Number of data points */
static int nTrace=10000;

static int nChange=1;

/* Locations of change points, last position is length of trace. */
static int k[MAX+1];

/* Open probabilities of segments */
static double prob[MAX+1];

/* Open current */
static double Io = -40;

/* Standard deviation of normally-distributed noise */
static double sigma_Io=5;

/* Prints test currents to file fp. */
void generateCurrents(FILE *fp, 
		      const gsl_rng *r,
		      const int *k, const double *prob, int nChange,
		      double Io, double sig_Io) {

  int i, j=0;
  for (i=0; i<=nChange+1; i++) {
    int k0=i-1>=0?k[i-1]:0;
    int k1=k[i];
    for(j=k0; j<k1; j++) {
      double current=gsl_ran_gaussian(r, sig_Io);
      if(gsl_rng_uniform(r) < prob[i]) current+=Io;

      fprintf(fp, "%g\n", current);
    } 
  }
}

#define OUT stdout
#define ERR stderr

/* reads at most max double values from file fp. */
int readDoubles(FILE *fp, double * v, int max) {
  int i;
  double d;

  while(  (i<max) && (fscanf(fp, "%lf", &d) != EOF) ) {
    v[i++] = d;
    fprintf(OUT, "%g\n", d);
  }
  return i;
}

/* reads at most max double values from file fp. */
int readInts(FILE *fp, int * v, int max) {
  int i;
  int d;

  while( (i<max) && (fscanf(fp, "%d", &d) != EOF) ) {
    v[i++] = d;
    fprintf(OUT, "%i\n", d);
  }
  return i;
}

int main(int argc, char **argv) {
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
  const int minPar=2, maxPar=7;
  const char * usage="generateTestCurrents kFile pFile [seed] [outfile] [openCurr] [stdDev]";
  char *kFn, *pFn, *outFn="-";
  FILE *kFp, *prFp, *outFp;
  int seed = 42;
  int nK, nP;
  int i;
  /* k[0]=3000, k[1]=5000, k[2]=nTrace; */
  /* prob[0]=0.4, prob[1]=0.6, prob[2]=0.8; */
  
  if((argc <= minPar) || (argc>maxPar)) fprintf(ERR, "%s\n", usage), exit(1);
  kFn=argv[1];
  kFp=fopen(kFn, "r");
  nK=readInts(kFp, k, MAX);
  fclose(kFp);

  pFn=argv[2];
  prFp=fopen(pFn, "r");
  nP=readDoubles(prFp, prob, MAX);
  fclose(prFp);

  if(nK != nP) fprintf(ERR, "Expected nP=%i probabilities for nK=%i changepoints\n", nK, nK-1), exit(1);

  if(argc>3) {
    seed=atoi(argv[3]);
  }

  if(argc>4) {
    outFn=argv[4];
  }
  if( strcmp(outFn, "-") != 0 ) {
    outFp=fopen(outFn, "w");
  }
  else {
    outFp=OUT;
  }

  if(argc > 5) {
    Io=atof(argv[5]);
  }

  if(argc > 6) {
    sigma_Io=atof(argv[6]);
  }

  fprintf(ERR, "Generating test data set with %i data points...\n", k[nK-1]);
  fprintf(ERR, "k = (");
  for(i=0; i<nK-2; i++) {
    fprintf(ERR, "%i, ", k[i]);
  }
  fprintf(ERR, "%i)\n", k[nK-2]);


  fprintf(ERR, "p = (");
  for(i=0; i<nP-1; i++) {
    fprintf(ERR, "%g, ", prob[i]);
  }
  fprintf(ERR, "%g)\n", prob[nP-1]);

  fprintf(ERR, "seed = %i\n", seed);
  fprintf(ERR, "Io=%g (%g)\n", Io, sigma_Io);
  fprintf(ERR, "Output written to: %s\n", outFn);

  generateCurrents(outFp, r, k, prob,  nChange, Io,  sigma_Io);

  return 0;  
}
 
