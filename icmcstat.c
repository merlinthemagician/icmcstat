/*
 * icmcstat.c
 *
 * MCMC method for calculating the distribution of a change point k
 * for a bernoulli-distributed data set whose parameter p changes
 * instantly.
 *
 *
 * Ivo Siekmann, 26/04/2012
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_randist.h>

#include "bernoulliDist.h"
#include "cp_proposal.h"
#include "cp_likelihood.h"
#include "cp_prior.h"
#include "mp_mcmc.h"
#include "mp_parameter.h"

#define OUT stdout
#define ERR stderr

/* Number of change points */
#define NCHANGE 250

/* Bernoulli-distributed data with change point */
static int *CDFdata=NULL/* , nCDFdata=1000 */;

static double * probs;/* [NCHANGE+1]; */
static double * probsProp;/* [NCHANGE+1]; */

static int * k;/* [NCHANGE+1]; */
static int * kProp;/* [NCHANGE+1]; */

/* Reads a list of doubles. */
long readDoubles(FILE *fp, double *v, int buf) {
  long k=0;
  double d;
  char buff[20];
  while( (k<buf) && (fgets( buff, 20, fp ))) {
    /* If buff has more than just \n */
    if((strlen(buff)>1) && (sscanf(buff,"%lf", &d)))
      v[k++]=d;
  }
  return k;
}

/*Number of changepoints and probabilities */
static int nProb=NCHANGE+1, nK=NCHANGE;

/*Seed of random number generator */
static int seed=42;

/* Number of iterations */
static int nIter=6000;

/* Step width for moving changepoint locations.*/
static int dK=5;

/* 50 % threshold of channel current */
static double threshold=-20;

/* Prefix for output files */
static char *prefix="RJtest";

/* Filename of data file. */
static char *dataFn="testdata/IPR2_10uMIP3_5mMATP_10nMCa_C4_T02.dat";

/* Random number generator */
static gsl_rng *r;


static parameters *p0/* [NCHANGE+1] */;
static intparameters *ip0/* [NCHANGE] */;

/* Read current data from dataFn, threshold currents and generate
   counts of open events. */
void processData(char *dataFn) {
  long nTr=60000*25;
  double * traceData=malloc(nTr*sizeof(double));
  long nTrace;
  FILE *data_fp=stdin;

  /* If data file is not stdin, open file. */
  if(strcmp(dataFn, "-")!=0) data_fp=fopen(dataFn,"r");

  /* Read and process input data */
  fprintf(OUT, "Processing data...\n");
  traceData = malloc(nTr*sizeof(double));
  if(!traceData) fprintf(ERR, "Out of memory\n"), exit(1);
  nTrace = readDoubles(data_fp, traceData,nTr);
  printf("Read %li lines...\n", nTrace);

  CDFdata=malloc(nTrace*sizeof(int));
  double2CDF(CDFdata, traceData, threshold, nTrace);

  fprintf(OUT, "Generating CDF:\n");
  setData(CDFdata, nTrace);
  fprintf(OUT, "done.\n");
  free(traceData);

  /* If data file is not stdin, close file. */
  if(strcmp(dataFn, "-")!=0) fclose(data_fp);
}

/* Initialise */
void initialise() {
  /* Allocate memory for parameters */
  p0=malloc(nProb*sizeof(parameters));
  ip0=malloc(nK*sizeof(intparameters));

  fprintf(OUT, "Setting parameters...\n");
  initialiseParameters(p0, 0.1, 0, 1, nProb);
  probs=malloc( nProb*sizeof(double));
  probsProp=malloc( nProb*sizeof(double));
  initialisePorg(p0, probs, nProb);
  initialisePprop(p0, probsProp, nProb);

  setParameterNames(p0, "p", nProb);
  /* printParameterNames(OUT, p0, nProb); */
  /* printParameters(OUT,p0, nProb); */

  /* Allocate memory for changepoints */
  k=malloc(nK*sizeof(int));
  kProp=malloc(nK*sizeof(int));
  initialiseIntPorg(ip0, k, nK);
  initialiseIntPprop(ip0, kProp, nK);

  /* Set name prefix and parameters for changepoints */
  setIntParameterNames(ip0, "k", nK);
  /* sample width dK - sample between p-dK and p+dK. */
  initialiseIntParameters(ip0, dK, 0, getNdata()-1, nK);

  /* Initialising random number generator */
  r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r,seed);
  initUniform(r, p0, nProb);
  /* Changepoint locations in ascending order. */
  initIntUniformAscending(r, ip0, 0, getNdata()-1, nK);
  
}


/* Process command line parameters:
 argv[1]: Name of data file or '-' for stdin
 argv[2]: Number of iterations
 argv[3]: seed for random number generator
 argv[4]: (optional) Prefix for output files
*/
void getArgs(int argc, char **argv) {
  char* endptr;
  const char* usage="icmcstat FILE ITERATIONS SEED [OUTPUTPREFIX]";

  if(argc<4) fprintf(ERR, "%s\n", usage), exit(1);
  dataFn=argv[1];  
  printf("Reading data from ");

  if(!strcmp(dataFn, "-")) printf("stdin\n");
  else printf("%s\n", dataFn);

  printf("%s iterations\n", argv[2]);
  nIter = strtol(argv[2], &endptr, 10);

  seed = strtol(argv[3], &endptr, 10);
  printf("Seed for random number generator: %i\n", seed);

  /* Optional parameter*/
  if(argc==5) {
    prefix=argv[4];
    printf("Prefix prepended to output files: %s\n", prefix); 
  }
}

/* Process command line parameters:
 argv[1]: Name of data file or '-' for stdin
 argv[2]: Number of changepoints
 argv[3]: Number of iterations
 argv[4]: seed for random number generator
 argv[5]: (optional) Prefix for output files
*/
void getFixedArgs(int argc, char **argv) {
  char* endptr;
  const char* usage="icmcstatFixed FILE NK ITERATIONS SEED [OUTPUTPREFIX]";

  if(argc<5) fprintf(ERR, "%s\n", usage), exit(1);
  dataFn=argv[1];  
  printf("Reading data from ");

  if(!strcmp(dataFn, "-")) printf("stdin\n");
  else printf("%s\n", dataFn);

  printf("%s nK\n", argv[2]);
  nK = strtol(argv[2], &endptr, 10);
  nProb=nK+1;

  printf("%s iterations\n", argv[3]);
  nIter = strtol(argv[3], &endptr, 10);

  seed = strtol(argv[4], &endptr, 10);
  printf("Seed for random number generator: %i\n", seed);

  /* Optional parameter*/
  if(argc==6) {
    prefix=argv[5];
    printf("Prefix prepended to output files: %s\n", prefix); 
  }
}

/* Dynamically allocates room for s plus additional n characters */
static char *strnsave(const char * s, size_t n) {
  char *p = malloc(strlen(s) + n + 1);
  if(!p) fputs("strsave: out of memory\n", ERR), exit(1);

  return strcpy(p, s);
}

/* Adds prefix to s */
static char* addPrefix(const char *s, const char *prefix) {
  char *out = strnsave(prefix, strlen(s));

  return  strcat(out, s);
}

int main(int argc, char **argv) {
  /* Initial number of changepoints */
  int nK0=1;
#ifdef FIXED
  FILE *pfp, *kfp, *Lfp;
  char *pfn="Probs.dat", *kfn="K.dat", *Lfn="Likelihood.dat";
  prefix=NULL;
#endif

#ifndef FIXED
  getArgs(argc, argv);
#else
  getFixedArgs(argc, argv);
#endif
  processData(dataFn);

  initialise();

#ifdef FIXED
  if(prefix) {
    pfn=addPrefix(pfn, prefix);
    kfn=addPrefix(kfn, prefix);
    Lfn=addPrefix(Lfn, prefix);
  }
  pfp=fopen(pfn,"w");
  kfp=fopen(kfn,"w");
  Lfp=fopen(Lfn,"w");

  setData(CDFdata, getNdata());
  setSample(sampleMixed);
  iterate(pfp, kfp, Lfp, /* dPriorNew */dPriorFixedChangepoints, dFixedLikelihood,p0,nProb, ip0, nK, seed, nIter);
#endif

#ifdef DEFAULT
  iterateRJ(prefix, prefix, prefix, prefix,
	    sampleBirthDeath,
	    dPriorChangepoints,dLikelihood,
	    p0, nK0+1,1, nProb,
	    ip0, nK0,1, nK, seed,  nIter);
#endif

#ifdef GEO
  iterateRJ(prefix, prefix, prefix, prefix,
	    sampleBirthDeath,
	    dPriorChangepointsGeo,dLikelihood,
	    p0, nK0+1,1, nProb,
	    ip0, nK0,1, nK, seed,  nIter);
#endif

#ifdef NEGATIVEBINK1
  iterateRJ(prefix, prefix, prefix, prefix,
	    sampleBirthDeath,
	    dPriorChangepointsNegativeBinK1,dLikelihood,
	    p0, nK0+1,1, nProb,
	    ip0, nK0,1, nK, seed,  nIter);
#endif
#ifdef NEGATIVEBINK2
  iterateRJ(prefix, prefix, prefix, prefix,
	    sampleBirthDeath,
	    dPriorChangepointsNegativeBinK2,dLikelihood,
	    p0, nK0+1,1, nProb,
	    ip0, nK0,1, nK, seed,  nIter);
#endif
#ifdef NJGEONEGATIVEBINK1
  iterateRJ(prefix, prefix, prefix, prefix,
	    sampleBirthDeath,
	    dPriorChangepoints_nJGeoNegativeBinK1,dLikelihood,
	    p0, nK0+1,1, nProb,
	    ip0, nK0,1, nK, seed,  nIter);
#endif
#ifdef NJGEONEGATIVEBINK2
  iterateRJ(prefix, prefix, prefix, prefix,
	    sampleBirthDeath,
	    dPriorChangepoints_nJGeoNegativeBinK2,dLikelihood,
	    p0, nK0+1,1, nProb,
	    ip0, nK0,1, nK, seed,  nIter);
#endif

  return 0;  
}
 
