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

#ifndef WINDOWS
/* For command line parameters */
#include <unistd.h>
#endif

#ifdef WINDOWS
#include "getopt.h"
#endif

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
#define NCHANGE 1000

/* Bernoulli-distributed data with change point */
static int *CDFdata=NULL/* , nCDFdata=1000 */;

static double * probs;/* [NCHANGE+1]; */
static double * probsProp;/* [NCHANGE+1]; */

static int * k;/* [NCHANGE+1]; */
static int * kProp;/* [NCHANGE+1]; */

#define MAX_LEN 10000

char *lastLine(FILE *fp, char *buf, size_t max_len) {
  char *last_newline; 
  char *last_line;
      
  fseek(fp, -max_len, SEEK_END);
  fread(buf, max_len-1, 1, fp);

  buf[max_len-1] = '\0';
  /* TODO: Check for NULL in case '\n' is not found. */
  last_newline= strrchr(buf, '\n');
  last_line = last_newline+1;

  /* printf("%s\n", last_line); */
  return last_line;
}

/* Convert s containing ints to an int vector*/
size_t str2ints(char* s, int *nums) {
  int i=0;
  char *p, *endP=NULL;
  int item;
  for(p=s; ;p=endP)  {
    item=strtol(p, &endP, 10);
    if(p == endP) break;
    nums[i++]=item;
  }
  return i;
}

/* Convert s containing ints to an int vector*/
size_t str2doubles(char* s, double *nums) {
  int i=0;
  char *p, *endP=NULL;
  double item;
  for(p=s; ;p=endP)  {
    item=strtod(p, &endP);
    if(p == endP) break;
    nums[i++]=item;
  }
  return i;
}

/* Reads a list of doubles. */
long readDoubles(FILE *fp, double *v, int buf) {
  long k=0;
  /* long double d; */
  double d;
  char buff[100];
  while( (k<buf) && (fgets( buff, 100, fp ))) {
    /* If buff has more than just \n */
    /* if((strlen(buff)>1) && (sscanf(buff,"%Lg", &d))) */
    if((strlen(buff)>1) && (sscanf(buff,"%lg", &d)))
      v[k++]=d;
  }
  return k;
}

/* Reads a list of ints. */
long readInts(FILE *fp, int *v, int buf) {
  long k=0;
  int i;
  char buff[20];
  while( (k<buf) && (fgets( buff, 20, fp ))) {
    /* If buff has more than just \n */
    if((strlen(buff)>1) && (sscanf(buff,"%i", &i)))
      v[k++]=i;
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
static double threshold=/* -20 */-20;

/* Prefix for output files */
static char *prefix="results/RJtest";

/* Filename of data file. */
static char *dataFn="testdata/IPR2_10uMIP3_5mMATP_10nMCa_C4_T02.dat";

/* Random number generator */
static gsl_rng *r;

/* Should run be restarted? */
static int restart=0;

/* Files containing restart data for k and p in last line */
static char* re_K;
static char* re_P;


static parameters *p0/* [NCHANGE+1] */;
static intparameters *ip0/* [NCHANGE] */;

/* Read current data from dataFn, threshold currents and generate
   counts of open events. */
void processData(char *dataFn) {
  long nTr=60000*25;
  double * traceData=malloc(nTr*sizeof(double));
  int *intData=malloc(nTr*sizeof(int));
  long nTrace;
  FILE *data_fp=stdin;

  /* If data file is not stdin, open file. */
  if(strcmp(dataFn, "-")!=0) data_fp=fopen(dataFn,"r");

  /* Read and process input data */
  fprintf(OUT, "Processing data...\n");

  #ifndef BERNOULLIDATA
  traceData = malloc(nTr*sizeof(double));
  if(!traceData) fprintf(ERR, "Out of memory\n"), exit(1);
  nTrace = readDoubles(data_fp, traceData,nTr);
  printf("Read %li lines...\n", nTrace);

  CDFdata=malloc(nTrace*sizeof(int));
  double2CDF(CDFdata, traceData, threshold, nTrace);

  fprintf(OUT, "Generating CDF:\n");
  free(traceData);
#else
  nTrace = readInts(data_fp, intData,nTr);
  printf("Read %i lines...\n", nTrace);

  CDFdata=malloc(nTrace*sizeof(int));
  bernoulliCDF(CDFdata, intData, nTrace);

  fprintf(OUT, "Generating CDF:\n");
  free(intData);
#endif

  setData(CDFdata, nTrace);
  fprintf(OUT, "done.\n");

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

static double (* prior)(const parameters *p, 
			int actP, int propP, int nP, 
			const intparameters *ip, 
			int actIP, int propIP, int nIP, 
			likelihood *L) = dPriorChangepoints;

static const char *prior_nJ_names[] = {"Poisson", "geometric"};

static const char *prior_Loc_names[] = {"Order statistics", "negative binomial (k=1)", "negative binomial (k=2)"};

/* Process command line parameters:
 argv[1]: Name of data file or '-' for stdin
 argv[2]: Number of iterations
 argv[3]: seed for random number generator
 argv[4]: (optional) Prefix for output files

 argv[5]: (optional) For RESTART: Initial changepoint locations
 argv[6]: (optional) For RESTART: Initial success probabilities
*/
void getArgs(int argc, char **argv) {
  const int nArgs=4;
  const int nOpt=3;
  
  char* endptr;
  const char* usage="icmcstat FILE ITERATIONS SEED [OUTPUTPREFIX] [K0FILE] [P0FILE]";
  double alpha=1, beta=1;

  int opt;

  int prior_nJ=0;  /* 0: Poisson, 1: geometric */
  int prior_Loc=0; /* 0: order statistic, 1: negative binomial (k=1), 2: negative binomal (k=2) */
  
  while( (opt=getopt(argc, argv, "gpq")) != -1 ){
      switch(opt) {
      case 'g': 
	prior_nJ=1;
	break;
      
      case 'p': 
	if(!prior_Loc) prior_Loc=1;
	else fprintf(ERR, "ERROR: Only either of the options -p or -q option can be used\n"), exit(1);
	break;
      
      case 'q': 
	if(!prior_Loc) prior_Loc=2;
	else fprintf(ERR, "ERROR: Only either of the options -p or -q option can be used\n"), exit(1);
	break;
      
	//default: fprintf(ERR, "ERROR: Unknown options\n"), exit(1);
      }
  }

  fprintf(ERR, "Priors selected:\n");
  switch(prior_nJ) {
  case 0: {
    switch(prior_Loc) {
    case 0: {
      prior=dPriorChangepoints;
      break;
    }
    case 1: {
      prior=dPriorChangepointsNegativeBinK1;
      break;
    }
    case 2: {
      prior=dPriorChangepointsNegativeBinK2;
      break;
    }
    }
    break;
  }
  case 1: {
    switch(prior_Loc) {
    case 0: {
      prior=dPriorChangepointsGeo;
      break;
    }
    case 1: {
      prior=dPriorChangepoints_nJGeoNegativeBinK1;
      break;
    }
    case 2: {
      prior=dPriorChangepoints_nJGeoNegativeBinK2;
      break;
    }
    }
    break; 
  }
  }

  fprintf(ERR, "Number of changepoints: %s, Changepoint locations: %s\n",
	  prior_nJ_names[prior_nJ], prior_Loc_names[prior_Loc]);

  fprintf(ERR, "optind=%d, argc=%d\n", optind,argc); 
  
  if( (argc<optind+nArgs-1)|| (argc>optind+nArgs+nOpt-1)) fprintf(ERR, "%s\n", usage), exit(1);
  dataFn=argv[optind];  
  printf("Reading data from ");

  if(!strcmp(dataFn, "-")) printf("stdin\n");
  else printf("%s\n", dataFn);

  printf("%s iterations\n", argv[optind+1]);
  nIter = strtol(argv[optind+1], &endptr, 10);

  seed = strtol(argv[optind+2], &endptr, 10);
  printf("Seed for random number generator: %i\n", seed);

  /* Optional parameters:*/
  /* Prefix */
  if(argc>=nArgs+1) {
    prefix=argv[optind+3];
    printf("Prefix prepended to output files: %s\n", prefix); 
  }

  /* Restart */
  if(argc==optind+nArgs+nOpt-1) {
    restart=1;
    re_K=argv[optind+4];
    re_P=argv[optind+5];
    printf("Restarting with parameters from:\n\t%s\n\t%s\n", re_K, re_P);
  }

  fprintf(ERR, "Setting alpha=%g, beta=%g\n", alpha, beta);
  setAlpha(alpha); setBeta(beta);

  fprintf(ERR, "Done.\n");
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
  double alpha=1, beta=1;

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

  fprintf(ERR, "Setting alpha=%g, beta=%g\n", alpha, beta);
  setAlpha(alpha); setBeta(beta);
}

#ifdef FIXED
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
#endif

/* Set int parameters to specific locations, estimate probabilities */
void initialiseLocations(intparameters *ip, 
			  const int *k0, const int *CDFdata, int nCDF,
			  parameters *p, int n) {
  int i;
  int succ, fail;

  for(i=0; i<n; i++) {
    int kPrev=(i-1>=0)?k0[i-1]:0;
    setIntParameter(ip, i, k0[i]);
    setIntProposal(ip, i, k0[i]);

    succ=succInterval(CDFdata, k0[i], kPrev);
    fail=k0[i]-kPrev-succ;
    setParameter(p, i, (double)succ/(succ+fail));
    setProposal(p, i, (double)succ/(succ+fail));
  }
  succ=succInterval(CDFdata, nCDF, k0[n-1]);
  fail=nCDF-k0[n-1]-succ;
  setParameter(p, n, (double)succ/(succ+fail));
  setProposal(p, n, (double)succ/(succ+fail));
}

/* Set int parameters to specific locations, estimate probabilities */
void initialiseLocationsAndProbs(intparameters *ip, parameters *p,
				 const int *k0, const double*p0,
				 int n) {
  int i;

  for(i=0; i<n; i++) {
    setIntParameter(ip, i, k0[i]);
    setIntProposal(ip, i, k0[i]);
    
    setParameter(p, i, p0[i]);
    setProposal(p, i, p0[i]);
  }

  setParameter(p, n, p0[n]);
  setProposal(p, n, p0[n]);
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

  if(!restart) {

    iterateRJ(prefix, prefix, prefix, prefix,
	      sampleBirthDeath,
	      /* dPriorChangepoints */prior,dLikelihood,
	      p0, nK0+1,1, nProb,
	      ip0, nK0,1, nK, seed,  nIter);
  }
  else {
  /* int k0[]={50000, 100000}; */
  /* nK0=2; */
  /* initialiseLocations(ip0, k0, CDFdata,getNdata()-1, */
  /* 		       p0, 2); */
  FILE *re_fp=fopen(re_K, "rb");
  char *buf=malloc(MAX_LEN+1);
  char *ll;
  int nReK=0, nReP=0;
  int i;
  
  int nums[NCHANGE+1];
  double dNums[NCHANGE+2];
  
  ll=lastLine(re_fp, buf, MAX_LEN);
  fclose(re_fp);
  
  nReK=str2ints(ll, nums);
  /* for(i=0; i<nReK; i++) { */
  /*   fprintf(OUT, "%d\t", nums[i]); */
  /* } */
  /* fprintf(OUT, "\nSuccessfully converted %d integers\n", nReK); */
  

  re_fp=fopen(re_P, "rb");
  ll=lastLine(re_fp, buf, MAX_LEN);
  free(buf);
  fclose(re_fp);
  
  nReP=str2doubles(ll, dNums);
  /* for(i=0; i<nReP; i++) { */
  /*   fprintf(OUT, "%f\t", dNums[i]); */
  /* } */
  /* fprintf(OUT, "\nSuccessfully converted %d doubles\n", nReP); */

  initialiseLocationsAndProbs(ip0, p0,
			      nums+1, dNums+1,
			      nReK-1);
  /* printIntParameters(OUT, ip0, nReK-1); */
  /* printParameters(OUT, p0, nReP-1); */
  /* exit(1); */

  iterateRJ(prefix, prefix, prefix, prefix,
	    sampleBirthDeath,
	    /* dPriorChangepoints */prior,dLikelihood,
	    p0, nReP-1,1, nProb,
	    ip0, nReK-1,1, nK, seed,  nIter);
  }
  return 0;  
}
 
