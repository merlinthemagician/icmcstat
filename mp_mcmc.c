/*
 *  mp_mcmc.c
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
#include <string.h>

#ifndef HPC
#include <gsl/gsl_rng.h>
#else
#include "gsl/gsl_rng.h"
#endif

#include "mp_parameter.h"
#include "mp_proposal.h"
#include "likelihood.h"

#define ERR stderr
#define OUT stdout

/* Global temperature for tempering */
static double globalTempPar=1, globalTempProp=1;

/* Set temperature for simulated tempering */
void setTemperature(double Tpar, double Tprop) {
  globalTempPar=Tpar; globalTempProp=Tprop;
}

double getTpar() {
  return globalTempPar;
}


double getTprop() {
  return globalTempProp;
}


/* Was the last move just accepted? */
static int justAccepted = 0;

/* Was the last move just accepted? */
int parameterUpdated() {
  return justAccepted;
}

/* sample positive parameters */
void samplePositive(const gsl_rng *r, 
		    parameters* p, int nPar, 
		    intparameters *ip, int nIPar) {
  if(p)  proposePositiveAB(r, p, 0, nPar);
  if(ip) proposeIntPositiveAB(r, ip, 0, nIPar);
}

/* sample parameters, constrained by max and min */
void sampleConstrained(const gsl_rng *r, 
		       parameters* p, int nPar,
		       intparameters *ip, int nIPar) {
  static int hello=1;
  if(hello) fprintf(OUT, "sampleConstrained()\n"), hello=0;
  if(p)  {proposeConstrainedAB(r, p, 0, nPar);}
  if(ip) proposeIntConstrainedAB(r, ip, 0, nIPar);
}

/* sample generation */
static void (*sample)(const gsl_rng *r, 
		      parameters *p, int,
		      intparameters *ip, int)=sampleConstrained;

/* Set new sample generator */
void setSample(void (*newSample)(const gsl_rng *r, 
		      parameters *p, int,
			      intparameters *ip, int)) {
  sample = newSample;
}


/* sample one randomly chosen parameter, constrained by max and min */
void sampleOneDoubleConstrained(const gsl_rng *r, 
				parameters* p, int nPar,
				intparameters *ip, int nIPar) {
  if(p)  {
    int i0=gsl_rng_uniform_int(r,nPar);
    proposeConstrainedAB(r, p, i0, 1);
  }
}

/* generic prior */
static double (*diffPrior)(const parameters *p, int nP,
			   const intparameters *ip, int nIP);

/* generic prior */
static double (*diffPosterior)(const parameters *p, int nP,
			       const intparameters *ip, int nIP);

/* outputs all parameters in file fp, starting from index k0 */
void generalOutputK0(FILE *fp, int iteration, 
		     const parameters *p, int k0, int nP) {
  fprintf(fp, "%i\t", iteration);
  printParametersK0(fp, p, k0, nP);
}

void generalOutput(FILE *fp, int iteration, 
		   const parameters *p, int nP) {
  generalOutputK0(fp, iteration, p, 0,  nP);
}

/* output of parameters */
static void (*output)(FILE *fp, int iteration,
		      const parameters *p, int nP)=generalOutput;


/* "Name" for iterations column*/
static const char *itheader="Iterations";

/*Output titles and initial values*/
void generalPrintTitles(FILE *p_fp, const parameters *p, int nP) {
  if(p_fp) { 
    fprintf(p_fp, "%s\t", itheader),printParameterNames(p_fp, p, nP);
  }
}

/* output of parameters */
static void (*printTitles)(FILE *fp, const parameters *p, int nP)=generalPrintTitles;

/* change output routine */
void setOutput(void (*newOutput)(FILE *fp, int iteration, 
				 const parameters *p, int nP)) {
  output=newOutput;
}

/* change output routine */
void setPrintTitles(void (*newPT)(FILE *fp,const parameters *p, int nP)) {
  printTitles=newPT;
}

/* compute likelihood of parameters */
void simpleDiffLikelihood(const parameters *p, int nP,
			  const intparameters *ip, int nIP,
			  likelihood *L) {
  double dPrior=0;
  double dPost=0;

  if(diffPrior) dPrior=diffPrior(p,nP, ip, nIP);
  if(diffPosterior) dPost=diffPosterior(p, nP, ip, nIP);
  /* fprintf(OUT, "simpleDiffLikelihood(): Setting likelihood"); */
  setDiffLikelihood(L, dPrior, dPost);
  /* fprintf(OUT, "simpleDiffLikelihood(): L ="); */
  /* printLikelihood(OUT, L); */
  /* exit(1); */
}

/* likelihood difference */
static void (*diffLikelihood)(const parameters *p, int nP,
			      const intparameters *ip, int nIP,
			      likelihood *L)=simpleDiffLikelihood;

/* accept proposal based upon likelihood L, simple Metropolis-Hastings */
int acceptProposal(const gsl_rng *r, const likelihood *L) {
  double diff=L->total[L->nTotal-1];
  int accept;
  if(diff >= 0) {accept = 1;}
  else {
    double Pproposal=exp(diff),samp;
    accept = ( (samp=gsl_rng_uniform(r)) < Pproposal);
    /* fprintf(OUT, "acceptProposal(): samp=%f, Pproposal = %f, accept = %i\n", */
    /* 	    samp, Pproposal, accept); */
  }

  return accept;
}

/* change parameters, compute likelihood and check if accepted */
int simpleMHstep(const gsl_rng *r, 
		 parameters *p, int nP,
		 intparameters *ip, int nIP,
		 likelihood *L){
  /* fprintf(OUT, "simpleMHstep(): Sample parameters...\n"); */
  sample(r, p, nP, ip, nIP);
  /* printFullParameters(OUT, p, nP); */
  /* printIntFullParameters(OUT, ip, nIP); */
  /* printIntProposal(OUT, ip, nIP); */
  /* exit(1); */
  /* fprintf(OUT, "simpleMHstep(): Likelihood...\n"); */
  diffLikelihood(p, nP, ip, nIP, L);
  /* fprintf(OUT, "simpleMHstep(): Computed likelihood.\n"); */
  /* printLikelihood(OUT, L); */
  return acceptProposal(r, L);
}

/* change parameters, compute likelihood and check if accepted */
int doubleMHstep(const gsl_rng *r, 
		 double (*dPrior)(const parameters *p, int nP, 
				  likelihood *L),
		 double (*dPosterior)(const parameters *p, int nP, 
				      likelihood *L),
		 parameters *p, int nP,
/* 		 intparameters *ip, int nIP, */
		 likelihood *L){

  sample(r, p, nP, NULL, 0);
  dPrior(p,nP,L);
  dPosterior(p,nP,L);
/*   diffLikelihood(p, nP, ip, nIP, L); */
  /* fprintf(OUT, "doubleMHstep(): Computed likelihood.\n"); */
  /* printLikelihood(OUT, L); */
  return acceptProposal(r, L);
}

/* change parameters, compute likelihood and check if accepted */
int MHstep(const gsl_rng *r, 
		 double (*dPrior)(const parameters *p, int nP, 
				  const intparameters *ip, int nIP, 
				  likelihood *L),
		 double (*dPosterior)(const parameters *p, int nP, 
				      const intparameters *ip, int nIP, 
				      likelihood *L),
		 parameters *p, int nP,
		 intparameters *ip, int nIP,
		 likelihood *L){

  /* sample both integer and parameter sets */
  sample(r, p, nP, ip, nIP);
  dPrior(p,nP,ip, nIP, L);
  dPosterior(p,nP,ip, nIP, L);
/*   diffLikelihood(p, nP, ip, nIP, L); */
  /* fprintf(OUT, "MHstep(): Computed likelihood.\n"); */
  /* printLikelihood(OUT, L); */
  return acceptProposal(r, L);
}

/* accept proposal based upon likelihood L, simple Metropolis-Hastings */
int acceptTemperedProposal(const gsl_rng *r, const likelihood *L, 
			   double Tpar, double Tprop) {
  double diff=L->total[L->nTotal-1];
  int accept;
  
  /* Normal move within chain (same temperature) */
  /* if(Tpar==Tprop) { */
  if(diff >= 0) {accept = 1;}
  else {
    double Pproposal=exp(diff/Tpar),samp;
    accept = ( (samp=gsl_rng_uniform(r)) < Pproposal);
    /* fprintf(OUT, "acceptTemperedProposal(): samp=%f, Pproposal = %f, accept = %i\n", */
    /* 	    samp, Pproposal, accept); */
  }
  /* Move between chains par <-> prop */
  if(accept && (Tpar !=Tprop )) {
    /* double  diff2=-diff/Tprop; */
    double alpha=exp(diff/Tpar+diff/Tprop), samp;
    /* Check if move to Tprop is accepted */
    if( (samp=gsl_rng_uniform(r)) < alpha) { 
      fprintf(OUT, "acceptTemperedProposal(): New temperature Told=%g, Tnew=%g\n", Tpar, Tprop);
      setTemperature(Tprop, Tprop);
      /* Tpar=Tprop; */
    } else { fprintf(OUT, "acceptTemperedProposal(): REJECTED: New temperature Told=%g, Tnew=%g\n", Tpar, Tprop);
      setTemperature(Tpar,Tpar);
/* Tprop=Tpar; */}
  }
    

  return accept;
}

/* Metropolis-Hastings sampler: output to files p_fp, ip_fp, diff
 * prior and posterior, initial values for p0 and i0,
 * iterations. Convention for prior/posterior: Both return the
 * difference and set values in likelihood which is given as an
 * argument.
 */
void iterateDiff(FILE *p_fp, FILE *ip_fp, FILE *like_fp,
		 double (*dPrior)(const parameters *p, int nP,
				  const intparameters *ip, int nIP),
		 double (*dPosterior)(const parameters *p, int nP,
				      const intparameters *ip, int nIP),
		 const parameters *p0, int nP0,
		 const intparameters *ip0, int nIP0,
		 int seed, int nIter) {
  parameters *p=NULL;
  intparameters *ip=NULL;
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
  likelihood *L=allocDiffLikelihood();

  int i=0, nAccepted=0;

  gsl_rng_set(r, seed);
  fprintf(OUT, "iterateDiff(): Random number generator initialised...\n");

  if(p0) p=malloc(nP0*sizeof(parameters));
  if(ip0) ip=malloc(nIP0*sizeof(intparameters));

  fprintf(OUT, "iterateDiff(): Memory allocated...\n");

  if(p0) {
    parameterMemcpy(p,p0, nP0);
    fprintf(OUT, "iterateDiff(): %i double parameters copied...\n", nP0);
  }

  if(ip0) {
    intParameterMemcpy(ip,ip0, nIP0);
    fprintf(OUT, "iterateDiff(): %i integer parameters copied...\n", nIP0);
  }

  /*set prior, posterior */
  diffPrior=dPrior;
  if(!diffPrior) fprintf(ERR, "iterateDiff(): Warning: No prior set!\n");
  diffPosterior=dPosterior;
  if(!diffPosterior) fprintf(ERR, "iterateDiff(): Warning: No posterior set!\n");

  fprintf(OUT, "iterateDiff(): Prior/posterior set...\n");

  /*Output titles and initial values*/
  if(p_fp) { 
    /* printParameterNames(OUT, p, nP0); */
    /* fprintf(OUT, "%s\n",p->names[0]); */
    /* exit(1); */
    fprintf(p_fp, "%s\t",itheader),printParameterNames(p_fp, p, nP0);
    fprintf(p_fp,"%i\t",i),printParameters(p_fp, p, nP0);
    fflush(p_fp);
  }
  if(ip_fp){
    fprintf(ip_fp, "%s\t", itheader), printIntParameterNames(ip_fp, ip, nIP0);
    fprintf(ip_fp,"%i\t",i),printIntParameters(ip_fp, ip, nIP0);
    fflush(ip_fp);
  }
  if (like_fp) { fprintf(like_fp, "%s\t", itheader),
      fprintf(like_fp, "diffPrior\tdiffPosterior\tdiffTotal\n");
    fflush(like_fp);
  }

  fprintf(OUT, "iterateDiff(): Initial conditions printed...\n");
  for(i=0; i<nIter; i++) {
    /* fprintf(OUT, "Step %i\n", i); */
    if(simpleMHstep(r, p, nP0, ip, nIP0, L)) {
      nAccepted++;
      if(like_fp) fprintf(like_fp, "%i\t", i),printLikelihood(like_fp, L);

      if(p) {
	copyProp2Org(p, nP0);
	fprintf(p_fp, "%i\t", i), printParameters(p_fp,p, nP0);
      }

      if(ip) {
	copyIntProp2Org(ip, nIP0);
	fprintf(ip_fp, "%i\t", i),printIntParameters(ip_fp,ip, nIP0);
      }
    }
  }
  fprintf(OUT, "iterateDiff(): Acceptance ratio %.2f %%\n", 
	  (float)nAccepted/nIter*100);
  gsl_rng_free(r);
}

/* copies values from proposal to original parameters */
static void (*copyDoublePar)(parameters *, int)= copyProp2Org;

/* sets routine for copying double parameters */
void setCopyPar(void (*copy)(parameters *, int)) {
  copyDoublePar=copy;
}

/* SHOULD BE PROGRAMMED AS A CALL TO THE MORE GENERAL iterateDiff() */
/* Metropolis-Hastings sampler: output to files p_fp, ip_fp, diff
   prior and posterior, initial values for p0 and i0, iterations */
void iterateDoubleDiff(FILE *p_fp, FILE *like_fp,
		       double (*dPrior)(const parameters *p, int nP,
				  const intparameters *ip, int nIP),
		       double (*dPosterior)(const parameters *p, int nP,
				  const intparameters *ip, int nIP),
		       parameters *p0, int nP0,
		       int seed, int nIter) {
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
  likelihood *L=allocDiffLikelihood();

  int i=0, nAccepted=0;

  gsl_rng_set(r, seed);
  fprintf(OUT, "iterateDoubleDiff(): Random number generator initialised...\n");

  /*set prior, posterior */
  diffPrior=dPrior;
  if(!diffPrior) fprintf(ERR, "iterateDoubleDiff(): Warning: No prior set!\n");
  diffPosterior=dPosterior;
  if(!diffPosterior) fprintf(ERR, "iterateDoubleDiff(): Warning: No posterior set!\n");

  fprintf(OUT, "iterateDoubleDiff(): Prior/posterior set...\n");

  /*Output titles and initial values*/
  if(p_fp) { 
    fprintf(p_fp, "%s\t", itheader),printParameterNames(p_fp, p0, nP0);
    fprintf(p_fp,"%i\t",i),printParameters(p_fp, p0, nP0);
    fflush(p_fp);
  }
  if (like_fp) { fprintf(like_fp, "%s\t", itheader),
      fprintf(like_fp, "diffPrior\tdiffPosterior\tdiffTotal\n");
    fflush(like_fp);
  }

  fprintf(OUT, "iterateDoubleDiff(): Initial conditions printed...\n");
  for(i=0; i<nIter; i++) {
    justAccepted=1;
    if(simpleMHstep(r, p0, nP0, NULL, 0, L)) {
      justAccepted=0;
      nAccepted++;
      if(like_fp) fprintf(like_fp, "%i\t", i),printLikelihood(like_fp, L);
      fflush(like_fp);

      if(p0) {
	copyDoublePar(p0,nP0);
	/* copyProp2Org(p0, nP0); */
	/* fprintf(p_fp, "%i\t", i), */ output(p_fp,i,p0, nP0), fflush(p_fp);
      }
    }
  }
  fprintf(OUT, "iterateDoubleDiff(): Acceptance ratio %.2f %%\n", 
	  (float)nAccepted/nIter*100);
  gsl_rng_free(r);
}

/* static const char* likelihoodTitles="oldPrior\tnewPrior\tdiffPrior\toldLikelihood\tnewLikelihood\tdiffLikelihood\toldPosterior\tnewPosterior\tdiffPosterior\n"; */

#define LTITLES "oldPrior\tnewPrior\tdiffPrior\toldLikelihood\tnewLikelihood\tdiffLikelihood\toldPosterior\tnewPosterior\tdiffPosterior\n"

/* SHOULD BE PROGRAMMED AS A CALL TO THE MORE GENERAL iterateDiff() */
/* Metropolis-Hastings sampler: output to files p_fp, prior and
   posterior, initial values for p0, iterations */
void iterateDouble(FILE *p_fp, FILE *like_fp,
		   double (*dPrior)(const parameters *p, int nP, 
				    likelihood *L),
		   double (*dPosterior)(const parameters *p, int nP, 
					likelihood *L),
		   parameters *p0, int nP0,
		   int seed, int nIter) {
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
  likelihood *L=allocLikelihood();

  int i=0, nAccepted=0;

  gsl_rng_set(r, seed);
  fprintf(OUT, "iterateDouble(): Random number generator initialised...\n");

  printTitles(p_fp, p0, nP0);
  output(p_fp,i,p0, nP0), fflush(p_fp);

  if (like_fp) { fprintf(like_fp, "%s\t", itheader),
      fprintf(like_fp, LTITLES);
    fflush(like_fp);
  }


  fprintf(OUT, "iterateDouble(): Initial conditions printed...\n");
  for(i=0; i<nIter; i++) {
    const int printDot=/* ((i*80)/nIter); */nIter/80;
    
    if(i%printDot==0) fprintf(OUT, "."), fflush(OUT);

    justAccepted=1;
    sample(r, p0, nP0, NULL, 0);
    dPrior(p0,nP0,L);
    dPosterior(p0,nP0,L);
    /* if(doubleMHstep(r, dPrior, dPosterior, p0, nP0, L)) { */
    if(acceptTemperedProposal(r,L, getTpar(), getTprop())) {
      /* fprintf(OUT, "iterateDouble(): accepted. Told=%g, Tnew=%g\n", getTpar(), getTprop()); */
      justAccepted=0;
      nAccepted++;
      if(like_fp && (getTpar()==1)) fprintf(like_fp, "%i\t", i),printLikelihood(like_fp, L);
      fflush(like_fp);

      if(p0) {
	/* fprintf(OUT, "iterateDouble(): Accepted!\n"); */
	/* printFullParameters(OUT, p0,nP0); */
	copyDoublePar(p0,nP0);
	/* copyProp2Org(p0, nP0); */
	if(getTpar()==1) {
	/* fprintf(p_fp, "%i\t", i), */ output(p_fp,i,p0, nP0), fflush(p_fp);
	}
      }

    }
  }
  fprintf(OUT, "iterateDouble(): Acceptance ratio %.2f %%\n", 
	  (float)nAccepted/nIter*100);
  gsl_rng_free(r);
}

static int tWalkSampledFirst=0;

/* Returns the index of the sampled parameter */
int getTWalkSampled() {
  return tWalkSampledFirst;
}

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
			int seed, int nIter) {
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
  likelihood *L=allocLikelihood();

  int i=0, nAccepted=0;

  gsl_rng_set(r, seed);
  fprintf(OUT, "iterateDouble(): Random number generator initialised...\n");

  /*Output titles and initial values*/
  if(p_fp) { 
    fprintf(p_fp, "%s\t", itheader),printParameterNames(p_fp, p1, nP1);
    fprintf(p_fp,"%i\t",i),printParameters(p_fp, p1, nP1);
    fflush(p_fp);
  }
  if (like_fp) { fprintf(like_fp, "%s\t", itheader),
      fprintf(like_fp, LTITLES);
    fflush(like_fp);
  }

  fprintf(OUT, "iterateDouble(): Initial conditions printed...\n");
  for(i=0; i<nIter; i++) {
    int sampledFirst=sample(r, p1, nP1, p2, nP2);
    parameters *sampledP=sampledFirst?p1:p2;
    int nP=sampledFirst?nP1:nP2;
    double propDensity=getLogTwalkScale();
    const int printDot=/* ((i*80)/nIter); */nIter/80;
    
    /* fprintf(OUT, "iterateDoubleTwalk(): iteration %i, printDot=%i (%i)\n", i, */
    /* 	    printDot, printDot%80); */
    if(i%printDot==0) fprintf(OUT, "."), fflush(OUT);

    tWalkSampledFirst=sampledFirst;
    justAccepted=1;
    prior(sampledP, nP, L);	
    likely(sampledP, nP, L);
    
    addToLikelihood(L, propDensity);

    /* if(acceptProposal(r,L)) { */
    if(acceptTemperedProposal(r,L, getTpar(), getTprop())) {
      /* if(doubleMHstep(r, dPrior, dPosterior, p0, nP0, L)) { */
      /* fprintf(OUT, "iterateDoubleTwalk(): accepted. Told=%g, Tnew=%g\n", getTpar(), getTprop()); */
      justAccepted=0;
      /* nAccepted++; */
      if(like_fp && (getTpar()==1)) {
	fprintf(like_fp, "%i\t", i),printLikelihood(like_fp, L);
	fflush(like_fp);
      } 

      copyDoublePar(sampledP,nP);
      /* if(p0) { */
      /* copyProp2Org(p0, nP0); */
      /* fprintf(p_fp, "%i\t", i), */ /* output(p_fp,i,p0, nP0), fflush(p_fp); */
      if(getTpar()==1) nAccepted++, output(p_fp,i,sampledP, nP), fflush(p_fp);
      /* } */

    }
  }
  fprintf(OUT,"\n");
  fprintf(OUT, "iterateDouble(): Acceptance ratio %.2f %%\n", 
	  (float)nAccepted/nIter*100);
  gsl_rng_free(r);
}

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
		 int seed, int nIter) {
  parameters *p=NULL;
  intparameters *ip=NULL;
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
  likelihood *L=allocLikelihood();

  int i=0, nAccepted=0;

  gsl_rng_set(r, seed);
  fprintf(OUT, "iterate(): Random number generator initialised...\n");

  if(p0) p=malloc(nP0*sizeof(parameters));
  if(ip0) ip=malloc(nIP0*sizeof(intparameters));

  fprintf(OUT, "iterate(): Memory allocated...\n");

  if(p0) {
    parameterMemcpy(p,p0, nP0);
    fprintf(OUT, "iterate(): %i double parameters copied...\n", nP0);
  }

  if(ip0) {
    intParameterMemcpy(ip,ip0, nIP0);
    fprintf(OUT, "iterate(): %i integer parameters copied...\n", nIP0);
  }

  /*set prior, posterior */
  /* NOT HERE! */

  /*Output titles and initial values*/
  if(p_fp) { 
    /* printParameterNames(OUT, p, nP0); */
    /* fprintf(OUT, "%s\n",p->names[0]); */
    /* exit(1); */
    fprintf(p_fp, "%s\t",itheader),printParameterNames(p_fp, p, nP0);
    fprintf(p_fp,"%i\t",i),printParameters(p_fp, p, nP0);
    fflush(p_fp);
  }
  if(ip_fp){
    fprintf(ip_fp, "%s\t", itheader), printIntParameterNames(ip_fp, ip, nIP0);
    fprintf(ip_fp,"%i\t",i),printIntParameters(ip_fp, ip, nIP0);
    fflush(ip_fp);
  }
  if (like_fp) { fprintf(like_fp, "%s\t", itheader),
      fprintf(like_fp, LTITLES);
    fflush(like_fp);
  }

  fprintf(OUT, "iterateDiff(): Initial conditions printed...\n");
  for(i=0; i<nIter; i++) {
    /* fprintf(OUT, "Step %i\n", i); */
    if(MHstep(r, dPrior, dPosterior,p,  nP0, ip, nIP0, L)) {
      nAccepted++;
      if(like_fp) fprintf(like_fp, "%i\t", i),printLikelihood(like_fp, L);

      if(p) {
	copyProp2Org(p, nP0);
	fprintf(p_fp, "%i\t", i), printParameters(p_fp,p, nP0);
      }

      if(ip) {
	copyIntProp2Org(ip, nIP0);
	fprintf(ip_fp, "%i\t", i),printIntParameters(ip_fp,ip, nIP0);
      }
    }
  }
  fprintf(OUT, "iterateDiff(): Acceptance ratio %.2f %%\n", 
	  (float)nAccepted/nIter*100);
  gsl_rng_free(r);
}

static char *allocFn(const char *pref, const char*mid, const char * suf) {
  return malloc(strlen(pref)+strlen(mid)+strlen(suf)+1);
}

/* Metropolis-Hastings sampler: output to files p_fp, ip_fp, diff
 * prior and posterior, initial values for p0 and i0, iterations. kFn
 * is the file name where "dimensions" of model will be saved.
 * Convention for prior/posterior: Both return the difference and set
 * values in likelihood which is given as an argument.
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
	       int seed, int nIter) {
  parameters *p=NULL;
  intparameters *ip=NULL;
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
  likelihood *L=allocLikelihood();

  int i=0, nAccepted=0;
  int actP=actP0,actIP=actIP0;
  const char *suf=".dat";
  const char *pMid="_p000";
  const char *ipMid="_ip000";
  const char *likeMid="_likelihood000";
  const char *kMid="_nK";
  char * pFmt="%s_p%03i.dat";
  char * ipFmt="%s_ip%03i.dat";
  char * pfn=allocFn(p_prefix, pMid, suf);
  char * ipfn=allocFn(ip_prefix, ipMid, suf);
  char * likefn=allocFn(like_prefix, likeMid, suf);
  char *kfn=allocFn(k_prefix, kMid, suf);
  FILE *like_fp, *k_fp;
  sprintf(likefn,"%s%s",like_prefix,"_likelihood.dat");
  like_fp=fopen(likefn,"w");

  sprintf(kfn,"%s%s",k_prefix,"_nK.dat");
  k_fp=fopen(kfn,"w");

  /*Initialise format strings for file names*/
  gsl_rng_set(r, seed);
  fprintf(OUT, "iterateRJ(): Random number generator initialised...\n");

  if(p0) p=malloc(nP0*sizeof(parameters));
  if(ip0) ip=malloc(nIP0*sizeof(intparameters));

  fprintf(OUT, "iterateRJ(): Memory allocated...\n");

  if(p0) {
    parameterMemcpy(p,p0, nP0);
    fprintf(OUT, "iterateRJ(): %i double parameters copied...\n", nP0);
  }

  if(ip0) {
    intParameterMemcpy(ip,ip0, nIP0);
    fprintf(OUT, "iterateRJ(): %i integer parameters copied...\n", nIP0);
  }

  /* Titles and initial values for actP and actIP */
  if(k_fp) {
    fprintf(k_fp, "%s",itheader);
    if(p) fprintf(k_fp, "\t%s", "\"n_K (d)\"");
    if(ip) fprintf(k_fp, "\t%s", "\"n_K\"");
    fprintf(k_fp, "\n");
    fprintf(k_fp, "%i", 0);
    if(p) fprintf(k_fp, "\t%i", actP0);
    if(ip) fprintf(k_fp, "\t%i", actIP0);
    fprintf(k_fp, "\n");
    fflush(k_fp);
  }

  /*Output titles and initial values*/
  if(p_prefix) { 
    FILE *p_fp;
    int k;
    /* fprintf(ERR, "minP = %i, nP0 = %i\n", minP, nP0); */
    /* Initialise from minP up to nP0: nP0+nP0-minIP */
    for(k=minP;k<nP0+1; k++) {
      sprintf(pfn,pFmt,p_prefix,k);
      /* fprintf(OUT, "%s\n", pfn); */
      p_fp=fopen(pfn, "w");
      fprintf(p_fp, "%s\t",itheader),printParameterNames(p_fp, p, k);
      fflush(p_fp);
      fclose(p_fp);
    }
    /* fprintf(p_fp,"%i\t",i),printParameters(p_fp, p, actP); */
  }
  if(ip_prefix){
    FILE *ip_fp;
    int k;
    /* Initialise from minIP up to nIP0: nIP0+nIP0-minIP */
    /* fprintf(ERR, "minIP = %i, nIP0 = %i\n", minIP, nIP0); */
    for(k=minIP; k<nIP0+1; k++) {
      sprintf(ipfn,ipFmt,ip_prefix,k);
      /* fprintf(OUT, "%s\n", ipfn); */
      ip_fp=fopen(ipfn, "w");
      fprintf(ip_fp, "%s\t", itheader), printIntParameterNames(ip_fp, ip, k);
      fflush(ip_fp), fclose(ip_fp);
    }
    /* fprintf(ip_fp,"%i\t",i),printIntParameters(ip_fp, ip, actIP); */
  }
  if (like_fp) {
    fprintf(like_fp, "%s\t", itheader),fprintf(like_fp, LTITLES);
    fflush(like_fp);
  }

  fprintf(OUT, "iterateRJ(): Initial conditions printed...\n");
  for(i=0; i<nIter; i++) {
    int propP=actP, propIP=actIP;
    /* fprintf(OUT, "Step %i, actP=%i, propP=%i\n", i, actP, propP); */
    const int printDot=/* ((i*80)/nIter); */nIter/75;
    
    /* if(i%printDot==0) fprintf(OUT, "."), fflush(OUT); */

    /* Next step */
    sample(r, p, actP, &propP, nP0, ip, actIP, &propIP,  nIP0);
    dPrior(p, actP, propP, nP0, ip, actIP, propIP, nIP0, L);
    dLikelihood(p, actP, propP, nP0, ip, actIP, propIP, nIP0, L);
    
    /* if(acceptProposal(r,L)) { */
    if(acceptTemperedProposal(r,L, globalTempPar, globalTempProp)) {
      /* fprintf(OUT, "iterateRJ(): accepted!\n"); */
      nAccepted++;
      actP=propP, actIP=propIP;
      /* fprintf(OUT, "iterateRJ(): Likelihood printed\n"); */

      /* if conditional is unnecessary */
      /* if(globalTempPar != globalTempProp) { globalTempPar=globalTempProp;} */

      /* only write to file if temperature is 1 */
      if(like_fp && (globalTempPar == 1) ) {
	if(i>0)
	fprintf(like_fp, "%i\t", i),printLikelihood(like_fp, L);
      }

      if(k_fp && (globalTempPar == 1) ) {
	fprintf(k_fp, "%i", i);
	if(p) fprintf(k_fp, "\t%i", actP);
	if(ip) fprintf(k_fp, "\t%i", actIP);
	fprintf(k_fp, "\n");
	fflush(k_fp);
      }

      if(p) {
	FILE *p_fp;
	/* Copy if accepted ... */
	copyProp2Org(p, actP);
	/* ... but only write to output if this number of parameters
	   exist */
	if((globalTempPar==1) && (actP>=minP) && (actP<=nP0)) {
	  sprintf(pfn, pFmt, p_prefix, actP);
	  p_fp=fopen(pfn,"a");
	  /* fprintf(OUT, "iterateRJ(): Copying %i parameters\n", actP); */
	  fprintf(p_fp, "%i\t", i), printParameters(p_fp,p, actP);
	  fflush(p_fp), fclose(p_fp);
	}
      }
      
      /* fprintf(OUT, "iterateRJ(): Parameters printed\n"); */

      if(ip) {
	FILE *ip_fp;
	/* Copy if accepted ... */
	copyIntProp2Org(ip, nIP0);
	/* ... but only write to output if this number of parameters
	   exist */
	if((globalTempPar==1) && (actIP>=minIP) && (actIP<=nIP0)) {
	  sprintf(ipfn, ipFmt, ip_prefix, actIP);
	  ip_fp=fopen(ipfn,"a");

	  fprintf(ip_fp, "%i\t", i),printIntParameters(ip_fp,ip, actIP);
	  /* fprintf(OUT, "%i\t", i),printIntParameters(OUT,ip, actIP); */
	  fflush(ip_fp), fclose(ip_fp);
	}
      }
      
    }
    /* else { */
    /*   globalTempProp=globalTempPar; */
    /* } */
  }
  fprintf(OUT, "iterateRJ(): Acceptance ratio %.2f %%\n", 
	  (float)nAccepted/nIter*100);
  gsl_rng_free(r);
}
