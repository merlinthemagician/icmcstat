/*
 * bernoulliOneChangePointMCMC.c
 *
 * MCMC method for calculating the distribution of a change point k
 * for a bernoulli-distributed data set whose parameter p changes
 * instantly.
 *
 * Ivo Siekmann, 2/11/2010
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
#include <gsl/gsl_sf.h>

#include <limits.h>

#include "bernoulliDist.h"
#include "cp_prior.h"
#include "mp_mcmc.h"
#include "mp_parameter.h"
#include "mp_proposal.h"

#define OUT stdout
#define ERR stderr

/* Number of successes observed between k0 and k1 */
int succInterval(const int *CDF, int k0, int k1) {
  int succIn1=(k0<=0)?CDF[k1-1]:CDF[k1-1]-CDF[k0-1];
  return succIn1;
}

/* Number of failures observed between k0 and k1 */
int failInterval(const int *CDF, int k0, int k1) {
  return k1-k0-succInterval(CDF, k0, k1);
}

/* Bernoulli-distributed data with change point */
static int *CDFdata=NULL, nCDFdata=100;

/* Set number of data points in trace */
void setNdata(int n) {
  nCDFdata=n;
}

/* Get number of data points in trace */
int getNdata() {
  return nCDFdata;
}

/* Set number of data points in trace */
const int *getData() {
  return CDFdata;
}

/* Get number of data points in trace */
void setData(int *CDF, int n) {
  CDFdata=CDF;
  nCDFdata=n;
}


/* For birth/death move */
static double proposalScale=1;

/* Get scaling factor for birth/death move */
double getProposalScale() {
  return proposalScale;
}

void setProposalScale(double p) {
  proposalScale=p;
}

/* Probability of moving changepoint location or open probability */
static double mixingRatio =0.5;


/* Sample probabilities according to successes and failures */
void sampleProbabilities(const gsl_rng *r,
			 const int *CDFdata, int nCDFdata,
			 parameters *p, int nPar, const intparameters *ip, int nIPar) {
  int k;
  for(k=0; k<nPar; k++) {
    int k0=(k<=0)?0:getIntProposal(ip,k-1);
    int k1=(k<nIPar)?getIntProposal(ip,k):nCDFdata-1;
    int success=succInterval(CDFdata, k0, k1)/* (k0<=0)?CDFdata[k1-1]:CDFdata[k1-1]-CDFdata[k0-1] */;
    int failures=k1-k0-success;
    double estP;
    double alpha=getAlpha(), beta=getBeta();

    estP=gsl_ran_beta(r, success+alpha+1, failures+beta+1);
    /* fprintf(OUT, "p=%g\n", estP); */
    setProposal(p, k, estP);
  }
}

/* Samples either one change point or one probability */
void sampleMixed(const gsl_rng *r,
		 parameters* p, int nPar,
		 intparameters *ip, int nIPar) {
  static int hello=1;
  double pr=gsl_rng_uniform(r);
  int moveK=pr<mixingRatio;
  const int *CDFdata=getData(), nCDFdata=getNdata();

  if(hello) fprintf(OUT, "sampleMixed()...\n"), hello=0;

  /* Set proposed changepoint locations and probabilities to current
     values.*/
  resetProposals(p, nPar, ip, nIPar);

  /* Assuming there is NO changepoint */
  if(nIPar<=0) moveK=0;
  if(moveK) {
    double pJump=0.1;
    /* Which k shall be moved? */
    int indK=gsl_rng_uniform_int(r, nIPar);
    int kDown=(indK-1>=0)?getIntParameter(ip,indK-1):0;
    int kUp=(indK+1<nIPar)?getIntParameter(ip,indK+1):getNdata();
    int kNew;
    int succ0, succ1;
    int fail0, fail1;
    double p0, p1;


    setIntMin(ip, indK, kDown);
    setIntMax(ip, indK, kUp);

    if(gsl_rng_uniform(r)>=pJump) {
      /*Move only one parameter by width */
      proposeIntConstrainedAB(r, ip, indK, 1);
    } else {
      proposeIntShiftInterval(r, ip, indK, 1);
    }
      /* Update probabilities */
      kNew=getIntProposal(ip,indK);

      succ0=succInterval(CDFdata, kDown, kNew);
      succ1=succInterval(CDFdata, kNew, kUp);
      fail0=kNew-kDown-succ0;
      fail1=kUp-kNew-succ1;
      p0=gsl_ran_beta(r, succ0, fail0);
      p1=gsl_ran_beta(r, succ1,fail1);

      setProposal(p, indK, p0);
      setProposal(p, indK+1, p1);
  }
  else {
    sampleProbabilities(r, CDFdata, nCDFdata, p, nPar, ip, nIPar);
  }
}


/* Selects changepoint position and adds new changepoint before or after */
void sampleBirth(const gsl_rng *r,
		 parameters* p, int actP, int nPar,
		 intparameters *ip, int actIP, int nIPar) {
  /* select parameter */
  int kIndex=gsl_rng_uniform_int(r, actIP+1);
  int kDown=(kIndex<=0)?0:getIntParameter(ip,kIndex-1);
  int kUp=(kIndex==actIP)?getNdata()-1:getIntParameter(ip,kIndex)-1;

  /*sample new change point*/
  int kNew;
  int i;
  static int hello=1;
  const int *CDF=getData(), nCDF=getNdata();
  int succ0, succ1, fail0, fail1;
  double p0, p1;
  double alpha=getAlpha(), beta=getBeta();

  /* Adding the first changepoint */
  if(actIP==0) {
    /* k must be between 1 and nCDF-1*/
    int kNew=gsl_rng_uniform_int(r,nCDF-1)+1;

    /* Successes and failures before and after new changepoint kNew
       are used for sampling p0 and p1... */
    succ0=succInterval(CDF, 0, kNew), fail0=failInterval(CDF, 0, kNew);
    succ1=succInterval(CDF, kNew, nCDF), fail1=failInterval(CDF, kNew, nCDF);

    /* ... from beta distributions. */
    p0=gsl_ran_beta(r, succ0+1+alpha, fail0+beta+1);
    p1=gsl_ran_beta(r, succ1+1+alpha, fail1+beta+1);

    /* Set new proposal */
    setIntProposal(ip, 0, kNew);
    setProposal(p,0,p0);
    setProposal(p,1,p1);

    /* Balance with corresponding death move */
    proposalScale=1/((double)(nCDF-1));
    return;
  }

  /* Find kUp and kDown that are more than one data point apart. */
  while(kUp-kDown<=1) {
    kIndex=gsl_rng_uniform_int(r, actIP+1);
    kDown=(kIndex<=0)?0:getIntParameter(ip,kIndex-1);
    kUp=(kIndex==actIP)?getNdata()-1:getIntParameter(ip,kIndex)-1;
  }

  /* Balance with corresponding death move */
  proposalScale=1/((double)(kUp-kDown));

  do {
    kNew=kDown+gsl_rng_uniform_int(r, kUp-kDown);
  }while (kNew==0);

  if(hello) fprintf(OUT, "sampleBirth()\n"), hello=0;

  /*Copy kIndex parameters unchanged */
  copyIntOrg2Prop(ip, kIndex);

  /* Insert new changepoint */
  setIntProposal(ip, kIndex, kNew);

  /* Shift changepoints after inserted changepoint. */
  for(i=kIndex+1; i<nIPar; i++) {
    setIntProposal(ip, i, getIntParameter(ip,i-1));
  }

  /* Now copy kIndex open probabilities. */
  if(kIndex-1>0) copyOrg2Prop(p, kIndex);

  /* Sample open probabilities for the segments 'left' and 'right' of
     the new changepoint. */
  succ0=succInterval(CDF, kDown, kNew), fail0=failInterval(CDF, kDown, kNew);
  succ1=succInterval(CDF, kNew, kUp), fail1=failInterval(CDF, kNew, kUp);
  p0=gsl_ran_beta(r, succ0+alpha+1, fail0+beta+1);
  p1=gsl_ran_beta(r, succ1+alpha+1, fail1+beta+1);

  setProposal(p, kIndex, p0);
  setProposal(p, kIndex+1, p1);

  /* Shift the remaining open probabilities. */
  for(i=kIndex+2; i<nPar; i++) {
    double pI=getParameter(p,i-1);
    setProposal(p, i, pI);
  }
}

/* Selects changepoint position and adds new changepoint before or after */
void sampleDeath(const gsl_rng *r,
		 parameters* p, int actP, int nPar,
		 intparameters *ip, int actIP, int nIPar) {
  /* select change point that will be deleted */
  int kIndex=gsl_rng_uniform_int(r, actIP);
  int kDown, kUp;
  int succ, fail;
  /* Open probabilities that replaces probabilities 'left' and 'right'
     of deleted changepoint.*/
  double pNew;
  int i;
  static int hello=1;
  const int*CDF=getData(), nCDF=getNdata();
  double alpha=getAlpha(), beta=getBeta();

  if(hello) fprintf(OUT, "sampleDeath()\n"), hello=0;

  /*Copy kIndex parameters unchanged */
  copyIntOrg2Prop(ip, kIndex);

  /* Omit kIndex */
  for(i=kIndex; i<nIPar-1; i++) {
    setIntProposal(ip, i, getIntParameter(ip,i+1));
  }

  /*Copy kIndex double parameters unchanged */
  copyOrg2Prop(p, kIndex);

  kDown=(kIndex<=0)?0:getIntParameter(ip, kIndex-1);
  kUp=(kIndex==actIP)?nCDF-1:getIntParameter(ip,kIndex)-1;
  succ=succInterval(CDF, kDown, kUp), fail=failInterval(CDF, kDown, kUp);
  pNew=gsl_ran_beta(r, succ+alpha+1,fail+beta+1);
  setProposal(p, kIndex, pNew);

  /* Balance with corresponding birth move */
  proposalScale=(double)kUp-kDown;

  for(i=kIndex+1; i<nPar-1; i++) {
    setProposal(p, i, getParameter(p,i+1));
  }
}

/* Birth/Death proposal or moving changepoints. */
void sampleBirthDeath(const gsl_rng *r,
		      parameters *p, int actP,
		      int *propP, int nPar,
		      intparameters *ip, int actIP,
		      int *propIP, int nIPar) {

  /* Probability for making a birth or a death move. */
  double birthDeath=0.5;

  if(gsl_rng_uniform(r) < birthDeath) {
    /* If no changepoints, only birth move is possible */
    if(actIP == 0) {
      sampleBirth(r, p, actP, nPar, ip,actIP, nIPar);
      *propIP=actIP+1;
      *propP=actP+1;
    }
    /* If maximum number of changepoints is reached, only death move
       is possible. */
    else if (actIP==nIPar-1){
      sampleDeath(r, p, actP, nPar, ip, actIP, nIPar);
      *propIP=actIP-1;
      *propP=actP-1;
    }
    else {
      /* Choose birth or death with 50 % */
      if(gsl_rng_uniform(r) < 0.5) {
	sampleBirth(r, p, actP, nPar, ip,actIP, nIPar);
	*propIP=actIP+1;
	*propP=actP+1;
      }
      else {
	sampleDeath(r, p, actP, nPar, ip, actIP, nIPar);
	*propIP=actIP-1;
	*propP=actP-1;
      }
    }
  }
  /* 'Normal move', no birth/death */
  else {
    sampleMixed(r, p, actP, ip, actIP);setProposalScale(1);
  }
}
