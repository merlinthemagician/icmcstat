/*
 * cp_prior.c
 *
 * Priors for number of changepoints, changepoint locations and
 * conjugate prior for open probabilities.
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
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

#include "mp_parameter.h"
#include "likelihood.h"
#include "cp_proposal.h"

#define OUT stdout
#define ERR stderr
 
/* Hyperparameters for conjugate prior */
static double alpha=1, beta = 1;

double getAlpha() {
  return alpha;
}

double getBeta() {
  return beta;
}

/* PRIOR */
/* returns changepoint m-1 or 0 or length of dataset. */
int getCPpar(const intparameters *ip, int m, int nIP) {
  if(m<=0) return 0;
  if(m>=nIP) return getNdata();
  return getIntParameter(ip, m-1);
}

/* returns changepoint m-1 or 0 or length of dataset. */
int getCPprop(const intparameters *ip, int m, int nIP) {
  if(m<=0) return 0;
  if(m>=nIP) return getNdata();
  return getIntProposal(ip, m-1);
}

/* Even-spaced order statistics as a prior for k*/
double logOrderStatistics(int km1, int k0, int kp1) {
  return log(kp1-k0) + log(k0-km1);
}


/* Scaling factor for prior for m changepoints in nData data
   points. The scaling factor is

   1 / Binomial(n, k) with n=nData-1, k = 2*m+1, see (Fearnhead, 2006)

   We know: (n+1)* B( n - k + 1, k + 1) = 1 / Binomial(n, k)

   where B is the Beta function. We use an implementation that
   calculates the logarithm 'logB' of the Beta function so that

   log( (n+1) * B( n - k + 1, k+1) ) 

   becomes 

   log( n+1 ) + logB( n-k+1, k+1).
     
*/
double logPriorScale(int nData, int m) {
  int status;
  gsl_sf_result result;
  int alpha=(nData-1) - (2*m+1) + 1;
  int beta=(2*m+1) + 1;
  
  status = gsl_sf_lnbeta_e ( alpha, beta, &result);
  if(status != GSL_SUCCESS) {
    fprintf(ERR, "Evaluation of beta function B(%i,%i) failed.\n",
	    alpha, beta);
    exit(1);
  }
  return  log(nData)+result.val;
}

/* Prior parameter for number of change points. lambda=3 is chosen for
   somewhat heuristic reasons by (Green, 1995). */
static int lambda=3;

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. */
double dPriorChangepoints(const parameters *p, 
				     int actP, int propP, int nP, 
				     const intparameters *ip, 
				     int actIP, int propIP, int nIP, 
				     likelihood *L) {
  int i;
  double logPar=0, logProp=0;
  int nData=getNdata();
  /* Prior for NUMBER of changepoints */
  double nPriorPar=gsl_ran_poisson_pdf(actIP, lambda);
  double nPriorProp=gsl_ran_poisson_pdf(propIP, lambda);

  /* Prior for number of changepoints */
  logPar+=nPriorPar;
  logProp+=nPriorProp;


  /* Prior for LOCATION. See bottom of page 207 in (Fearnhead,
     2006) */
  for(i=0; i<=actIP; i++){
    int k0=getCPpar(ip, i, actIP), k1=getCPpar(ip,i+1, actIP);
    logPar+=(k1>k0)?log(k1-k0):0;
  }

  for(i=0; i<=propIP; i++){
    int k0P=getCPprop(ip, i, propIP), k1P=getCPprop(ip, i+1, propIP);
    logProp+=(k1P>k0P)?log(k1P-k0P):0;
  }

  /* 'Divide' by scaling factors */
  logPar+=logPriorScale(nData, actIP);
  logProp+=logPriorScale(nData, propIP);

  setPrior(L, logPar, logProp);

  return logProp-logPar;
}


/* END: PRIOR */
