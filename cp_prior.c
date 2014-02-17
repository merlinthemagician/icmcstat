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

/* Routines for Fearnhead's suggested negative binomial process prior,
   see Fearnhead (2006) */

/* Calculates binomial coefficient Binomial(n, k) from beta function. */
double binomialFromBeta(double n, double k) {
  int status;
  gsl_sf_result result;
  /* fprintf(ERR, "Calculating B(%i,%i).\n", */
  /* 	  (int)(n-k+1), (int)k+1); */
  if(n-k+1<k+1) {fprintf(ERR, "Setting B(%i,%i)=0.\n", 
			     (int)(n-k+1), (int)k+1);return 0;}
  status = gsl_sf_beta_e ( n-k+1, k+1, &result);
  if(status != GSL_SUCCESS) {
    fprintf(ERR, "Evaluation of beta function B(%i,%i) failed.\n",
	    (int)(n-k+1), (int)k+1);
    exit(1);
  }
  return 1.0/result.val/(n+1);
}

static int nb_k=1;
static double nb_p=0.001;

/* g according to Fearnhead (2006) */
double negativeBin_g(int t, int k) {
  double pT=pow(nb_p,t), qTminK=pow(1-nb_p,t-k);
  double result=binomialFromBeta(t-k, k-1)*pT*qTminK;
  /* fprintf(ERR, "negativeBin_g(): %g\n", result); */
  return result;
}

/* g according to Fearnhead (2006) */
double negativeBin_g0(int t) {
  int i;
  double sum=0;
  for(i=1; i<=nb_k; i++) {
    sum+=negativeBin_g(t, i);
  }
  return sum/nb_k;
}

/* G according to Fearnhead (2006) */
double negativeBin_G(int t)  {
  int s;
  double sum=0;
  for(s=nb_k; s<=t; s++) {
    sum+=negativeBin_g(s, nb_k);
  }
  return sum;
}

/* Negative binomial log-prior ratio according to Fearnhead (2006) */
double logNegativeBinomialPriorDiff(const intparameters *ip, 
				    int actIP, int propIP, int nIP) {
  int i;
  double logPar=0, logProp=0;

  int j0,j1;

  if(actIP > 0 ) {
    /* First changepoint */
    j0=getCPpar(ip, 1, actIP);
    fprintf(ERR, "Par:\n j1=%i\n", j0);

    logPar += log(negativeBin_g0(j0));
    fprintf(ERR, "log-g0=%g\n", logPar);
    
    for(i=2;i<=actIP; i++) {
      double factor;
      j1=getCPpar(ip, i, actIP);
      factor=log(negativeBin_g(j1-j0, nb_k));
      fprintf(ERR, "g(%i)=%g\n", j1, factor);
      logPar+=log(factor);
      j0=j1;
    }

    /* Last changepoint */
    j0=j1;
    j1=getNdata();
    logPar+=log(1-negativeBin_G(j1-j0));
  }

  if(propIP>0) {
    /* First changepoint */
    j0=getCPprop(ip, 1, propIP);
    logProp += log(negativeBin_g0(j0));

    for(i=2;i<=propIP; i++) {
      j1=getCPprop(ip, i, propIP);
      logProp+=log(negativeBin_g(j1-j0, nb_k));
      j0=j1;
    }

    /* Last changepoint */
    j1=getNdata();
    logProp+=log(1-negativeBin_G(j1-j0));
  }
  
  fprintf(ERR, "logNegativeBinomialPriorDiff(): prop=%g, par=%g, diff=%g\n", logProp,logPar, logProp-logPar);
  return logProp-logPar;
}

/* log-Prior for a fixed number of consecutive changepoints calculated
   for current sample and proposal. The conjugate prior for the
   probabilities is accounted for in the likelihood. */
double dPriorFixedChangepoints(const parameters *p, int nP, 
			  const intparameters *ip, int nIP, 
			  likelihood *L) {
  int i;
  double logPar=0, logProp=0;
  int nData=getNdata();

  /* Prior for LOCATION. See bottom of page 207 in (Fearnhead,
     2006) */
  for(i=0; i<=nIP; i++){
    int k0=getCPpar(ip, i, nIP), k1=getCPpar(ip,i+1, nIP);
    logPar+=(k1>k0)?log(k1-k0):0;
  }

  for(i=0; i<=nIP; i++){
    int k0P=getCPprop(ip, i, nIP), k1P=getCPprop(ip, i+1, nIP);
    logProp+=(k1P>k0P)?log(k1P-k0P):0;
  }

  /* 'Divide' by scaling factors */
  logPar+=logPriorScale(nData, nIP);
  logProp+=logPriorScale(nData, nIP);

  setPrior(L, logPar, logProp);

  return logProp-logPar;
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

/* Geometric prior for number of change points. */
static double geoP=0.001;

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. Uses a geometric
   prior for the number of changepoints whereas the prior for the
   locations is based upon the uniform order statistics.*/
double dPriorChangepointsGeo(const parameters *p, 
				     int actP, int propP, int nP, 
				     const intparameters *ip, 
				     int actIP, int propIP, int nIP, 
				     likelihood *L) {
  int i;
  double logPar=0, logProp=0;
  int nData=getNdata();
  /* Geometric prior for NUMBER of changepoints: in contrast to GSL definition: 
     k=0, ..., p is SUCCESS probability */
  double nPriorPar=gsl_ran_geometric_pdf(1+actIP, 1-geoP);
  double nPriorProp=gsl_ran_geometric_pdf(1+propIP, 1-geoP);

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

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. The uniform order
   statistics prior is replaced by a negative binomial process
   prior (Fearnhead, 2006). */
double dPriorChangepointsNegativeBin(const parameters *p, 
				      int actP, int propP, int nP, 
				      const intparameters *ip, 
				      int actIP, int propIP, int nIP, 
				      likelihood *L) {
  double logPar=0, logProp=0;
  /* Prior for NUMBER of changepoints */
  double nPriorPar=gsl_ran_poisson_pdf(actIP, lambda);
  double nPriorProp=gsl_ran_poisson_pdf(propIP, lambda);

  /* Prior for number of changepoints */
  logPar+=nPriorPar;
  logProp+=nPriorProp;


  /* Prior for LOCATION. See page 205 in (Fearnhead,
     2006) */
  if((actIP >= 1) && (propIP >= 1)) {
    logProp+=logNegativeBinomialPriorDiff(ip, actIP, propIP, nIP);
  }

  setPrior(L, logPar, logProp);

  return logProp-logPar;
}

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. Uses a Poisson
   prior for the number of changepoint and a negative binomial prior
   for the changepoint locations with k=1 (geometric prior).*/
double dPriorChangepointsNegativeBinK1(const parameters *p, 
				       int actP, int propP, int nP, 
				       const intparameters *ip, 
				       int actIP, int propIP, int nIP, 
				       likelihood *L) {
  int i;
  double logPar=0, logProp=0;
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
    logPar+=(k1>k0+1)?(k1-k0-1)*log(1-nb_p):0;
  }

  for(i=0; i<=propIP; i++){
    int k0P=getCPprop(ip, i, propIP), k1P=getCPprop(ip, i+1, propIP);
    logProp+=(k1P>k0P+1)?(k1P-k0P-1)*log(1-nb_p):0;
  }

  /* 'Divide' by scaling factors */
  logPar+=actIP*log(nb_p);
  logProp+=propIP*log(nb_p);

  setPrior(L, logPar, logProp);

  return logProp-logPar;
}

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. Uses a geometric
   prior for the number of changepoint and a negative binomial prior
   for the changepoint locations with k=1 (geometric prior).*/
double dPriorChangepoints_nJGeoNegativeBinK1(const parameters *p, 
				       int actP, int propP, int nP, 
				       const intparameters *ip, 
				       int actIP, int propIP, int nIP, 
				       likelihood *L) {
  int i;
  double logPar=0, logProp=0;
  /* Prior for NUMBER of changepoints */
  /* Geometric prior for NUMBER of changepoints: in contrast to GSL definition: 
     k=0, ..., p is SUCCESS probability */
  double nPriorPar=gsl_ran_geometric_pdf(1+actIP, 1-geoP);
  double nPriorProp=gsl_ran_geometric_pdf(1+propIP, 1-geoP);

  /* Prior for number of changepoints */
  logPar+=nPriorPar;
  logProp+=nPriorProp;


  /* Prior for LOCATION. See bottom of page 207 in (Fearnhead,
     2006) */
  for(i=0; i<=actIP; i++){
    int k0=getCPpar(ip, i, actIP), k1=getCPpar(ip,i+1, actIP);
    logPar+=(k1>k0+1)?(k1-k0-1)*log(1-nb_p):0;
  }

  for(i=0; i<=propIP; i++){
    int k0P=getCPprop(ip, i, propIP), k1P=getCPprop(ip, i+1, propIP);
    logProp+=(k1P>k0P+1)?(k1P-k0P-1)*log(1-nb_p):0;
  }

  /* 'Divide' by scaling factors */
  logPar+=actIP*log(nb_p);
  logProp+=propIP*log(nb_p);

  setPrior(L, logPar, logProp);

  return logProp-logPar;
}

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. Uses a Poisson
   prior for the number of changepoint and a negative binomial prior
   for the changepoint locations with k=2.*/
double dPriorChangepointsNegativeBinK2(const parameters *p, 
				       int actP, int propP, int nP, 
				       const intparameters *ip, 
				       int actIP, int propIP, int nIP, 
				       likelihood *L) {
  int i;
  double logPar=0, logProp=0;
  /* Prior for NUMBER of changepoints */
  double nPriorPar=gsl_ran_poisson_pdf(actIP, lambda);
  double nPriorProp=gsl_ran_poisson_pdf(propIP, lambda);

  int k0, k1, k0P, k1P;

  /* Prior for number of changepoints */
  logPar+=nPriorPar;
  logProp+=nPriorProp;


  /* Prior for LOCATION. See bottom of page 207 in (Fearnhead,
     2006) */

  /* g0 */
  k0=getCPpar(ip, 1, actIP);
  logPar += log(nb_p) + (k0-1)*log(1-nb_p) + log(1+(k0-2)*nb_p*(1-nb_p));

  k0P=getCPprop(ip, 1, propIP);
  logProp += log(nb_p) + (k0P-1)*log(1-nb_p) + log(1+(k0P-2)*nb_p*(1-nb_p));
    
  /*g*/
  for(i=1; i<actIP; i++){
    k0=getCPpar(ip, i, actIP), k1=getCPpar(ip,i+1, actIP);
    logPar+=(k1>k0+2)?log(k1-k0-2) + 2*log(nb_p)  + (k1-k0-2)*log(1-nb_p):0;
  }

  for(i=1; i<propIP; i++){
    k0P=getCPprop(ip, i, propIP), k1P=getCPprop(ip, i+1, propIP);
    logProp+=(k1P>k0P+2)?log(k1P-k0P-2) + 2*log(nb_p)  + (k1P-k0P-2)*log(1-nb_p):0;
  }

  /*1-G*/
  k0=getCPpar(ip, actIP-1, actIP), k1=getCPpar(ip,actIP, actIP);
  logPar+=log(1-negativeBin_G(k1-k0));

  k0P=getCPprop(ip, propIP-1, propIP), k1P=getCPprop(ip,propIP, propIP);
  logProp+=log(1-negativeBin_G(k1P-k0P));

  setPrior(L, logPar, logProp);

  return logProp-logPar;
}

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. Uses a geometric
   prior for the number of changepoint and a negative binomial prior
   for the changepoint locations with k=2.*/
double dPriorChangepoints_nJGeoNegativeBinK2(const parameters *p, 
				       int actP, int propP, int nP, 
				       const intparameters *ip, 
				       int actIP, int propIP, int nIP, 
				       likelihood *L) {
  int i;
  double logPar=0, logProp=0;

  /* Prior for NUMBER of changepoints */
  /* Geometric prior for NUMBER of changepoints: in contrast to GSL definition: 
     k=0, ..., p is SUCCESS probability */
  double nPriorPar=gsl_ran_geometric_pdf(1+actIP, 1-geoP);
  double nPriorProp=gsl_ran_geometric_pdf(1+propIP, 1-geoP);

  int k0, k1, k0P, k1P;

  /* Prior for number of changepoints */
  logPar+=nPriorPar;
  logProp+=nPriorProp;

  /* Prior for LOCATION. See bottom of page 207 in (Fearnhead,
     2006) */

  /* g0 */
  k0=getCPpar(ip, 1, actIP);
  logPar += log(nb_p) + (k0-1)*log(1-nb_p) + log(1+(k0-2)*nb_p*(1-nb_p));

  k0P=getCPprop(ip, 1, propIP);
  logProp += log(nb_p) + (k0P-1)*log(1-nb_p) + log(1+(k0P-2)*nb_p*(1-nb_p));
    
  /*g*/
  for(i=1; i<actIP; i++){
    k0=getCPpar(ip, i, actIP), k1=getCPpar(ip,i+1, actIP);
    logPar+=(k1>k0+2)?log(k1-k0-2) + 2*log(nb_p)  + (k1-k0-2)*log(1-nb_p):0;
  }

  for(i=1; i<propIP; i++){
    k0P=getCPprop(ip, i, propIP), k1P=getCPprop(ip, i+1, propIP);
    logProp+=(k1P>k0P+2)?log(k1P-k0P-2) + 2*log(nb_p)  + (k1P-k0P-2)*log(1-nb_p):0;
  }

  /*1-G*/
  k0=getCPpar(ip, actIP-1, actIP), k1=getCPpar(ip,actIP, actIP);
  logPar+=log(1-negativeBin_G(k1-k0));

  k0P=getCPprop(ip, propIP-1, propIP), k1P=getCPprop(ip,propIP, propIP);
  logProp+=log(1-negativeBin_G(k1P-k0P));

  setPrior(L, logPar, logProp);

  return logProp-logPar;
}
