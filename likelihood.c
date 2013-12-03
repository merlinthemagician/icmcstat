/*
 *  likelihood.c
 *
 *
 * Data structure for likelihood representation
 * 
 * 
 *
 * Ivo Siekmann, 15/10/2009
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include "likelihood.h"

#define ERR stderr

/* allocates diff likelihood */
static likelihood* generalNalloc(int n) {
  likelihood* out=malloc(sizeof(likelihood));
  out->prior=malloc(n*sizeof(double)), out->nPrior=n;
  out->posterior=malloc(n*sizeof(double)), out->nPosterior=n;
  out->total=malloc(n*sizeof(double)), out->nTotal=n;

  return out;
}

/* frees memory allocated for L */
void freeLikelihood(likelihood *L) {
  free(L->prior);
  free(L->posterior);
  free(L->total);
  free(L);
}

/* allocates likelihood */
likelihood* allocLikelihood() {
  return generalNalloc(3);
}

/* allocates diff likelihood */
likelihood* allocDiffLikelihood() {
  return generalNalloc(1);
}

static void printComponent(FILE *fp, const double *v, int n) {
  if(n==1) {
    fprintf(fp, "%f", v[0]);
  } else if(n==3) {
    fprintf(fp, "%f\t%f\t%f", v[0], v[1], v[2]);
  }
  else {
    fprintf(ERR, "Warning: Invalid likelihood component\n");
  }
}

/*prints prior likelihood only*/
void printPriorL(FILE *fp, const likelihood* L) {
  printComponent(fp, L->prior, L->nPrior);
}

/*prints posterior likelihood only*/
void printPosteriorL(FILE *fp, const likelihood* L) {
  printComponent(fp, L->posterior, L->nPosterior);
}

/*prints total likelihood only*/
void printTotalL(FILE *fp, const likelihood* L) {
  printComponent(fp, L->total, L->nTotal);
}


/* prints likelihood */
void printLikelihood(FILE *fp, const likelihood *L) {
/*   fprintf(fp, "Prior:\t"); */
  printPriorL(fp, L);
  fprintf(fp, "\t");

/*   fprintf(fp, "Posterior:\t"); */
  printPosteriorL(fp, L);
  fprintf(fp, "\t");

/*   fprintf(fp, "total:\t"); */
  printTotalL(fp, L);
  fprintf(fp, "\t");

  fprintf(fp, "\n");
}


/* Sets component v to values in newV */
static void setComponent(double *v, double d, double dP) {
  /* fprintf(stdout, "setComponent(): Setting component:\n"); */
  /* fprintf(stdout, "setComponent(): Pointer %p:\n", v); */
  v[0] = d, v[1]=dP;
  v[2]=v[1]-v[0];
}

static void printLengthWarning(int length) {
  fprintf(ERR, "Length above %i expected\n", length);
  exit(1);
}

static double totalLikelihood(likelihood *L) {
  int i;
/*   fprintf(stdout, "nPrior=%i, nPost=%i, nTotal = %i\n",  */
/* 	  L->nPrior, L->nPosterior, L->nTotal); */
  for(i=0; i< (L->nTotal); i++) {
    if( (i >= L->nPrior) || (i >= L->nPosterior) ) printLengthWarning(i);
    L->total[i]=(L->prior[i])+(L->posterior[i]);
  }
  return L->total[L->nTotal-1];
}

/* sets likelihoods for current state and proposal, computes
   differences and total */
void setLikelihood(likelihood *L, double prior, double priorP,
		   double post, double postP) {
  if(L->nPrior != 3) printLengthWarning(3);
  setComponent(L->prior, prior, priorP);

  if(L->nPosterior != 3) printLengthWarning(3);
  setComponent(L->posterior, post, postP);

  if(L->nTotal != 3) printLengthWarning(3);
  totalLikelihood(L);
}

/* returns posterior likelihood diff */
double getPostDiff(const likelihood * L) {
  return L->posterior[0];
}

/*sets posterior likelihood */
void setPrior(likelihood *L, double priorSample, double priorProp) {
/*   L->posterior[0] = postSample; */
/*   L->posterior[1] = postProp; */
/*   L->posterior[2] = L->posterior[1]-L->posterior[0]; */
  setComponent(L->prior, priorSample, priorProp);
  totalLikelihood(L);
}

/*sets posterior likelihood */
double setPosterior(likelihood *L, double postSample, double postProp) {
/*   L->posterior[0] = postSample; */
/*   L->posterior[1] = postProp; */
/*   L->posterior[2] = L->posterior[1]-L->posterior[0]; */
  /* fprintf(stdout, "setPosterioir(): Setting likelihood\n"); */
  setComponent(L->posterior, postSample, postProp);
  /* fprintf(stdout, "setPosterior(): done\n"); */
  return totalLikelihood(L);
  /* fprintf(stdout, "setPosterior(): Setting likelihood\n"); */
}

/* returns likelihood for proposal */
double getLikelihoodProposal(const likelihood *L) {
  if(L->nPosterior!=3)  printLengthWarning(3), exit(1);

  return L->posterior[1];
}

/* returns likelihood for proposal */
double getLikelihoodParameter(const likelihood *L) {
  if(L->nPosterior!=3)  printLengthWarning(3), exit(1);

  return L->posterior[0];
}

/* returns likelihood for proposal */
void setLikelihoodProposal(const likelihood *L, double d) {
  if(L->nPosterior!=3)  printLengthWarning(3), exit(1);

  L->posterior[1]=d;
}

/* returns likelihood for proposal */
void setLikelihoodParameter(const likelihood *L, double d) {
  if(L->nPosterior!=3)  printLengthWarning(3), exit(1);

  L->posterior[0]=d;
}

/* Adds offsets dLold and dLprop to the likelihood */
double addToLikelihoodL1L2(likelihood *L, double dLold, double dLprop) {
  double Lold, Lprop, Ltot=0;
  if(L->nPosterior==3) {
    Lold=getLikelihoodParameter(L);
    Lprop=getLikelihoodProposal(L);
    Ltot = setPosterior(L,Lold+dLold, Lprop+dLprop);
  }
  return Ltot;
}


/* Adds deltaL to the likelihood */
void addToLikelihood(likelihood *L, double deltaL) {
  double Lold;
  if(L->nPosterior==3) {
    Lold=getLikelihoodProposal(L);
    setLikelihoodProposal(L,Lold+deltaL);
  }
  totalLikelihood(L);  
}

/*sets posterior likelihood diff */
void setPostDiff(likelihood *L, double postDiff) {
  L->posterior[0] = postDiff;
  totalLikelihood(L);
}

/* sets likelihoods for current state and proposal, computes
   differences and total */
void setDiffLikelihood(likelihood *L, double priorDiff, double postDiff) {
  L->prior[0]=priorDiff;
  L->posterior[0]=postDiff;
  totalLikelihood(L);
}
