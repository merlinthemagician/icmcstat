#ifndef LIKELIHOOD
#define LIKELIHOOD
/*
 *  likelihood.h
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

#include <stdio.h>
#include <stdlib.h>

/* Represents total likelihood as well as likelihoods of prior and
   posterior */
typedef struct likelihood {
  double *prior;
  int nPrior; /* 1 or 3*/
  double *posterior;
  int nPosterior; /* 1 or 3*/
  double *total;
  int nTotal; /* 1 or 3*/
} likelihood;


/* allocates likelihood */
likelihood* allocLikelihood();

/* allocates diff likelihood */
likelihood* allocDiffLikelihood();


/* frees memory allocated for L */
void freeLikelihood(likelihood *L);

/* prints likelihood */
void printLikelihood(FILE *fp, const likelihood *L);

/*prints total likelihood only*/
void printTotalL(FILE *fp, const likelihood* L);

/*prints prior likelihood only*/
void printPriorL(FILE *fp, const likelihood* L);

/*prints posterior likelihood only*/
void printPosteriorL(FILE *fp, const likelihood* L);

/* sets likelihoods for current state and proposal, computes
   differences and total */
void setLikelihood(likelihood *L, double prior, double priorP,
		   double post, double postP);

/* sets likelihoods for current state and proposal, computes
   differences and total */
void setDiffLikelihood(likelihood *L, double priorDiff, double postDiff);

/*sets posterior likelihood diff */
void setPostDiff(likelihood *L, double postDiff);

/*sets pior likelihood */
void setPrior(likelihood *L, double priorSample, double priorProp);

/*sets posterior likelihood */
double setPosterior(likelihood *L, double postSample, double postProp);

/* returns posterior likelihood diff */
double getPostDiff(const likelihood * L);

/* Adds offsets dLold and dLprop to the likelihood */
double addToLikelihoodL1L2(likelihood *L, double dLold, double dLprop);

/* Adds deltaL to the likelihood */
void addToLikelihood(likelihood *L, double deltaL);
#endif

