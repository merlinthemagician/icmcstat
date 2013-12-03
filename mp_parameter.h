#ifndef PARAMETER
#define PARAMETER
/*
 * mp_parameter.h
 *
 * MCMC parameter data type
 * 
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

#include <stdio.h>
#include <stdlib.h> 

#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

/* #define NRATES 20 */

/* Data structure for double parameters: original parameters and
   proposal are both saved as an array of double pointers.  Each
   parameter can have a name (string). */
typedef struct parameters{
  char *name;
  double  *org;
  double *prop;
  double width;
  double min;
  double max;
} parameters;

/* Data structure for int parameters: original parameters and
   proposal are both saved as an array of double pointers.  Each
   parameter can have a name (string). */
typedef struct intparameters{
  char *name;
  int  *org;
  int *prop;
  int width;
  double min;
  double max;
} intparameters;

/* Returns uniformly distributed integers from {i0, ..., i0+n-1} */
int intUniform(const gsl_rng *r, int i0, int n);

/*initialises parameter with default values for width, min and max */
void initialiseParameters(parameters *p, double width, double min, double max, int nPar);

/* gives i-th sampling width of the parameters */
double getWidth(const parameters *p, int i);

/* set max of i-th parameter to d */
void setMax(parameters *p, int i, double d);

/* gives minimum constraint of parameter i */
double getMin(const parameters *p, int i);

/* gives maximum constraint of parameter i */
double getMax(const parameters *p, int i);

/* Copies paramerter new to dest */
void parameterMemcpy(parameters *dest, const parameters *new, int nPar);

/* allocates memory for n parameters */
intparameters *allocIntParameters(int n);

/*initialises parameter with default values for width, min and max */
void initialiseIntParameters(intparameters *p, int width, int min, int max, int nPar);

/* Copies parameter new to dest */
void intParameterMemcpy(intparameters *dest, const intparameters *new, int nPar);

/* initialise parameters of a continuous-time Markov model from the
   matrix Q */
/* parameters *initMarkovParameters(gsl_matrix *Q); */

/* saves parameters in file fp, starting from index k0 */
void printParametersK0(FILE *fp, const parameters *p, int k0, int nPar);

/* saves parameters in file fp */
void printParameters(FILE *fp, const parameters *p, int nPar);

/* saves parameters in file fp with all additional information */
void printFullParameters(FILE *fp, const parameters *p, int nPar);

/* saves parameters in file fp */
void printIntParameters(FILE *fp, const intparameters *p, int nPar);

/* copy proposal after accepting */
void copyAcceptedProposal(parameters *p);

/* set parameters to components of d */
void setParameters(parameters *p, int pn0, 
		   double* const d, int dn0, int nD);

/* set proposals to components of d */
void setProposals(parameters *p, int pn0,
		  double* const d, int dn0, int nD);

/* initialise parameters with random numbers, uniformly distributed
   over [0,1] */
void initParUniform(parameters *p);

/* check if proposal is positive */
int isProposalPositive();

/*prints parameter names to file fp */
void printParameterNames(FILE *fp, const parameters *p, int nPar);

/*prints parameter names to file fp */
void printIntParameterNames(FILE *fp, const intparameters *p, int nPar);

/* copies values from proposal to original parameters */
void copyProp2Org(parameters *p, int nPar);

/* copies values from proposal to original parameters */
void copyProp2OrgK0(parameters *p, int k0, int nPar);

/* copies values from proposal to original parameters */
void copyIntProp2Org(intparameters *ip, int nPar);

/* gives i-th of the proposals */
double getParameter(const parameters *p, int i);

/* gives i-th of the proposal */
double getProposal(const parameters *p, int i);

/* set i-th sampling width of the parameters */
void setWidth(parameters *p, int i, double d);

/* set i-th sampling width of the parameters */
void setIntWidth(intparameters *p, int i, int d);

/* set i-th parameter to d */
void setParameter(parameters *p, int i, double d);

/* set i-th parameter to d */
void setParameterPtr(parameters *p, int i, double *d);

/* Swap parameters i and j */
void swapParameters(parameters *p, int i, int j);

/* set i-th proposal to d */
void setProposal(parameters *p, int i, double d);

/* set i-th parameter to d */
void setProposalPtr(parameters *p, int i, double *d);

/* sets name of i-th parameter */
void setParameterName(parameters *p, int i, char * const name);

/* sets names of parameters, starting index at k0 */
void setParameterNamesK0(parameters *p, char * const prefix, int k0, int n);

/* sets name of i-th parameter */
void setParameterNames(parameters *p, char * const prefix, int n);

/* Initialises Parameters with pointers to elements in double array */
void initialisePorg(parameters *p, double *d, int n);

/* Initialises Parameter proposals with pointers to elements in double array */
void initialisePprop(parameters *p, double *d, int n);

/* Initialises Parameters with pointers to elements in double array */
void initialiseIntPorg(intparameters *ip,  int *k, int n);

/* Initialises parameter proposals with pointers to elements in double array */
void initialiseIntPprop(intparameters *ip, int *k, int n);

/* sets name of i-th parameter */
void setIntParameterName(intparameters *p, int i, char *name);

/* sets names of parameters */
void setIntParameterNames(intparameters *ip, char * const prefix, int n);

/* set i-th parameter to d */
void setIntParameter(intparameters *p, int i, int d);

/* set i-th parameter to d */
void setIntParameterPtr(intparameters *p, int i, int *d);

/* set i-th parameter to d */
void setIntProposalPtr(intparameters *p, int i, int *d);

/* gives i-th of the proposals */
int getIntParameter(const intparameters *p, int i);

/* gives i-th of the proposal */
int getIntProposal(const intparameters *p, int i);

/* gives i-th sampling width of the parameters */
int getIntMin(const intparameters *p, int i);

/* gives i-th sampling width of the parameters */
int getIntMax(const intparameters *p, int i);

/* Prints constraints for each int parameter: width of random walk,
   min, max */
void printIntConstraints(FILE *fp, const intparameters *p, int nPar);

/* Does proposal match constraints ? */
int matchConstraintsK0(const parameters *p, int k0, int nPar);

/* Does proposal match constraints ? */
int matchConstraints(const parameters *p, int nPar);

/* resets all proposals to original parameter values */
void resetProposals(parameters *p, int nPar,
		    intparameters *ip, int nIPar);


/* initialise parameters with random numbers, uniformly distributed
   over [0,1] */
void initUniform(const gsl_rng *r, parameters *par, int nPar);

/* initialise parameters with random numbers, uniformly distributed
   over [0,scale] */
void initUniformScale(const gsl_rng *r, parameters *par, 
		      double scale, int nPar);

/* initialise parameters with random numbers, uniformly distributed
   over {i0,..., n-1}, sort in ascending order */
void initUniformAscending(const gsl_rng *r, parameters *par, int nPar);

/* gives width of parameter i */
int getIntWidth(const intparameters *p, int i);

/* initialise parameters with random numbers, uniformly distributed
   over {i0,..., n-1} */
void initIntUniform(const gsl_rng *r, intparameters *par, int i0, int n, int nPar);

/* initialise parameters with random numbers, uniformly distributed
   over {i0,..., n-1}, sort in ascending order */
void initIntUniformAscending(const gsl_rng *r, intparameters *par, int i0, int n, int nPar);

/* copies values from original parameters to proposal */
void copyOrg2Prop(parameters *p, int nPar);

/* copies values from original parameters to proposal */
void copyIntOrg2Prop(intparameters *p, int nIPar);

/* set i-th proposal to d */
void setIntProposal(intparameters *p, int i, int d);

/* proposal, uniformly distributed over {i0,..., n-1}, sort in
   ascending order */
void proposeIntUniformAscending(const gsl_rng *r, intparameters *par, int i0, int n, int nPar);

/* gives i-th sampling width of the parameters */
int getIntMin(const intparameters *p, int i);

/* set i-th sampling width of the parameters */
void setIntMin(intparameters *p, int i, int d);

/* gives i-th sampling width of the parameters */
int getIntMax(const intparameters *p, int i);

/* set i-th sampling width of the parameters */
void setIntMax(intparameters *p, int i, int d);

#endif
