/*
 * mp_proposal.c
 *
 * MCMC proposals
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
#include <string.h>
#include <math.h>

#ifndef HPC
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#else
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#endif

#include "mp_parameter.h"

#define EPS 1e-12
#define OUT stdout
#define ERR stderr

/*Returns next value for uniform random walk*/
double walk(const gsl_rng *r, double width) {
  return (gsl_rng_uniform(r)-0.5)*width;
}

/* Walk until positive */
double walkPositive(const gsl_rng *r, 
		    const double par, const double delta) {
  double w;
  do {
    w = walk(r, delta);
  } while(par + w <= 0);
  return w;
}

static void changePositive(const gsl_rng *r, parameters *p, 
			   int i, double width) {
  double par=getParameter(p, i);
  double w=walkPositive(r, par, width);
  setProposal(p, i, par+w);
}

/* Propose parameters */
void proposePositiveAB(const gsl_rng *r, parameters *p, int i0, int steps) {
  int i, n=i0+steps;

  for(i=i0; i<n; i++) {
    changePositive(r, p, i,  getWidth(p,i));
  }
}

/* check if proposal is in the right range */
int isProposalMinMax(const parameters *p, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    double prop = getProposal(p, i);
    if( (prop < getMin(p,i)) || (prop > getMax(p,i)) ) {
/*       fprintf(ERR, "negative proposal found:\n"); */
      return 0;
    }
  }
  return 1;
}


static void changeConstrained(const gsl_rng *r, parameters *p, 
			      int i) {
  double par=getParameter(p, i);
  double w;
  double min = p[i].min/* getMin(p,i) */, max = p[i].max/* getMax(p,i) */;
  do {
    double width = p[i].width/* getWidth(&p[i]); *//* getWidth(p,i) */;
    w = walk(r, width);
    /* fprintf(OUT, "changeConstrained(): i = %i, width = %f, par = %f, w = %f, min=%f, max=%f\n", */
    /* 	    i, width, par, w, min, max); */
    /* printFullParameters(OUT, &(p[i]), 1); */
  } while( (par + w < min) || (par + w > max) );
  setProposal(p, i, par+w);
  /* fprintf(OUT, "Proposal=%f\n",getProposal(p,i)); */
}

/* Propose parameters */
void proposeConstrainedAB(const gsl_rng *r, parameters *p, int i0, int steps) {
  int i;
  for(i=i0; i<i0+steps; i++) {
    changeConstrained(r, p, i);
  }
  /* fprintf(OUT, "proposeConstrainedAB():\n"); */
  /* printParameters(OUT, p, steps); */
  /* fprintf(OUT, "proposeConstrainedAB():\n"); */
  /* printProposal(OUT, p, steps); */
}

/* Propose parameters: Ratio of new and old is log(uNew/uOld)=u, where
   u is uniformly distributed over (-0.5, 0.5). */
void proposeExpAB(const gsl_rng *r, parameters *p, int i0, int steps) {
  int i;
  for(i=i0; i<i0+steps; i++) {
    double u=gsl_rng_uniform(r)-0.5;
    double oldP=getParameter(p,i);
    /* exp of uniform random number in interval ]-0.5,0.5 [*/
    /* fprintf(OUT, "proposeExpAB(): oldP = %f, newP=%f\n", oldP, oldP*exp(u)); */
    setProposal(p, i, oldP*exp(u));
    /* changeConstrained(r, p, i); */
  }
}



/* Does proposal match constraints ? */
int matchConstraintsK0(const parameters *p, int k0, int nPar) {
  int i;
  for(i=k0; i<k0+nPar; i++) {
    double prop = getProposal(p, i);
    double min = getMin(p, i);
    double max = getMax(p, i);
    if( (prop < min) || (prop > max)) {
      return 0;
    }
  }
  return 1;
}

/* Does proposal match constraints ? */
int matchConstraints(const parameters *p, int nPar) {
  return matchConstraintsK0(p, 0, nPar);
}

/* check if proposal is positive */
int isProposalPositive(const parameters *p, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    double prop = getProposal(p, i);
    if( prop < EPS) {
/*       fprintf(ERR, "negative proposal found:\n"); */
      return 0;
    }
  }
  return 1;
}

/*Returns next value for uniform random walk*/
static int intSymmWalk(const gsl_rng *r, int width) {
  return intUniform(r, -width, 2*width+1);
}

/* Walk until positive */
static int intSymmWalkPositive(const gsl_rng *r, 
				  const int par, const int delta) {
  int w;
  do {
    w = intSymmWalk(r, delta);
  } while(par + w <= 0);
  return w;
}

static void changeIntPositive(const gsl_rng *r, intparameters *p, 
			      int i, int width) {
  int par=getIntParameter(p, i);
  int w=intSymmWalkPositive(r, par, width);
  setIntProposal(p, i, par+w);
}

/* Propose parameters */
void proposeIntPositiveAB(const gsl_rng *r, intparameters *p, 
			 int i0, int steps) {
  int i, n=i0+steps;
  for(i=i0; i<n; i++) {
    changeIntPositive(r, p, i,  getIntWidth(p,i));
  }
}

/* Samples between the current values of the parameters and the maxima
   or the minima and the currrent values. */
void proposeIntIntervalAB(const gsl_rng *r, intparameters *p, 
			  int i0, int steps) {
  int i, n=i0+steps;
  for(i=i0; i<n; i++) {
    int par=getIntParameter(p, i);
    int newPar;
    int up=gsl_rng_uniform(r)<0.5;
    if(up) {
      int max=getIntMax(p, i);
      newPar=intUniform(r, par+1, max-(par+1)+1);
      setIntProposal(p, i, newPar);
      /* fprintf(OUT, ""); */
    }
    else {
      int min=getIntMin(p, i);
      newPar=intUniform(r, min, (par-min)+1);
      setIntProposal(p, i, newPar);
    }
  }
}

/* Samples between the current values of the parameters and the maxima
   or the minima and the currrent values. */
void proposeIntShiftInterval(const gsl_rng *r, intparameters *p, 
			     int i0, int nIPar) {
  /* Moves any parameter starting from i0... probably usually i0=0 */
  /* int par=getIntParameter(p, i0); */
  int newPar;
  int min=(i0>0)?getIntParameter(p,i0-1):getIntMin(p,i0);
  int max =(i0+1<nIPar)?getIntParameter(p,i0+1):getIntMax(p, i0);
  newPar=intUniform(r, min+1, max-(min+1)-1);
  /* fprintf(OUT, "proposeIntShiftInterval(): i0=%i, min=%i, max=%i, new=%i\n", */
  /* 	  i0, min+1, max-(min+1)-1, newPar); */
  setIntProposal(p, i0, newPar);
  /* int up=gsl_rng_uniform(r)<0.5; */
  /* if(up) { */
  /*   /\* Sample between current value of this and next parameter or */
  /*      between current value and maximal value *\/ */
  /*   int max=(i0+1<nIPar)?getIntParameter(p,i0+1):getIntMax(p, i0); */
  /*   /\* fprintf(OUT, "k%i = %i, k%i = %i, max=%i\n", i0, par,i0+1, getIntParameter(p,i0+1),max); *\/ */
  /*   if(par+1<max-(par+1)+1) newPar=intUniform(r, par+1, max-(par+1)+1);  */
  /*   else newPar=par; */
  /*   fprintf(OUT, "k%i: up: oldK = %i, newK = %i, upper = %i\n", i0, par, newPar, max); */
  /*   setIntProposal(p, i0,newPar);  */
  /* } else {  */
  /*   /\* Sample between preceding and current value or between minimum */
  /*      and current value. *\/ */
  /*   int min=(i0-1>=0)?getIntParameter(p,i0-1):getIntMin(p, i0);  */
  /*   if(min<(par-min)) newPar=intUniform(r, min, (par-min));  */
  /*   else newPar=par; */
  /*   fprintf(OUT, "k%i, down: oldK = %i, newK = %i,  lower = %i\n", i0, par, newPar, min); */
  /*   setIntProposal(p, i0, newPar);  */
  /* }  */
}


static void changeIntConstrained(const gsl_rng *r, intparameters *p, 
				 int i) {
  int par=getIntParameter(p, i), w;
  int width = p[i].width/* getIntWidth(p,i) */;
  int min=p[i].min/* getIntMin(p,i) */;
  int max=p[i].max/* getIntMax(p,i) */;

  do {
    w = intSymmWalk(r, width);
/*     w = intUpDown(r, width); */
    /* fprintf(OUT, "changeIntConstrained(): width = %i, par =%i, w = %i, min = %i, max = %i\n", */
    /* 	    width, par, w, min, max); */
  } while( (par + w < min) || (par + w > max) );
  setIntProposal(p, i, par+w);

/*   do { */
/*     sample = intUniform(r, par-width,2*width+1); */
/*   } while( (sample < min) || (sample > max) ); */
/*   setIntProposal(p, i, sample); */
}

/* Propose parameters */
void proposeIntConstrainedAB(const gsl_rng *r, intparameters *p, int i0, int steps) {
  int i;
  for(i=i0; i<i0+steps; i++) {
    changeIntConstrained(r, p, i);
  }
}


/* check if proposal is positive */
int isIntProposalPositive(const intparameters *p, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    int prop = getIntProposal(p, i);
    if( prop < EPS) {
/*       fprintf(ERR, "negative proposal found:\n"); */
      return 0;
    }
  }
  return 1;
}

/****************** The whole thing again... only simpler *************/
/* General form:
   propose(gsl_rng *r, const double *par, double *prop, likelihood *L, int nP);
 */

/* Generates proposals prop for all parameters par using a uniform
   random walk with width delta. Checks that proposals are within
   constraints min/max. */
void proposeAllUniformDouble(gsl_rng *r, 
			     const double *par, double *prop, 
			     double delta, 
			     int nP) {
  int i;

  for(i=0; i<nP; i++) {
    double w=walk(r, delta);
    prop[i]=par[i]+w;
  }
}

/* Generates proposals prop for all parameters par using a uniform
   random walk with width delta. Checks that proposals are within
   constraints min/max. */
void proposeOneUniformDouble(gsl_rng *r, 
			     const double *par, double *prop, 
			     int i0,
			     double delta, 
			     int nP) {
  double  w;

  if(i0>=nP) {
    fprintf(ERR, "proposeOneUniformDoubleMinMax(): index %i out of bounds (nP=%i\n", i0, nP);
    exit(1);
  }
    
  w=walk(r, delta);
  prop[i0]=par[i0]+w;
}


/* Generates proposals prop for all parameters par using a uniform
   random walk with width delta. Checks that proposals are within
   constraints min/max. */
void proposeAllUniformDoubleMinMax(gsl_rng *r, 
				   const double *par, double *prop, 
				   double delta, 
				   double min, double max, int nP) {
  int i;

  for(i=0; i<nP; i++) {
    double w=walk(r, delta);
    while( (par[i]+w <min) || (par[i]+w>max)) w=walk(r, delta);
    prop[i]=par[i]+w;
  }
}

/* Generates proposals prop[i0] for parameter par[i0] using a uniform
   random walk with width delta. Checks that proposals are within
   constraints min/max. */
void proposeOneUniformDoubleMinMax(gsl_rng *r, 
				   const double *par, double *prop, 
				   int i0,
				   double delta, 
				   double min, double max, int nP) {
  double  w;

  if(i0>=nP) {
    fprintf(ERR, "proposeOneUniformDoubleMinMax(): index %i out of bounds (nP=%i\n", i0, nP);
    exit(1);
  }
    
  w=walk(r, delta);
  while( (par[i0]+w <min) || (par[i0]+w>max)) w=walk(r, delta);
  prop[i0]=par[i0]+w;

}

/************************* t-Walk **************************************/
/* upper bound for moving parameters */
#define N1PHI 5

/* How many coordinates have been moved?*/
static int nTr=0;

/* Has one parameter been set to a value below zero ? */
static int zeroPar = 0;


/* probability for moving coordinate j*/
static double pMove=1;

/* for computing the probability that a coordinate is moved */
static double laphibound(double x) { return (-1.0/(x-1.0))*log((N1PHI-x)/(N1PHI-1.0)); }

/* Initialises the probability for moving a coordinate j*/
void initpMove(int nPar) {
  double nPhi= (N1PHI-(N1PHI-1)*exp(-laphibound(2)*(nPar-1)));
  pMove = nPhi/((double)nPar);
  fprintf(OUT, "pMove set to %f\n", pMove);
}

/************************* END of: initialisers ************************/

static double twalkScale = 1;

/* returns scaling factor for the acceptance ratio for last t-walk move */
double getTwalkScale() {
  return twalkScale;
}

/* returns scaling factor for the acceptance ratio for last t-walk move */
double getLogTwalkScale() {
  return log(twalkScale);
}


/* for walk move */
static double walkDist(const gsl_rng *r, double a) {
  double u = gsl_rng_uniform(r);
  return a/(a+1)*( -1 + 2*u + a*u*u);
}

static double hWalk(const gsl_rng *r, double x, double xP) {
  double out =  x + (x - xP)*walkDist(r, 0.5);
/*   fprintf(OUT, "hwalk(): x = %f, xP = %f, out = %f\n", x, xP, out); */
/*   exit(1); */
  return out;
}

static double hWalkPositive(const gsl_rng *r, double x, double xP) {
  double out = hWalk(r, x, xP);
  while ( out <= 0 ) {
    out = hWalk(r, x, xP);
  }
  return out;
}

static double hWalkConstrained(const gsl_rng *r, double x, double xP,
			       double xmin, double xmax) {
  double out = hWalk(r, x, xP);
  while ( (x <= xmin) || (x>=xmax) ) {
    out = hWalk(r, x, xP);
  }
  return out;
}

/* Applies move to a list of parameters */
static void moveDoubles(const gsl_rng *r, 
			double (*move)(const gsl_rng *, double, double),
			parameters *par, const parameters * parP, 
			int i0, int steps) {
  int i;

/*   printParam(OUT, par), printParam(OUT, parP); */
/*   fprintf(OUT, "par[%i] = %f\n", i0, getParameter(par, i0)); */

  for(i=i0; i<i0+steps; i++) {
    double x = getParameter(par, i),xP = getParameter(parP, i);
    double prop = move(r, x, xP);

    setProposal(par, i, prop);
  }
}

/* Applies move to a list of parameters */
static void moveDoublesPart(const gsl_rng *r, 
			    double (*move)(const gsl_rng *, double, double),
			    parameters *par, const parameters * parP, 
			    int i0, int steps) {
  int i;

  /*   printParam(OUT, par), printParam(OUT, parP); */
  /*   fprintf(OUT, "par[%i] = %f\n", i0, getParameter(par, i0)); */
  nTr=0;
  zeroPar=0;
  for(i=i0; i<i0+steps; i++) {
    double x = getParameter(par, i),xP = getParameter(parP, i);
    double prop = x;
    
    if(gsl_rng_uniform(r) < pMove) {
      /* fprintf(OUT, "moveDoublesPart(): Move parameter %i\n", i); */
      nTr++;
      prop=move(r, x, xP);
    }
    setProposal(par, i, prop);
  }
}

/* Applies move to a list of parameters */
static void moveDoublesConstrained(const gsl_rng *r, 
				   double (*move)(const gsl_rng *, double, double, double, double),
				   parameters *par, const parameters * parP, 
				   int i0, int steps) {
  int i;

  /*   printParam(OUT, par), printParam(OUT, parP); */
  /*   fprintf(OUT, "par[%i] = %f\n", i0, getParameter(par, i0)); */
  nTr=0;
  zeroPar=0;
  for(i=i0; i<i0+steps; i++) {
    double x = getParameter(par, i),xP = getParameter(parP, i);
    double prop = x;
    
    if(gsl_rng_uniform(r) < pMove) {
      /* fprintf(OUT, "moveDoublesPart(): Move parameter %i\n", i); */
      nTr++;
      prop=move(r, x, xP, par[i].min, par[i].max);
    }
    setProposal(par, i, prop);
  }
}

/* Implementation of the t-walk move 'walk' */
int sampleWalk(const gsl_rng *r, 
	       parameters *p1, int nP1,
	       parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleWalk()\n"), hello=0;
  moveDoubles(r,  hWalk,  p, prop, 
	      0, nP);
  twalkScale=1;
  return sampleFirst;
}


/* Implementation of the t-walk move 'walk' */
int sampleWalkPositive(const gsl_rng *r, 
		       parameters *p1, int nP1,
		       parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleWalk()\n"), hello=0;
  moveDoublesConstrained(r,  hWalkConstrained,  p, prop, 
			 0, nP);
  twalkScale=1;
  return sampleFirst;
}

/* Implementation of the t-walk move 'walk' */
int sampleWalkConstrained(const gsl_rng *r, 
			  parameters *p1, int nP1,
			  parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleWalkConstrained()\n"), hello=0;
  moveDoublesConstrained(r,  hWalkConstrained,  p, prop, 
			 0, nP);
  twalkScale=1;
  return sampleFirst;
}

/* Implementation of the t-walk move 'walk' - with constraints. */
int sampleWalkPositivePart(const gsl_rng *r, 
			   parameters *p1, int nP1,
			   parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleWalk()\n"), hello=0;
  moveDoublesPart(r,  hWalkPositive,  p, prop, 
		  0, nP);
  twalkScale=1;
  return sampleFirst;
}

/*************** traverse move ************************/

static double caseTraverse(double a) {
  return (a-1.0)/(2.0*a);
}

static double powxy(double x, double y) {
  if(x == 0) {
    if(y > 0) return 0;
    else {
      fprintf(ERR, "powxy(): Division by zero\n");
      exit(1);
    }
  }

  return exp(y*log(x));
}

/* for traverse move */
static double traverseDist(const gsl_rng *r, double a) {
  int case0 = ( gsl_rng_uniform(r) <= caseTraverse(a) ); 
  double u = gsl_rng_uniform(r);
  double expo, beta;

  if(case0) expo = 1.0/(a+1.0);
  else expo = 1.0/(1.0-a);

  beta = powxy(u, expo);
/*   fprintf(OUT, "traverseDist(): u = %f, expo = %f, u^expo = %f\n", */
/* 	  u, expo, beta); */
/*   exit(1); */

/*   twalkScale = pow(beta, nPar-2); */
/*   fprintf(OUT, "traverseDist(): nPar-2 = %i, scale = %f\n", */
/* 	  nPar-2,twalkScale); */


  return beta;
}

static double traverseBeta;

static double hTraverse(const gsl_rng *r, double x, double xP) {
  double out;

  out =  xP + (xP - x)*traverseBeta;

  if(out < 0) zeroPar = 1;

  return out;
}

/* changes parP according to walk move */
static void traverseMove(const gsl_rng *r, 
			 parameters *par, const parameters * parP, 
			 int i0, int steps) {
  traverseBeta = traverseDist(r, 4);
  moveDoublesPart(r, hTraverse, par, parP, i0, steps);
  twalkScale = pow(traverseBeta, nTr-2);
}

/* NEVER DO THAT! Traverse gets stuck when trying to force it! */
/* changes parP according to traverse move */
static void traverseMovePositivePart(const gsl_rng *r, 
				     parameters *par, const parameters * parP, 
				     int i0, int steps) {
  /* NEVER try to loop - traverse gets stuck! */
  do {
  /* zeroPar = 0; */
  /* nTr=0; */
    traverseBeta = traverseDist(r, 4);
    moveDoublesPart(r, hTraverse, par, parP, i0, steps);
    /* traverseMove(r, par, parP, i0, steps); */
  } while(zeroPar);
  /* if(zeroPar) twalkScale = 0, fprintf(ERR, "traverseMovePositive(): Negative sample...\n"); */
}

/* Implementation of the t-walk move 'traverse' */
int sampleTraversePositivePart(const gsl_rng *r, 
			   parameters *p1, int nP1,
			   parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleTraverse()\n"), hello=0;

  traverseMovePositivePart(r, p, prop, 0, nP);
  return sampleFirst;
}

/* Implementation of the t-walk move 'traverse' */
int sampleTraverse(const gsl_rng *r, 
			   parameters *p1, int nP1,
			   parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleTraverse()\n"), hello=0;

  traverseMove(r, p, prop, 0, nP);
  return sampleFirst;
}


/*************** END OF: traverse move ************************/

/* Implementation of the t-walk moves 'traverse' and 'walk' */
int sampleWalkTraversePositivePart(const gsl_rng *r, 
				   parameters *p1, int nP1,
				   parameters *p2, int nP2) {
  static int hello=1;

  if (hello) fprintf(ERR, "sampleWalkTraverse()\n"), hello=0;

  if(gsl_rng_uniform(r) < 0.5) return sampleWalkPositivePart(r, p1, nP1, p2, nP2);
  else return sampleTraversePositivePart(r, p1, nP1, p2, nP2);
}

/* Implementation of the t-walk moves 'traverse' and 'walk' */
int sampleWalkTraverse(const gsl_rng *r, 
		       parameters *p1, int nP1,
		       parameters *p2, int nP2) {
  static int hello=1;

  if (hello) fprintf(ERR, "sampleWalkTraverse()\n"), hello=0;

  if(gsl_rng_uniform(r) < 0.5) return sampleWalk(r, p1, nP1, p2, nP2);
  else return sampleTraverse(r, p1, nP1, p2, nP2);
}

/************* hop move *************************************/
/***************** helper routines for hop and blow moves ************/

/* changes parP according to walk move */
void hopMove(const gsl_rng *r, 
	     parameters *par,
	     parameters * const parP, int nP) {
  int i, nTr=0;
  double parMaxDist=0, propMaxDist=0;
  double sumPropSqrDeviates=0;

  /* Calculate max deviation between x and xP and initialise proposal
     with normal random variables or zeroes. */
  for (i=0; i<nP; i++) {
    /* Check if parameter shall move */
    if(gsl_rng_uniform(r) < pMove){
      double x=getParameter(par, i);
      double xP=getParameter(parP, i);
      double d=fabs(x-xP);
      if (d > parMaxDist) parMaxDist=d;
      /* Set par to N(0,1) */
      setProposal(par, i, gsl_ran_ugaussian (r));

      nTr++;
    }
    else setProposal(par, i, 0);
  }
  
  for(i=0; i< nP; i++) {
    double x=getParameter(par, i);
    double xOffs=getProposal(par, i)*parMaxDist/3.0;
    double xOffsSqr=xOffs*xOffs;
    double absXoffs=fabs(xOffs);
    double prop=x+xOffs;
    if(prop<0) zeroPar=1;
    if(absXoffs > propMaxDist) propMaxDist=absXoffs;

    sumPropSqrDeviates+=xOffsSqr;
    
    setProposal(par, i, prop);
  }

  /* Positive sign in exp due to dividing by g(y|x,x') */
  twalkScale = exp(9.0/(2.0*parMaxDist*parMaxDist)*sumPropSqrDeviates);
  twalkScale *= pow(parMaxDist/propMaxDist, nTr);
  /* fprintf(OUT, "hopMove(): twalkScale = %f\n", twalkScale); */

/*   fprintf(OUT, "maxParDist = %f, scalePar = %g\n", */
/* 	  maxParDist, scalePar); */
/*   fprintf(OUT, "maxPropDist = %f, propDist=%f, scaleProp = %g, twalkScale = %g, log(twalkScale)=%g\n", */
/* 	  maxPropDist, propDist, scaleProp, twalkScale, log(twalkScale)); */

/*   exit(1); */

}

/* Implementation of the t-walk move 'hop' */
int sampleHopMove(const gsl_rng *r, 
		  parameters *p1, int nP1,
		  parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleHop()\n"), hello=0;

  hopMove(r, p, prop, nP);
  return sampleFirst;
}

/* changes parP according to walk move */
void hopMovePositive(const gsl_rng *r, 
		     parameters *par, parameters * const parP, int nP) {
  do {
    zeroPar=0;
    hopMove(r, par, parP, nP);
  } while(zeroPar);
}

/* Implementation of the t-walk move 'traverse' */
int sampleHopMovePositive(const gsl_rng *r, 
			  parameters *p1, int nP1,
			  parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleHop()\n"), hello=0;

  hopMovePositive(r, p, prop, nP);
  return sampleFirst;
}

/* changes parP according to walk move */
void blowMove(const gsl_rng *r, 
	     parameters *par,
	     parameters * const parP, int nP) {
  int i, nTr=0;
  double parMaxDist=0, propMaxDist=0;
  double sumPropSqrDeviates=0;
  double sumParParSqrDeviates=0;

  /* Calculate max deviation between x and xP and initialise proposal
     with normal random variables or zeroes. */
  for (i=0; i<nP; i++) {
    /* Check if parameter shall move */
    if(gsl_rng_uniform(r) < pMove){
      double x=getParameter(par, i);
      double xP=getParameter(parP, i);
      double d=fabs(x-xP);
      if (d > parMaxDist) parMaxDist=d;
      sumParParSqrDeviates+=(x-xP)*(x-xP);
      /* Set par to N(0,1) */
      setProposal(par, i, gsl_ran_ugaussian (r));

      nTr++;
    }
    else setProposal(par, i, 0);
  }
  
  for(i=0; i< nP; i++) {
    double x=getParameter(par, i);
    double xP=getParameter(parP, i);
    double xOffs=getProposal(par, i)*parMaxDist;
    double xOffsSqr=xOffs*xOffs;
    double absXoffs=fabs(xOffs);
    /*Hoping that normal random variables are not too small */
    int movedI=fabs(getProposal(par,i))>=EPS;
    double prop=movedI?xP+xOffs:x;
    if(prop<0) zeroPar=1;
    if(absXoffs > propMaxDist) propMaxDist=absXoffs;

    sumPropSqrDeviates+=movedI?xOffsSqr:0;
    
    setProposal(par, i, prop);
  }

  /* Positive sign in exp due to dividing by g(y|x,x') */
  twalkScale = exp(1.0/(parMaxDist*parMaxDist)*sumPropSqrDeviates);
  twalkScale *= exp(-1.0/(propMaxDist*propMaxDist)*sumParParSqrDeviates);
  twalkScale *= pow(parMaxDist/propMaxDist, nTr);
  /* fprintf(OUT, "blowMove(): twalkScale = %f\n", twalkScale); */
}

/* Implementation of the t-walk move 'traverse' */
int sampleBlowMove(const gsl_rng *r, 
		   parameters *p1, int nP1,
		   parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleHop()\n"), hello=0;

  blowMove(r, p, prop, nP);
  return sampleFirst;
}

/* changes parP according to walk move */
void blowMovePositive(const gsl_rng *r, 
		      parameters *par, parameters * const parP, int nP) {
  do {
    zeroPar=0;
    blowMove(r, par, parP, nP);
  } while(zeroPar);
}

/* Implementation of the t-walk move 'traverse' */
int sampleBlowMovePositive(const gsl_rng *r, 
			   parameters *p1, int nP1,
			   parameters *p2, int nP2) {
  int sampleFirst=gsl_rng_uniform(r) < 0.5;
  parameters *p=sampleFirst?p1:p2;
  const parameters*prop=sampleFirst?p2:p1;
  int nP=sampleFirst?nP1:nP2;
  static int hello=1;

  if (hello) fprintf(ERR, "sampleHop()\n"), hello=0;

  blowMovePositive(r, p, prop, nP);
  return sampleFirst;
}

/* Implementation of the t-walk moves 'traverse' and 'walk' */
int sampleHopBlowPositive(const gsl_rng *r, 
			  parameters *p1, int nP1,
			  parameters *p2, int nP2) {
  static int hello=1;

  if (hello) fprintf(ERR, "sampleHopBlow()\n"), hello=0;

  if(gsl_rng_uniform(r) < 0.5) return sampleHopMovePositive(r, p1, nP1, p2, nP2);
  else return sampleBlowMovePositive(r, p1, nP1, p2, nP2);
}

/* Implementation of the t-walk moves 'traverse' and 'walk' */
int sampleHopBlow(const gsl_rng *r, 
		  parameters *p1, int nP1,
		  parameters *p2, int nP2) {
  static int hello=1;

  if (hello) fprintf(ERR, "sampleHopBlow()\n"), hello=0;

  if(gsl_rng_uniform(r) < 0.5) return sampleHopMove(r, p1, nP1, p2, nP2);
  else return sampleBlowMove(r, p1, nP1, p2, nP2);
}


/* Implementation of the t-walk move 'traverse' */
int sampleTwalkPositive(const gsl_rng *r, 
			parameters *p1, int nP1,
			parameters *p2, int nP2) {
  /* As in Christen and Fox (2010) */
  double pHopBlow=1.0/122.0;
  static int hello=0;
  if (hello) fprintf(ERR, "sampleTwalk()\n"), hello=0;
  
  if(gsl_rng_uniform(r) < pHopBlow)
    return sampleWalkTraversePositivePart(r,p1,nP1,p2,nP2);
  
  return sampleHopBlowPositive(r,p1,nP1,p2,nP2);
}

/* Implementation of the t-walk move 'traverse' */
int sampleTwalk(const gsl_rng *r, 
		parameters *p1, int nP1,
		parameters *p2, int nP2) {
  /* As in Christen and Fox (2010) */
  double pHopBlow=1.0/122.0;
  static int hello=0;
  if (hello) fprintf(ERR, "sampleTwalk()\n"), hello=0;
  
  if(gsl_rng_uniform(r) < pHopBlow)
    return sampleHopBlow(r,p1,nP1,p2,nP2);

  return sampleWalkTraverse(r,p1,nP1,p2,nP2);

}


/* /\* distance of proposal of par to parameter parP *\/ */
/* static double maxParParDistAB(parameters * const par, parameters * const parP, int i0, int steps) { */
/*   int i; */
/*   double max=0; */
/*   for(i=i0; i<i0+steps; i++) { */
/*     double prop=getProposal(par, i), xP = getParameter(parP, i); */
/*     double dist = fabs(prop-xP); */
/*     if(dist > max) max = dist; */
/*   } */
/*   return max; */
/* } */


/* /\* difference between parameter and proposal *\/ */
/* static double sqrPropParDistAB(parameters * const par, int i0, int steps){ */
/*   int i; */
/*   double sum=0; */

/*   for(i=i0; i<i0+steps; i++) { */
/*     double x = getProposal(par, i), xP = getParameter(par, i); */
/*     sum += (x-xP)*(x-xP); */
/*   } */

/*   return sum; */
/* } */

/* /\* difference between parameter and proposal *\/ */
/* static double sqrPropParDist(parameters * const par, int nP){ */
/*   return sqrPropParDistAB(par, 0, nP); */
/* } */

/* /\* difference between parameter and proposal *\/ */
/* static double sqrPropOtherParDistAB(parameters * const par, */
/* 				    parameters * const parP, */
/* 				    int i0, int steps){ */
/*   int i; */
/*   double sum=0; */

/*   for(i=i0; i<i0+steps; i++) { */
/*     double x = getProposal(par, i), xP = getParameter(parP, i); */
/*     sum += (x-xP)*(x-xP); */
/*   } */

/*   return sum; */
/* } */

/* /\* difference between parameter and proposal *\/ */
/* static double sqrPropOtherParDist(parameters * const par,  */
/* 				  parameters * const parP, int nP){ */
/*   return sqrPropOtherParDistAB(par, parP, 0, nP); */
/* } */

/* /\* difference between proposal and auxiliary parameter *\/ */
/* static double sqrParParDistAB(parameters * const par, parameters * const parP, */
/* 			      int i0, int steps){ */
/*   int i; */
/*   double sum=0; */

/*   for(i=i0; i<i0+steps; i++) { */
/*     double x = getParameter(par, i), xP = getParameter(parP, i); */
/*     sum += (x-xP)*(x-xP); */
/*   } */

/*   return sum; */
/* } */

/* /\* difference between proposal and auxiliary parameter *\/ */
/* static double sqrParParDist(parameters * const par, parameters * const parP, int nP){ */
/*   return sqrParParDistAB(par, parP, 0, nP); */
/* } */

/* /\* returns the normal deviates for a given standard deviation stdDev and */
/*    meanDist  */
/*    NO SCALING FACTORS (2 PI) INCLUDED !!! *\/ */
/* static double multiNormDist(double meanDist, double stdDev, int n) { */
/*   double stdDevN=pow(stdDev, n), var = stdDev*stdDev; */
/*   double meanOvar=meanDist/(2*var); */
/*   double normDist = exp(-meanOvar)/stdDevN; */
/* /\*   fprintf(OUT, "multiNormDist(): meanDist = %f, stdDevN=%g, var = %f, meanOvar = %f, exp(-meanOvar) = %g, normDist=%f\n", *\/ */
/* /\* 	  meanDist, stdDevN, var, meanOvar, exp(-meanOvar), normDist); *\/ */
/*   return normDist; */
/* } */

/* /\* for each parameter this vector contains *\/ */
/* static double *ZhopNblow; */

/* /\* allocate memory for hop and blow move *\/ */
/* static void inithopNblow(int nrates) { */
/*   ZhopNblow = malloc(nrates*sizeof(*ZhopNblow)); */
/* } */


/* /\* generate random numbers scaled by maximum distance which is */
/*    returned *\/ */
/* static double generateZAB(const gsl_rng *r,  */
/* 			  const parameters *par, const parameters *parP, */
/* 			  double *z, */
/* 			  int i0, int steps) { */
/*   int i; */
/*   double max=0; */

/*   /\* nTr=0; *\/ */
/*   /\* first run: find maximum and select z *\/ */
/*   for(i=i0; i<i0+steps; i++) { */
/*     double x = getParameter(par, i), xP=getParameter(parP, i); */
/*     double dist=fabs(x - xP); */
/* /\*     fprintf(OUT, "generateZ(): x = %f, xP = %f, dist = %f\n", *\/ */
/* /\* 	    x, xP, dist); *\/ */
/*     /\* decide if parameter is moved *\/ */
/*     if(gsl_rng_uniform(r) <= pMove) { */
/*       if(dist > max) max = dist; */

/*       ZhopNblow[i] = 1; */
/*       /\* nTr++; *\/ */
/*     } else ZhopNblow[i] = 0; */
/*   } */

/*   for(i=i0; i<i0+steps; i++) { */
/*     if(ZhopNblow[i]) ZhopNblow[i] = gsl_ran_gaussian (r, max); */
/* /\*     fprintf(OUT, "%f\t", ZhopNblow[i]); *\/ */
/*   } */
/* /\*   fprintf(OUT, "\nExiting...\n"); *\/ */
/* /\*   exit(1); *\/ */
/*   return max; */
/* } */

/* /\* generate random numbers scaled by maximum distance which is */
/*    returned *\/ */
/* static double generateZ(const gsl_rng *r,  */
/* 			const parameters *par, const parameters *parP, */
/* 			double *z, int nP) { */
/*   return generateZAB(r, par, parP,z, 0, nP); */
/* } */

/* /\* distance of proposal of par to parameter parP *\/ */
/* static double maxPropOtherParDistAB(parameters * const par, parameters * const parP, int i0, int steps) { */
/*   int i; */
/*   double max=0; */
/*   for(i=i0; i<i0+steps; i++) { */
/*     double prop=getProposal(par, i), xP = getParameter(parP, i); */
/*     double dist = fabs(prop-xP); */
/*     if(dist > max) max = dist; */
/*   } */
/*   return max; */
/* } */

/* /\* distance of proposal of par to parameter parP *\/ */
/* static double maxPropOtherParDist(parameters * const par, parameters * const parP, int nP) { */
/*   return maxPropOtherParDistAB(par, parP, 0, nP); */
/* } */

/* /\* changes parP according to walk move *\/ */
/* static void hopMove(const gsl_rng *r,  */
/* 		    parameters *par, parameters * const parP,  */
/* 		    int i0, int steps) { */
/*   int i; */
/*   double propDist; */
/*   double scaleProp, scalePar, maxParDist, maxPropDist; */

/*   maxParDist=generateZAB(r, par, parP, ZhopNblow, i0, steps); */

/*   for(i=i0; i<i0+steps; i++) { */
/*     double x=getParameter(par, i); */
/*     double prop=x + ZhopNblow[i]/3.0; */

/*     setProposal(par, i, prop); */
/*   } */
/*   maxPropDist=maxPropOtherParDistAB(par, parP, i0, steps); */
/*   propDist=sqrPropParDistAB(par, i0, steps); */
/* /\*   parDist=sqrParParDist(par, parP); *\/ */

/*   /\* scalePar=multiNormDist(propDist, maxParDist/3.0, nToMove); *\/ */
/*   /\* scaleProp=multiNormDist(propDist, maxPropDist/3.0, nToMove); *\/ */


/*   twalkScale = scaleProp/scalePar; */

/* /\*   fprintf(OUT, "maxParDist = %f, scalePar = %g\n", *\/ */
/* /\* 	  maxParDist, scalePar); *\/ */
/* /\*   fprintf(OUT, "maxPropDist = %f, propDist=%f, scaleProp = %g, twalkScale = %g, log(twalkScale)=%g\n", *\/ */
/* /\* 	  maxPropDist, propDist, scaleProp, twalkScale, log(twalkScale)); *\/ */

/* /\*   exit(1); *\/ */

/* } */

/* /\* changes parP according to walk move *\/ */
/* static void blowMove(const gsl_rng *r,  */
/* 		    parameters *par, parameters * const parP,  */
/* 		    int i0, int steps) { */
/*   int i; */
/*   double propDist, parDist; */
/*   double scaleProp, scalePar, maxParDist, maxPropDist; */

/*   maxParDist=generateZAB(r, par, parP, ZhopNblow, i0, steps); */

/*   for(i=i0; i<i0+steps; i++) { */
/*     double xP=getParameter(parP, i); */
/*     double prop=xP + ZhopNblow[i]; */

/*     setProposal(par, i, prop); */
/*   } */
/*   maxPropDist=maxPropOtherParDistAB(par, parP, i0, steps); */
/*   propDist=sqrPropOtherParDistAB(par, parP, i0, steps); */
/*   parDist=sqrParParDistAB(par, parP, i0, steps); */

/*   /\* scalePar=multiNormDist(propDist, maxParDist, nToMove); *\/ */
/*   /\* scaleProp=multiNormDist(parDist, maxPropDist, nToMove); *\/ */


/*   twalkScale = scaleProp/scalePar; */

/* /\*   fprintf(OUT, "maxParDist = %f, scalePar = %g\n", *\/ */
/* /\* 	  maxParDist, scalePar); *\/ */
/* /\*   fprintf(OUT, "maxPropDist = %f, propDist=%f, scaleProp = %g, twalkScale = %g, log(twalkScale)=%g\n", *\/ */
/* /\* 	  maxPropDist, propDist, scaleProp, twalkScale, log(twalkScale)); *\/ */

/* /\*   exit(1); *\/ */

/* } */
/************************* END: hop/blow move **************************/

/************************* int Parameter *******************************/
/* Generates proposals prop for all parameters par using a uniform
   random walk with width delta. Checks that proposals are within
   constraints min/max. */
void proposeAllUniformInt(gsl_rng *r, 
			  const int *par, int *prop, 
			  int delta, 
			  int nP) {
  int i;

  for(i=0; i<nP; i++) {
    double w=walk(r, delta);
    prop[i]=par[i]+w;
  }
}

/* Generates proposals prop for all parameters par using a uniform
   random walk with width delta. Checks that proposals are within
   constraints min/max. */
/* void proposeOneUniformDouble(gsl_rng *r,  */
/* 			     const double *par, double *prop,  */
/* 			     int i0, */
/* 			     int delta,  */
/* 			     int nP) { */
/*   double  w; */

/*   if(i0>=nP) { */
/*     fprintf(ERR, "proposeOneUniformDoubleMinMax(): index %i out of bounds (nP=%i\n", i0, nP); */
/*     exit(1); */
/*   } */
    
/*   w=walk(r, delta); */
/*   prop[i0]=par[i0]+w; */
/* } */


/* Generates proposals prop for all parameters par using a uniform
   random walk with width delta. Checks that proposals are within
   constraints min/max. */
/* void proposeAllUniformDoubleMinMax(gsl_rng *r,  */
/* 				   const double *par, double *prop,  */
/* 				   double delta,  */
/* 				   double min, double max, int nP) { */
/*   int i; */

/*   for(i=0; i<nP; i++) { */
/*     double w=walk(r, delta); */
/*     while( (par[i]+w <min) || (par[i]+w>max)) w=walk(r, delta); */
/*     prop[i]=par[i]+w; */
/*   } */
/* } */

/* Generates proposals prop[i0] for parameter par[i0] using a uniform
   random walk with width delta. Checks that proposals are within
   constraints min/max. */
/* void proposeOneUniformDoubleMinMax(gsl_rng *r,  */
/* 				   const double *par, double *prop,  */
/* 				   int i0, */
/* 				   double delta,  */
/* 				   double min, double max, int nP) { */
/*   double  w; */

/*   if(i0>=nP) { */
/*     fprintf(ERR, "proposeOneUniformDoubleMinMax(): index %i out of bounds (nP=%i\n", i0, nP); */
/*     exit(1); */
/*   } */
    
/*   w=walk(r, delta); */
/*   while( (par[i0]+w <min) || (par[i0]+w>max)) w=walk(r, delta); */
/*   prop[i0]=par[i0]+w; */

/* } */

