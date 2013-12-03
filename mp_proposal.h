#ifndef PROPOSAL
#define PROPOSAL
/*
 * mp_parameter.c
 *
 * MCMC parameter data type
 * 
 * 
 *
 * Ivo Siekmann, 10/01/2012
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */
/*Returns next value for uniform random walk*/
double walk(const gsl_rng *r, double width);

/* Walk until positive */
double walkPositive(const gsl_rng *r, 
		    const double par, const double delta);

/* Propose parameters */
void proposePositiveAB(const gsl_rng *r, parameters *p, int i0, int steps);

/* Propose parameters */
void proposeIntPositiveAB(const gsl_rng *r, intparameters *p, 
			  int i0, int steps);

/* Propose parameters */
void proposeConstrainedAB(const gsl_rng *r, parameters *p, int i0, int steps);

/* Propose parameters: Ratio of new and old is log(uNew/uOld)=u, where
   u is uniformly distributed over (-0.5, 0.5). */
void proposeExpAB(const gsl_rng *r, parameters *p, int i0, int steps);

/* Propose parameters */
void proposeIntConstrainedAB(const gsl_rng *r, intparameters *p, int i0, int steps);

/* Samples between the current values of the parameters and the maxima
   or the minima and the currrent values. */
void proposeIntIntervalAB(const gsl_rng *r, intparameters *p, 
			  int i0, int steps);

/* Samples between the current values of the parameters and the maxima
   or the minima and the currrent values. */
void proposeIntShiftInterval(const gsl_rng *r, intparameters *p, 
			     int i0, int steps);

/* returns scaling factor for the acceptance ratio for last t-walk move */
double getLogTwalkScale();

/* Implementation of the t-walk move 'walk' */
int sampleWalkPositive(const gsl_rng *r, 
		       parameters *p1, int nP1,
		       parameters *p2, int nP2);

/* Implementation of the t-walk move 'walk' */
int sampleWalk(const gsl_rng *r, 
	       parameters *p1, int nP1,
	       parameters *p2, int nP2);

int sampleWalkPositivePart(const gsl_rng *r, 
			   parameters *p1, int nP1,
			   parameters *p2, int nP2);

/* Implementation of the t-walk move 'traverse' */
int sampleTraversePositive(const gsl_rng *r, 
			   parameters *p1, int nP1,
			   parameters *p2, int nP2);

/* Implementation of the t-walk move 'traverse' */
int sampleWalkTraversePositive(const gsl_rng *r, 
			   parameters *p1, int nP1,
			       parameters *p2, int nP2);

/* Implementation of the t-walk move 'traverse' */
int sampleTraversePositivePart(const gsl_rng *r, 
			   parameters *p1, int nP1,
			       parameters *p2, int nP2);

/* Implementation of the t-walk move 'traverse' */
int sampleTraverse(const gsl_rng *r, 
		   parameters *p1, int nP1,
		   parameters *p2, int nP2);

/* Implementation of the t-walk moves 'traverse' and 'walk' */
int sampleWalkTraversePositivePart(const gsl_rng *r, 
			   parameters *p1, int nP1,
				   parameters *p2, int nP2);

/* Implementation of the t-walk moves 'traverse' and 'walk' */
int sampleWalkTraverse(const gsl_rng *r, 
		       parameters *p1, int nP1,
		       parameters *p2, int nP2);

/* Initialises the probability for moving a coordinate j*/
void initpMove(int nPar);

/* changes parP according to walk move */
void hopMove(const gsl_rng *r, 
	     parameters *par, int nP1,
	     parameters * const parP, int nP2);

/* Implementation of the t-walk move 'traverse' */
int sampleHopMove(const gsl_rng *r, 
		  parameters *p1, int nP1,
		  parameters *p2, int nP2);

/* Implementation of the t-walk move 'traverse' */
int sampleHopMovePositive(const gsl_rng *r, 
			  parameters *p1, int nP1,
			  parameters *p2, int nP2);

/* Implementation of the t-walk move 'blow' */
int sampleBlowMove(const gsl_rng *r, 
		   parameters *p1, int nP1,
		   parameters *p2, int nP2);

/* Implementation of the t-walk moves 'traverse' and 'walk' */
int sampleHopBlow(const gsl_rng *r, 
		  parameters *p1, int nP1,
		  parameters *p2, int nP2);

/* Implementation of the t-walk move 'blow' */
int sampleBlowMovePositive(const gsl_rng *r, 
			   parameters *p1, int nP1,
			   parameters *p2, int nP2);

/* Implementation of the t-walk move 'traverse' */
int sampleTwalkPositive(const gsl_rng *r, 
			parameters *p1, int nP1,
			parameters *p2, int nP2);

/* Implementation of the t-walk move 'traverse' */
int sampleTwalk(const gsl_rng *r, 
		parameters *p1, int nP1,
		parameters *p2, int nP2);
#endif
