#ifndef CP_PRIOR
#define CP_PRIOR
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
#include "mp_parameter.h"
#include "likelihood.h"

/* Conjugate prior for open probabilities: Get hyperparameter alpha */
double getAlpha();

/* Conjugate prior for open probabilities: Get hyperparameter beta */
double getBeta();

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. */
double dPriorChangepoints(const parameters *p, 
			  int actP, int propP, int nP, 
			  const intparameters *ip, 
			  int actIP, int propIP, int nIP, 
			  likelihood *L);

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. */
double dPriorChangepointsGeo(const parameters *p, 
			     int actP, int propP, int nP, 
			     const intparameters *ip, 
			     int actIP, int propIP, int nIP, 
			     likelihood *L);

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
				     likelihood *L);

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. The uniform order
   statistics prior is replaced by a negative binomial process
   prior (Fearnhead, 2006). */
double dPriorChangepointsNegativeBinK1(const parameters *p, 
				      int actP, int propP, int nP, 
				      const intparameters *ip, 
				      int actIP, int propIP, int nIP, 
				     likelihood *L);

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. The uniform order
   statistics prior is replaced by a negative binomial process
   prior (Fearnhead, 2006). */
double dPriorChangepointsNegativeBinK1(const parameters *p, 
				      int actP, int propP, int nP, 
				      const intparameters *ip, 
				      int actIP, int propIP, int nIP, 
				     likelihood *L);

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. */
double dPriorChangepointsNegativeBinK2(const parameters *p, 
				       int actP, int propP, int nP, 
				       const intparameters *ip, 
				       int actIP, int propIP, int nIP, 
				       likelihood *L);

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. */
double dPriorChangepoints_nJGeoNegativeBinK1(const parameters *p, 
				       int actP, int propP, int nP, 
				       const intparameters *ip, 
				       int actIP, int propIP, int nIP, 
					     likelihood *L);

/* log-Prior for consecutive changepoints calculated for current
   sample and proposal. Priors for the NUMBER of changepoints and the
   changepoint LOCATIONS are added. The conjugate prior for the
   probabilities is accounted for in the likelihood. */
double dPriorChangepoints_nJGeoNegativeBinK2(const parameters *p, 
				       int actP, int propP, int nP, 
				       const intparameters *ip, 
				       int actIP, int propIP, int nIP, 
					     likelihood *L);
#endif
