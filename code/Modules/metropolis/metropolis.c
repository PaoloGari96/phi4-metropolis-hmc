
/*
 *  File metropolic.c
 *  Metropolis algorithm for the phi4 theory
 *
 *  double action(act_params_t *apars)
 *      This routine computes the action S[phi] for the global field phi in
 *      lattice.h and the parameters kappa and lambda from apars.
 *      S = Sum_x [ -2*kappa*sum_mu phi_x phi_{x+mu}+phi_x^2+lambda(phi_x^2-1)^2 ]
 */

#include "metropolis.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


double action(act_params_t *apars)
{
   int i,mu;
   double phin,S,phi2;
   double kappa =apars->kappa;
   double lambda=apars->lambda;

   S=0;

   /* loop over all sites */
   for (i=0;i<V;i++)
   {
      /*sum over neighbors in positive direction*/
      phin=0;
      for (mu=0;mu<D;mu++) phin+=phi[hop[i][mu]];

      phi2=phi[i]*phi[i];
      S+=-2*kappa*phin*phi[i]+phi2+lambda*(phi2-1.0)*(phi2-1.0);
   }

   return S;
}

/********************************************************************** 
 * metropolis()
 * Does ntherm+nsweep sweeps over the whole lattice of the Metropolis
 * algorithm. Measurement after each sweep>=ntherm, the averaged measured 
 * values are printed out in fixed intervals, controlled by naccu.
 **********************************************************************/

double metropolis(act_params_t *apars, metro_params_t *mpars)
{
   int i,acc;
   static double rand[2];
   double S, S_prime, deltaS, tmp;
   double delta=mpars->delta;
   acc=0;

   for(i=0;i<V;i++)
   {
      ranlxd(rand,2); 
      S = action(apars); 
      tmp = phi[i];
      phi[i]+=delta*(1.*rand[0]-0.5);
      S_prime = action(apars);
      deltaS = S_prime - S;

      if (rand[1] < exp(-deltaS))   {acc++;}
      else {phi[i] = tmp;}
   }

   return acc;
}


