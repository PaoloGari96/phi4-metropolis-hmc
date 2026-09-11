
/*
 *  File hmc.c
 *  Hybrid Montecarlo algorithm for the phi4 theory
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

double hamiltonian(act_params_t *apars)
{
   int i;
   double sum;
   sum=0;
   for (i=0;i<V;i++) sum+=(mom[i]*mom[i]);
   sum*=0.5;
   sum+=action(apars);

   return sum;
}

void move_pi(act_params_t *apars, double deltat)
{
   double F[V];
   double kappa=apars->kappa;
   double lambda=apars->lambda;
   double phin,phi2;
   int i,mu;

   for(i=0;i<V;i++){

      phin=0;
      for (mu=0;mu<D;mu++) phin+=(phi[hop[i][mu]] + phi[hop[i][D+mu]]);
      phi2=phi[i]*phi[i];
      F[i]=-2*kappa*phin + 2*phi[i] + 4*lambda*(phi2-1.0)*phi[i];
   }

   for(i=0;i<V;i++){
      mom[i]+=(-1)*deltat*F[i];
   }

}

void move_phi(act_params_t *apars, double deltat)
{
   int i;

   for(i=0;i<V;i++){
      phi[i]+=deltat*mom[i];
   }

}

void Leapfrog(act_params_t *apars, double deltat, int N)
{
   int i;

   for (i=0;i<N;i++){
      move_pi(apars,0.5*deltat);
      move_phi(apars,deltat);
      move_pi(apars,0.5*deltat);
   }

}

/********************************************************************** 
 * hybrid montecarlo
 **********************************************************************/

double hmc(act_params_t *apars, metro_params_t *mpars, double deltat, int N)
{
   int i;
   static double rand[1];
   double H, H_prime,deltaH;
   double phi_old[V];
   double x1[V],x2[V];
   double acc;

   acc=0;
   ranlxd(rand,1);
   ranlxd(x1,V);
   ranlxd(x2,V);

   for(i=0;i<V;i++){
      mom[i] = (sqrt(-2*log(1-x1[i]))) * cos(2*M_PI*(1-x2[i]));
   }
   
   H = hamiltonian(apars);

   for(i=0;i<V;i++)
   {
      phi_old[i] = phi[i];
   }

   Leapfrog(apars,deltat,N);

   H_prime = hamiltonian(apars);
   deltaH = H_prime - H;

	if (rand[0]<exp(-deltaH)) {acc++;}
	else 
   {
      for(i=0;i<V;i++)  {phi[i] = phi_old[i];}
   }

   return deltaH;
}


