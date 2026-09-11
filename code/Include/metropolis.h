#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "phi4.h"

/* METROPOLIS_C */
extern double action(act_params_t *apars);
extern double hamiltonian(act_params_t *apars);
extern void move_pi(act_params_t *apars, double deltat);
extern void move_phi(act_params_t *apars, double deltat);
extern void Leapfrog(act_params_t *apars, double deltat, int N);
extern double metropolis(act_params_t *apars, metro_params_t *mpars);
extern double hmc(act_params_t *apars, metro_params_t *mpars, double deltat, int N);

#endif
