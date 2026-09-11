# Metropolis and HMC algorithms

Metropolis and Hybrid Monte Carlo simulations of the lattice phi^4 scalar
field theory (magnetization, susceptibility, Binder cumulant, autocorrelation
/ jackknife error analysis, spontaneous symmetry breaking). Written for the
Fisica Computazionale course.

Full writeup, derivations and plots: `RelazioneMetropolis&HMC.pdf`.

## Layout

```
code/Include/          headers (lattice geometry, RNG, parameter structs)
code/Modules/          shared building blocks
  geometry/hopping.c        nearest-neighbor table for the lattice
  random/ranlxd.c           ranlxd RNG (M. Luescher)
  extras/observables.c      magnetization observables (M, |M|, M^2, M^4)
  metropolis/metropolis.c   local Metropolis update + action
  hmc/hmc.c                 HMC update (leapfrog + accept/reject) + action
code/MainPrograms/     one .c file per executable (see below)
```

Note: `metropolis.c` and `hmc.c` both define `action()`, so a given binary
links against exactly one of the two — never both.

| Binary               | Source                                                | Purpose |
|----------------------|--------------------------------------------------------|---------|
| `metropolis`          | `MainPrograms/Metropolis/MainMetropolis.c`             | Metropolis sampling of the phi^4 theory |
| `ssb`                 | `MainPrograms/HybridMonteCarlo/NumericalSimulationsProgram/ssb.c` | HMC sampling, spontaneous symmetry breaking study |
| `exact-algorithm`      | `.../AlgorithmCheck/ExactAlgorithm.c`                  | Check ⟨exp(-ΔH)⟩ = 1 |
| `inexact-algorithm`    | `.../AlgorithmCheck/InexactAlgorithm.c`                | HMC without the accept/reject step |
| `hcons`                | `.../AlgorithmCheck/hcons.c`                           | Hamiltonian conservation along a trajectory |
| `phase-space`          | `.../AlgorithmCheck/PhaseSpace.c`                      | Phase-space measure preservation |
| `reversibility`        | `.../AlgorithmCheck/Reversibility.c`                   | Leapfrog time-reversibility check |

## Build

```sh
make                  # builds every program into bin/
make bin/metropolis    # build a single program
make clean
```

Requires a C compiler (`gcc` by default, see `Makefile`) and `libm`.

## Run

Every program takes one input file with whitespace-separated `key value`
pairs and writes its results to stdout:

```sh
./bin/metropolis input.in > out.dat
```

`metropolis` needs: `kappa lambda ntherm nsweep delta seed naccu`.
The HMC-based programs (`ssb`, `exact-algorithm`, `inexact-algorithm`,
`hcons`, `phase-space`, `reversibility`) additionally need: `ntraj N`.

Example `input.in` (values used in the report for the first Metropolis
study):

```
kappa 0.15
lambda 1.145
ntherm 100000
nsweep 100000000
delta 0.25
seed 42
naccu 10000
```

For an HMC program, add e.g. `ntraj 1000` and `N 10` (number of leapfrog
steps per trajectory).

## Bug fixes applied

- `ssb.c` read the `ntraj` and `N` input values into `metro_params.ntraj`
  and `metro_params.N`, which are not fields of `metro_params_t` — this
  did not compile. Fixed to write into `hmc_params.ntraj` and
  `leapfrog_params.N`, matching every other HMC main program.
- `InexactAlgorithm.c` read `ntraj`/`N` into `metro_params.ntherm`/
  `metro_params.nsweep`, silently overwriting the thermalization/sweep
  counts and leaving the leapfrog step count `N` at 0 (so `dt = tzero/N`
  evaluated to infinity). Fixed the same way as above.
