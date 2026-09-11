# Builds the Metropolis and HMC programs for the phi^4 lattice theory.
# See Relazione/RelazioneMetropolis&HMC.pdf for the physics writeup.
#
# Each main program links against either the Metropolis module or the HMC
# module (both define an `action()` function, so they must not be linked
# together into the same binary).
#
# Usage:
#   make                 build every program into bin/
#   make bin/metropolis  build a single program
#   make clean           remove bin/

CC      := gcc
CFLAGS  := -O2 -Wall -ICodici/Include
LDLIBS  := -lm

BIN_DIR := bin
MAIN_DIR := Codici/MainPrograms
MOD_DIR  := Codici/Modules

COMMON_SRCS     := $(MOD_DIR)/geometry/hopping.c $(MOD_DIR)/random/ranlxd.c $(MOD_DIR)/extras/observables.c
METROPOLIS_SRCS := $(MOD_DIR)/metropolis/metropolis.c
HMC_SRCS        := $(MOD_DIR)/hmc/hmc.c

.PHONY: all clean

all: $(BIN_DIR)/metropolis \
     $(BIN_DIR)/ssb \
     $(BIN_DIR)/exact-algorithm \
     $(BIN_DIR)/inexact-algorithm \
     $(BIN_DIR)/hcons \
     $(BIN_DIR)/phase-space \
     $(BIN_DIR)/reversibility

$(BIN_DIR)/metropolis: $(MAIN_DIR)/Metropolis/MainMetropolis.c $(COMMON_SRCS) $(METROPOLIS_SRCS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BIN_DIR)/ssb: $(MAIN_DIR)/HybridMonteCarlo/NumericalSimulationsProgram/ssb.c $(COMMON_SRCS) $(HMC_SRCS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BIN_DIR)/exact-algorithm: $(MAIN_DIR)/HybridMonteCarlo/AlgorithmCheck/ExactAlgorithm.c $(COMMON_SRCS) $(HMC_SRCS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BIN_DIR)/inexact-algorithm: $(MAIN_DIR)/HybridMonteCarlo/AlgorithmCheck/InexactAlgorithm.c $(COMMON_SRCS) $(HMC_SRCS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BIN_DIR)/hcons: $(MAIN_DIR)/HybridMonteCarlo/AlgorithmCheck/hcons.c $(COMMON_SRCS) $(HMC_SRCS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BIN_DIR)/phase-space: $(MAIN_DIR)/HybridMonteCarlo/AlgorithmCheck/PhaseSpace.c $(COMMON_SRCS) $(HMC_SRCS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BIN_DIR)/reversibility: $(MAIN_DIR)/HybridMonteCarlo/AlgorithmCheck/Reversibility.c $(COMMON_SRCS) $(HMC_SRCS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR)
