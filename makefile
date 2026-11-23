# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
#
# MLFMM2HRTF: builds main_HRTF_Ray and main_HRTF_Diff.
#   make                                   both solvers
#   make Ray | make Diff                   one solver
#   make ARCH=-march=native                processor older than ~2013 or ARM
#   make LAPACK="-lopenblas"               another LAPACK/BLAS
#   make clean                             before rebuilding with other options

.SUFFIXES:
.DEFAULT_GOAL := all

FC = gfortran
# -march=x86-64-v3 (AVX2 + FMA, processors since ~2013-2015): -31 % computing time
ARCH ?= -march=x86-64-v3
FCFLAGS = -fopenmp -g -O3 $(ARCH)
LAPACK ?= -llapack -lblas

SRC_DIRECTORY = ./mlfmm
LDFLAGS = $(LAPACK) -I$(SRC_DIRECTORY)

OBJS_F08 := $(patsubst %.f08, %.o, $(wildcard $(SRC_DIRECTORY)/*.f08))
OBJS_F90 := $(patsubst %.f90, %.o, $(wildcard $(SRC_DIRECTORY)/*.f90))

# The drivers "use" the modules of ./mlfmm: rebuild them when a module source changes
MODULE_SRCS := $(wildcard $(SRC_DIRECTORY)/*.f08) $(wildcard $(SRC_DIRECTORY)/*.f90)
TOP_OBJS := main_HRTF_Diff.o main_HRTF_Ray.o routines_Fait_tout.o routines_gmres.o
$(TOP_OBJS) : $(MODULE_SRCS)

PROGRAMS = main_HRTF_Diff main_HRTF_Ray

all: modules $(PROGRAMS)

Ray: modules main_HRTF_Ray

Diff: modules main_HRTF_Diff

main_HRTF_Ray : main_HRTF_Ray.o routines_Fait_tout.o routines_gmres.o $(OBJS_F08) $(OBJS_F90)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

main_HRTF_Diff : main_HRTF_Diff.o routines_Fait_tout.o routines_gmres.o $(OBJS_F08) $(OBJS_F90)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

modules:
	$(MAKE) -C $(SRC_DIRECTORY) ARCH="$(ARCH)" LAPACK="$(LAPACK)" all

%.o: %.f08
	$(FC) $(FCFLAGS) -fallow-argument-mismatch -c $< $(LDFLAGS)

clean:
	-rm -f *.mod *.o $(PROGRAMS)
	$(MAKE) -C $(SRC_DIRECTORY) clean

.PHONY: all clean Ray Diff modules
