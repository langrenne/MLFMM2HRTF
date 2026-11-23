# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
#
# mlfmm: MLFMM/BEM modules of MLFMM2HRTF (called by ../makefile).

FC = gfortran
ARCH ?= -march=x86-64-v3
FCFLAGS = -fopenmp -g -O3 $(ARCH)
FFLAGS = -Wall
LAPACK ?= -llapack -lblas
LDFLAGS = $(LAPACK)

all : module__Precision.o module__Constants.o sparse_kit.o module_gmres.o module_fft.o \
	module_sort.o parametres.o module_Utils.o module_Acoustiques.o module_Octree.o \
	module_BEM.o module_GMSH.o module_MLFMM.o

# Dependencies between modules (Fortran compilation order)
module__Constants.o  : module__Precision.o
parametres.o         : module__Precision.o module__Constants.o
# FFT: modern_fftpack (fftpack_*.f90, MIT + UCAR BSD-type, LICENSE_fftpack.md) through module_fft
fftpack_auxiliary_routines.o : fftpack_precision.o
fftpack_error_handlers.o : fftpack_auxiliary_routines.o
fftpack_complex_transform_routines.o : fftpack_precision.o fftpack_auxiliary_routines.o
fftpack_complex_initialization_routines.o fftpack_multiple_complex_forward.o fftpack_multiple_complex_backward.o : fftpack_complex_transform_routines.o
module_fft.o         : fftpack_complex_transform_routines.o fftpack_complex_initialization_routines.o \
	fftpack_multiple_complex_forward.o fftpack_multiple_complex_backward.o fftpack_error_handlers.o
module_Utils.o       : parametres.o
module_Acoustiques.o : module__Precision.o module__Constants.o module_Utils.o
module_Octree.o      : module__Precision.o
module_BEM.o         : module__Precision.o module_Octree.o sparse_kit.o module_Utils.o parametres.o
module_GMSH.o        : module_BEM.o module_Utils.o module_sort.o
module_MLFMM.o       : parametres.o module_Acoustiques.o module_Octree.o module_BEM.o module_fft.o
module_gmres.o       : sparse_kit.o

%.o : %.f08
	$(FC) $(FCFLAGS) $(FFLAGS) -c -Wno-tabs $< $(LDFLAGS) -o $@
%.o : %.f90
	$(FC) $(FCFLAGS) $(FFLAGS) -fallow-argument-mismatch -c $< $(LDFLAGS) -o $@

clean :
	rm -f *.o *.mod *.smod

.PHONY : all clean
