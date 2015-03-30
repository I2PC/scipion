FC=gfortran
FFLAGS=-O3

all:nma_diag_arpack

nma_diag_arpack:diag_arpack.o
	$(FC) $(FFLAGS) -m64 -o $@ diag_arpack.o -larpack -lblas -llapack
.f.o:
	$(FC) $(FFLAGS) -c $*.f -m64
