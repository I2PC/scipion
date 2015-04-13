FC=gfortran
FFLAGS=-O3

all:nma_check_modes nma_diagrtb nma_elnemo_pdbmat

nma_check_modes:check_modes.o
	$(FC) $(FFLAGS) -m64 -o $@ check_modes.o
nma_diagrtb:diagrtb.o
	$(FC) $(FFLAGS) -m64 -o $@ diagrtb.o
nma_elnemo_pdbmat:pdbmat_aux.o
# grep -v '^DoC' pdbmat.f > pdbmat_aux.f
	$(FC) $(FFLAGS) -m64 -o $@ pdbmat_aux.o
.f.o:
	$(FC) $(FFLAGS) -c $*.f -m64
