#!/bin/sh
PROCS=$1
shift
mpirun -np $PROCS $SCIPION_HOME/software/em/xmipp/bin/xmipp_mpi_classify_CLTomo_prog "$@"
