#!/bin/sh
PROCS=$1
shift       
PYTHONPATH=$SCIPION_HOME/software/lib/python2.7/site-packages/sh_alignment:$PYTHONPATH
mpirun -np $PROCS $SCIPION_HOME/software/em/xmipp/bin/xmipp_mpi_classify_CLTomo_prog "$@"
