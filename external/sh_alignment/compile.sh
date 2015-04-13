#!/bin/sh

swig -python frm/swig/frm.i 

CCFLAGS="-O3 -fPIC -I. -c"
PYTHONDIR=../python/Python-2.7.2
PYTHONLIB=python2.7
PYTHONPACKAGES=../../lib/python2.7/site-packages
NUMPYDIR=$PYTHONPACKAGES/numpy/core

gcc $CCFLAGS ./frm/src/lib_vec.c
gcc $CCFLAGS ./frm/src/lib_std.c
gcc $CCFLAGS ./frm/src/lib_vio.c
gcc $CCFLAGS ./frm/src/lib_pwk.c
gcc $CCFLAGS ./frm/src/lib_vwk.c
gcc $CCFLAGS ./frm/src/lib_tim.c
gcc $CCFLAGS ./frm/src/lib_pio.c
gcc $CCFLAGS ./frm/src/lib_eul.c
gcc $CCFLAGS ./SpharmonicKit27/OURmods.c
gcc $CCFLAGS ./SpharmonicKit27/primitive_FST.c
gcc $CCFLAGS ./SpharmonicKit27/MathFace.c
gcc $CCFLAGS ./SpharmonicKit27/fft_grids.c
gcc $CCFLAGS ./SpharmonicKit27/permroots.c
gcc $CCFLAGS ./SpharmonicKit27/weights.c
gcc $CCFLAGS ./SpharmonicKit27/newFCT.c
gcc $CCFLAGS ./SpharmonicKit27/FST_semi_memo.c
gcc $CCFLAGS ./SpharmonicKit27/indextables.c
gcc $CCFLAGS ./SpharmonicKit27/FFTcode.c
gcc $CCFLAGS ./SpharmonicKit27/naive_synthesis.c
gcc $CCFLAGS ./SpharmonicKit27/OURperms.c
gcc $CCFLAGS ./SpharmonicKit27/seminaive.c
gcc $CCFLAGS ./SpharmonicKit27/primitive.c
gcc $CCFLAGS ./SpharmonicKit27/cospmls.c
gcc $CCFLAGS ./SpharmonicKit27/oddweights.c
gcc $CCFLAGS ./SpharmonicKit27/csecond.c

CCFLAGS="-O3 -fPIC -I$PYTHONDIR/Include -I$NUMPYDIR/include -c -I$PYTHONDIR"
gcc $CCFLAGS ./frm/swig/frm_wrap.c
CCFLAGS="-O3 -fPIC -I. -Ifrm/src -I$PYTHONDIR/Include -ISpharmonicKit27 -I$NUMPYDIR/include -c -L../../lib "
gcc $CCFLAGS ./frm/swig/frm.c

gcc -shared -o _swig_frm.so frm.o permroots.o csecond.o frm_wrap.o indextables.o weights.o FST_semi_memo.o oddweights.o lib_pwk.o newFCT.o lib_eul.o lib_vio.o OURmods.o lib_vec.o lib_pio.o seminaive.o MathFace.o FFTcode.o OURperms.o primitive.o fft_grids.o naive_synthesis.o lib_std.o cospmls.o lib_vwk.o lib_tim.o primitive_FST.o ../../lib/libfftw3.so
rm *.o

if [ ! -d $PYTHONPACKAGES/sh_alignment ]
then
    mkdir $PYTHONPACKAGES/sh_alignment
fi
cp *.py $PYTHONPACKAGES/sh_alignment
cp ./frm/swig/swig_frm.py _swig_frm.so $PYTHONPACKAGES/sh_alignment

cp -r tompy $PYTHONPACKAGES/sh_alignment
