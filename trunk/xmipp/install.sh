#!/bin/bash

#Some flags variables

DO_UNTAR=false
DO_SQLITE=false
DO_TCLTK=false
DO_PYTHON=true
DO_FFTW=false
DO_TIFF=false
DO_ARPACK=true
DO_JAVA=false

DO_CLEAN=false
DO_STATIC=false
DO_DOWNLOAD=false

CPU=$@



# Some other vars

GREEN="\033[32m"
RED="\033[31m"
ENDC="\033[0m"

#################### PARSING PARAMETERS ###########################

#for p in $@; do
#   echo "param: $p" 
#done
#exit 1

#Some path variables
XMIPP_HOME=$PWD
EXT_PATH=$XMIPP_HOME/external
BUILD_PATH=$XMIPP_HOME/build

#External libraries versions
VSQLITE=sqlite-3.6.23
VTCLTK=8.5.10
VPYTHON=Python-2.7.2
VFFTW=fftw-3.2.2
VTIFF=tiff-3.9.4
VARPACK=arpack++-2.3
VNUMPY=numpy-1.6.1
VMATLIBPLOT=matplotlib-1.0.1

################# HELPER FUNCTIONS ##################
TIMESTAMP=""
tic()
{
   TIMESTAMP="$(date +%s)"
}
toc()
{
   NOW="$(date +%s)"
   ELAPSED="$(expr $NOW - $TIMESTAMP)"
   echo "*** Elapsed time: $ELAPSED seconds"
}

compile_library()
{

   tic
   LIB=$1
   PREFIX_PATH=$2
   SUFFIX_PATH=$3
   CONFIGFLAGS=$4
   LIBS_PATH=$5
   _PATH=$EXT_PATH/$PREFIX_PATH/$LIB/$SUFFIX_PATH
  echo
  echo -e "$GREEN*** Compiling $LIB ...$ENDC"
  echo "--> cd $_PATH"
  cd $_PATH

 if ! $DO_STATIC; then
	echo "--> Enabling shared libraries..."
	CONFIGFLAGS="--enable-shared $CONFIGFLAGS"
 fi

  if $DO_CLEAN; then
    echo "--> make distclean > /dev/null 2>&1"
    make distclean > /dev/null 2>&1
  fi

  echo "--> ./configure $CONFIGFLAGS >$BUILD_PATH/${LIB}_configure.log 2>&1"
  ./configure $CONFIGFLAGS >$BUILD_PATH/${LIB}_configure.log 2>&1
  echo "--> make $CPU >$BUILD_PATH/${LIB}_make.log 2>&1"
  make $CPU >$BUILD_PATH/${LIB}_make.log 2>&1
  toc
}

compile_pymodule()
{
   MOD=$1
   _PATH=$EXT_PATH/python/$MOD
   _PYTHON=$EXT_PATH/python/$VPYTHON/python
   echo "--> cd $_PATH"
   cd $_PATH
   echo "--> $_PYTHON setup.py install --prefix $XMIPP_HOME >$BUILD_PATH/${MOD}_setup_install.log 2>&1"
   $_PYTHON setup.py install --prefix $XMIPP_HOME >$BUILD_PATH/${MOD}_setup_install.log 2>&1 
   
}

#This function should be called from XMIPP_HOME
install_libs()
{
  cd $XMIPP_HOME
  LIBPATH=../external/$1; shift
  COMMON=$1; shift
  SUFFIXES=$@
  for suffix in $SUFFIXES; do
     LIBNAME=$COMMON$suffix
     echo "--> ln -sf $LIBPATH/$LIBNAME lib/$LIBNAME"
     ln -sf $LIBPATH/$LIBNAME lib/$LIBNAME  
  done
}

install_bin()
{
  cd $XMIPP_HOME
  BINPATH=../external/$1
  LINKNAME=bin/$2
  echo "--> ln -sf $BINPATH $LINKNAME"
  ln -sf $BINPATH $LINKNAME
}

create_dir()
{
  DIR=$1
  if [ -d $DIR ]; then 
      echo "--> Dir $DIR exists."
  else
    echo "--> mkdir $DIR"
    mkdir $DIR
  fi
}

#################### NEEDED FOLDERS: bin lib build ##############
echo -e "$GREEN*** Checking needed folders ...$ENDC"
create_dir build
create_dir bin
create_dir lib

#################### DECOMPRESSING EXTERNAL LIBRARIES ###########################
if $DO_UNTAR; then  
  tic
  dirs=". python"
  echo
  echo -e "$GREEN*** Decompressing external libraries ...$ENDC"
  #Enter to external dir
  echo "--> cd $EXT_PATH"
  cd $EXT_PATH
  for d in $dirs; do
    echo "--> cd $d"
    cd $d 
    for file in $(ls *.tgz); do
      echo "--> tar -xvzf $file > /dev/null"
      tar -xvzf $file > /dev/null
    done
    echo "--> cd -"
    cd - > /dev/null
  done
  toc
fi

#################### SQLITE ###########################
if $DO_SQLITE; then
  compile_library $VSQLITE "." "." "CPPFLAGS=-w CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1" ".libs"
  install_bin $VSQLITE/sqlite3 xmipp_sqlite3
  install_libs $VSQLITE/.libs libsqlite3. a la so so.0
fi

#################### TCL/TK ###########################
if $DO_TCLTK; then
  compile_library tcl$VTCLTK python unix ""
  compile_library tk$VTCLTK python unix ""
install_libs python/tcl$VTCLTK/unix libtcl8.5. so
install_libs python/tk$VTCLTK/unix libtk8.5. so
fi

#################### PYTHON ###########################
if $DO_PYTHON; then
	STATIC_BACKUP=$DO_STATIC
	DO_STATIC=true
 export CPPFLAGS="-I$EXT_PATH/$VSQLITE/ -I$EXT_PATH/python/tk$VTCLTK/generic -I$EXT_PATH/python/tcl$VTCLTK/generic"
 compile_library $VPYTHON python "." ""
  install_bin python/$VPYTHON/python xmipp_python
DO_STATIC=$STATIC_BACKUP
#  install_libs python/$VPYTHON libpython2.7. a so so.1.0
  compile_pymodule $VNUMPY
  compile_pymodule $VMATLIBPLOT
fi

#################### FFTW ###########################
if $DO_FFTW; then
  compile_library $VFFTW "." "." "--enable-threads"
  install_libs $VFFTW/.libs libfftw3. a la so so.3
  install_libs $VFFTW/threads/.libs libfftw3_threads. a la so so.3
fi

#################### TIFF ###########################
if $DO_TIFF; then
  compile_library $VTIFF "." "." "CPPFLAGS=-w"
  install_libs $VTIFF/libtiff/.libs libtiff. a la so so.3
fi

#################### ARPACK ###########################
if $DO_ARPACK; then
  compile_library $VARPACK "." "." ""
  install_libs $VARPACK/src/.libs libarpack++. a la so so.2
fi
exit 0
#################### JAVA ###########################
if $DO_JAVA; then
  JAVAC=$(readlink -f `which javac`)
  if [ -e $JAVAC ]; then 
      echo "Java found at: $JAVAC"
  else
    echo "Java jdk folder not found"
    read -p "Do you want to install java-jdk(inside xmipp, doesn't require admin privileges)? (Y/n)" $
  fi
fi
