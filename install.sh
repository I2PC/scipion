#!/bin/sh

#Some flags variables

DO_UNTAR=true
DO_SQLITE=true
DO_TCLTK=true
DO_PYTHON=true
DO_PYMOD=true
DO_FFTW=true
DO_TIFF=true
DO_JPEG=true
DO_HDF5=true
DO_ARPACK=false
DO_FRM=false
DO_NMA=true

DO_CLEAN=true
DO_STATIC=false
DO_UPDATE=false
DO_SETUP=true
DO_GUI=true
DO_UNATTENDED=false
DO_MATLAB=false

export NUMBER_OF_CPU=1

# Some other vars

#################### PARSING PARAMETERS ###########################
TAKE_CPU=false
TAKE_CONFIGURE=false
TAKE_COMPILE=false
CONFIGURE_ARGS=""
COMPILE_ARGS=""
GUI_ARGS="gui"

OS_TYPE=`uname`
IS_MAC=false
IS_CYGWIN=false
IS_MINGW=false
IS_LINUX=false
echo "The OS is $OS_TYPE"
case "$OS_TYPE" in
  Darwin)
    IS_MAC=true
    CONFIGURE_ARGS="mpi=True MPI_CXX=mpic++ MPI_LINKERFORPROGRAMS=mpic++"
    ;;
  CYGWIN*)
    IS_CYGWIN=true
    ;;
  MINGW*)
    IS_MINGW=true
    ;;
  *)
    IS_LINUX=true
    ;;
esac

for param in $@; do
 if $TAKE_CPU; then
    NUMBER_OF_CPU=$param
    TAKE_CPU=false
 elif $TAKE_CONFIGURE && [ "$param" != "compile" ]; then
     echo "param: $param"
     CONFIGURE_ARGS="$CONFIGURE_ARGS $param"
 elif $TAKE_COMPILE && [ "$param" != "configure" ]; then
     COMPILE_ARGS="$COMPILE_ARGS $param"
 else
    case $param in
        "disable_all")
			    DO_UNTAR=false
			    DO_SQLITE=false
			    DO_TCLTK=false
			    DO_PYTHON=false
			    DO_PYMOD=false
			    DO_FFTW=false
			    DO_TIFF=false
			    DO_JPEG=false
			    DO_HDF5=false
			    DO_ARPACK=false
			    DO_NMA=false
			    DO_MATLAB=false
			    DO_SETUP=false;;			
        "-j")               TAKE_CPU=true;;
        "untar=true")       DO_UNTAR=true;;
        "untar=false")      DO_UNTAR=false;;
        "sqlite=true")      DO_SQLITE=true;;
        "sqlite=false")     DO_SQLITE=false;;
        "tcltk=true")       DO_TCLTK=true;;
        "tcltk=false")      DO_TCLTK=false;;
        "python=true")      DO_PYTHON=true;DO_PYMOD=true;;
        "python=false")     DO_PYTHON=false;DO_PYMOD=false;;
        "pymodules=true")   DO_PYMOD=true;;
        "pymodules=false")  DO_PYMOD=false;;
        "fftw=true")        DO_FFTW=true;;
        "fftw=false")       DO_FFTW=false;;
        "tiff=true")        DO_TIFF=true;;
        "tiff=false")       DO_TIFF=false;;
        "jpeg=true")        DO_JPEG=true;;
        "jpeg=false")       DO_JPEG=false;;
        "hdf5=true")        DO_HDF5=true;;
        "hdf5=false")       DO_HDF5=false;;
        "arpack=true")      DO_ARPACK=true;;
        "arpack=false")     DO_ARPACK=false;;
        "frm=true")         DO_FRM=true;;
        "frm=false")        DO_FRM=false;;
        "matlab=true")      DO_MATLAB=true;;
        "matlab=false")     DO_MATLAB=false;;
        "nma=true")         DO_NMA=true;;
        "nma=false")        DO_NMA=false;;
        "clean=true")       DO_CLEAN=true;;
        "clean=false")      DO_CLEAN=false;;
        "static=true")      DO_STATIC=true;;
        "static=false")     DO_STATIC=false;;
        "gui=false")        GUI_ARGS="";;
	"setup=true")       DO_SETUP=true;;
        # This two if passed should be at the end and 
        # will setup arguments for configure and compilation steps
        "configure")        TAKE_CONFIGURE=true;
                            TAKE_COMPILE=false;;
        "compile")          TAKE_CONFIGURE=false;
                            TAKE_COMPILE=true;;
        "unattended=true")  DO_UNATTENDED=true;;
        "unattended=false") DO_UNATTENDED=false;;
         *)          echo "Unrecognized option $param, exiting..."; exit 1
    esac
 fi 
done

if $DO_UNATTENDED; then
  CONFIGURE_ARGS="$CONFIGURE_ARGS unattended"
fi

#Some path variables
export XMIPP_HOME=$PWD
export PATH=$XMIPP_HOME/bin:$PATH
export LD_LIBRARY_PATH=$XMIPP_HOME/lib:$LD_LIBRARY_PATH
if $IS_MAC; then
	export DYLD_FALLBACK_LIBRARY_PATH=$XMIPP_HOME/lib:$DYLD_FALLBACK_LIBRARY_PATH
fi

# Create file to include from BASH this Xmipp installation
INC_FILE=.xmipp.bashrc
echo "export XMIPP_HOME=$PWD" > $INC_FILE
echo 'export PATH=$XMIPP_HOME/bin:$PATH' >> $INC_FILE
echo 'export LD_LIBRARY_PATH=$XMIPP_HOME/lib:$LD_LIBRARY_PATH' >> $INC_FILE
echo 'if [ "$BASH" != "" ]; then' >> $INC_FILE
echo '# Load global autocomplete file ' >> $INC_FILE
echo 'test -s $XMIPP_HOME/.xmipp.autocomplete && . $XMIPP_HOME/.xmipp.autocomplete || true' >> $INC_FILE
echo '# Load programs autocomplete file ' >> $INC_FILE
echo 'test -s $XMIPP_HOME/.xmipp_programs.autocomplete && . $XMIPP_HOME/.xmipp_programs.autocomplete || true' >> $INC_FILE
echo 'fi' >> $INC_FILE

if $IS_MAC; then
	echo 'export DYLD_FALLBACK_LIBRARY_PATH=$XMIPP_HOME/lib:$DYLD_FALLBACK_LIBRARY_PATH' >> $INC_FILE
fi
echo " "    >> $INC_FILE
echo " "    >> $INC_FILE

echo "# Xmipp Aliases 						 "    >> $INC_FILE
echo "## Setup ##                        "    >> $INC_FILE
echo "alias xconfigure='./setup.py -j $NUMBER_OF_CPU configure compile ' " >> $INC_FILE
echo "alias xcompile='./setup.py -j $NUMBER_OF_CPU compile ' "             >> $INC_FILE
echo "alias xupdate='./setup.py -j $NUMBER_OF_CPU update compile ' "       >> $INC_FILE
echo "## Interface ##                        "    >> $INC_FILE
echo "alias xa='xmipp_apropos'               "    >> $INC_FILE
echo "alias xb='xmipp_browser'               "    >> $INC_FILE
echo "alias xp='xmipp_protocols'             "    >> $INC_FILE
echo "alias xmipp='xmipp_protocols'          "    >> $INC_FILE
echo "alias xs='xmipp_showj'                 "    >> $INC_FILE
echo "alias xmipp_show='xmipp_showj'         "    >> $INC_FILE
echo "alias xsj='xmipp_showj'                "    >> $INC_FILE
echo "alias xij='xmipp_imagej'               "    >> $INC_FILE
echo "## Image ##                            "    >> $INC_FILE
echo "alias xic='xmipp_image_convert'        "    >> $INC_FILE
echo "alias xih='xmipp_image_header'         "    >> $INC_FILE
echo "alias xio='xmipp_image_operate'        "    >> $INC_FILE
echo "alias xis='xmipp_image_statistics'     "    >> $INC_FILE
echo "## Metadata ##                         "    >> $INC_FILE
echo "alias xmu='xmipp_metadata_utilities'   "    >> $INC_FILE
echo "alias xmp='xmipp_metadata_plot'        "    >> $INC_FILE
echo "## Transformation ##                   "    >> $INC_FILE
echo "alias xtg='xmipp_transform_geometry'   "    >> $INC_FILE
echo "alias xtf='xmipp_transform_filter'     "    >> $INC_FILE
echo "alias xtn='xmipp_transform_normalize'  "    >> $INC_FILE
echo "## Other ##                            "    >> $INC_FILE
echo "alias xrf='xmipp_resolution_fsc'       "    >> $INC_FILE
echo "alias xrs='xmipp_resolution_ssnr'      "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "## Configuration ##                                                          "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "# This file will serve to customize some settings of you Xmipp installation  "    >> $INC_FILE
echo "# Each setting will be in the form o a shell variable set to some value      "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "#---------- GUI ----------                                                   "    >> $INC_FILE
echo "# If you set to 1 the value of this variable, by default the programs        "    >> $INC_FILE
echo "# will launch the gui when call without argments, default is print the help  "    >> $INC_FILE
echo "export XMIPP_GUI_DEFAULT=0                                                   "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "# If you set to 0 the value of this variable the script generated            "    >> $INC_FILE
echo "# by programs gui will not be erased and you can use the same parameters     "    >> $INC_FILE
echo "export XMIPP_GUI_CLEAN=1                                                     "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "#---------- Parallel ----------                                              "    >> $INC_FILE
echo "# This variable will point to your job submition template file               "    >> $INC_FILE
echo "export XMIPP_PARALLEL_LAUNCH=config_launch.py                                "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "# If you have .xmipp.cfg in your home folder it will override                "    >> $INC_FILE
echo "# this configurations                                                        "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "test -s ~/.xmipp.cfg && . ~/.xmipp.cfg || true                               "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE

chmod a+x $INC_FILE


# for CSH or TCSH
INC_FILE=.xmipp.csh
echo "setenv XMIPP_HOME $PWD" > $INC_FILE
echo 'setenv PATH $XMIPP_HOME/bin:$PATH' >> $INC_FILE
echo 'if($?LD_LIBRARY_PATH) then' >> $INC_FILE
echo '  setenv LD_LIBRARY_PATH $XMIPP_HOME/lib:$LD_LIBRARY_PATH' >> $INC_FILE
echo 'else' >> $INC_FILE
echo '  setenv LD_LIBRARY_PATH $XMIPP_HOME/lib' >> $INC_FILE
echo 'endif' >> $INC_FILE

if $IS_MAC; then
  echo 'if($?DYLD_FALLBACK_LIBRARY_PATH) then' >> $INC_FILE
  echo '  setenv DYLD_FALLBACK_LIBRARY_PATH $XMIPP_HOME/lib:$DYLD_FALLBACK_LIBRARY_PATH' >> $INC_FILE
  echo 'else' >> $INC_FILE
  echo '  setenv DYLD_FALLBACK_LIBRARY_PATH $XMIPP_HOME/lib' >> $INC_FILE
  echo 'endif' >> $INC_FILE
fi
echo 'test -s $XMIPP_HOME/.xmipp.alias && source $XMIPP_HOME/.xmipp.alias || true' >> $INC_FILE

echo " "    >> $INC_FILE
echo " "    >> $INC_FILE

echo "# Xmipp Aliases 						 "    >> $INC_FILE
echo "## Setup ##                        "    >> $INC_FILE
echo "alias xconfigure './setup.py -j $NUMBER_OF_CPU configure compile ' " >> $INC_FILE
echo "alias xcompile './setup.py -j $NUMBER_OF_CPU compile ' "             >> $INC_FILE
echo "alias xupdate './setup.py -j $NUMBER_OF_CPU update compile ' "       >> $INC_FILE
echo "## Interface ##                        "    >> $INC_FILE
echo "alias xa 'xmipp_apropos'               "    >> $INC_FILE
echo "alias xb 'xmipp_browser'               "    >> $INC_FILE
echo "alias xp 'xmipp_protocols'             "    >> $INC_FILE
echo "alias xmipp 'xmipp_protocols'          "    >> $INC_FILE
echo "alias xs 'xmipp_showj'                 "    >> $INC_FILE
echo "alias xmipp_show 'xmipp_showj'         "    >> $INC_FILE
echo "alias xsj 'xmipp_showj'                "    >> $INC_FILE
echo "## Image ##                            "    >> $INC_FILE
echo "alias xic 'xmipp_image_convert'        "    >> $INC_FILE
echo "alias xih 'xmipp_image_header'         "    >> $INC_FILE
echo "alias xio 'xmipp_image_operate'        "    >> $INC_FILE
echo "alias xis 'xmipp_image_statistics'     "    >> $INC_FILE
echo "## Metadata ##                         "    >> $INC_FILE
echo "alias xmu 'xmipp_metadata_utilities'   "    >> $INC_FILE
echo "alias xmp 'xmipp_metadata_plot'        "    >> $INC_FILE
echo "## Transformation ##                   "    >> $INC_FILE
echo "alias xtg 'xmipp_transform_geometry'   "    >> $INC_FILE
echo "alias xtf 'xmipp_transform_filter'     "    >> $INC_FILE
echo "alias xtn 'xmipp_transform_normalize'  "    >> $INC_FILE
echo "## Other ##                            "    >> $INC_FILE
echo "alias xrf 'xmipp_resolution_fsc'       "    >> $INC_FILE
echo "alias xrs 'xmipp_resolution_ssnr'      "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "## Configuration ##                                                          "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "# This file will serve to customize some settings of you Xmipp installation  "    >> $INC_FILE
echo "# Each setting will be in the form o a shell variable set to some value      "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "#---------- GUI ----------                                                   "    >> $INC_FILE
echo "# If you set to 1 the value of this variable, by default the programs        "    >> $INC_FILE
echo "# will launch the gui when call without argments, default is print the help  "    >> $INC_FILE
echo "setenv XMIPP_GUI_DEFAULT 0                                                   "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "# If you set to 0 the value of this variable the script generated            "    >> $INC_FILE
echo "# by programs gui will not be erased and you can use the same parameters     "    >> $INC_FILE
echo "setenv XMIPP_GUI_CLEAN 1                                                     "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "#---------- Parallel ----------                                              "    >> $INC_FILE
echo "# This variable will point to your job submition template file               "    >> $INC_FILE
echo "setenv XMIPP_PARALLEL_LAUNCH config_launch.py                                "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "# If you have .xmipp.cfg in your home folder it will override                "    >> $INC_FILE
echo "# this configurations                                                        "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE
echo "test -s ~/.xmipp.cfg && source ~/.xmipp.cfg || true                          "    >> $INC_FILE
echo "                                                                             "    >> $INC_FILE

chmod a+x $INC_FILE



EXT_PATH=$XMIPP_HOME/external
BUILD_PATH=$XMIPP_HOME/build

#External libraries versions
VSQLITE=sqlite-3.6.23
VSQLITE_EXT=sqliteExt
VTCLTK=8.5.10
VPYTHON=Python-2.7.2
VFFTW=fftw-3.3.3
VTIFF=tiff-3.9.4
VJPEG=jpeg-8c
VHDF5=hdf5-1.8.10
VARPACK=arpack++-2.3
VNUMPY=numpy-1.6.1
VSCIPY=scipy-0.12.0
VMATLIBPLOT=matplotlib-1.1.0
VPYMPI=mpi4py-1.2.2
#read star files
VPYCIFRW=PyCifRW-3.3
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

GREEN="\033[32m"
RED="\033[31m"
ENDC="\033[0m"

# Print a green msg using terminal escaped color sequence
echoGreen()
{
    printf "$GREEN %b $ENDC\n" "$1"
}

# Execute and print the sequence
echoExec()
{
 COMMAND="$@"
 echo '-->' $COMMAND
 $COMMAND
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
  echoGreen "*** Compiling $LIB ..."
  echoExec "cd $_PATH"

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
  echo "--> make -j $NUMBER_OF_CPU >$BUILD_PATH/${LIB}_make.log 2>&1"
  make -j $NUMBER_OF_CPU >$BUILD_PATH/${LIB}_make.log 2>&1
  toc
}

compile_pymodule()
{
   MOD=$1
   _PATH=$EXT_PATH/python/$MOD
   #_PYTHON=$EXT_PATH/python/$VPYTHON/python
   echo "--> cd $_PATH"
   cd $_PATH
   echo "--> xmipp_python setup.py install --prefix $XMIPP_HOME >$BUILD_PATH/${MOD}_setup_install.log 2>&1"
   xmipp_python setup.py install --prefix $XMIPP_HOME >$BUILD_PATH/${MOD}_setup_install.log 2>&1 
   
}

#This function should be called from XMIPP_HOME
# Parameter: Library_Path Library_name Lib_Version_Number 
install_libs()
{
  cd $XMIPP_HOME
  LIBPATH=external/$1; shift
  COMMON="$1"; shift
  VERSION=$1; shift
  COPY=$1
  SUFFIXES=".a .la "
  if $IS_MAC; then
	SUFFIXES="$SUFFIXES .dylib .$VERSION.dylib"
  elif $IS_MINGW; then
        SUFFIXES="$SUFFIXES .dll.a -$VERSION.dll"
  else 
	SUFFIXES="$SUFFIXES .so .so.$VERSION"
  fi
  
  for suffix in $SUFFIXES; do
     LIBNAME=$COMMON$suffix
     if $COPY; then
     	     if [ -e lib/$LIBNAME ]; then
	         rm -f lib/$LIBNAME
             fi
	     echo "--> cp -f $LIBPATH/$LIBNAME lib/$LIBNAME"
	     cp -f $LIBPATH/$LIBNAME lib/$LIBNAME  
     else
	     echo "--> ln -sf ../$LIBPATH/$LIBNAME lib/$LIBNAME"
	     ln -sf ../$LIBPATH/$LIBNAME lib/$LIBNAME  
     fi
  done
}

install_bin()
{
  cd $XMIPP_HOME
  BINPATH=../external/$1
  LINKNAME=bin/$2
  echoExec "ln -sf $BINPATH $LINKNAME"
}

create_dir()
{
  DIR=$1
  if [ -d $DIR ]; then 
      echo "--> Dir $DIR exists."
  else
    echoExec "mkdir $DIR"
  fi
}

#################### NEEDED FOLDERS: bin lib build ##############
echoGreen "*** Checking needed folders ..."
create_dir build
create_dir bin
create_dir lib


#################### DECOMPRESSING EXTERNAL LIBRARIES ###########################
if $DO_UNTAR; then  
  tic
  dirs=". python"
  echo
  echoGreen "*** Decompressing external libraries ..."
  #Enter to external dir
  echo "--> cd $EXT_PATH"
  cd $EXT_PATH
  for d in $dirs; do
    echo "--> cd $d"
    cd $d 
    for file in $(ls *.tgz); do
      echo "--> tar -xvzf $file > /dev/null"
      tar -xvzf $file > /dev/null 2>&1
    done
    echo "--> cd -"
    cd - > /dev/null 2>&1
  done
  toc
fi

#################### SQLITE ###########################
if $DO_SQLITE; then
  if $IS_MAC; then
    #compile_library $VSQLITE "." "." "CPPFLAGS=-w CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1 -I/opt/local/include -I/opt/local/lib -I/sw/include -I/sw/lib -lsqlite3" ".libs"
    compile_library $VSQLITE "." "." "CPPFLAGS=-w CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1" ".libs"
  else
    compile_library $VSQLITE "." "." "CPPFLAGS=-w CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1" ".libs"
  fi
  #execute sqlite to avoid relinking in the future
  echo "select 1+1 ;" | $VSQLITE/sqlite3
  install_bin $VSQLITE/sqlite3 xmipp_sqlite3
  install_libs $VSQLITE/.libs libsqlite3 0 false
  #compile math library for sqlite
  cd $EXT_PATH/$VSQLITE_EXT
  if $IS_MINGW; then
    gcc -shared -I. -o libsqlitefunctions.dll extension-functions.c
    cp libsqlitefunctions.dll $XMIPP_HOME/lib/libXmippSqliteExt.dll
  elif $IS_MAC; then
    gcc -fno-common -dynamiclib extension-functions.c -o libsqlitefunctions.dylib
    cp libsqlitefunctions.dylib $XMIPP_HOME/lib/libXmippSqliteExt.dylib
  else  
    gcc -fPIC -lm -shared  extension-functions.c -o libsqlitefunctions.so
    cp libsqlitefunctions.so $XMIPP_HOME/lib/libXmippSqliteExt.so
  fi
fi

#################### FFTW ###########################
if $DO_FFTW; then
  if $IS_MINGW; then
    FFTWFLAGS=" CPPFLAGS=-I/c/MinGW/include CFLAGS=-I/c/MinGW/include"
  else
    FFTWFLAGS=""
  fi
  FLAGS="$FFTWFLAGS --enable-threads"
  compile_library $VFFTW "." "." $FLAGS
  install_libs $VFFTW/.libs libfftw3 3 true
  install_libs $VFFTW/threads/.libs libfftw3_threads 3 true

  FLAGS="$FFTWFLAGS --enable-float"
  compile_library $VFFTW "." "." $FLAGS
  install_libs $VFFTW/.libs libfftw3f 3 true
fi

#################### JPEG ###########################
if $DO_JPEG; then
compile_library $VJPEG "." "." "CPPFLAGS=-w"
install_libs $VJPEG/.libs libjpeg 8 false
fi

#################### TIFF ###########################
if $DO_TIFF; then
  compile_library $VTIFF "." "." "CPPFLAGS=-w --with-jpeg-include-dir=$EXT_PATH/$VJPEG --with-jpeg-lib-dir=$XMIPP_HOME/lib"
  install_libs $VTIFF/libtiff/.libs libtiff 3 false
fi

#################### HDF5 ###########################
if $DO_HDF5; then
compile_library $VHDF5 "." "." "CPPFLAGS=-w --enable-cxx"
install_libs $VHDF5/src/.libs libhdf5 7 false
install_libs $VHDF5/c++/src/.libs libhdf5_cpp 7 false
fi

#################### TCL/TK ###########################
if $DO_TCLTK; then
  if $IS_MAC; then
    compile_library tcl$VTCLTK python macosx "--disable-xft"
    compile_library tk$VTCLTK python macosx "--disable-xft"
  elif $IS_MINGW; then
    compile_library tcl$VTCLTK python win "--disable-xft CFLAGS=-I/c/MinGW/include CPPFLAGS=-I/c/MinGW/include"
    compile_library tk$VTCLTK python win "--disable-xft --with-tcl=../../tcl$VTCLTK/win CFLAGS=-I/c/MinGW/include CPPFLAGS=-I/c/MinGW/include"
  else
    compile_library tcl$VTCLTK python unix "--disable-xft --enable-threads"
    compile_library tk$VTCLTK python unix "--disable-xft --enable-threads"
  fi
fi

#################### PYTHON ###########################
EXT_PYTHON=$EXT_PATH/python

if $DO_PYTHON; then
  echoGreen "PYTHON SETUP"
  export CPPFLAGS="-I$EXT_PATH/$VSQLITE/ -I$EXT_PYTHON/tk$VTCLTK/generic -I$EXT_PYTHON/tcl$VTCLTK/generic"
  if $IS_CYGWIN; then
    export CPPFLAGS="-I/usr/include -I/usr/include/ncurses $CPPFLAGS"
    export LDFLAGS="-L$EXT_PYTHON/$VPYTHON -L$XMIPP_HOME/lib -L$EXT_PYTHON/tk$VTCLTK/unix -L$EXT_PYTHON/tcl$VTCLTK/unix"
    export LD_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/unix:$EXT_PYTHON/tcl$VTCLTK/unix:$LD_LIBRARY_PATH"
    echo "--> export CPPFLAGS=$CPPFLAGS"
    echo "--> export LDFLAGS=$LDFLAGS"
    echo "--> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH"	 
  elif $IS_MAC; then
    export LDFLAGS="-L$EXT_PYTHON/$VPYTHON -L$XMIPP_HOME/lib -L$EXT_PYTHON/tk$VTCLTK/macosx -L$EXT_PYTHON/tcl$VTCLTK/macosx"
    export LD_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/macosx:$EXT_PYTHON/tcl$VTCLTK/macosx:$LD_LIBRARY_PATH"
    export DYLD_FALLBACK_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/macosx:$EXT_PYTHON/tcl$VTCLTK/macosx:$DYLD_FALLBACK_LIBRARY_PATH"
    echo "--> export CPPFLAGS=$CPPFLAGS"
    echo "--> export LDFLAGS=$LDFLAGS"
    echo "--> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
    echo "--> export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH"
  elif $IS_MINGW; then
    export CPPFLAGS="-I/usr/include -I/usr/include/ncurses -I/c/MinGW/include $CPPFLAGS -D__MINGW32__ -I$EXT_PYTHON/$VPYTHON/Include -I$EXT_PYTHON/$VPYTHON/PC "
    export CFLAGS="$CPPFLAGS"
    export LDFLAGS="-L$EXT_PYTHON/$VPYTHON -L$XMIPP_HOME/lib -L$EXT_PYTHON/tk$VTCLTK/win -L$EXT_PYTHON/tcl$VTCLTK/win"
    export LD_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/win:$EXT_PYTHON/tcl$VTCLTK/win:$LD_LIBRARY_PATH"
    echo "--> export CPPFLAGS=$CPPFLAGS"
    echo "--> export LDFLAGS=$LDFLAGS"
    echo "--> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
  else
    export LDFLAGS="-L$EXT_PYTHON/$VPYTHON -L$XMIPP_HOME/lib -L$EXT_PYTHON/tk$VTCLTK/unix -L$EXT_PYTHON/tcl$VTCLTK/unix"
    export LD_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/unix:$EXT_PYTHON/tcl$VTCLTK/unix:$LD_LIBRARY_PATH"
    echo "--> export CPPFLAGS=$CPPFLAGS"
    echo "--> export LDFLAGS=$LDFLAGS"
    echo "--> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH"     
  fi
  echoGreen "Copying our custom python files ..."
  echoExec "cd $EXT_PYTHON"
  echoExec "cp ./xmipp_setup.py $VPYTHON/setup.py"
  echoExec "chmod a+x $VPYTHON/setup.py"
  #cp ./xmipp_setup.py $VPYTHON/setup.py
  #I thick these two are not needed
  #cp ./xmipp__iomodule.h $VPYTHON/Modules/_io/_iomodule.h
  #echo "--> cp ./xmipp__iomodule.h $VPYTHON/Modules/_io/_iomodule.h"
    
  compile_library $VPYTHON python "." ""

  # Create the python launch script with necessary environment variable settings
  PYTHON_BIN=$XMIPP_HOME/bin/xmipp_python
  echo "--> Creating python launch script $PYTHON_BIN ..."
  printf "#!/bin/sh\n\n" > $PYTHON_BIN
  if $IS_CYGWIN; then
    printf 'VPYTHON=%b \n' "$VPYTHON" >> $PYTHON_BIN
    printf 'VTCLTK=%b \n\n' "$VTCLTK" >> $PYTHON_BIN
    printf 'EXT_PYTHON=$XMIPP_HOME/external/python \n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tcl$VTCLTK/unix:$EXT_PYTHON/tk$VTCLTK/unix:$LD_LIBRARY_PATH \n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$XMIPP_HOME/applications/tests/pythonlib:$XMIPP_HOME/lib/python2.7/site-packages:$PYTHONPATH \n' >> $PYTHON_BIN
    printf 'export TCL_LIBRARY=$EXT_PYTHON/tcl$VTCLTK/library \n' >> $PYTHON_BIN
    printf 'export TK_LIBRARY=$EXT_PYTHON/tk$VTCLTK/library \n\n' >> $PYTHON_BIN
    printf 'PYTHONCYGWINLIB=`find $EXT_PYTHON/$VPYTHON/build -name "lib.cygwin*" -type d`\n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$PYTHONCYGWINLIB:$LD_LIBRARY_PATH\n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$PYTHONCYGWINLIB:$PYTHONPATH\n' >> $PYTHON_BIN
    printf '$EXT_PYTHON/$VPYTHON/python.exe "$@"\n' >> $PYTHON_BIN
  elif $IS_MINGW; then
    printf 'VPYTHON=Python27 \n' >> $PYTHON_BIN
    printf 'VTCLTK=8.5 \n\n' >> $PYTHON_BIN
    printf 'EXT_PYTHON=/c \n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/$VPYTHON/tcl/tcl$VTCLTK:$EXT_PYTHON/$VPYTHON/tcl/tk$VTCLTK:$LD_LIBRARY_PATH \n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$XMIPP_HOME/applications/tests/pythonlib:$EXT_PYTHON/$VPYTHON:$XMIPP_HOME/lib:$PYTHONPATH \n' >> $PYTHON_BIN
    printf 'export TCL_LIBRARY=$EXT_PYTHON/$VPYTHON/tcl/tcl$VTCLTK \n' >> $PYTHON_BIN
    printf 'export TK_LIBRARY=$EXT_PYTHON/$VPYTHON/tcl/tk$VTCLTK \n\n' >> $PYTHON_BIN
    printf '$EXT_PYTHON/$VPYTHON/python.exe "$@"\n' >> $PYTHON_BIN
  elif $IS_MAC; then
    printf 'VPYTHON=%b \n' "$VPYTHON" >> $PYTHON_BIN
    printf 'VTCLTK=%b \n\n' "$VTCLTK" >> $PYTHON_BIN
    printf 'EXT_PYTHON=$XMIPP_HOME/external/python \n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tcl$VTCLTK/macosx:$EXT_PYTHON/tk$VTCLTK/macosx:$LD_LIBRARY_PATH \n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$XMIPP_HOME/applications/tests/pythonlib:$XMIPP_HOME/lib/python2.7/site-packages:$PYTHONPATH \n' >> $PYTHON_BIN
    printf 'export TCL_LIBRARY=$EXT_PYTHON/tcl$VTCLTK/library \n' >> $PYTHON_BIN
    printf 'export TK_LIBRARY=$EXT_PYTHON/tk$VTCLTK/library \n\n' >> $PYTHON_BIN
    printf 'export DYLD_FALLBACK_LIBRARY_PATH=$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tcl$VTCLTK/macosx:$EXT_PYTHON/tk$VTCLTK/macosx:$DYLD_FALLBACK_LIBRARY_PATH \n' >> $PYTHON_BIN	
    printf '$EXT_PYTHON/$VPYTHON/python.exe "$@"\n' >> $PYTHON_BIN
  else
    printf 'VPYTHON=%b \n' "$VPYTHON" >> $PYTHON_BIN
    printf 'VTCLTK=%b \n\n' "$VTCLTK" >> $PYTHON_BIN
    printf 'EXT_PYTHON=$XMIPP_HOME/external/python \n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tcl$VTCLTK/unix:$EXT_PYTHON/tk$VTCLTK/unix:$LD_LIBRARY_PATH \n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$XMIPP_HOME/applications/tests/pythonlib:$XMIPP_HOME/lib/python2.7/site-packages:$PYTHONPATH \n' >> $PYTHON_BIN
    printf 'export TCL_LIBRARY=$EXT_PYTHON/tcl$VTCLTK/library \n' >> $PYTHON_BIN
    printf 'export TK_LIBRARY=$EXT_PYTHON/tk$VTCLTK/library \n\n' >> $PYTHON_BIN

    printf '$EXT_PYTHON/$VPYTHON/python "$@"\n' >> $PYTHON_BIN
  fi
  echoExec "chmod a+x $PYTHON_BIN"
  #make python directory accesible by anybody
  echoExec "chmod -R a+x $XMIPP_HOME/external/python/Python-2.7.2"

fi    

#################### PYTHON MODULES ###########################
if $DO_PYMOD; then
  compile_pymodule $VNUMPY
  export CPPFLAGS="-I$EXT_PATH/$VSQLITE/ -I$EXT_PYTHON/tk$VTCLTK/generic -I$EXT_PYTHON/tcl$VTCLTK/generic"
  if $IS_MAC; then
    export LDFLAGS="-L$EXT_PYTHON/$VPYTHON -L$XMIPP_HOME/lib -L$EXT_PYTHON/tk$VTCLTK/macosx -L$EXT_PYTHON/tcl$VTCLTK/macosx"
    export LD_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/macosx:$EXT_PYTHON/tcl$VTCLTK/macosx:$LD_LIBRARY_PATH"
    export DYLD_FALLBACK_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/macosx:$EXT_PYTHON/tcl$VTCLTK/macosx:$DYLD_FALLBACK_LIBRARY_PATH"
    echoExec "ln -s $XMIPP_HOME/bin/xmipp_python $XMIPP_HOME/bin/python2.7"
    echoExec "cd $EXT_PYTHON/$VMATLIBPLOT"
    echoExec "ln -s $XMIPP_HOME/bin/xmipp_python $XMIPP_HOME/bin/pythonXmipp" 
    echoExec "make -f make.osx clean"
    echoExec "make -f make.osx PREFIX=$XMIPP_HOME PYVERSION=Xmipp fetch deps mpl_install"
    echoExec "rm $XMIPP_HOME/bin/pythonXmipp"
    echoExec "rm $XMIPP_HOME/bin/python2.7"
  elif $IS_MINGW; then
    export LDFLAGS="-L$EXT_PYTHON/$VPYTHON -L$XMIPP_HOME/lib -L$EXT_PYTHON/tk$VTCLTK/win -L$EXT_PYTHON/tcl$VTCLTK/win"
    export LD_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/win:$EXT_PYTHON/tcl$VTCLTK/win:$LD_LIBRARY_PATH"
    echoExec "ln -s $XMIPP_HOME/bin/xmipp_python $XMIPP_HOME/bin/python2.7"
    echoExec "cd $EXT_PYTHON/$VMATLIBPLOT"
    echoExec "ln -s $XMIPP_HOME/bin/xmipp_python $XMIPP_HOME/bin/pythonXmipp"
  else
    export LDFLAGS="-L$EXT_PYTHON/$VPYTHON -L$XMIPP_HOME/lib -L$EXT_PYTHON/tk$VTCLTK/unix -L$EXT_PYTHON/tcl$VTCLTK/unix"
    export LD_LIBRARY_PATH="$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tk$VTCLTK/unix:$EXT_PYTHON/tcl$VTCLTK/unix:$LD_LIBRARY_PATH"
    echoExec "cp $EXT_PYTHON/matplotlib_setupext.py $EXT_PYTHON/$VMATLIBPLOT/setupext.py"
    #The following is needed from matplotlib to works
    echoExec "cd $EXT_PYTHON/tk$VTCLTK/unix/"
    echoExec "ln -sf libtk8.5.so  libtk.so"
    echoExec "cd $EXT_PYTHON/tcl$VTCLTK/unix/"
    echoExec "ln -sf libtcl8.5.so  libtcl.so"
  fi
  compile_pymodule $VMATLIBPLOT
  compile_pymodule $VPYMPI
  compile_pymodule $VPYCIFRW
  
  if $DO_FRM; then
    # Fast Rotational Matching
    export LDFLAGS="-shared $LDFLAGS"
    compile_pymodule $VSCIPY
    cd $EXT_PATH/sh_alignment
    ./compile.sh
  fi
fi

#################### ARPACK ###########################
if $DO_ARPACK; then
  compile_library $VARPACK "." "." ""
  install_libs $VARPACK/src/.libs libarpack++ 2 false
fi

# Launch the configure/compile python script 
cd $XMIPP_HOME

#echoGreen "Compiling XMIPP ..."
#echoGreen "CONFIGURE: $CONFIGURE_ARGS"
#echoGreen "COMPILE: $COMPILE_ARGS"
#echoGreen "GUI: $GUI_ARGS"

if $DO_SETUP; then
	echoExec "./setup.py -j $NUMBER_OF_CPU configure $CONFIGURE_ARGS compile $COMPILE_ARGS $GUI_ARGS install"
fi

#################### NMA ###########################
if $DO_NMA; then
    echoExec "cd $XMIPP_HOME/external/NMA/ElNemo"
    echoExec "make" 
    echoExec "cp nma_* $XMIPP_HOME/bin"
    echoExec "cd $XMIPP_HOME/external/NMA/NMA_cart"
    echoExec "make" 
    echoExec "cp nma_* $XMIPP_HOME/bin"
    echoExec "cd $XMIPP_HOME"
    echoExec "cp $XMIPP_HOME/external/NMA/nma_* $XMIPP_HOME/bin"
    echoExec "cp $XMIPP_HOME/external/NMA/m_inout_Bfact.py $XMIPP_HOME/bin"
fi

exit 0

#################### JAVA ###########################
install_jdk()
{
  LINUX="http://download.oracle.com/otn-pub/java/jdk/6u27-b07/jdk-6u27-linux-i586.bin"
  LINUX64="http://download.oracle.com/otn-pub/java/jdk/6u27-b07/jdk-6u27-linux-x64.bin"
  MACOSX="PENDING"

  # Which is our OS?
  case "$(uname -s)" in
  Darwin)
	  JDK_URL=$MACOSX;;
  Linux)
	  case "$(uname -m)" in
		  x86_64) JDK_URL=$LINUX64;;
		  *) JDK_URL=$LINUX;;
	  esac;;
  esac
  # Download jdk from Oracle site
  echo "wget $JDK_URL -o /dev/null -O $EXT_PATH/java/jdk.bin"
  wget $JDK_URL -o /dev/null -O $EXT_PATH/java/jdk.bin
  # Install jdk
  $EXT_PATH/java/jdk.bin
}

exit 0
if $DO_JAVA; then
  echo -e "$GREEN*** Checking jdk ...$ENDC"  
  JAVA_HOME=$EXT_PATH/java/jvm
  
  if [ `which java`) ]; then 
      JAVAC=$(readlink -f `which javac`)
      echo "Java found at: $JAVAC"
      JDK_PATH=$(dirname `dirname $JAVAC`)
  else
    echo "Java jdk folder not found"
    read -p "Do you want to install java-jdk(inside xmipp, doesn't require admin privileges)? (Y/n)" ANSWER
    if [ $ANSWER -eq 'y' ]; then
       install_jdk
       JDK_PATH=$EXT_PATH/java/jdk1.6.0_27
    fi
    ln -sf $JDK_PATH $JAVA_HOME 
  fi
fi
