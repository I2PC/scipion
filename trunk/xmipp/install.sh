#!/bin/sh

#Some flags variables

DO_UNTAR=true
DO_SQLITE=true
DO_TCLTK=true
DO_PYTHON=true
DO_FFTW=true
DO_TIFF=true
DO_ARPACK=false

DO_CLEAN=true
DO_STATIC=false
DO_GUI=true

export NUMBER_OF_CPU=1



# Some other vars


#################### PARSING PARAMETERS ###########################
TAKE_CPU=false
TAKE_CONFIGURE=false
TAKE_COMPILE=false
CONFIGURE_ARGS=""
COMPILE_ARGS=""
GUI_ARGS="gui"

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
			DO_FFTW=false
			DO_TIFF=false
			DO_ARPACK=false;;        
        "-j")             TAKE_CPU=true;;
        "untar=true")   DO_UNTAR=true;;
        "untar=false")   DO_UNTAR=false;;
        "sqlite=true")   DO_SQLITE=true;;
        "sqlite=false")   DO_SQLITE=false;;
        "tcltk=true")   DO_TCLTK=true;;
        "tcltk=false")   DO_TCLTK=false;;
        "python=true")   DO_PYTHON=true;;
        "python=false")   DO_PYTHON=false;;
        "fftw=true")   DO_FFTW=true;;
        "fftw=false")   DO_FFTW=false;;
        "tiff=true")   DO_TIFF=true;;
        "tiff=false")   DO_TIFF=false;;
        "arpack=true")   DO_ARPACK=true;;
        "arpack=false")   DO_ARPACK=false;;
        "clean=true")   DO_CLEAN=true;;
        "clean=false")   DO_CLEAN=false;;
        "static=true")   DO_STATIC=true;;
        "static=false")   DO_STATIC=false;;
        "gui=false")   GUI_ARGS="";;
        # This two if passed should be at the end and 
        # will setup arguments for configure and compilation steps
        "configure") TAKE_CONFIGURE=true;
                     TAKE_COMPILE=false;;
        "compile")   TAKE_CONFIGURE=false;
                     TAKE_COMPILE=true;;
         *)          echo "Unrecognized option $param, exiting..."; exit 1
    esac
 fi 
done

#Some path variables
export XMIPP_HOME=$PWD
export PATH=$XMIPP_HOME/bin:$PATH
export LD_LIBRARY_PATH=$XMIPP_HOME/lib:$LD_LIBRARY_PATH

#create file to include from BASH this Xmipp installation
INC_FILE=.xmipp.bashrc
echo "export XMIPP_HOME=$PWD" > $INC_FILE
echo 'export PATH=$XMIPP_HOME/bin:$PATH' >> $INC_FILE
echo 'export LD_LIBRARY_PATH=$XMIPP_HOME/lib:$LD_LIBRARY_PATH' >> $INC_FILE
chmod u+x $INC_FILE

# for CSH or TCSH
INC_FILE=.xmipp.csh
echo "setenv XMIPP_HOME $PWD" > $INC_FILE
echo 'setenv PATH $XMIPP_HOME/bin:$PATH' >> $INC_FILE
echo 'setenv LD_LIBRARY_PATH $XMIPP_HOME/lib:$LD_LIBRARY_PATH' >> $INC_FILE
chmod u+x $INC_FILE

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


GREEN="\033[32m"
RED="\033[31m"
ENDC="\033[0m"

# Print a green msg using terminal escaped color sequence
echoGreen()
{
    printf "$GREEN %b $ENDC\n" "$1"
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
fi

#################### PYTHON ###########################
if $DO_PYTHON; then
    echoGreen "PYTHON SETUP"
    EXT_PYTHON=$EXT_PATH/python
    export CPPFLAGS="-I$EXT_PATH/$VSQLITE/ -I$EXT_PYTHON/tk$VTCLTK/generic -I$EXT_PYTHON/tcl$VTCLTK/generic"
    export LDFLAGS="-L$XMIPP_HOME/lib -L$EXT_PYTHON/tk$VTCLTK/unix -L$EXT_PYTHON/tcl$VTCLTK/unix"
    export LD_LIBRARY_PATH="$EXT_PYTHON/tk$VTCLTK/unix:$EXT_PYTHON/tcl$VTCLTK/unix:$LD_LIBRARY_PATH"
    echo "--> export CPPFLAGS=$CPPFLAGS"
    echo "--> export LDFLAGS=$LDFLAGS"
    # Copy our custom python files:
    cd $EXT_PYTHON
    echo "-->  cd $EXT_PYTHON"
    cp ./xmipp_setup.py $VPYTHON/setup.py
    echo "--> cp ./xmipp_setup.py $VPYTHON/setup.py"
    compile_library $VPYTHON python "." ""

    # Create the python launch script with necessary environment variable settings
    PYTHON_BIN=$XMIPP_HOME/bin/xmipp_python
    printf "#!/bin/sh\n\n" > $PYTHON_BIN
    printf 'VPYTHON=%b \n' "$VPYTHON" >> $PYTHON_BIN
    printf 'VTCLTK=%b \n\n' "$VTCLTK" >> $PYTHON_BIN
    printf 'EXT_PYTHON=$XMIPP_HOME/external/python \n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$EXT_PYTHON/$VPYTHON:$EXT_PYTHON/tcl$VTCLTK/unix:$EXT_PYTHON/tk$VTCLTK/unix:$LD_LIBRARY_PATH \n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$XMIPP_HOME/applications/tests/pythonlib:$XMIPP_HOME/lib/python2.7/site-packages:$PYTHONPATH \n' >> $PYTHON_BIN
    printf 'export TCL_LIBRARY=$EXT_PYTHON/tcl$VTCLTK/library \n' >> $PYTHON_BIN
    printf 'export TK_LIBRARY=$EXT_PYTHON/tk$VTCLTK/library \n\n' >> $PYTHON_BIN
    printf '$EXT_PYTHON/$VPYTHON/python "$@" \n ' >> $PYTHON_BIN
    chmod u+x $PYTHON_BIN
    
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

# Launch the configure/compile python script 
cd $XMIPP_HOME

#echoGreen "Compiling XMIPP ..."
#echoGreen "CONFIGURE: $CONFIGURE_ARGS"
#echoGreen "COMPILE: $COMPILE_ARGS"
#echoGreen "GUI: $GUI_ARGS"

echo "--> ./xmipp -j $NUMBER_OF_CPU configure $CONFIGURE_ARGS compile $COMPILE_ARGS $GUI_ARGS" install
./xmipp -j $NUMBER_OF_CPU configure $CONFIGURE_ARGS compile $COMPILE_ARGS $GUI_ARGS install

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
