#!/bin/sh
#set -x #uncomment for debugging

#############################
# XMIPP installation script #
# --------------------------################################################################################
# X-Window-based Microscopy Image Processing Package                                                       #
# Xmipp is a suite of image processing programs, primarily aimed at single-particle 3D electron microscopy #
# More info: http://xmipp.cnb.csic.es - xmipp@cnb.csic.es                                                  #
# Instruct Image Processing Center - National Center of Biotechnology - CSIC. Spain                        #
############################################################################################################


##################################################################################
#################### DEFINITIONS #################################################
##################################################################################

#Some interface definitions
TITLE="Xmipp installation script"
BLACK='\033[30m'
WHITE='\033[37m'
YELLOW='\033[33m'
RED='\033[31m'
BLUE='\033[36m'
GREEN="\033[32m"
RED="\033[31m"
ENDC="\033[0m"

export NUMBER_OF_CPU=1
GLOB_STATE=0
GLOB_COMMAND=""
TIMESTAMP=""
EXT_PATH=
BUILD_PATH=
OS_TYPE=$(uname)
IS_MAC=0
IS_MINGW=0
IS_LINUX=0
DELETE_ANSWER="n"
INTERARRAY=
INTERELEMARRAY=

#Some flags variables
DO_UNTAR=1
DO_COMPILE=1
DO_CONFIGURE=1
DO_CLEAN=1
DO_SETUP=1
DO_GUI=1
DO_UNATTENDED=0

DO_CLTOMO=0 #for scipy
CONFIGURE_ARGS=""
COMPILE_ARGS=""
GUI_ARGS="gui"


#External libraries definitions
VALGLIB=3.8.0
ALGLIB_FOLDER="alglib-${VALGLIB}.cpp"
ALGLIB_TAR="${ALGLIB_FOLDER}.tgz"
DO_ALGLIB=0

VBILIB=0.0
BILIB_FOLDER=bilib
BILIB_TAR="${BILIB_FOLDER}.tgz"
DO_BILIB=0

VCONDOR=0.0
CONDOR_FOLDER=condor
CONDOR_TAR="${CONDOR_FOLDER}.tgz"
DO_CONDOR=0

VFFTW=3.3.3
FFTW_FOLDER="fftw-${VFFTW}"
FFTW_TAR="${FFTW_FOLDER}.tgz"
DO_FFTW=0

VGTEST=1.6.0
GTEST_FOLDER="gtest-${VGTEST}"
GTEST_TAR="${GTEST_FOLDER}.tgz"
DO_GTEST=0

VHDF5=1.8.10
HDF5_FOLDER="hdf5-${VHDF5}"
HDF5_TAR="${HDF5_FOLDER}.tgz"
DO_HDF5=0

VIMAGEJ=1.45g
IMAGEJ_FOLDER="imagej"
IMAGEJ_TAR="${IMAGEJ_FOLDER}.tgz"
DO_IMAGEJ=0

VJPEG=8c
JPEG_FOLDER="jpeg-${VJPEG}"
JPEG_TAR="jpegsrc.v${VJPEG}.tgz"
DO_JPEG=0

VSCONS=1.2.0
SCONS_FOLDER="scons"
SCONS_TAR="${SCONS_FOLDER}.tgz"
DO_SCONS=0

VSQLITE=3.6.23
SQLITE_FOLDER="sqlite-${VSQLITE}"
SQLITE_TAR="${SQLITE_FOLDER}.tgz"
SQLITE_EXT_FOLDER=sqliteExt
DO_SQLITE=0

VTIFF=3.9.4
TIFF_FOLDER="tiff-${VTIFF}"
TIFF_TAR="${TIFF_FOLDER}.tgz"
DO_TIFF=0

VNMA=0.0
NMA_FOLDER="NMA"
NMA_TAR="${NMA_FOLDER}.tgz"
DO_NMA=0

VSHALIGNMENT=0.0
SHALIGNMENT_FOLDER="sh_alignment"
SHALIGNMENT_TAR="${SHALIGNMENT_FOLDER}.tgz"
DO_SHALIGNMENT=0

#External libraries arrays. For adding a new external library, define his decompress folder and tar names, put them in EXTERNAL_LIBRARIES and EXTERNAL_LIBRARIES_FILES arrays, and put 1 whther it has to be installed by default, 0 otherwise in EXTERNAL_LIBRARIES_DEFAULT in the appropiate positions. You will also need to stablish a DO* variable in the *DO array, to let the script work
EXTERNAL_LIBRARIES="         $ALGLIB_FOLDER $BILIB_FOLDER $CONDOR_FOLDER $FFTW_FOLDER $GTEST_FOLDER $HDF5_FOLDER $IMAGEJ_FOLDER $JPEG_FOLDER $SCONS_FOLDER $SQLITE_FOLDER $TIFF_FOLDER $NMA_FOLDER $SHALIGNMENT_FOLDER "
EXTERNAL_LIBRARIES_FILES="   $ALGLIB_TAR    $BILIB_TAR    $CONDOR_TAR    $FFTW_TAR    $GTEST_TAR    $HDF5_TAR    $IMAGEJ_TAR    $JPEG_TAR    $SCONS_TAR    $SQLITE_TAR    $TIFF_TAR    $NMA_TAR    $SHALIGNMENT_TAR    "
EXTERNAL_LIBRARIES_DO="      $DO_ALGLIB     $DO_BILIB     $DO_CONDOR     $DO_FFTW     $DO_GTEST     $DO_HDF5     $DO_IMAGEJ     $DO_JPEG     $DO_SCONS     $DO_SQLITE     $DO_TIFF     $DO_NMA     $DO_SHALIGNMENT     "
EXTERNAL_LIBRARIES_DEFAULT="        1             1              1            1             1            1              1            1             1              1            1           0                   0       "


#Python definitions
#Python
VPYTHON=2.7.2
PYTHON_FOLDER="Python-${VPYTHON}"
PYTHON_TAR="${PYTHON_FOLDER}.tgz"
DO_PYTHON=1

#Python modules
DO_PYMOD=0

VMATLIBPLOT=1.1.0
MATLIBPLOT_FOLDER="matplotlib-${VMATLIBPLOT}"
MATLIBPLOT_TAR="${MATLIBPLOT_FOLDER}.tgz"
DO_MATLIBPLOT=0

VPYMPI=1.2.2
PYMPI_FOLDER="mpi4py-${VPYMPI}"
PYMPI_TAR="${PYMPI_FOLDER}.tgz"
DO_PYMPI=0

VNUMPY=1.6.1
NUMPY_FOLDER="numpy-${VNUMPY}"
NUMPY_TAR="${NUMPY_FOLDER}.tgz"
DO_NUMPY=0

VSCIPY=0.12.0
SCIPY_FOLDER="scipy-${VSCIPY}"
SCIPY_TAR="${SCIPY_FOLDER}.tgz"
DO_SCIPY=0

VPSUTIL=0.7.1
PSUTIL_FOLDER="psutil-${VPSUTIL}"
PSUTIL_TAR="${PSUTIL_FOLDER}.tgz"
DO_PSUTIL=1

VTCLTK=8.5.10
DO_TCLTK=1
TCL_FOLDER="tcl${VTCLTK}"
TCL_TAR="${TCL_FOLDER}.tgz"
DO_TCL=0
TK_FOLDER="tk${VTCLTK}"
TK_TAR="${TK_FOLDER}.tgz"
DO_TK=0

#Python modules arrays. For adding a new python module, define his decompress folder and tar names, put them in PYTHON_MODULES and PYTHON_MODULES_FILES arrays, and put 1 whther it has to be installed by default, 0 otherwise in PYTHON_MODULES_DEFAULT in the appropiate position. You will also need to stablish a DO* variable in the *DO array, to let the script work
PYTHON_MODULES="        $MATLIBPLOT_FOLDER $PYMPI_FOLDER $NUMPY_FOLDER $SCIPY_FOLDER $TCL_FOLDER $TK_FOLDER $PSUTIL_FOLDER "
PYTHON_MODULES_FILES="  $MATLIBPLOT_TAR    $PYMPI_TAR    $NUMPY_TAR    $SCIPY_TAR    $TCL_TAR    $TK_TAR    $PSUTIL_TAR    "
PYTHON_MODULES_DO="     $DO_MATLIBPLOT     $DO_PYMPI     $DO_NUMPY     $DO_SCIPY     $DO_TCL     $DO_TK     $DO_PSUTIL     "
PYTHON_MODULES_DEFAULT="           1             1             1             0           1          1          1           "


##################################################################################
#################### FUNCTIONS ###################################################
##################################################################################

#function that searchs in an array for an element and returns the position of that element in the array
indexOf(){
  local i=1 S=$1 found=0; shift
    for elem in $@; do
      if [ -n $elem ]; then
      if [ "$S" != "$elem" ]; then
        i=$(expr $i + 1)
      else
        found=$i
      fi
      fi
    done
  return $found
}

setElem()
{
  local i=1 elemS=$1 value=$2; shift 2
  array=$@
  found=1
  INTERARRAY=''
  for elem in $@; do
    if [ $i -eq $elemS ]; then
      INTERARRAY="${INTERARRAY} ${value}"
      found=0
    else
      INTERARRAY="${INTERARRAY} ${elem}"
    fi
    i=$(expr $i + 1)
  done
  return $found
}

elemAt()
{
  local i=0 S=$1 found=0; shift
  INTERELEMARRAY=''
  for elem in $@; do
    if [ $i -eq $S ]; then
      INTERELEMARRAY=${elem}
      found=$i
    fi
    i=$(expr $i + 1)
  done
  return $found
}

len()
{
  local i=0
  for elem in $@; do
    i=$(expr $i + 1)
  done
  return $i
}

decideOS()
{
  echo "The OS is $OS_TYPE"
  case "$OS_TYPE" in
    Darwin)
      IS_MAC=1
      CONFIGURE_ARGS="mpi=True MPI_CXX=mpic++ MPI_LINKERFORPROGRAMS=mpic++"
      ;;
    MINGW*)
      IS_MINGW=1
      ;;
    *)
      IS_LINUX=1
      ;;
  esac
}

# Helping function to get the timestamp for measuring the time spent at any part of the code
tic()
{
   TIMESTAMP="$(date +%s)"
}

# Helping function to get the second timestamp for measuring the difference with the first and then know the time spent at any part of the code
toc()
{
   NOW="$(date +%s)"
   ELAPSED="$(expr $NOW - $TIMESTAMP)"
   echo "*** Elapsed time: $ELAPSED seconds"
}

# Print a green msg using terminal escaped color sequence
echoGreen()
{
    printf "$GREEN ${1} $ENDC\n"
}

# Pinrt a red msg using terminal escaped color sequence
echoRed()
{
    printf "$RED ${1} $ENDC\n"
}

# Check last return status by checking GLOB_STATE and GLOB_COMMAND vars. It can receive a parameter, If 1 is given means the program has to exit if non-zero return state is detected
check_state()
{
  if [ $GLOB_STATE -ne 0 ]; then
    echoRed "WARNING: command returned a non-zero status (${GLOB_STATE})"
    echoRed "COMMAND: $GLOB_COMMAND"
    case $1 in
      1)
        exitGracefully $GLOB_STATE "$GLOB_COMMAND"
      ;;
    esac
  fi
  return $GLOB_STATE
}

# Function that echoes the provided command passed as first argument and redirect it to file passed as second one checking its returning state. If 1 is given as third argument, then a non-zero exit status will result in a program exit.
echoExec()
{
  if [ $# -gt 3 ]; then
    echoRed "Error: bad parameter number on echoExec function. Exiting"
    exitGracefully
  fi
  COMMAND="$1"
  if [ $# -ne 1 ]; then
    REDIRECTION="$2"
  fi
  GLOB_COMMAND=${COMMAND}
  if ([ "${REDIRECTION}" = "/dev/null" ] || [ $# -eq 1 ]); then
    echo '-->' $COMMAND 
  else
    echo '-->' $COMMAND '>' $REDIRECTION '2>&1'
  fi
  if [ $# -eq 1 ]; then
    ${COMMAND}
  else
    ${COMMAND} > ${REDIRECTION} 2>&1
  fi
  GLOB_STATE=$?
  if [ $# -eq 2 ]; then
    check_state
  elif [ $# -eq 3 ]; then
    check_state 1
    if [ $3 -eq 2 ]; then
      if [ ${GLOB_STATE} -ne 0 ]; then
        echoRed "Printing Log:"
	cat ${REDIRECTION}
      fi
    fi
  fi
  return $GLOB_STATE
}

# Function that retreives from *DO arrays whether the asked external library or python module is enabled or not. It receives exactly 2 parameters:
# The first one can be "library" for external libraries or "pymodule" for python module
# The second one is the *TAR name of the library/module to search for
shouldIDoIt()
{
  ANS=0
  if [ $# -ne 2 ]; then
    echoRed "Error: bad parameter number on shouldIDoIt function. Exiting"
    exitGracefully
  fi

  case $1 in
    library)
      indexOf "$2" "${EXTERNAL_LIBRARIES_FILES}"
      ind=$?
      if [ ${ind} -lt 1 ]; then
        echoRed "Error: bad parameter ($2) on shouldIDoIt. File not found. Exiting"
        exitGracefully
      fi
      ind=$(expr $ind - 1)
      elemAt $ind "${EXTERNAL_LIBRARIES_DO}"
      ANS=${INTERELEMARRAY}
      ;;
    pymodule)
      indexOf "$2" "${PYTHON_MODULES_FILES}"
      ind=$?
      if [ ${ind} -lt 1 ]; then
        echoRed "Error: bad parameter ($2) on shouldIDoIt. File not found. Exiting"
        exitGracefully
      fi
      ind=$(expr $ind - 1)
      elemAt $ind "${PYTHON_MODULES_DO}"
      ANS=${INTERELEMARRAY}
      ;;
    *)
      echoRed "Error: bad parameter ($1) on shouldIDoIt function. Exiting"
      exitGracefully
      ;;
  esac

  return $ANS
}

# Function that set a value in the *DO arrays. It receives exactly 3 parameters:
# The first one can be "library" for external libraries or "pymodule" for python module
# The second one is the *TAR name of the library/module to search
# The third value is the value to set. It can be true or false
doIt()
{
  ANS=0
  if [ $# -ne 3 ]; then
    echoRed "Error: bad parameter number on doIt function. Exiting"
    exitGracefully
  fi

  case $1 in
    library)
      indexOf "$2" "${EXTERNAL_LIBRARIES_FILES}"
      ind=$?
      if [ ${ind} -lt 1 ]; then
        echoRed "Error: bad parameter ($2) on doIt. File not found. Exiting"
        exitGracefully
      fi
#      ind=$(expr $ind - 1)
      setElem $ind $3 "${EXTERNAL_LIBRARIES_DO}"
      EXTERNAL_LIBRARIES_DO=${INTERARRAY}
      ;;
    pymodule)
      indexOf "$2" "${PYTHON_MODULES_FILES}"
      ind=$?
      if [ ${ind} -lt 1 ]; then
        echoRed "Error: bad parameter ($2) on shouldIDoIt. File not found. Exiting"
        exitGracefully
      fi
 #     ind=$(expr $ind - 1)
      setElem $ind $3 "${PYTHON_MODULES_DO}"
      PYTHON_MODULES_DO=${INTERARRAY}
      ;;
    *)
      echoRed "Error: bad parameter ($1) on shouldIDoIt function. Exiting"
      exitGracefully
      ;;
  esac

  return $ANS
}

welcomeMessage()
{
  printf "${BLACK}0000000000000000000000000000000000000000000000000001\n"
  printf "${BLACK}0000000000000000P!!00000!!!!!00000!!0000000000000001\n"
  printf "${BLACK}000000000000P'  ${RED}.:==.           ,=;:.  ${BLACK}\"400000000001\n"
  printf "${BLACK}0000000000!  ${RED}.;=::.::,         .=:.-:=,.  ${BLACK}!000000001\n"
  printf "${BLACK}0000000P' ${RED}.=:-......::=      ::;.......-:=. ${BLACK}\"0000001\n"
  printf "${BLACK}0000000,${RED}.==-.........::;.   .;::.-.......-=:${BLACK}.a#00001\n"
  printf "${BLACK}0000000'${RED}.=;......--:.:.:=:.==::.:--:.:...,=- ${BLACK}!000001\n"
  printf "${BLACK}000000    ${RED}-=...:.:.:-:::::=;::.::.:.:..-:=-    ${BLACK}00001\n"
  printf "${BLACK}0000P       ${RED}==:.:-::. ${YELLOW}.aa.${RED} :::::::::.::=:      ${BLACK}\"4001\n"
  printf "${BLACK}00001        ${RED}:;:::::..${YELLOW}:#0:${RED} =::::::::::::        ${BLACK}j#01\n"
  printf "${BLACK}0001   ${YELLOW}aa _aas  _aa_  .aa, _a__aa_.  .a__aas,    ${BLACK}j#1\n"
  printf "${BLACK}"'0001   '"${YELLOW}"'4WU*!4#gW9!##i .##; 3##P!9#Ai .##U!!Q#_   '"${BLACK}j#1\n"
  printf "${BLACK}"'001    '"${YELLOW}"'3O.  :xO.  ]Oi .XX: ]O( .  X2 .XC;  :xX.   '"${BLACK}01\n"
  printf "${BLACK}"'001    '"${YELLOW}"'dU.  :jU.  ]Ui .WW: ]UL,..aXf .Ums  jd*    '"${BLACK}01\n"
  printf "${BLACK}"'0001   '"${YELLOW}"'4W.  :dW. .]W1 :WW: ]WVXNWO~  .#U*#WV!     '"${BLACK}01\n"
  printf "${BLACK}"'0001        '"${RED}"'.............. '"${YELLOW}"'[#1  '"${RED}"'.... '"${YELLOW}"'.#A)        '"${BLACK}j01\n"
  printf "${BLACK}"'00001      '"${RED}"':=::-:::::;;;;: '"${YELLOW}"'301 '"${RED}"'::::: '"${YELLOW}"'.0A)       '"${BLACK}j#01\n"
  printf "${BLACK}0000L    ${RED}.;::.--::::::::;: ,,..:::::... .::.   ${BLACK}_d001\n"
  printf "${BLACK}000000  ${RED}:;:...-.-.:.:::=;   =;:::---.:....::,  ${BLACK}00001\n"
  printf "${BLACK}000000!${RED}:::.....-.:.::::-     :=:-.:.-......:= ${BLACK}!00001\n"
  printf "${BLACK}000000a  ${RED}=;.........:;        .:;........,;;  ${BLACK}a00001\n"
  printf "${BLACK}"'0000000La '"${RED}"'--:_....:=-           :=:...:_:--'"${BLACK}_a#000001\n"
  printf "${BLACK}"'0000000000a  '"${RED}"'-=;,:=              -;::=:-  '"${BLACK}a000000001\n"
  printf "${BLACK}"'00000000000Laaa '"${RED}"'-\aa                - '"${BLACK}aaa00000000001\n"
  printf "${BLACK}00000000000000#00000000aaaaaad000000#00#0#0000000001"
  echo ""
  printf "${WHITE}Welcome to $TITLE ${ENDC}\n"
  echo ""
  return 0;
}

#function that prints the script usage help
helpMessage()
{
  printf "\n"
  printf "###################################\n"
  printf '# XMIPP INSTALLATION SCRIPT USAGE #\n'
  printf "###################################\n"
  printf "${RED}NAME\n"
  printf "${WHITE}  install.sh - XMIPP installation script. \n"
  printf "\n"
  printf "${RED}SYNOPSIS\n"
  printf "${WHITE}  ./install.sh ${BLUE}[OPERATIONS] [OPTIONS]${WHITE}\n"
  printf "\n"
  printf "${RED}DESCRIPTION\n"
  printf "${WHITE}  Script that automates the XMIPP compilation process. When this script is executed, the compilation sequence starts, depending on the selected options. If no option is given, the script follows the sequence:\n"
  printf "1- untar the external libraries, xmipp_python and xmipp_python modules.\n"
  printf "2- Configure and compile one by one every external library\n"
  printf "3- Configure and compile xmipp_python\n"
  printf "4- Compile and install python modules in xmipp_python\n"
  printf "5- Launch SConscript to compile Xmipp\n"
  printf "6- Run the unitary tests\n"
  printf "\n"
  printf "No option is mandatory. The default behaviour (without any given option) is untar, configure and compile a subset of the external libraries, python and launch the Xmipp compilation with GUI. Following options are accepted:\n"
  printf "\n"
  printf "\n"
  printf "GENERAL OPTIONS:\n"
  printf "\n"
  printf "${BLUE}--disable-all${WHITE},${BLUE} -d${WHITE}\n"
  printf "    Disable every option in order to let the user to enable just some of them.\n"
  printf "${BLUE}--num-cpus=${YELLOW}<NUMCPU>${WHITE},${BLUE} -j ${YELLOW}<NUMCPU>${WHITE}\n"
  printf "    Provides a number of CPUs for the compilation process. Default value is 2.\n"
  printf "${BLUE}--configure-args=${YELLOW}<ARGS>${WHITE}\n"
  printf "    Store the arguments the user wants to provide to configure process before compilation start.\n"
  printf "${BLUE}--compile-args=${YELLOW}<ARGS>${WHITE}\n"
  printf "    Store the arguments the user wants to provide to compilation process.\n"
  printf "${BLUE}--unattended=${YELLOW}[true|false]${WHITE}\n"
  printf "    Tells the script to asume default where no option is given on every question and don't show any GUI. This is intended to be used when no human will be watching the screen (i.e. other scripts).\n"
  printf "${BLUE}--help,${BLUE} -h${WHITE}\n"
  printf "    Shows this help message\n"
  printf "\n"
  printf "\n"
  printf "OPERATIONS:\n"
  printf "\n"
  printf "${BLUE}--clean${WHITE},${BLUE} -u${WHITE}\n"
  printf "    Clean the list of libraries provided to the script (or default libraries array if not provided).\n"
  printf "${BLUE}--untar${WHITE},${BLUE} -u${WHITE}\n"
  printf "    Untar the list of libraries provided to the script (or default libraries array if not provided).\n"
  printf "${BLUE}--configure${WHITE},${BLUE} -f${WHITE}\n"
  printf "    Configure the list of libraries provided to the script (or default libraries array if not provided).\n"
  printf "${BLUE}--compile${WHITE},${BLUE} -c${WHITE}\n"
  printf "    Compile the list of libraries provided to the script (or default libraries array if not provided).\n"
  printf "\n"
  printf "    NOTE: If you don't provide any operation, the script will asume you want all of them. This is used as a shortcut for recompiling from zero a library/module. For example, the command ./install.sh --disable-all --hdf5 will be the same as ./install.sh --disable-all --untar --clean --configure --compile --hdf5.\n"
  printf "\n"
  printf "\n"
  printf "EXTERNAL LIBRARIES OPTIONS:\n"
  printf "\n"
  printf "${BLUE}--alglib=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over alglib library. When just --alglib is given, true is asumed.\n"
  printf "${BLUE}--bilib=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over bilib library. When just --bilib is given, true is asumed.\n"
  printf "${BLUE}--condor=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over condor library. When just --condor is given, true is asumed.\n"
  printf "${BLUE}--fftw=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over fftw library. When just --fftw is given, true is asumed.\n"
  printf "${BLUE}--gtest=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over gtest library. When just --gtest is given, true is asumed.\n"
  printf "${BLUE}--hdf5=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over hdf5 library. When just --hdf5 is given, true is asumed.\n"
  printf "${BLUE}--imagej=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over imagej library. When just --hdf5 is given, true is asumed.\n"
  printf "${BLUE}--jpeg=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over jpeg library. When just --jpeg is given, true is asumed.\n"
  printf "${BLUE}--scons=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over scons library. When just --scons is given, true is asumed.\n"
  printf "${BLUE}--sqlite=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over sqlite library. When just --sqlite is given, true is asumed.\n"
  printf "${BLUE}--tiff=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over tiff library. When just --tiff is given, true is asumed.\n"
  printf "${BLUE}--nma=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over nma library. When just --nma is given, true is asumed. NMA uses fortran compiler. If FC environment variable is set, installer will use that fortran compiler, gfortran will be used otherwise. If FFLAGS variable is set those flags will be passed to the fortran compilation, otherwise -O3 will be asumed.\n"
  printf "${BLUE}--cltomo=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over cltomo library. When just --cltomo is given, true is asumed.\n"
  printf "\n"
  printf "\n"
  printf "PYTHON-RELATED OPTIONS:\n"
  printf "\n"
  printf "${BLUE}--python=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over xmipp_python. Just --python is equal to --python=true.\n"
  printf "${BLUE}--matplotlib=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over matplotlib module. --matplotlib is equivalent to --matplotlib=true.\n"
  printf "${BLUE}--mpi4py=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over mpi4py module. --mpi4py is equivalent to --mpi4py=true.\n"
  printf "${BLUE}--numpy=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over numpy module. --numpy is equivalent to --numpy=true.\n"
  printf "${BLUE}--psutil=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over psutil module. --scipy is equivalent to --psutil=true.\n"
  printf "${BLUE}--tcl-tk=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over tcl and tk libraries. --tcl-tk is equivalent to --tcl-tk=true.\n"
  printf "${BLUE}--tcl=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over tcl library. If --tcl is given, --tcl=true will be understood.\n"
  printf "${BLUE}--tk=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over tk library. Where --tk means --tk=true.\n"
  printf "${BLUE}--pymodules=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute or not selected operation over every python modules. When just --pymodules is given, true is asumed.\n"
  printf "\n"
  printf "\n"
  printf "XMIPP-RELATED OPTIONS:\n"
  printf "\n"
  printf "${BLUE}--xmipp=${YELLOW}[true|false]${WHITE}\n"
  printf "    Execute selected operation over xmipp. --xmipp is equivalent to --xmipp=true.\n"
  printf "${BLUE}--gui=${YELLOW}[true|false]${WHITE}\n"
  printf "    When launching scons compilation, select true whether you want the Xmipp compilation GUI or false otherwise. Where --gui means --gui=true.\n"
  printf "\n"
  printf "\n"
  return 0;
}

takeArguments()
{
  while [ "$1" ]; do
    case $1 in
      --disable-all|-d)
        DO_UNTAR=0
        DO_COMPILE=0
        DO_CONFIGURE=0
	DO_CLEAN=0
        for lib in ${EXTERNAL_LIBRARIES_FILES}; do
          doIt library ${lib} 0
        done
        for mod in ${PYTHON_MODULES_FILES}; do
          doIt pymodule ${mod} 0
        done
        DO_PYTHON=0
        DO_TCLTK=0
        DO_PYMOD=0
        DO_CLTOMO=0
        DO_SETUP=0
        ;;
      --num-cpus=*)
        ;;
      -j)
        NUMBER_OF_CPU=$2
        shift
        ;;
      -j*)
        NUMBER_OF_CPU=$(echo "$1"|sed -e 's/-j//g')
        ;;
      --configure-args=*)
        CONFIGURE_ARGS=$(echo "$1"|cut -d '=' -f2)
        ;;
      --compile-args=*)
        COMPILE_ARGS=$(echo "$1"|cut -d '=' -f2)
        ;;
      --unattended=*) #[true|false]
        WITH_UNATTENDED=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_UNATTENDED}" = "true" ]; then
          DO_UNATTENDED=1
        elif [ "${WITH_UNATTENDED}" = "false" ];then 
          DO_UNATTENDED=0
        else
          echoRed "Parameter --unattended only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --help|-h)
        helpMessage
        exitGracefully
        ;;
      --clean|-l)
        DO_CLEAN=1
        ;;
      --untar|-u)
        DO_UNTAR=1
        ;;
      --configure|-f)
        DO_CONFIGURE=1
        ;;
      --compile|-c)
        DO_COMPILE=1
        ;;
      --alglib)
        doIt library ${ALGLIB_TAR} 1
        ;;
      --alglib=*)
        WITH_ALGLIB=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_ALGLIB}" = "true" ]; then
          doIt library ${ALGLIB_TAR} 1
        elif [ "${WITH_ALGLIB}" = "false" ]; then
          doIt library ${ALGLIB_TAR} 0
        else
          echoRed "Parameter --alglib only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --bilib)
        doIt library ${BILIB_TAR} 1
        ;;
      --bilib=*)
        WITH_BILIB=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_BILIB}" = "true" ]; then
          doIt library ${BILIB_TAR} 1
        elif [ "${WITH_BILIB}" = "false" ]; then
          doIt library ${BILIB_TAR} 0
        else
          echoRed "Parameter --bilib only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --condor)
        doIt library ${CONDOR_TAR} 1
        ;;
      --condor=*)
        WITH_CONDOR=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_CONDOR}" = "true" ]; then
          doIt library ${CONDOR_TAR} 1
        elif [ "${WITH_CONDOR}" = "false" ]; then
          doIt library ${CONDOR_TAR} 0
        else
          echoRed "Parameter --condor only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --fftw)
        doIt library ${FFTW_TAR} 1
        ;;
      --fftw=*)
        WITH_FFTW=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_FFTW}" = "true" ]; then
          doIt library ${FFTW_TAR} 1
        elif [ "${WITH_FFTW}" = "false" ]; then
          doIt library ${FFTW_TAR} 0
        else
          echoRed "Parameter --fftw only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --gtest)
        doIt library ${GTEST_TAR} 1
        ;;
      --gtest=*)
        WITH_GTEST=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_GTEST}" = "true" ]; then
          doIt library ${GTEST_TAR} 1
        elif [ "${WITH_GTEST}" = "false" ]; then
          doIt library ${GTEST_TAR} 0
        else
          echoRed "Parameter --gtest only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --hdf5)
        doIt library ${HDF5_TAR} 1
        ;;
      --hdf5=*)
        WITH_HDF5=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_HDF5}" = "true" ]; then
          doIt library ${HDF5_TAR} 1
        elif [ "${WITH_HDF5}" = "false" ]; then
          doIt library ${HDF5_TAR} 0
        else
          echoRed "Parameter --hdf5 only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --imagej)
        doIt library ${IMAGEJ_TAR} 1
        ;;
      --imagej=*)
        WITH_IMAGEJ=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_IMAGEJ}" = "true" ]; then
          doIt library ${IMAGEJ_TAR} 1
        elif [ "${WITH_IMAGEJ}" = "false" ]; then
          doIt library ${IMAGEJ_TAR} 0
        else
          echoRed "Parameter --imagej only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --jpeg)
        doIt library ${JPEG_TAR} 1
        ;;
      --jpeg=*)
        WITH_JPEG=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_JPEG}" = "true" ]; then
          doIt library ${JPEG_TAR} 1
        elif [ "${WITH_JPEG}" = "false" ]; then
          doIt library ${JPEG_TAR} 0
        else
          echoRed "Parameter --jpeg only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --scons)
        doIt library ${SCONS_TAR} 1
        ;;
      --scons=*)
        WITH_SCONS=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_SCONS}" = "true" ]; then
          doIt library ${SCONS_TAR} 1
        elif [ "${WITH_SCONS}" = "false" ]; then
          doIt library ${SCONS_TAR} 0
        else
          echoRed "Parameter --scons only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --sqlite)
        doIt library ${SQLITE_TAR} 1
        ;;
      --sqlite=*)
        WITH_SQLITE=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_SQLITE}" = "true" ]; then
          doIt library ${SQLITE_TAR} 1
        elif [ "${WITH_SQLITE}" = "false" ]; then
          doIt library ${SQLITE_TAR} 0
        else
          echoRed "Parameter --sqlite only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --tiff)
        doIt library ${TIFF_TAR} 1
        ;;
      --tiff=*)
        WITH_TIFF=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_TIFF}" = "true" ]; then
          doIt library ${TIFF_TAR} 1
        elif [ "${WITH_TIFF}" = "false" ]; then
          doIt library ${TIFF_TAR} 0
        else
          echoRed "Parameter --tiff only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --nma)
        doIt library ${NMA_TAR} 1
        ;;
      --nma=*)
        WITH_NMA=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_NMA}" = "true" ]; then
          doIt library ${NMA_TAR} 1
        elif [ "${WITH_NMA}" = "false" ]; then
          doIt library ${NMA_TAR} 0
        else
          echoRed "Parameter --nma only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --cltomo)
        DO_CLTOMO=1
        doIt library ${SHALIGNMENT_TAR} 1
        doIt pymodule ${SCIPY_TAR} 1
        ;;
      --cltomo=*)
        WITH_CLTOMO=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_CLTOMO}" = "true" ]; then
          DO_CLTOMO=1
          doIt library ${SHALIGNMENT_TAR} 1
          doIt pymodule ${SCIPY_TAR} 1
        elif [ "${WITH_CLTOMO}" = "false" ]; then
          DO_CLTOMO=0
          doIt library ${SHALIGNMENT_TAR} 0
          doIt pymodule ${SCIPY_TAR} 0
        else
          echoRed "Parameter --cltomo only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --python)
        DO_PYTHON=1
        DO_PYMOD=1
        ;;
      --python=*)
        WITH_PYTHON=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_PYTHON}" = "true" ]; then
          DO_PYTHON=1
          DO_PYMOD=1
        elif [ "${WITH_PYTHON}" = "false" ]; then
          DO_PYTHON=1
          DO_PYMOD=1
        else
          echoRed "Parameter --python only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --matplotlib)
        doIt pymodule ${MATLIBPLOT_TAR} 1
        ;;
      --matplotlib=*)
        WITH_MATLIBPLOT=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_MATLIBPLOT}" = "true" ]; then
          doIt pymodule ${MATLIBPLOT_TAR} 1
        elif [ "${WITH_MATLIBPLOT}" = "false" ]; then
          doIt pymodule ${MATLIBPLOT_TAR} 0
        else
          echoRed "Parameter --matplotlib only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --mpi4py)
        doIt pymodule ${PYMPI_TAR} 1
        ;;
      --mpi4py=*)
        WITH_PYMPI=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_PYMPI}" = "true" ]; then
          doIt pymodule ${PYMPI_TAR} 1
        elif [ "${WITH_PYMPI}" = "false" ]; then
          doIt pymodule ${PYMPI_TAR} 0
        else
          echoRed "Parameter --mpi4py only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --numpy)
        doIt pymodule ${NUMPY_TAR} 1
        ;;
      --numpy=*)
        WITH_NUMPY=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_NUMPY}" = "true" ]; then
          doIt pymodule ${NUMPY_TAR} 1
        elif [ "${WITH_NUMPY}" = "false" ]; then
          doIt pymodule ${NUMPY_TAR} 0
        else
          echoRed "Parameter --numpy only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --psutil)
        doIt pymodule ${PSUTIL_TAR} 1
        ;;
      --psutil=*)
        WITH_PSUTIL=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_PSUTIL}" = "true" ]; then
          doIt pymodule ${PSUTIL_TAR} 1
        elif [ "${WITH_PSUTIL}" = "false" ]; then
          doIt pymodule ${PSUTIL_TAR} 0
        else
          echoRed "Parameter --psutil only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --tcl-tk)
        DO_TCLTK=1
        ;;
      --tcl-tk=*)
        WITH_TCLTK=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_TCLTK}" = "true" ]; then
          DO_TCLTK=1
        elif [ "${WITH_TCLTK}" = "false" ]; then
          DO_TCLTK=0
        else
          echoRed "Parameter --tcl-tk only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --tcl)
        doIt pymodule ${TCL_TAR} 1
        ;;
      --tcl=*)
        WITH_TCL=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_TCL}" = "true" ]; then
          doIt pymodule ${TCL_TAR} 1
        elif [ "${WITH_TCL}" = "false" ]; then
          doIt pymodule ${TCL_TAR} 0
        else
          echoRed "Parameter --tcl only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --tk)
        doIt pymodule ${TK_TAR} 1
        ;;
      --tk=*)
        WITH_TK=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_TK}" = "true" ]; then
          doIt pymodule ${TK_TAR} 1
        elif [ "${WITH_TK}" = "false" ]; then
          doIt pymodule ${TK_TAR} 0
        else
          echoRed "Parameter --tk only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --pymodules)
        DO_PYMOD=1
        ;;
      --pymodules=*)
        WITH_PYMOD=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_PYMOD}" = "true" ]; then
          DO_PYMOD=1
        elif [ "${WITH_PYMOD}" = "false" ]; then
          DO_PYMOD=0
        else
          echoRed "Parameter --pymodules only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      --xmipp)
        DO_SETUP=1
        ;;
      --xmipp=*)
        WITH_XMIPP=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_XMIPP}" = "true" ]; then
          DO_SETUP=1
        elif [ "${WITH_XMIPP}" = "false" ]; then
          DO_SETUP=0
        else
          echoRed "Parameter --xmipp only accept true or false values. Ignored and assuming default value."
        fi
        ;;

      --gui=*)
        WITH_GUI=$(echo "$1"|cut -d '=' -f2)
        if [ "${WITH_GUI}" = "true" ]; then
          WITH_GUI=1
        elif [ "${WITH_GUI}" = "false" ]; then
          WITH_GUI=0
        else
          echoRed "Parameter --gui only accept true or false values. Ignored and assuming default value."
        fi
        ;;
      *)
        echoRed "Error: Unrecognized option $1, exiting..."
        helpMessage
        exitGracefully
        ;;
    esac
    shift
  done

  # Some decisions to make after arguments are taken
  # to allow a user do everything on a library/module by just typing its name...
  if ([ ${DO_CLEAN} -eq 0 ] && [ ${DO_UNTAR} -eq 0 ] && [ ${DO_CONFIGURE} -eq 0 ] && [ ${DO_COMPILE} -eq 0 ] ); then
    DO_UNTAR=1
    DO_CLEAN=1
    DO_CONFIGURE=1
    DO_COMPILE=1
  fi
}

takeDefaults()
{
# Update *_DO arrays to represent what the *_DEFAULT arrays say we should do by default
  lib=0
  len "${EXTERNAL_LIBRARIES_FILES}"
  length=$?
  while [ $lib -lt $length ]; do
    elemAt $lib "${EXTERNAL_LIBRARIES_DEFAULT}"
    aux=${INTERELEMARRAY}
    aux2=$(expr $lib + 1)
    setElem $aux2 $aux "${EXTERNAL_LIBRARIES_DO}"
    EXTERNAL_LIBRARIES_DO=${INTERARRAY}
    lib=$(expr $lib + 1)
  done

  mod=0
  len "${PYTHON_MODULES_FILES}"
  length=$?
  while [ $mod -lt $length ]; do
    elemAt $mod "${PYTHON_MODULES_DEFAULT}"
    aux=${INTERELEMARRAY}
    aux2=$(expr $mod + 1)
    setElem $aux2 $aux "${PYTHON_MODULES_DO}"
    PYTHON_MODULES_DO=${INTERARRAY}
    mod=$(expr $mod + 1)
  done
}

create_bashrc_file()
{
  INC_FILE=$1
  echo "export XMIPP_HOME=$PWD" > $INC_FILE
  echo 'export PATH=$XMIPP_HOME/bin:$PATH' >> $INC_FILE
  echo 'export LD_LIBRARY_PATH=$XMIPP_HOME/lib:$LD_LIBRARY_PATH' >> $INC_FILE
  echo 'export XMIPP_FONT_NAME=Verdana' >> $INC_FILE
  echo 'export XMIPP_FONT_SIZE=10' >> $INC_FILE
  echo 'if [ "$BASH" != "" ]; then' >> $INC_FILE
  echo '# Load global autocomplete file ' >> $INC_FILE
  echo 'test -s $XMIPP_HOME/.xmipp.autocomplete && . $XMIPP_HOME/.xmipp.autocomplete || true' >> $INC_FILE
  echo '# Load programs autocomplete file ' >> $INC_FILE
  echo 'test -s $XMIPP_HOME/.xmipp_programs.autocomplete && . $XMIPP_HOME/.xmipp_programs.autocomplete || true' >> $INC_FILE
  echo 'fi' >> $INC_FILE
  
  if [ $IS_MAC -eq 1 ]; then
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
}

create_tcsh_file()
{
  INC_FILE=$1
  echo "setenv XMIPP_HOME $PWD" > $INC_FILE
  echo 'setenv PATH $XMIPP_HOME/bin:$PATH' >> $INC_FILE
  echo 'if($?LD_LIBRARY_PATH) then' >> $INC_FILE
  echo '  setenv LD_LIBRARY_PATH $XMIPP_HOME/lib:$LD_LIBRARY_PATH' >> $INC_FILE
  echo 'else' >> $INC_FILE
  echo '  setenv LD_LIBRARY_PATH $XMIPP_HOME/lib' >> $INC_FILE
  echo 'endif' >> $INC_FILE
  echo 'setenv XMIPP_FONT_NAME Verdana' >> $INC_FILE
  echo 'setenv XMIPP_FONT_SIZE 10' >> $INC_FILE
  
  if [ $IS_MAC -eq 1 ]; then
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
}
 
compile_library()
  {
  tic
  LIB=$1
  PREFIX_PATH=$2
  SUFFIX_PATH=$3
  CONFIGFLAGS=$4
  LIBS_PATH=$5
  _PATH=${EXT_PATH}/${PREFIX_PATH}/${LIB}/${SUFFIX_PATH}
  echo
  echoGreen "*** Compiling ${LIB} ..."
  echoExec "cd ${_PATH}" "/dev/null" 1
  echoExec "make -j $NUMBER_OF_CPU" "$BUILD_PATH/${LIB}_make.log" 2
  echoExec "cd -" "/dev/null" 1
  toc
  echo
}

configure_library()
{
  tic
  LIB=$1
  PREFIX_PATH=$2
  SUFFIX_PATH=$3
  CONFIGFLAGS=$4
  LIBS_PATH=$5
  _PATH=${EXT_PATH}/${PREFIX_PATH}/${LIB}/${SUFFIX_PATH}
  echo
  echoGreen "*** Configuring ${LIB} ..."
  echoExec "cd ${_PATH}" "/dev/null" 1

  echo "--> Enabling shared libraries..."
  CONFIGFLAGS="--enable-shared ${CONFIGFLAGS}"
  
  echoExec "./configure ${CONFIGFLAGS}" "${BUILD_PATH}/${LIB}_configure.log" 2
  echoExec "cd -" "/dev/null" 1
  toc
  echo
}

clean_library()
{
  tic
  LIB=$1
  PREFIX_PATH=$2
  SUFFIX_PATH=$3
  CONFIGFLAGS=$4
  LIBS_PATH=$5
  _PATH=${EXT_PATH}/${PREFIX_PATH}/${LIB}/${SUFFIX_PATH}
  echo
  echoGreen "*** Cleaning ${LIB} ..."
  echoExec "cd ${_PATH}" "/dev/null" 1
  echoExec "make distclean" "/dev/null"
  echoExec "cd -" "/dev/null" 1
  toc
  echo
}

compile_pymodule()
{
  tic
  MOD=$1
  _PATH=$EXT_PATH/python/$MOD
  echo
  echoGreen "*** Compiling ${MOD} ..."
  #_PYTHON=$EXT_PATH/python/$PYTHON_FOLDER/python
  echoExec "cd $_PATH" "/dev/null" 1
  echoExec "xmipp_python setup.py install --prefix $XMIPP_HOME" "$BUILD_PATH/${MOD}_setup_install.log" 2
  echoExec "cd -" "/dev/null" 1
  toc
  echo
}

#This function should be called from XMIPP_HOME
# Parameter: Library_Path Library_name Lib_Version_Number 
install_libs()
{
  tic
  echoExec "cd $XMIPP_HOME" "/dev/null" 1
  LIBPATH=external/$1; shift
  COMMON="$1"; shift
  VERSION=$1; shift
  COPY=$1
  SUFFIXES=".a .la "
  echo
  echoGreen "*** Installing lib ${LIBPATH} ..."
  if [ $IS_MAC -eq 1 ]; then
	SUFFIXES="$SUFFIXES .dylib .$VERSION.dylib"
  elif [ $IS_MINGW -eq 1 ]; then
        SUFFIXES="$SUFFIXES .dll.a -$VERSION.dll"
  else 
	SUFFIXES="$SUFFIXES .so .so.$VERSION"
  fi
  
  for suffix in $SUFFIXES; do
     LIBNAME=$COMMON$suffix
     if $COPY; then
	     echoExec "cp -f $LIBPATH/$LIBNAME lib/$LIBNAME" "/dev/null" 1
     else
	     echoExec "ln -sf ../$LIBPATH/$LIBNAME lib/$LIBNAME " "/dev/null" 1
     fi
  done
  echoExec "cd -" "/dev/null" 1
  toc
  echo
}

#This function should be called from XMIPP_HOME
# Parameter: Library_Path Library_name Lib_Version_Number 
uninstall_libs()
{
  tic
  echoExec "cd $XMIPP_HOME" "/dev/null" 1
  LIBPATH=external/$1; shift
  COMMON="$1"; shift
  VERSION=$1; shift
  COPY=$1
  SUFFIXES=".a .la "
  echo
  echoGreen "*** Uninstalling lib ${LIBPATH} ..."
  if [ $IS_MAC -eq 1 ]; then
	SUFFIXES="$SUFFIXES .dylib .$VERSION.dylib"
  elif [ $IS_MINGW -eq 1 ]; then
        SUFFIXES="$SUFFIXES .dll.a -$VERSION.dll"
  else 
	SUFFIXES="$SUFFIXES .so .so.$VERSION"
  fi
  
  for suffix in $SUFFIXES; do
     LIBNAME=$COMMON$suffix
     if [ -e lib/$LIBNAME ]; then
       echoExec "rm -f lib/$LIBNAME" "/dev/null"
     fi
  done
  echoExec "cd -" "/dev/null" 1
  toc
  echo
}

install_bin()
{
  tic
  echoExec "cd $XMIPP_HOME" "/dev/null" 1
  BINPATH=../external/$1
  LINKNAME=bin/$2
  echo
  echoGreen "*** Installing bin ${LINKNAME} ..."
  echoExec "ln -sf $BINPATH $LINKNAME" "/dev/null" 1
  echoExec "cd -" "/dev/null" 1
  toc
  echo
}

uninstall_bin()
{
  echoExec "cd $XMIPP_HOME" "/dev/null" 1
  BINPATH=../external/$1
  LINKNAME=bin/$2
  echo
  echoGreen "*** Uninstalling bin ${LINKNAME} ..."
  echoExec "rm $LINKNAME" "/dev/null"
  echoExec "cd -" "/dev/null" 1
  echo
}

create_dir()
{
  DIR=$1
  if [ -d $DIR ]; then 
    echoRed "--> Dir $DIR exists."
  else
    echoExec "mkdir $DIR" "/dev/null" 1
  fi
}

initial_definitions()
{
  export XMIPP_HOME=$PWD
  export PATH=$XMIPP_HOME/bin:$PATH
  export LD_LIBRARY_PATH=$XMIPP_HOME/lib:$LD_LIBRARY_PATH
  if [ $IS_MAC -eq 1 ]; then
    export DYLD_FALLBACK_LIBRARY_PATH=$XMIPP_HOME/lib:$DYLD_FALLBACK_LIBRARY_PATH
  fi
  EXT_PATH=$XMIPP_HOME/external
  BUILD_PATH=$XMIPP_HOME/build
}

decompressExternals()
{
  tic
  echo
  echoExec "cd ${EXT_PATH}" "/dev/null" 1
  echoGreen "*** Decompressing external libraries ..."
  lib=0
  len "${EXTERNAL_LIBRARIES}"
  length=$?
  while [ $lib -lt $length ]; do
    elemAt $lib "${EXTERNAL_LIBRARIES_DO}"
    if [ ${INTERELEMARRAY} -eq 1 ]; then
      elemAt $lib "${EXTERNAL_LIBRARIES}"
      if [ -d ${INTERELEMARRAY} ]; then
        if ([ $DO_UNATTENDED -eq 0 ] && [ ${DELETE_ANSWER} != "Y" ] && [ ${DELETE_ANSWER} != "N" ]); then
          echo "${INTERELEMARRAY} folder exists, do you want to permanently remove it? (y)es/(n)o/(Y)es-to-all/(N)o-to-all"
          read DELETE_ANSWER
        else
          DELETE_ANSWER="Y"
        fi
	if ([ ${DELETE_ANSWER} = "y" ] || [ ${DELETE_ANSWER} = "Y" ]); then
          echoExec "rm -rf ${INTERELEMARRAY}" "/dev/null"
        else
          echoRed "Library ${INTERELEMARRAY} folder remains untouched."
        fi
      fi
      elemAt $lib "${EXTERNAL_LIBRARIES_FILES}"
      echoExec "tar -xvzf ${INTERELEMARRAY}" "/dev/null" 1
    fi
    lib=$(expr $lib + 1)
  done
  echoExec "cd -" "/dev/null" 1
  toc
}

decompressPython()
{
  tic
  echo
  echoExec "cd ${EXT_PATH}/python" "/dev/null" 1
  echoGreen "*** Decompressing Python ***"
  if [ ${DO_PYTHON} -eq 1 ]; then
    if [ -d ${PYTHON_FOLDER} ]; then
      if ([ $DO_UNATTENDED -eq 0 ] && [ ${DELETE_ANSWER} != "Y" ] && [ ${DELETE_ANSWER} != "N" ]); then
        echo "${PYTHON_FOLDER} folder exists, do you want to permanently remove it? (y)es/(n)o/(Y)es-to-all/(N)o-to-all"
        read DELETE_ANSWER
      else
        DELETE_ANSWER="Y"
      fi
      if ([ ${DELETE_ANSWER} = "y" ] || [ ${DELETE_ANSWER} = "Y" ]); then
        echoExec "rm -rf ${PYTHON_FOLDER}" "/dev/null"
      else
        echoRed "${PYTHON_FOLDER} folder remains untouched."
      fi   
    fi
    echoExec "tar -xvzf ${PYTHON_TAR}" "/dev/null" 1
  fi
  echoExec "cd -" "/dev/null" 1
  toc
}

decompressPythonModules()
{
  tic
  echo
  echoExec "cd ${EXT_PATH}/python" "/dev/null" 1
  echoGreen "*** Decompressing Python modules ..."
  lib=0
  len "${PYTHON_MODULES}"
  length=$?
  while [ $lib -lt $length ]; do
    elemAt $lib "${PYTHON_MODULES_DO}"
    if [ ${INTERELEMARRAY} -eq 1 ]; then
      elemAt $lib "${PYTHON_MODULES}"
      if [ -d ${INTERELEMARRAY} ]; then
        if ([ $DO_UNATTENDED -eq 0 ] && [ ${DELETE_ANSWER} != "Y" ] && [ ${DELETE_ANSWER} != "N" ]); then
          elemAt $lib "${PYTHON_MODULES}"
          echo "${INTERELEMARRAY} folder exists, do you want to permanently remove it? (y)es/(n)o/(Y)es-to-all/(N)o-to-all"
          read DELETE_ANSWER
        else
          DELETE_ANSWER="Y"
        fi
	if ([ ${DELETE_ANSWER} = "y" ] || [ ${DELETE_ANSWER} = "Y" ]); then
          echoExec "rm -rf ${INTERELEMARRAY}" "/dev/null"
        else
          echoRed "Library ${INTERELEMARRAY} folder remains untouched."
        fi
      fi
      elemAt $lib "${PYTHON_MODULES_FILES}"
      echoExec "tar -xvzf ${INTERELEMARRAY}" "/dev/null" 1
    fi
    lib=$(expr $lib + 1)
  done
  echoExec "cd -" "/dev/null" 1
  toc
}

preparePythonEnvironment()
{
  export CPPFLAGS="-I${EXT_PATH}/${SQLITE_FOLDER} -I${EXT_PYTHON}/${TK_FOLDER}/generic -I${EXT_PYTHON}/${TCL_FOLDER}/generic"
  if [ $IS_MAC -eq 1 ]; then
    export LDFLAGS="-L${EXT_PYTHON}/${PYTHON_FOLDER} -L${XMIPP_HOME}/lib -L${EXT_PYTHON}/${TK_FOLDER}/macosx -L${EXT_PYTHON}/${TCL_FOLDER}/macosx"
    export LD_LIBRARY_PATH="${EXT_PYTHON}/${PYTHON_FOLDER}:${EXT_PYTHON}/${TK_FOLDER}/macosx:${EXT_PYTHON}/${TCL_FOLDER}/macosx:${LD_LIBRARY_PATH}"
    export DYLD_FALLBACK_LIBRARY_PATH="${EXT_PYTHON}/${PYTHON_FOLDER}:${EXT_PYTHON}/${TK_FOLDER}/macosx:${EXT_PYTHON}/${TCL_FOLDER}/macosx:${DYLD_FALLBACK_LIBRARY_PATH}"
    echoExec "ln -s ${XMIPP_HOME}/bin/xmipp_python ${XMIPP_HOME}/bin/python2.7" "/dev/null"
    echoExec "cd ${EXT_PYTHON}/${MATLIBPLOT_FOLDER}" "/dev/null" 1
    echoExec "ln -s ${XMIPP_HOME}/bin/xmipp_python ${XMIPP_HOME}/bin/pythonXmipp" "/dev/null" 
    echoExec "make -f make.osx clean" "/dev/null" 1
    echoExec "make -f make.osx PREFIX=${XMIPP_HOME} PYVERSION=Xmipp fetch deps mpl_install" "/dev/null" 1
    echoExec "rm ${XMIPP_HOME}/bin/pythonXmipp" "/dev/null"
    echoExec "rm ${XMIPP_HOME}/bin/python2.7" "/dev/null"
  elif [ $IS_MINGW -eq 1 ]; then
    export LDFLAGS="-L${EXT_PYTHON}/${PYTHON_FOLDER} -L${XMIPP_HOME}/lib -L${EXT_PYTHON}/${TK_FOLDER}/win -L${EXT_PYTHON}/${TCL_FOLDER}/win"
    export LD_LIBRARY_PATH="${EXT_PYTHON}/${PYTHON_FOLDER}:${EXT_PYTHON}/${TK_FOLDER}/win:${EXT_PYTHON}/${TCL_FOLDER}/win:${LD_LIBRARY_PATH}"
    echoExec "ln -s ${XMIPP_HOME}/bin/xmipp_python ${XMIPP_HOME}/bin/python2.7" "/dev/null"
    echoExec "cd ${EXT_PYTHON}/${MATLIBPLOT_FOLDER}" "/dev/null" 1
    echoExec "ln -s ${XMIPP_HOME}/bin/xmipp_python ${XMIPP_HOME}/bin/pythonXmipp" "/dev/null"
  else
    export LDFLAGS="-L${EXT_PYTHON}/${PYTHON_FOLDER} -L${XMIPP_HOME}/lib -L${EXT_PYTHON}/${TK_FOLDER}/unix -L${EXT_PYTHON}/${TCL_FOLDER}/unix"
    export LD_LIBRARY_PATH="${EXT_PYTHON}/${PYTHON_FOLDER}:${EXT_PYTHON}/${TK_FOLDER}/unix:${EXT_PYTHON}/${TCL_FOLDER}/unix:${LD_LIBRARY_PATH}"

    shouldIDoIt pymodule ${MATLIBPLOT_TAR}
    if [ $? -eq 1 ]; then
      echoExec "cp ${EXT_PYTHON}/matplotlib_setupext.py ${EXT_PYTHON}/${MATLIBPLOT_FOLDER}/setupext.py" "/dev/null" 1
    fi
    #The following is needed from matplotlib to works
    shouldIDoIt pymodule ${TK_TAR}
    if [ $? -eq 1 ]; then
      echoExec "cd ${EXT_PYTHON}/${TK_FOLDER}/unix/" "/dev/null" 1
      echoExec "ln -sf libtk8.5.so  libtk.so" "/dev/null"
      echoExec "cd -" "/dev/null" 1
    fi
    shouldIDoIt pymodule ${TCL_TAR}
    if [ $? -eq 1 ]; then
      echoExec "cd ${EXT_PYTHON}/${TCL_FOLDER}/unix/" "/dev/null" 1
      echoExec "ln -sf libtcl8.5.so  libtcl.so" "/dev/null"
      echoExec "cd -" "/dev/null" 1
    fi
  fi
}

exitGracefully()
{
  if [ $# -eq 1 ]; then
    GLOB_STATE=$1
    GLOB_COMMAND="Unknown"
  elif [ $# -eq 2 ]; then
    GLOB_STATE=$1
    GLOB_COMMAND="$2"
  elif [ $# -gt 2 ]; then
    echoRed "Error: too much parameters passed to exitGracefully function."
  fi

  if [ ${GLOB_STATE} -eq 0 ]; then
    echoGreen "Program exited succesfully."
  else
    echoRed "Program exited with non-zero status (${GLOB_STATE}), produced by command ${GLOB_COMMAND}."
  fi
  exit $GLOB_STATE
}

##################################################################################
#################### INITIAL TASKS ###############################################
##################################################################################

welcomeMessage
initial_definitions
takeDefaults
takeArguments $@
decideOS
create_bashrc_file .xmipp.bashrc # Create file to include from BASH this Xmipp installation
create_tcsh_file .xmipp.csh      # for CSH or TCSH

if [ $DO_UNATTENDED -eq 1 ]; then
  CONFIGURE_ARGS="$CONFIGURE_ARGS unattended "
fi


##################################################################################
#################### NEEDED FOLDERS: bin lib build ###############################
##################################################################################

echoGreen "*** Checking needed folders ..."
create_dir build
create_dir bin
create_dir lib


##################################################################################
#################### DECOMPRESSING EVERYTHING ####################################
##################################################################################

if [ $DO_UNTAR -eq 1 ]; then 
  decompressExternals
  decompressPython
  decompressPythonModules
fi

##################################################################################
#################### COMPILING EXTERNAL LIBRARIES ################################
##################################################################################


#################### ALGLIB ###########################
shouldIDoIt library ${ALGLIB_TAR}
if ([ $? -eq 1 ] && ([ $DO_COMPILE -eq 1 ] || [ $DO_CONFIGURE -eq 1 ])); then
  echo
  echoGreen "Note: alglib configure and compile is delegated into Scons. For compiling it xcompile must be used."
  echo
fi

#################### BILIB ###########################
shouldIDoIt library ${BILIB_TAR}
if ([ $? -eq 1 ] && ([ $DO_COMPILE -eq 1 ] || [ $DO_CONFIGURE -eq 1 ])); then
  echo
  echoGreen "Note: bilib configure and compile is delegated into Scons. For compiling it xcompile must be used."
  echo
fi

#################### CONDOR ###########################
shouldIDoIt library ${CONDOR_TAR}
if ([ $? -eq 1 ] && ([ $DO_COMPILE -eq 1 ] || [ $DO_CONFIGURE -eq 1 ])); then
  echo
  echoGreen "Note: condor configure and compile is delegated into Scons. For compiling it xcompile must be used."
  echo
fi

#################### FFTW ###########################
shouldIDoIt library ${FFTW_TAR}
if [ $? -eq 1 ]; then
  if [ $DO_CLEAN  -eq 1 ]; then
    uninstall_libs ${FFTW_FOLDER}/.libs libfftw3 3 true
    uninstall_libs ${FFTW_FOLDER}/threads/.libs libfftw3_threads 3 true
    uninstall_libs ${FFTW_FOLDER}/.libs libfftw3f 3 true
  fi
  if [ $DO_CONFIGURE  -eq 1 ]; then
    if [ $IS_MINGW -eq 1 ]; then
      FFTWFLAGS=" CPPFLAGS=-I/c/MinGW/include CFLAGS=-I/c/MinGW/include"
    else
      FFTWFLAGS=""
    fi
    FLAGS="${FFTWFLAGS} --enable-threads"
    configure_library ${FFTW_FOLDER} "." "." ${FLAGS}
  fi
  if [ $DO_COMPILE  -eq 1 ]; then
    compile_library ${FFTW_FOLDER} "." "." ${FLAGS}
    install_libs ${FFTW_FOLDER}/.libs libfftw3 3 true
    install_libs ${FFTW_FOLDER}/threads/.libs libfftw3_threads 3 true
  fi

  if [ $DO_CONFIGURE  -eq 1 ]; then
    FLAGS="${FFTWFLAGS} --enable-float"
    configure_library ${FFTW_FOLDER} "." "." ${FLAGS}
  fi
  if [ $DO_COMPILE  -eq 1 ]; then
    compile_library ${FFTW_FOLDER} "." "." ${FLAGS}
    install_libs ${FFTW_FOLDER}/.libs libfftw3f 3 true
  fi
fi

#################### GTEST ###########################
shouldIDoIt library ${GTEST_TAR}
if ([ $? -eq 1 ] && ([ $DO_COMPILE -eq 1 ] || [ $DO_CONFIGURE -eq 1 ])); then
  echo
  echoGreen "Note: gtest configure and compile is delegated into Scons. For compiling it xcompile must be used."
  echo
fi

#################### HDF5 ###########################
shouldIDoIt library ${HDF5_TAR}
if [ $? -eq 1 ]; then
  if [ $DO_CLEAN  -eq 1 ]; then
    uninstall_libs ${HDF5_FOLDER}/src/.libs libhdf5 7 false
    uninstall_libs ${HDF5_FOLDER}/c++/src/.libs libhdf5_cpp 7 false
  fi
  if [ $DO_COMPILE  -eq 1 ]; then
    configure_library ${HDF5_FOLDER} "." "." "CPPFLAGS=-w --enable-cxx"
    compile_library ${HDF5_FOLDER} "." "." "CPPFLAGS=-w --enable-cxx"
    install_libs ${HDF5_FOLDER}/src/.libs libhdf5 7 false
    install_libs ${HDF5_FOLDER}/c++/src/.libs libhdf5_cpp 7 false
  fi
fi

#################### IMAGEJ ###########################
shouldIDoIt library ${IMAGEJ_TAR}
if ([ $? -eq 1 ] && ([ $DO_COMPILE -eq 1 ] || [ $DO_CONFIGURE -eq 1 ])); then
  echo
  echoGreen "Note: imagej configure and compile is delegated into Scons. For compiling it xcompile must be used."
  echo
fi

#################### JPEG ###########################
shouldIDoIt library ${JPEG_TAR}
if [ $? -eq 1 ]; then
  if [ $DO_CLEAN  -eq 1 ]; then
    uninstall_libs ${JPEG_FOLDER}/.libs libjpeg 8 false
  fi
  if [ $DO_CONFIGURE  -eq 1 ]; then
    configure_library ${JPEG_FOLDER} "." "." "CPPFLAGS=-w"
  fi
  if [ $DO_COMPILE  -eq 1 ]; then
    compile_library ${JPEG_FOLDER} "." "." "CPPFLAGS=-w"
    install_libs ${JPEG_FOLDER}/.libs libjpeg 8 false
  fi
fi

#################### IMAGEJ ###########################
shouldIDoIt library ${SCONS_TAR}
if ([ $? -eq 1 ] && ([ $DO_COMPILE -eq 1 ] || [ $DO_CONFIGURE -eq 1 ])); then
  echo
  echoGreen "Note: scons doesn't need to be configured nor compiled. "
  echo
fi

#################### SQLITE ###########################
shouldIDoIt library ${SQLITE_TAR}
if [ $? -eq 1 ]; then
  if [ $DO_CLEAN  -eq 1 ]; then
    clean_library ${SQLITE_FOLDER} "." "." "CPPFLAGS=-w CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1" ".libs"
    if [ $IS_MINGW -eq 1 ]; then
      rm ${EXT_PATH}/${SQLITE_EXT_FOLDER}/libsqlitefunctions.dll ${XMIPP_HOME}/lib/libXmippSqliteExt.dll
    elif [ $IS_MAC -eq 1 ]; then
      rm ${EXT_PATH}/${SQLITE_EXT_FOLDER}/libsqlitefunctions.dylib ${XMIPP_HOME}/lib/libXmippSqliteExt.dylib
    else  
      rm ${EXT_PATH}/${SQLITE_EXT_FOLDER}/libsqlitefunctions.so ${XMIPP_HOME}/lib/libXmippSqliteExt.so
    fi
    uninstall_bin ${SQLITE_FOLDER}/sqlite3 xmipp_sqlite3
    uninstall_libs ${SQLITE_FOLDER}/.libs libsqlite3 0 false
  fi
  if [ $DO_CONFIGURE -eq 1 ]; then
    configure_library ${SQLITE_FOLDER} "." "." "CPPFLAGS=-w CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1" ".libs"
  fi
  if [ $DO_COMPILE -eq 1 ]; then
    compile_library ${SQLITE_FOLDER} "." "." "CPPFLAGS=-w CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1" ".libs"
    #execute sqlite to avoid relinking in the future
    echo "select 1+1 ;" | ${XMIPP_HOME}/external/${SQLITE_FOLDER}/sqlite3
    install_bin ${SQLITE_FOLDER}/sqlite3 xmipp_sqlite3
    install_libs ${SQLITE_FOLDER}/.libs libsqlite3 0 false
    #compile math library for sqlite
    cd ${EXT_PATH}/${SQLITE_EXT_FOLDER}
    if [ $IS_MINGW -eq 1 ]; then
      gcc -shared -I. -o libsqlitefunctions.dll extension-functions.c
      cp libsqlitefunctions.dll ${XMIPP_HOME}/lib/libXmippSqliteExt.dll
    elif [ $IS_MAC -eq 1 ]; then
      gcc -fno-common -dynamiclib extension-functions.c -o libsqlitefunctions.dylib
      cp libsqlitefunctions.dylib ${XMIPP_HOME}/lib/libXmippSqliteExt.dylib
    else  
      gcc -fPIC  -shared  extension-functions.c -o libsqlitefunctions.so -lm
      cp libsqlitefunctions.so ${XMIPP_HOME}/lib/libXmippSqliteExt.so
    fi
  fi
fi

#################### TIFF ###########################
shouldIDoIt library ${TIFF_TAR}
if [ $? -eq 1 ]; then
  if [ $DO_CLEAN  -eq 1 ]; then
    uninstall_libs ${TIFF_FOLDER}/libtiff/.libs libtiff 3 false
  fi
  if [ $DO_CONFIGURE  -eq 1 ]; then
    configure_library ${TIFF_FOLDER} "." "." "CPPFLAGS=-w --with-jpeg-include-dir=${EXT_PATH}/${JPEG_FOLDER} --with-jpeg-lib-dir=${XMIPP_HOME}/lib"
  fi
  if [ $DO_COMPILE  -eq 1 ]; then
    compile_library ${TIFF_FOLDER} "." "." "CPPFLAGS=-w --with-jpeg-include-dir=${EXT_PATH}/${JPEG_FOLDER} --with-jpeg-lib-dir=${XMIPP_HOME}/lib"
    install_libs ${TIFF_FOLDER}/libtiff/.libs libtiff 3 false
  fi
fi

#################### NMA ###########################
shouldIDoIt library ${NMA_TAR}
if [ $? -eq 1 ]; then
  echoExec "cd ${XMIPP_HOME}/external/NMA/ElNemo" "${XMIPP_HOME}/build/make_NMA.log" 1
  if [ $DO_CLEAN -eq 1 ]; then
    echoExec "make clean" "${XMIPP_HOME}/build/make_NMA.log"
  fi
  if [ $DO_COMPILE  -eq 1 ]; then
    FC=${FC:-gfortran}
    FFLAGS=${fc:-'-O3'}
    echoExec "make FC=${FC} FFLAGS=${FFLAGS}" "${XMIPP_HOME}/build/make_NMA.log" 1
    echoExec "cp nma_* ${XMIPP_HOME}/bin" "${XMIPP_HOME}/build/make_NMA.log"
    echoExec "cd -" "${XMIPP_HOME}/build/make_NMA.log"  1
    echoExec "cd ${XMIPP_HOME}/external/NMA/NMA_cart" "${XMIPP_HOME}/build/make_NMA.log" 1
    echoExec "make" "${XMIPP_HOME}/build/make_NMA.log" 1
    echoExec "cp nma_* ${XMIPP_HOME}/bin" "${XMIPP_HOME}/build/make_NMA.log"
    echoExec "cd -" "${XMIPP_HOME}/build/make_NMA.log"  1
    echoExec "cd ${XMIPP_HOME}" "${XMIPP_HOME}/build/make_NMA.log" 1
    echoExec "cp ${XMIPP_HOME}/external/NMA/nma_* ${XMIPP_HOME}/bin" "${XMIPP_HOME}/build/make_NMA.log"
    echoExec "cp ${XMIPP_HOME}/external/NMA/m_inout_Bfact.py ${XMIPP_HOME}/bin" "${XMIPP_HOME}/build/make_NMA.log"
  fi
  echoExec "cd -" "${XMIPP_HOME}/build/make_NMA.log"  1
fi

#################### TCL/TK ###########################
if [ $DO_TCLTK -eq 1 ]; then
  if [ $IS_MAC -eq 1 ]; then
    shouldIDoIt pymodule ${TCL_TAR}
    if [ $? -eq 1 ]; then
      if [ $DO_CONFIGURE  -eq 1 ]; then
        configure_library ${TCL_FOLDER} python macosx "--disable-xft"
      fi
      if [ $DO_COMPILE  -eq 1 ]; then
        compile_library ${TCL_FOLDER} python macosx "--disable-xft"
      fi
    fi
    shouldIDoIt pymodule ${TK_TAR}
    if [ $? -eq 1 ]; then
      if [ $DO_CONFIGURE  -eq 1 ]; then
        configure_library ${TK_FOLDER} python macosx "--disable-xft"
      fi
      if [ $DO_COMPILE  -eq 1 ]; then
        compile_library ${TK_FOLDER} python macosx "--disable-xft"
      fi
    fi
  elif [ $IS_MINGW -eq 1 ]; then
    shouldIDoIt pymodule ${TCL_TAR}
    if [ $? -eq 1 ]; then
      if [ $DO_CONFIGURE  -eq 1 ]; then
        configure_library ${TCL_FOLDER} python win "--disable-xft CFLAGS=-I/c/MinGW/include CPPFLAGS=-I/c/MinGW/include"
      fi
      if [ $DO_COMPILE  -eq 1 ]; then
        compile_library ${TCL_FOLDER} python win "--disable-xft CFLAGS=-I/c/MinGW/include CPPFLAGS=-I/c/MinGW/include"
      fi
    fi
    shouldIDoIt pymodule ${TK_TAR}
    if [ $? -eq 1 ]; then
      if [ $DO_CONFIGURE  -eq 1 ]; then
        configure_library ${TK_FOLDER} python win "--disable-xft --with-tcl=../../${TCL_FOLDER}/win CFLAGS=-I/c/MinGW/include CPPFLAGS=-I/c/MinGW/include"
      fi
      if [ $DO_COMPILE  -eq 1 ]; then
        compile_library ${TK_FOLDER} python win "--disable-xft --with-tcl=../../${TCL_FOLDER}/win CFLAGS=-I/c/MinGW/include CPPFLAGS=-I/c/MinGW/include"
      fi
    fi
  else
    shouldIDoIt pymodule ${TCL_TAR}
    if [ $? -eq 1 ]; then
      if [ $DO_CONFIGURE  -eq 1 ]; then
        configure_library ${TCL_FOLDER} python unix "--enable-threads"
      fi
      if [ $DO_COMPILE  -eq 1 ]; then
        compile_library ${TCL_FOLDER} python unix "--enable-threads"
      fi
    fi
    shouldIDoIt pymodule ${TK_TAR}
    if [ $? -eq 1 ]; then
      if [ $DO_CONFIGURE  -eq 1 ]; then
        configure_library ${TK_FOLDER} python unix "--enable-threads"
      fi
      if [ $DO_COMPILE  -eq 1 ]; then
        compile_library ${TK_FOLDER} python unix "--enable-threads"
      fi
    fi
  fi
fi

##################################################################################
#################### COMPILING PYTHON ############################################
##################################################################################

EXT_PYTHON=${EXT_PATH}/python

if [ $DO_PYTHON -eq 1 ]; then
  echoGreen "PYTHON SETUP"
  export CPPFLAGS="-I${EXT_PATH}/${SQLITE_FOLDER} -I${EXT_PYTHON}/${TK_FOLDER}/generic -I${EXT_PYTHON}/${TCL_FOLDER}/generic"
  if [ $IS_MAC -eq 1 ]; then
    export LDFLAGS="-L${EXT_PYTHON}/${PYTHON_FOLDER} -L${XMIPP_HOME}/lib -L${EXT_PYTHON}/${TK_FOLDER}/macosx -L${EXT_PYTHON}/${TCL_FOLDER}/macosx"
    export LD_LIBRARY_PATH="${EXT_PYTHON}/${PYTHON_FOLDER}:${EXT_PYTHON}/${TK_FOLDER}/macosx:${EXT_PYTHON}/${TCL_FOLDER}/macosx:${LD_LIBRARY_PATH}"
    export DYLD_FALLBACK_LIBRARY_PATH="${EXT_PYTHON}/${PYTHON_FOLDER}:${EXT_PYTHON}/${TK_FOLDER}/macosx:${EXT_PYTHON}/${TCL_FOLDER}/macosx:${DYLD_FALLBACK_LIBRARY_PATH}"
    echo "--> export CPPFLAGS=${CPPFLAGS}"
    echo "--> export LDFLAGS=${LDFLAGS}"
    echo "--> export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
    echo "--> export DYLD_FALLBACK_LIBRARY_PATH=${DYLD_FALLBACK_LIBRARY_PATH}"
  elif [ $IS_MINGW -eq 1 ]; then
    export CPPFLAGS="-I/usr/include -I/usr/include/ncurses -I/c/MinGW/include ${CPPFLAGS} -D__MINGW32__ -I${EXT_PYTHON}/${PYTHON_FOLDER}/Include -I${EXT_PYTHON}/${PYTHON_FOLDER}/PC "
    export CFLAGS="${CPPFLAGS}"
    export LDFLAGS="-L${EXT_PYTHON}/${PYTHON_FOLDER} -L${XMIPP_HOME}/lib -L${EXT_PYTHON}/${TK_FOLDER}/win -L${EXT_PYTHON}/${TCL_FOLDER}/win"
    export LD_LIBRARY_PATH="${EXT_PYTHON}/${PYTHON_FOLDER}:${EXT_PYTHON}/${TK_FOLDER}/win:${EXT_PYTHON}/${TCL_FOLDER}/win:${LD_LIBRARY_PATH}"
    echo "--> export CPPFLAGS=${CPPFLAGS}"
    echo "--> export LDFLAGS=${LDFLAGS}"
    echo "--> export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
  else
    export LDFLAGS="-L${EXT_PYTHON}/${PYTHON_FOLDER} -L${XMIPP_HOME}/lib -L${EXT_PYTHON}/${TK_FOLDER}/unix -L${EXT_PYTHON}/${TCL_FOLDER}/unix"
    export LD_LIBRARY_PATH="${EXT_PYTHON}/${PYTHON_FOLDER}:${EXT_PYTHON}/${TK_FOLDER}/unix:${EXT_PYTHON}/${TCL_FOLDER}/unix:${LD_LIBRARY_PATH}"
    echo "--> export CPPFLAGS=${CPPFLAGS}"
    echo "--> export LDFLAGS=${LDFLAGS}"
    echo "--> export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
  fi
  echoGreen "Copying our custom python files ..."
  echoExec "cd ${EXT_PYTHON}" "/dev/null" 1
  echoExec "cp ./xmipp_setup.py ${PYTHON_FOLDER}/setup.py" "/dev/null" 
  echoExec "chmod a+x ${PYTHON_FOLDER}/setup.py" "/dev/null" 1
  #cp ./xmipp_setup.py $PYTHON_FOLDER/setup.py
  #I thick these two are not needed
  #cp ./xmipp__iomodule.h $PYTHON_FOLDER/Modules/_io/_iomodule.h
  #echo "--> cp ./xmipp__iomodule.h $PYTHON_FOLDER/Modules/_io/_iomodule.h"
    
  configure_library ${PYTHON_FOLDER} python "." ""
  compile_library ${PYTHON_FOLDER} python "." ""

  # Create the python launch script with necessary environment variable settings
  PYTHON_BIN=${XMIPP_HOME}/bin/xmipp_python
  echo "--> Creating python launch script $PYTHON_BIN ..."
  printf "#!/bin/sh\n\n" > $PYTHON_BIN
  if [ $IS_MINGW -eq 1 ]; then
    printf 'PYTHON_FOLDER=Python27 \n' >> $PYTHON_BIN
    printf 'VTCLTK=8.5 \n\n' >> $PYTHON_BIN
    printf 'EXT_PYTHON=/c \n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$EXT_PYTHON/$PYTHON_FOLDER:$EXT_PYTHON/$PYTHON_FOLDER/tcl/tcl$VTCLTK:$EXT_PYTHON/$PYTHON_FOLDER/tcl/tk$VTCLTK:$LD_LIBRARY_PATH \n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$XMIPP_HOME/applications/tests/pythonlib:$EXT_PYTHON/$PYTHON_FOLDER:$XMIPP_HOME/lib:$PYTHONPATH \n' >> $PYTHON_BIN
    printf 'export TCL_LIBRARY=$EXT_PYTHON/$PYTHON_FOLDER/tcl/tcl$VTCLTK \n' >> $PYTHON_BIN
    printf 'export TK_LIBRARY=$EXT_PYTHON/$PYTHON_FOLDER/tcl/tk$VTCLTK \n\n' >> $PYTHON_BIN
    printf '$EXT_PYTHON/$PYTHON_FOLDER/python.exe "$@"\n' >> $PYTHON_BIN
  elif [ $IS_MAC -eq 1 ]; then
    printf "PYTHON_FOLDER=${PYTHON_FOLDER} \n" >> $PYTHON_BIN
    printf "VTCLTK=${VTCLTK} \n\n" >> $PYTHON_BIN
    printf 'EXT_PYTHON=$XMIPP_HOME/external/python \n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$EXT_PYTHON/$PYTHON_FOLDER:$EXT_PYTHON/tcl$VTCLTK/macosx:$EXT_PYTHON/tk$VTCLTK/macosx:$LD_LIBRARY_PATH \n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$XMIPP_HOME/applications/tests/pythonlib:$XMIPP_HOME/lib/python2.7/site-packages:$PYTHONPATH \n' >> $PYTHON_BIN
    printf 'export TCL_LIBRARY=$EXT_PYTHON/tcl$VTCLTK/library \n' >> $PYTHON_BIN
    printf 'export TK_LIBRARY=$EXT_PYTHON/tk$VTCLTK/library \n\n' >> $PYTHON_BIN
    printf 'export DYLD_FALLBACK_LIBRARY_PATH=$EXT_PYTHON/$PYTHON_FOLDER:$EXT_PYTHON/tcl$VTCLTK/macosx:$EXT_PYTHON/tk$VTCLTK/macosx:$DYLD_FALLBACK_LIBRARY_PATH \n' >> $PYTHON_BIN	
    printf '$EXT_PYTHON/$PYTHON_FOLDER/python.exe "$@"\n' >> $PYTHON_BIN
  else
    printf "PYTHON_FOLDER=${PYTHON_FOLDER} \n" >> $PYTHON_BIN
    printf "VTCLTK=${VTCLTK} \n\n" >> $PYTHON_BIN
    printf 'EXT_PYTHON=$XMIPP_HOME/external/python \n' >> $PYTHON_BIN
    printf 'export LD_LIBRARY_PATH=$EXT_PYTHON/$PYTHON_FOLDER:$EXT_PYTHON/tcl$VTCLTK/unix:$EXT_PYTHON/tk$VTCLTK/unix:$LD_LIBRARY_PATH \n' >> $PYTHON_BIN
    printf 'export PYTHONPATH=$XMIPP_HOME/lib:$XMIPP_HOME/protocols:$XMIPP_HOME/applications/tests/pythonlib:$XMIPP_HOME/lib/python2.7/site-packages:$PYTHONPATH \n' >> $PYTHON_BIN
    printf 'export TCL_LIBRARY=$EXT_PYTHON/tcl$VTCLTK/library \n' >> $PYTHON_BIN
    printf 'export TK_LIBRARY=$EXT_PYTHON/tk$VTCLTK/library \n\n' >> $PYTHON_BIN
#    printf 'source ${XMIPP_HOME}/bin/activate' >> $PYTHON_BIN
    printf '$EXT_PYTHON/$PYTHON_FOLDER/python "$@"\n' >> $PYTHON_BIN
  fi
  echoExec "chmod a+x ${PYTHON_BIN}" 1
  #make python directory accesible by anybody
  echoExec "chmod -R a+x ${XMIPP_HOME}/external/python/Python-2.7.2" 1

fi


##################################################################################
#################### COMPILING PYTHON MODULES ####################################
##################################################################################

preparePythonEnvironment

if [ $DO_PYMOD -eq 1 ]; then
  doIt pymodule ${NUMPY_TAR} 1
  doIt pymodule ${MATLIBPLOT_TAR} 1
  doIt pymodule ${PYMPI_TAR} 1
  doIt pymodule ${SCIPY_TAR} 1
  doIt pymodule ${PSUTIL_TAR} 1
fi

shouldIDoIt pymodule ${NUMPY_TAR}
if [ $? -eq 1 ]; then
  compile_pymodule ${NUMPY_FOLDER}
fi
shouldIDoIt pymodule ${MATLIBPLOT_TAR}
if [ $? -eq 1 ]; then
  compile_pymodule ${MATLIBPLOT_FOLDER}
fi
shouldIDoIt pymodule ${PYMPI_TAR}
if [ $? -eq 1 ]; then
  compile_pymodule ${PYMPI_FOLDER}
fi

##################################################################################
#################### COMPILING SHALIGNMENT AND SCIPY #############################
##################################################################################

if [ $DO_CLTOMO -eq 1 ]; then
  # Fast Rotational Matching
  export LDFLAGS="-shared ${LDFLAGS}"
  shouldIDoIt pymodule ${SCIPY_TAR}
  if [ $? -eq 1 ]; then
    compile_pymodule ${SCIPY_FOLDER}
  fi
  shouldIDoIt library ${SHALIGNMENT_TAR}
  if [ $? -eq 1 ]; then
    echoExec "cd ${EXT_PATH}/${SHALIGNMENT_FOLDER}" "${XMIPP_HOME}/build/cltomo.log" 1
    echoExec "./compile.sh" "${XMIPP_HOME}/build/cltomo.log" 1
  fi
fi

shouldIDoIt pymodule ${PSUTIL_TAR}
if [ $? -eq 1 ]; then
  compile_pymodule ${PSUTIL_FOLDER}
fi

# Launch the configure/compile python script 

#echoGreen "CONFIGURE: $CONFIGURE_ARGS"
#echoGreen "COMPILE: $COMPILE_ARGS"
#echoGreen "GUI: $GUI_ARGS"

if [ $DO_SETUP -eq 1 ]; then
  echoGreen "Compiling XMIPP ..."
  echoExec "cd ${XMIPP_HOME}" "/dev/null" 1
  echoExec "./setup.py -j ${NUMBER_OF_CPU} configure ${CONFIGURE_ARGS} compile ${COMPILE_ARGS} ${GUI_ARGS} install"
fi

exit ${GLOB_STATE}

