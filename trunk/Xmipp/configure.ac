#############################################################################
##									   ##
##  This file generates configure scripts. It is only for development      ##
##  processes. It should not be distributed with working versions          ##
##									   ##
#############################################################################

# Initialization of configure parameters and checks for common tools
# like C compiler, libtool tool, C flags, etc..

AC_INIT
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE(Xmipp,0.9.1)
AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_LANG(C++)


# If not --prefix is passed to the configure script, Xmipp programs and
# libraries will be installed under Xmipp dir by default.

AC_PREFIX_DEFAULT(`pwd`)

# xmipp_ string is prepended by default to installed xmipp programs. It can be 
# changed using --program-prefix=your_string

if test "x$program_prefix" = xNONE; then
	program_transform_name="s,^,xmipp_,"
fi

XMIPPDIR=`pwd`
AC_SUBST(XMIPPDIR)

# If --enable-debug , -g compiler flag
# will be enabled.

AC_ARG_ENABLE([debug],
  AC_HELP_STRING([--enable-debug], [Enable the use of debugging capabilities.]),
  [use_debug=$enableval], [use_debug=no])

if test $use_debug = no; then
	CFLAGS=`echo $CFLAGS | sed 's/\-g //g'`
	CXXFLAGS=`echo $CXXFLAGS | sed 's/\-g //g'`
fi

if test $use_debug = yes; then
	CFLAGS=$CFLAGS" -g"
	CXXFLAGS=$CXXFLAGS" -g"
fi

# If --enable-profiling , -p compiler flag
# will be enabled.

AC_ARG_ENABLE([profiling],
  AC_HELP_STRING([--enable-profiling], [Enable the use of profiling capabilities.]),
  [use_profiling=$enableval], [use_profiling=no])

if test $use_profiling = no; then
	CFLAGS=`echo $CFLAGS | sed 's/\-pg //pg'`
	CXXFLAGS=`echo $CXXFLAGS | sed 's/\-pg //pg'`
fi

if test $use_debug = yes; then
	CFLAGS=$CFLAGS" -pg"
	CXXFLAGS=$CXXFLAGS" -pg"
fi

AC_ARG_ENABLE([static-exec],
  AC_HELP_STRING([--enable-static-exec], [Enables static executables building]),
  [build_static_exec=$enableval], [build_static_exec=no])

if test $build_static_exec = yes; then
	EXEC_TYPE="-all-static"
else
	EXEC_TYPE=
fi

AC_SUBST(EXEC_TYPE)

# Checking for PERL, csh and sh interpreters
# This is needed for Xmipp scripts. If not found scripts wont be installed.

AC_PATH_PROG( PERL_PATH, perl, "not found" )
AC_PATH_PROG( CSH_PATH, csh, "not found" )
AC_PATH_PROG( SH_PATH, sh, "not found" )

#Checking for TIFF library as well as some of its capabilities
AC_CHECK_LIBTIFF

TIFF_LIBS=
if test x"$have_tiff" = xyes; then
	if test x"$have_tiffjpg" = xyes; then
		if test x"$have_tiffzip" = xyes ; then
			TIFF_LIBS="-lm -ljpeg -lz "
			use_tiff=yes
		fi
	fi
fi

AM_CONDITIONAL( BUILD_TIFF, test x"$use_tiff" = xyes )
AC_SUBST(TIFF_LIBS)

#Checking for the presence and usability of MPI library
ACX_MPI(USE_MPI=true)
AM_CONDITIONAL(BUILD_MPI, test x"$USE_MPI" = xtrue )

# checking for Qt
# If Qt is present then set HAVE_QT to 1 in config.h file
# else print error and exit 
BNV_HAVE_QT
AM_CONDITIONAL(BUILD_QT, test x"$have_qt" = xyes )

if have_qt=yes
then
	AC_DEFINE([HAVE_QT],1,Define HAVE_QT to 1 if Qt is present and usable)
else
	AC_MSG_RESULT([WARNING!!!: Qt Not found, Qt capabilities will be disabled.])
fi

#checking for OpenGL
MDL_HAVE_OPENGL

if test -z "$GL_LIBS"
then
	AC_MSG_RESULT([WARNING!!!: GL libs not found, OpenGL disabled.])
else
	if test -z "$GL_CFLAGS"
	then
		AC_MSG_RESULT([WARNING!!!: GL headers not found, OpenGL disabled.])
	else
		AC_MSG_RESULT([OpenGL flags: $GL_CFLAGS])
		AC_MSG_RESULT([OpenGL libs: $GL_LIBS])
	fi
fi

AC_ARG_ENABLE([GL],
  AC_HELP_STRING([--disable-GL], [Disable GL programs]),
  [case "${enableval}" in
     yes) USE_GL=true ;;
      no) USE_GL=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --disable-GL) ;;
   esac], [USE_GL=true])

AM_CONDITIONAL(BUILD_OGL, test -n "$GL_LIBS")
AM_CONDITIONAL(BUILD_OGL, test -n "$GL_CFLAGS")
AM_CONDITIONAL(BUILD_OGL, test x"$USE_GL" = xtrue)

# Checks for the presence of some system libs 

AC_CHECK_LIB(Xext, main)
AC_CHECK_LIB(Xt, main)
AC_CHECK_LIB(X11, main)
AC_CHECK_LIB(pthread, main)

# Generates Makefiles for the different source directories.

AC_OUTPUT( Other_packages/Makefile \
	Other_packages/LibINRIA/Makefile \
	Lib/Makefile \
	Lib/XmippData/Makefile \
	Lib/Reconstruction/Makefile \
	Lib/Classification/Makefile \
	Lib/XmippInterface/Makefile \
	Lib/OldXmipp/Makefile \
	Lib/XmippGraphics/Makefile \
	Applications/Makefile \
	Applications/OldXmipp/Makefile \
	Applications/OldXmipp/MakeSpectra/Makefile \
	Applications/OldXmipp/FindSCenter/Makefile \
	Applications/Batches/Makefile \
	Applications/Batches/Do_Selfile/Makefile \
	Applications/Batches/Gnuplot_mplot/Makefile \
	Applications/Batches/Gnuplot_mprint/Makefile \
	Applications/Batches/Iv/Makefile \
	Applications/Batches/Look_at_FFT/Makefile \
	Applications/Batches/PDB2descr/Makefile \
	Applications/Batches/Pdb2surface/Makefile \
	Applications/Batches/PDBsidechain/Makefile \
	Applications/Batches/Showsel/Makefile \
	Applications/Batches/ShowSOM/Makefile \
	Applications/Batches/ShowSpectra/Makefile \
	Applications/Batches/ShowSpectraSOM/Makefile \
	Applications/Src/Makefile \
	Applications/Src/2DInterpolation/Makefile \
	Applications/Src/3DInterpolation/Makefile \
	Applications/Src/Adapt_for_Spider/Makefile \
	Applications/Src/Add_noise/Makefile \
	Applications/Src/Adjust_CTF/Makefile \
	Applications/Src/Adjust_volume/Makefile \
	Applications/Src/Adjust_surface/Makefile \
	Applications/Src/Align2D/Makefile \
	Applications/Src/Align3D/Makefile \
	Applications/Src/Ang_distribution/Makefile \
	Applications/Src/Angular_distance/Makefile \
	Applications/Src/Angular_predict/Makefile \
	Applications/Src/Angular_predict_continuous/Makefile \
	Applications/Src/Angular_refinement/Makefile \
	Applications/Src/ApplyGeo/Makefile \
	Applications/Src/AppendDocFile/Makefile \
	Applications/Src/Are_different/Makefile \
	Applications/Src/Art/Makefile \
	Applications/Src/Assign_angles/Makefile \
	Applications/Src/Assign_CTF/Makefile \
	Applications/Src/Assign_multiple_CTFs/Makefile \
	Applications/Src/BatchSOM/Makefile \
	Applications/Src/Break_sym/Makefile \
	Applications/Src/Centilt/Makefile \
	Applications/Src/CodeStatis/Makefile \
	Applications/Src/Combine_stats/Makefile \
	Applications/Src/Compare_Selfiles/Makefile \
	Applications/Src/CorrectPhase/Makefile \
	Applications/Src/Correlation/Makefile \
	Applications/Src/CpSel/Makefile Applications/Src/Create_phantom/Makefile \
	Applications/Src/CTF_profile/Makefile \
	Applications/Src/deconvoluteCTF/Makefile \
	Applications/Src/Data2Img/Makefile \
	Applications/Src/Data2Vol/Makefile \
	Applications/Src/Deform/Makefile \
	Applications/Src/Denoising/Makefile \
	Applications/Src/Denoising_micrograph/Makefile \
	Applications/Src/DocFile_Histogram/Makefile \
	Applications/Src/DownSample/Makefile \
	Applications/Src/Draw_surface/Makefile \
	Applications/Src/DWT/Makefile \
	Applications/Src/Edit/Makefile \
        Applications/Src/EM2Angles/Makefile \
	Applications/Src/Enhance_PSD/Makefile \
	Applications/Src/Euler/Makefile \
	Applications/Src/Evaluate/Makefile \
	Applications/Src/Evaluate_FSCs/Makefile \
	Applications/Src/Exp_Shift_Computation/Makefile \
	Applications/Src/FCMeans/Makefile \
	Applications/Src/FFT/Makefile \
	Applications/Src/FindCenter3D/Makefile \
	Applications/Src/FKCN/Makefile \
	Applications/Src/Flip/Makefile \
	Applications/Src/FourierFilter/Makefile \
	Applications/Src/FSOM/Makefile \
	Applications/Src/Headerinfo/Makefile \
	Applications/Src/Histogram/Makefile \
	Applications/Src/IDR_ART/Makefile \
	Applications/Src/IDWT/Makefile \
	Applications/Src/IFFT/Makefile \
	Applications/Src/Img2Data/Makefile \
	Applications/Src/Infogeo/Makefile \
	Applications/Src/KCMeans/Makefile \
	Applications/Src/KerDenSOM/Makefile \
	Applications/Src/Lattice_vectors/Makefile \
	Applications/Src/MakeSpectra/Makefile \
	Applications/Src/Mark/Makefile \
	Applications/Src/Markhan/Makefile \
	Applications/Src/Mask/Makefile \
	Applications/Src/Microscope/Makefile \
	Applications/Src/MLalign2D/Makefile \
	Applications/Src/MLalign2D_combine/Makefile \
	Applications/Src/MLrefine3D/Makefile \
	Applications/Src/Morphology/Makefile \
	Applications/Src/MPIArt/Makefile \
	Applications/Src/MPI_WBP/Makefile \
	Applications/Src/MPI_MLalign2D/Makefile \
	Applications/Src/MPI_MLrefine3D/Makefile \
	Applications/Src/MPI_Projection_matching/Makefile \
	Applications/Src/MvSel/Makefile \
	Applications/Src/Normalize/Makefile \
	Applications/Src/Operate/Makefile \
	Applications/Src/PCA/Makefile \
	Applications/Src/PDBPhantom/Makefile \
	Applications/Src/Phan3D/Makefile \
	Applications/Src/Phantom_transform/Makefile \
	Applications/Src/Project/Makefile \
	Applications/Src/Projection_Neighbourhood/Makefile \
	Applications/Src/Projection_matching/Makefile \
	Applications/Src/ProjectPCA/Makefile \
	Applications/Src/Pyramid/Makefile \
	Applications/Src/Radon_Transform/Makefile \
	Applications/Src/RandomMapping/Makefile \
	Applications/Src/Random_phantom/Makefile \
	Applications/Src/Range_adjust/Makefile \
	Applications/Src/Raw22spi/Makefile \
	Applications/Src/Recons_test/Makefile \
	Applications/Src/Reverse_endian/Makefile \
	Applications/Src/Resolution/Makefile \
	Applications/Src/RmSel/Makefile \
	Applications/Src/Rotate/Makefile \
	Applications/Src/RotateDescr/Makefile \
	Applications/Src/Sammon/Makefile \
	Applications/Src/Sample_size/Makefile \
	Applications/Src/Scale/Makefile \
	Applications/Src/Scissor/Makefile \
	Applications/Src/Segment/Makefile \
	Applications/Src/Select_Images/Makefile \
	Applications/Src/SelStats/Makefile \
	Applications/Src/Separate_Objects/Makefile \
	Applications/Src/Shift/Makefile \
	Applications/Src/Show/Makefile \
	Applications/Src/Skew/Makefile \
	Applications/Src/SOM/Makefile \
	Applications/Src/SpAR/Makefile \
	Applications/Src/SpARMA/Makefile \
	Applications/Src/Spider22CCP4/Makefile \
	Applications/Src/Spots2RealSpace2D/Makefile \
	Applications/Src/Split_selfile/Makefile \
	Applications/Src/SSNR/Makefile \
	Applications/Src/Statis/Makefile \
	Applications/Src/Surface/Makefile \
	Applications/Src/Symmetrize/Makefile \
	Applications/Src/Taxa_tilt/Makefile \
	Applications/Src/Threshold/Makefile \
	Applications/Src/Tiff2Raw/Makefile \
	Applications/Src/ToPDB/Makefile \
	Applications/Src/ToPhantom/Makefile \
	Applications/Src/Unbend/Makefile \
	Applications/Src/Visualize_FFT/Makefile \
	Applications/Src/Vol2Data/Makefile \
	Applications/Src/VolSampler/Makefile \
	Applications/Src/Voxels22Blobs/Makefile \
	Applications/Src/Wbp/Makefile \
	Applications/Src/Window/Makefile \
	Applications/Src/Window_micrograph/Makefile \
	Applications/Src/XMask/Makefile \
	Makefile \
	m4/Makefile )
		
