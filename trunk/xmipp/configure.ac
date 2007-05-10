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
AM_INIT_AUTOMAKE(Xmipp,1.2)
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

# Fast mathematics

if test $GCC = yes ; then
	CFLAGS=$CFLAGS" -O3 -fomit-frame-pointer -ffast-math -finline-functions -funroll-loops"
	CXXFLAGS=$CXXFLAGS" -O3 -fomit-frame-pointer -ffast-math -finline-functions -funroll-loops"
else
	CFLAGS=$CFLAGS" -O"
	CXXFLAGS=$CXXFLAGS" -O"
fi

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
#QTLIBEXT=so
#AM_PATH_QT([3.3.0], $QTLIBEXT, have_qt=yes, have_qt=no)
AM_CONDITIONAL(BUILD_QT, test x"$have_qt" = xyes )

if have_qt=yes
then
	AC_DEFINE([HAVE_QT],1,Define HAVE_QT to 1 if Qt is present and usable)
else
	AC_MSG_RESULT([WARNING!!!: Qt Not found, Qt capabilities will be disabled.])
fi

# Checks for the presence of some system libs

AC_CHECK_LIB(Xext, main)
AC_CHECK_LIB(Xt, main)
AC_CHECK_LIB(X11, main)
AC_CHECK_LIB(pthread, main)

# Generates Makefiles for the different source directories.

AC_OUTPUT( \
    external/Makefile \
	external/bilib/Makefile \
	external/cuba/Makefile \
	external/inria/Makefile \
	libraries/Makefile \
	libraries/data/Makefile \
	libraries/reconstruction/Makefile \
	libraries/reconstruction_mpi/Makefile \
	libraries/classification/Makefile \
	libraries/interface/Makefile \
	libraries/legacy/Makefile \
	libraries/graphics/Makefile \
	applications/Makefile \
	applications/legacy/Makefile \
	applications/legacy/find_center2d/Makefile \
	applications/scripts/Makefile \
    applications/scripts/backwards/Makefile \
	applications/scripts/selfile_create/Makefile \
	applications/scripts/convert_pdb2descr/Makefile \
	applications/scripts/convert_pdb2surface/Makefile \
	applications/scripts/extract_sidechain_from_pdb/Makefile \
	applications/scripts/mlf_align2d/Makefile \
	applications/scripts/mpi_mlf_align2d/Makefile \
	applications/scripts/mlf_refine3d/Makefile \
	applications/scripts/mpi_mlf_refine3d/Makefile \
	applications/scripts/protocols/Makefile \
	applications/programs/Makefile \
	applications/programs/add_noise/Makefile \
	applications/programs/adjust_volume_grey_levels/Makefile \
	applications/programs/align2d/Makefile \
	applications/programs/align_tilt_pairs/Makefile \
	applications/programs/align_volumes/Makefile \
	applications/programs/angular_assign_for_tomogram/Makefile \
	applications/programs/angular_continuous_assign/Makefile \
	applications/programs/angular_discrete_assign/Makefile \
	applications/programs/angular_distance/Makefile \
	applications/programs/angular_distribution_show/Makefile \
	applications/programs/angular_neighbourhood/Makefile \
	applications/programs/angular_projection_matching/Makefile \
	applications/programs/average/Makefile \
	applications/programs/break_symmetry/Makefile \
	applications/programs/classify_batch_som/Makefile \
	applications/programs/classify_fcmeans/Makefile \
	applications/programs/classify_fkcn/Makefile \
	applications/programs/classify_fsom/Makefile \
	applications/programs/classify_kcmeans/Makefile \
	applications/programs/classify_kerdensom/Makefile \
	applications/programs/classify_pca/Makefile \
	applications/programs/classify_project_pca/Makefile \
	applications/programs/classify_sammon_projection/Makefile \
	applications/programs/classify_som/Makefile \
	applications/programs/convert_data2descr/Makefile \
	applications/programs/convert_data2img/Makefile \
	applications/programs/convert_data2pdb/Makefile \
	applications/programs/convert_data2vol/Makefile \
	applications/programs/convert_img2data/Makefile \
	applications/programs/convert_pdb2vol/Makefile \
	applications/programs/convert_raw22spi/Makefile \
	applications/programs/convert_spi22ccp4/Makefile \
	applications/programs/convert_spi22em/Makefile \
	applications/programs/convert_tiff2raw/Makefile \
	applications/programs/convert_vol2data/Makefile \
	applications/programs/convert_vol2pseudo/Makefile \
	applications/programs/convert_voxels22blobs/Makefile \
	applications/programs/correlation/Makefile \
	applications/programs/crystal_aph2img/Makefile \
	applications/programs/crystal_create_surface/Makefile \
	applications/programs/crystal_euler2mrc/Makefile \
	applications/programs/crystal_lattice_vectors/Makefile \
	applications/programs/crystal_shear/Makefile \
	applications/programs/crystal_skew/Makefile \
	applications/programs/crystal_unbend/Makefile \
	applications/programs/ctf_correct_idr/Makefile \
	applications/programs/ctf_correct_phase/Makefile \
	applications/programs/ctf_estimate_from_micrograph/Makefile \
	applications/programs/ctf_profile/Makefile \
	applications/programs/ctf_view/Makefile \
	applications/programs/denoise/Makefile \
	applications/programs/docfile_append/Makefile \
	applications/programs/docfile_histogram/Makefile \
	applications/programs/edit/Makefile \
	applications/programs/find_center3d/Makefile \
	applications/programs/foms_evaluate/Makefile \
	applications/programs/fourier_filter/Makefile \
	applications/programs/fourier_transform/Makefile \
	applications/programs/header_apply/Makefile \
	applications/programs/header_assign/Makefile \
	applications/programs/header_extract/Makefile \
	applications/programs/header_print/Makefile \
	applications/programs/header_reset/Makefile \
	applications/programs/statistics/Makefile \
	applications/programs/histogram/Makefile \
	applications/programs/make_spectra/Makefile \
	applications/programs/mask/Makefile \
	applications/programs/mask_design/Makefile \
	applications/programs/micrograph_denoise/Makefile \
	applications/programs/micrograph_downsample/Makefile \
	applications/programs/micrograph_mark/Makefile \
	applications/programs/micrograph_scissor/Makefile \
	applications/programs/micrograph_window/Makefile \
	applications/programs/mirror/Makefile \
	applications/programs/ml_align2d/Makefile \
	applications/programs/ml_align2d_combine/Makefile \
	applications/programs/ml_refine3d/Makefile \
	applications/programs/morphology/Makefile \
	applications/programs/mpi_angular_continuous_assign/Makefile \
	applications/programs/mpi_angular_discrete_assign/Makefile \
	applications/programs/mpi_angular_projection_matching/Makefile \
	applications/programs/mpi_art/Makefile \
	applications/programs/mpi_ml_align2d/Makefile \
	applications/programs/mpi_ml_refine3d/Makefile \
	applications/programs/mpi_reconstruct_wbp/Makefile \
	applications/programs/mpi_run/Makefile \
	applications/programs/normalize/Makefile \
	applications/programs/operate/Makefile \
	applications/programs/phantom_create/Makefile \
	applications/programs/phantom_create_random/Makefile \
	applications/programs/phantom_simulate_microscope/Makefile \
	applications/programs/phantom_test_reconstruction/Makefile \
	applications/programs/phantom_transform/Makefile \
	applications/programs/precompute_sampling/Makefile \
	applications/programs/project/Makefile \
	applications/programs/psd_enhance/Makefile \
	applications/programs/range_adjust/Makefile \
	applications/programs/reconstruct_art/Makefile \
	applications/programs/reconstruct_wbp/Makefile \
	applications/programs/resolution_fsc/Makefile \
	applications/programs/resolution_ssnr/Makefile \
	applications/programs/reverse_endian/Makefile \
	applications/programs/rotate/Makefile \
	applications/programs/scale/Makefile \
	applications/programs/scale_pyramid/Makefile \
	applications/programs/selfile_compare/Makefile \
	applications/programs/selfile_copy/Makefile \
	applications/programs/selfile_delete/Makefile \
	applications/programs/selfile_move/Makefile \
	applications/programs/selfile_select/Makefile \
	applications/programs/selfile_split/Makefile \
	applications/programs/selfile_statistics/Makefile \
	applications/programs/separate_objects/Makefile \
	applications/programs/shift/Makefile \
	applications/programs/show/Makefile \
	applications/programs/sort_by_statistics/Makefile \
	applications/programs/symmetrize/Makefile \
	applications/programs/threshold/Makefile \
	applications/programs/try_symmetry/Makefile \
	applications/programs/volume_segment/Makefile \
	applications/programs/window/Makefile \
	Makefile \
	m4/Makefile )

