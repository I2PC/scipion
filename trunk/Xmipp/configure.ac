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
AM_INIT_AUTOMAKE(Xmipp,0.9)
AC_PROG_LIBTOOL
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_LANG(C++)

# Checking for --devel flag. If present, programs are buit with prefix devel, if not
# they will use prefix xmipp

AC_ARG_WITH([devel],
[  --devel                 build development versions],AC_SUBST(XMIPP_PREFIX,[devel_]),AC_SUBST(XMIPP_PREFIX,[xmipp_]))

# checking for the existence of VTK on regular PATH
#AC_MSG_CHECKING("VTK library dir")
#AC_ARG_WITH(vtk-libdir,[  --with-vtk-libdir=DIR           location of the 
#vtk lib directory],
#	[LDFLAGS="${LDFLAGS} -L${withval} -lvtkImaging  -lvtkGraphics 
#-lvtkFiltering -lvtkIO -lvtkjpeg -lvtkpng -lvtkzlib -lvtkCommon -dl" ])

#AC_ARG_WITH(vtk-incdir,[  --with-vtk-incdir=DIR           location of the 
#vtk include directory],
#	[CPPFLAGS="${CPPFLAGS} -I${withval}"])

#AC_CHECK_LIB(vtkCommon, main,,
#	[AC_MSG_ERROR([VTK not found; Check your installation or try 
#./configure --with-vtk-libdir=/libdir/for/vtk])])
#AC_CHECK_LIB(vtkImaging, main,,
#	[AC_MSG_ERROR([VTK not found; Check your installation or try 
#./configure --with-vtk-libdir=/libdir/for/vtk])])
#AC_CHECK_LIB(vtkGraphics, main,,
#	[AC_MSG_ERROR([VTK not found; Check your installation or try 
#./configure --with-vtk-libdir=/libdir/for/vtk])])
#AC_CHECK_HEADERS([vtkImageFFT.h])


#Checking for TIFF library as well as some of its capabilities
AC_CHECK_LIBTIFF


#Checking for the presence and usability of MPI library
ACX_MPI

#checking for Qt
BNV_HAVE_QT

# If Qt is present then set HAVE_QT to 1 in config.h file
# else print error and exit 

if have_qt=yes
then
	AC_DEFINE([HAVE_QT],1,Define HAVE_QT to 1 if Qt is present and usable)
else
	AC_MSG_ERROR([Qt Not found, check your installation.])
fi

#checking for OpenGL
MDL_HAVE_OPENGL

if test -z "$GL_LIBS"
then
	AC_MSG_RESULT([GL libs not found, OpenGL disabled.])
else
	if test -z "$GL_CFLAGS"
	then
		AC_MSG_RESULT([GL headers not found, OpenGL disabled.])
	else
		AC_MSG_RESULT([OpenGL flags: $GL_CFLAGS])
		AC_MSG_RESULT([OpenGL libs: $GL_LIBS])
		PHAN3D="phan3D"
		AC_SUBST([PHAN3D])
	fi
fi


# Checks for the presence of some system libs 

AC_CHECK_LIB(Xext, main)
AC_CHECK_LIB(Xt, main)
AC_CHECK_LIB(X11, main)
AC_CHECK_LIB(pthread, main)

# Generates Makefiles for the different source directories.

AC_OUTPUT( Lib/Makefile \
	Lib/XmippData/Makefile \
	Lib/Reconstruction/Makefile \
	Lib/Classification/Makefile \
	Lib/XmippInterface/Makefile \
	Lib/OldXmipp/Makefile \
	Lib/XmippGraphics/Makefile \
	Applications/Makefile \
	Applications/Src/Makefile \
	Applications/Src/2DInterpolation/Makefile \
	Applications/Src/3DInterpolation/Makefile \
	Applications/Src/Adapt_for_Spider/Makefile \
	Applications/Src/Add_noise/Makefile \
	Applications/Src/Adjust_CTF/Makefile \
	Applications/Src/Adjust_surface/Makefile \
	Applications/Src/Align2D/Makefile \
	Applications/Src/Align3D/Makefile \
	Applications/Src/Ang_distribution/Makefile \
	Applications/Src/Angular_distance/Makefile \
	Applications/Src/Angular_predict/Makefile \
	Applications/Src/Angular_refinement/Makefile \
	Applications/Src/Are_different/Makefile \
	Applications/Src/Art/Makefile \
	Applications/Src/Assign_angles/Makefile \
	Applications/Src/Assign_CTF/Makefile \
	Applications/Src/Assign_multiple_CTFs/Makefile \
	Applications/Src/BatchSOM/Makefile \
	Applications/Src/CCP4Spider/Makefile \
	Applications/Src/CodeStatis/Makefile \
	Applications/Src/Combine_stats/Makefile \
	Applications/Src/CorrectPhase/Makefile \
	Applications/Src/Correlation/Makefile \
	Applications/Src/CpSel/Makefile Applications/Src/Create_phantom/Makefile \
	Applications/Src/CTF_profile/Makefile \
	Applications/Src/Data2Img/Makefile \
	Applications/Src/Data2Vol/Makefile \
	Applications/Src/Deform/Makefile \
	Applications/Src/Denoising/Makefile \
	Applications/Src/Denoising_micrograph/Makefile \
	Applications/Src/DocFile_Histogram/Makefile \
	Applications/Src/DownSample/Makefile \
	Applications/Src/Draw_surface/Makefile \
	Applications/Src/DWT/Makefile \
	Applications/Src/Euler/Makefile \
	Applications/Src/Evaluate/Makefile \
	Applications/Src/Evaluate_FSCs/Makefile \
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
	Applications/Src/Mark/Makefile \
	Applications/Src/Mask/Makefile \
	Applications/Src/Microscope/Makefile \
	Applications/Src/Morphology/Makefile \
	Applications/Src/MPIArt/Makefile \
	Applications/Src/MPI_WBP/Makefile \
	Applications/Src/MvSel/Makefile \
	Applications/Src/Normalize/Makefile \
	Applications/Src/Operate/Makefile \
	Applications/Src/PCA/Makefile \
	Applications/Src/Phan3D/Makefile \
	Applications/Src/Phantom_transform/Makefile \
	Applications/Src/Project/Makefile \
	Applications/Src/Projection_Neighbourhood/Makefile \
	Applications/Src/ProjectPCA/Makefile \
	Applications/Src/Radon_Transform/Makefile \
	Applications/Src/RandomMapping/Makefile \
	Applications/Src/Random_phantom/Makefile \
	Applications/Src/Raw22spi/Makefile \
	Applications/Src/Recons_test/Makefile \
	Applications/Src/Reverse_endian/Makefile \
	Applications/Src/Resolution/Makefile \
	Applications/Src/RmSel/Makefile \
	Applications/Src/Rotate/Makefile \
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
	Applications/Src/Spots2RealSpace2D/Makefile \
	Applications/Src/SSNR/Makefile \
	Applications/Src/Statis/Makefile \
	Applications/Src/Surface/Makefile \
	Applications/Src/Symmetrize/Makefile \
	Applications/Src/Taxa_tilt/Makefile \
	Applications/Src/Threshold/Makefile \
	Applications/Src/Tiff2Raw/Makefile \
	Applications/Src/ToPDB/Makefile \
	Applications/Src/ToPhantom/Makefile \
	Applications/Src/Visualize_FFT/Makefile \
	Applications/Src/Vol2Data/Makefile \
	Applications/Src/VolSampler/Makefile \
	Applications/Src/Voxels22Blobs/Makefile \
	Applications/Src/Window/Makefile \
	Applications/Src/Window_micrograph/Makefile \
	Applications/Src/XMask/Makefile \
	Makefile \
	m4/Makefile \
	po/Makefile )
		
