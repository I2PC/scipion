dnl @synopsis AC_CHECK_LIBTIFF
dnl
dnl Tests for presence of TIFF library. If libtiff has been found, 
dnl `HAVE_LIBTIFF' will be defined. If the libtiff requires additional
dnl libraries (usualy zlib for deflate and libjpeg for JPEG compression),
dnl they will be added to LIBS.
dnl
dnl The test also checks the presence of different codecs in libtiff:
dnl HAVE_TIFFINITZIP tells deflate codec has been compiled in,
dnl HAVE_TIFFINITJPEG tells JPEG codec has been compiled in,
dnl HAVE_TIFFINITSGILOG tells SGI logLuv codec has been compuiled in,
dnl HAVE_TIFFLZWENCODER tells the LZW encoder is present.
dnl
dnl @version $Id$
dnl @author Jan Prikryl <prikryl@acm.org>
AC_DEFUN(AC_CHECK_LIBTIFF, [
dnl zlib and libbjpeg might be needed later - check for them now; if they
dnl are not available, nothing happens - maybe libtiff does not need them 
  olibs=$LIBS
  AC_CHECK_LIB(jpeg, jpeg_start_compress)
  AC_CHECK_LIB(z, zlibVersion)
  zjlibs=$LIBS
dnl Now check for libtiff; try it without zlib and libjpeg first,
dnl in case this fails try it once more with zlib and libjpeg. We shall
dnl probably check for linking with zlib and libjpeg separately in future.
  have_tiff=no
  with_zj=no
  AC_MSG_CHECKING([for TIFFOpen in -ltiff])
  AC_CACHE_VAL(ac_cv_have_libtiff, [
    LIBS="-ltiff $olibs"
    AC_TRY_LINK_FUNC(TIFFOpen, [ have_tiff=yes ],
    [ if test "$olibs" != "$zjlibs" ; then
        LIBS="-ltiff $zjlibs"
        AC_TRY_LINK_FUNC(TIFFOpen, [ have_tiff=yes; with_zj=yes ])
      fi ])
    ac_cv_have_libtiff="have_tiff=$have_tiff with_zj=$with_zj"])
dnl retrieve the cached value
  eval "$ac_cv_have_libtiff"
dnl set the environment
  if test $have_tiff = yes ; then
    if test $with_zj = yes ; then
      LIBS="-ltiff $zjlibs"
dnl   libs have a trailing space which gets lost with echo 
      oes=`echo $olibs`
      res="yes, requires "`echo $zjlibs | sed "s%$oes%%g"`
    else
      LIBS="-ltiff $olibs"
      res="yes"
    fi
  else
    LIBS="$olibs"
    res="no"
  fi
  AC_MSG_RESULT([$res])
  AC_DEFINE(HAVE_LIBTIFF, 1, [Define to 1 if translation of program messages to the user's native language is requested.])
  AM_CONDITIONAL(TIFF, test $have_tiff = yes)
dnl check for codecs; we assume that these can both compress and decompress
  if test $have_tiff = yes ; then   
    AC_MSG_CHECKING([for available nonstandard codecs in -ltiff])
    AC_CACHE_VAL(ac_cv_tiff_codecs, [
      AC_TRY_LINK_FUNC(TIFFInitZIP,    [have_tiffzip=yes], 
[have_tiffzip=no ])
      ac_cv_tiff_codecs="have_tiffzip=$have_tiffzip "
      AC_TRY_LINK_FUNC(TIFFInitJPEG,   [have_tiffjpg=yes], 
[have_tiffjpg=no ])
      ac_cv_tiff_codecs=$ac_cv_tiff_codecs"have_tiffjpg=$have_tiffjpg "
      AC_TRY_LINK_FUNC(TIFFInitSGILog, [have_tiffsgi=yes], 
[have_tiffsgi=no ])
      ac_cv_tiff_codecs=$ac_cv_tiff_codecs"have_tiffsgi=$have_tiffsgi " ])
dnl retrieve the cached value
    eval "$ac_cv_tiff_codecs"
dnl print the result and set the variables
    res=""
    if test $have_tiffzip = yes ; then
      res=$res"zip "
      AC_DEFINE(HAVE_TIFFINITZIP, 1, [Define to 1 if translation of program messages to the user's native language is requested.])
    fi
    if test $have_tiffjpg = yes ; then
      res=$res"jpeg "
      AC_DEFINE(HAVE_TIFFINITJPEG, 1, [Define to 1 if translation of program messages to the user's native language is requested.])
    fi
    if test $have_tiffsgi = yes ; then
      res=$res"sgilog "
      AC_DEFINE(HAVE_TIFFINITSGILOG, 1, [Define to 1 if translation of program messages to the user's native language is requested.])
    fi
    if test "x$res" = x ; then
      res="none"
    fi
    AC_MSG_RESULT([$res])
dnl check for LZW encoder; LZW decoder shall be compiled in by default
    AC_MSG_CHECKING([for LZW encoder in -ltiff])
    AC_CACHE_VAL(ac_cv_tiff_lzwenc, [
      AC_TRY_RUN([#include <tiffio.h>
  	int main(void)
  	{
  	  uint16 compress;
  	  TIFF *tif = TIFFOpen("conftest.tif", "w");
  	  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
  	  TIFFGetField(tif, TIFFTAG_COMPRESSION, &compress);
  	  if (compress != COMPRESSION_LZW)
  	    exit(-1);
  	  exit(0);
  	  return 0;    
  	}], [ ac_cv_tiff_lzwenc=yes ], [ ac_cv_tiff_lzwenc=no ], [ 
        ac_cv_tiff_lzwenc=no ])
      ])
dnl output
    if test $ac_cv_tiff_lzwenc = yes ; then
      AC_DEFINE(HAVE_TIFFLZWENCODER, 1, [Define to 1 if translation of program messages to the user's native language is requested.])
    fi
    AC_MSG_RESULT([$ac_cv_tiff_lzwenc])
  fi
])
