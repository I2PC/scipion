dnl qt.m4
dnl Adapted to GePhex by Georg Seidel <georg.seidel@web.de>
dnl Changes made: 
dnl    - added support for Darwin (made shared library extension a
dnl      new parameter
dnl    - added check for libqt-mt
dnl    - added minimum version check
dnl    - replaced AC_ERROR with AC_MSG_RESULT
dnl    - moved evaluation of ACTION-IF-FOUND and ACTION-IF-NOT-FOUND
dnl      to the end
dnl    - added a lot more guess dirs
dnl    - adapted to qt4 only
dnl    - adapted to mac os x frameworks
dnl     
dnl Original version from Rik Hemsley:
dnl   Copyright (C) 2001 Rik Hemsley (rikkus) <rik@kde.org>


dnl AM_PATH_QT(MINIMUM_VERSION, LIBEXT, 
dnl            [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_QT],
[
AC_CHECKING([for Qt ...])

AC_LANG_SAVE
AC_LANG_CPLUSPLUS

saved_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
saved_LIBRARY_PATH="$LIBRARY_PATH"
saved_CXXFLAGS="$CXXFLAGS"
saved_LDFLAGS="$LDFLAGS"
saved_LIBS="$LIBS"

# the test looks inside the following files to find the qt headers, libs
# and binaries
GUESS_QT_INC_DIRS="$QTDIR/include $QTDIR/include/qt /usr/include /usr/include/qt  /usr/include/qt4 /usr/local/include /usr/local/include/qt /usr/local/include/qt4 /usr/X11R6/include/ /usr/X11R6/include/qt /usr/X11R6/include/X11/qt /usr/X11R6/include/qt4 /usr/lib/qt/include /usr/lib/qt4/include /usr/lib/qt-4/include /usr/lib/qt-4.0/include /usr/lib/qt-4.1/include /usr/lib/qt-4.2/include /sw/include/qt"

GUESS_QT_LIB_DIRS="$QTDIR/lib /usr/lib /usr/local/lib /usr/X11R6/lib /usr/local/qt/lib /usr/lib/qt/lib /usr/lib/qt4/lib /usr/lib/qt-4/lib /usr/lib/qt-4.0/lib /usr/lib/qt-4.1/lib /usr/lib/qt-4.2/lib /usr/lib/qt-4.3/lib /usr/lib/qt-4.4/lib /sw/lib"

GUESS_QT_BIN_DIRS="$QTDIR/bin /usr/bin /usr/local/bin /usr/local/bin/qt4 /usr/X11R6/bin /usr/lib/qt/bin /usr/lib/qt4/bin /usr/lib/qt-4/bin /usr/lib/qt-340/bin /usr/lib/qt-4.1/bin /usr/lib/qt-4.2/bin /usr/lib/qt-4.3/bin /usr/lib/qt-4.4/bin /sw/bin"

HAVE_QT=no
min_qt_version=ifelse([$1], ,4.0.0, $1)

AC_ARG_WITH([qt-libdir],
  [AC_HELP_STRING([--with-qt-libdir=PFX],
                  [Prefix where Qt library is installed (optional)])],
  qt_libdir="$withval",
  qt_libdir=""
)

AC_ARG_WITH([qt-incdir],
  [AC_HELP_STRING([--with-qt-incdir=PFX],
                  [Prefix where Qt includes are installed (optional)])],
  qt_incdir="$withval",
  qt_incdir=""
)

AC_ARG_WITH([qt-bindir],
  [AC_HELP_STRING([--with-qt-bindir=PFX],
                  [Prefix where moc and uic are installed (optional)])],
  qt_bindir="$withval",
  qt_bindir=""
)

AC_MSG_CHECKING([include path])

dnl If we weren't given qt_incdir, have a guess.

if test "x$qt_incdir" != "x"
then
  AC_MSG_RESULT([specified as $qt_incdir])
else

  for dir in $GUESS_QT_INC_DIRS
  do
    if test -e $dir/qobject.h || test -e $dir/Qt/qobject.h
    then
      qt_incdir="$dir -I$dir/QtCore -I$dir/Qt"
      AC_MSG_RESULT([assuming $dir])
      break
    fi

    if test -e $dir/../lib/QtCore.framework/Headers/qobject.h
    then
      qt_incdir="$dir/../lib/QtCore.framework/Headers -F$dir/../lib"
      AC_MSG_RESULT([using frameworks, assuming $dir])
      break
    fi

  done

  if test "x$qt_incdir" = "x"
  then
    AC_MSG_RESULT([not found])
  fi

fi

dnl If we weren't given qt_libdir, have a guess.

AC_MSG_CHECKING([library path])

if test "x$qt_libdir" != "x"
then
  AC_MSG_RESULT([specified as $qt_libdir])
else
  for dir in $GUESS_QT_LIB_DIRS
  do
    if test -e $dir/libQtCore.$2
    then
      qt_libdir=$dir
      qt_ld_flag="-lQtCore -lQtGui"
      AC_MSG_RESULT([assuming $dir])
      break
    fi

    if test -e $dir/libQtCore4.$2
    then
      qt_libdir=$dir
      qt_ld_flag="-lQtCore4 -lQtGui4"
      AC_MSG_RESULT([assuming $dir])
      break
    fi

    dnl Look for frameworks

    if test -e $dir/QtCore.framework
    then
      qt_libdir=$dir
      qt_ld_flag="-F$dir -framework QtCore -framework QtGui"
      AC_MSG_RESULT([using frameworks, assuming $dir])
      break
    fi
  done

  if test "x$qt_ld_flag" = "x"
  then
    AC_MSG_RESULT([not found])
  fi

fi


dnl If we weren't given qt_bindir, have a guess.

AC_MSG_CHECKING([binary directory])

if test "x$qt_bindir" != "x"
then
  AC_MSG_RESULT([specified as $qt_bindir])
else

  for dir in $GUESS_QT_BIN_DIRS
  do
    if test -x $dir/moc -a -x $dir/uic
    then
      qt_bindir=$dir
      AC_MSG_RESULT([assuming $dir])
      break
    fi
  done

  if test "x$qt_bindir" = "x"
  then
    AC_MSG_RESULT([not found])
  fi
fi

dnl ifelse is talked about in m4 docs

if test "x$qt_incdir" = "x"
then
  AC_MSG_RESULT([Can't find includes])
elif test "x$qt_ld_flag" = "x"
then
  AC_MSG_RESULT([Can't find library])
elif test "x$qt_bindir" = "x"
then
  AC_MSG_RESULT([Can't find moc and/or uic])
else
  HAVE_QT=yes
fi

LDFLAGS="$LDFLAGS -L$qt_libdir"
LIBS="$LIBS $qt_ld_flag"

CXXFLAGS="$CXXFLAGS -I$qt_incdir"

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$qt_libdir"

dnl we got this far, now start checking if we have the right version
if test "x$HAVE_QT" = "xyes"
then
  AC_MSG_CHECKING(for qt - version >= $min_qt_version)
	dnl now run a short C app that tells us if the version is ok or not
        rm -f conf.qttest
	AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <QtCore/qglobal.h>

int
main ()
{
  int major, minor, micro;
  char ver[50];

  system ("touch conf.qttest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  strncpy(ver, "$min_qt_version", sizeof(ver) - 1);
  if (sscanf(ver, "%d.%d.%d", &major, &minor, &micro) != 3)
  {
    printf ("%s, bad version string\n", "$min_qt_version");
    exit (1);
  }

  if ( QT_VERSION >= major*100 + minor*10 + micro )
  {
    printf("%s - ", QT_VERSION_STR);
    return 0;
  }
  else
  {
    printf("\n*** Version of QT found is too old: QT_VERSION is %i.\n",
            QT_VERSION);
    printf("*** Upgrade QT and remove the file config.cache (if it exists)\n");
    printf("*** before re-running configure\n");
    return 1;
  }
}
    ],
    [
	AC_MSG_RESULT(yes)
	HAVE_QT=yes
    ],
    [
	AC_MSG_RESULT(no)
        HAVE_QT=no
    ])
fi

found_qt="no"

if test "x$HAVE_QT" = "xyes"
then
  if test -f conf.qttest ; then
      found_qt="yes"
  else
      AC_MSG_CHECKING([Could not run QT test program, checking if a Qt program links...])

      AC_TRY_LINK([
       #include <qstring.h>
      ],
      [
       QString s("Hello, world!");
       qDebug(s.latin1());
      ],
      found_qt="yes"
      AC_MSG_RESULT([ok]),
      AC_MSG_RESULT([failed - check config.log for details])
      )
  fi

  if test "x$found_qt" = "xyes"
  then
    QT_CXXFLAGS="-I$qt_incdir"

    QT_LIBS="-L$qt_libdir $qt_ld_flag"
    MOC="$qt_bindir/moc"
    UIC="$qt_bindir/uic"
    HAVE_QT=yes
  else
    HAVE_QT=no
  fi

  if test "x$HAVE_QT" = "xyes"
  then
    ifelse([$3], , :, [$3])
  else
    ifelse([$4], , :, [$4])
  fi
  
fi

AC_SUBST(QT_CXXFLAGS)
AC_SUBST(QT_LIBS)
AC_SUBST(MOC)
AC_SUBST(UIC)

AC_LANG_RESTORE()

LD_LIBRARY_PATH="$saved_LD_LIBRARY_PATH"
LIBRARY_PATH="$saved_LIBRARY_PATH"
CXXFLAGS="$saved_CXXFLAGS"
LDFLAGS="$saved_LDFLAGS"
LIBS="$saved_LIBS"
rm -f conf.qttest
])
