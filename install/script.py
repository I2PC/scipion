# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************

import sys
from install.funcs import Environment, progInPath



env = Environment(args=sys.argv)


#  ************************************************************************
#  *                                                                      *
#  *                              Libraries                               *
#  *                                                                      *
#  ************************************************************************

# In order to get both the float and double libraries of fftw
# we need to execute ./configure; make; make install twice
# see: http://www.fftw.org/fftw2_doc/fftw_6.html
env.addLibrary(
    'fftw3',
    tar='fftw-3.3.4.tgz',
    flags=['--enable-threads', '--enable-shared'],
    clean=True) # We need to clean to configure again with --enable-float
    
env.addLibrary(
    'fftw3f',
    tar='fftw-3.3.4.tgz',
    flags=['--enable-threads', '--enable-shared', '--enable-float'])

osBuildDir = 'tcl8.6.1/unix'
osFlags = ['--enable-threads']

tcl = env.addLibrary(
    'tcl',
    tar='tcl8.6.1-src.tgz',
    buildDir=osBuildDir,
    targets=[env.getLib('tcl8.6')],
    flags=osFlags)

osBuildDir = 'tk8.6.1/unix'
osFlags = ['--enable-threads']

tk = env.addLibrary(
    'tk',
    tar='tk8.6.1-src.tgz',
    buildDir=osBuildDir,
    targets=[env.getLib('tk8.6')],
    libChecks=['xft'],
    flags=osFlags,
    deps=[tcl])


# Special case: tk does not make the link automatically, go figure.
tk_wish = env.addTarget('tk_wish')
tk_wish.addCommand('ln -v -s wish8.6 wish',
                   targets='software/bin/wish',
                   cwd='software/bin')

zlib = env.addLibrary(
    'zlib',
    targets=[env.getLib('z')],
    tar='zlib-1.2.8.tgz',
    configTarget='zlib.pc')

jpeg = env.addLibrary(
    'jpeg',
    tar='libjpeg-turbo-1.3.1.tgz',
    flags=['--without-simd'])
    #flags=([] if progInPath('nasm') else ['--without-simd']))

png = env.addLibrary(
    'png',
    tar='libpng-1.6.16.tgz',
    deps=[zlib])

tiff = env.addLibrary(
     'tiff',
     tar='tiff-3.9.4.tgz',
     deps=[zlib, jpeg])

sqlite = env.addLibrary(
    'sqlite3',
    tar='sqlite-3.6.23.tgz',
    flags=['CPPFLAGS=-w',
           'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1'])

hdf5 = env.addLibrary(
     'hdf5',
     tar='hdf5-1.8.14.tgz',
     flags=['--enable-cxx', '--enable-shared'],
     targets=[env.getLib('hdf5'), env.getLib('hdf5_cpp')],
     configAlways=True,
     deps=[zlib])

python = env.addLibrary(
    'python',
    tar='Python-2.7.8.tgz',
    targets=[env.getLib('python2.7'), env.getBin('python')],
    flags=['--enable-shared'],
    deps=[sqlite, tk, zlib])

lapack = env.addLibrary(
    'lapack',
    tar='lapack-3.5.0.tgz',
    flags=['-DBUILD_SHARED_LIBS:BOOL=ON',
           '-DLAPACKE:BOOL=ON'],
    cmake=True, 
    default=False)

opencv = env.addLibrary(
    'opencv',
    tar='opencv-2.4.9.tgz',
    targets=[env.getLib('opencv_core')],
    cmake=True,
    default=False)

# ---------- Libraries required by PyTom 
pcre = env.addLibrary(
    'pcre',
    tar='pcre-8.36.tgz',
    targets=[env.getBin('pcretest')],
    default=False)

swig = env.addLibrary(
    'swig',
    tar='swig-3.0.2.tgz',
    targets=[env.getBin('swig')],
    makeTargets=['Source/Swig/tree.o'],
    deps=[pcre],
    default=False)

boost = env.addLibrary(
    'boost',
    tar='boost_1_56_0.tgz',
    commands=[('cp -rf software/tmp/boost_1_56_0/boost software/include/', 
               'software/include/boost')],
    default=False)


#  ************************************************************************
#  *                                                                      *
#  *                           Python Modules                             *
#  *                                                                      *
#  ************************************************************************

# The flag '--old-and-unmanageable' used in some modules avoids
# creating a single Python egg. That way the modules create a full
# directory with the name of package, and we use that as a target.

setuptools = env.addModule(
    'setuptools',
    tar='setuptools-5.4.1.tgz',
    targets=['setuptools.pth'])

scons = env.addModule(
    'scons',
    targets=[env.getBin('scons')],
    tar='scons-2.3.4.tgz')

numpy = env.addModule(
    'numpy',
    tar='numpy-1.8.1.tgz')

six = env.addModule(
    'six',
    tar='six-1.7.3.tgz',
    targets=['six.py'],
    flags=['--old-and-unmanageable'])

dateutil = env.addModule(
    'dateutil',
    tar='python-dateutil-1.5.tgz',
    flags=['--old-and-unmanageable'],
    deps=[setuptools, six])

pyparsing = env.addModule(
    'pyparsing',
    targets=['pyparsing.py'],
    tar='pyparsing-2.0.2.tgz')

matplotlib = env.addModule(
    'matplotlib',
    tar='matplotlib-1.3.1.tgz',
    flags=['--old-and-unmanageable'],
    deps=[numpy, png, dateutil, pyparsing])

env.addModule(
    'psutil',
    tar='psutil-2.1.1.tgz',
    flags=['--old-and-unmanageable'])

env.addModule(
    'mpi4py',
    tar='mpi4py-1.3.1.tgz')

env.addModule(
    'scipy',
    tar='scipy-0.14.0.tgz',
    default=False,
    deps=[lapack, numpy, matplotlib])

env.addModule(
    'bibtexparser',
    tar='bibtexparser-0.5.tgz')

django = env.addModule(
    'django',
    tar='Django-1.5.5.tgz')

env.addModule(
    'paramiko',
    tar='paramiko-1.14.0.tgz',
    default=False)

env.addModule(
    'Pillow',
    tar='Pillow-2.5.1.tgz',
    targets=['PIL'],
    flags=['--old-and-unmanageable'],
    deps=[setuptools, jpeg])

env.addModule(
    'winpdb',
    tar='winpdb-1.4.8.tgz',
    default=False)

pyzmq = env.addModule(
    'pyzmq',
    tar='pyzmq-2.2.0.1.tar.gz',
    default=False)

jinja2 = env.addModule(
    'jinja2',
    tar='Jinja2-2.7.3.tar.gz',
    default=False)

tornado = env.addModule(
    'tornado',
    tar='tornado-4.0.2.tar.gz',
    default=False)

lxml = env.addModule(
    'lxml',
    tar='lxml-3.4.1.tgz',
    libChecks=['libxml-2.0', 'libxslt'],
    deps=[], # libxml2, libxslt],
    default=False)
# libxml2 and libxslt are checked instead of compiled because
# they are so hard to compile right.

env.addModule(
    'ipython',
    tar='ipython-2.1.0.tar.gz',
    deps=[pyzmq, jinja2, tornado],
    default=False)



#  ************************************************************************
#  *                                                                      *
#  *                       External (EM) Packages                         *
#  *                                                                      *
#  ************************************************************************

# extraActions is a list of (target, command) to run after installation.

env.addPackage2('xmipp',
               tar='xmipp_scipion.tgz',
               buildDir='xmipp_scipion',
               reqs={'mpi': 'cxx',
                     'freetype': 'cxx',
                     'X11': 'cxx',
                     'png': 'cxx',
                     'ncurses': 'cxx',
                     'ssl': 'cxx',
                     'readline': 'cxx'},
               #deps=[opencv],
               default=False)
# In case you want to install an older version of Xmipp, you can use
# the extraActions parameter instead of using its own SConscript, like this:
# 
#               extraActions=[('xmipp.bashrc',
#                             './install.sh --unattended=true --gui=false -j %s'
#                              % GetOption('num_jobs'))],

env.addPackage('bsoft',
               tar='bsoft1_8_8_Fedora_12.tgz',
               default=False)

env.addPackage('ctffind',
               tar='ctffind_V3.5.tgz',
               default=False)

env.addPackage('ctffind4',
               tar='ctffind_V4.0.13.tgz',
               default=False)

env.addPackage('eman',
               tar='eman2.1.linux64.tgz',
               commands=[('./eman2-installer', 
                          'software/em/eman/eman2.bashrc')],
               default=False)

env.addPackage('frealign',
               tar='frealign_v9.07.tgz',
               default=False)


# env.addPackage2('pytom',
#                tar='pytom_develop0.962.tgz',
#                extraActions=[('pytomc/libs/libtomc/libs/libtomc.%s' % libSuffix,
#                               'PATH=%s/software/bin:%s '
#                               'LD_LIBRARY_PATH=%s/software/lib:%s '
#                               'MPILIBDIR=%s MPIINCLUDEDIR=%s SCIPION_HOME=%s '
#                               './scipion_installer'
#                               % (shome, environ.get('PATH', ''),
#                                  shome, environ.get('LD_LIBRARY_PATH', ''),
#                                  env['MPI_LIBDIR'], env['MPI_INCLUDE'], shome))],
#                deps=[boost, 'fftw3', 'fftw3f', swig, lxml],
#                default=False)

env.addPackage2('relion',
               tar='relion-1.3.tgz',
               extraActions=[
                   ('relion_build.log', './INSTALL.sh -j 5'
                    )],
               default=False)

env.addPackage('resmap',
               tar='resmap-1.1.5-scipion.tgz',
               deps=['scipy'],
               default=False)

env.addPackage2('spider',
               tar='spider-web-21.13.tgz',
               neededProgs=['csh'],
               default=False)

env.addPackage2('motioncorr',
               tar='motioncorr_v2.1.tgz',
               default=False)

env.addPackage2('simple',
               tar='simple2.tgz',
               default=False)


env.execute()
