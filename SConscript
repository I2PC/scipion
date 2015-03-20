#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:     I. Foche Perez (ifoche@cnb.csic.es)
# *              J. Burguet Castell (jburguet@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************


# First import the environment (this comes from SConstruct)
Import('env')

from os import environ

#  ************************************************************************
#  *                                                                      *
#  *                              Libraries                               *
#  *                                                                      *
#  ************************************************************************

# First, we check the compilers
if not env.GetOption('clean'):
    env = env.CompilerConfig()

# We might want to add freetype and make tcl depend on it. That would be:
# freetype = env.AddLibrary(
#     'freetype',
#     tar='freetype-2.5.3.tgz',
#     autoConfigTarget='config.mk')
# But because freetype's compilation is a pain, it's better to use whatever
# version is in the system.

# fftw = env.AddLibrary(
#     'fftw3',
#     tar='fftw-3.3.4.tgz',
#     targets=[File('#software/lib/libfftw3.so').abspath],
#     flags=['--enable-threads', '--enable-shared'],
#     clean=[Dir('#software/tmp/fftw-3.3.4')])

# I needed to add two compilations of fftw to get float extensions
fftw1 = env.AddLibrary(
    'fftw1',
    tar='fftw-3.3.4.tgz',
    targets=[File('#software/lib/libfftw3.so').abspath],
    flags=['--enable-threads', '--enable-shared'],
    clean=[Dir('#software/tmp/fftw-3.3.4')])
    
fftw2 = env.AddLibrary(
    'fftw2',
    tar='fftw-3.3.4f.tgz',
    targets=[File('#software/lib/libfftw3f.so').abspath],
    flags=['--enable-threads', '--enable-shared', '--enable-float'],
    clean=[Dir('#software/tmp/fftw-3.3.4f')])

fftw = env.Alias('fftw3', [fftw1, fftw2])
env.Default(fftw)

tcl = env.AddLibrary(
    'tcl',
    tar='tcl8.6.1-src.tgz',
    buildDir='tcl8.6.1/unix',
    targets=[File('#software/lib/libtcl8.6.so').abspath],
    flags=['--enable-threads'],
    clean=[Dir('#software/tmp/tcl8.6.1').abspath])

tk = env.AddLibrary(
    'tk',
    tar='tk8.6.1-src.tgz',
    buildDir='tk8.6.1/unix',
    targets=[File('#software/lib/libtk8.6.so').abspath],
    libChecks=['xft'],
    flags=['--enable-threads'],
    deps=[tcl],
    clean=[Dir('#software/tmp/tk8.6.1').abspath])

tk_wish = env.Command(
    'software/bin/wish',
    'software/bin/wish8.6',
    Action('ln -v -s wish8.6 wish', 'Linking wish8.6 to wish in software/bin',
           chdir='software/bin'))
env.Default(tk_wish)
env.Depends(tk_wish, tk)
# Special case: tk does not make the link automatically, go figure.

zlib = env.AddLibrary(
    'zlib',
    tar='zlib-1.2.8.tgz',
    targets=[File('#software/lib/libz.so').abspath],
    addPath=False,
    autoConfigTargets='zlib.pc')

jpeg = env.AddLibrary(
    'jpeg',
    tar='libjpeg-turbo-1.3.1.tgz',
    targets=[File('#software/lib/libjpeg.so').abspath],
    flags=['--without-simd'])
    #flags=([] if env.ProgInPath('nasm') else ['--without-simd']))

png = env.AddLibrary(
    'png',
    tar='libpng-1.6.16.tgz',
    deps=[zlib])

tiff = env.AddLibrary(
     'tiff',
     tar='tiff-3.9.4.tgz',
     targets=[File('#software/lib/libtiff.so').abspath],
     deps=[zlib])

sqlite = env.AddLibrary(
    'sqlite3',
    tar='sqlite-3.6.23.tgz',
    targets=[File('#software/lib/libsqlite3.so').abspath],
    flags=['CPPFLAGS=-w',
           'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1'])

hdf5 = env.AddLibrary(
     'hdf5',
     tar='hdf5-1.8.14.tgz',
     flags=['--enable-cxx'],
     targets=[File('#software/lib/libhdf5.so').abspath, 
              File('#software/lib/libhdf5_cpp.so').abspath],
     deps=[zlib])

python = env.AddLibrary(
    'python',
    tar='Python-2.7.8.tgz',
    targets=[File('#software/lib/libpython2.7.so').abspath,
             File('#software/bin/python').abspath],
    flags=['--enable-shared'],
    deps=[sqlite, tk, zlib])

libxml2 = env.AddLibrary(
    'libxml2',
    tar='libxml2-2.9.2.tgz',
    targets=['lib/libxml2.so'],
    deps=[python],
    default=False)

libxslt = env.AddLibrary(
    'libxslt',
    tar='libxslt-1.1.28.tgz',
    targets=['lib/libxslt.so'],
    deps=[libxml2],
    default=False)
# This library is pretty complicated to compile right. For the moment,
# we will require it to be in the system (with development headers)
# for anyone that wants to use it, instead of compiling it automatically.

pcre = env.AddLibrary(
    'pcre',
    tar='pcre-8.36.tgz',
    targets=[File('#software/bin/pcretest').abspath],
    default=False)

swig = env.AddLibrary(
    'swig',
    tar='swig-3.0.2.tgz',
    targets=[File('#software/bin/swig').abspath],
    makeTargets='swig install',
    deps=[pcre],
    default=False)
# We have to add the "makeTargets" part because swig needs to call
# "make" before "make install". Horrible.

env.AddLibrary(
    'parallel',
    tar='parallel-20140922.tgz',
    targets=[File('#software/bin/parallel').abspath],
    deps=[zlib])

shome = env['SCIPION_HOME']  # short notation, we use it quite a lot
boost_headers_only = env.ManualInstall(
    'boost_headers_only',
    tar='boost_1_56_0.tgz',
    extraActions=[
        ('%s/software/include/boost' % shome,
         'cp -rf boost %s/software/include' % shome)],
    default=False)

lapack = env.ManualInstall(
    'lapack',
    tar='lapack-3.5.0.tgz',
    neededProgs=['cmake'],
    extraActions=[
        ('%s/software/tmp/lapack-3.5.0/Makefile' % shome,
         'cmake -DBUILD_SHARED_LIBS:BOOL=ON -DLAPACKE:BOOL=ON '
         '-DCMAKE_INSTALL_PREFIX:PATH=%s/software .' % shome),
        ('%s/software/lib/liblapack.so' % shome,
         'make install')],
    default=False)

opencv = env.ManualInstall(
    'opencv',
    tar='opencv-2.4.9.tgz',
    neededProgs=['cmake'],
    extraActions=[
        ('%s/software/tmp/opencv-2.4.9/Makefile' % shome,
         'cmake -DCMAKE_INSTALL_PREFIX:PATH=%s/software .' % shome),
        ('%s/software/lib/libopencv_core.so' % shome,
         'make install')],
    default=False)


#  ************************************************************************
#  *                                                                      *
#  *                           Python Modules                             *
#  *                                                                      *
#  ************************************************************************

# Helper function to include the python dependency automatically.
def addModule(*args, **kwargs):
    kwargs['deps'] = kwargs.get('deps', []) + python
    return env.AddModule(*args, **kwargs)


# The flag '--old-and-unmanageable' used in some modules avoids
# creating a single Python egg. That way the modules create a full
# directory with the name of package, and we use that as a target.

setuptools = addModule(
    'setuptools',
    tar='setuptools-5.4.1.tgz',
    targets=['setuptools.pth'])

numpy = addModule(
    'numpy',
    tar='numpy-1.8.1.tgz')

six = addModule(
    'six',
    tar='six-1.7.3.tgz',
    targets=['six.py'],
    flags=['--old-and-unmanageable'])

dateutil = addModule(
    'dateutil',
    tar='python-dateutil-1.5.tgz',
    flags=['--old-and-unmanageable'],
    deps=[setuptools, six])

pyparsing = addModule(
    'pyparsing',
    targets=['pyparsing.py'],
    tar='pyparsing-2.0.2.tgz')

matplotlib = addModule(
    'matplotlib',
    tar='matplotlib-1.3.1.tgz',
    flags=['--old-and-unmanageable'],
    deps=[numpy, png, dateutil, pyparsing])

addModule(
    'psutil',
    tar='psutil-2.1.1.tgz',
    flags=['--old-and-unmanageable'])

addModule(
    'mpi4py',
    tar='mpi4py-1.3.1.tgz')

addModule(
    'scipy',
    tar='scipy-0.14.0.tgz',
    default=False,
    deps=[lapack, numpy, matplotlib])

addModule(
    'bibtexparser',
    tar='bibtexparser-0.5.tgz')

django = addModule(
    'django',
    tar='Django-1.5.5.tgz')

addModule(
    'paramiko',
    tar='paramiko-1.14.0.tgz',
    default=False)

addModule(
    'Pillow',
    tar='Pillow-2.5.1.tgz',
    targets=['PIL'],
    flags=['--old-and-unmanageable'],
    deps=[setuptools, jpeg])

addModule(
    'winpdb',
    tar='winpdb-1.4.8.tgz',
    default=False)

pyzmq = addModule(
    'pyzmq',
    tar='pyzmq-2.2.0.1.tar.gz',
    default=False)

jinja2 = addModule(
    'jinja2',
    tar='Jinja2-2.7.3.tar.gz',
    default=False)

tornado = addModule(
    'tornado',
    tar='tornado-4.0.2.tar.gz',
    default=False)

lxml = addModule(
    'lxml',
    tar='lxml-3.4.1.tgz',
    libChecks=['libxml-2.0', 'libxslt'],
    deps=[], # libxml2, libxslt],
    default=False)
# libxml2 and libxslt are checked instead of compiled because
# they are so hard to compile right.

addModule(
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

env.AddPackage('xmipp',
               tar='xmipp_scipion.tgz',
               buildDir='xmipp_scipion',
               reqs={'mpi': 'cxx',
                     'freetype': 'cxx',
                     'X11': 'cxx',
                     'png': 'cxx',
                     'ncurses': 'cxx',
                     'ssl': 'cxx',
                     'readline': 'cxx'},
               default=False)
# In case you want to install an older version of Xmipp, you can use
# the extraActions parameter instead of using its own SConscript, like this:
# 
#               extraActions=[('xmipp.bashrc',
#                             './install.sh --unattended=true --gui=false -j %s'
#                              % GetOption('num_jobs'))],

env.AddPackage('bsoft',
               tar='bsoft1_8_8_Fedora_12.tgz',
               default=False)

env.AddPackage('ctffind',
               tar='ctffind_V3.5.tgz',
               default=False)

env.AddPackage('ctffind4',
               tar='ctffind_V4.0.13.tgz',
               default=False)

env.AddPackage('eman',
               tar='eman2.1.linux64.tgz',
               extraActions=[('eman2.bashrc', './eman2-installer')],
               default=False)

env.AddPackage('frealign',
               tar='frealign_v9.07.tgz',
               default=False)

env.AddPackage('pytom',
               tar='pytom_develop0.962.tgz',
               extraActions=[('pytomc/libs/libtomc/libs/libtomc.so',
                              'PATH=%s/software/bin:%s '
                              'LD_LIBRARY_PATH=%s/software/lib:%s '
                              'MPILIBDIR=%s MPIINCLUDEDIR=%s SCIPION_HOME=%s '
                              './scipion_installer'
                              % (shome, environ.get('PATH', ''),
                                 shome, environ.get('LD_LIBRARY_PATH', ''),
                                 env['MPI_LIBDIR'], env['MPI_INCLUDE'], shome))],
               deps=[boost_headers_only, fftw, swig, lxml],
               default=False)

env.AddPackage('relion',
               tar='relion-1.3.tgz',
               extraActions=[
                   ('relion_build.log', './INSTALL.sh -j %s'
                    % GetOption('num_jobs'))],
               default=False)

env.AddPackage('spider',
               tar='spider-web-21.13.tgz',
               neededProgs=['csh'],
               default=False)

env.AddPackage('motioncorr',
               tar='motioncorr_v2.1.tgz',
               default=False)

env.AddPackage('simple',
               tar='simple2.tgz',
               default=False)
