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
import os
import sys
from install.funcs import Environment, progInPath

get = lambda x: os.environ.get(x, 'y').lower() in ['true', 'yes', 'y', '1']


env = Environment(args=sys.argv)

noOpencv = '--no-opencv' in sys.argv or not get('OPENCV')
noScipy = '--no-scipy' in sys.argv or not get('SCIPY')


#  ************************************************************************
#  *                                                                      *
#  *                              Libraries                               *
#  *                                                                      *
#  ************************************************************************

cmake = env.addLibrary(
    'cmake',
    tar='cmake-3.2.2.tgz',
    targets=[env.getBin('cmake')],
    commands=[('cd software/tmp/cmake-3.2.2; '
               './bootstrap --prefix=../.. --parallel=%d' % env.getProcessors(),
               'software/tmp/cmake-3.2.2/Makefile'),
              ('cd software/tmp/cmake-3.2.2; make install -j %d'
               % env.getProcessors(), 'software/bin/cmake')],
    default=False)

# In order to get both the float and double libraries of fftw
# we need to execute ./configure; make; make install twice
# see: http://www.fftw.org/fftw2_doc/fftw_6.html
fftw3 = env.addLibrary(
    'fftw3',
    tar='fftw-3.3.4.tgz',
    flags=['--enable-threads', '--enable-shared'],
    clean=True) # We need to clean to configure again with --enable-float
    
fftw3f = env.addLibrary(
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

sh_alignment = env.addLibrary(
    'sh_alignment',
    tar='sh_alignment.tgz',
    commands=[('cd software/tmp/sh_alignment; make install',
               'software/lib/python2.7/site-packages/sh_alignment/frm.py')],
    deps=[python, swig],
    default=False)

lapack = env.addLibrary(
    'lapack',
    tar='lapack-3.5.0.tgz',
    flags=['-DBUILD_SHARED_LIBS:BOOL=ON',
           '-DLAPACKE:BOOL=ON'],
    cmake=True,
    neededProgs=['gfortran'],
    default=False)

arpack = env.addLibrary(
    'arpack',
    tar='arpack-96.tgz',
    neededProgs=['gfortran'],
    commands=[('cd software/bin; ln -s $(which gfortran) f77',
               'software/bin/f77'),
              ('cd software/tmp/arpack-96; make all',
               'software/lib/libarpack.a')],
    default=False)
# See http://modb.oce.ulg.ac.be/mediawiki/index.php/How_to_compile_ARPACK

if get('CUDA'):
    opencvFlags = ['-DWITH_CUDA:BOOL=ON']
else:
    opencvFlags = ['-DWITH_CUDA:BOOL=OFF']
opencv = env.addLibrary(
    'opencv',
    tar='opencv-2.4.9.tgz',
    targets=[env.getLib('opencv_core')],
    flags=opencvFlags,
    cmake=True,
    default=not noOpencv)

# ---------- Libraries required by PyTom 

boost = env.addLibrary(
    'boost',
    tar='boost_1_56_0.tgz',
    commands=[('cp -rf software/tmp/boost_1_56_0/boost software/include/', 
               'software/include/boost')],
    default=False)

nfft3 = env.addLibrary(
    'nfft3',
    tar='nfft-3.2.3.tgz',
    deps=[fftw3],
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
    tar='numpy-1.8.1.tgz',
    deps=[lapack])

six = env.addModule(
    'six',
    tar='six-1.7.3.tgz',
    targets=['six-1.7.3*'])

dateutil = env.addModule(
    'dateutil',
    tar='python-dateutil-1.5.tgz',
    targets=['python_dateutil-1.5*'],
    deps=[setuptools, six])

pyparsing = env.addModule(
    'pyparsing',
    tar='pyparsing-2.0.2.tgz',
    targets=['pyparsing.py'])

matplotlib = env.addModule(
    'matplotlib',
    tar='matplotlib-1.3.1.tgz',
    targets=['matplotlib-1.3.1*'],
    numpyIncludes=True,
    deps=[numpy, png, dateutil, pyparsing])

psutil = env.addModule(
    'psutil',
    targets=['psutil-2.1.1*'],
    tar='psutil-2.1.1.tgz')

mpi4py = env.addModule(
    'mpi4py',
    tar='mpi4py-1.3.1.tgz')

scipy = env.addModule(
    'scipy',
    tar='scipy-0.14.0.tgz',
    default=not noScipy,
    deps=[lapack, numpy, matplotlib])

bibtexparser = env.addModule(
    'bibtexparser',
    tar='bibtexparser-0.5.tgz')

django = env.addModule(
    'django',
    tar='Django-1.5.5.tgz')

paramiko = env.addModule(
    'paramiko',
    tar='paramiko-1.14.0.tgz',
    default=False)

pillow = env.addModule(
    'Pillow',
    tar='Pillow-2.5.1.tgz',
    targets=['Pillow-2.5.1*'],
    deps=[setuptools, jpeg])

winpdb = env.addModule(
    'winpdb',
    tar='winpdb-1.4.8.tgz',
    targets=[env.getBin('winpdb')],
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
    targets=['lxml-3.4.1*'],
    libChecks=['libxml-2.0', 'libxslt'],
    deps=[], # libxml2, libxslt],
    incs=['/usr/include/libxml2'],
    default=False)
# libxml2 and libxslt are checked instead of compiled because
# they are so hard to compile right.

ipython = env.addModule(
    'ipython',
    tar='ipython-2.1.0.tar.gz',
    deps=[pyzmq, jinja2, tornado],
    default=False)

cython = env.addModule(
    'cython',
    tar='Cython-0.22.tgz',
    targets=['Cython-0.22*'],
    default=False)

cythongsl = env.addModule(
    'cythongsl',
    tar='CythonGSL-0.2.1.tgz',
    targets=['CythonGSL-0.2.1*'],
    default=False,
    deps=[cython])
# TODO: add checks for dependencies: GSL


#  ************************************************************************
#  *                                                                      *
#  *                       External (EM) Packages                         *
#  *                                                                      *
#  ************************************************************************

# 'commands' is a list of (command, [targets]) to run after installation.


env.addPackage('bsoft-1.8.8',
               tar='bsoft1_8_8_Fedora_12.tgz',
               default=False)

env.addPackage('bsoft-1.9.0',
               tar='bsoft1_9_0_Fedora_20.tgz',
               default=False)

env.addPackage('ctffind',
               tar='ctffind_V3.5.tgz',
               default=False)

env.addPackage('ctffind4',
               tar='ctffind_V4.0.15.tgz',
               default=False)

env.addPackage('summovie',
               tar='summovie_1.0.2.tgz',
               default=False)

env.addPackage('unblur',
               tar='unblur_1.0_150529.tgz',
               default=False)

env.addPackage('eman2.11',
               tar='eman2.11.linux64.tgz',
               commands=[('./eman2-installer', 
                          'eman2.bashrc')],
               default=False)

env.addPackage('eman2.12',
               tar='eman2.12.linux64.tgz',
               commands=[('./eman2-installer', 
                          'eman2.bashrc')],
               default=False)

env.addPackage('frealign',
               tar='frealign_v9.07.tgz',
               default=False)

libSuffix = env.getLibSuffix()

env.addPackage('pytom',
               tar='pytom-0.963beta-scipion.tgz',
               commands=[('./scipion_installer', 
                          ['pytomc/libs/libtomc/libs/libtomc.%s' % libSuffix] + 
                          ['pytomc/swigModules/_pytom_%s.%s' % (s, libSuffix) 
                           for s in ['mpi', 'freqweight', 'volume', 'fftplan', 'numpy']])],
               deps=[boost, fftw3, fftw3f, nfft3,
                     swig, lxml, numpy, scipy,
                     matplotlib, mpi4py, pillow],
               default=False)

env.addPackage('relion-1.4',
               tar='relion-1.4.tgz',
               commands=[('./INSTALL.sh -j %d' % env.getProcessors(),
                          ['relion_build.log',
                           'bin/relion'])],
               default=False)

env.addPackage('relion-1.4_float',
               tar='relion-1.4_float.tgz',
               commands=[('./INSTALL.sh -j %d' % env.getProcessors(),
                          ['relion_build.log',
                           'bin/relion'])],
               default=False)

env.addPackage('relion-1.3',
               tar='relion-1.3.tgz',
               commands=[('./INSTALL.sh -j %d' % env.getProcessors(),
                          ['relion_build.log',
                           'bin/relion'])],
               default=False)

env.addPackage('resmap',
               tar='resmap-1.1.5-scipion.tgz',
               deps=['scipy'],
               default=False)

env.addPackage('spider',
               tar='spider-web-21.13.tgz',
               neededProgs=['csh'],
               default=False)

env.addPackage('motioncorr',
               tar='motioncorr_v2.1.tgz',
               default=False)

env.addPackage('simple',
               tar='simple2.tgz',
               default=False)

env.addPackage('chimera',
               tar='chimera-1.10.1-linux_x86_64.tgz',
               targetDir='chimera-1.10.1',
               commands=[('./scipion_installer','bin/chimera')],
               default=False)

env.addPackage('dogpicker',
               tar='dogpicker-0.2.1.tgz',
               default=False)

env.addPackage('nma',
               tar='nma.tgz',
               commands=[('cd ElNemo; make; mv nma_* ..', 'nma_elnemo_pdbmat'),
                         ('cd NMA_cart; LDFLAGS=-L%s/lib make; mv nma_* ..' %
                          os.environ['SCIPION_SOFTWARE'], 'nma_diag_arpack')],
               deps=['arpack'],
               default=False)

cryoem = env.addPackage(
                'cryoem',
                tar='cryoem-1.0.tgz',
                default=False,
                pythonMod=True,
                deps=[numpy, scipy, matplotlib, cythongsl])

env.addPackage('gEMpicker_v1.1',
               tar='gEMpicker_v1.1.tgz',
               default=False)


env.execute()
