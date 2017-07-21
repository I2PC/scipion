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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
import sys
from install.funcs import Environment, progInPath

get = lambda x: os.environ.get(x, 'y').lower() in ['true', 'yes', 'y', '1']


env = Environment(args=sys.argv)

noOpencv = '--no-opencv' in sys.argv or not get('OPENCV')
noScipy = '--no-scipy' in sys.argv or not get('SCIPY')


#  *******************************
#  *  PATHS
#  *******************************
# GET the real path where scipion is installed
SCIPION = env._args[0]
SCIPION = os.path.realpath(SCIPION)
SCIPION = os.path.dirname(SCIPION)
SCIPION = os.path.abspath(SCIPION)

SW = SCIPION + '/software'
SW_BIN = SW + '/bin'
SW_LIB = SW + '/lib'
SW_INC = SW + '/include'

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
    opencvFlags = ['-DWITH_FFMPEG=OFF -DWITH_CUDA:BOOL=ON']
else:
    opencvFlags = ['-DWITH_FFMPEG=OFF -DWITH_CUDA:BOOL=OFF']
opencv = env.addLibrary(
    'opencv',
    tar='opencv-2.4.13.tgz',
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

sklearn = env.addModule(
    'sklearn',
    tar='scikit-learn-0.17.tar.gz',
    default=False,
    deps=[scipy, numpy, cython])



#  ************************************************************************
#  *                                                                      *
#  *                       External (EM) Packages                         *
#  *                                                                      *
#  ************************************************************************

# 'commands' is a list of (command, [targets]) to run after installation.


env.addPackage('bsoft', version='1.8.8',
               tar='bsoft1_8_8_Fedora_12.tgz')

env.addPackage('bsoft', version='1.9.0',
               tar='bsoft1_9_0_Fedora_20.tgz')

env.addPackage('ctffind', version='3.6',
               tar='ctffind_V3.5.tgz')

env.addPackage('ctffind4', version='4.0.15',
               tar='ctffind_V4.0.15.tgz')

env.addPackage('ctffind4', version='4.1.5',
               tar='ctffind_V4.1.5.tgz')

env.addPackage('ctffind4', version='4.1.8',
               tar='ctffind_V4.1.8.tgz')

env.addPackage('summovie', version='1.0.2',
               tar='summovie_1.0.2.tgz')

env.addPackage('unblur', version='1.0.15',
               tar='unblur_1.0_150529.tgz')

env.addPackage('unblur', version='1.0.2',
               tar='unblur_1.0.2.tgz')

eman2_commands = [('./eman2-installer',
                   'eman2.*rc')]

env.addPackage('eman', version='2.11',
               tar='eman2.11.linux64.tgz',
               commands=eman2_commands)

env.addPackage('eman', version='2.12',
               tar='eman2.12.linux64.tgz',
               commands=eman2_commands)

env.addPackage('frealign', version='9.07',
               tar='frealign_v9.07.tgz')

relion_commands = [('./INSTALL.sh -j %d' % env.getProcessors(),
                          ['relion_build.log',
                           'bin/relion_refine'])]

env.addPackage('relion', version='1.3',
               tar='relion-1.3.tgz',
               commands=relion_commands)

env.addPackage('relion', version='1.4',
               tar='relion-1.4.tgz',
               commands=relion_commands)

env.addPackage('relion', version='1.4f',
               tar='relion-1.4_float.tgz',
               commands=relion_commands)

# Define FFTW3 path variables
relion_vars = [('FFTW_LIB', SW_LIB),
               ('FFTW_INCLUDE', SW_INC)]

relion2_commands = [('cmake -DGUI=OFF -DCMAKE_INSTALL_PREFIX=./ .', []),
                    ('make -j %d' % env.getProcessors(), ['bin/relion_refine'])]

env.addPackage('relion', version='2.0',
               tar='relion-2.0.4.tgz',
               commands=relion2_commands,
               updateCuda=True,
               vars=relion_vars)

env.addPackage('localrec', version='1.1.0',
               tar='localrec-1.1.0.tgz')

env.addPackage('localrec', version='1.2.0',
               tar='localrec-1.2.0.tgz')

env.addPackage('resmap', version='1.1.5s2',
               tar='resmap-1.1.5-s2.tgz',
               deps=['scipy'])

env.addPackage('spider', version='21.13',
               tar='spider-web-21.13.tgz',
               neededProgs=['csh'])

env.addPackage('motioncorr', version='2.1',
               tar='motioncorr_v2.1.tgz')

env.addPackage('motioncor2', version='16.03.16',
               tar='motioncor2_03162016.tgz')

env.addPackage('motioncor2', version='16.10.19',
               tar='motioncor2_10192016.tgz')

env.addPackage('motioncor2', version='17.01.30',
               tar='motioncor2_01302017.tgz')

env.addPackage('simple', version='2.1',
               tar='simple2.tgz')

env.addPackage('chimera', version='1.10.1',
               tar='chimera-1.10.1-linux_x86_64.tgz',
               targetDir='chimera-1.10.1',
               commands=[('./scipion_installer','bin/chimera')])

env.addPackage('dogpicker', version='0.2.1',
               tar='dogpicker-0.2.1.tgz')

env.addPackage('nma',
               tar='nma.tgz',
               commands=[('cd ElNemo; make; mv nma_* ..', 'nma_elnemo_pdbmat'),
                         ('cd NMA_cart; LDFLAGS=-L%s/lib make; mv nma_* ..' %
                          os.environ['SCIPION_SOFTWARE'], 'nma_diag_arpack')],
               deps=['arpack'])

env.addPackage('cryoem', version='1.0',
                tar='cryoem-1.0.tgz',
                pythonMod=True, default=False,
                numpyIncludes=True,
                deps=[numpy, scipy, matplotlib, cythongsl])

env.addPackage('gEMpicker', version='1.1',
               tar='gEMpicker_v1.1.tgz')

env.addPackage('Gctf', version='0.50',
               tar='Gctf_v0.50.tgz')

env.addPackage('Gctf', version='1.06',
               tar='Gctf_v1.06.tgz')

env.addPackage('Gautomatch', version='0.53',
               tar='Gautomatch_v0.53.tgz')

env.addPackage('mag_distortion', version='1.0.1',
               tar='mag_distortion-1.0.1.tgz')

env.addPackage('ethan', version='1.2',
               tar='ethan-1.2.tgz',
               commands=[('make', 'ethan')])

env.execute()
