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
import glob
import os
import sys

import shutil

from install.funcs import Environment

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


SW_BIN = env.getBinFolder()
SW_LIB = env.getLibFolder()
SW_INC = env.getIncludeFolder()
SW_TMP = env.getTmpFolder()
SW_PYT_PACK = env.getPythonPackagesFolder()

#  *******************************
#  *  DETECT CURRENT INSTALLATION
#  *******************************
# Try to detect current installation and correct it if necessary


def clean_python_2_7_8_installation():
    # Detects installations where python 2.7.8 was installed.
    # In those installations we where using sqlite 3.6.23 and matplotlib-1.3.1
    # A bit of a hack but we will check based on matplotlib path!
    # Also this is not an exhaustive clean that might be more detailed
    # but enough to trigger the proper installation of the new versions.

    oldMatplotLibPath = Environment.getPythonPackagesFolder() + '/matplotlib-1.3.1*'
    print'Matplot lib path at %s' % oldMatplotLibPath

    def removeByPattern(pattern):
        for f in glob.glob(pattern):
            os.remove(f)

    # If old matplot lib exists
    if len(glob.glob(oldMatplotLibPath)) != 0:
        print "OLD Installation identified: removing Python and sqlite"

        # remove sqlite3 3.6.23
        sqliteLibs = Environment.getLibFolder() + "/libsqlite3*"
        removeByPattern(sqliteLibs)

        sqliteInc = Environment.getIncludeFolder() + "/sqlite3*"
        removeByPattern(sqliteInc)

        os.remove(Environment.getBinFolder() + "/sqlite3")

        # remove python installation
        pythonBinaries = Environment.getBinFolder() + "/python*"
        removeByPattern(pythonBinaries)

        # Python at include
        pythonIncludes = Environment.getIncludeFolder() + "/python2.7"
        shutil.rmtree(pythonIncludes)

        # Python at lib folder
        shutil.rmtree(Environment.getPythonFolder())

        return


clean_python_2_7_8_installation()

#  ************************************************************************
#  *                                                                      *
#  *                              Libraries                               *
#  *                                                                      *
#  ************************************************************************

cmake = env.addLibrary(
    'cmake',
    tar='cmake-3.2.2.tgz',
    targets=[env.getBin('cmake')],
    commands=[('cd ' + SW_TMP + '/cmake-3.2.2; '
               './bootstrap --prefix=../.. --parallel=%d' % env.getProcessors(),
               SW_TMP + '/cmake-3.2.2/Makefile'),
              ('cd ' + SW_TMP + '/cmake-3.2.2; make install -j %d'
               % env.getProcessors(), SW_BIN + '/cmake')],
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
                   targets=SW_BIN + '/wish',
                   cwd= SW_BIN)

zlib = env.addLibrary(
    'zlib',
    targets=[env.getLib('z')],
    tar='zlib-1.2.8.tgz',
    configTarget='zlib.pc')

jpeg = env.addLibrary(
    'jpeg',
    tar='libjpeg-turbo-1.3.1.tgz',
    flags=['--without-simd'])

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
    tar='SQLite-1a584e49.tgz',
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
    tar='Python-2.7.14.tgz',
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
    commands=[('cd ' + SW_TMP + '/sh_alignment; make install',
               SW_PYT_PACK + '/sh_alignment/frm.py')],
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
    commands=[('cd ' + SW_BIN + '; ln -s $(which gfortran) f77',
               SW_BIN + '/f77'),
              ('cd ' + SW_TMP + '/arpack-96; make all',
               SW_LIB +'libarpack.a')],
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
    commands=[('cp -rf ' + SW_TMP + '/boost_1_56_0/boost ' + SW_INC + '/',
               SW_INC + '/boost')],
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

# Add pip to our python
pip = env.addTarget('pip')
pip.addCommand('python scripts/get-pip.py', targets=SW_PYT_PACK + '/pip',
               default=True, final=True)

# Scons has a different pattern: it is expected to be in bin..TODO
scons = env.addModule(
     'scons',
     targets=[env.getBin('scons')],
     tar='scons-2.3.4.tgz')
#env.addPipModule('scons','2.3.6', target='scons-2.3.6')

# Required python modules
env.addPipModule('setuptools', '38.2.5')
numpy = env.addPipModule('numpy','1.14.1')
matplotlib = env.addPipModule('matplotlib', '1.5.3', target='matplotlib-1.5.3*')
env.addPipModule('psutil', '2.1.1', target='psutil-2.1.1*')
env.addPipModule('mpi4py', '1.3.1')
scipy = env.addPipModule('scipy', '0.14.0',
                     default=not noScipy, deps=[lapack, matplotlib])
env.addPipModule('bibtexparser', '0.6.2')
env.addPipModule('django', '1.5.5')
env.addPipModule('Pillow', '2.5.1', target='Pillow-2.5.1*',
    deps=[jpeg])


# Optional python modules
env.addPipModule('paramiko', '1.14.0', default=False)
# 1.4.8 could not be found ! Using latest available
env.addPipModule('winpdb', '1.3.6', default=False)

env.addPipModule('lxml', '3.4.1', target='lxml-3.4.1*', default=False)
env.addPipModule('requests', '2.18.4', default=False)

# These were dependencies of iPython
env.addPipModule('pyzmq', '2.2.0.1', target='pyzmq*', default=False)
env.addPipModule('jinja2', '2.7.3', default=False)
env.addPipModule('tornado', '4.0.2', default=False)
env.addPipModule('ipython', '2.1.0', target='IPython', default=False)
cython = env.addPipModule('cython', '0.22', target='Cython-0.22*', default=False)
cythongsl = env.addPipModule('cythongsl','0.2.1',
                         target='CythonGSL-0.2.1*',
                         default=False, deps=[cython])
env.addPipModule('scikit-learn', '0.17', target='scikit_learn*',
             default=False, deps=[scipy, cython])


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

env.addPackage('relion', version='2.1',
              tar='relion-2.1.tgz',
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

env.addPackage('motioncor2', version='17.01.30',
               tar='motioncor2_01302017.tgz')

env.addPackage('motioncor2', version='1.0.0',
               tar='motioncor2_1.0.0.tgz')

env.addPackage('motioncor2', version='1.0.2',
               tar='motioncor2-1.0.2.tgz')

env.addPackage('motioncor2', version='1.0.4',
               tar='motioncor2-1.0.4.tgz')

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
                         ('cd NMA_cart; LDFLAGS=-L%s make; mv nma_* ..' %
                          Environment.getLibFolder(), 'nma_diag_arpack')],
               deps=['arpack'])

env.addPackage('cryoem', version='1.0',
                tar='cryoem-1.0.tgz',
                pythonMod=True, default=False,
                deps=[scipy, matplotlib, cythongsl])

env.addPackage('powerfit', version='2.0',
                tar='powerfit.tgz',
                targets=['powerfit-2.0*'],
                pythonMod=True, default=False,
                deps=[numpy, scipy, fftw3])

env.addPackage('gEMpicker', version='1.1',
               tar='gEMpicker_v1.1.tgz')

env.addPackage('gctf', version='0.50',
               tar='Gctf_v0.50.tgz')

env.addPackage('gctf', version='1.06',
               tar='Gctf_v1.06.tgz')

env.addPackage('gautomatch', version='0.53',
               tar='Gautomatch_v0.53.tgz')

env.addPackage('mag_distortion', version='1.0.1',
               tar='mag_distortion-1.0.1.tgz')

env.addPackage('ethan', version='1.2',
               tar='ethan-1.2.tgz',
               commands=[('make', 'ethan')])

fsc_commands = [('conda env create -f environment.yml && touch IS_INSTALLED',
                 'IS_INSTALLED')]

env.addPackage('nysbc-3DFSC', version='2.5',
               tar='nysbc-3DFSC_2.5.tgz',
               commands=fsc_commands,
               neededProgs=['conda'])

env.addPackage('cryoEF', version='1.1.0',
               tar='cryoEF_v1.1.0.tgz')

env.execute()
