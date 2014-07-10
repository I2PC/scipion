#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:     I. Foche Perez (ifoche@cnb.csic.es)
# *              J. Burguet Castell (jburguet@cnb.csic.es)
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
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************


# First import the environment (this comes from SConstruct)
Import('env')


#  ************************************************************************
#  *                                                                      *
#  *                              Libraries                               *
#  *                                                                      *
#  ************************************************************************

# Tcl/Tk
tcl = env.AddLibrary(
    'tcl',
    tar='tcl8.6.1-src.tgz',
    buildDir='tcl8.6.1/unix',
    libs=['lib/libtcl8.6.so'],
    flags=['--enable-threads'])

tk = env.AddLibrary(
    'tk',
    tar='tk8.6.1-src.tgz',
    buildDir='tk8.6.1/unix',
    libs=['lib/libtk8.6.so'],
    flags=['--enable-threads'],
    deps=[tcl])

# sqlite
sqlite = env.AddLibrary(
    'sqlite',
    tar='sqlite-3.6.23.tgz',
    libs=['lib/libsqlite3.so'],
    flags=['CPPFLAGS=-w',
           'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1'])

# zlib
zlib = env.AddLibrary(
    'zlib',
    tar='zlib-1.2.8.tgz',
    buildDir='zlib-1.2.8',
    libs=['lib/libz.so'],
    autoConfigTarget='configure.log')

# Python
python = env.AddLibrary(
    'python',
    tar='Python-2.7.8.tgz',
    libs=['lib/libpython2.7.a', 'bin/python'],
    flags=['--enable-shared',
           'CFLAGS=-I/usr/include/ncurses'],
    deps=[sqlite, tk])


#  ************************************************************************
#  *                                                                      *
#  *                           Python Modules                             *
#  *                                                                      *
#  ************************************************************************

# Helper function to include the python dependency automatically.
def addModule(*args, **kwargs):
    kwargs['deps'] = kwargs.get('deps', []) + python
    return env.AddModule(*args, **kwargs)


# numpy
addModule('numpy',
          tar='numpy-1.8.1.tgz')

# matplotlib
addModule('matplotlib',
          tar='matplotlib-1.3.1.tgz')

# psutil
addModule('psutil',
          tar='psutil-2.1.1.tgz')

# mpi4py
addModule('mpi4py',
          tar='mpi4py-1.3.1.tgz')

# scipy
addModule('scipy',
          tar='scipy-0.14.0.tgz',
          default=False)

# bibtex
addModule('bibtexparser',
          tar='bibtexparser-0.5.tgz')

# django
addModule('django',
          tar='Django-1.5.5.tgz')

# paramiko
addModule('paramiko',
          tar='paramiko-1.14.0.tgz',
          default=False)

# PIL
addModule('PIL',
          tar='Imaging-1.1.7.tgz',
          deps=[zlib])



# ##########################
# # EXECUTING COMPILATIONS #
# ##########################

# # EXTERNAL LIBRARIES

# t= env.AutoConfig(source=Dir('software/tmp/%s' % SCIPION['LIBS']['sqlite'][DIR]),
#                AutoConfigTarget='Makefile',
#                AutoConfigSource='configure',
#                AutoConfigParams=['CPPFLAGS=-w',
#                                  'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1',
#                                  '--prefix=%s' % abspath(SCIPION['FOLDERS'][SOFTWARE_FOLDER])])

# # env.AutoConfig(source=Dir('software/tmp/%s' % SCIPION['LIBS']['sqlite'][DIR]),
# #                target='makefile-sqlite',
# #                AutoConfigTarget='Makefile',
# #                AutoConfigSource='configure',
# #                AutoConfigParams=['CPPFLAGS=-w',
# #                                  'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1',
# #                                  '--prefix=%s' % abspath(SCIPION['FOLDERS'][SOFTWARE_FOLDER])])
# env['CROSS_BUILD'] = False
# env.Make(source=t,
#          target='software/lib/libsqlite3.so',
#          MakePath='software/tmp/%s' % SCIPION['LIBS']['sqlite'][DIR],
#          MakeEnv=os.environ,
#          MakeTargets='install',
#          MakeStdOut='aaaaa')

# # env.CompileLibrary('sqlite',
# #                    source=sqliteUntar,
# #                    flags=['CPPFLAGS=-w',
# #                           'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1',
# #                           '--prefix=%s' % join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)],
# #                    target='libsqlite3.so',
# #                    autoSource='configure') #SCIPION['LIBS']['sqlite'][DIR])

# # env.CompileLibrary('tcl',
# #                    source=tclUntar,
# #                    flags=['--enable-threads',
# #                           '--prefix=%s' % join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)],
# #                    target='libtcl.so',
# #                    autoSource=join('unix','Makefile.in'),
# #                    autoTarget=join('unix','Makefile'),
# #                    makePath='unix')

# # env.CompileLibrary('tk',
# #                    source=tkUntar,
# #                    flags=['--enable-threads',
# # #                          '--with-tcl="%s"' % join(Dir(join(SCIPION['FOLDERS'][SOFTWARE_FOLDER], 'lib64')).abspath),
# #                           '--prefix=%s' % join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)],
# #                    target='libtk.so',
# #                    autoSource=join('unix','Makefile.in'),
# #                    autoTarget=join('unix','Makefile'),
# #                    makePath='unix')

# # env.CompileLibrary('zlib',
# #                    source=zlibUntar,
# #                    autoSource='configure',
# #                    autoTarget='configure.log',
# #                    flags=['--prefix=%s' % join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath)],
# #                    target='libz.so')

# # # PYTHON

# # pythonMake = env.CompileLibrary('python',
# #                                 source=pythonUntar,
# #                                 flags=['--prefix=%s' % join(Dir(SCIPION['FOLDERS'][SOFTWARE_FOLDER]).abspath),
# #                                        'CFLAGS=-I/usr/include/ncurses'],
# #                                 target='libpython2.7.so',
# #                                 autoSource='Makefile.pre.in')

# # EMPackagesDeps = env.CompileWithSetupPy('python', deps=pythonMake)


# # # PYTHON MODULES

# # EMPackagesDeps += env.CompileWithSetupPy('numpy')
# # EMPackagesDeps += env.CompileWithSetupPy('matplotlib')
# # EMPackagesDeps += env.CompileWithSetupPy('psutil')
# # EMPackagesDeps += env.CompileWithSetupPy('mpi4py')
# # EMPackagesDeps += env.CompileWithSetupPy('scipy')
# # EMPackagesDeps += env.CompileWithSetupPy('bibtexparser')
# # EMPackagesDeps += env.CompileWithSetupPy('django')
# # EMPackagesDeps += env.CompileWithSetupPy('paramiko')
# # EMPackagesDeps += env.CompileWithSetupPy('pil')


# # # EM PACKAGES

# # # Xmipp3.1
# # env.AddPackage('xmipp',
# #                dft=False,
# #                tar='Xmipp-3.1-src.tgz',
# #                dir='xmipp',
# #                url='http://xmipp.cnb.csic.es/Downloads/Xmipp-3.1-src.tgz')

# # # Bsoft
# # env.AddPackage('bsoft',
# #                dft=False,
# #                tar='bsoft1_8_8_Fedora_12.tgz',
# #                dir='bsoft')

# # # CtfFind
# # env.AddPackage('ctffind',
# #                dft=False,
# #                tar='ctffind_V3.5.tgz',
# #                dir='ctf')

# # # EMAN2
# # env.AddPackage('eman2',
# #                dft=False,
# #                tar='eman2.1beta3.linux64.tgz',
# #                dir='EMAN2')

# # # frealign
# # env.AddPackage('frealign',
# #                dft=False,
# #                tar='frealign_v9.07.tgz')

# # # relion
# # env.AddPackage('relion',
# #                dft=False,
# #                tar='relion-1.2.tgz')
# # # spider
# # env.AddPackage('spider',
# #                dft=False,
# #                tar='spider-web-21.13.tgz',
# #                dir='spider-web')

# # if not GetOption('clean'):
# #     env.Download('xmipp', type='EMPackage')
# #     env.Download('bsoft', type='EMPackage')
# #     env.Download('ctffind', type='EMPackage')
# #     env.Download('eman2', type='EMPackage')
# #     env.Download('frealign', type='EMPackage')
# #     env.Download('relion', type='EMPackage')
# #     env.Download('spider', type='EMPackage')

# # # Purge option
# # if GetOption('purge'):
# #     env.RemoveInstallation()
