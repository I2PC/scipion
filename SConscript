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


#  ************************************************************************
#  *                                                                      *
#  *                              Libraries                               *
#  *                                                                      *
#  ************************************************************************

# We might want to add freetype and make tcl depend on it. That would be:
# freetype = env.AddLibrary(
#     'freetype',
#     tar='freetype-2.5.3.tgz',
#     autoConfigTarget='config.mk')
# But because freetype's compilation is a pain, it's better to use whatever
# version is in the system.

fftw = env.AddLibrary(
    'fftw',
    tar='fftw-3.3.4.tgz',
    targets=['lib/libfftw3.so'],
    flags=['--enable-threads', '--enable-shared'],
    default=False)

tcl = env.AddLibrary(
    'tcl',
    tar='tcl8.6.1-src.tgz',
    buildDir='tcl8.6.1/unix',
    targets=['lib/libtcl8.6.so'],
    flags=['--enable-threads'],
    clean=['software/tmp/tcl8.6.1'])

tk = env.AddLibrary(
    'tk',
    tar='tk8.6.1-src.tgz',
    buildDir='tk8.6.1/unix',
    targets=['lib/libtk8.6.so'],
    libChecks=['xft'],
    flags=['--enable-threads'],
    deps=[tcl],
    clean=['software/tmp/tk8.6.1'])

zlib = env.AddLibrary(
    'zlib',
    tar='zlib-1.2.8.tgz',
    targets=['lib/libz.so'],
    addPath=False,
    autoConfigTarget='zlib.pc')

jpeg = env.AddLibrary(
    'jpeg',
    tar='libjpeg-turbo-1.3.1.tgz',
    flags=([] if env.ProgInPath('nasm') else ['--without-simd']))

sqlite = env.AddLibrary(
    'sqlite',
    tar='sqlite-3.6.23.tgz',
    targets=['lib/libsqlite3.so'],
    flags=['CPPFLAGS=-w',
           'CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1'])

python = env.AddLibrary(
    'python',
    tar='Python-2.7.8.tgz',
    targets=['lib/libpython2.7.so', 'bin/python'],
    flags=['--enable-shared'],
    deps=[sqlite, tk, zlib])

env.AddLibrary(
    'parallel',
    tar='parallel-20140922.tgz',
    targets=['bin/parallel'],
    deps=[zlib])

boost_headers_only = env.ManualInstall(
    'boost_headers_only',
    tar='boost_1_56_0.tgz',
    extraActions=[
        ('%s/software/include/boost' % env['SCIPION_HOME'],
         'cp -rf boost %s/software/include' % env['SCIPION_HOME'])],
    default=False)

swig = env.AddLibrary(
    'swig',
    tar='swig-3.0.2.tgz',
    flags=['--without-pcre'],
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
    deps=[numpy, dateutil, pyparsing])

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
    deps=[numpy, matplotlib])

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
               extraActions=[('xmipp.bashrc',
                             './install.sh --unattended=true --gui=false -j %s'
                              % GetOption('num_jobs'))],
               default=False)

env.AddPackage('bsoft',
               tar='bsoft1_8_8_Fedora_12.tgz',
               default=False)

env.AddPackage('ctffind',
               tar='ctffind_V3.5.tgz',
               default=False)

env.AddPackage('eman',
               tar='eman2.1beta3.linux64.tgz',
               extraActions=[('eman2.bashrc', './eman2-installer')],
               default=False)

env.AddPackage('frealign',
               tar='frealign_v9.07.tgz',
               default=False)

env.AddPackage('pytom',
               tar='pytom_develop0.962.tgz',
               extraActions=[('pytomc/libs/libtomc/libs/libtomc.so',
                             'MPILIBDIR=%s MPIINCLUDEDIR=%s SCIPION_HOME=%s ./scipion_installer'
                              % (env['MPI_LIBDIR'],env['MPI_INCLUDE'],env['SCIPION_HOME']))],
               deps=[boost_headers_only, fftw, swig],
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

# This last one already contains the binary.


# TODO: check if we have to use the "purge" option below:

# # Purge option
# if GetOption('purge'):
#     env.RemoveInstallation()
