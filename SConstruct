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

# Builders and pseudobuilders used be SConscript to install things.


import os
from os.path import join, abspath, splitext
from glob import glob
import tarfile
import fnmatch
import platform
import subprocess
import SCons.Script
import SCons.SConf


# OS boolean vars
MACOSX = (platform.system() == 'Darwin')
WINDOWS = (platform.system() == 'Windows')
LINUX = (platform.system() == 'Linux')

# URL where we have most of our tgz files for libraries, modules and packages.
URL_BASE = 'http://scipionwiki.cnb.csic.es/files/scipion/software'

# Define our builders.
Download = Builder(action='wget -nv $SOURCE -c -O $TARGET')
Untar = Builder(action='tar -C $cdir --recursive-unlink -xzf $SOURCE')

# Create the environment the whole build will use.
env = Environment(ENV=os.environ,
                  BUILDERS=Environment()['BUILDERS'],
                  tools=['Make', 'AutoConfig'],
                  toolpath=[join('software', 'install', 'scons-tools')])
#Export('env')
# TODO: BUILDERS var added from the tricky creation of a new environment.
# If not, they lose default builders like "Program", which are needed later
# (by CheckLib and so on). See http://www.scons.org/doc/2.0.1/HTML/scons-user/x3516.html
# See how to change it into a cleaner way (not doing BUILDERS=Environment()['BUILDERS']!)

# Message from autoconf and make, so we don't see all its verbosity.
env['AUTOCONFIGCOMSTR'] = "Configuring $TARGET from $SOURCES"
env['MAKECOMSTR'] = "Compiling & installing $TARGET from $SOURCES "


def progInPath(env, prog):
    "Is program prog in PATH?"
    return any(os.path.exists('%s/%s' % (base, prog)) for base in
               os.environ.get('PATH', '').split(os.pathsep))


# This functions was previously used for checking the presence of a
# library, with the additional functionality of letting the user decide to
# continue or not. Now we're using SCons stuff (its CheckLib function). We
# may want to delete this function in the future
def checkConfigLib(target, source, env):
    "See if we have library <name> from software/log/lib_<name>.log"

    # This is used to create the CheckConfigLib builder.

    # source is ignored and must be ''
    # target must look like 'software/log/lib_<name>.log'

    for tg in map(str, target):  # "target" is a list of SCons.Node.FS.File
        name = tg[len('software/log/lib_'):-len('.log')]
        try:
            subprocess.check_call(['pkg-config', '--cflags', '--libs', name],
                                  stdout=open(os.devnull, 'w'),
                                  stderr=subprocess.STDOUT)
        except (subprocess.CalledProcessError, OSError) as e:
            try:
                subprocess.check_call(['%s-config' % name, '--cflags'])
            except (subprocess.CalledProcessError, OSError) as e:
                print """
  ************************************************************************
    Warning: %s not found. Please consider installing it first.
  ************************************************************************

Continue anyway? (y/n)""" % name
                if raw_input().upper() != 'Y':
                    Exit(2)

                # If we continue without the lib, we write the target log
                # so scons knows not to ask again.
                open(tg, 'w').write('%s off\n' % name)
        # If we do have it, write the tg so we do not ask again either.
        open(tg, 'w').write('%s on\n' % name)


CheckConfigLib = Builder(action=checkConfigLib)


# Add the path to dynamic libraries so the linker can find them.
if LINUX:
    env.AppendUnique(LIBPATH=os.environ.get('LD_LIBRARY_PATH', ''))
elif MACOSX:
    print "OS not tested yet"
    env.AppendUnique(LIBPATH=os.environ.get('DYLD_FALLBACK_LIBRARY_PATH', ''))
elif WINDOWS:
    print "OS not tested yet"

# Python and SCons versions are fixed
env.EnsurePythonVersion(2,7)
env.EnsureSConsVersion(2,3,2)
# TODO: see after all is clean and crispy if we can avoid fixing the versions.
# We can specify a range of valid version after we check it works with them.


#  ************************************************************************
#  *                                                                      *
#  *                       Auxiliar functions                             *
#  *                                                                      *
#  ************************************************************************


def appendUnique(elist, element):
    'Add element to a list only if it doesnt previously exist'
    if element not in elist:
        if not isinstance(element, basestring):
            elist.extend(element)
        else:
            elist.append(element)


def CheckMPI(context, mpi_inc, mpi_libpath, mpi_lib, mpi_cc, mpi_cxx, mpi_link, replace):
    context.Message('* Checking for MPI ... ')

    lastLIBS = context.env.get('LIBS', '')
    lastLIBPATH = context.env['LIBPATH']
    lastCPPPATH = context.env['CPPPATH']
    lastCC = context.env['CC']
    lastCXX = context.env['CXX']

    # TODO Replace() also here?
    context.env.Append(LIBS=mpi_lib, LIBPATH=mpi_libpath,
                       CPPPATH=mpi_inc)
    context.env.Replace(LINK=mpi_link)
    context.env.Replace(CC=mpi_cc, CXX=mpi_cxx)

    # Test only C++ mpi compiler
    ret = context.TryLink('''
    #include <mpi.h>
    int main(int argc, char** argv)
    {
        MPI_Init(0, 0);
        MPI_Finalize();
        return 0;
    }
    ''', '.cpp')

    # NOTE: We don't want MPI flags for not-mpi programs (always revert)
    # env['mpi'] remains 1 so those can be enabled again when needed
    if not replace:
        context.env.Replace(LIBS=lastLIBS)
        context.env.Replace(LIBPATH=lastLIBPATH)
        context.env.Replace(CPPPATH=lastCPPPATH)
        context.env.Replace(CC=lastCC)
        context.env.Replace(CXX=lastCXX)

    context.Result(ret)
    return ret


#  ************************************************************************
#  *                                                                      *
#  *                           Pseudobuilders                             *
#  *                                                                      *
#  ************************************************************************

# We have 6 "Pseudo-Builders" http://www.scons.org/doc/HTML/scons-user/ch20.html
#
# They are:
#   addLibrary        - install a library
#   addModule         - install a Python module
#   addPackage        - install an EM package
#   addPackageLibrary - install a EM package library
#   addJavaLibrary    - install a java jar
#   addProgram        - install a EM package program
#   manualInstall     - install by manually running commands
#
# Their structure is similar:
#   * Define reasonable defaults
#   * Set --with-<name> option as appropriate
#   * Concatenate builders
#
# For the last step we concatenate the builders this way:
#   target1 = Builder1(env, target1, source1)
#   SideEffect('dummy', target1)
#   target2 = Builder2(env, target2, source2=target1)
#   SideEffect('dummy', target2)
#   ...
#
# So each target becomes the source of the next builder, and the
# dependency is solved. Also, we use SideEffect('dummy', ...) to
# ensure it works in a parallel build (see
# http://www.scons.org/wiki/SConsMethods/SideEffect), and does not try
# to do one step while the previous one is still running in the background.


def addLibrary(env, name, tar=None, buildDir=None, configDir=None, 
               targets=None, libChecks=[], url=None, flags=[], addPath=True,
               autoConfigTargets='Makefile', deps=[], clean=[], default=True):
    """Add library <name> to the construction process.

    This pseudobuilder checks that the needed programs are in PATH,
    downloads the given url, untars the resulting tar file, configures
    the library with the given flags, compiles it (in the given
    buildDir) and installs it. It also tells SCons about the proper
    dependencies (deps).

    If addPath=False, we will not pass the variables PATH and
    LD_LIBRARY_PATH pointing to our local installation directory.

    If default=False, the library will not be built unless the option
    --with-<name> is used.

    Returns the final targets, the ones that Make will create.

    """
    # Use reasonable defaults.
    tar = tar or ('%s.tgz' % name)
    url = url or ('%s/external/%s' % (URL_BASE, tar))
    buildDir = buildDir or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]
    configDir = configDir or buildDir
    targets = targets or [File('#software/lib/lib%s.so' % name).abspath]
    if isinstance(buildDir, basestring):
        buildDir = [buildDir]
        configDir = [configDir]
        flags = [flags]
        autoConfigTargets = [autoConfigTargets]
        targets = [targets]
    if len(buildDir) != len(configDir) != len(flags):
        print >> sys.stderr, 'ERROR: buildDir, configDir and flags length must be equal. Exiting...'
        Exit(1)

    # Add "software/lib" and "software/bin" to LD_LIBRARY_PATH and PATH.
    def pathAppend(var, value):
        valueOld = os.environ.get(var, '')
        i = 0  # so if flags is empty, we put it at the beginning too
        for flag in flags:
            for i in range(len(flag)):
                if flag[i].startswith('%s=' % var):
                    valueOld = flag.pop(i).split('=', 1)[1] + ':' + valueOld
                    break
            flag.insert(i, '%s=%s:%s' % (var, value, valueOld))

    if addPath:
        pathAppend('LD_LIBRARY_PATH', Dir('#software/lib').abspath)
        pathAppend('PATH', Dir('#software/bin').abspath)

    # Install everything in the appropriate place.
    for flag in flags:
        flag += ['--prefix=%s' % Dir('#software').abspath,
                 '--libdir=%s' % Dir('#software/lib').abspath]  # not lib64

    # Add the option --with-name, so the user can call SCons with this
    # to activate the library even if it is not on by default.
    if not default:
        AddOption('--with-%s' % name, dest=name, action='store_true',
                  help='Activate library %s' % name)

    # Add all the checks
    for lib in libChecks:
        libraryTest(env, lib)
    # Alternatively we could use CheckConfigLib() (as in commit cac431d117)

    # Create and concatenate the builders.
    tDownload = Download(env, File('#software/tmp/%s' % tar).abspath, Value(url))
    SideEffect('dummy', tDownload)  # so it works fine in parallel builds
    tUntar = Untar(env, File('#software/tmp/%s/configure' % configDir[0]).abspath, tDownload,
                   cdir=Dir('#software/tmp').abspath)
    SideEffect('dummy', tUntar)  # so it works fine in parallel builds
    Clean(tUntar, Dir('#software/tmp/%s' % buildDir).abspath)
    tConfig = tMake = []
    toReturn = []
    for x, folder in enumerate(buildDir):
        tConfig.append(*env.AutoConfig(
            source=Dir('#software/tmp/%s' % configDir[x]),
            AutoConfigTarget=autoConfigTargets[x],
            AutoConfigSource='configure',
            AutoConfigParams=flags[x],
            AutoConfigStdOut=File('#software/log/%s_config_%s.log' % (name, x)).abspath))
        SideEffect('dummy', tConfig[x])  # so it works fine in parallel builds
        env.Depends(tConfig[x], tUntar)
        make = env.Make(
            source=tConfig[x],
            target=targets[x],
            MakePath=Dir('#software/tmp/%s' % buildDir[x]).abspath,
            MakeEnv=os.environ,
            MakeTargets='install',
            MakeStdOut=File('#software/log/%s_make_%s.log' % (name, x)).abspath)
        if not isinstance(make, basestring):
            tMake+=make
        else:
            tMake.append(make)
        tMake.append(make)
        SideEffect('dummy', tMake[x]) # so it works fine in parallel builds
        for cFile in clean:
            Clean(tMake[x], cFile)
        # Clean the special generated files
        # Add the dependencies.
        for dep in deps:
            env.Depends(tConfig[x], dep)
        # Make library by default
        if default or GetOption(name):
            env.Default(tMake[x])
        toReturn += [tMake[x]]
    return toReturn


def addPackageLibrary(env, name, dirs=None, tars=None, untarTargets=None, patterns=None, incs=None, libs=None, prefix=None, suffix=None, installDir=None, libpath=['lib'], deps=[], mpi=False, cuda=False, default=True):
    """Add self-made and compiled shared library to the compilation process
    
    This pseudobuilder access given directory, compiles it
    and installs it. It also tells SCons about it dependencies.

    If default=False, the library will not be built unless the option
    --with-<name> is used.

    Returns the final targets, the ones that Make will create.
    """
    libs = libs or []
    dirs = dirs or []
    tars = tars or []
    untarTargets = untarTargets or ['configure']
    lastTarget = deps
    incs = incs or []
    prefix = 'lib' if prefix is None else prefix
    suffix = '.so' if suffix is None else suffix
    
    basedir = 'lib'
    fullname = prefix + name
    patterns = patterns or []
    sources = []
    untars = []
    libpath.append(Dir('#software/lib').abspath)
    
    if tars:
        for x, dir in enumerate(dirs):
            tarDestination, tarName = splitext(tars[x])
            tUntar = Untar(env, File(join(dirs[x], untarTargets[x])).abspath, 
                           File(tars[x]),
                           cdir=Dir(dirs[x]).abspath)
            SideEffect('dummy', tUntar)
            env.Depends(tUntar, deps)
            untars += tUntar
        lastTarget = untars

    for x, p in enumerate(patterns):
        sourcesRel = join(dirs[x], p)
        if not sourcesRel.startswith(basedir):
            sourcesRel = join(basedir, p)
#        if p.endswith(".cu"):
#            cudaFiles = True
        # select sources inside tarfile but we only get those files that match the pattern
        if tars: 
            tarfiles = tarfile.open(tars[x]).getmembers() 
            sources += [join(dirs[x], tarred.name) for tarred in tarfiles if fnmatch.fnmatch(tarred.name, p)]
        else:
            sources += glob(join(dirs[x], patterns[x]))

    mpiArgs = {}
    if mpi:
        libpath.append(env['MPI_LIBDIR'])
        libs.append(env['MPI_LIB'])
        mpiArgs = {'CC': env['MPI_CC'],
                   'CXX': env['MPI_CXX']}
        conf = Configure(env, custom_tests = {'CheckMPI': CheckMPI})
        if not conf.CheckMPI(env['MPI_INCLUDE'], env['MPI_LIBDIR'], env['MPI_LIB'], env['MPI_CC'], env['MPI_CXX'], env['MPI_LINKERFORPROGRAMS'], False):
            print >> sys.stderr, 'ERROR: MPI is not properly working. Exiting...'
            Exit(1)
        env = conf.Finish()
        env.PrependENVPath('PATH', env['MPI_BINDIR'])
    
    # FIXME: There must be a key in env dictionary that breaks the compilation. Please find it to make it more beautiful
    env2 = Environment()
    env2['ENV']['PATH'] = env['ENV']['PATH']

    
    #print env2.Dump()
    library = env2.SharedLibrary(
              target=join(basedir, fullname),
              source=sources,
              CPPPATH=incs + [env['CPPPATH']] + [Dir('#software/include').abspath],
              LIBPATH=libpath,
              LIBS=libs,
              SHLIBPREFIX=prefix,
              SHLIBSUFFIX=suffix,
              **mpiArgs
              )
    SideEffect('dummy', library)

    for previous in lastTarget:
        env.Depends(library, previous)
    
    if installDir:
        install = env.Install(installDir, library)
        SideEffect('dummy', install)
        lastTarget = install
    else:
        lastTarget = library
    env.Default(lastTarget)
    
    env.Depends(lastTarget, deps)

    return lastTarget


def symLink(env, target, source):
    #As the link will be in bin/ directory we need to move up
    sources = source
    if isinstance(target, list):
        link = str(target[0])
        sources = str(source[0])
    else:
        link = target
    sources = os.path.relpath(sources, os.path.split(link)[0])
    if os.path.lexists(link):
        os.remove(link)
    print 'Linking to %s from %s' % (sources, link)
    os.symlink(sources, link)
    return target



def addJavaLibrary(env, name, jar=None, dirs=None, patterns=None, installDir=None, buildDir=None, classDir=None, sourcePath=None, deps=[], default=True):
    """Add self-made and compiled java library to the compilation process
    
    This pseudobuilder access given directory, compiles it
    and installs it. It also tells SCons about it dependencies.

    If default=False, the library will not be built unless the option
    --with-<name> is used.

    Returns the final targets, the ones that Make will create.
    """
    jar = jar or '%s.jar' % name
    dirs = dirs or ['java/src']
    patterns = patterns or 'unset'
    buildDir = buildDir or 'java/build'
    installDir = installDir or 'java/lib'
    classDir = classDir or 'java/build'
    sourcePath = sourcePath or 'java/src'
    sources = []
    if patterns == 'unset':
        patterns = []
        for x in range(len(dirs)):
            patterns += ['*.java']

    
    if not default:
        return None
    
    deps = [join(installDir, '%s.jar' % name) for name in deps]
    for x, source in enumerate(dirs):
        sources += glob(join(dirs[x], patterns[x]))
    
    env2 = Environment()
    env2['ENV']['PATH'] = env['ENV']['PATH']
    env2['ENV']['JAVA_ROOT'] = env['ENV']['JAVA_ROOT']
    env2['ENV']['JAVA_HOME'] = env['ENV']['JAVA_HOME']
    env2['ENV']['JAVA_HOME'] = env['ENV']['JAVA_BINDIR']
    env2['ENV']['JRE_HOME'] = env['ENV']['JRE_HOME']
    env2.AppendUnique(JAVACLASSPATH=":".join(glob(join(Dir(installDir).abspath,'*.jar'))))
    env2.AppendUnique(JAVASOURCEPATH=Dir(sourcePath).abspath)

    jarCreation = env2.Jar(target=join(buildDir, jar), source=sources)
    SideEffect('dummy', jarCreation)
    lastTarget = jarCreation
    env.Default(jarCreation)
    for dep in deps:
        env.Depends(jarCreation, File(dep).abspath)
    
    install = env.Install(installDir, jarCreation)
    SideEffect('dummy', install)
    lastTarget = install

    env.Default(lastTarget)
    
    return lastTarget
    


def addJavaTest(env, name, source, installDir=None, default=True):
    """Add java test to the compilation process
    
    This pseudobuilder executes a java test using the Command builder.

    If default=False, the test will not be done unless the option
    --with-java-tests is used.

    Returns the final targets, the ones that Command returns.
    """
    if not default and not GetOption('run_java_tests'):
        return ''
    installDir = installDir or 'java/lib'
    classPath = ":".join(glob(join(Dir(installDir).abspath,'*.jar')))
    cmd = '%s -cp %s org.junit.runner.JUnitCore xmipp.test.%s' % (join(env['JAVA_BINDIR'], 'java'), classPath, name)
    runTest = env.Command(name, join(installDir, source), cmd)
    env.Alias('run_java_tests', runTest)
    Default(runTest)


def addProgram(env, name, src=None, pattern=None, installDir=None, libPaths=[], incs=[], libs=[], cxxflags=[], linkflags=[], deps=[], mpi=False, cuda=False, default=True):
    """Add, compile and install a program to the compilation process
    
    This pseudobuilder compiles a C++ program using CXX compiler and linker.
    
    This is designed to compile the different parts of a EM software

    If default=False, the program will not be compiled unless the 
    --with-<program-name> is used.

    Returns the final targets, the ones that Command returns.
    
    """

    if not default and not GetOption(name):
        AddOption('--with-%s' % name, dest=name, action='store_true',
                  help='Add the program %s to the compilation' % name)
        return ''
    src = src or ['src']
    pattern = pattern or ['*.cpp']
    installDir = installDir or 'bin'
    libs = libs or []
    libPathsCopy = libPaths + ['lib', '#software/lib']
    incs = incs or []
    incs += ['libraries', '#software/include', '#software/include/python2.7']
    if cuda:
        libs += ['cudart', 'cublas', 'cufft', 'curand', 'cusparse', 'npp', 'nvToolsExt', 'opencv_gpu']
        incs += [join(env['CUDA_SDK_PATH'], "CUDALibraries","common","inc"),
                 join(env['CUDA_SDK_PATH'], "shared","inc"),
                 join(env['CUDA_SDK_PATH'],"CUDALibraries","common","lib","linux"),
                 join("/usr","local","cuda","lib64"),
                 env['CUDA_LIB_PATH']]
    sources = []
    for x, dir in enumerate(src):
        sources += glob(join(dir, pattern[x]))
    
    ccCopy = env['MPI_CC'] if mpi else env['CC']
    cxxCopy = env['MPI_CXX'] if mpi else env['CXX']
    linkCopy = env['MPI_LINKERFORPROGRAMS'] if mpi else env['LINKERFORPROGRAMS']
    incsCopy = incs + env['CPPPATH']
    libsCopy = libs
    cxxflagsCopy = cxxflags + env['CXXFLAGS']
    linkflagsCopy = linkflags + env['LINKFLAGS']
    ldLibraryPathCopy = [env['LIBPATH']]
    appendUnique(libPathsCopy, env.get('LIBPATH', ''))
    env2 = Environment()
    if mpi: 
        appendUnique(incsCopy, env['MPI_INCLUDE'])
        appendUnique(libPathsCopy, env['MPI_LIBDIR'])
        appendUnique(libsCopy, env['MPI_LIB'])
        appendUnique(ldLibraryPathCopy, env['MPI_LIBDIR'])
	env2['ENV']['LD_LIBRARY_PATH'] = env['ENV'].get('LD_LIBRARY_PATH', '')
	env2['ENV']['PATH'] = env['ENV']['PATH']


    program = env2.Program(
                          File(join(installDir, name)).abspath,
                          source=sources,
                          CC=ccCopy,
                          CXX=cxxCopy,
                          CPPPATH=incsCopy,
                          LIBPATH=libPathsCopy,
                          LIBS=libsCopy,
                          CXXFLAGS=cxxflagsCopy,
                          LINKFLAGS=linkflagsCopy,
                          LINK=linkCopy,
                          LD_LIBRARY_PATH=ldLibraryPathCopy
                          )
    env2.Default(program)
    
    if deps: 
        env2.Depends(program, [dep for dep in deps])
    
    return program


def addModule(env, name, tar=None, buildDir=None, targets=None,
              url=None, flags=[], deps=[], clean=[], default=True):
    """Add Python module <name> to the construction process.

    This pseudobuilder downloads the given url, untars the resulting
    tar file, configures the module with the given flags, compiles it
    (in the given buildDir) and installs it. It also tells SCons about
    the proper dependencies (deps).

    If default=False, the module will not be built unless the option
    --with-<name> is used.

    Returns the final target (software/lib/python2.7/site-packages/<name>).

    """
    # Use reasonable defaults.
    tar = tar or ('%s.tgz' % name)
    url = url or ('%s/python/%s' % (URL_BASE, tar))
    buildDir = buildDir or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]
    targets = targets or [name]
    flags += ['--prefix=%s' % Dir('#software').abspath]

    # Add the option --with-name, so the user can call SCons with this
    # to activate the module even if it is not on by default.
    if not default:
        AddOption('--with-%s' % name, dest=name, action='store_true',
                  help='Activate module %s' % name)

    # Create and concatenate the builders.
    tDownload = Download(env, 'software/tmp/%s' % tar, Value(url))
    SideEffect('dummy', tDownload)  # so it works fine in parallel builds
    tUntar = Untar(env, 'software/tmp/%s/setup.py' % buildDir, tDownload,
                   cdir='software/tmp')
    SideEffect('dummy', tUntar)  # so it works fine in parallel builds
    Clean(tUntar, 'software/tmp/%s' % buildDir)
    tInstall = env.Command(
        ['software/lib/python2.7/site-packages/%s' % t for t in targets],
        tUntar,
        Action('PYTHONHOME="%(root)s" LD_LIBRARY_PATH="%(root)s/lib" '
               'PATH="%(root)s/bin:%(PATH)s" '
               '%(root)s/bin/python setup.py install %(flags)s > '
               '%(root)s/log/%(name)s.log 2>&1' % {'root': Dir('#software').abspath,
                                                   'PATH': os.environ['PATH'],
                                                   'flags': ' '.join(flags),
                                                   'name': name},
               'Compiling & installing %s > software/log/%s.log' % (name, name),
               chdir=Dir('#software/tmp/%s' % buildDir).abspath))
    SideEffect('dummy', tInstall)  # so it works fine in parallel builds

    # Clean the special generated files
    for cFile in clean:
        Clean(lastTarget, cFile)

    # Add the dependencies.
    for dep in deps:
        env.Depends(tInstall, dep)

    if default or GetOption(name):
        env.Default(tInstall)

    return tInstall


def compilerConfig(env):
    """Check the good state of the C and C++ compilers and return the proper env."""

    conf = Configure(env)
    # ---- check for environment variables
    if 'CC' in os.environ:
        conf.env.Replace(CC=os.environ['CC'])
    else:
        conf.env.Replace(CC='gcc')
    print(">> Using C compiler: " + conf.env.get('CC'))

    if 'CFLAGS' in os.environ:
        conf.env.Replace(CFLAGS=os.environ['CFLAGS'])
        print(">> Using custom C build flags")

    if 'CXX' in os.environ:
        conf.env.Replace(CXX=os.environ['CXX'])
    else:
        conf.env.Replace(CXX='g++')
    print(">> Using C++ compiler: " + conf.env.get('CXX'))

    if 'CXXFLAGS' in os.environ:
        conf.env.Append(CPPFLAGS=os.environ['CXXFLAGS'])
        print(">> Appending custom C++ build flags : " + os.environ['CXXFLAGS'])

    if 'LDFLAGS' in os.environ:
        conf.env.Append(LINKFLAGS=os.environ['LDFLAGS'])
        print(">> Appending custom link flags : " + os.environ['LDFLAGS'])

    conf.CheckCC()
    conf.CheckCXX()
    env = conf.Finish()
    return env


def libraryTest(env, name, lang='c'):
    """Check the existence of a concrete C/C++ library."""
    env2 = Environment(LIBS=env.get('LIBS',''))
    conf = Configure(env2)
    conf.CheckLib(name, language=lang)
    env2 = conf.Finish()
    # conf.Finish() returns the environment it used, and we may want to use it,
    # like:  return conf.Finish()  but we don't do that so we keep our env clean :)



def addPackage(env, name, tar=None, buildDir=None, url=None, neededProgs=[],
               extraActions=[], deps=[], clean=[], reqs=[], default=True):
    """Add external (EM) package <name> to the construction process.

    This pseudobuilder downloads the given url, untars the resulting
    tar file and copies its content from buildDir into the
    installation directory <name>. It also tells SCons about the
    proper dependencies (deps).

    extraActions is a list of (target, command) that should be
    executed after the package is properly installed.

    If default=False, the package will not be built unless the option
    --with-<name> or --with-all-packages is used.

    Returns the final target (software/em/<name>).

    """
    # Use reasonable defaults.
    tar = tar or ('%s.tgz' % name)
    url = url or ('%s/em/%s' % (URL_BASE, tar))
    buildDir = buildDir or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]
    confPath = File('#software/cfg/%s.cfg' % name).abspath
    

    # Minimum requirements must be accomplished. To check them, we use
    # the req list, iterating with SConf CheckLib on it
    for req in reqs:
        libraryTest(env, req, reqs[req])

    # Add the option --with-<name>, so the user can call SCons with this
    # to get the package even if it is not on by default.
    AddOption('--with-%s' % name, dest=name, metavar='%s_HOME' % name.upper(),
              nargs='?', const='unset',
              help=("Get package %s. With no argument, download and "
                    "install it. To use an existing installation, pass "
                    "the package's directory." % name))
    # So GetOption(name) will be...
    #   None      if we did *not* pass --with-<name>
    #   'unset'   if we passed --with-<name> (nargs=0)
    #   PKG_HOME  if we passed --with-<name>=PKG_HOME (nargs=1)

    # See if we have used the --with-<package> option and exit if appropriate.
    if GetOption('withAllPackages'):
        defaultPackageHome = 'unset'
        # we asked to install all packages, so it is at least as if we
        # also did --with-<name>
    else:
        defaultPackageHome = None
        # by default it is as if we did not use --with-<name>

    packageHome = GetOption(name) or defaultPackageHome
    
    if not (default or packageHome):
        return ''
    
    lastTarget = tLink = None
    # If we do have a local installation, link to it and exit.
    if packageHome != 'unset':  # default value when calling only --with-package
        #FIXME: only for operating while programming
        if not os.path.exists(packageHome):
            # If it's a completed local installation. Just link to it and do nothing more.
            return env.Command(
                Dir('software/em/%s/bin' % name),
                Dir(packageHome),
                Action('rm -rf %s && ln -v -s %s %s' % (name, packageHome, name),
                       'Linking package %s to software/em/%s' % (name, name),
                       chdir='software/em'))

    # Check that all needed programs are there.
    for p in neededProgs:
        if not progInPath(env, p):
            print """
  ************************************************************************
    Warning: Cannot find program "%s" needed by %s
  ************************************************************************

Continue anyway? (y/n)""" % (p, name)
            if raw_input().upper() != 'Y':
                Exit(2)

    # Donload, untar, link to it and execute any extra actions.
    tDownload = Download(env, 'software/tmp/%s' % tar, Value(url))
    SideEffect('dummy', tDownload)  # so it works fine in parallel builds
    tUntar = Untar(env, Dir('software/em/%s/bin' % buildDir), tDownload,
                   cdir='software/em')
    SideEffect('dummy', tUntar)  # so it works fine in parallel builds
    Clean(tUntar, 'software/em/%s' % buildDir)
    if buildDir != name:
        # Yep, some packages untar to the same directory as the package
        # name (hello Xmipp), and that is not so great. No link to it.
        tLink = env.Command(
            Dir('#software/em/%s/bin').abspath % name,  # TODO: find smtg better than "/bin"
            Dir('#software/em/%s' % buildDir),
            Action('rm -rf %s && ln -v -s %s %s' % (name, buildDir, name),
                   'Linking package %s to software/em/%s' % (name, name),
                   chdir='software/em'))
    else:
        tLink = tUntar  # just so the targets are properly connected later on
    SideEffect('dummy', tLink)  # so it works fine in parallel builds
    lastTarget = tLink
    
    # Load Package vars
    # First we search this in cfg folder and otherwise in package home
    if not os.path.exists(confPath):
        confPath = join(packageHome, '%s.cfg' % name)
    if not os.path.exists(confPath):
        confPath = 'unset'
    opts = Variables(confPath)
    opts.Add('PACKAGE_SCRIPT')
    opts.Add('PRIVATE_KEYS')
    opts.Add('BOOL_PRIVATE_KEYS')
    opts.Update(env) 
    for var in env.get('PRIVATE_KEYS'):
        opts.Add(var)
        opts.Update(env)
        opts.Add(*(env.get(var)))
    for var in env.get('BOOL_PRIVATE_KEYS'):
        opts.Add(var)
        opts.Update(env)
        opts.Add(BoolVariable(*env.get(var)))
    opts.Update(env)
    opts.Save(confPath + "_tmp", env)
    Help(opts.GenerateHelpText(env, sort=cmp))
#    print opts.keys()
    scriptPath = env.get('PACKAGE_SCRIPT')
    altScriptPath = join(packageHome, 'SConscript')
    # FIXME: this scriptPath wont be reachable uless Xmipp is already downloaded
    env.Replace(packageDeps=lastTarget)
    if scriptPath and os.path.exists(scriptPath):
        print "Reading SCons script file at %s" % scriptPath
        lastTarget = env.SConscript(scriptPath, exports='env')
    elif os.path.exists(altScriptPath):
        print "Config file not present. Trying SConscript in package home %s" % altScriptPath
        lastTarget = env.SConscript(altScriptPath, exports='env')
    else:
        # If we can't find the SConscript, then only extra actions can be done
        for target, command in extraActions:
            lastTarget = env.Command('software/em/%s/%s' % (name, target),
                                     lastTarget,
                                     Action(command, chdir='software/em/%s' % name))
            SideEffect('dummy', lastTarget)  # so it works fine in parallel builds

    # Clean the special generated files
    for cFile in clean:
        Clean(lastTarget, cFile)
    # Add the dependencies. Do it to the "link target" (tLink), so any
    # extra actions (like setup scripts) have everything in place.
    for dep in deps:
        env.Depends(lastTarget, dep)

    if (default or packageHome):
        for lastT in lastTarget:
            if lastT is not None:
                env.Default(lastT)

    return lastTarget


def manualInstall(env, name, tar=None, buildDir=None, url=None,
                  extraActions=[], deps=[], clean=[], default=True):
    """Just download and run extraActions.

    This pseudobuilder downloads the given url, untars the resulting
    tar file and runs extraActions on it. It also tells SCons about
    the proper dependencies (deps).

    extraActions is a list of (target, command) that should be
    executed after the package is properly installed.

    If default=False, the package will not be built unless the option
    --with-<name> is used.

    Returns the final target in extraActions.

    """
    # Use reasonable defaults.
    tar = tar or ('%s.tgz' % name)
    url = url or ('%s/external/%s' % (URL_BASE, tar))
    buildDir = buildDir or tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0]

    # Add the option --with-name, so the user can call SCons with this
    # to activate it even if it is not on by default.
    if not default:
        AddOption('--with-%s' % name, dest=name, action='store_true',
                  help='Activate %s' % name)

    # Donload, untar, and execute any extra actions.
    tDownload = Download(env, 'software/tmp/%s' % tar, Value(url))
    SideEffect('dummy', tDownload)  # so it works fine in parallel builds
    tUntar = Untar(env, 'software/tmp/%s/README' % buildDir, tDownload,
                   cdir='software/tmp')
    SideEffect('dummy', tUntar)  # so it works fine in parallel builds
    Clean(tUntar, 'software/tmp/%s' % buildDir)
    lastTarget = tUntar

    for target, command in extraActions:
        lastTarget = env.Command(
            target,
            lastTarget,
            Action(command, chdir='software/tmp/%s' % buildDir))
        SideEffect('dummy', lastTarget)  # so it works fine in parallel builds

    # Clean the special generated files.
    for cFile in clean:
        Clean(lastTarget, cFile)

    # Add the dependencies.
    for dep in deps:
        env.Depends(tUntar, dep)

    if default or GetOption(name):
        env.Default(lastTarget)

    return lastTarget


# Add methods so SConscript can call them.
env.AddMethod(addLibrary, 'AddLibrary')
env.AddMethod(addModule, 'AddModule')
env.AddMethod(addPackage, 'AddPackage')
env.AddMethod(manualInstall, 'ManualInstall')
env.AddMethod(compilerConfig, 'CompilerConfig')
env.AddMethod(addPackageLibrary, 'AddPackageLibrary')
env.AddMethod(addJavaLibrary, 'AddJavaLibrary')
env.AddMethod(symLink, 'SymLink')
env.AddMethod(addJavaTest, 'AddJavaTest')
env.AddMethod(addProgram, 'AddProgram')
env.AddMethod(progInPath, 'ProgInPath')


#  ************************************************************************
#  *                                                                      *
#  *                            Extra options                             *
#  *                                                                      *
#  ************************************************************************


opts = Variables(None, ARGUMENTS)

opts.Add('SCIPION_HOME', 'Scipion base directory', abspath('.'))
opts.Add('JAVAC', 'Java compiler', 'javac')
opts.Add('MPI_CC', 'MPI C compiler', 'mpicc')
opts.Add('MPI_CXX', 'MPI C++ compiler', 'mpiCC')
opts.Add('MPI_LINKERFORPROGRAMS', 'MPI Linker for programs', 'mpiCC')
opts.Add('MPI_INCLUDE', 'MPI headers dir ', '/usr/include')
opts.Add('MPI_LIBDIR', 'MPI libraries dir ', '/usr/lib')
opts.Add('MPI_LIB', 'MPI library', 'mpi')
opts.Add('MPI_BINDIR', 'MPI binaries', '/usr/bin')

opts.Update(env)

Help('\nVariables that can be set:\n')
Help(opts.GenerateHelpText(env))
Help('\n')




AddOption('--with-all-packages', dest='withAllPackages', action='store_true',
          help='Get all EM packages')

env.SConscript('SConscript', exports='env')

# Add original help (the one that we would have if we didn't use
# Help() before). But remove the "usage:" part (first line).
phelp = SCons.Script.Main.OptionsParser.format_help().split('\n')
Help('\n'.join(phelp[1:]))
# This is kind of a hack, because the #@!^ scons doesn't give you easy
# access to the original help message.