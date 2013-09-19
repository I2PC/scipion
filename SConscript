#!/usr/bin/env python

Import('env')

# Required for custom functions
import os
from os.path import join
from types import *

FFTWDir = "external/fftw-3.3.1"
FFTWLibs = ['fftw3', 'fftw3_threads']
TIFFDir = "external/tiff-3.9.4"
TIFFLibs = ['tiff']
JPEGDir = "external/jpeg-8c"
JPEGLibs = ['jpeg']
HDF5Dir = "external/hdf5-1.8.10/src/"
HDF5Libs = ['hdf5', 'hdf5_cpp']
CYGWIN = env['PLATFORM'] == 'cygwin'
MACOSX = env['PLATFORM'] == 'darwin'
MINGW = env['PLATFORM'] == 'win32'
if CYGWIN:
	TIFFLibs.append('z')
ARPACKppDir = "external/arpack++-2.3"
ARPACKppLibs = ['arpack++']
SQliteDir = "external/sqlite-3.6.23"
SQLiteLibs = ['sqlite3']

#PYSQliteDir = "external/pysqlite-2.6.3"

PythonDir = "external/python/Python-2.7.2"
PythonLibs = []

PluginLibs = {}
PluginResources = {}
javaEnumDict = {'ImageWriteMode': ['libraries/data/xmipp_image_base.h', 'WRITE_'],
            'CastWriteMode': ['libraries/data/xmipp_image_base.h', 'CW_'],
            'MDLabel': ['libraries/data/metadata_label.h', 'MDL_'],
            'XmippError': ['libraries/data/xmipp_error.h', 'ERR_']}

copyJar = None

def AddMatchingFiles((pattern, blacklist, sources), directory, files):
    ''' Callback, adds all matching files in dir '''
    import fnmatch
    for file in fnmatch.filter(files, pattern):
        if file not in blacklist:
           # DBG print 'Adding ' + join(directory, file)
           sources.append(join(directory, file))

def Glob(dir, pattern, blacklist):
    ''' Custom made globbing '''
    
    sources = []
    os.path.walk(dir, AddMatchingFiles, (pattern, blacklist, sources))
    return sources

def AddLastSlash(string):
    ''' Low trick for correct parsing of paths '''
    str = string.strip();
    if len(str) == 0:
        return "";
    if not str.endswith('/'):
        str = str + '/'
    return str

def SymLink(target, source, env=None):
    #As the link will be in bin/ directory we need to move up
    if isinstance(target, list):
        link = str(target[0])
        source = str(source[0])
    else:
        link = target
    
    if CYGWIN:
    	from shutil import copyfile
        copyfile(source, link)
    else:    
        source = os.path.relpath(source, os.path.split(link)[0])
        if os.path.lexists(link):
            os.remove(link)
        os.symlink(source, link)

def AddProgramLink(target, source, PrependXmipp=True):
    binprefix = join(env['prefix'], 'bin')
    bintarget = join(binprefix, env['prepend'] + target)
    if PrependXmipp:
        binsource = join(binprefix, env['prepend'] + source)
    else:
        binsource = source
    command = env.Command(bintarget, binsource,
        [Chmod('$SOURCE', 0755), SymLink])
    alias = env.Alias(target, command)
    env.Default(alias)
    return command

def AddBatch(name, basedir, extension=''):
    # setup
    basedir = AddLastSlash(basedir)
    fullname = env['prepend'] + name
    binprefix = join(env['prefix'], 'bin')

    # action
    source = join(basedir, 'batch_' + name + extension)
    target = join(binprefix, fullname)
    command = env.Command(target, source,
        [Chmod('$SOURCE', 0755), SymLink])

    # alias
    alias = env.Alias(fullname, command)
    # Consider batch also as xmipp_programs
    env.Alias('xmipp_programs', command)
    #install = env.Install(binprefix, command)
    #env.Alias(fullname, install)
    env.Default(alias)
    return command

def AddProgram(name, basedir, sources_pattern='*.cpp', skip_list=[],
               includes=[], libpath=[], libs=[], cxxflags=[],
               linkflags=[]):
    ''' add a new program to the build list '''
    # setup
    basedir = AddLastSlash(basedir)
    fullname = env['prepend'] + name
    sources = Glob(basedir, sources_pattern, skip_list)
    binprefix = join(env['prefix'], 'bin')

    # FIXME fix for static executables
    if env['static']:
        cxxflags += [env['STATIC_FLAG']]
        linkflags += [env['STATIC_FLAG']]
    Decider('MD5-timestamp')
    # action
    program = env.Program(
        join(basedir, fullname),
        sources,
        CPPPATH=includes + [env['CPPPATH']],
        LIBPATH=libpath + [env['LIBPATH']],
        LIBS=libs + [env['LIBS']],
        CXXFLAGS=cxxflags + [env['CXXFLAGS']],
        LINKFLAGS=linkflags + [env['LINKFLAGS']],
        LINK=env['LINKERFORPROGRAMS']
        )

    install = env.Install(binprefix, program)
    alias = env.Alias(fullname, install)
    env.Default(alias)
    return alias

def AddMPIProgram(name, basedir, sources_pattern='*.cpp', skip_list=[],
                  includes=[], libpath=[], libs=[], cxxflags=[],
                  linkflags=[]):

    # setup
    basedir = AddLastSlash(basedir)
    fullname = env['prepend'] + name
    sources = Glob(basedir, sources_pattern, skip_list)
    binprefix = join(env['prefix'], 'bin')

    # FIXME fix for static executables
    if env['static']:
        cxxflags += [env['STATIC_FLAG']]
        linkflags += [env['STATIC_FLAG']]

    # action
    program = env.Program(
        join(basedir, fullname),
        sources,
        CC=env['MPI_CC'],
        CXX=env['MPI_CXX'],
        CPPPATH=includes + [env['CPPPATH']] + [env['MPI_INCLUDE']],
        LIBPATH=libpath + [env['LIBPATH']] + [env['MPI_LIBDIR']],
        LIBS=libs + [env['LIBS']] + [env['MPI_LIB']],
        CXXFLAGS=cxxflags + [env['CXXFLAGS']],
        LINKFLAGS=linkflags + [env['LINKFLAGS']],
        LINK=env['MPI_LINKERFORPROGRAMS'],
        LD_LIBRARY_PATH=[env['LIBPATH']] + [env['MPI_LIBDIR']]
        )

    # alias
    alias = env.Alias(fullname, program)
    install = env.Install(binprefix, program)
    env.Alias(fullname, install)
    env.Default(alias)

# Add a program integrated in the Xmipp structure
def AddXmippProgram(name, libs=[], folder='programs', incPaths=[], libPaths=[],
                    useCudaEnvironment=False):
    finalLibPath = ['lib']
    finalLibPath.append(libPaths)
    finalIncludePath = ['libraries', '#', '#'+HDF5Dir]
    finalIncludePath.append(incPaths)
    finalLibs = libs + ['XmippData', 'XmippExternal'] + FFTWLibs + SQLiteLibs + TIFFLibs + JPEGLibs + HDF5Libs
    if useCudaEnvironment:
    	finalLibs += ['cudart', 'cutil', 'shrutil_x86_64' ]
    	finalIncludePath += [env['CUDA_SDK_PATH'] + "/CUDALibraries/common/inc",
                           env['CUDA_SDK_PATH'] + "/shared/inc"]
    	finalLibPath += [env['CUDA_SDK_PATH'] + "/CUDALibraries/common/lib",
                       env['CUDA_SDK_PATH'] + "/shared/lib",
                       env['CUDA_SDK_PATH'] + "/CUDALibraries/common/lib/linux",
		       "/usr/local/cuda/lib64",
                       env['CUDA_LIB_PATH']]
    if 'XmippRecons' in finalLibs and not 'XmippClassif' in finalLibs:
        finalLibs.append('XmippClassif')
    if 'XmippRecons' in finalLibs and int(env['cuda']):
	finalLibs.append("XmippReconsCuda");
    if int(env["arpack"]):
        finalLibs += ['arpack++', 'arpack', 'lapack', 'blas']
    program = AddProgram(name, 'applications/%s/%s' % (folder, name), '*.cpp', [],
        finalIncludePath, finalLibPath, finalLibs, [], [])
    env.Alias('xmipp_programs', program)
    return program
	
# Add a program integrated in the Xmipp structure
def AddXmippTest(name, testprog, command):
    #testprog = AddXmippProgram(name, ['gtest'], 'tests')
    testname = 'xmipp_' + name
    xmlFileName = 'applications/tests/OUTPUT/' + testname + ".xml"
    if  os.path.exists(xmlFileName):
       os.remove(xmlFileName)
    testcase = env.Alias('run_' + name , env.Command(xmlFileName, testname, command))
    env.Depends(testcase, testprog)
    test = env.Alias('run_tests', testcase)
    AlwaysBuild(testcase)
    return testcase


def AddXmippCTest(name):
    testprog = AddXmippProgram(name, ['gtest', 'XmippRecons','XmippDimred'], 'tests')
    AddXmippTest(name, testprog, "$SOURCE --gtest_output=xml:$TARGET")

def AddXmippPythonTest(name):
    #print "Adding python test: ", name
    #FIXME ROB
    testprog = AddBatch(name, 'applications/tests/' + name, '.py')
    test = AddXmippTest(name, testprog, "$SOURCE $TARGET")
    return test


def AddXmippJavaTest(name):
    #javac xmipp/ImageGeneric_Test.java -classpath ~/xmipp_svn/java/lib/XmippJNI.jar:/usr/share/java/junit4.jar
    #mv xmipp/ImageGeneric_Test.class xmipp/jni/.
    #java  -classpath :/usr/share/java/junit4.jar:xmipp:/home/roberto/xmipp_svn/java/lib/XmippJNI.jar   org.junit.runner.JUnitCore xmipp.jni.ImageGeneric_Test  
    pass

    #env.Default(test)

def AddXmippMPIProgram(name, libs=[]):
    finalLibPath = ['lib']
    finalIncludePath = ['libraries', '#', '#'+HDF5Dir]
    finalLibs = libs + ['XmippData', 'XmippExternal', 'XmippParallel'] + FFTWLibs + SQLiteLibs + TIFFLibs + JPEGLibs + HDF5Libs
    if int(env["arpack"]):
        finalLibs += ['arpack++', 'arpack', 'lapack', 'blas']
    if 'XmippRecons' in finalLibs and not 'XmippClassif' in finalLibs:
        finalLibs.append('XmippClassif')
    if 'XmippRecons' in finalLibs and int(env['cuda']):
		finalLibs.append("XmippReconsCuda");
    for i in range(len(libs)):
       if libs[i] == 'XmippRecons':
          finalLibPath += ['libraries/reconstruction']
       elif libs[i] == 'XmippInterface':
          finalLibPath += ['libraries/interface']
       elif libs[i] == 'XmippRecons_Interface':
          finalLibPath += ['libraries/interface']
          finalLibs.insert(i + 1, 'XmippInterface')
       elif libs[i] == 'XmippReconsMPI':
          finalLibPath += ['libraries/reconstruction_mpi']
       elif libs[i] == 'XmippClassif':
          finalLibPath += ['libraries/classification']
    AddMPIProgram(name, 'applications/programs/' + name, '*.cpp', [],
        finalIncludePath, finalLibPath, finalLibs, [], [])

# For Roberto's new lib
def AddMPILibrary(name, basedir, sources, includes, libpath=[], libs=[]):
    # setup
    basedir = AddLastSlash(basedir)
    libprefix = join(env['prefix'], 'lib')
    #for x in sources:
    #    sources[sources.index(x)] = basedir + x

    # separate local and global includes
    for x in includes:
        if x[0] != '#':
            includes[includes.index(x)] = basedir + x

    # action
    # FIXME Exclusive may not be what users want
    if int(env['static']):
        library = env.StaticLibrary(
            basedir + name,
            sources,
            CPPPATH=includes + [env['CPPPATH']] + [env['MPI_INCLUDE']],
            CC=env['MPI_CC'],
            CXX=env['MPI_CXX'],
            LIBPATH=[env['MPI_LIBDIR']] + libpath,
            LIBS=[env['MPI_LIB']] + libs
            )
    else:
        library = env.SharedLibrary(
            basedir + name,
            sources,
            CPPPATH=includes + [env['CPPPATH']] + [env['MPI_INCLUDE']],
            CC=env['MPI_CC'],
            CXX=env['MPI_CXX'],
            LIBPATH=[env['MPI_LIBDIR']] + libpath,
            LIBS=[env['MPI_LIB']] + libs
            )

    # alias
    alias = env.Alias(name, library)
    install = env.Install(libprefix, library)
    env.Alias(name, install)
    env.Default(alias)

def AddLibrary(name, basedir, sources, includes, libpath=[], libs=[],
    shlibprefix='lib', shlibsuffix=env['SHLIBSUFFIX'], useCudaEnvironment=False):
    # setup
    basedir = AddLastSlash(basedir)
    libprefix = join(env['prefix'], 'lib')
    cudaFiles = False
    for x in sources:
        if x.find(basedir) == -1:
            sources[sources.index(x)] = basedir + x
	if x.endswith(".cu"):
	    cudaFiles = True

    # separate local and global includes
    for x in includes:
        if type(x) is ListType: 
            y=x
            includes.remove(x)
            for j in y:
                includes.append(j)
        else:
            if x[0] != '#' and x[0] != '/':
                includes[includes.index(x)] = basedir + x

    if useCudaEnvironment:
        envToUse = env.Clone()
        envToUse.Tool('cuda')
        envToUse.Append(CXXFLAGS=['-DWITH_CUDA'])
	if cudaFiles:
            envToUse.Append(NVCCFLAGS="-shared --compiler-options '-fPIC'")
        libs += ['cutil_x86_64', 'shrutil_x86_64', 'cudart']
    	includes += [env['CUDA_SDK_PATH'] + "/CUDALibraries/common/inc",
                   env['CUDA_SDK_PATH'] + "/shared/inc"]
    	libpath += [env['CUDA_SDK_PATH'] + "/CUDALibraries/common/lib",
                       env['CUDA_SDK_PATH'] + "/shared/lib",
                       env['CUDA_SDK_PATH'] + "/CUDALibraries/common/lib/linux",
		       "/usr/local/cuda/lib64",
                       env['CUDA_LIB_PATH']]
    else:
        envToUse = env

    if int(envToUse['static']):
        library = envToUse.StaticLibrary(
            basedir + name,
            sources,
            CPPPATH=includes + [envToUse['CPPPATH']],
            LIBPATH=libpath,
            LIBS=libs
            )
    else:
        library = envToUse.SharedLibrary(
            basedir + name,
            sources,
            CPPPATH=includes + [envToUse['CPPPATH']],
            LIBPATH=libpath,
            LIBS=libs,
            SHLIBPREFIX=shlibprefix,
            SHLIBSUFFIX=shlibsuffix
            )
    for lib in libs:
        for ext in ['.a', '.so', '.dll', '.dylib']:
            rootname = 'lib/lib' + lib
            if os.path.exists(rootname + ext):
                Depends(library, rootname + ext)

    install = envToUse.Install(libprefix, library)
    alias = envToUse.Alias(name, install)
    envToUse.Default(alias)
    return alias

def CreateFileList(path, pattern, filename, root='', root2=''):
    fOut = open(filename, 'w+')
    files = [f.replace(root, root2) + '\n' for f in Glob(path, pattern, [])]
#    print '************************************************'
#    print 'path', path
#    print 'pattern', pattern
#    print 'filename', filename
#    print 'files', files
#    print '************************************************'
    fOut.writelines(files)
    fOut.close()
    
def Cmd(cmd):
    print cmd
    os.system(cmd)
    
def javaBuildName(*args):
    return join(env['JAVADIR'], 'build', *args)
 
def javaLibName(*args):   
    return join(env['JAVADIR'], 'lib', *args)

def javaSrcName(*args):
    return join(env['JAVADIR'], 'src', *args)

def CompileJava(name, dependencies=[]):
    ''' name parameter is expected without .java extension '''
    source = javaSrcName(name + '.java')
    target = javaBuildName(name + '.class')
    buildDir = javaBuildName()
    cmd = env['JAVAC'] + ' -cp "java/lib/*"'
    compileCmd = env.Command(target, source, '%(cmd)s -d %(buildDir)s %(source)s' % locals())
    dependencies.append('XmippJNI') # XmippJNI dependency is added by default
    deps = [javaLibName(name + '.jar') for name in dependencies]
    for lib in deps:
        env.Depends(compileCmd, lib)
    env.Default(compileCmd)
    return compileCmd

def CompileJavaJar(target, source, env):    
    srcDir = str(source[0])
    buildDir = javaBuildName()
    jarfile = str(target[0])
    name = os.path.basename(jarfile)
    listfile = javaBuildName(name + '_source.txt')
    classfile = javaBuildName(name + '_classes.txt')
    CreateFileList(srcDir, '*.java', listfile)
    Cmd(env['JAVAC'] + ' -cp "java/lib/*" -d %(buildDir)s -sourcepath %(srcDir)s @%(listfile)s' % locals())
    
    classDir = join(buildDir, os.path.relpath(srcDir, javaSrcName()))
    # This is needed for compiling IJ plugins
    # where the file 'plugins.config' need to be include in the final .jar file
    configFile = 'plugins.config'
    pluginDest = ''
    if os.path.exists(join(srcDir, configFile)):
    	pluginDest = join(classDir, configFile)
    	Cmd('cp %s %s' % (join(srcDir, configFile), pluginDest))
    # We need to create a txt file with the list of all .class files 
    # with the following format:
    # -C java/build xmipp/package/A.class
    # -C java/build xmipp/package/B.class
    # ...
    cmd = 'jar ' + env['JARFLAGS']
    CreateFileList(classDir, '*.class', classfile, buildDir + '/', '-C %(buildDir)s ' % locals())
    Cmd('%(cmd)s %(jarfile)s @%(classfile)s %(pluginDest)s' % locals())
    
def AddJavaLibrary(name, src, dependencies=[]):#, sourceList, includes, libpath=[], libs=[]):
    # all libraries are assumed to be inside xmipp
    srcDir = javaSrcName('xmipp', src)
#    listfile = javaBuildName(name + '_source.txt')
#    classfile = javaBuildName(name + '_classes.txt')    
    jarfile = javaLibName(name + '.jar')
    if name != 'XmippJNI': # XmippJNI dependency is added by default
        dependencies.append('XmippJNI')
    deps = [javaLibName(name + '.jar') for name in dependencies]
    buildJar = env.Command(jarfile, [srcDir] + deps, CompileJavaJar)
    #Add sources dependencies, not handled now by Scons
    sources = Glob(srcDir, '*.java', [])
    for s in sources:
        env.Depends(buildJar, s)
    for lib in deps:
        env.Depends(buildJar, lib)
    env.Alias(name + '.jar', buildJar)
    # Group all jar targets under 'java' alias
    env.Alias('java', buildJar)
    env.Default(buildJar)
    return buildJar

def AddJavaTest(testName):
    cmd = 'java -cp "java/lib/*" org.junit.runner.JUnitCore xmipp.test.' + testName
    runTest = env.Command(testName, javaLibName('XmippTest.jar'), cmd)
    env.Alias('run_java_tests', runTest)

def AddJavaIJPlugin(name, src, dependencies=[]):
    ''' this function does the same of AddJavaLibrary and 
    add a link to the .jar from ij/plugins'''
    jarfile = name + '.jar'
    buildJar = AddJavaLibrary(name, src, dependencies)
    copyIJPlugin = env.Command(join('external/imagej/plugins', jarfile), buildJar, SymLink)
    env.Alias('java', copyIJPlugin)
    env.Default(copyIJPlugin)    
    
def removeAll(basedir, regexp):
	import glob
	files = glob.glob(basedir + regexp)
	for i in range(len(files)):
		os.remove(files[i])

def AddJavaApp(name, appDir, outDir):
    import shutil
    sourceDir = join(appDir, 'src')
    resourcesDir = join(sourceDir, 'resources')
    libsDir = join(appDir , 'libs')
    metainfDir = join(sourceDir, 'META-INF')
    buildDir = join(appDir , 'build')
    buildResourcesDir = join(buildDir, 'resources')
    buildMetainf = join(buildDir, 'META-INF')
    outLibsDir = join(outDir, 'libs')
    # Clears "buildDir"
    if not os.path.exists(buildDir):
        os.mkdir(buildDir)
        #os.mkdir(buildResourcesDir)

    env.Append(JAVACLASSPATH=join(libsDir, '*'))
    buildClasses = env.Java(buildDir, sourceDir, JAVAVERSION='1.6')
    # Copies MANIFEST
    copyMetainf = env.Install(buildMetainf, Glob(metainfDir, "MANIFEST.MF", []))
    buildJar = env.Jar(join(outDir, name + '.jar'), buildDir, JARCHDIR=buildDir)

    appLibs = {}
    jarList = Glob(libsDir, "*.jar", [])
    for l in jarList:
        lname = os.path.basename(l)
        if not appLibs.has_key(lname):
            install = env.Install(outLibsDir, l)
            appLibs[lname] = install
        env.Depends(buildClasses, appLibs[lname])

    copyPluginResources = env.Install(buildResourcesDir, Glob(resourcesDir, "*?.???", []))
    env.Depends(buildJar, [copyMetainf, copyPluginResources])#, buildClasses])
    pluginAlias = env.Alias(name, buildJar)
    env.Depends(buildClasses, copyJar)

    return pluginAlias

def AddImageJPlugin(name, pluginDir, outDir, requiredPlugins=[]):
    import shutil
    sourceDir = join(pluginDir, 'src')
    resourcesDir = join(sourceDir, 'resources')
    libsDir = join(pluginDir , 'libs')
    macrosDir = join(pluginDir , 'macros')
    buildDir = join(pluginDir , 'build')
    buildResourcesDir = join(buildDir, 'resources')
    # Clears "buildDir"
    if not os.path.exists(buildDir):
        os.mkdir(buildDir)
        os.mkdir(buildResourcesDir)
    # Copies 'plugins.config' to build dir...
    buildClasses = env.Java(buildDir, sourceDir, JAVAVERSION='1.6')
    env.Depends(buildClasses, requiredPlugins)
    copyPluginConfig = env.Install(buildDir, join(sourceDir, 'plugins.config'))#Glob(sourceDir, 'plugins.config', []))
    buildJar = env.Jar(join(outDir, name + '.jar'), buildDir, JARCHDIR=buildDir)

    jarList = Glob(libsDir, "*.jar", [])
    for l in jarList:
        lname = os.path.basename(l)
        if not PluginLibs.has_key(lname):
            install = env.Install(outDir, l)
            PluginLibs[lname] = install
        env.Depends(buildClasses, PluginLibs[lname])

    copyPluginResources = env.Install(buildResourcesDir, Glob(resourcesDir, "*?.???", []))
    env.Depends(buildJar, [copyPluginConfig, copyPluginResources, buildClasses])
    pluginAlias = env.Alias(name, buildJar)
    env.Depends(buildClasses, copyJar)

    outMacrosDir = join(outDir, "..", "macros")
    env.Alias('copyMacros', env.Install(outMacrosDir, Glob(macrosDir, "*?.???", [])))
    return pluginAlias

# Gtest
if int(env['gtest']):
    DataSources = Glob('external/gtest-1.6.0/fused-src/gtest', 'gtest-all.cc', [])
    AddLibrary('gtest', 'external/gtest-1.6.0/fused-src/gtest', DataSources, ['#'],
               [], [])

# Bilib
BilibSources = Glob('external/bilib/sources', '*.cc', [])

# INRIA
INRIASources = Glob('external/inria', '*.cc', [])

# Condor
CondorSources = Glob('external/condor', '*.cpp', [])

# AlgLib
AlglibSources = Glob('external/alglib/src', '*.cpp', [])

AddLibrary('XmippExternal', 'external',
   INRIASources + BilibSources + CondorSources + AlglibSources,
   ['bilib', 'bilib/headers', 'bilib/types'])

# sqliteExt
SqliteExtSources = Glob('external/sqliteExt', '*.c', [])
AddLibrary('XmippSqliteExt', 'external',
   SqliteExtSources,
   ['#'],[''],['-lm'],'lib','.so')

# XmippData
DataSources = Glob('libraries/data', '*.cpp', [])

libraries=['#', '#'+HDF5Dir]
if MINGW:
    import sys
    sys.setrecursionlimit(22500)
    libraries.append(env['MINGW_PATHS'])
    AddLibrary('XmippData', '', DataSources, libraries, 
               ['lib'], ['XmippExternal','regex','rt'] + FFTWLibs + TIFFLibs + JPEGLibs + HDF5Libs + SQLiteLibs)
elif MACOSX:
    AddLibrary('XmippData', 'libraries/data', DataSources, libraries,
               ['lib'], ['XmippExternal'] + FFTWLibs + TIFFLibs + JPEGLibs + HDF5Libs + SQLiteLibs)  
else:
    AddLibrary('XmippData', 'libraries/data', DataSources, libraries,
               ['lib'], ['XmippExternal','rt'] + FFTWLibs + TIFFLibs + JPEGLibs + HDF5Libs + SQLiteLibs)  


#Xmipp Python Extension
PyExtSources = Glob('libraries/bindings/python', '*.cpp', [])
#import distutils.sysconfig
pythonLibName = 'xmipp'
pythonIncludes = ["#" + join(PythonDir, dir) for dir in [".", "Include"]]
pythonIncludes.append("#lib/python2.7/site-packages/numpy/core/include") 

libpath = ['lib']
libraries = ['XmippData', 'XmippRecons', 'XmippExternal'] + FFTWLibs + TIFFLibs + JPEGLibs + HDF5Libs + SQLiteLibs
if CYGWIN or MACOSX:
    libpath.append(PythonDir)
    libraries.append("python2.7")

pythonbinding = AddLibrary(pythonLibName, 'libraries/bindings/python', PyExtSources,
           ['#libraries', "#", '#'+HDF5Dir] + pythonIncludes,
           libpath, libraries, '')

# in MACOSX Python requires module libraries as .so instead of .dylib
if MACOSX:
	command = env.Command('lib/' + pythonLibName + '.so', 'lib/' + pythonLibName + '.dylib', SymLink)
	env.Alias(pythonLibName, command)

# Reconstruction
ReconsSources = Glob('libraries/reconstruction', '*.cpp', ["angular_gcar.cpp"])
ReconsLib = ['XmippExternal', 'XmippData', 'pthread', 'XmippClassif'] + FFTWLibs + TIFFLibs + JPEGLibs + HDF5Libs + SQLiteLibs
ReconsIncDir = ['#libraries', '#external', '#', '#'+HDF5Dir]
ReconsLibDir = ['lib']
if int(env['arpack']):
    ReconsSources.append("angular_gcar.cpp")
    ReconsLib += ['arpack', 'lapack', 'blas']
AddLibrary('XmippRecons', 'libraries/reconstruction', ReconsSources,
           ReconsIncDir, ReconsLibDir, ReconsLib, useCudaEnvironment=int(env['cuda']))

if int(env['cuda']):
    ReconsCudaSources = Glob('libraries/reconstruction', '*.cu', [])
    AddLibrary('XmippReconsCuda', 'libraries/reconstruction', ReconsCudaSources,
           ['#'], ['lib'], [], useCudaEnvironment=True)

# Classification
ClassificationSources = Glob('libraries/classification', '*.cpp', [])
AddLibrary('XmippClassif', 'libraries/classification', ClassificationSources,
    ['#libraries', '#', '#'+HDF5Dir], ['lib'], ['XmippExternal', 'XmippData'])

# Dimensionality reduction
DimRedSources = Glob('libraries/dimred', '*.cpp', [])
AddLibrary('XmippDimred', 'libraries/dimred', DimRedSources,
    ['#libraries', '#'], ['lib'], ['XmippExternal', 'XmippData'])

# XmippParallel
ParallelSources = Glob('libraries/parallel', '*.cpp', []);
AddMPILibrary("XmippParallel", 'libraries/parallel', ParallelSources, ["#", "#libraries", "#external", '#'+HDF5Dir],
              ['lib'], ['XmippExternal', 'XmippData', 'XmippRecons', 'XmippClassif'] + FFTWLibs + TIFFLibs + JPEGLibs + HDF5Libs + SQLiteLibs)

# Interface
InterfaceSources = Glob('libraries/interface', '*.cpp', [])
AddLibrary('XmippInterface', 'libraries/interface', InterfaceSources,
    ['#libraries', '#external', '#', '#'+HDF5Dir], ['lib'], ['XmippExternal', 'XmippData', 'pthread'])

# Recons Interface
#AddLibrary('XmippRecons_Interface', '#libraries/reconstruction',
#           ReconsInterfaceSources, ['#libraries', '#external', '#'],['lib'],['XmippExternal','XmippData','XmippRecons','XmippInterface'])

def WriteJavaEnum(class_name, header_file, pattern, log):
    java_file = "java/src/xmipp/jni/%s.java" % class_name
    env.Depends(java_file, header_file)
    f = open(header_file)
    fOut = open(java_file, 'w+')
    counter = 0;
    last_label_pattern = pattern + "LAST_LABEL"
    fOut.write("package xmipp.jni; \n")
    fOut.write("public class " + class_name + " {\n")

    for line in f:
        l = line.strip();
        if l.startswith(pattern):
            if '///' in l:
                l, comment = l.split('///')
            else:
                comment = ''
            if l.startswith(last_label_pattern):
                l = l.replace(last_label_pattern, last_label_pattern + " = " + str(counter) + ";")
            if (l.find("=") == -1):
                l = l.replace(",", " = %d;" % counter)
                counter = counter + 1;
            else:
                l = l.replace(",", ";")

            fOut.write("   public static final int %s ///%s\n" % (l, comment))
    fOut.write("}\n")
    fOut.close()
    f.close()
    # Write log file
    if log:
        from datetime import datetime
        d = str(datetime.now())
        #d = date.today();
        log.write("Java file '%s' successful generated at %s\n" % (java_file, d))

def ExtractEnumFromHeader(source, target, env):
    # this is very ugly, we still need more scons knowledge

    log = open(str(target[0]), 'w+')
    for (class_name, list) in javaEnumDict.iteritems():
        WriteJavaEnum(class_name, list[0], list[1], log)

    log.close()
    return None



if int(env['java']):
    libDir = 'java/lib'
    env.Append(JAVACLASSPATH=libDir)
    env['JAVABUILDPATH'] = 'java/build'
    env['JAVADIR'] = 'java'
    env['ENV']['LANG'] = 'en_GB.UTF-8'
    env['JARFLAGS'] = '-Mcf'	# Default "cf". "M" = Do not add a manifest file.
    buildDir = env['JAVABUILDPATH']
    # Create the build dir if not exist
    if not os.path.exists(buildDir): # create classes dir if not exists
        Execute(Mkdir(buildDir))
    # Set -g debug options if debugging
    if env['debug'] == True:
        env['JAVAC'] = 'javac -g'

    # Update enums from c++ headers, if not exist, generate it
    for (class_name, class_list) in javaEnumDict.iteritems():
        java_file = "java/src/xmipp/jni/%s.java" % class_name
        #if not os.path.exists(java_file):
        WriteJavaEnum(class_name, class_list[0], class_list[1], None)

    env.Alias('javaEnums', env.Command("libraries/bindings/java/src/xmipp/jni/enums.changelog",
                            ["libraries/data/xmipp_image_base.h", "libraries/data/metadata_label.h" ], ExtractEnumFromHeader))

    JavaInterfaceSources = Glob('libraries/bindings/java', '*.cpp', [])
    JavaDependLibraries = ['XmippData', 'pthread', 'XmippRecons', 'XmippClassif', 'XmippExternal']
    if int(env['arpack']):
        JavaDependLibraries += ['arpack++', 'arpack', 'lapack', 'blas']
    # Compilation of the c code needed for java jni binding
    javaJniC = AddLibrary('XmippJNI', 'libraries/bindings/java', JavaInterfaceSources, ['#libraries', '#external', '#', '#'+HDF5Dir] + env['JNI_CPPPATH'], ['lib'],JavaDependLibraries)

    # Create some jar links
    cmd = env.Command(join(libDir, 'ij.jar'), 'external/imagej/ij.jar', SymLink)
    env.Default(cmd)
    javaJni = AddJavaLibrary('XmippJNI', 'jni', ['ij'])
    env.Depends(javaJni, javaJniC)
    AddJavaLibrary('XmippUtils', 'utils', ['ij', 'commons-cli-1.1'])
    AddJavaLibrary('XmippIJ', 'ij/commons', ['XmippUtils'])
    AddJavaLibrary('XmippViewer', 'viewer', ['XmippIJ', 'XmippUtils', 'ij', 'jfreechart-1.0.13'])
    AddJavaLibrary('XmippTest', 'test', ['XmippViewer', 'junit4-4.8.2', 'core-1.1'])
    AddJavaIJPlugin('XmippIJPlugin_MasksToolbar', 'ij/plugins/maskstoolbar', [])
    #Add java tests
    AddJavaTest('FilenameTest')
    AddJavaTest('ImageGenericTest')
    AddJavaTest('MetadataTest')
    #env.Default('run_java_tests')
    # Compile the HandleExtraFileTypes for ImageJ I/O
    compileHEFT = CompileJava('HandleExtraFileTypes', ['XmippViewer'])
    copyHandleExtraFileTypes = env.Install('external/imagej/plugins/', javaBuildName('HandleExtraFileTypes.class'))
    env.Default(copyHandleExtraFileTypes)
    # with the alias sometimes does not perform the copy?
    # env.Alias('java', copyHandleExtraFileTypes) # Add copy of this class as a Java dependency


# --- Programs

# Src (apps)

if not int(env['release']):
    AddXmippProgram('angular_commonline', ['XmippRecons'])
AddXmippProgram('angular_continuous_assign', ['XmippRecons'])
AddXmippProgram('angular_discrete_assign', ['XmippRecons'])
AddXmippProgram('angular_distance', ['XmippRecons'])
AddXmippProgram('angular_distribution_show', ['XmippInterface'])
AddXmippProgram('angular_neighbourhood', ['XmippRecons'])
AddXmippProgram('angular_projection_matching', ['XmippRecons'])
AddXmippProgram('angular_project_library', ['XmippRecons'])
AddXmippProgram('angular_rotate')
AddXmippProgram('classify_analyze_cluster', ['XmippClassif'])
AddXmippProgram('classify_compare_classes', ['XmippRecons'])
AddXmippProgram('classify_evaluate_classes', ['XmippRecons'])
AddXmippProgram('classify_kerdensom', ['XmippClassif'])
AddXmippProgram('ctf_correct_wiener3d', ['XmippRecons'])
AddXmippProgram('ctf_correct_idr', ['XmippRecons'])
AddXmippProgram('ctf_create_ctfdat', ['XmippRecons'])
AddXmippProgram('ctf_enhance_psd', ['XmippRecons'])
AddXmippProgram('ctf_estimate_from_micrograph', ['XmippRecons'])
AddXmippProgram('ctf_estimate_from_psd', ['XmippRecons'])
AddXmippProgram('ctf_group', ['XmippRecons'])
AddXmippProgram('ctf_phase_flip', ['XmippRecons'])
AddXmippProgram('ctf_show', ['XmippRecons'])
AddXmippProgram('ctf_sort_psds', ['XmippRecons'])
#AddXmippProgram('denoise', ['XmippRecons'])
if not int(env['release']):
    AddXmippProgram('idr_xray_tomo', ['XmippRecons'])
AddXmippProgram('image_align', ['XmippRecons'])
AddXmippProgram('image_align_tilt_pairs', ['XmippRecons'])
AddXmippProgram('image_common_lines', ['XmippRecons'])
AddXmippProgram('image_convert')
AddXmippProgram('image_find_center')
AddXmippProgram('image_header')
AddXmippProgram('image_histogram')
AddXmippProgram('image_operate')
AddXmippProgram('image_rotational_pca', ['XmippRecons'])
AddXmippProgram('image_resize', ['XmippRecons'])
AddXmippProgram('image_rotational_spectra', ['XmippRecons'])
AddXmippProgram('image_sort_by_statistics', ['XmippRecons'])
AddXmippProgram('image_separate_objects')
AddXmippProgram('image_statistics')
AddXmippProgram('image_vectorize')
AddXmippProgram('matrix_dimred', ['XmippDimred'])
#AddXmippProgram('mean_shift')
AddXmippProgram('metadata_convert_to_spider', ['XmippInterface'])
AddXmippProgram('metadata_histogram')
AddXmippProgram('metadata_import')
AddXmippProgram('metadata_split')
AddXmippProgram('metadata_utilities')
AddXmippProgram('metadata_xml')
AddXmippProgram('micrograph_scissor'),
AddXmippProgram('micrograph_automatic_picking', ['XmippRecons'])
AddXmippProgram('ml_align2d', ['XmippRecons'])
AddXmippProgram('mlf_align2d', ['XmippRecons'])
AddXmippProgram('ml_refine3d', ['XmippRecons'])
AddXmippProgram('mlf_refine3d', ['XmippRecons'])
AddXmippProgram('ml_tomo', ['XmippRecons'])
AddXmippProgram('mrc_create_metadata')
AddXmippProgram('nma_alignment', ['XmippRecons'])
AddXmippProgram('flexible_alignment', ['XmippRecons'])
AddXmippProgram('pdb_nma_deform', ['XmippRecons'])
AddXmippProgram('pdb_analysis', ['XmippRecons'])
AddXmippProgram('phantom_create', ['XmippRecons'])
AddXmippProgram('phantom_project', ['XmippRecons', 'XmippInterface'])
AddXmippProgram('phantom_simulate_microscope', ['XmippRecons'])
AddXmippProgram('phantom_transform', ['XmippRecons', 'XmippInterface'])
AddXmippProgram('reconstruct_art', ['XmippRecons'])
AddXmippProgram('reconstruct_art_pseudo', ['XmippRecons'])
if not int(env['release']):
    AddXmippProgram('reconstruct_art_xray', ['XmippRecons'])
AddXmippProgram('reconstruct_wbp', ['XmippRecons'])
AddXmippProgram('reconstruct_fourier', ['XmippRecons'])
AddXmippProgram('resolution_fsc')
if not int(env['release']):
    AddXmippProgram('resolution_ibw', ['XmippRecons'])
AddXmippProgram('resolution_ssnr', ['XmippRecons'])
AddXmippProgram('transform_add_noise')
AddXmippProgram('transform_adjust_volume_grey_levels', ['XmippRecons'])
AddXmippProgram('transform_center_image')
AddXmippProgram('transform_dimred', ['XmippDimred'])
AddXmippProgram('transform_downsample', ['XmippRecons'])
AddXmippProgram('transform_filter', ['XmippRecons'])
AddXmippProgram('transform_geometry')
AddXmippProgram('transform_mask')
AddXmippProgram('transform_mirror')
AddXmippProgram('transform_morphology')
AddXmippProgram('transform_normalize')
AddXmippProgram('transform_randomize_phases')
AddXmippProgram('transform_range_adjust')
AddXmippProgram('transform_symmetrize', ['XmippRecons'])
AddXmippProgram('transform_threshold', ['XmippRecons'])
AddXmippProgram('transform_window')
#AddXmippProgram('fourier_projection', ['XmippRecons'])
#AddXmippProgram('test_sql')
#AddXmippProgram('template_threads')
if not int(env['release']):
	AddXmippProgram('tomo_align_dual_tilt_series', ['XmippRecons'])
	AddXmippProgram('tomo_align_refinement', ['XmippRecons'])
AddXmippProgram('tomo_align_tilt_series', ['XmippRecons' ], useCudaEnvironment=int(env['cuda']))
AddXmippProgram('tomo_detect_missing_wedge', ['XmippRecons'])
AddXmippProgram('tomo_project', ['XmippRecons'])
AddXmippProgram('tomo_remove_fluctuations', ['XmippRecons'])
AddXmippProgram('tomo_extract_subvolume', ['XmippRecons'])
AddXmippProgram('volume_align')
AddXmippProgram('volume_center')
AddXmippProgram('volume_correct_bfactor', ['XmippRecons'])
AddXmippProgram('volume_enhance_contrast', ['XmippRecons'])
AddXmippProgram('volume_find_symmetry')
AddXmippProgram('volume_from_pdb', ['XmippRecons'])
AddXmippProgram('volume_reslice')
AddXmippProgram('volume_segment', ['XmippRecons'])
AddXmippProgram('volume_structure_factor')
AddXmippProgram('volume_to_pseudoatoms', ['XmippRecons'])
AddXmippProgram('volume_to_web')
AddXmippProgram('xray_import', ['XmippRecons'])
AddXmippProgram('xray_psf_create')
AddXmippProgram('xray_project', ['XmippRecons'])

#if not int(env['release']):
	#AddXmippProgram('xray_volume_correct', ['XmippRecons'])	

if int(env['arpack']):
    AddXmippProgram('angular_gcar', ['XmippRecons'])


# --- Scripts

# Python Batches (apps)
#
AddBatch('apropos', 'applications/scripts/apropos', '.py')
AddBatch('compile', 'applications/scripts/compile', '.py')
AddBatch('export_emx', 'applications/scripts/export_emx', '.py')
AddBatch('import_box', 'applications/scripts/import_box', '.py')
AddBatch('import_ctfparam', 'applications/scripts/import_ctfparam', '.py')
AddBatch('import_ctfdat', 'applications/scripts/import_ctfdat', '.py')
AddBatch('import_emx', 'applications/scripts/import_emx', '.py')
#AddBatch('metadata_operate', 'applications/scripts/metadata_operate','.py')
AddBatch('metadata_plot', 'applications/scripts/metadata_plot', '.py')
AddBatch('metadata_selfile_create', 'applications/scripts/metadata_selfile_create', '.py')
protocols_main = AddBatch('protocols', 'protocols', '.py')
env.Alias('protocols', protocols_main)
AddBatch('browser', 'applications/scripts/browser', '.py')
#AddBatch('browserj', 'applications/scripts/browserj', '.py')
AddBatch('micrograph_particle_picking', 'applications/scripts/micrograph_particle_picking', '.py')
AddBatch('chimera_client', 'applications/scripts/chimera_client', '.py')
#AddBatch('metadata_showj', 'applications/scripts/metadata_showj', '.py')
AddBatch('micrograph_tiltpair_picking', 'applications/scripts/micrograph_tiltpair_picking', '.py')
AddBatch('mpi_steps_runner', 'protocols', '.py')
AddBatch('projections_explorerj', 'applications/scripts/projections_explorerj', '.py')
#AddBatch('rot_spectraj', 'applications/scripts/rot_spectraj', '.py')
AddBatch('showj', 'applications/scripts/showj', '.py')
#AddBatch('stitchingj', 'applications/scripts/stitchingj', '.py')
AddBatch('tomoj', 'applications/scripts/tomoj', '.py')
AddBatch('visualize_preprocessing_micrographj', 'applications/scripts/visualize_preprocessing_micrograph', '.py')

# Shell script files
#
SymLink('bin/xmipp_imagej', 'external/runImageJ')

# Shell script files
#
SymLink('bin/xmipp_imagej', 'external/runImageJ')

# MPI
AddXmippMPIProgram('mpi_angular_class_average', ['XmippRecons'])
AddXmippMPIProgram('mpi_angular_continuous_assign', ['XmippRecons'])
if not int(env['release']):
    AddXmippMPIProgram('mpi_angular_gcar_commonlines', ['XmippRecons'])
AddXmippMPIProgram('mpi_angular_projection_matching', ['XmippRecons'])
AddXmippMPIProgram('mpi_angular_project_library', ['XmippRecons'])
AddXmippMPIProgram('mpi_classify_CL2D', ['XmippRecons'])
AddProgramLink('classify_CL2D', 'mpi_classify_CL2D')
if not int(env['release']):
    AddXmippMPIProgram('mpi_classify_CL3D', ['XmippRecons'])
    AddProgramLink('classify_CL3D', 'mpi_classify_CL3D')
AddXmippMPIProgram('mpi_classify_CL2D_core_analysis', ['XmippRecons'])
if not int(env['release']):
    AddXmippMPIProgram('mpi_classify_FTTRI', ['XmippRecons'])
AddXmippMPIProgram('mpi_ctf_correct_idr', ['XmippRecons'])
AddXmippMPIProgram('mpi_ctf_sort_psds', ['XmippRecons'])
AddXmippMPIProgram('mpi_image_operate')
AddXmippMPIProgram('mpi_image_rotational_pca', ['XmippRecons'])
AddXmippMPIProgram('mpi_performance_test', ['XmippRecons'])
AddXmippMPIProgram('mpi_image_resize', ['XmippRecons'])
AddXmippMPIProgram('mpi_image_sort', ['XmippRecons'])
AddProgramLink('image_sort', 'mpi_image_sort')
AddXmippMPIProgram('mpi_ml_align2d', ['XmippRecons'])
AddXmippMPIProgram('mpi_ml_tomo', ['XmippRecons'])
AddXmippMPIProgram('mpi_mlf_align2d', ['XmippRecons'])
AddXmippMPIProgram('mpi_ml_refine3d', ['XmippRecons'])
AddXmippMPIProgram('mpi_mlf_refine3d', ['XmippRecons'])
AddXmippMPIProgram('mpi_nma_alignment', ['XmippRecons'])
AddXmippMPIProgram('mpi_xray_project', ['XmippRecons'])
AddXmippMPIProgram('mpi_reconstruct_art', ['XmippRecons'])
AddXmippMPIProgram('mpi_reconstruct_wbp', ['XmippRecons'])
AddXmippMPIProgram('mpi_reconstruct_fourier', ['XmippRecons'])
AddXmippMPIProgram('mpi_run', ['XmippRecons'])
AddXmippMPIProgram('mpi_tomo_extract_subvolume', ['XmippRecons'])
AddXmippMPIProgram('mpi_transform_filter', ['XmippRecons'])
AddXmippMPIProgram('mpi_transform_symmetrize', ['XmippRecons'])
AddXmippMPIProgram('mpi_transform_geometry', ['XmippRecons'])
AddXmippMPIProgram('mpi_transform_mask', ['XmippRecons'])
AddXmippMPIProgram('mpi_transform_normalize', ['XmippRecons'])
if not int(env['release']):
    AddXmippMPIProgram('mpi_write_test', ['XmippRecons'])
#    AddXmippMPIProgram('template_threads', ['XmippRecons'])
#    AddXmippMPIProgram('template_mpi', ['XmippRecons'])

#---- Tests
if int(env['gtest']):
     AddXmippCTest('test_ctf')
     AddXmippCTest('test_euler')
     AddXmippCTest('test_fftw')
     AddXmippCTest('test_filters')
     AddXmippCTest('test_fringe_processing')          
     AddXmippCTest('test_funcs')
     AddXmippCTest('test_geometry')
     AddXmippCTest('test_image')
     AddXmippCTest('test_image_generic')
     AddXmippCTest('test_matrix')
     AddXmippCTest('test_metadata')
     AddXmippCTest('test_multidim')
     AddXmippCTest('test_polar')
     AddXmippCTest('test_polynomials')          
     AddXmippCTest('test_sampling')
     AddXmippCTest('test_symmetries')
     AddXmippCTest('test_transformation')
     AddXmippCTest('test_dimred')
     AddXmippCTest('test_wavelets')
     #env.Depends('run_tests', [fftw, tiff, sqlite])
     #python tests
     test = AddXmippPythonTest('test_pythoninterface')
     #AddXmippPythonTest('test_projectionmatching')
     AddXmippPythonTest('test_pysqlite')
     AddXmippPythonTest('test_emx')
     env.Depends(test, pythonbinding)
     env.Depends('run_tests', 'xmipp_programs')
     #env.Default('run_tests'     )

if int(env['matlab']):
    def CompileMatlab(name, dependencies=[]):
        ''' name parameter is expected without .java extension '''
        source = 'libraries/bindings/matlab/'+name+".cpp"
        target = 'libraries/bindings/matlab/'+name+".mexa64"
        command = env['MATLAB_DIR'] + '/bin/mex -O -outdir libraries/bindings/matlab -I. -Ilibraries -Llib -lXmippRecons -lXmippData -lXmippExternal '+source
        compileCmd = env.Command(target, source, command)
        env.Default(compileCmd)
        return compileCmd
    
    bindings = ['adjust_ctf', 'align2d', 'ctf_correct_phase',
        'mask', 'mirror', 'morphology', 'normalize', 'psd_enhance',
        'resolution', 'rotate', 'scale', 'scale_pyramid', 'volume_segment']
    for i in range(len(bindings)):
       CompileMatlab("tom_xmipp_"+bindings[i]+"_wrapper")
    CompileMatlab('xmipp_read')
    CompileMatlab('xmipp_nma_read_alignment')
    CompileMatlab('xmipp_nma_save_cluster')
    CompileMatlab('xmipp_read_structure_factor')

# Clean
# Configuration or cleaning
#if env.GetOption('clean'):
#    print '* Cleaning  ...'
#    os.system("( cd external/fftw-3.2.2    ; make clean >& /dev/null )");
#    os.system("( cd external/sqlite-3.6.23 ; make clean >& /dev/null )");
