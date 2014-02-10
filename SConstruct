#!/usr/bin/env python

# basic setup, import all environment and custom tools
import os
import platform 
import SCons.Script
CONFIG = '.xmipp_scons.options'
BASHRC = '.xmipp.bashrc'
CSH = '.xmipp.csh'
MPI_LIBDIR = ''
MPI_BINDIR = ''

if platform.system() == 'Windows':
    env = Environment(tools = ['mingw'], ENV = os.environ)
    env['ENV']['JAVA_HOME'] = "/c/Java/jdk1.6.0_34"
    env.PrependENVPath('PATH', 'C:\\MinGW\\bin')
    env.PrependENVPath('LIB', 'C:\\MinGW\\lib') 
else:
    env = Environment(ENV=os.environ,
          tools=['default', 'disttar'],
	  toolpath=['external/scons/ToolsFromWiki'])
    env.AppendUnique(LIBPATH=os.environ['LD_LIBRARY_PATH'])
    env.AppendUnique(LIBPATH=['/usr/lib64/openmpi/lib','/usr/lib64/mpi/gcc/openmpi/lib64','/usr/lib/openmpi'])

# avoid cruft in top dir
base_dir = 'build'
if not os.path.exists(base_dir):
    Execute(Mkdir(base_dir))
base_dir += '/'

# use only one signature file
env.SConsignFile(base_dir + 'SCons.dblite')

# read -or not- the cached -non default- options
if (ARGUMENTS['mode'] == 'configure'):
    opts = Variables(None, ARGUMENTS)
else:
    opts = Variables('.xmipp_scons.options', ARGUMENTS)

if (ARGUMENTS['mode'] == 'dependencies'):
    conf = Configure(env)
    checking = {}
    if not conf.CheckLib('mpi', None, None, 'cxx'):
        checking['mpi'] = "not found"
    if not conf.CheckLib('freetype', None, None, 'cxx'):
        checking['freetype'] = "not found"
    if not conf.CheckLib('X11', None, None, 'cxx'):
        checking['X11'] = "not found"
    if not conf.CheckLib('png', None, None, 'cxx'): 
        checking['png'] = "not found"
    if not conf.CheckLib('ncurses', None, None, 'cxx'): 
        checking['ncurses'] = "not found"
    if not conf.CheckLib('ssl', None, None, 'cxx'): 
        checking['ssl'] = "not found"
    if not conf.CheckLib('readline', None, None, 'cxx'): 
        checking['readline'] = "not found"
    if checking == {}:
        print 'All dependencies satisfied, proceeding with compilation'
    else:
        print 'Some dependencies unsatisfied, please check the following list and install them all:'
        for k, v in checking.items():
            print u'{0}: {1}'.format(k, v)
        ans = ""
        if 'unattended' in ARGUMENTS:
            if ARGUMENTS['unattended'] == 'yes':
                print "Unattended compilation selected, proceeding with the compilation."
        else:
            while ans != ("y" or "Y" or "n" or "N"):
                ans = raw_input("Do you still want to proceed with the compilation? (y/n):")
                if ans == "n" or ans == "N":
                    import os
                    print "You've choosen to abort the installation. Aborting..."
                    os._exit(os.EX_OSERR)
                if ans != "y" and ans != "Y":
                    import os
                    print "Unknown option, please answer y/Y/n/N"
    env = conf.Finish()


#print 'en opts', opts['MPI_LINKERFORPROGRAMS']

opts.Add('CC', 'The C compiler', 'gcc')
opts.Add('CXX', 'The C++ compiler', 'g++')

# Hack, some architectures required this
opts.Add('LINKERFORPROGRAMS', 'Linker for programs', 'g++')

if platform.system()=='Windows':
    opts.Add('CCFLAGS', 'The C compiler flags', '-fpermissive -I/c/MinGW/include')
    opts.Add('CXXFLAGS', 'The C++ compiler flags', '-fpermissive -I/c/MinGW/include')
    opts.Add(BoolVariable('release', 'Release mode', 'yes'))
else:
    opts.Add('CCFLAGS', 'The C compiler flags', '-std=c99')
    opts.Add('CXXFLAGS', 'The C++ compiler flags', None)
    opts.Add(BoolVariable('release', 'Release mode', 'yes'))

opts.Add(BoolVariable('debug', 'Build debug version?', 'no'))
#Profile version implies debug and then it will be ignored
opts.Add(BoolVariable('profile', 'Build profile version?', 'no'))
opts.Add(BoolVariable('warn', 'Show warnings?', 'no'))
opts.Add(BoolVariable('fast', 'Fast?', 'no'))
opts.Add(BoolVariable('static', 'Prevent dynamic linking?', 'no'))

opts.Add('prepend', 'What to prepend to executable names', 'xmipp')
opts.Add(BoolVariable('quiet', 'Hide command line?', 'yes'))

opts.Add(BoolVariable('java', 'Build the java programs?', 'yes'))

opts.Add(BoolVariable('gtest', 'Build tests?', 'yes'))

#opts.Add(BoolVariable('mpi', 'Build the MPI programs?', 'yes'))

opts.Add('MPI_CC', 'MPI C compiler', 'mpicc')
opts.Add('MPI_CXX', 'MPI C++ compiler', 'mpiCC')
opts.Add('MPI_LINKERFORPROGRAMS', 'MPI Linker for programs', 'mpiCC')
opts.Add('MPI_INCLUDE', 'MPI headers dir ', '/usr/include')
opts.Add('MPI_LIBDIR', 'MPI libraries dir ', '/usr/lib')
opts.Add('MPI_LIB', 'MPI library', 'mpi')
opts.Add('MPI_BINDIR', 'MPI binaries', '/usr/bin')

#MINGW 
opts.Add('MINGW_PATHS', 'Include path for MinGW', '')

opts.Add('prefix', 'Base installation directory', Dir('.').abspath)

opts.Add(BoolVariable('matlab', 'Build the Matlab bindings?', 'no'))
opts.Add('MATLAB_DIR', 'Matlab installation dir', '/usr/local/MATLAB/R2011a')

opts.Add(BoolVariable('cuda', 'Build GPU stuff?', 'no'))
opts.Add('CUDA_SDK_PATH', 'CUDA SDK dir', '/root/NVIDIA_GPU_Computing_SDK')
opts.Add('CUDA_LIB_PATH', 'CUDA RunTimeLib dir', '/usr/local/cuda/lib64')


#opts.Add(BoolVariable('verbose_fftw', 'Verbose configuring of FFTW libraries?', 'no'))
#opts.Add('FFTWFLAGS', 'Additional flags for FFTW configure', '--enable-threads')
#
#opts.Add(BoolVariable('verbose_tiff', 'Verbose configuring of TIFF libraries?', 'no'))
#opts.Add('TIFFFLAGS', 'Additional flags for TIFF configure', 'CPPFLAGS=-w')
#
#opts.Add(BoolVariable('verbose_sqlite', 'Verbose configuring of SQLite libraries?', 'no'))
#opts.Add('SQLITEFLAGS', 'Additional flags for SQLite configure', 'CPPFLAGS=-w CFLAGS=-DSQLITE_ENABLE_UPDATE_DELETE_LIMIT=1')
#
#opts.Add(BoolVariable('verbose_python', 'Verbose configuring of Python compilation?', 'no'))
#opts.Add('PYTHONFLAGS', 'Additional flags for Python configure', '')

opts.Add('JAVAC', 'Java compiler', 'javac')
opts.Add('JAVA_HOME', 'Java installation directory', '')
opts.Add('JNI_CPPPATH', 'Directory of jni.h', '')
#print 'en opts-2', opts['MPI_LINKERFORPROGRAMS']

opts.Update(env)
# generate help for options
Help(opts.GenerateHelpText(env, sort=cmp))

# FIXME Hack, for several flags in command-line
env['CCFLAGS'] = Split(env['CCFLAGS'])
env['CXXFLAGS'] = Split(env['CXXFLAGS'])
env['JARFLAGS'] = '-Mcf'    # Default "cf". "M" = Do not add a manifest file.

# These defaults are needed for both the custom tests and the compilation
env.SetDefault(LIBS='')
env.SetDefault(LIBPATH='')
env.SetDefault(CPPPATH='')

def checkMpiInBashrc(mpiPath, type, shellType='bash'):
    fileToOpen = foundPath = stringToSearch = separator = ''
    exists = different = False
    if shellType == 'bash':
        fileToOpen = BASHRC
        separator = '='
    elif shellType == 'csh':
        fileToOpen = CSH
        separator = ' '
    if type == 'lib':
        stringToSearch = 'XMIPP_MPI_LIBDIR'
    elif type == 'bin':
        stringToSearch = 'XMIPP_MPI_BINDIR'
    if os.path.exists(fileToOpen):
        for line in open(fileToOpen):
            parts = line.split(separator)
            if len(parts) == 2:
                if parts[0].strip() == ('export '+stringToSearch):
                    exists = True
                    foundPath = parts[1]
                    different = (foundPath.strip() != mpiPath)
            elif len(parts) == 3:
                if parts[1].strip() == (stringToSearch):
                    exists = True
                    foundPath = parts[2]
                    different = (foundPath.strip() != mpiPath)
                    break #In CSH file, we have to find only the first occurrence             
    return exists, different, foundPath

def addMpiToBashrc(mpiPath, type, shellType='bash', replace=False):
    fileToOpen = preserv = stringToSearch = stringToAppend = ''
    separator = '='
    exportation = 'export'
    if shellType == 'bash':
        fileToOpen = BASHRC
        separator = '='
        exportation = 'export'
    elif shellType == 'csh':
        fileToOpen = CSH
        separator = ' '
        exportation = 'setenv'
    if os.path.exists(fileToOpen):
        lines = open(fileToOpen, 'r').readlines()
        filew = open(fileToOpen, 'w')
        found = -1
        if type == 'lib':
            stringToSearch = 'LD_LIBRARY_PATH'
            stringToAppend = 'XMIPP_MPI_LIBDIR'
        elif type == 'bin':
            stringToSearch = 'PATH'
            stringToAppend = 'XMIPP_MPI_BINDIR'
        if replace:
            stringToSearch = (exportation + ' ' + stringToAppend)
        else:
            stringToSearch = (exportation+ ' ' + stringToSearch)
        for index, line in enumerate(lines):
            parts = line.strip().split(separator)
            if len(parts) == 2:
                if parts[0].strip() == stringToSearch:
                    found = index
                    preserv = parts[1].strip()
            if len(parts) == 3:
                if (parts[0].strip()+' '+parts[1].strip()) == stringToSearch:
                    found = index
                    preserv = parts[2].strip()
                    break
        if found != -1:
            if not replace:
                lines[found] = (stringToSearch + separator + preserv + ':$' + stringToAppend + "\n")
                lines[found] = (exportation + ' ' + stringToAppend + separator + mpiPath + "\n") + lines[found]
            else:
                lines[found] = (stringToSearch+separator+mpiPath+"\n")
            filew.writelines(lines)


def AppendIfNotExists(**args):
    append = True
    for k, v in args.iteritems():
        if v in env[k]:
            append = False
    if append:
        env.Append(**args)
    
# This is required for both modes
env['STATIC_FLAG'] = '-static'

if (ARGUMENTS['mode'] == 'configure'):
    # --- This is the configure mode

    # Custom tests
    def CheckMPI(context, mpi_inc, mpi_libpath, mpi_lib, mpi_cc, mpi_cxx, mpi_link):
        context.Message('* Checking for MPI ... ')

        lastLIBS = context.env['LIBS']
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

        context.env.Replace(LIBS=lastLIBS)
        context.env.Replace(LIBPATH=lastLIBPATH)
        context.env.Replace(CPPPATH=lastCPPPATH)
        context.env.Replace(CC=lastCC)
        context.env.Replace(CXX=lastCXX)

        context.Result(ret)
        return ret

    # Configuration or cleaning
    if env.GetOption('clean'):
        print '* Cleaning  ...'
        if 'distclean' in COMMAND_LINE_TARGETS:
            print '* Deleting configuration ...'
            Execute(Delete(base_dir))
            Execute(Delete(env['prefix']))
    else:
        print '* Configuring  ...'
        config_dir = base_dir + 'config.tests'
        config_log = base_dir + 'config.log'

    # release?
    # This is done in compile mode
#    if int(env['release']):
#        AppendIfNotExists(CCFLAGS='-DRELEASE_MODE')
#        AppendIfNotExists(CXXFLAGS='-DRELEASE_MODE')

    # static?
    if int(env['static']):
        AppendIfNotExists(CCFLAGS='$STATIC_FLAG')
        AppendIfNotExists(LINKFLAGS='$STATIC_FLAG')

    # osx?
#    if env['PLATFORM'] == 'darwin':
#        AppendIfNotExists(CCFLAGS='-m64')
#        AppendIfNotExists(CXXFLAGS='-m64')
#        AppendIfNotExists(LINKFLAGS='-m64')

    # mingw?
    if platform.system() == 'Windows':
         AppendIfNotExists(CCFLAGS='-f permissive -Ilibraries/data -I/c/MinGW/include')
         AppendIfNotExists(CXXFLAGS='-f permissive -Ilibraries/data -I/c/MinGW/include')
#         AppendIfNotExists(LINKFLAGS='-pthread')

    # Non-GUI configuration environment
    conf = Configure(env, {'CheckMPI' : CheckMPI}, config_dir, config_log)

    # MPI
    if not conf.CheckMPI(env['MPI_INCLUDE'], env['MPI_LIBDIR'],
                             env['MPI_LIB'], env['MPI_CC'], env['MPI_CXX'], env['MPI_LINKERFORPROGRAMS']):
        import os
        print '* ERROR: Did not find MPI library. MPI support is mandatory.'
        os._exit(os.EX_OSERR)

    # Java
    import sys
    sys.path.append("external/scons/ToolsFromWiki")
    import ConfigureJNI
    if not ConfigureJNI.ConfigureJNI(env):
        import os
        print '* ERROR: Did not find JNI header. Java support is mandatory.'
        os._exit(os.EX_OSERR)
    else:
        print '* JNI header found: ', env['JNI_CPPPATH']

    # MATLAB
    if int(env['matlab']):
        print '* Checking for Matlab ... ',
        if not os.path.exists(env['MATLAB_DIR'] + '/bin/matlab'):
            print 'no'
            print '* Did not find Matlab. Disabling ...'
            env['matlab'] = 0
        else:
            print 'yes'

#    def ConfigureExternalLibrary(libName, libFlags, libpath, verbose):
#            
#        libpath = os.path.join('external', libpath)
#        if not GetOption("help"):
#            if os.path.exists(os.path.join(libpath, 'Makefile')):
#                command = 'cd %s ; make distclean > /dev/null' % libpath
#                print command
#                os.system(command)
#    
#        if not int(env['static']):
#            libFlags += " --enable-shared"
#    
#        if env['PLATFORM'] == 'darwin':
#            libFlags += " CFLAGS='-m64'"
#        
#        msg = "* Configuring " + libName
#        cmd = "cd %s;./configure %s" % (libpath, libFlags)
#        
#        logFile = os.path.join("..", "..", "build", libName + "_configure.log")
#        if not verbose:
#            msg += "(see log file '%s' for details)" % logFile
#            cmd += " > " + logFile
#            
#        print msg, "..."
#        print "   command:  ", cmd
#        os.system(cmd)        
    
#    ConfigureExternalLibrary('sqlite', env['SQLITEFLAGS'], 
#                             'sqlite-3.6.23', int(env['verbose_sqlite']))
#    ConfigureExternalLibrary('python', env['PYTHONFLAGS'], 
#                             'Python-2.7.2', int(env['verbose_python']))
#    ConfigureExternalLibrary('fftw', env['FFTWFLAGS'], 
#                             'fftw-3.2.2', int(env['verbose_fftw']))
#    ConfigureExternalLibrary('tiff', env['TIFFFLAGS'], 
#                             'tiff-3.9.4', int(env['verbose_tiff']))

    # Finish configuration
    env = conf.Finish()
    CONFIG = '.xmipp_scons.options'
    opts.Save(CONFIG, env)
    if os.path.exists(CONFIG):
        for line in open(CONFIG):
            parts = line.split('=')
            if len(parts) == 2:
                assign = '%s = "%s"' % (parts[0], parts[1])
                if parts[0].strip() == 'MPI_LIBDIR':
                    parts[1] = parts[1].replace("'", "").strip()
                    os.environ['LD_LIBRARY_PATH'] += os.pathsep + parts[1]
                    os.environ['MPI_LIBDIR'] = os.pathsep + parts[1]
                    MPI_LIBDIR = os.environ['MPI_LIBDIR']
                if parts[0].strip() == 'MPI_BINDIR':
                    parts[1] = parts[1].replace("'", "").strip()
                    os.environ['PATH'] += os.pathsep + parts[1]
                    os.environ['MPI_BINDIR'] = os.pathsep + parts[1]
                    MPI_BINDIR = os.environ['MPI_BINDIR']
    
    # When we are going to install, we need to put some vars in bashrc file
    if not 'unattended' in ARGUMENTS or (ARGUMENTS('unattended') != 'yes'):
        for shell in ['bash', 'csh']:
            try:
                mpiLibDirExists, mpiLibDirIsDifferent, mpiLibDirPath = checkMpiInBashrc(MPI_LIBDIR, 'lib', shell)
                if not mpiLibDirExists:
                    addMpiToBashrc(MPI_LIBDIR, 'lib', shell)
                elif mpiLibDirIsDifferent:
                    addMpiToBashrc(MPI_LIBDIR, 'lib', shell, True)
                else:
                    print "MPI_LIBDIR untouched. Already set in "+BASHRC+" file"
            except NameError:
                print "Be careful, your MPI_LIBDIR has not been set!"
            try:
                mpiBinDirExists, mpiBinDirIsDifferent, mpiBinDirPath = checkMpiInBashrc(MPI_BINDIR, 'bin', shell)
                if not mpiBinDirExists:
                    addMpiToBashrc(MPI_BINDIR, 'bin', shell)
                elif mpiBinDirIsDifferent:
                    addMpiToBashrc(MPI_BINDIR, 'bin', shell, True)
                else:
                    print "MPI_BINDIR untouched. Already set in "+BASHRC+" file"
            except NameError:
                print "Be careful, your MPI_BINDIR has not been set! You may need to manually set it in your bashrc"

elif (ARGUMENTS['mode'] == 'compile'):
    # --- This is the compilation mode

    # add separator to prepend (not shown in help)
    env['prepend'] = env['prepend'] + '_'

    dp = False #debug or profile?
    dpFlags = ""
    if int(env['profile']):
        dp = True
        dpFlags = ['-g', '-pg']
    elif int(env['debug']):
        dp = True
        dpFlags = ['-g']
    # profile or debug?
    if dp:
        AppendIfNotExists(CXXFLAGS=['-D__XMIPP_DEBUG__'])
        AppendIfNotExists(CXXFLAGS=dpFlags)
        AppendIfNotExists(CCFLAGS=dpFlags)
        AppendIfNotExists(LINKFLAGS=dpFlags)
        #required for print stack trace
    
    # Activate release version when compiling
    if int(env['release']):
        AppendIfNotExists(CCFLAGS=['-DRELEASE_MODE'])
        AppendIfNotExists(CXXFLAGS=['-DRELEASE_MODE'])
        
    if not int(env['cuda']):
	if env['PLATFORM'] != 'cygwin' and env['PLATFORM'] != 'win32':
            AppendIfNotExists(CXXFLAGS=['-rdynamic'])
        AppendIfNotExists(CXXFLAGS=['-O0'])
    else:
        # "safe" optimization
        AppendIfNotExists(CXXFLAGS=['-O2'])
        AppendIfNotExists(CCFLAGS=['-O2'])
        AppendIfNotExists(LINKFLAGS=['-s'])

    if env['PLATFORM'] == 'darwin':
#        env.Append(CXXFLAGS=['-m64'])
#        env.Append(CCFLAGS=['-m64'])
#        env.Append(LINKFLAGS=['-m64'])
        AppendIfNotExists(CXXFLAGS=['-I/usr/include/malloc'])
        AppendIfNotExists(CCFLAGS=['-I/usr/include/malloc'])
#        AppendIfNotExists(LINKFLAGS=['-m64'])


# stdout_handle = os.popen('svnversion -n .')
#    svnver = stdout_handle.read()
#    env.Append(CXXFLAGS="-D'SVN_REV=\""+ svnver +"\"'")
#    env.Append(CCFLAGS="-D'SVN_REV=\""+ svnver +"\"'")

    # Add threads
    env.Append(LINKFLAGS=['-lpthread'])
    #env.Append(CXXFLAGS=['-lpthread'])
    #env.Append(CCFLAGS=['-lpthread'])

    # warnings?
    if int(env['warn']) or int(env['debug']):
        env.Append(CXXFLAGS=['-Wall','-pedantic','-Wno-variadic-macros','-Wno-long-long','-Wno-deprecated'])
    else:
        env.Append(CXXFLAGS=['-w'])
        # TODO suppress linker warnings too... what's the flag?
        # env.Append(LINKFLAGS = ['-Wl,???'])

    # fast?
    # COSS' work. I dont like this. Classic debug vs release (asolano)
    if int(env['fast']):
        env.Append(CXXFLAGS=['-O3', '-fomit-frame-pointer', '-ffast-math',
                   '-finline-functions', '-funroll-loops'])

    # verbosity (for $SCONS_VERSION < 0.96.2 option has no effect)
    if int(env['quiet']):
        env['CCCOMSTR'] = 'Compiling $SOURCE'
        env['CXXCOMSTR'] = 'Compiling $SOURCE'
        env['SHCXXCOMSTR'] = 'Compiling $SOURCE'
        env['SHCCCOMSTR'] = 'Compiling $SOURCE'
        env['LINKCOMSTR'] = 'Linking $SOURCE'
        env['ARCOMSTR'] = 'Archiving $TARGET'
        env['SHLINKCOMSTR'] = 'Linking $TARGET'
        env['RANLIBCOMSTR'] = 'Indexing $TARGET'
        env['TARCOMSTR'] = 'Archiving $TARGET'
        env['INSTALLSTR'] = 'Installing $TARGET'

    Export('env')
    env.SConscript('SConscript')

elif (ARGUMENTS['mode'] == 'docs'):
    action = env.Action("doxygen")
    env.Execute(action)

# TODO: make specific modes for generation of dist

# distribution
""" FIXME Testing, not ready for production
env['DISTTAR_FORMAT'] = 'bz2'
env.Append(DISTTAR_EXCLUDEEXTS = ['.o', '.os', '.so', '.a', '.dll', '.cc',
    '.cache', '.pyc', '.cvsignore', '.dblite', '.log', '.bz2'],
DISTTAR_EXCLUDEDIRS = ['CVS', '.svn', '.sconf_temp', 'dist']
)
disttar = env.DistTar(os.path.join(base_dir, 'Xmipp-1.1'), [env.Dir('#')])
env.Alias('dist', disttar)
"""

