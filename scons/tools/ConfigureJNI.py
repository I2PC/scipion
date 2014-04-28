import os
import sys
from SCons.Script import *
import SCons.Builder

def walkDirs(path):
    """helper function to get a list of all subdirectories"""
    def addDirs(pathlist, dirname, names):
        """internal function to pass to os.path.walk"""
        for n in names:
            f = os.path.join(dirname, n)
            if os.path.isdir(f):
                pathlist.append(f)
    pathlist = [path]
    os.path.walk(path, addDirs, pathlist)
    return pathlist

def ConfigureJNI(env):
    """Configure the given environment for compiling Java Native Interface
       c or c++ language files."""

    if not env.get('JAVAC'):
        print "The Java compiler must be installed and in the current path."
        return 0

    # first look for a shell variable called JAVA_HOME
    java_base = os.environ.get('JAVA_HOME')
    if not java_base:
        if sys.platform == 'darwin':
            # Apple's OS X has its own special java base directory
            java_base = '/System/Library/Frameworks/JavaVM.framework'
        else:
            # Search for the java compiler
            print "JAVA_HOME environment variable is not set. Searching for java... ",
            tmpJavac=env.WhereIs('javac')
	    if tmpJavac:
	        jcdir = os.path.dirname(tmpJavac)
            else:
	        jcdir=None
	    if not jcdir:
                print "not found."
                return 0
            # assuming the compiler found is in some directory like
            # /usr/jdkX.X/bin/javac, java's home directory is /usr/jdkX.X
            java_base = os.path.join(jcdir, "..")
            print "found."

    if sys.platform == 'cygwin':
        # Cygwin and Sun Java have different ideas of how path names
        # are defined. Use cygpath to convert the windows path to
        # a cygwin path. i.e. C:\jdkX.X to /cygdrive/c/jdkX.X
        java_base = os.popen("cygpath -up '"+java_base+"'").read().replace( \
                 '\n', '')

    if sys.platform == 'darwin':
        # Apple does not use Sun's naming convention
        # java_headers = [os.path.join(java_base, 'Headers')]
        # But it has a Home/include dir to follow this conventions too
        java_headers = [os.path.join(java_base, 'include')]
        # java_libs = [os.path.join(java_base, 'Libraries')]
        java_libs = [os.path.join(java_base, 'lib')]
    else:
        # windows and linux
        java_headers = [os.path.join(java_base, 'include')]
        java_libs = [os.path.join(java_base, 'lib')]
        # Sun's windows and linux JDKs keep system-specific header
        # files in a sub-directory of include
        if java_base == '/usr' or java_base == '/usr/local':
            # too many possible subdirectories. Just use defaults
            java_headers.append(os.path.join(java_headers[0], 'win32'))
            java_headers.append(os.path.join(java_headers[0], 'linux'))
            java_headers.append(os.path.join(java_headers[0], 'solaris'))
        else:
            # add all subdirs of 'include'. The system specific headers
            # should be in there somewhere
            java_headers = walkDirs(java_headers[0])

    # Add extra potentially useful environment variables
    env['JAVA_HOME'] = java_base
    env['JNI_CPPPATH'] = java_headers
    env['JNI_LIBPATH'] = java_libs
    return 1
