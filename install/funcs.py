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

import platform
import os
import sys
import time


# Then we get some OS vars
MACOSX = (platform.system() == 'Darwin')
WINDOWS = (platform.system() == 'Windows')
LINUX = (platform.system() == 'Linux')

#SCIPION_URL_SOFTWARE = os.environ['SCIPION_URL_SOFTWARE']
SCIPION_URL_SOFTWARE = 'http://scipionwiki.cnb.csic.es/files/scipion/software'
SCIPION_PROCS = 5


def greenText(text):
    return '\x1b[%sm%s\x1b[0m' % ('32', text)

def printText(text, green=False):
    if not green:
        t = text
    else:
        t = greenText(text)
        
    print t
    

class Command():
    def __init__(self, env, cmd, targets=None,  **kwargs):
        self._env = env
        self._cmd = cmd
        
        if targets is None:
            self._targets = []
        elif isinstance(targets, list):
            self._targets = targets
        else:
            self._targets = [targets]
        
        self._environ = kwargs.get('environ', None)
        self._cwd = kwargs.get('cwd', None)
        self._out = kwargs.get('out', None)
        self._always = kwargs.get('always', False)
        
    def _existsAll(self):
        """ Return True if all targets exists. """
        for t in self._targets:
            if not os.path.exists(t):
                return False
        return True
    
    def execute(self):
        if not self._always and self._targets and self._existsAll():
            printText("  Skipping command: %s" % self._cmd)
            printText("  All targets exists.")
        else:         
            cwd = os.getcwd()
            if self._cwd is not None:
                if self._env.doExecute():
                    os.chdir(self._cwd)
                printText("  cd %s" % self._cwd)
                
            cmd = self._cmd
            if self._out is not None:
                cmd += ' 1> %s 2>&1' % self._out
            printText("  %s " % cmd)
            if self._env.doExecute():
                os.system(cmd)
            # Return to working directory, useful
            # when changing dir before executing command
            os.chdir(cwd)
            
    def __str__(self):
        return "Command: %s, targets: %s" % (self._cmd, self._targets)
        
        
class Target():
    def __init__(self, env, name, *commands, **kwargs):
        self._env = env
        self._name = name
        self._default = kwargs.get('default', False)
        self._commandList = []
        for c in commands:
            self._commandList.append(c)
        
    def addCommand(self, *args, **kwargs):
        c = Command(self._env, *args, **kwargs)
        self._commandList.append(c)
        return c
        
    def _existsAll(self):
        for command in self._commandList:
            #print ">>> Checking %s" % command
            if not command._existsAll():
                #print "   not _existsAll()"
                return False
        return True
    
    def isDefault(self):
        return self._default
    
    def getName(self):
        return self._name
        
    def execute(self):
        t1 = time.time()
        
        printText("> Building '%s'" % self._name, green=True)
        if self._existsAll():
            printText("  All targets exists, skipping.")
        else:
            for command in self._commandList:
                command.execute()
                
        t2 = time.time()
        
        if self._env.doExecute():
            e = int(t2 - t1)
            mins = e / 60
            secs = e % 60
            if mins > 0:
                minStr = '%d min ' % mins
            else:
                minStr = ''
            
            printText('< Elapsed: %s%d secs' % (minStr, secs))
        
        
class Environment():
    
    def __init__(self, **kwargs):
        self._targetList = []
        self._args = kwargs.get('args', [])
        self._doExecute = not '--show' in self._args
        
        if LINUX:
            self._libSuffix = 'so' # Shared libraries extension name
        else:
            self._libSuffix = 'dylib'
            
        self._downloadCmd = 'wget -nv -c -O %s %s'
        self._tarCmd = 'tar --recursive-unlink -xzf %s'
        
        self.AddLibrary = self.addLibrary
    
    def doExecute(self):
        return self._doExecute
    
    def getLib(self, name):
        return 'software/lib/lib%s.%s' % (name, self._libSuffix)
    
    def getBin(self, name):
        return 'software/bin/%s' % name
    
    def createTarget(self, name, *commands, **kwargs):
        t = Target(self, name, *commands, **kwargs)
        self._targetList.append(t)        
        return t
        
    def addLibrary(self, name, **kwargs): 
                   
#                    tar=None, buildDir=None, configDir=None, 
#                    targets=[], makeTargets=None, libChecks=[], url=None, flags=[], addPath=True,
#                    autoConfigTargets='Makefile', deps=[], clean=[], default=True):
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
        tar = kwargs.get('tar', '%s.tgz' % name)
        url = kwargs.get('url', '%s/external/%s' % (SCIPION_URL_SOFTWARE, tar))
        buildDir = kwargs.get('buildDir', 
                              tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0])
        
        configDir = kwargs.get('configDir', buildDir)
        configTarget = kwargs.get('configTarget', 'Makefile')
        configAlways = kwargs.get('configAlways', False)
        
        flags = kwargs.get('flags', [])
        targets = kwargs.get('targets', [self.getLib(name)])    
        clean = kwargs.get('clean', False) # Execute make clean at the end??
        cmake = kwargs.get('cmake', False) # Use cmake instead of configure??
        
                
        # Download library tgz
        tarFile = 'software/tmp/%s' % tar
        buildPath = 'software/tmp/%s' % buildDir
        configPath = 'software/tmp/%s' % configDir
        makeFile = '%s/%s' % (configPath, configTarget)
        
        t = self.createTarget(name, default=kwargs.get('default', True))
        
        t.addCommand(self._downloadCmd % (tarFile, url),
                     targets=tarFile)
        t.addCommand(self._tarCmd % tar,
                     targets=buildPath,
                     cwd='software/tmp')
        
        prefixPath = os.path.abspath('software')
        
        if not cmake:
            flags.append('--prefix=%s' % prefixPath)
            flags.append('--libdir=%s/lib' % prefixPath)
        
            t.addCommand('./configure %s' % ' '.join(flags),
                         targets=makeFile,
                         cwd=configPath, 
                         out='%s/log/%s_configure.log' % (prefixPath, name),
                         always=configAlways)
        else:
            flags.append('-DCMAKE_INSTALL_PREFIX:PATH=%s .' % prefixPath)
            t.addCommand('cmake %s' % ' '.join(flags),
                         targets=makeFile,
                         cwd=configPath, 
                         out='%s/log/%s_cmake.log' % (prefixPath, name))
        
        t.addCommand('make -j %d' % SCIPION_PROCS,
                     cwd=buildPath,
                     out='%s/log/%s_make.log' % (prefixPath, name))
        
        t.addCommand('make install',
                     targets=targets,
                     cwd=buildPath, 
                     out='%s/log/%s_make_install.log' % (prefixPath, name))
        
        if clean:
            t.addCommand('make clean',
                         cwd=buildPath, 
                         out='%s/log/%s_make_clean.log' % (prefixPath, name))
            t.addCommand('rm %s' % makeFile)
            
    def addModule(self, name, **kwargs): 
                  #tar=None, buildDir=None, targets=None, libChecks=[],
                  #url=None, flags=[], deps=[], clean=[], default=True):
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
        tar = kwargs.get('tar', '%s.tgz' % name)
        url = kwargs.get('url', '%s/python/%s' % (SCIPION_URL_SOFTWARE, tar))
        buildDir = kwargs.get('buildDir', 
                              tar.rsplit('.tar.gz', 1)[0].rsplit('.tgz', 1)[0])
        
        targets = kwargs.get('targets', [name])
        flags = kwargs.get('flags', [])
        
        tarFile = 'software/tmp/%s' % tar
        buildPath = 'software/tmp/%s' % buildDir
        
        prefixPath = os.path.abspath('software')
        flags.append('--prefix=%s' % prefixPath)
    
        t = self.createTarget(name, default=kwargs.get('default', True))
        
        t.addCommand(self._downloadCmd % (tarFile, url),
                     targets=tarFile)
        t.addCommand(self._tarCmd % tar,
                     targets='%s/setup.py' % buildPath,
                     cwd='software/tmp')

        t.addCommand('PYTHONHOME="%(root)s" LD_LIBRARY_PATH="%(root)s/lib" '
                     'PATH="%(root)s/bin:%(PATH)s" '
    #               'CFLAGS="-I%(root)s/include" LDFLAGS="-L%(root)s/lib" '
    # The CFLAGS line is commented out because even if it is needed for modules
    # like libxml2, it causes problems for others like numpy and scipy (see for
    # example http://mail.scipy.org/pipermail/scipy-user/2007-January/010773.html)
    # TODO: maybe add an argument to the function to chose if we want them?
                     '%(root)s/bin/python setup.py install %(flags)s > '
                     '%(root)s/log/%(name)s.log 2>&1' % {'root': prefixPath,
                                                       'PATH': os.environ['PATH'],
                                                       'flags': ' '.join(flags),
                                                       'name': name},
                   targets=['software/lib/python2.7/site-packages/%s' % tg for tg in targets],
                   cwd=buildPath)
        
    def execute(self):
        for target in self._targetList:
            if (target.isDefault() or 
                target.getName() in self._args): 
                target.execute()
            
        sys.exit(0)