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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This modules serve to define some Configuration classes
mainly for project GUI
"""

import os
from os.path import join, exists

import pyworkflow as pw
from pyworkflow.object import *
from pyworkflow.mapper import SqliteMapper, XmlMapper

PATH = os.path.dirname(__file__)
SETTINGS = join(pw.HOME, 'settings')

def getConfigPath(filename):
    """Return a configuration filename from settings folder"""
    return join(SETTINGS, filename)


class ProjectConfig(OrderedObject):
    """A simple base class to store ordered parameters"""
    def __init__(self, menu=None, protocols=None, **args):
        OrderedObject.__init__(self, **args)
        self.menu = String(menu)
        self.protocols = String(protocols)
        self.icon = String('scipion_bn.xbm')
        self.logo = String('scipion_logo.gif')


class MenuConfig(OrderedObject):
    """Menu configuration in a tree fashion.
    Each menu can contains submenus.
    Leaf elements can contain actions"""
    def __init__(self, text=None, value=None, 
                 icon=None, tag=None, **args):
        """Constructor for the Menu config item.
        Arguments:
          text: text to be displayed
          value: internal value associated with the item.
          icon: display an icon with the item
          tag: put some tags to items
        **args: pass other options to base class.
        """
        OrderedObject.__init__(self, **args)
        #List.__init__(self, **args)
        self.text = String(text)
        self.value = String(value)
        self.icon = String(icon)
        self.tag = String(tag)
        self.childCount = 0
        
    def _getSubMenuName(self):
        self.childCount += 1
        return '_child_%03d' % self.childCount
        
    def addSubMenu(self, text, value=None, **args):
        subMenu = type(self)(text, value, **args)
        name = self._getSubMenuName()
        setattr(self, name, subMenu)
        return subMenu
    
    def __iter__(self):
        for v in self._getChilds():
            yield v        
                
    def __setattr__(self, name, value):
        if len(name) == 0:
            name = self._getSubMenuName()
        OrderedObject.__setattr__(self, name, value)
                
    def _getStr(self, prefix):
        s = prefix + "%s text = %s, icon = %s\n" % (self.getClassName(), self.text.get(), self.icon.get())
        for sub in self:
            s += sub._getStr(prefix + "  ")
        return s
    
    def _getChilds(self):
        return [v for k, v in self.__dict__.iteritems() if k.startswith('_child_')]
    
    def __len__(self):
        return len(self._getChilds())
    
    def isEmpty(self):
        return len(self._getChilds()) == 0
            
    def __str__(self):
        return self._getStr(' ')
    
    
class ProtocolConfig(MenuConfig):
    """Store protocols configuration """
    pass
    
    
class ExecutionHostConfig(OrderedObject):
    """ Configuration information for execution hosts. """
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self.label = String()
        self.hostname = String()
        self.mpiCommand = String()
        self.queueSystem = QueueSystemConfig()
        
class QueueSystemConfig(OrderedObject):
    def __init__(self, **args):
        OrderedObject.__init__(self, **args) 
        self.name = String()
        self.mandatory = Boolean()
        self.queues = List() # List for queue configurations
        self.submitCommand = String()
        self.submitTemplate = String()
        self.checkCommand = String()
        self.cancelCommand = String()
        
    def hasValue(self):
        return self.name.hasValue() and len(self.queues)
        
        
class QueueConfig(OrderedObject):
    def __init__(self, **args):
        OrderedObject.__init__(self, **args) 
        self.name = String('default')
        self.maxCores = Integer()
        self.allowMPI = Boolean()
        self.allowThreads = Boolean()
        self.maxHours = Integer()
        
        
class ExecutionHostMapper(XmlMapper):
    def __init__(self, filename, dictClasses=None, **args):
        if dictClasses is None:
            dictClasses = globals()
        XmlMapper.__init__(self, filename, dictClasses, **args)
        #self.setClassTag('ExecutionHostConfig.QueueSystemConfig', 'class_name')
        self.setClassTag('ExecutionHostConfig.String', 'attribute')
        self.setClassTag('ExecutionHostConfig.mpiCommand', 'name_only')
        self.setClassTag('QueueSystemConfig.name', 'attribute')
        self.setClassTag('QueueSystemConfig.mandatory', 'attribute')  
        self.setClassTag('QueueConfig.ALL', 'attribute')   
        self.setClassTag('List.QueueConfig', 'class_only')           
        #self.setClassTag('ProtocolConfig.ProtocolConfig', 'class_only')
    
class ConfigMapper(XmlMapper):
    """Sub-class of XmlMapper to store configurations"""
    def __init__(self, filename, dictClasses, **args):
        XmlMapper.__init__(self, filename, dictClasses, **args)
        self.setClassTag('MenuConfig.MenuConfig', 'class_only')
        self.setClassTag('MenuConfig.String', 'attribute')
        self.setClassTag('ProtocolConfig.ProtocolConfig', 'class_only')
        self.setClassTag('ProtocolConfig.String', 'attribute')
        
    def getConfig(self):
        return self.selectAll()[0]


def writeConfig(config, fn, mapperClass=ConfigMapper):
    fn = getConfigPath(fn)
    print "config file: ", fn
    if exists(fn):
        os.remove(fn)
    mapper = mapperClass(fn, globals())
    mapper.insert(config)
    mapper.commit()
    
    
def writeMenus():
    """Write default configuration files"""
    # Write menu configuration
    menu = MenuConfig()
    projMenu = menu.addSubMenu('Project')
    projMenu.addSubMenu('Browse files', 'browse', icon='folderopen.gif')
    projMenu.addSubMenu('Remove temporary files', 'delete', icon='delete.gif')
    projMenu.addSubMenu('Clean project', 'clean')
    projMenu.addSubMenu('Exit', 'exit')
    
    helpMenu = menu.addSubMenu('Help')
    helpMenu.addSubMenu('Online help', 'online_help', icon='online_help.gif')
    helpMenu.addSubMenu('About', 'about')
    
    writeConfig(menu, 'menu_default.xml')
    
    # Write another test menu
    menu = MenuConfig()
    m1 = menu.addSubMenu('Test')
    m1.addSubMenu('KK', icon='tree.gif')
    m1.addSubMenu('PP', icon='folderopen.gif')
    writeConfig(menu, 'menu_test.xml')
    
def writeProtocols():
    """ Write protocols configuration. """
    menu = ProtocolConfig()
    m1 = menu.addSubMenu('Micrographs', tag='section')
    
    #m2 = m1.addSubMenu('Micrographs')
    m1.addSubMenu(' Import', value='ProtImportMicrographs', 
                  tag='protocol', icon='bookmark.png')
    m1.addSubMenu('Preprocess', value='ProtPreprocessMicrographs',
                  tag='protocol_base')
    m1.addSubMenu('CTF estimation', value='ProtCTFMicrographs',
                  tag='protocol_base')
    
    m1 = menu.addSubMenu('Particles', tag='section')
    m1.addSubMenu('Import', value='ProtImportParticles', 
                  tag='protocol', icon='bookmark.png')
    m1.addSubMenu('Picking', value='ProtParticlePicking',
                  tag='protocol_base')
    m1.addSubMenu('Extract', value='ProtExtractParticles',
                  tag='protocol_base')    

    m1 = menu.addSubMenu('2D', tag='section')
    
    m1.addSubMenu('Align', value='ProtAlign',
                  tag = 'protocol_base')
    m1.addSubMenu('Classify', value='ProtClassify',
                  tag = 'protocol_base')
    m1.addSubMenu('Align+Classify', value='ProtAlignClassify',
                  tag = 'protocol_base')

    writeConfig(menu, 'protocols_default.xml')
    
    
def writeHosts():
    host = ExecutionHostConfig()
    host.label.set('localhost')
    host.hostname.set('localhost')
    host.mpiCommand.set('mpirun -np %(nodes)d -bynode %(command)s')
    
    queueSys = host.queueSystem
    #queueSys = QueueSystemConfig()
    queueSys.name.set('PBS/TORQUE')
    queueSys.mandatory.set(False)
    queueSys.submitCommand.set('qsub %(script)s')
    queueSys.submitTemplate.set("""
#!/bin/bash
### Inherit all current environment variables
#PBS -V
### Job name
#PBS -N %(jobName)s
### Queue name
###PBS -q %(queueName)s
### Standard output and standard error messages
#PBS -k eo
### Specify the number of nodes and thread (ppn) for your job.
#PBS -l nodes=%(nodes)d:ppn=%(threads)d
### Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=%(hours)d:00:00
# Use as working dir the path where qsub was launched
WORKDIR=$PBS_O_WORKDIR
#################################
### Set environment varible to know running mode is non interactive
export XMIPP_IN_QUEUE=1
### Switch to the working directory;
cd $WORKDIR
# Make a copy of PBS_NODEFILE 
cp $PBS_NODEFILE %(nodesfileBackup)s
# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`
# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`
### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Working directory is `pwd`
echo Using ${NPROCS} processors across ${NNODES} nodes
echo PBS_NODEFILE:
cat $PBS_NODEFILE
#################################

%(command)s
""")
    queueSys.cancelCommand.set('canceljob %(jobid)d')
    queueSys.checkCommand.set('qstat %(jobid)d')
    
    queue = QueueConfig()
    queue.maxCores.set(4)
    queue.allowMPI.set(True)
    queue.allowThreads.set(True)
    
    queueSys.queues.append(queue)
    
    
    writeConfig(host, 'execution_hosts.xml', mapperClass=ExecutionHostMapper)
    
def writeDefaults():
    writeMenus()
    writeProtocols()
    writeHosts() 
    # Write global configuration
    config = ProjectConfig(menu='menu_default.xml',
                           protocols='protocols_default.xml')
    writeConfig(config, 'configuration.xml')
    
    
    
if __name__ == '__main__':
    #Write default configurations
    writeDefaults()
    fn = getConfigPath('execution_hosts.xml')
    mapper = ExecutionHostMapper(fn, globals())
    l = mapper.selectAll()[0]
    print l.queueSystem.submitTemplate.get()
    