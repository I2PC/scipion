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
from pyworkflow.utils.path import getHomePath, makeFilePath, cleanPath
import pyworkflow as pw
from pyworkflow.object import *
from pyworkflow.hosts import *
from pyworkflow.mapper import SqliteMapper, XmlMapper

PATH = os.path.dirname(__file__)
SETTINGS = join(pw.HOME,'..','settings')

def getConfigPath(filename):
    """Return a configuration filename from settings folder"""
    return join(SETTINGS, filename)

        
def loadSettings(dbPath):
    """ Load a ProjectSettings from dbPath. """
    mapper = SqliteMapper(dbPath, globals())
    settingList = mapper.selectByClass('ProjectSettings')
    n = len(settingList)
    if n == 0:
        raise Exception("Can't load ProjectSettings from %s" % dbPath)
    elif n > 1:
        raise Exception("Only one ProjectSettings is expected in db, found %d in %s" % (n, dbPath))
    
    settings = settingList[0]
    settings.mapper = mapper
    
    return settings


class SettingList(List):
    """ Basically a list that also store an index of the last selection. """
    def __init__(self, **args):
        List.__init__(self, **args)
        self.currentIndex = Integer(0)
        
    def getIndex(self):
        return self.currentIndex.get()
    
    def setIndex(self, i):
        self.currentIndex.set(i)
        
    def getItem(self):
        """ Get the item corresponding to current index. """
        return self[self.getIndex()]
    
    
class ProjectSettings(OrderedObject):
    """ This class will store settings related to a project. """
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self.config = ProjectConfig()
        self.hostList = SettingList() # List to store different hosts configurations
        self.menuList = SettingList() # Store different menus
        self.protMenuList = SettingList() # Store different protocol configurations
        self.mapper = None # This should be set when load, or write 
        self.graphView = Boolean(False)
        
    def commit(self):
        """ Commit changes made. """
        self.mapper.commit()
        
    def addHost(self, hostConfig):
        self.hostList.append(hostConfig)
        
    def getHosts(self):
        return self.hostList
    
    def getHostById(self, hostId):
        return self.mapper.selectById(hostId)
    
    def getHostByLabel(self, hostLabel):
        for host in self.hostList:
            if host.label == hostLabel:
                return host
        return None
    
    def saveHost(self, host, commit=False):
        """ Save a host for project settings.
            If the hosts exists it is updated, else it is created.
        params:
            host: The host to update or create.
        """
        if not host in self.hostList:
            self.addHost(host)
        self.mapper.store(host)
        if commit:
            self.commit()
    
    def deleteHost(self, host, commit=False):
        """ Delete a host of project settings.
        params:
            hostId: The host id to delete.
        """
        if not host in self.hostList:
            raise Exception('Deleting host not from host list.')
        self.hostList.remove(host)
        self.mapper.delete(host)
        if commit:
            self.commit()
        
    def addMenu(self, menuConfig):
        self.menuList.append(menuConfig)
        
    def addProtocolMenu(self, protMenuConfig):
        self.protMenuList.append(protMenuConfig)
        
    def getConfig(self):
        return self.config
    
    def getCurrentMenu(self):
        """ Now by default return element at index 0,
        later should be stored the current index. 
        """
        return self.menuList.getItem()
    
    def getCurrentProtocolMenu(self):
        return self.protMenuList.getItem()
    
    def setCurrentProtocolMenu(self, index):
        """ Set the new protocol Menu given its index.
        The new ProtocolMenu will be returned.
        """
        self.protMenuList.setIndex(index)
        return self.getCurrentProtocolMenu()
    
    def write(self, dbPath=None):
        self.setName('ProjectSettings')
        if dbPath is not None:
            self.mapper = SqliteMapper(dbPath, globals())
            self.mapper.deleteAll()
            self.mapper.insert(self)
        else:
            if self.mapper is None:
                raise Exception("Can't write ProjectSettings without mapper or dbPath")
            self.mapper.store(self)
        
        self.mapper.commit()


class ProjectConfig(OrderedObject):
    """A simple base class to store ordered parameters"""
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
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
        self.childs = List()
        self.openItem = Boolean(args.get('openItem', False))
        
    def addSubMenu(self, text, value=None, **args):
        subMenu = type(self)(text, value, **args)
        self.childs.append(subMenu)
        return subMenu
    
    def __iter__(self):
        for v in self.childs:
            yield v
                
    def __len__(self):
        return len(self.childs)
    
    def isEmpty(self):
        return len(self.childs) == 0
            
    
class ProtocolConfig(MenuConfig):
    """Store protocols configuration """
    def __init__(self, text=None, value=None, **args):
        MenuConfig.__init__(self, text, value, **args)
        if 'openItem' not in args:
            self.openItem.set(self.tag.get() != 'protocol_base')
            
    
def addMenus(settings):
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
    
    #writeConfig(menu, 'menu_default.xml')
    settings.addMenu(menu)
    
    # Write another test menu
    menu = MenuConfig()
    m1 = menu.addSubMenu('Test')
    m1.addSubMenu('KK', icon='tree.gif')
    m1.addSubMenu('PP', icon='folderopen.gif')
    
    #writeConfig(menu, 'menu_test.xml')
    settings.addMenu(menu)
    
def addProtocols(settings):
    """ Write protocols configuration. """
    menu = ProtocolConfig("Protocols SPA")
    
    # ------------------- Micrographs ----------------------------
    m1 = menu.addSubMenu('Micrographs', tag='section')
    
    m1.addSubMenu(' Import', value='ProtImportMicrographs', 
                  tag='protocol', icon='bookmark.png')
    m1.addSubMenu('Preprocess', value='ProtPreprocessMicrographs',
                  tag='protocol_base')
    m1.addSubMenu('CTF estimation', value='ProtCTFMicrographs',
                  tag='protocol_base')
    
    # ------------------- Particles ----------------------------
    m1 = menu.addSubMenu('Particles', tag='section')
    
    m1.addSubMenu('Import', value='ProtImportParticles', 
                  tag='protocol', icon='bookmark.png')
    m1.addSubMenu('Picking', value='ProtParticlePicking',
                  tag='protocol_base')
    m1.addSubMenu('Extract', value='ProtExtractParticles',
                  tag='protocol_base')    
    m1.addSubMenu('Process', value='ProtProcessParticles',
                  tag='protocol_base')   
    
    # ------------------- 2D ----------------------------
    m1 = menu.addSubMenu('2D', tag='section')
    
    m1.addSubMenu('Align', value='ProtAlign',
                  tag = 'protocol_base')
    m1.addSubMenu('Classify', value='ProtClassify',
                  tag = 'protocol_base')
    m1.addSubMenu('Align+Classify', value='ProtAlignClassify',
                  tag = 'protocol_base')
    
    # ------------------- 3D ----------------------------
    m1 = menu.addSubMenu('3D', tag='section')
    
    m1.addSubMenu(' Import', value='ProtImportVolumes', 
                 tag='protocol', icon='bookmark.png')   
    m1.addSubMenu('Preprocess', value='ProtPreprocessVolumes',
                  tag='protocol_base')
    m1.addSubMenu('Refine', value='ProtRefine3D',
                  tag='protocol_base')
    m1.addSubMenu('Classify', value='ProtClassify3D',
                  tag='protocol_base')
    m1.addSubMenu('Validate', value='ProtValidate3D',
                  tag='protocol_base')
    m1.addSubMenu('Mask', value='ProtCreateMask3D',
                  tag='protocol_base')   
    
    settings.addProtocolMenu(menu)
    
    addSpiderMDAProtocols(settings)
    
    
def addSpiderMDAProtocols(settings):
    """ Write protocols related to Spider MDA workflow. """
    menu = ProtocolConfig("MDA workflow")
    
    # ------------------- Particles ----------------------------
    m1 = menu.addSubMenu('MDA workflow', tag='section')
    
    m1.addSubMenu(' Import particles', tag='protocol', icon='bookmark.png',
                  value='ProtImportParticles')
    m1.addSubMenu(' Filter (optional)', tag='protocol',
                  value='SpiderProtFilter')
    m1.addSubMenu(' Align', tag='protocol_base', openItem=True, 
                  value='ProtAlign')
    m1.addSubMenu(' Create mask (optional)', tag='protocol',
                  value='SpiderProtCustomMask')
    m1.addSubMenu(' Dimension reduction', tag='protocol',  
                  value='SpiderProtCAPCA')
    m1.addSubMenu(' Classification', tag='protocol',
                  value='SpiderProtClassifyWard')
    
    m2 = menu.addSubMenu('Protocol MDA', tag='protocol',
                         value='SpiderWfMDA')
            
    settings.addProtocolMenu(menu)
       
    
def getScipionHome(userHome):
    """ Returns default SCIPION_HOME from HOME. """
    return join(userHome, 'Scipion')
  
def setQueueSystem(host, maxCores):
    host.mpiCommand.set('mpirun -np %(JOB_NODES)d -bynode %(COMMAND)s')
    host.queueSystem = QueueSystemConfig()
    queueSys = host.queueSystem
    #queueSys = QueueSystemConfig()
    queueSys.name.set('PBS/TORQUE')
    queueSys.mandatory.set(False)
    queueSys.submitCommand.set('qsub %(JOB_SCRIPT)s')
    queueSys.submitTemplate.set("""
#!/bin/bash
### Inherit all current environment variables
#PBS -V
### Job name
#PBS -N %(JOB_NAME)s
### Queue name
###PBS -q %(JOB_QUEUE)s
### Standard output and standard error messages
#PBS -k eo
### Specify the number of nodes and thread (ppn) for your job.
#PBS -l nodes=%(JOB_NODES)d:ppn=%(JOB_THREADS)d
### Tell PBS the anticipated run-time for your job, where walltime=HH:MM:SS
#PBS -l walltime=%(JOB_HOURS)d:00:00
# Use as working dir the path where qsub was launched
WORKDIR=$PBS_O_WORKDIR
#################################
### Set environment varible to know running mode is non interactive
export XMIPP_IN_QUEUE=1
### Switch to the working directory;
cd $WORKDIR
# Make a copy of PBS_NODEFILE 
cp $PBS_NODEFILE %(JOB_NODEFILE)s
# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`
# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`
### Display the job context
echo Running on host `hostName`
echo Time is `date`
echo Working directory is `pwd`
echo Using ${NPROCS} processors across ${NNODES} nodes
echo PBS_NODEFILE:
cat $PBS_NODEFILE
#################################

%(JOB_COMMAND)s
""")
    queueSys.cancelCommand.set('canceljob %(JOB_ID)d')
    queueSys.checkCommand.set('qstat %(JOB_ID)d')
    
    queue = QueueConfig()
    queue.maxCores.set(maxCores)
    queue.allowMPI.set(True)
    queue.allowThreads.set(True)
    
    queueSys.queues = List()
    queueSys.queues.append(queue)
    
    
    
    
def addHosts(settings):
    host = HostConfig()
    host.label.set('localhost')
    host.hostName.set('localhost')
    host.hostPath.set(pw.SCIPION_HOME)   
    setQueueSystem(host, maxCores=4)
    
    #writeConfig(host, dbPath, mapperClass=HostMapper, clean=True)
    settings.addHost(host)
    
    host = HostConfig()
    host.label.set('glassfishdev')
    host.hostName.set('glassfishdev.cnb.csic.es')
    host.userName.set('apoza')
    host.password.set('BF6fYiFYiYD')
    host.hostPath.set(getScipionHome('/home/apoza'))
    
    #writeConfig(host, dbPath, mapperClass=HostMapper, clean=False)
    settings.addHost(host)
    
    host = HostConfig()
    host.label.set('crunchy')
    host.hostName.set('crunchy.cnb.csic.es')
    host.userName.set('apoza')
    host.password.set('nerfyeb4f1v')
    host.hostPath.set(getScipionHome('/gpfs/fs1/home/bioinfo/apoza'))
    setQueueSystem(host, maxCores=64)
    
    #writeConfig(host, dbPath, mapperClass=HostMapper, clean=False)
    settings.addHost(host)
    
def writeDefaults():
    settings = ProjectSettings()
    addMenus(settings)
    addProtocols(settings)
    addHosts(settings)
    
    #writeHosts(join(getHomePath(), SCIPION_PATH, SETTINGS_PATH) ) 
    
    
    #writeConfig(config, 'configuration.xml')
    dbPath = pw.SETTINGS
    print "Writing default settings to: ", dbPath
    if exists(dbPath):
        print "File exist, deleting..."
        cleanPath(dbPath)
    settings.write(dbPath)
    
    
    
if __name__ == '__main__':
    #Write default configurations
    writeDefaults()
    # Read and print to check
    settings = loadSettings(pw.SETTINGS)
    settings.printAll()
    
