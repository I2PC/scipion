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

import sys
import os
from os.path import join, exists
import time
from ConfigParser import ConfigParser
import json

import pyworkflow as pw
from pyworkflow.object import Boolean, Integer, String, List, OrderedObject, CsvList
from pyworkflow.hosts import HostConfig, QueueConfig, QueueSystemConfig # we need this to be retrieved by mapper
from pyworkflow.mapper import SqliteMapper

PATH = os.path.dirname(__file__)


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
    def __init__(self, confs={}, **kwargs):
        OrderedObject.__init__(self, **kwargs)
        self.config = ProjectConfig()
        self.hostList = SettingList() # List to store different hosts configurations
        self.menuList = SettingList() # Store different menus
        self.protMenuList = SettingList() # Store different protocol configurations
        self.mapper = None # This should be set when load, or write
        self.graphView = Boolean(False)
        self.runSelection = CsvList(int) # Store selected runs
        
    def loadConfig(self, confs={}):
        """ Load values from configuration files.
        confs can contains the files for configuration .conf files. 
        """
        # Load configuration
        self.addMenus(confs.get('menus', None))
        print "protocols config: ", confs.get('protocols', None)
        self.addProtocols(confs.get('protocols', None))
        self.addHosts(confs.get('hosts', None))

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
    
    def getGraphView(self):
        return self.graphView.get()
    
    def setGraphView(self, value):
        self.graphView.set(value)

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
        else:
            if self.mapper is None:
                raise Exception("Can't write ProjectSettings without mapper or dbPath")

        self.mapper.deleteAll()
        self.mapper.insert(self)
        self.mapper.commit()
        
    def addProtocolMenu(self, protMenuConfig):
        self.protMenuList.append(protMenuConfig)
        
    def addProtocols(self, protocolsConf=None):
        """ Read the protocol configuration from a .conf
        file similar of the one in ~/.config/scipion/menu.conf,
        which is the default one when no file is passed.
        """
    
        # Helper function to recursively add items to a menu.
        def add(menu, item):
            "Add item (a dictionary that can contain more dictionaries) to menu"
            children = item.pop('children', [])
            subMenu = menu.addSubMenu(**item)  # we expect item={'text': ...}
            for child in children:
                add(subMenu, child)  # add recursively to sub-menu
    
        # Read menus from users' config file.
        cp = ConfigParser()
        cp.optionxform = str  # keep case (stackoverflow.com/questions/1611799)
        SCIPION_MENU = protocolsConf or os.environ['SCIPION_MENU']
        # Also mentioned in /scipion . Maybe we could do better.
    
        try:
            assert cp.read(SCIPION_MENU) != [], 'Missing file %s' % SCIPION_MENU
    
            # Populate the protocol menu from the config file.
            for menuName in cp.options('PROTOCOLS'):
                menu = ProtocolConfig(menuName)
                children = json.loads(cp.get('PROTOCOLS', menuName))
                for child in children:
                    add(menu, child)
                self.addProtocolMenu(menu)
        except Exception as e:
            sys.exit('Failed to read settings. The reported error was:\n  %s\n'
                     'To solve it, delete %s and run again.' % (e, SCIPION_MENU))

    def addHosts(self, hostConf=None):
        #TODO: parse several hosts and include in the settings
        host = HostConfig()
        host.label.set('localhost')
        host.hostName.set('localhost')
        host.hostPath.set(pw.SCIPION_USER_DATA)
        #TODO: parse the host conf specific for each host
        host.addQueueSystem(hostConf)
        self.addHost(host)
        
    def addMenus(self, menusConf=None):
        """ Add the menu to project windows. """
        #TODO: read this from a .conf file
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
        self.addMenu(menu)
    
        # Write another test menu
        menu = MenuConfig()
        m1 = menu.addSubMenu('Test')
        m1.addSubMenu('KK', icon='tree.gif')
        m1.addSubMenu('PP', icon='folderopen.gif')
    
        #writeConfig(menu, 'menu_test.xml')
        self.addMenu(menu)



class ProjectConfig(OrderedObject):
    """A simple base class to store ordered parameters"""
    def __init__(self, **args):
        OrderedObject.__init__(self, **args)
        self.icon = String('scipion_bn.xbm')
        self.logo = String('scipion_logo_small.png')


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

    def addSubMenu(self, text, value=None, **args):
        if 'icon' not in args:
            tag = args.get('tag', None)
            if tag == 'protocol':
                args['icon'] = 'python_file.gif'
            elif tag == 'protocol_base':
                args['icon'] = 'class_obj.gif'
        return MenuConfig.addSubMenu(self, text, value, **args)

