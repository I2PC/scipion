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

class MenuConfig(List):
    """Menu configuration in a tree fashion.
    Each menu can contains submenus.
    Leaf elements can contain actions"""
    def __init__(self, text=None, icon=None, key=None, 
                 action=None, **args):
        List.__init__(self, **args)
        self.text = String(text)
        self.icon = String(icon)
        self.key = String(key)
        self.action = String(action)
        
    def addSubMenu(self, text, icon=None, **args):
        subMenu = type(self)(text, icon, **args)
        self.append(subMenu)
        return subMenu
    
    def _getStr(self, prefix):
        s = prefix + "%s text = %s, icon = %s\n" % (self.getClassName(), self.text.get(), self.icon.get())
        for sub in self:
            s += sub._getStr(prefix + "  ")
        return s
            
    def __str__(self):
        return self._getStr(' ')
    

class ProtocolConfig(MenuConfig):
    """Store protocols configuration """
    pass
    
class ConfigXmlMapper(XmlMapper):
    """Sub-class of XmlMapper to store configurations"""
    def __init__(self, filename, dictClasses=None, **args):
        if not exists(filename):
            raise Exception('Config file: %s does not exists' % filename)
        XmlMapper.__init__(self, filename, dictClasses, **args)
        self.setClassTag('MenuConfig.MenuConfig', 'class_only')
        self.setClassTag('MenuConfig.String', 'attribute')
        self.setClassTag('ProtocolConfig.ProtocolConfig', 'class_only')
        self.setClassTag('ProtocolConfig.String', 'attribute')
        
    def getConfig(self):
        return self.getAll()[0]

def writeConfig(config, fn):
    fn = getConfigPath(fn)
    if exists(fn):
        os.remove(fn)
    mapper = ConfigXmlMapper(fn, globals())
    mapper.insert(config)
    mapper.commit()
    
def writeDefaults():
    """Write default configuration files"""
    # Write menu configuration
    menu = MenuConfig()
    projMenu = menu.addSubMenu('Project')
    projMenu.addSubMenu('Browse files', 'folderopen.gif', action='browse')
    projMenu.addSubMenu('Remove temporary files', 'delete.gif')
    projMenu.addSubMenu('Clean project')
    projMenu.addSubMenu('Exit')
    
    helpMenu = menu.addSubMenu('Help')
    helpMenu.addSubMenu('Online help', 'online_help.gif')
    helpMenu.addSubMenu('About')
    
    writeConfig(menu, 'menu_default.xml')
    
    # Write another test menu
    menu = MenuConfig()
    m1 = menu.addSubMenu('Test')
    m1.addSubMenu('KK', 'tree.gif')
    m1.addSubMenu('PP', 'folderopen.gif')
    writeConfig(menu, 'menu_test.xml')
    
    # Write protocols configuration
    menu = ProtocolConfig()
    m1 = menu.addSubMenu('Preprocessing')
    
    m2 = m1.addSubMenu('Micrographs')
    m2.addSubMenu('Import', action='ProtImportMicrographs')
    m2.addSubMenu('Screen', action='ProtScreenMicrographs')
    m2.addSubMenu('Downsample', action='ProtDownsampleMicrographs')
    
    m2 = m1.addSubMenu('Particle Picking')
    m2.addSubMenu('Manual')
    m2.addSubMenu('Supervised')
    m2.addSubMenu('Automatic')
    
    m1 = menu.addSubMenu('2D')
    
    m1.addSubMenu('Align', action='ProtAlign')
    m1.addSubMenu('Classify', action='ProtClassify')
    m1.addSubMenu('Align+Classify', action='ProtAlignClassify')

    
    writeConfig(menu, 'protocols_default.xml')
    
    # Write global configuration
    config = ProjectConfig(menu='menu_default.xml',
                           protocols='protocols_default.xml')
    writeConfig(config, 'configuration.xml')
    
if __name__ == '__main__':
    writeDefaults()
    #Write default configurations
    