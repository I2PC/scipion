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

from pyworkflow.object import *
from pyworkflow.mapper import SqliteMapper, XmlMapper

PATH = os.path.dirname(__file__)
SETTINGS = join(PATH, 'settings')

def getConfigPath(filename):
    """Return a configuration filename from settings folder"""
    return join(SETTINGS, filename)

class Configuration(OrderedObject):
    """A simple base class to store ordered parameters"""
    pass

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
        subMenu = MenuConfig(text, icon, **args)
        self.append(subMenu)
        return subMenu
    
    def hasSubMenu(self):
        return len(self) == 0
    
    def _getStr(self, prefix):
        s = prefix + "MenuConfig text = %s, icon = %s\n" % (self.text.get(), self.icon.get())
        for sub in self:
            s += sub._getStr(prefix + "  ")
        return s
            
    def __str__(self):
        return self._getStr(' ')
    
class ConfigXmlMapper(XmlMapper):
    """Sub-class of XmlMapper to store configurations"""
    def __init__(self, filename, dictClasses=None, **args):
        XmlMapper.__init__(self, filename, dictClasses, **args)
        self.setClassTag('MenuConfig.MenuConfig', 'class_only')
        self.setClassTag('MenuConfig.String', 'attribute')

def writeDefaults():
    """Write default configuration files"""
    # Write menu configuration
    menu = MenuConfig()
    projMenu = menu.addSubMenu('Project')
    projMenu.addSubMenu('Browse files', 'folderopen.gif')
    projMenu.addSubMenu('Remove temporary files', 'delete.gif')
    projMenu.addSubMenu('Clean project')
    projMenu.addSubMenu('Exit')
    
    helpMenu = menu.addSubMenu('Help')
    helpMenu.addSubMenu('Online help', 'online_help.gif')
    helpMenu.addSubMenu('About')
    
    menuFn = join(SETTINGS, 'menu_default.xml')
    if exists(menuFn):
        os.remove(menuFn)
#    mapper = ConfigXmlMapper(menuFn, globals())
#    mapper.insert(menu)
    #mapper.commit()
    
    menu2 = MenuConfig()
    m1 = menu2.addSubMenu('Test')
    m2 = m1.addSubMenu('KK', 'tree.gif')
    m3 = m1.addSubMenu('PP', 'folderopen.gif')
    
    #print menu2
    print "creating mapper"
    mapper2 = ConfigXmlMapper(menuFn, globals())
    mapper2.insert(menu2)
    print "before commit"
    mapper2.commit()
    
    #print menu
    
    # Write protocols configuration
    
    # Write 
    
if __name__ == '__main__':
    writeDefaults()
    #Write default configurations
    