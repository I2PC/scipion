#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import sys
import os
import json
import datetime as dt
from collections import OrderedDict
from ConfigParser import ConfigParser  # FIXME Does not work in Python3

import pyworkflow as pw
import pyworkflow.object as pwobj
from pyworkflow.mapper import SqliteMapper


class ProjectSettings(pwobj.OrderedObject):
    """ Store settings related to a project. """

    COLOR_MODE_STATUS = 0
    COLOR_MODE_LABELS = 1
    COLOR_MODE_AGE = 2
    COLOR_MODES = (COLOR_MODE_STATUS, COLOR_MODE_LABELS, COLOR_MODE_AGE)

    def __init__(self, confs={}, **kwargs):
        pwobj.OrderedObject.__init__(self, **kwargs)
        self.config = ProjectConfig()
        # Store the current view selected by the user
        self.currentProtocolsView = pwobj.String()
        # Store the color mode: 0= Status, 1=Labels, ...
        self.colorMode = pwobj.Integer(ProjectSettings.COLOR_MODE_STATUS)
        self.nodeList = NodeConfigList()  # Store graph nodes positions and other info
        self.labelsList = LabelsList()  # Label list
        self.mapper = None  # This should be set when load, or write
        self.runsView = pwobj.Integer(1) # by default the graph view
        self.readOnly = pwobj.Boolean(False)
        self.runSelection = pwobj.CsvList(int) # Store selected runs
        self.dataSelection = pwobj.CsvList(int)  # Store selected runs
        # Some extra settings stored, now mainly used
        # from the webtools
        # Time when the project was created
        self.creationTime = pwobj.String(dt.datetime.now())
        # Number of days that this project is active
        # if None, the project will not expire
        # This is used in webtools where a limited time
        # is allowed for each project
        self.lifeTime = pwobj.Integer()
        # Set a disk quota for the project (in Gb)
        # if None, quota is unlimited
        self.diskQuota = pwobj.Integer()

    def commit(self):
        """ Commit changes made. """
        self.mapper.commit()

    def getRunsView(self):
        return self.runsView.get()

    def setRunsView(self, value):
        self.runsView.set(value)

    def getReadOnly(self):
        return self.readOnly.get()

    def setReadOnly(self, value):
        self.readOnly.set(value)

    def getCreationTime(self):
        return self.creationTime.datetime()

    def setCreationTime(self, value):
        self.creationTime.set(value)

    def getLifeTime(self):
        return self.lifeTime.get()

    def setLifeTime(self, value):
        self.lifeTime.set(value)

    def getConfig(self):
        return self.config

    def getProtocolView(self):
        return self.currentProtocolsView.get()

    def setProtocolView(self, protocolView):
        """ Set the new protocol Menu given its index.
        The new ProtocolMenu will be returned.
        """
        self.currentProtocolsView.set(protocolView)

    def getColorMode(self):
        return self.colorMode.get()

    def setColorMode(self, colorMode):
        """ Set the color mode to use when drawing the graph.
        """
        self.colorMode.set(colorMode)

    def statusColorMode(self):
        return self.getColorMode() == self.COLOR_MODE_STATUS

    def labelsColorMode(self):
        return self.getColorMode() == self.COLOR_MODE_LABELS

    def ageColorMode(self):
        return self.getColorMode() == self.COLOR_MODE_AGE

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

    def getNodes(self):
        return self.nodeList

    def getNodeById(self, nodeId):
        return self.nodeList.getNode(nodeId)

    def addNode(self, nodeId, **kwargs):
        return self.nodeList.addNode(nodeId, **kwargs)

    def getLabels(self):
        return self.labelsList

    @classmethod
    def load(cls, dbPath):
        """ Load a ProjectSettings from dbPath. """
        classDict = dict(globals())
        classDict.update(pwobj.__dict__)
        mapper = SqliteMapper(dbPath, classDict)
        settingList = mapper.selectByClass('ProjectSettings')
        n = len(settingList)

        if n == 0:
            raise Exception("Can't load ProjectSettings from %s" % dbPath)
        elif n > 1:
            raise Exception("Only one ProjectSettings is expected in db, found %d in %s" % (n, dbPath))

        settings = settingList[0]
        settings.mapper = mapper

        return settings


class ProjectConfig(pwobj.OrderedObject):
    """A simple base class to store ordered parameters"""
    def __init__(self, **args):
        pwobj.OrderedObject.__init__(self, **args)
        self.logo = pwobj.String('scipion_logo_small.png')
        # Do not store this object, unless we implement some kind of
        # icon customization
        self._objDoStore = False

class MenuConfig(object):
    """Menu configuration in a tree fashion.
    Each menu can contains submenus.
    Leaf elements can contain actions"""
    def __init__(self, text=None, value=None,
                 icon=None, tag=None, **kwargs):
        """Constructor for the Menu config item.
        Arguments:
          text: text to be displayed
          value: internal value associated with the item.
          icon: display an icon with the item
          tag: put some tags to items
        **args: pass other options to base class.
        """
        self.text = pwobj.String(text)
        self.value = pwobj.String(value)
        self.icon = pwobj.String(icon)
        self.tag = pwobj.String(tag)
        self.shortCut = pwobj.String(kwargs.get('shortCut', None))
        self.childs = pwobj.List()
        self.openItem = pwobj.Boolean(kwargs.get('openItem', False))

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

    def addSubMenu(self, text, value=None, shortCut=None, **args):
        if 'icon' not in args:
            tag = args.get('tag', None)
            if tag == 'protocol':
                args['icon'] = 'python_file.gif'
            elif tag == 'protocol_base':
                args['icon'] = 'class_obj.gif'

        args['shortCut'] = shortCut
        return MenuConfig.addSubMenu(self, text, value, **args)


class NodeConfig(pwobj.Scalar):
    """ Store Graph node information such as x, y. """
    
    def __init__(self, nodeId=0, x=None, y=None, selected=False, expanded=True):
        pwobj.Scalar.__init__(self)
        # Special node id 0 for project node
        self._values = {'id': nodeId, 
                        'x': pwobj.Integer(x).get(0), 
                        'y': pwobj.Integer(y).get(0), 
                        'selected': selected, 
                        'expanded': expanded,
                        'labels': []}
        
    def _convertValue(self, value):
        """Value should be a str with comma separated values
        or a list.
        """
        self._values = json.loads(value)
            
    def getObjValue(self):
        self._objValue = json.dumps(self._values)
        return self._objValue
    
    def get(self):
        return self.getObjValue()
        
    def getId(self):
        return self._values['id']
    
    def setX(self, x):
        self._values['x'] = x
        
    def getX(self):
        return self._values['x']
    
    def setY(self, y):
        self._values['y'] = y
        
    def getY(self):
        return self._values['y']
    
    def setPosition(self, x, y):
        self.setX(x)
        self.setY(y)
        
    def getPosition(self):
        return self.getX(), self.getY()
        
    def setSelected(self, selected):
        self._values['selected'] = selected
        
    def isSelected(self):
        return self._values['selected']
    
    def setExpanded(self, expanded):
        self._values['expanded'] = expanded
        
    def isExpanded(self):
        return self._values['expanded']

    def setLabels(self, labels):
        self._values['labels'] = labels

    def getLabels(self):
        return self._values.get('labels', None)
    
    def __str__(self):
        return 'NodeConfig: %s' % self._values
    
    
class NodeConfigList(pwobj.List):
    """ Store all nodes information items and 
    also store a dictionary for quick access
    to nodes query.
    """
    def __init__(self):
        self._nodesDict = {}
        pwobj.List.__init__(self)
        
    def getNode(self, nodeId):
        return self._nodesDict.get(nodeId, None)
    
    def addNode(self, nodeId, **kwargs):
        node = NodeConfig(nodeId, **kwargs)
        self._nodesDict[node.getId()] = node
        self.append(node)
        return node
        
    def updateDict(self):
        self._nodesDict.clear()
        for node in self:
            self._nodesDict[node.getId()] = node
            
    def clear(self):
        pwobj.List.clear(self)
        self._nodesDict.clear()
        

class Label(pwobj.Scalar):
    """ Store Label information """

    def __init__(self, labelId=None, name='', color=None):
        pwobj.Scalar.__init__(self)
        # Special node id 0 for project node
        self._values = {'id': labelId,
                        'name': name,
                        'color': color}

    def _convertValue(self, value):
        """Value should be a str with comma separated values
        or a list.
        """
        self._values = json.loads(value)

    def getObjValue(self):
        self._objValue = json.dumps(self._values)
        return self._objValue

    def get(self):
        return self.getObjValue()

    def getId(self):
        return self._values['id']

    def getName(self):
        return self._values['name']

    def setName(self, newName):
        self._values['name'] = newName

    def setColor(self, color):
        self._values['color'] = color

    def getColor(self):
        return self._values.get('color', None)

    def __str__(self):
        return 'Label: %s' % self._values

    def __eq__(self, other):
        return self.getName() == other.getName()


class LabelsList(pwobj.List):
    """ Store all labels information"""
    def __init__(self):
        self._labelsDict = {}
        pwobj.List.__init__(self)

    def getLabel(self, name):
        return self._labelsDict.get(name, None)

    def addLabel(self, label):
        self._labelsDict[label.getName()] = label
        self.append(label)
        return label

    def updateDict(self):
        self._labelsDict.clear()
        for label in self:
            self._labelsDict[label.getName()] = label

    def deleteLabel(self, label):
        self._labelsDict.pop(label.getName())
        self.remove(label)

    def clear(self):
        pwobj.List.clear(self)
        self._labelDict.clear()


class ProtocolTreeConfig:
    """ Handler class that groups functions and constants
    related to the protocols tree configuration.
    """
    ALL_PROTOCOLS = "All"
    TAG_PROTOCOL_DISABLED = 'protocol-disabled'
    TAG_PROTOCOL = 'protocol'
    TAG_SECTION = 'section'
    TAG_PROTOCOL_GROUP = 'protocol_group'
    PLUGIN_CONFIG_PROTOCOLS = 'protocols.conf'

    @classmethod
    def getProtocolTag(cls, isInstalled):
        """ Return the proper tag depending if the protocol is installed or not.
        """
        return cls.TAG_PROTOCOL if isInstalled else cls.TAG_PROTOCOL_DISABLED

    @classmethod
    def isAFinalProtocol(cls, v, k):
        if (issubclass(v, pw.viewer.ProtocolViewer) or
            v.isBase() or v.isDisabled()):
            return False

        # To remove duplicated protocol, ProtMovieAlignment turns into OF:
        # ProtMovieAlignment = XmippProtOFAlignment
        return v.__name__ == k

    @classmethod
    def __addToTree(cls, menu, item, checkFunction=None):
        """ Helper function to recursively add items to a menu.
        Add item (a dictionary that can contain more dictionaries) to menu
        If check function is added will use it to check if the value must be added.
        """
        children = item.pop('children', [])

        if checkFunction is not None:
            add = checkFunction(item)
            if not add:
                return
        subMenu = menu.addSubMenu(**item)  # we expect item={'text': ...}
        for child in children:
            cls.__addToTree(subMenu, child, checkFunction)  # add recursively to sub-menu

        return subMenu

    @classmethod
    def __inSubMenu(cls, child, subMenu):
        """
        Return True if child belongs to subMenu
        """
        for ch in subMenu:
            if child['tag'] == cls.TAG_PROTOCOL:
                if ch.value == child['value']:
                    return ch
            elif ch.text == child['text']:
                return ch
        return None

    @classmethod
    def _orderSubMenu(cls, session):
        """
        Order all children of a given session:
        The protocols first, then the sessions(the 'more' session at the end)
        """
        lengthSession = len(session.childs)
        if lengthSession > 1:
            childs = session.childs
            lastChildPos = lengthSession - 1
            if childs[lastChildPos].tag == cls.TAG_PROTOCOL:
                for i in range(lastChildPos - 1, -1, -1):
                    if childs[i].tag == cls.TAG_PROTOCOL:
                        break
                    else:
                        tmp = childs[i+1]
                        childs[i+1] = childs[i]
                        childs[i] = tmp
            else:
                for i in range(lastChildPos - 1, -1, -1):
                    if childs[i].tag == cls.TAG_PROTOCOL:
                        break
                    elif 'more' in str(childs[i].text).lower():
                        tmp = childs[i+1]
                        childs[i+1] = childs[i]
                        childs[i] = tmp

    @classmethod
    def __findTreeLocation(cls, subMenu, children, parent):
        """
        Locate the protocol position in the given view
        """
        for child in children:
            sm = cls.__inSubMenu(child, subMenu)
            if sm is None:
                cls.__addToTree(parent, child, cls.__checkItem)
                cls._orderSubMenu(parent)
            elif child['tag'] == cls.TAG_PROTOCOL_GROUP or child['tag'] == cls.TAG_SECTION:
                cls.__findTreeLocation(sm.childs, child['children'], sm)

    @classmethod
    def __checkItem(cls, item):
        """ Function to check if the protocol has to be added or not.
        Params:
            item: {"tag": "protocol", "value": "ProtImportMovies",
                   "text": "import movies"}
        """
        if item["tag"] != cls.TAG_PROTOCOL:
            return True

        # It is a protocol as this point, get the class name and
        # check if it is disabled
        protClassName = item["value"]
        protClass = pw.em.Domain.getProtocols().get(protClassName)

        return False if protClass is None else not protClass.isDisabled()

    @classmethod
    def __addAllProtocols(cls, protocols):
        # Add all protocols
        # FIXME: Check why this import is here
        allProts = pw.em.Domain.getProtocols()

        # Sort the dictionary
        allProtsSorted = OrderedDict(sorted(allProts.items(),
                                            key=lambda e: e[1].getClassLabel()))
        allProtMenu = ProtocolConfig(cls.ALL_PROTOCOLS)
        packages = {}

        # Group protocols by package name
        for k, v in allProtsSorted.iteritems():
            if cls.isAFinalProtocol(v, k):
                packageName = v.getClassPackageName()
                # Get the package submenu
                packageMenu = packages.get(packageName)

                # If no package menu available
                if packageMenu is None:
                    # Add it to the menu ...
                    packageLine = {"tag": "package", "value": packageName,
                                   "text": packageName}
                    packageMenu = cls.__addToTree(allProtMenu, packageLine)

                    # Store it in the dict
                    packages[packageName] = packageMenu

                # Add the protocol
                tag = cls.getProtocolTag(v.isInstalled())

                protLine = {"tag": tag, "value": k,
                            "text": v.getClassLabel(prependPackageName=False)}

                # If it's a new protocol
                if v.isNew() and v.isInstalled():
                    # add the new icon
                    protLine["icon"] = "newProt.png"

                cls.__addToTree(packageMenu, protLine)

        protocols[cls.ALL_PROTOCOLS] = allProtMenu

    @classmethod
    def __addProtocolsFromConf(cls, protocols, protocolsConfPath):
        """
        Load the protocols in the tree from a given protocols.conf file,
        either the global one in Scipion or defined in a plugin.
        """
        # Populate the protocols menu from the plugin config file.
        if os.path.exists(protocolsConfPath):
            cp = ConfigParser()
            cp.optionxform = str  # keep case
            cp.read(protocolsConfPath)
            #  Ensure that the protocols section exists
            if cp.has_section('PROTOCOLS'):
                for menuName in cp.options('PROTOCOLS'):
                    if menuName not in protocols:  # The view has not been inserted
                        menu = ProtocolConfig(menuName)
                        children = json.loads(cp.get('PROTOCOLS', menuName))
                        for child in children:
                            cls.__addToTree(menu, child, cls.__checkItem)
                        protocols[menuName] = menu
                    else:  # The view has been inserted
                        menu = protocols.get(menuName)
                        children = json.loads(cp.get('PROTOCOLS',
                                                     menuName))
                        cls.__findTreeLocation(menu.childs, children, menu)

    @classmethod
    def load(cls, protocolsConf):
        """ Read the protocol configuration from a .conf file similar to the
        one in scipion/config/protocols.conf,
        which is the default one when no file is passed.
        """

        protocols = OrderedDict()
        # Read the protocols.conf from Scipion (base) and create an initial
        # tree view
        cls.__addProtocolsFromConf(protocols, protocolsConf)

        # Read the protocols.conf of any installed plugin
        pluginDict = pw.em.Domain.getPlugins()
        pluginList = pluginDict.keys()
        for pluginName in pluginList:
            try:
                # Locate the plugin protocols.conf file
                protocolsConfPath = os.path.join(pluginDict[pluginName].__path__[0],
                                                 cls.PLUGIN_CONFIG_PROTOCOLS)
                cls.__addProtocolsFromConf(protocols, protocolsConfPath)
            except Exception as e:
                print('Failed to read settings. The reported error was:\n  %s\n'
                      'To solve it, fix %s and run again.' % (
                            e, protocolsConfPath))

            # Add all protocols to All view
        cls.__addAllProtocols(protocols)

        return protocols