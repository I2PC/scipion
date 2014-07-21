# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
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

#===============================================================================
# PROTOCOL TREE
#===============================================================================

def loadProtTree(project):
    protCfg = project.getSettings().getCurrentProtocolMenu()
    root = TreeItem('root', 'root', '', '')
    populateProtTree(root, protCfg)
    return root

class TreeItem():
    def __init__(self, name, tag, icon, openItem, protClassName=None, protClass=None):
        if protClass is None:
            self.name = name
        else:
            self.name = protClass.getClassLabel()
        self.tag = tag
        self.icon = icon
        self.openItem = openItem
        self.protClass = protClassName
        self.protRealClass = name
        self.childs = []
        
        
def populateProtTree(tree, obj):    
    from pyworkflow.em import emProtocolsDict
    
    for sub in obj:
        text = sub.text.get()
        value = sub.value.get(text)
        tag = sub.tag.get('')
        icon = sub.icon.get('')
        openItem = sub.openItem.get()
        item = TreeItem(text, tag, icon, openItem)
        tree.childs.append(item)
        # If have tag 'protocol_base', fill dynamically with protocol sub-classes
        protClassName = value.split('.')[-1]  # Take last part
        if sub.value.hasValue() and tag == 'protocol_base':
            prot = emProtocolsDict.get(protClassName, None)
            if prot is not None:
                for k, v in emProtocolsDict.iteritems():
                    if not v is prot and issubclass(v, prot):
                        protItem = TreeItem(k, 'protocol_class', 'python_file.gif', None, protClassName, v)
                        item.childs.append(protItem)
        else:
            item.protClass = protClassName
            populateProtTree(item, sub)


#===============================================================================
# OBJECT TREE
#===============================================================================

def loadObjTree(project):
    root = TreeItem('root', 'root', '', '')
    g = project.getSourceGraph()
    
    listObjs = []
    for node in g.getNodes():
        id = node.getName()
        if id != 'PROJECT':
            obj = project.mapper.selectById(id)
            listObjs.append(obj)
#            className = obj.__class__.__name__
            
    populateObjTree(root, listObjs)
    
    return root

def populateObjTree(tree, listObjs):
    
    for obj in listObjs:
        item = TreeItem(text, tag, icon, openItem)

        tree.childs.append(item)
    
    pass


    