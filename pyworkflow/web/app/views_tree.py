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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.web.pages import settings as django_settings
import pyworkflow.em as em

#===============================================================================
# PROTOCOL TREE
#===============================================================================

def loadProtTree(project):
    protCfg = project.getCurrentProtocolView()
    root = TreeItem('root', 'root', '', '')
    populateProtTree(root, protCfg)
#    return root
    
    html = convertProtTree(root)
    return html

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
    emProtocolsDict = em.getProtocols()
    
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
                        protItem = TreeItem(k, 'protocol_class', '/python_file.gif', None, protClassName, v)
                        item.childs.append(protItem)
        else:
            item.protClass = protClassName
            populateProtTree(item, sub)


def convertProtTree(tree):
    sections = tree.childs
    html = ''
    for section in sections:
        if section.tag == 'section':
            html += '<li><span class="section">'+ section.name +'</span><ul>'        
            html += convertProtTree(section)
            html += '</ul></li>'
        else:
            html += getProtChildrens(section)
    return html
        
        
def getProtChildrens(tree):
    html = ''
        
    if tree.tag == 'protocol':
        html += '<li><span class="protocol">'
        if tree.icon != None:
            html += '<img src="' + django_settings.ABSOLUTE_URL + '/resources/'+ tree.icon +'"/>'
        function = 'javascript:popup("/form/?protocolClass='+ tree.protClass +'")'
        html += ' <a href=' + function + '>'+ tree.name +'</a>'
        html += '</span></li>'
    
    elif tree.tag == 'protocol_class':
        html += '<li><span class="protocol_class">'
        if tree.icon != None:
            html += '<img src="' + django_settings.ABSOLUTE_URL + '/resources/'+ tree.icon +'"/>'
        function = 'javascript:popup("/form/?protocolClass='+ tree.protRealClass +'")'
        html += ' <a href=' + function + '>'+ tree.name +'</a>'
        html += '</span></li>'
    
    elif tree.tag == 'protocol_base':
        openItem = ""
        if tree.openItem != True:
            openItem = "closed"
        html += '<li class="'+ openItem +'"><span class="protocol_base">'
        if tree.icon != None:
            html += '<img src="' + django_settings.ABSOLUTE_URL + '/resources/' + tree.icon +'"/>'
        html += tree.name + '</span>'
    
        html += '<ul>'
        childrens = tree.childs
        for child in childrens:
            html += getProtChildrens(child)
        html += '</ul>'
    
    elif tree.tag == 'url':
        html += '<li><span class="protocol">'
        if tree.icon != None:
            html += '<img src="' + django_settings.ABSOLUTE_URL + '/resources/'+ tree.icon +'"/>'
        function = 'javascript:popup("'+ tree.protClass +'")'
        html += ' <a href=' + function + '>'+ tree.name +'</a>'
        html += '</span></li>'
    
    
    return html

#===============================================================================
# OBJECT TREE
#===============================================================================

def getGraphClassesNode(project):
    from pyworkflow.utils.graph import Graph
    classesGraph = Graph()
    
    # Method to create class nodes
    def createClassNode(classObj):
        """ Add the object class to hierarchy and 
        any needed subclass. """
        className = classObj.__name__
        classNode = classesGraph.getNode(className)
        
        if not classNode:
            classNode = classesGraph.createNode(className)
            if className != 'EMObject' and classObj.__bases__:
                baseClass = classObj.__bases__[0]
                for b in classObj.__bases__:
                    if b.__name__ == 'EMObject':
                        baseClass = b
                        break
                parent = createClassNode(baseClass)
                parent.addChild(classNode)
            classNode.count = 0
        return classNode
    
    g = project.getSourceGraph()

    for node in g.getNodes():
        id = node.getName()
        if id != 'PROJECT':
            obj = project.mapper.selectById(id)
            classNode = createClassNode(obj.getClass())
            classNode.count += 1
    
    return classesGraph
    
def populateObjTree(tree, elements):
    for node in elements:
        if node.getName() != "ROOT":
#            print "HERE->", node.getName()
#            print "CHILDS: ", [i.getName() for i in node.getChilds()]
            item = ObjItem(node.getName(), node.count)
            tree.childs.append(item)
            if len(node.getChilds()) > 0:
                populateObjTree(item, node.getChilds())
         
class ObjItem():
    def __init__(self, name, count=None):
        self.name = name
        self.count = count
        self.openItem = True
        self.childs = []
        

def convertObjTree(tree):
    sections = tree.childs
    for section in sections:
        html = getChildrens(section)
    return html
        
        
def getChildrens(tree):
    hasChilds = False
    
    html = '<span>'
    if len(tree.childs) > 0:
        html += tree.name
        hasChilds = True
    else:
        html = '<strong>'+tree.name +' ('+str(tree.count)+')</strong>'
    html += '</span>'
        
    if hasChilds:
        html += '<ul>'
        childrens = tree.childs
        for child in childrens:
            html += '<li>'
            html += getChildrens(child)
            html += '</li>'
        html += '</ul>'
    
    return html

