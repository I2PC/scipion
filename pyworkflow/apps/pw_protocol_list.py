#!/usr/bin/env python
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
List all existing protocols within Scipion
"""

import sys
from pyworkflow.viewer import Viewer
import pyworkflow.em as em


def getFirstLine(doc):
    """ Get the first non empty line from doc. """
    if doc:
        for lines in doc.split('\n'):
            l = lines.strip()
            if l:
                return l
    return ''


if __name__ == '__main__':
    count = 0
    withDoc = '--with-doc' in sys.argv
    asciidoc = '--asciidoc' in sys.argv
    
    emProtocolsDict = em.getProtocols()
    
    emCategories = [('Imports', em.ProtImport, []),
                    ('Micrographs', em.ProtMicrographs, []),
                    ('Particles', em.ProtParticles, []),
                    ('2D', em.Prot2D, []),
                    ('3D', em.Prot3D, [])]
    
    protDict = {}
    
    # Group protocols by package name
    for k, v in emProtocolsDict.iteritems():
        packageName = v.getClassPackageName()
        
        if packageName not in protDict:
            protDict[packageName] = []
        
        if not issubclass(v, Viewer) and not v.isBase():
            protDict[packageName].append((k, v))
            
            for c in emCategories:
                if issubclass(v, c[1]):
                    c[2].append((k, v))
           
    
    def iterGroups(protDict):
        groups = protDict.keys()
        groups.sort(key=lambda x: 1000-len(protDict[x]))
        
        for g in groups:
            yield g, protDict[g]
            
    def printProtocols(prots):
        protList = [(k, v, v.getClassLabel()) for k, v in prots]
        protList.sort(key=lambda x: x[2])
        
        for k, v, l in protList:
            doc = getFirstLine(v.__doc__) if withDoc else ''
            print "* link:%s[%s]: %s" % (k, l, doc)
        
        
    if asciidoc:
        print ":toc:\n:toc-placement!:\n\ntoc::[]\n"
        
        print "\n== By Categories\n"
        for c in emCategories:
            print "\n=== %s\n" % c[0]
            printProtocols(c[2])
        
        print "\n== By Packages\n"
        for group, prots in iterGroups(protDict):
            print "\n=== ", group, "(%d protocols)\n" % len(prots)
            printProtocols(prots)
        
    else:
        for group, prots in iterGroups(protDict):
            print "-" * 100
            print "Package: ", group, "(%d protocols)" % len(prots)
            for k, v in prots:
                print "   %s ( %s )" % (k, v.getClassLabel())
                if withDoc:
                    print "      doc: ", v.__doc__
                #count += 1
    
    
            
