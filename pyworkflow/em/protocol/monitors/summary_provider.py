# -*- coding: utf-8 -*-
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow.object as pwobj
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.protocol import getUpdatedProtocol


class SummaryProvider(TreeProvider):
    """Create the tree elements for a Protocol run"""
    def __init__(self, protocol):
        TreeProvider.__init__(self)
        self.protocol = protocol
        self.getColumns = lambda: [('Name', 300), ('Output', 150),
                                   ('Number', 100)]
        self._parentDict = {}
        self.acquisition = []
        self.refreshObjects()

    def getObjects(self):
        return self._objects

    def refreshObjects(self):
        objects = []
        objIds = []  # need to store ids too to avoid duplication in runs table

        def addObj(objId, name, output='', size='', parent=None):
            if objId not in objIds:
                obj = pwobj.Object(objId=objId)
                obj.name = name
                obj.output = output
                obj.outSize = size
                obj._objParent = parent
                objIds.append(objId)
                objects.append(obj)
                return obj
            else:
                return None

        inputProts = self.protocol.getInputProtocols()
        runs = [getUpdatedProtocol(p) for p in inputProts]
        g = self.protocol.getProject().getGraphFromRuns(runs)

        nodes = g.getRoot().iterChildsBreadth()

        for n in nodes:
            prot = n.run
            pobj = addObj(prot.getObjId(),
                          '%s (id=%s)' % (prot.getRunName(), prot.strId()))

            for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
                outSet.load()
                outSet.loadAllProperties()
                addObj(outSet.getObjId(), '', outName, outSet.getSize(), pobj)
                outSet.close()
                # Store acquisition parameters in case of the import protocol
                from pyworkflow.em import ProtImportImages
                # NOTE by Yaiza: we force the string containing the Å to be unicode
                # because this is the encoding used when generating report in report_html.py
                if isinstance(prot, ProtImportImages):
                    self.acquisition = [("Microscope Voltage: ",
                                         prot.voltage.get()),
                                        ("Spherical aberration: ",
                                         prot.sphericalAberration.get()),
                                        ("Magnification: ",
                                         prot.magnification.get()),
                                        ("Pixel Size (Å/px): ",
                                         outSet.getSamplingRate())
                                        ]

        self._objects = objects

    def getObjectInfo(self, obj):
        info = {'key': obj.strId(),
                'parent': obj._objParent,
                'text': obj.name,
                'values': (obj.output, obj.outSize),
                'open': True
               }

        return info
