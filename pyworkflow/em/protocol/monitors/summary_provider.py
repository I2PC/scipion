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

import datetime

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

        def addObj(objId, name, output='', size='', rate='', trend='', parent=None):
            if objId not in objIds:
                obj = pwobj.Object(objId=objId)
                obj.name = name
                obj.output = output
                obj.outSize = size
                obj.rate = rate
                obj.trend = trend
                obj._objParent = parent
                objIds.append(objId)
                objects.append(obj)
                return obj
            else:
                return None

        prots = [getUpdatedProtocol(p) for p in self.protocol.getInputProtocols()]

        for prot in prots:
            pobj = addObj(prot.getObjId(),
                          '%s (id=%s)' % (prot.getRunName(), prot.strId()))
            for outName, outSet in prot.iterOutputAttributes(pwobj.Set):
                outSet.load()
                outSet.loadAllProperties()
                # outSetId needs to be compound id to avoid duplicate ids
                outSetId = '%s.%s' % (outSet.getObjId(), prot.getObjId())
                rate, trend = self.getOutputStats(outSet)
                addObj(outSetId, '', outName, outSet.getSize(), rate, trend, pobj)
                outSet.close()
                # Store acquisition parameters in case of the import protocol
                from pyworkflow.em import ProtImportImages
                # NOTE by Yaiza: we force the string containing the Å to be unicode
                # because this is the encoding used when generating report in report_html.py
                if isinstance(prot, ProtImportImages):
                    self.acquisition = [("Microscope Voltage (kV): ",
                                         prot.voltage.get()),
                                        ("Spherical aberration (mm): ",
                                         prot.sphericalAberration.get()),
                                        ("Magnification: ",
                                         prot.magnification.get()),
                                        (u"Pixel Size (Å/px): ",
                                         round(outSet.getSamplingRate(),2))
                                        ]
                    if prot.dosePerFrame.get() is not None:
                        self.acquisition.append((u"Dose per frame (e/Å²):",
                                                 prot.dosePerFrame.get()))

        self._objects = objects

    def getOutputStats(self, outSet):

        TREND_ICONS_ROOT_URL = "https://raw.githubusercontent.com/I2PC/scipion/dm_ratesHTMLreport/pyworkflow/resources/"
        MORE_INCREASING_TREND_ICON = TREND_ICONS_ROOT_URL + "Trend-MoreIncreasing.png"
        INCREASING_TREND_ICON = TREND_ICONS_ROOT_URL + "Trend-Increasing.png"
        STABLE_TREND_ICON = TREND_ICONS_ROOT_URL + "Trend-Stable.png"
        DECREASING_TREND_ICON = TREND_ICONS_ROOT_URL + "Trend-Decreasing.png"
        MORE_DECREASING_TREND_ICON = TREND_ICONS_ROOT_URL + "Trend-MoreDecreasing.png"
        STOPPED_TREND_ICON = TREND_ICONS_ROOT_URL + "Trend-Stopped.png"
        NODATA_TREND_ICON = TREND_ICONS_ROOT_URL + "Trend-NoData.png"

        DELTA_TIME = 0.1  # in hours (it can be a fraction)

        time0 = outSet.getFirstItem().getObjCreationAsDate()
        finished = False
        if outSet.isStreamClosed():
            timeF = datetime.datetime.utcnow()  # just in case that is empty
            for p in outSet.iterItems(orderBy='creation',
                                      direction='DESC'):
                timeF = p.getObjCreationAsDate()
                break
            finished = True
        else:
            timeF = datetime.datetime.utcnow()

        time = (timeF-time0).total_seconds() / 3600.0  # in hours

        rate = outSet.getSize() / time  # items/hour

        if time < DELTA_TIME*2 or finished:
            return '%.0f' % rate, NODATA_TREND_ICON

        timeI = timeF - datetime.timedelta(hours=DELTA_TIME)
        timeIstr = timeI.strftime(pwobj.String.DATETIME_FORMAT)
        # This loop is just over the last minutes to count items
        whereStr = 'creation>"%s"' % timeIstr
        bunchSize = len([x for x in outSet.iterItems(orderBy='creation',
                                                     direction='DESC',
                                                     where=whereStr)])

        timeBunch = (timeF - timeI).total_seconds() / 3600.0  # in hours
        newRate = bunchSize/timeBunch

        trendValue = newRate/rate  # maybe it could be (newRate-rate)/rate

        highestThres = 2
        highThres = 1.3
        lowThres = 0.8
        lowestThres = 0.3

        if trendValue > highestThres:
            trend = MORE_INCREASING_TREND_ICON
        elif trendValue > highThres:
            trend = INCREASING_TREND_ICON
        elif trendValue > lowThres:
            trend = STABLE_TREND_ICON
        elif trendValue > lowestThres:
            trend = DECREASING_TREND_ICON
        elif trendValue > 0:
            trend = MORE_DECREASING_TREND_ICON
        else:
            trend = STOPPED_TREND_ICON

        return '%.0f' % rate, trend

    def getObjectInfo(self, obj):
        info = {'key': obj.strId(),
                'parent': obj._objParent,
                'text': obj.name,
                'values': (obj.output, obj.outSize),
                'open': True
               }

        return info
