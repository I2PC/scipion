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

import pyworkflow.viewer as pwviewer
from pyworkflow.em.constants import COLOR_OTHER, COLOR_CHOICES
from pyworkflow.em.convert import ImageHandler


class LocalResolutionViewer(pwviewer.ProtocolViewer):
    """
    Visualization tools for local resolution results.

    """
    binaryCondition = ('(colorMap == %d) ' % (COLOR_OTHER))

    def __init__(self, *args, **kwargs):
        pwviewer.ProtocolViewer.__init__(self, *args, **kwargs)

    def getImgData(self, imgFile):
        import numpy as np
        img = ImageHandler().read(imgFile)
        imgData = img.getData()

        maxRes = np.amax(imgData)
        imgData2 = np.ma.masked_where(imgData < 0.1, imgData, copy=True)
        minRes = np.amin(imgData2)

        return imgData2, minRes, maxRes

    def getSlice(self, index, volumeData):
        return int(index*volumeData.shape[0] / 9)

    def getSliceImage(self, volumeData, sliceNumber, dataAxis):
        if dataAxis == 'y':
            imgSlice = volumeData[:, sliceNumber, :]
        elif dataAxis == 'x':
            imgSlice = volumeData[:, :, sliceNumber]
        else:
            imgSlice = volumeData[sliceNumber, :, :]

        return imgSlice

    def createChimeraScript(self, scriptFile, fnResVol, fnOrigMap, sampRate):
        import pyworkflow.gui.plotter as plotter
        import os
        from itertools import izip
        fhCmd = open(scriptFile, 'w')
        imageFile = os.path.abspath(fnResVol)

        _, minRes, maxRes = self.getImgData(imageFile)

        stepColors = self._getStepColors(minRes, maxRes)
        colorList = plotter.getHexColorList(stepColors, self._getColorName())

        fnVol = os.path.abspath(fnOrigMap)

        fhCmd.write("background solid white\n")

        fhCmd.write("open %s\n" % fnVol)
        fhCmd.write("open %s\n" % (imageFile))

        fhCmd.write("volume #0 voxelSize %s\n" % (str(sampRate)))
        fhCmd.write("volume #1 voxelSize %s\n" % (str(sampRate)))
        fhCmd.write("volume #1 hide\n")

        scolorStr = ''
        for step, color in izip(stepColors, colorList):
            scolorStr += '%s,%s:' % (step, color)
        scolorStr = scolorStr[:-1]
        line = ("scolor #0 volume #1 perPixel false cmap " + scolorStr + "\n")
        fhCmd.write(line)

        scolorStr2 = ''
        for step, color in izip(stepColors, colorList):
            indx = stepColors.index(step)
            if ((indx % 4) != 0):
                scolorStr2 += '" " %s ' % color
            else:
                scolorStr2 += '%s %s ' % (step, color)
        line = ("colorkey 0.01,0.05 0.02,0.95 labelColor None "
                + scolorStr2 + " \n")
        fhCmd.write(line)
        fhCmd.close()

    def _getStepColors(self, minRes, maxRes, numberOfColors=13):
        inter = (maxRes - minRes) / (numberOfColors - 1)
        rangeList = []
        for step in range(0, numberOfColors):
            rangeList.append(round(minRes + step * inter, 2))
        return rangeList

    def _getColorName(self):
        if self.colorMap.get() != COLOR_OTHER:
            return COLOR_CHOICES[self.colorMap.get()]
        else:
            return self.otherColorMap.get()


