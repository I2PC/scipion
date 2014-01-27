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
This module implement the classes to create plots on xmipp.
"""
#from pyworkflow.gui.plotter import Plotter
from pyworkflow.em.plotter import EmPlotter

class XmippPlotter(EmPlotter):
    ''' Class to create several plots with Xmipp utilities'''
    
    def plotMdAngularDistribution(self, title, md, color='blue'):
        '''Create an special type of subplot, representing the angular
        distribution of weight projections. A metadata should be provided containing
        labels: MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_WEIGHT '''
        from math import radians
        from xmipp import MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_WEIGHT
        
        rot = [radians(md.getValue(MDL_ANGLE_ROT, objId)) for objId in md]
        tilt = [md.getValue(MDL_ANGLE_TILT, objId) for objId in md]
        weight = [md.getValue(MDL_WEIGHT, objId) for objId in md]
        
        self.plotAngularDistribution(title, rot, tilt, weight)
    
    def plotMd(self, md, mdLabelX, mdLabelY, color='g',**args):
        """ plot metadata columns mdLabelX and mdLabelY
            if nbins is in args then and histogram over y data is made
        """
        if mdLabelX:
            xx = []
        else:
            xx = range(1, md.size() + 1)
        yy = []
        for objId in md:
            if mdLabelX:
                xx.append(md.getValue(mdLabelX, objId))
            yy.append(md.getValue(mdLabelY, objId))
        
        self.plotHist(xx, yy, color, **args)
        
    def plotMdFile(self, mdFilename, mdLabelX, mdLabelY, color='g', **args):
        """ plot metadataFile columns mdLabelX and mdLabelY
            if nbins is in args then and histogram over y data is made
        """
        from xmipp import MetaData
        md = MetaData(mdFilename)
        self.plotMd(md, mdLabelX, mdLabelY, color='g',**args)
      
