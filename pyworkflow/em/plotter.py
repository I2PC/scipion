# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This module implement the classes to create plots on xmipp.
"""
from itertools import izip
import matplotlib.pyplot as plt 

import pyworkflow.em.metadata as md
from pyworkflow.gui.plotter import Plotter



class EmPlotter(Plotter):
    ''' Class to create several plots'''
    def __init__(self, x=1, y=1, mainTitle="", **kwargs):
        Plotter.__init__(self, x, y, mainTitle, **kwargs)

    def plotAngularDistribution(self, title, rot, 
                                tilt, weight=[], max_p=40, 
                                min_p=5, max_w=2, min_w=1, color='blue'):
        '''Create an special type of subplot, representing the angular
        distribution of weight projections. '''
        if weight:
            max_w = max(weight)
            min_w = min(weight)
            a = self.createSubPlot(title, 'Min weight=%(min_w).2f, Max weight=%(max_w).2f' % locals(), '', projection='polar')
            for r, t, w in izip(rot, tilt, weight):
                pointsize = int((w - min_w)/(max_w - min_w + 0.001) * (max_p - min_p) + min_p)
                a.plot(r, t, markerfacecolor=color, marker='.', markersize=pointsize)
        else:
            a = self.createSubPlot(title, 'Empty plot', '', projection='polar')
            for r, t in izip(rot, tilt):
                a.plot(r, t, markerfacecolor=color, marker='.', markersize=10)
                
    def plotAngularDistributionFromMd(self, mdFile, title, **kwargs):
        """ Read the values of rot, tilt and weights from
        the medata and plot the angular distribution.
        In the metadata:
            rot: MDL_ANGLE_ROT
            tilt: MDL_ANGLE_TILT
            weight: MDL_WEIGHT
        """
        angMd = md.MetaData(mdFile)
        rot = []
        tilt = []
        weight = []
        
        for row in md.iterRows(angMd):
            rot.append(row.getValue(md.MDL_ANGLE_ROT))
            tilt.append(row.getValue(md.MDL_ANGLE_TILT))
            weight.append(row.getValue(md.MDL_WEIGHT))
            
        return self.plotAngularDistribution(title, rot, tilt, weight, **kwargs)
        
    def plotHist(self, yValues, nbins, color='blue', **kwargs):
        """ Create an histogram. """
        self.hist(yValues, nbins, facecolor=color, **kwargs)
        
    def plotMatrix(self,_matrix,cmap='Greens'
                       , xticksLablesMajor=None
                       , yticksLablesMajor=None
                       , rotationX=90.
                       , rotationY=0.):
        im = plt.imshow(_matrix, interpolation="none", cmap=cmap)
        if (xticksLablesMajor is not None):       
            plt.xticks(range(len(xticksLablesMajor)), 
                                 xticksLablesMajor[:len(xticksLablesMajor)],
                                 rotation=rotationX)
        if (yticksLablesMajor is not None):       
            plt.yticks(range(len(yticksLablesMajor)),
                                 yticksLablesMajor[:len(yticksLablesMajor)],
                                 rotation=rotationY)
        cax = plt.colorbar(im)
        #im.cmap.set_over('g')#outbound values

    def plotData(self, xValues, yValues, color='blue', **kwargs):
        """ Shortcut function to plot some values.
        Params:
            xValues: list of values to show in x-axis
            yValues: list of values to show as values in y-axis
            color: color for the plot.
            **kwargs: keyword arguments that accepts:
                marker, linestyle
        """

        self.plot(xValues, yValues, color, **kwargs)