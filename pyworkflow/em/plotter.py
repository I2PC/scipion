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
from pyworkflow.gui.plotter import Plotter

class EmPlotter(Plotter):
    ''' Class to create several plots'''
    
    def plotAngularDistribution(self, title, rot, tilt, weight=[], max_p=40, min_p=5, max_w=2, min_w=1, color='blue'):
        '''Create an special type of subplot, representing the angular
        distribution of weight projections. '''
#        from math import radians
        
        if (len(weight) > 0):
            max_w = max(weight)
            min_w = min(weight)
            a = self.createSubPlot(title, 'Min weight=%(min_w)f, Max weight=%(max_w)f' % locals(), '', projection='polar')
        else:
            a = self.createSubPlot(title, 'Empty plot', '', projection='polar')
      
        for i, rotItem in enumerate(rot):
            if (len(weight) > 0):
                pointsize = int((weight[i] - min_w)/(max_w - min_w + 0.001) * (max_p - min_p) + min_p)
            else:
                pointsize = 1
            a.plot(rot[i], tilt[i], markerfacecolor=color, marker='.', markersize=pointsize)
    
    def plotHist(self, xCol, yCol, color, **args):
        """ plot columns xCol and yCol
            if nbins is in args then and histogram over y data is made
        """
        nbins = args.pop('nbins', None)
        marker = args.pop('marker', None)
        linestyle = args.pop('linestyle', None)
        if nbins is None:
            if not marker is None:
                args['marker'] = marker     
            if not linestyle is None:
                args['linestyle'] = linestyle
            self.plot(xCol, yCol, color, **args) #no histogram
        else:
            self.hist(yCol,nbins, facecolor=color, **args)
    
