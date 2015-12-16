# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              C.O.S. Sorzano (coss@cnb.csic.es)
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

import os
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from protocol_structure_mapping import XmippProtStructureMapping
import pylab
import matplotlib.pyplot as plt
import numpy as np
from sklearn import manifold
from mpl_toolkits.mplot3d import proj3d
import pyworkflow.protocol.params as params






class XmippProtStructureMappingViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    
    _label = 'viewer validate_overfitting'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtStructureMapping]
        
    def _defineParams(self, form):
        form.addSection(label='Show StructMap')
        form.addParam('numberOfDimensions', params.IntParam, default=2,
                      label="Number of dimensions",
                      help='In normal practice, it should be 1, 2 or, at most, 3.')
        form.addParam('doShowPlot', params.LabelParam,
                      label="Display the StructMap")
    
    def _getVisualizeDict(self):
        return {'doShowPlot': self._visualize}
    
        
    def _visualize(self, e=None):
        fnOutput = self.protocol._defineResultsName()
        
        if not os.path.exists(fnOutput):
            return [self.errorMessage('The necessary metadata was not produced\n'
                                      'Execute again the protocol\n',
                                      title='Missing result file')]
        
        
        volList = [vol.clone() for vol in self.protocol._iterInputVolumes()]
        
        data = []
        fnFreqOpen = open(fnOutput, "r")
        for line in fnFreqOpen:
            fields = line.split()
            rowdata = map(float, fields)
            data.extend(rowdata)
                
        count = 0
        distance = [[0 for i in volList] for i in volList]
        nVoli = 1
        for i in volList:
            nVolj = 1
            for j in volList:
                distance[(nVoli-1)][(nVolj-1)] = data[count]
                count += 1
                nVolj += 1
            nVoli += 1    
        
        count = 0
        labels = []
        for voli in volList:
            labels.append("%s"%voli.getObjLabel())
            count += 1
        
        nComponent = self.numberOfDimensions.get()
        
        if nComponent == 1:
            mds = manifold.MDS(n_components=1, metric=True, max_iter=3000, eps=1e-9, dissimilarity="precomputed", n_jobs=1)
            embed3d = mds.fit(distance).embedding_ 
            
            AX = plt.subplot(111)
            val = 0
            plot = plt.plot(embed3d[:, 0], np.zeros_like(embed3d[:, 0]) + val, 'o', c='g')
            plt.xlabel('Dimension 1', fontsize=11)
            AX.set_yticks([1])
            plt.title('StructMap')
            
            for label, x, y in zip(labels, embed3d[:, 0], np.zeros_like(embed3d[:, 0]) + val):
                plt.annotate(
                             label, 
                             xy = (x, y), xytext = (-8, 8),
                             textcoords = 'offset points', ha = 'right', va = 'bottom',fontsize=9,
                             bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3))
            plt.show()
                   
        elif nComponent == 2:
            mds = manifold.MDS(n_components=2, metric=True, max_iter=3000, eps=1e-9, dissimilarity="precomputed", n_jobs=1)
            embed3d = mds.fit(distance).embedding_ 
            plot = plt.scatter(embed3d[:, 0], embed3d[:, 1], marker='o', c='g')
            plt.xlabel('Dimension 1', fontsize=11)
            plt.ylabel('Dimension 2', fontsize=11)
            plt.title('StructMap')
            
            for label, x, y in zip(labels, embed3d[:, 0], embed3d[:, 1]):
                plt.annotate(
                            label, 
                            xy = (x, y), xytext = (-8, 8),
                            textcoords = 'offset points', ha = 'right', va = 'bottom',fontsize=9,
                            bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3))
                            #arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            
            plt.grid(True)
            plt.show()
               
        else: 
                        
            fig = pylab.figure()
            ax = fig.add_subplot(111, projection = '3d')
            #ax = plt.subplot(111, projection = '3d')
            mds = manifold.MDS(n_components=3, metric=True, max_iter=3000, eps=1e-9, dissimilarity="precomputed", n_jobs=1)
            embed3d = mds.fit(distance).embedding_ 
            #plot = ax.plot(embed3d[:, 0], embed3d[:, 1], embed3d[:, 2], 'o', c='g')
            ax.scatter(embed3d[:, 0], embed3d[:, 1], embed3d[:, 2], marker = 'o', c='g')
            ax.set_xlabel('Dimension 1', fontsize=11)
            ax.set_ylabel('Dimension 2', fontsize=11)
            ax.set_zlabel('Dimension 3', fontsize=11)
            ax.text2D(0.05, 0.95, "StructMap", transform=ax.transAxes)
                        
            tX, tY, _ = proj3d.proj_transform(embed3d[:, 0], embed3d[:, 1], embed3d[:, 2], ax.get_proj())
            Labels = []
            for i in range(len(embed3d[:, 0])):
                text = labels[i]
                label = ax.annotate(text,
                                    xycoords='data',
                                    xy = (tX[i], tY[i]), xytext = (-8, 8),
                                    textcoords = 'offset points', ha = 'right', va = 'bottom', fontsize=9,
                                    bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5))
                                    
                Labels.append(label)
            def update_position(e):
                x2, y2, _ = proj3d.proj_transform(embed3d[:, 0], embed3d[:, 1], embed3d[:, 2], ax.get_proj())
                for i in range(len(embed3d[:, 0])):
                    label = Labels[i]
                    label.xy = x2[i],y2[i]
                    label.update_positions(fig.canvas.renderer)
                fig.canvas.draw()
            fig.canvas.mpl_connect('button_release_event', update_position)
            pylab.show()
        
        return plot
        
    def _validate(self):
        errors = []
        
        numberOfDimensions=self.numberOfDimensions.get()
        if numberOfDimensions > 3 or numberOfDimensions < 1:
            errors.append("The number of dimensions should be 1, 2 or, at most, 3.")
        
        return errors
    