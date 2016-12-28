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
from pyworkflow.gui.plotter import plt
import numpy as np
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
        fnOutput1 = self.protocol._defineResultsName1()
        fnOutput2 = self.protocol._defineResultsName2()
        fnOutput3 = self.protocol._defineResultsName3()
        
        if not os.path.exists(fnOutput3):
            return [self.errorMessage('The necessary metadata was not produced\n'
                                      'Execute again the protocol\n',
                                      title='Missing result file')]
        
        
        volList = [vol.clone() for vol in self.protocol._iterInputVolumes()]
        
        # 1 Dimension
        data1 = []
        fnCoordinate1 = open(fnOutput1, "r")
        
        for line in fnCoordinate1:
            fields = line.split()
            rowdata = map(float, fields)
            data1.extend(rowdata)
        
        data1N = [[0 for i in range(1)] for i in volList]
        count = 0
        
        for i in volList:
            data1N[count][0] = data1[count]
            count += 1
        data1N = np.array (data1N)    
        
        # 2 Dimensions      
        data2 = []
        fnCoordinate2 = open(fnOutput2, "r")
        
        for line in fnCoordinate2:
            fields = line.split()
            rowdata = map(float, fields)
            data2.extend(rowdata)
            
        count = 0
        data2N = [[0 for i in range(2)] for i in volList]
        nVolj = 1
        
        for j in range(2):
            nVoli = 1
            for i in volList:
                data2N[(nVoli-1)][(nVolj-1)] = data2[count]
                count += 1
                nVoli += 1
            nVolj += 1 
        data2N = np.array (data2N)   
                 
        # 3 Dimensions    
        data3 = []
        fnCoordinate3 = open(fnOutput3, "r")
        
        for line in fnCoordinate3:
            fields = line.split()
            rowdata = map(float, fields)
            data3.extend(rowdata)
       
        
        count = 0
        data3N = [[0 for i in range(3)] for i in volList]
        nVolj = 1
        
        for j in range(3):
            nVoli = 1
            for i in volList:
                data3N[(nVoli-1)][(nVolj-1)] = data3[count]
                count += 1
                nVoli += 1
            nVolj += 1 
        data3N = np.array (data3N)
                        
        count = 0
        labels = []
        
        for voli in volList:
            labels.append("%s"%voli.getObjLabel())
            count += 1
        
        nComponent = self.numberOfDimensions.get()
        
        if nComponent == 1:
            AX = plt.subplot(111)
            val = 0
            plot = plt.plot(data1N[:, 0], np.zeros_like(data1N[:, 0]) + val, 'o', c='g')
            plt.xlabel('Dimension 1', fontsize=11)
            AX.set_yticks([1])
            plt.title('StructMap')
            
            for label, x, y in zip(labels, data1N[:, 0], np.zeros_like(data1N[:, 0 ]) + val):
                plt.annotate(
                             label, 
                             xy = (x, y), xytext = (-8, 8),
                             textcoords = 'offset points', ha = 'right', va = 'bottom',fontsize=9,
                             bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3))
            plt.show()
                        
        elif nComponent == 2:
            plot = plt.scatter(data2N[:, 0], data2N[:, 1], marker='o', c='g')
            plt.xlabel('Dimension 1', fontsize=11)
            plt.ylabel('Dimension 2', fontsize=11)
            plt.title('StructMap')
            
            for label, x, y in zip(labels, data2N[:, 0], data2N[:, 1]):
                plt.annotate(
                            label, 
                            xy = (x, y), xytext = (-8, 8),
                            textcoords = 'offset points', ha = 'right', va = 'bottom',fontsize=9,
                            bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3))
                                        
            plt.grid(True)
            plt.show()
               
        else: 
                        
            fig = plt.figure()
            ax = fig.add_subplot(111, projection = '3d')
           
            ax.scatter(data3N[:, 0], data3N[:, 1], data3N[:, 2], marker = 'o', c='g')
            ax.set_xlabel('Dimension 1', fontsize=11)
            ax.set_ylabel('Dimension 2', fontsize=11)
            ax.set_zlabel('Dimension 3', fontsize=11)
            ax.text2D(0.05, 0.95, "StructMap", transform=ax.transAxes)
                        
            tX, tY, _ = proj3d.proj_transform(data3N[:, 0], data3N[:, 1],
                                               data3N[:, 2], ax.get_proj())
            Labels = []
            
            for i in range(len(data3N[:, 0])):
                text = labels[i]
                label = ax.annotate(text,
                                    xycoords='data',
                                    xy = (tX[i], tY[i]), xytext = (-8, 8),
                                    textcoords = 'offset points', ha = 'right',
                                     va = 'bottom', fontsize=9,
                                    bbox = dict(boxstyle = 'round,pad=0.5',
                                                 fc = 'yellow', alpha = 0.5))
                                    
                Labels.append(label)
            def update_position(e):
                x2, y2, _ = proj3d.proj_transform(data3N[:, 0], data3N[:, 1],
                                                  data3N[:, 2], ax.get_proj())
                
                for i in range(len(data3N[:, 0])):
                    label = Labels[i]
                    label.xy = x2[i],y2[i]
                    label.update_positions(fig.canvas.renderer)
                fig.canvas.draw()
            fig.canvas.mpl_connect('button_release_event', update_position)
            plt.show()
        
        return plot
        
    def _validate(self):
        errors = []
        
        numberOfDimensions=self.numberOfDimensions.get()
        if numberOfDimensions > 3 or numberOfDimensions < 1:
            errors.append("The number of dimensions should be 1, 2 or, at most, 3.")
        
        return errors
    