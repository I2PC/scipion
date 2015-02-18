# **************************************************************************
# *
# * Authors:     L. del Cano (ldelcano@cnb.csic.es)
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
from pyworkflow.em.packages.xmipp3.convert import readSetOfClasses
"""
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from protocol_reconstruct_significant import XmippProtReconstructSignificant
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
from pyworkflow.protocol.params import LabelParam
from convert import readSetOfParticles
import glob

        
class XmippReconstructSignificantViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer reconstruct significant'
    _targets = [XmippProtReconstructSignificant]
    _environments = [DESKTOP_TKINTER]
    
    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('volumeToVisualize', IntParam, default=0,
                      label="Volume to visualize")      
        form.addParam('doShowAngles', LabelParam, default=False,
                      label="Visualize Significant Angular assignment")      
        form.addParam('doShowImagesSignificant', LabelParam, default=False,
                      label="Visualize Optimal Angular assignment")
        if self.protocol.keepIntermediate:
            form.addParam('iterationToVisualize', IntParam, default=-1,
                          label="Iteration to visualize", help="-1 for last iteration")
    
    def _getVisualizeDict(self):
        return {'doShowAngles': self._viewAngles,
                'doShowImagesSignificant': self._viewImagesSignificant}        

    def _getMetadataSqlite(self, fn, clean=False):
        """ Read the classes from Xmipp metadata and write as sqlite file. """
        fnSqlite = fn + '.sqlite'
        if clean:
            cleanPath(fnSqlite)
            
        if not os.path.exists(fnSqlite):
            imagesSet = SetOfParticles(filename=fnSqlite)
            readSetOfParticles(fn+".xmd",imagesSet)
            imagesSet.write()
            imagesSet.close()
        return fnSqlite

    def _getIteration(self, Nvolumes):
        lastIteration=self.protocol.getLastIteration(Nvolumes)
        if self.protocol.keepIntermediate:
            if self.iterationToVisualize==-1:
                return lastIteration
            else:
                if self.iterationToVisualize<=lastIteration:
                    return self.iterationToVisualize.get()
                else:
                    return -1
        else:
            return lastIteration

    def _viewAngles(self, e=None):
        Nvolumes=self.protocol.getNumberOfVolumes()
        if self.volumeToVisualize>=Nvolumes:
            return [self.errorMessage("Requested volume does not exist, the volume number is larger than the number of calculated volumes. First volumes is number 0")]
        iteration=self._getIteration(Nvolumes)
        if iteration!=-1:
            fnAngles=self.protocol._getExtraPath('angles_iter%03d_%02d'%(iteration,self.volumeToVisualize.get()))
            fnSqlite=self._getMetadataSqlite(fnAngles)
            return [ObjectView(self._project, self.protocol.strId(), fnSqlite)]
        else:
            return [self.errorMessage("Requested iteration does not exist, the iteration number is larger than the number of calculated iterations. First iteration is number 0")]

    def _viewImagesSignificant(self, e=None):
        Nvolumes=self.protocol.getNumberOfVolumes()
        if self.volumeToVisualize>=Nvolumes:
            return [self.errorMessage("Requested volume does not exist, the volume number is larger than the number of calculated volumes. First volumes is number 0")]
        iteration=self._getIteration(Nvolumes)
        if iteration!=-1:
            fnImages=self.protocol._getExtraPath('images_significant_iter%03d_%02d'%(iteration,self.volumeToVisualize.get()))
            fnSqlite=self._getMetadataSqlite(fnImages)
            return [ObjectView(self._project, self.protocol.strId(), fnSqlite)]
        else:
            return [self.errorMessage("Requested iteration does not exist, the iteration number is larger than the number of calculated iterations. First iteration is number 0")]
