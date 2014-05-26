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
This module implement the wrappers around xmipp_showj
visualization program.
"""

import os
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.data import *
from pyworkflow.em.protocol import *
from xmipp3 import getXmippPath
from protocol_preprocess_micrographs import XmippProtPreprocessMicrographs
from protocol_extract_particles import XmippProtExtractParticles, ProtImportParticles
from protocol_screen_particles import XmippProtScreenParticles
from protocol_cl2d_align import XmippProtCL2DAlign
from protocol_cl2d import XmippProtCL2D
from protocol_kerdensom import XmippProtKerdensom
from protocol_rotational_spectra import XmippProtRotSpectra
from protocol_screen_classes import XmippProtScreenClasses
from protocol_helical_parameters import XmippProtHelicalParameters
from protocol_convert_to_pseudoatoms import XmippProtConvertToPseudoAtoms
from protocol_identify_outliers import XmippProtIdentifyOutliers
from protocol_preprocess import XmippProtPreprocessVolumes
from convert import *
from os.path import dirname, join
from pyworkflow.utils import makePath, runJob, copyTree
import pyworkflow as pw
import xmipp



class XmippViewer(Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [Image, SetOfImages, SetOfCoordinates, SetOfClasses2D, SetOfClasses3D, 
                SetOfMovies, ProtExtractParticles, XmippProtScreenParticles, 
                XmippProtKerdensom, XmippProtRotSpectra,  
                SetOfCTF, NormalModes, XmippProtScreenClasses,
                XmippProtConvertToPseudoAtoms, XmippProtIdentifyOutliers]
    
    def __init__(self, **args):
        Viewer.__init__(self, **args)
        self._views = []   
        
        
    def visualize(self, obj, **args):
        self._visualize(obj, **args)
        
        for v in self._views:
            v.show()
            
    def _visualize(self, obj, **args):
        cls = type(obj)

        if issubclass(cls, Volume):
            fn = getImageLocation(obj)
            if fn.endswith('.mrc'):
                fn += ":mrc"
            self._views.append(DataView(fn))
                 
        elif issubclass(cls, Image):
            fn = getImageLocation(obj)
            self._views.append(DataView(fn))
            
        elif issubclass(cls, NormalModes):
            self._views.append(DataView(fn))
              
        elif issubclass(cls, SetOfMicrographs):
            
            mdFn = getattr(obj, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:
                fn = obj.getFileName()
            self._views.append(ObjectView(fn, "Micrographs", self._project.getName(), obj.strId(), obj.strId(), **args))
            
        elif issubclass(cls, SetOfMovies):
            fn = self._getTmpPath(obj.getName() + '_movies.xmd')
            writeSetOfMovies(obj, fn)
            self._views.append(ObjectView(fn, "Micrographs", self._project.getName(), obj.strId(), obj.strId(), **args))    
            
                
        elif issubclass(cls, SetOfCoordinates):
            micSet = obj.getMicrographs()#accessing mics to provide metadata file
            if micSet is None:
                raise Exception('visualize: SetOfCoordinates has no micrographs set.')
            
            mdFn = getattr(micSet, '_xmippMd', None)
            tmpDir = self._getTmpPath(obj.getName()) 
            makePath(tmpDir)
            
            if mdFn:
                fn = mdFn.get()
            else: # happens if protocol is not an xmipp one
                fn = self._getTmpPath(micSet.getName() + '_micrographs.xmd')
                writeSetOfMicrographs(micSet, fn)
            posDir = getattr(obj, '_xmippMd', None)#extra dir istead of md file for SetOfCoordinates
            if posDir:
                copyTree(posDir.get(), tmpDir)
            else:
                writeSetOfCoordinates(tmpDir, obj)   
                           
            self._views.append(CoordinatesObjectView(fn, tmpDir, 'review', self._project.getName(), obj.strId()))
        
        elif issubclass(cls, SetOfParticles):
#            mdFn = getattr(obj, '_xmippMd', None)
#            if mdFn:
#                fn = mdFn.get()
#
#            else:
            fn = obj.getFileName()
            self._views.append(ObjectView(fn, "Particles", self._project.getName(), obj.strId(), obj.strId()))
               
                    
        elif issubclass(cls, SetOfVolumes):
            mdFn = getattr(obj, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:
                fn = obj.getFileName()
            self._views.append(DataView(fn))
        
        elif issubclass(cls, SetOfClasses2D):
#            mdFn = getattr(obj, '_xmippMd', None)
#            if mdFn:
#                fn = mdFn.get()
#            else:
            fn = obj.getFileName()
            self._views.append(ObjectView(fn, "Classes2D", self._project.getName(), obj.strId(), obj.getImages().strId()))  
            
        elif issubclass(cls, SetOfClasses3D):
#            mdFn = getattr(obj, '_xmippMd', None)
#            if mdFn:
#                fn = mdFn.get()
#            else:
            fn = obj.getFileName()

            self._views.append(ObjectView(fn, "Classes3D", self._project.getName(), obj.strId(), obj.getImages().strId(), extraParams=args.get('extraParams', '')))
              
        elif issubclass(cls, SetOfCTF):
#            mdFn = getattr(obj, '_xmippMd', None)
#            if mdFn:
#                fn = mdFn.get()
##            else:
##                fn = self._getTmpPath(obj.getName() + '_ctfs.xmd')
##                writeSetOfCTFs(obj, fn)
#            else:
            fn = obj.getFileName()
            self._views.append(DataView(fn, viewParams={MODE: 'metadata', RENDER: 'psd psdEnhanced image1 image2', ORDER:'id psd psdEnhanced image1 image2', ZOOM: 50}))  
         
        elif issubclass(cls, XmippProtExtractParticles) or issubclass(cls, XmippProtScreenParticles):
            self._visualize(obj.outputParticles)
            fn = getattr(obj.outputParticles, '_xmippMd', None)
            if fn:
                md = xmipp.MetaData(fn) 
                # If Zscore on output images plot Zscore particle sorting
                if md.containsLabel(xmipp.MDL_ZSCORE):
                    from plotter import XmippPlotter
                    xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
                    xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
                    xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_ZSCORE)
                    self._views.append(xplotter)
    
        elif issubclass(cls, XmippProtRotSpectra):
            self._visualize(obj.outputClasses, extraParams='--mode rotspectra --columns %d' % obj.SomXdim.get())
        
        elif issubclass(cls, XmippProtKerdensom):
            self._visualize(obj.outputClasses, extraParams='--columns %d' % obj.SomXdim.get())

        elif issubclass(cls, XmippProtScreenClasses) or issubclass(cls, XmippProtIdentifyOutliers):
            self._views.append(DataView(obj.getVisualizeInfo().get(), viewParams={'mode': 'metadata'}))

        elif issubclass(cls, XmippProtConvertToPseudoAtoms):
            self._views.append(CommandView(obj._getPath('chimera.cmd')))

        return self._views
        
