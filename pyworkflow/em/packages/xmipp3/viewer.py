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
from pyworkflow.em.data_tiltpairs import MicrographsTiltPair, ParticlesTiltPair, CoordinatesTiltPair
from convert import *
from os.path import dirname, join
from pyworkflow.utils import makePath, runJob, copyTree
import pyworkflow as pw
import xmipp

from protocol_cl2d_align import XmippProtCL2DAlign
from protocol_cl2d import XmippProtCL2D
from protocol_convert_to_pseudoatoms import XmippProtConvertToPseudoAtoms
from protocol_ctf_discrepancy import XmippProtCTFDiscrepancy
from protocol_extract_particles import XmippProtExtractParticles
from protocol_extract_particles_pairs import XmippProtExtractParticlesPairs
from protocol_helical_parameters import XmippProtHelicalParameters
from protocol_identify_outliers import XmippProtIdentifyOutliers
from protocol_kerdensom import XmippProtKerdensom
from protocol_particle_pick_automatic import XmippParticlePickingAutomatic
from protocol_particle_pick import XmippProtParticlePicking
from protocol_particle_pick_pairs import XmippProtParticlePickingPairs
from protocol_preprocess import XmippProtPreprocessVolumes
from protocol_preprocess_micrographs import XmippProtPreprocessMicrographs
from protocol_rotational_spectra import XmippProtRotSpectra
from protocol_screen_classes import XmippProtScreenClasses
from protocol_screen_particles import XmippProtScreenParticles


class XmippViewer(Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [  
                  CoordinatesTiltPair 
                , Image
                , MicrographsTiltPair
                , NormalModes
                , ParticlesTiltPair
                , ProtExtractParticles
                , SetOfClasses2D
                , SetOfClasses3D
                , SetOfCoordinates
                , SetOfCTF
                , SetOfImages
                , SetOfMovies
                , XmippParticlePickingAutomatic
                , XmippProtConvertToPseudoAtoms
                , XmippProtExtractParticlesPairs
                , XmippProtIdentifyOutliers
                , XmippProtKerdensom
                , XmippProtParticlePicking
                , XmippProtParticlePickingPairs
                , XmippProtRotSpectra
                , XmippProtScreenClasses
                , XmippProtScreenParticles
                ]
    
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
            self._views.append(DataView(fn, viewParams={RENDER: 'image'}))
                 
        elif issubclass(cls, Image):
            fn = getImageLocation(obj)
            self._views.append(DataView(fn))
            
        elif issubclass(cls, NormalModes):
            self._views.append(DataView(fn))
              
        elif issubclass(cls, SetOfMicrographs):
            
            fn = obj.getFileName()
            self._views.append(ObjectView(self._project.getName(), obj.strId(), fn, viewParams={MODE: MODE_MD}, **args))
            
        elif issubclass(cls, SetOfMovies):
            fn = self._getTmpPath(obj.getName() + '_movies.xmd')
            writeSetOfMovies(obj, fn)
            self._views.append(ObjectView(self._project.getName(), obj.strId(), fn, **args))    


        elif issubclass(cls, MicrographsTiltPair):          
#             fnU = obj.getUntilted().getFileName()
#             fnT = obj.getTilted().getFileName()
#             self._views.append(ObjectView(self._project.getName(), obj.strId(), fnU, **args))            
#             self._views.append(ObjectView(self._project.getName(), obj.strId(), fnT, **args))
            labels = 'id enabled _untilted._filename _tilted._filename'
            self._views.append(ObjectView(self._project.getName(), obj.strId(), obj.getFileName(), 
                                          viewParams={ORDER: labels, 
                                                      VISIBLE: labels, 
                                                      MODE: MODE_MD}))
            
        elif issubclass(cls, ParticlesTiltPair):          
            labels = 'id enabled _untilted._filename _tilted._filename'
            self._views.append(ObjectView(self._project.getName(), obj.strId(), obj.getFileName(), 
                                          viewParams={ORDER: labels, 
                                                      VISIBLE: labels, RENDER:'_untilted._filename _tilted._filename',
                                                      MODE: MODE_MD}))

                            
        elif issubclass(cls, SetOfCoordinates):
            micSet = obj.getMicrographs()  # accessing mics to provide metadata file
            if micSet is None:
                raise Exception('visualize: SetOfCoordinates has no micrographs set.')
            
            mdFn = getattr(micSet, '_xmippMd', None)
            tmpDir = self._getTmpPath(obj.getName()) 
            makePath(tmpDir)
            
            if mdFn:
                fn = mdFn.get()
            else:  # happens if protocol is not an xmipp one
                fn = self._getTmpPath(micSet.getName() + '_micrographs.xmd')
                writeSetOfMicrographs(micSet, fn)
            posDir = getattr(obj, '_xmippMd', None)  # extra dir istead of md file for SetOfCoordinates
            if posDir:
                copyTree(posDir.get(), tmpDir)
            else:
                writeSetOfCoordinates(tmpDir, obj)   
            self._views.append(CoordinatesObjectView(fn, tmpDir))

        elif issubclass(cls, SetOfParticles):
            fn = obj.getFileName()
            labels = 'id enabled _index _filename _xmipp_zScore _sampling '
            labels += '_ctfModel._defocusU _ctfModel._defocusV _ctfModel._defocusAngle'
            self._views.append(ObjectView(self._project.getName(), obj.strId(), fn,
                                          viewParams={ORDER: labels, 
                                                      VISIBLE: labels, 
                                                      'sortby': '_ctfModel._xmipp_zScore desc', RENDER:'_filename'}))
               
        elif issubclass(cls, SetOfVolumes):
            fn = obj.getFileName()
            self._views.append(ObjectView(self._project.getName(), obj.strId(), fn))
        
        elif issubclass(cls, SetOfClasses2D):
            fn = obj.getFileName()
            self._views.append(ClassesView(self._project.getName(), obj.strId(), fn))
            
        elif issubclass(cls, SetOfClasses3D):
            fn = obj.getFileName()
            self._views.append(Classes3DView(self._project.getName(), obj.strId(), fn))            
              
        elif issubclass(cls, SetOfCTF):
            fn = obj.getFileName()
#            self._views.append(DataView(fn, viewParams={MODE: 'metadata'}))
            psdLabels = '_psdFile _xmipp_enhanced_psd _xmipp_ctfmodel_quadrant _xmipp_ctfmodel_halfplane'
            labels = 'id enabled comment %s _defocusU _defocusV _defocusAngle _xmippCTFCritFirstZero _micObj._filename' % psdLabels 
            self._views.append(ObjectView(self._project.getName(), obj.strId(), fn, 
                                          viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels, ZOOM: 50, RENDER: psdLabels}))    

        elif issubclass(cls, CoordinatesTiltPair):
            tmpDir = self._getTmpPath(obj.getName()) 
            makePath(tmpDir)

            mdFn = join(tmpDir, 'input_micrographs.xmd')
            writeSetOfMicrographsPairs(obj.getUntilted().getMicrographs(),
                                        obj.getTilted().getMicrographs(), 
                                        mdFn) 
            parentProtId = obj.getObjParentId()
            parentProt = self.getProject().mapper.selectById(parentProtId)
            extraDir = parentProt._getExtraPath()
            
            extraDir = parentProt._getExtraPath()
            #TODO: Review this if ever a non Xmipp CoordinatesTiltPair is available
#             writeSetOfCoordinates(tmpDir, obj.getUntilted()) 
#             writeSetOfCoordinates(tmpDir, obj.getTilted()) 
            
            scipion =  "%s %s \"%s\" %s" % ( pw.PYTHON, pw.join('apps', 'pw_create_coords.py'), self.getProject().getDbPath(), obj.strId())
            app = "xmipp.viewer.particlepicker.tiltpair.TiltPairPickerRunner"
            args = " --input %(mdFn)s --output %(extraDir)s --mode readonly --scipion %(scipion)s"%locals()
        
            runJavaIJapp("%dg"%(2), app, args, True)
         
        elif issubclass(cls, XmippProtExtractParticles) or issubclass(cls, XmippProtScreenParticles):
            particles = obj.outputParticles
            self._visualize(particles)
            
            fn = obj._getPath('images.xmd')
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
            
        elif issubclass(cls, XmippProtParticlePicking):
            self._visualize(obj.getCoords())
            
        elif issubclass(cls, XmippParticlePickingAutomatic):
        	micSet = obj.getInputMicrographs()
        	mdFn = getattr(micSet, '_xmippMd', None)
        	if mdFn:
        		micsfn = mdFn.get()
        	else:  # happens if protocol is not an xmipp one
        		micsfn = self._getTmpPath(micSet.getName() + '_micrographs.xmd')
                writeSetOfMicrographs(micSet, micsfn)
                
           	posDir = getattr(obj.getCoords(), '_xmippMd').get()  # extra dir istead of md file for SetOfCoordinates
           	scipion = "%s %s \"%s\" %s" % (pw.PYTHON, pw.join('apps', 'pw_create_coords.py'), self.getProject().getDbPath(), obj.strId())
        	app = "xmipp.viewer.particlepicker.training.SupervisedPickerRunner"# Launch the particle picking GUI
        	args = "--input %(micsfn)s --output %(posDir)s --mode review  --scipion %(scipion)s" % locals()
        	runJavaIJapp("%dg" % (2), app, args, True)
            # self._views.append(CoordinatesObjectView(fn, tmpDir, 'review', self._project.getName(), obj.strId()))

        elif issubclass(cls, XmippProtParticlePickingPairs):
            tmpDir = self._getTmpPath(obj.getName()) 
            makePath(tmpDir)

            mdFn = join(tmpDir, 'input_micrographs.xmd')
            writeSetOfMicrographsPairs(obj.outputCoordinatesTiltPair.getUntilted().getMicrographs(),
                                        obj.outputCoordinatesTiltPair.getTilted().getMicrographs(), 
                                        mdFn) 
            extraDir = obj._getExtraPath()
            scipion =  "%s %s \"%s\" %s" % ( pw.PYTHON, pw.join('apps', 'pw_create_coords.py'), self.getProject().getDbPath(), obj.strId())
            app = "xmipp.viewer.particlepicker.tiltpair.TiltPairPickerRunner"
            args = " --input %(mdFn)s --output %(extraDir)s --mode readonly --scipion %(scipion)s"%locals()
        
            runJavaIJapp("%dg"%(obj.memory.get()), app, args, True)
            
        elif issubclass(cls, XmippProtExtractParticlesPairs):
            self._visualize(obj.outputParticlesTiltPair)
            
        return self._views
        
