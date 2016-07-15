# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
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

import math
from glob import glob
import numpy

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import atsas
from pyworkflow.utils.path import createLink

# TODO: Move to 3D Tools
class AtsasProtConvertPdbToSAXS(ProtPreprocessVolumes):
    """ Protocol for converting a PDB file (true atoms or pseudoatoms) into a
    SAXS curve.
    
    This is actually a wrapper to the program Crysol from Atsas.
    See documentation at:
       http://www.embl-hamburg.de/biosaxs/manuals/crysol.html
    """
    _label = 'convert PDB to SAXS curve'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass='PdbFile',
                      label="Input structure", important=True)
        form.addParam('numberOfSamples', IntParam, default=256,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of samples',
                      help='Number of samples of the curve (parameter ns)') 
        form.addParam('maximumFrequency', FloatParam, default=0.3,
                      expertLevel=LEVEL_ADVANCED,
                      label='Maximum frequency (1/A)',
                      help='Maximum frequency of the curve (parameter sm)') 
        form.addParam('numberOfHarmonics', IntParam, default=18,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of harmonics',
                      help='Number of harmonics to generate the curve '
                           '(parameter lm)')
        form.addParam('experimentalSAXS', FileParam, filter="*.dat", default='',
                      label='Experimental SAXS curve (optional)',
                      help="This parameter is optional. If it is given the "
                           "simulated SAXS curve will be compared to the "
                           "experimental one")
        form.addParam('otherCrysol', StringParam, default='', 
                      label='Other parameters for Crysol',
                      help='See http://www.embl-hamburg.de/biosaxs/manuals/crysol.html') 
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('crysolWrapper')
        
    #--------------------------- STEPS functions -------------------------------
    def crysolWrapper(self):
        experimentalSAXS=""
        if self.experimentalSAXS.get()!="":
            experimentalSAXS=os.path.abspath(self.experimentalSAXS.get())
        inputStructure=os.path.abspath(self.inputStructure.get().getFileName())
        self._enterWorkingDir()
        if experimentalSAXS!="":
            createLink(experimentalSAXS,'experimental_SAXS_curve.dat')
            experimentalSAXS='experimental_SAXS_curve.dat'
        createLink(inputStructure,'pseudoatoms.pdb')
        self.runJob("crysol",
                    "pseudoatoms.pdb %s /lm %d /sm %f /ns %d %s"
                    %(experimentalSAXS, self.numberOfHarmonics,
                      self.maximumFrequency, self.numberOfSamples,
                      self.otherCrysol))
        self.runJob("mv","*log *txt extra")
        if experimentalSAXS=="":
            self.runJob("mv","*alm extra")
        self._leaveWorkingDir()       
        
    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append('Number of samples: %d' % self.numberOfSamples)
        summary.append('Maximum frequency: %d' % self.maximumFrequency)
        summary.append('Number of Harmonics: %d' % self.numberOfHarmonics)

        if not self.experimentalSAXS.empty():
            summary.append('Experimental SAXS curve: %s' % self.experimentalSAXS)

        if not self.otherCrysol.empty():
            summary.append('Other crysol parameters: %d' % self.otherCrysol)

        # Goodness of fit
        fnInt = self._getPath("pseudoatoms00.fit")

        if exists(fnInt):
            x = numpy.loadtxt(fnInt,skiprows=1)
            diff = numpy.log(x[:,1])-numpy.log(x[:,2])
            idx = numpy.isfinite(diff)
            RMS = math.sqrt(1.0/numpy.sum(idx)*numpy.dot(diff[idx],diff[idx]))
            summary.append("RMS=%f" % RMS)
        return summary

    def _methods(self):
        summary = []
        summary.append('We computed the SAXS curve of the model %s.'
                       % self.inputStructure.get().getNameId())
        summary.append('We simulated the SAXS curve using the method described '
                       'in [Svergun1995] up to a frequency of %f (1/Angstroms),'
                       ' with %d harmonics and %d samples.'
                       % (self.maximumFrequency, self.numberOfHarmonics,
                          self.numberOfSamples))
        return summary

    def _citations(self):
        return ['Svergun1995']
    
    def _validate(self):
        errors = []
        if which('crysol') is '':
            errors.append('You should have the program crysol in the PATH')
        return errors
        