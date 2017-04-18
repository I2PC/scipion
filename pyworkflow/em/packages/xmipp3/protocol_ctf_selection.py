# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Amaya Jimenez    (ajimenez@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import collections
from itertools import izip

from pyworkflow.utils.path import cleanPath, removeBaseExt
from pyworkflow.object import Set, Float, String, Object
import pyworkflow.protocol.params as params
import pyworkflow.em as em
from pyworkflow.em.metadata import Row, MetaData

import convert 
import xmipp



class XmippProtCTFSelection(em.ProtCTFMicrographs):
    """
    Protocol to make a selection of meaningful CTFs in basis of the defocus values,
    the astigmatism, and the resolution.
    """
    _label = 'ctf selection'
    
    def __init__(self, **args):
        em.ProtCTFMicrographs.__init__(self, **args)
        self._freqResol = {}

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        # Read a ctf estimation
        form.addParam('inputCTF', params.PointerParam, pointerClass='SetOfCTF',
                      important=True,
                      label="Estimated CTF",
                      help='Estimated CTF to evaluate.')
        form.addParam('maxDefocus', params.FloatParam, default=1,
                      label='Maximum defocus (A)',
                      help='Maximum value for defocus in Angstroms. '
                           'If the evaluated CTF does not fulfill this requirement, it will be discarded.')
        form.addParam('minDefocus', params.FloatParam, default=1,
                      label='Mimimum defocus (A)',
                      help='Minimum value for defocus in Angstroms. '
                           'If the evaluated CTF does not fulfill this requirement, it will be discarded.')
        form.addParam('astigmatism', params.FloatParam, default=1,
                      label='Astigmatism (A)',
                      help='Maximum value allowed for astigmatism in Angstroms. '
                           'If the evaluated CTF does not fulfill this requirement, it will be discarded.')
        form.addParam('resolution', params.FloatParam, default=1,
                  label='Resolution (A)',
                  help='Minimum value for resolution in Angstroms. '
                       'If the evaluated CTF does not fulfill this requirement, it will be discarded.')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self._insertFunctionStep('readCTF')


    # --------------------------- STEPS functions -------------------------------
    def readCTF(self):
        inputCTF = self.inputCTF.get()

        for ctf in inputCTF:
            ctfName = ctf.getMicrograph().getMicName()
            print "ctf: ", ctfName
            defocusU = ctf.getDefocusU()
            print "defocusU = ", defocusU
            defocusV = ctf.getDefocusV()
            print "defocusV = ", defocusV
            astigm = defocusU - defocusV
            #astigm = ctf.getDefocusRatio()
            print "astigm = ", astigm
            resol = ctf._ctffind4_ctfResolution.get()
            print "resol = ", resol


    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        message = []
        return message
    
    def _methods(self):
        pass#nothing here
    
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        message = [ ]
        fileCTF = self.inputCTF.get()
        if fileCTF == None:
            message.append("You must specify a set of CTFs.")
        return message

    def _stepsCheck(self):
        # Just to avoid the stream checking inherited from ProtCTFMicrographs
        pass
    
        
