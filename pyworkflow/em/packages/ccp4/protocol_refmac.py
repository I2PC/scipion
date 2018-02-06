# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
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

import os
import stat

import pyworkflow.protocol.constants as const
from pyworkflow import VERSION_1_2
from pyworkflow.em import Volume, PdbFile
from pyworkflow.em.packages.ccp4.refmac_template_ifft import template_ifft
from pyworkflow.em.packages.ccp4.refmac_template_mask import \
    template_mask
from pyworkflow.em.packages.ccp4.refmac_template_refine import template_refine
from pyworkflow.em.pdb_handler import fixCRYSrecordToPDBFile
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.convert_header.CCP4.convert import (
    adaptFileToCCP4, runCCP4Program, ORIGIN, getProgram)
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, \
    BooleanParam


class CCP4ProtRunRefmac(EMProtocol):
    """ generates files for volumes and FSCs to submit structures to EMDB
    """
    _label = 'refmac'
    _program = ""
    _version = VERSION_1_2
    refmacMaskScriptFileName = "mask_refmac.sh"
    refmacIfftScriptFileName = "ifft_refmac.sh"
    refmacRefineScriptFileName = "refine_refmac.sh"
    OutPdbFileName = "%s-refined.pdb"
    createMaskLogFileName="mask.log"
    refineLogFileName = "refine.log"
    fftLogFileName = "ifft.log"
    maskedMapFileName = "masked_fs"
    refmacShiftsNames=["pdbin_cell", "pdbin_shifts", "pdbout_cell",
                       "pdbout_shifts"]
    REFMAC='refmac5'

    def __init__(self, **kwargs):
            EMProtocol.__init__(self, **kwargs)

    #--------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, label="Input Volume", important=True,
                      pointerClass='Volume',
                      help='This is the unit cell volume.')#que pasa si la extension no es mrc?
        form.addParam('inputStructure', PointerParam, label="Input PDB file", important=True,
                      pointerClass='PdbFile', help='Specify a PDB object.')
        form.addParam('maxResolution', FloatParam, default=5,
                      label='Max. Resolution (A):', help="Max resolution used in the refinement (Angstroms).")
        form.addParam('minResolution', FloatParam, default=200,
                      label='Min. Resolution (A):', help="Min resolution used in the refinement (Angstroms).")
        form.addParam('nRefCycle', IntParam, default=30, expertLevel=const.LEVEL_ADVANCED,
                      label='Number of refinement iterations:',
                      help='Specify the number of cycles of refinement.\n')
        form.addParam('weightMatrix', FloatParam, default=0.01, expertLevel=const.LEVEL_ADVANCED, label= 'Matrix refinement weight:',
                      help='Weight between the density map and the chemical constraints. Smaller means less weight for the EM map.\n')
        form.addParam('generateMaskedVolume', BooleanParam, default=True, label="Generate masked volume",
                      expertLevel=const.LEVEL_ADVANCED, important=True, help='If set to True, the masked volume will be generated')
        form.addParam('SFCALCmapradius', FloatParam, default=3, expertLevel=const.LEVEL_ADVANCED,
                      label='SFCALC mapradius:', help='Specify how much around molecule should be cut (Angstroms)')
        form.addParam('SFCALCmradius', FloatParam, default=3, expertLevel=const.LEVEL_ADVANCED,
                      label='SFCALC mradius:', help='Specify the radius (Angstroms)to calculate the mask around molecule')
        form.addParam('BFactorSet', FloatParam, default=0, expertLevel=const.LEVEL_ADVANCED,
                      label='B Factor:', help='Specify the B factor value prior to refinement')
        form.addParam('RefiSharpen', FloatParam, default=0, expertLevel=const.LEVEL_ADVANCED,
                      label='Map sharpening:', help='Specify the map sharpening to be used during refinement')


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fixCRYSrecordToPDBFileStep')
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('createDataDictStep')
        self._insertFunctionStep('createMaskScriptFileStep')
        self._insertFunctionStep('executeMaskRefmacStep')
        self._insertFunctionStep('createIfftScriptFileStep')
        self._insertFunctionStep('executeIfftStep')
        self._insertFunctionStep('createIfftOutputStep') # create masked
                                                          # volume
        self._insertFunctionStep('createRefineScriptFileStep')
        self._insertFunctionStep('executeRefineRefmacStep')
        self._insertFunctionStep('createRefmacOutputStep') # create output
                                                           # pdb file
        self._insertFunctionStep('writeFinalResultsTableStep') #Print output
        # results

    # --------------------------- STEPS functions --------------------------------------------
    def fixCRYSrecordToPDBFileStep(self):
        fnVol = self._getInputVolume()
        self.fixedPDBFileName = fixCRYSrecordToPDBFile(
        self.inputStructure.get().getFileName(),
        self._getTmpPath(),
            x=fnVol.getDim()[0],
            y=fnVol.getDim()[1],
            z=fnVol.getDim()[2],
            alpha=90., beta=90., gamma=90.)

    def convertInputStep(self):
        """ convert 3Dmaps to MRC '.mrc' format
        """
        # TODO: IF NO VOLUME NAME USE THE VALUE ASSOCIATED TO THE PDB FILE
        # get input 3D map filename
        fnVol = self._getInputVolume()
        inFileName = fnVol.getFileName()
        # create local copy of 3Dmap (tmp3DMapFile.mrc)
        localInFileName = self._getVolumeFileName()
        origin = fnVol.getOrigin(returnInitIfNone=True).getShifts()
        sampling = fnVol.getSamplingRate()
        adaptFileToCCP4(inFileName, localInFileName, origin, sampling,
                        ORIGIN)

    def createDataDictStep(self):
        self.dict = {}
        self.dict['CCP4_HOME'] = os.environ['CCP4_HOME']
        self.dict['REFMAC_BIN'] = getProgram(self.REFMAC)
        # dict['PDBDIR'] =  os.path.dirname(self.inputStructure.get().getFileName())
        self.dict['PDBDIR'] = os.path.dirname(self.fixedPDBFileName)
        pdfileName = os.path.splitext(
            os.path.basename(self.inputStructure.get().getFileName()))[0]
        self.dict['PDBFILE'] = pdfileName

        self.dict['MAPFILE'] = self._getVolumeFileName()
        #self.inputVolume.get(
        # ).getFileName().replace(     ':mrc', '')

        #       dict['CHIMERA_BIN'] = os.path.join(os.environ['CHIMERA_HOME'],"bin","chimera")
        self.dict['RESOMIN'] = self.minResolution.get()
        self.dict['RESOMAX'] = self.maxResolution.get()
        self.dict['NCYCLE'] = self.nRefCycle.get()
        self.dict['WEIGHT MATRIX'] = self.weightMatrix.get()
        self.dict['OUTPUTDIR'] = self._getExtraPath('')
        self.dict['MASKED_VOLUME'] = self.generateMaskedVolume.get()
        self.dict['SFCALC_mapradius'] = self.SFCALCmapradius.get()
        self.dict['SFCALC_mradius'] = self.SFCALCmradius.get()
        #self.dict['XDIM'] = self.inputVolume.get().getDim()[0]
        #self.dict['YDIM'] = self.inputVolume.get().getDim()[1]
        #self.dict['ZDIM'] = self.inputVolume.get().getDim()[2]
        if self.BFactorSet.get() == 0:
            self.dict['BFACTOR_SET'] = "#BFACtor SET 0"
        else:
            self.dict['BFACTOR_SET'] = "BFACtor SET %f" % self.BFactorSet.get()
        if self.RefiSharpen.get() == 0:
            self.dict['REFI_SHARPEN'] = "#REFI sharpen SET 0"
        else:
            self.dict[
                'REFI_SHARPEN'] = "REFI sharpen %f" % self.RefiSharpen.get()

    def createMaskScriptFileStep(self):
        data_mask = template_mask % self.dict
        f_mask = open(self._getMaskScriptFileName(), "w")
        f_mask.write(data_mask)
        f_mask.close()
        os.chmod(self._getMaskScriptFileName(), stat.S_IEXEC | stat.S_IREAD |
                 stat.S_IWRITE)

    def executeMaskRefmacStep(self):
        # Generic is a env variable that coot uses as base dir for some
        # but not all files. "" force a trailing slash
        runCCP4Program(self._getMaskScriptFileName(),"",
                       {'GENERIC':self._getExtraPath("")})

    def createIfftScriptFileStep(self):
        # sampling
        fnVol = self._getInputVolume()
        sampling = fnVol.getSamplingRate()
        # parse refmac shifts
        refmacShiftDict = self.parseRefmacShiftFile()
        shiftDict = refmacShiftDict[self.refmacShiftsNames[0]]
        self.dict['XDIM']=int(shiftDict[0] / sampling)
        self.dict['YDIM']=int(shiftDict[1] / sampling)
        self.dict['ZDIM']=int(shiftDict[2] / sampling)
        # create template file
        data_ifft = template_ifft % self.dict
        f_ifft = open(self._getIfftScriptFileName(), "w")
        f_ifft.write(data_ifft)
        f_ifft.close()
        os.chmod(self._getIfftScriptFileName(), stat.S_IEXEC | stat.S_IREAD |
                 stat.S_IWRITE)

    def executeIfftStep(self):
        # Generic is a env variable that coot uses as base dir for some
        # but not all files. "" force a trailing slash
        runCCP4Program(self._getIfftScriptFileName(),"",
                       {'GENERIC':self._getExtraPath("")})

    def createIfftOutputStep(self):
        vol = Volume()
        volLocation = vol.setLocation(self._getExtraPath('%s.map'
                                             %self.maskedMapFileName)) #
        # ifft volume
        fnVol = self._getInputVolume()
        sampling = fnVol.getSamplingRate()
        vol.setSamplingRate(sampling)
        #
        """
        ccp4header = Ccp4Header(volLocation, readHeader=True)
        t = Transform()
        x, y, z = ccp4header.getOffset()  # origin output vol coordinates
        fnVol = self._getInputVolume()
        # x, y, z origin input vol coordinates
        x_origin = fnVol.getOrigin().getShifts()[0]
        y_origin = fnVol.getOrigin().getShifts()[1]
        z_origin = fnVol.getOrigin().getShifts()[2]
        # x, y, z origin output vol coordinates
        x += fnVol.getDim()[0] / 2. - x_origin
        y += fnVol.getDim()[1] / 2. - y_origin
        z += fnVol.getDim()[2] / 2. - z_origin
        t.setShifts(-x, -y, -z)  # we follow chimera convention no MRC
        vol.setOrigin(t)
        """
        #
        self._defineOutputs(outputMaskedVolume=vol)
        self._defineSourceRelation(fnVol, self.outputMaskedVolume)
        self._defineSourceRelation(self.inputStructure, self.outputMaskedVolume)

    def createRefineScriptFileStep(self):
        data_refine = template_refine % self.dict
        f_refine = open(self._getRefineScriptFileName(), "w")
        f_refine.write(data_refine)
        f_refine.close()
        os.chmod(self._getRefineScriptFileName(), stat.S_IEXEC | stat.S_IREAD |
                 stat.S_IWRITE)

    def executeRefineRefmacStep(self):
        # Generic is a env variable that coot uses as base dir for some
        # but not all files. "" force a trailing slash
        runCCP4Program(self._getRefineScriptFileName(),"",
                       {'GENERIC':self._getExtraPath("")})

    def createRefmacOutputStep(self):
        pdb = PdbFile()
        pdb.setFileName(self._getOutPdbFileName())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure, self.outputPdb)
        fnVol = self._getInputVolume()
        self._defineSourceRelation(fnVol, self.outputPdb)

    def writeFinalResultsTableStep(self):
        with open(self._getlogFileName()) as input_data:
            for line in input_data:
                if line.strip() == '$TEXT:Result: $$ Final results $$':
                    break
            for line in input_data:
                if line.strip() == '$$':
                    break
                print line

    # --------------------------- UTLIS functions --------------------------------------------
    def _getInputVolume(self):
        if self.inputVolume.get() is None:
            fnVol = self.inputStructure.get().getVolume()
        else:
            fnVol = self.inputVolume.get()
        return fnVol
    def _getOutPdbFileName(self):
        pdfileName = os.path.splitext(os.path.basename(
            self.inputStructure.get().getFileName()))[0]
        return self._getExtraPath(self.OutPdbFileName%pdfileName)

    #def _getVolName(self):
        #return self.Volume.get()
        #pass
    def _getMaskScriptFileName(self):
        return self._getTmpPath(self.refmacMaskScriptFileName)

    def _getIfftScriptFileName(self):
        return self._getTmpPath(self.refmacIfftScriptFileName)

    def _getRefineScriptFileName(self):
        return self._getTmpPath(self.refmacRefineScriptFileName)

    def _getlogFileName(self):
        return self._getExtraPath(self.refineLogFileName)

    def _getVolumeFileName(self, baseFileName="tmp3DMapFile.mrc"):
        return self._getExtraPath(baseFileName)

    def parseRefmacShiftFile(self):
        _shiftFileName =  self._getExtraPath("_shifts.txt")
        f = open(_shiftFileName, 'r')
        refmacShiftDict={}
        for i in range (4):
            l = []
            words = f.readline().split()
            for d in words[2:]:
                l.append(float(d))
            refmacShiftDict[self.refmacShiftsNames[i]] = l
        ##
        s = ""
        for k, v in refmacShiftDict.iteritems():
            s += "%s: %s\n"%(str(k), str(v))
        print s
        ##
        return refmacShiftDict








