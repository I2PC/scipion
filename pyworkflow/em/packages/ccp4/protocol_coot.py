# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *
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

import pyworkflow.utils as pwutils
from pyworkflow import VERSION_1_2
from pyworkflow.em import PdbFile
from pyworkflow.em import Volume
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import EMObject
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.utils.ccp4_utilities.convert import (getProgram,
                                                        copyMRCHeader,
                                                        runCCP4Program,
                                                        cootPdbTemplateFileName,
                                                        cootScriptFileName,
                                                        Ccp4Header)
from pyworkflow.protocol.params import MultiPointerParam, PointerParam, \
    BooleanParam
from pyworkflow.utils.properties import Message
from pyworkflow.em.data import Transform


#TODO: viewer

class CootRefine(EMProtocol):
    """Coot is an interactive graphical application for
macromolecular model building, model completion
and validation. IMPORTANT: press "w" in coot to transfer
the pdb file from coot  to scipion '
"""
    _label = 'coot refinement'
    _program = ""
    _version = VERSION_1_2

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolumes', MultiPointerParam, pointerClass="Volume",
                      label='Input Volume/s', allowsNull=True,
                      help="Set of volumes to process")
        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize', important=True,
                      help='If set to True, particles will be normalized in the'
                           'way COOT prefers it. It is recommended to '
                           '*always normalize your particles* if the maximum value is higher than 1.')
        form.addParam('pdbFileToBeRefined', PointerParam, pointerClass="PdbFile",
                      label='PDB to be refined',
                      help="PDB file to be refined. This PDB object, after refinement, will be saved")
        form.addParam('inputPdbFiles', MultiPointerParam, pointerClass="PdbFile",
                      label='Other referece PDBs',
                      help="Other PDB files used as reference. These PDB "
                           "objects will not be saved")
        form.addSection(label='Help')
        form.addLine('Press "w" in coot to transfer the pdb file from coot  to scipion')
        form.addLine("You may also excute (from script -> python) the command scipion_write(imol)")
        form.addLine("where imol is the PDB id")
        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        #test loop over inputVol
        self.inVolumes = []
        self.norVolumesNames = []
        if self.inputVolumes.get() is None:
            vol = self.pdbFileToBeRefined.get().getVolume()
            inFileName = vol.getFileName()
            self.inVolumes.append(vol)
            self.norVolumesNames.append(self._getVolumeFileName(inFileName))
        else:
            for vol in self.inputVolumes.get():
                inFileName = vol.getFileName()
                self.inVolumes.append(vol)
                self.norVolumesNames.append(self._getVolumeFileName(inFileName))

        convertId = self._insertFunctionStep('convertInputStep', self.inVolumes,
                                             self.norVolumesNames)
        self._insertFunctionStep('runCootStep', self.inVolumes, self.norVolumesNames,
                                 prerequisites=[convertId],
                                 interactive=True)
        #self._insertFunctionStep('createOutputStep', inVolumes,
        #                         norVolumesNames)

    # --------------------------- STEPS functions ---------------------------------------------------

    def convertInputStep(self, inVolumes, norVolumesNames):
        """ convert 3D maps to MRC '.mrc' format
        """

        for inVol, norVolName in zip(inVolumes, norVolumesNames):
            inVolName  = inVol.getFileName()
            if inVolName.endswith("mrc"): inVolName += ":mrc"
            if norVolName.endswith("mrc"): norVolName += ":mrc"
            if not os.path.exists(norVolName):
                if self.doNormalize:
                    img = ImageHandler()._img
                    img.read(inVolName)
                    mean, dev, min, max = img.computeStats()
                    img.inplaceMultiply(1./max)
                    img.write(norVolName)
                else:
                    ImageHandler().convert(inVolName, norVolName)
                copyMRCHeader(inVolName, norVolName, inVol.getOrigin(
                returnInitIfNone=True).getShifts(),
                              inVol.getSamplingRate())

        createScriptFile(0,#imol
                         self._getTmpPath(cootScriptFileName),
                         self._getExtraPath(cootPdbTemplateFileName))

    def runCootStep(self, inVolumes, norVolumesNames):
        #find last created PDB output file
        template = self._getExtraPath(cootPdbTemplateFileName)
        counter=1
        while os.path.isfile(template%counter):
            counter += 1

        #if there is not previous output use pdb file form
        #otherwise use last created pdb file
        if counter == 1:
            pdbFileToBeRefined = self.pdbFileToBeRefined.get().getFileName()
        else:
            pdbFileToBeRefined = template%(counter-1)
            self._log.info("Using last created PDB file named=%s", pdbFileToBeRefined)
        args = ""
        args +=  " --pdb " + pdbFileToBeRefined
        for pdb in self.inputPdbFiles:
            args += " --pdb " + pdb.get().getFileName() # other pdb files
        args +=  " --script " + self._getTmpPath(cootScriptFileName) # script wit auxiliary files
        for volName in norVolumesNames:
            args += " --map " + volName
        #_envDict['COOT_PDB_TEMPLATE_FILE_NAME'] = self._getExtraPath(cootPdbTemplateFileName)
        self._log.info('Launching: ' + getProgram(os.environ['COOT']) + ' ' + args)

        #run in the background
        runCCP4Program(getProgram(os.environ['COOT']), args)

        self.createOutputStep(   inVolumes,
                                 norVolumesNames,
                                 counter)

    def createOutputStep(self, inVolumes,
                                 norVolumesNames,
                                 init_counter=1):
        """ Copy the PDB structure and register the output object.
        """
        template = self._getExtraPath(cootPdbTemplateFileName)
        counter=init_counter
        while os.path.isfile(template%counter):
            #counter -=1
            pdb = PdbFile()
            pdb.setFileName(template%counter)

            outputs = {"outputPdb_%04d"%counter: pdb}
            self._defineOutputs(**outputs)

            #self._defineOutputs(outputPdb=pdb)
            self._defineSourceRelation(self.inputPdbFiles, pdb)
            #self._defineSourceRelation(self.inputVolumes, self.outputPdb)

            for vol in inVolumes:
                self._defineSourceRelation(vol, pdb)
            counter += 1

        if self.doNormalize:
            counter=init_counter
            for inVol, norVolName in zip(inVolumes,norVolumesNames):
                if os.path.exists(norVolName):
                    #break
                    outVol = Volume()
                    sampling = inVol.getSamplingRate()
                    origin = inVol.getOrigin(
                        returnInitIfNone=True).getShifts()
                    outVol.setSamplingRate(sampling)
                    outVol.setOrigin(origin)

                inFileName  = vol.getFileName()
                if inFileName.endswith('.mrc'):
                    inFileName = inFileName + ":mrc"
                outFileName = self._getVolumeFileName(inFileName)
                outVol.setFileName(outFileName)
                outputs = {"output3DMap_%04d"%counter: outVol}
                self._defineOutputs(**outputs)
                self._defineSourceRelation(vol, outVol)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        program = getProgram(os.environ['COOT'])
        if program is None:
            errors.append("Missing variables COOT and/or CCP4_HOME")
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: ~/.config/scipion/scipion.conf")
            errors.append("and set COOT and CCP4_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("CCP4_HOME = %s" % os.environ['CCP4_HOME'])
                errors.append("COOT = %s" % os.environ['COOT'])

        # Check that the input volume exist
        if (not self.pdbFileToBeRefined.get().hasVolume()) \
                and (self.inputVolumes.get() is None):
            errors.append("Error: You should provide a volume.\n")

        return errors

    def _summary(self):
        #Think on how to uprrrrdate this summary with created PDB
        summary = []
        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes(EMObject):
                summary.append("*%s:* \n %s " % (key, output.getObjComment()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")

        return methodsMsgs

    def _citations(self):
        return ['Emsley_2004']

    # --------------------------- UTILS functions --------------------------------------------------

    def _getVolumeFileName(self, inFileName):
        return os.path.join(self._getExtraPath(''),
                            pwutils.replaceBaseExt(inFileName, 'mrc'))

    def replace_at_index(self, tup, ix, val):
       return tup[:ix] + (val,) + tup[ix+1:]

cootScriptHeader='''import ConfigParser
import os
from subprocess import call
mydict={}
mydict['imol']=%d
mydict['aa_main_chain']="B"
mydict['aa_auxiliary_chain']="BB"
mydict['aaNumber']=37
mydict['step']=5
mydict['outfile']='%s'
'''

cootScriptBody='''

def beep(time):
   """I simply do not know how to create a portable beep sound.
      This system call seems to work pretty well if you have sox
      installed"""
   try:
      command = "play --no-show-progress -n synth %f sin 880"%time
      print command
      os.system(command)
   except:
      pass

def _change_chain_id(signStep):
    """move a few aminoacid between chains"""
    global mydict
    dic = dict(mydict)
    if signStep < 0:
        dic['fromAaNumber'] = mydict['aaNumber'] - dic['step'] +1
        dic['toAaNumber']   = mydict['aaNumber']
        dic['fromAaChain']  = mydict['aa_auxiliary_chain']
        dic['toAaChain']    = mydict['aa_main_chain']
    else:
        dic['fromAaNumber'] = mydict['aaNumber']
        dic['toAaNumber']   = mydict['aaNumber'] + dic['step'] -1
        dic['fromAaChain']  = mydict['aa_main_chain']
        dic['toAaChain']    = mydict['aa_auxiliary_chain']
    mydict['aaNumber'] = mydict['aaNumber'] + (dic['step'] * signStep)
    command = "change_chain_id(%(imol)d, '%(fromAaChain)s', '%(toAaChain)s', 1, %(fromAaNumber)d, %(toAaNumber)d)"%dic
    doIt(command)

def _refine_zone(signStep):
    """Execute the refine command"""
    global  mydict
    dic = dict(mydict)
    if signStep <0:
        dic['fromAaNumber'] = mydict['aaNumber'] - dic['step']
        dic['toAaNumber']   = mydict['aaNumber'] + 2
        mydict['aaNumber']  = mydict['aaNumber'] - dic['step']
    else:
        dic['fromAaNumber'] = mydict['aaNumber'] - 2
        dic['toAaNumber']   = mydict['aaNumber'] + dic['step']
        mydict['aaNumber']  = mydict['aaNumber'] + dic['step']
    command = 'refine_zone(%(imol)s, "%(aa_main_chain)s", %(fromAaNumber)d, %(toAaNumber)d, "")'%dic
    doIt(command)

def _updateMol():
    """update global variable using a file as
    [myvars]
    imol: 0
    aa_main_chain: A
    aa_auxiliary_chain: AA
    aaNumber: 82
    called /tmp/config.ini"""
    global mydict
    config = ConfigParser.ConfigParser()
    config.read(os.environ.get('COOT_INI',"/tmp/coot.ini"))
    try:
        mydict['imol']               = int(config.get("myvars", "imol"))
        mydict['aa_main_chain']      = config.get("myvars", "aa_main_chain")
        mydict['aa_auxiliary_chain'] = config.get("myvars", "aa_auxiliary_chain")
        mydict['aaNumber']           = int(config.get("myvars", "aaNumber"))
        mydict['step']               = int(config.get("myvars", "step"))
        mydict['outfile']            = int(config.get("myvars", "outfile"))
    except ConfigParser.NoOptionError:
        pass
    print ("reading:", "/tmp/coot.ini", mydict)
    beep(0.1)

def getOutPutFileName(template):
    """get name based on template that does not exists
    %04d will be incremented untill it does not exists"""
    counter=1
    if "%04d" in template:
        while os.path.isfile(template%counter):
             counter += 1

    return template%counter

def _write():
    """write pdb file, default names
       can be overwritted using coot.ini"""
    #imol = getOutPutFileName(mydict['imol'])
    #outFile = getOutPutFileName(mydict['outfile'])
    dic = dict(mydict)
    dic['outfile']=getOutPutFileName(dic['outfile'])
    command = "write_pdb_file(%(imol)s,'%(outfile)s')"%dic
    doIt(command)
    beep(0.1)

def scipion_write(imol=0):
    """scipion utility for writting files
    args: model number, 0 by default"""
    global mydict
    mydict['imol']=imol
    _write()

def doIt(command):
    """launch command"""
    print "********", command
    eval(command)
    #beep(0.1)

def _printEnv():
    for key in os.environ.keys():
       print "%30s %s \\n" % (key,os.environ[key])

#change chain id
add_key_binding("change_chain_id_down","x", lambda: _change_chain_id(-1))
add_key_binding("change_chain_id_down","X", lambda: _change_chain_id(1))

#refine aminoacid segment
add_key_binding("refine zone m","z", lambda: _refine_zone(1))
add_key_binding("refine zone m","Z", lambda: _refine_zone(-1))

#update global variables
add_key_binding("init global variables","U", lambda: _updateMol())

#write file
add_key_binding("write pdb file","w", lambda: _write())

#print environ
add_key_binding("print enviroment","e", lambda: _printEnv())

'''

def createScriptFile(imol, scriptFile, pdbFile):
    f = open(scriptFile,"w")
    f.write(cootScriptHeader%(imol, pdbFile))
    f.write(cootScriptBody)
    f.close()

