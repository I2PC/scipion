#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from os.path import join, exists
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
from xmipp import MetaData, Image, MDL_IMAGE, MDL_ITER, MDL_LL, AGGR_SUM, MDL_REF3D, MDL_WEIGHT, \
getBlocksInMetaDataFile, MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDValueEQ, MDL_SAMPLINGRATE, MDL_CTF_MODEL
from protlib_utils import runShowJ, getListFromVector, getListFromRangeString
from protlib_parser import ProtocolParser
from protlib_xmipp import redStr, cyanStr
from protlib_gui_ext import showWarning
from protlib_filesystem import xmippExists, findAcquisitionInfo, moveFile, replaceBasenameExt
from protocol_ml2d import lastIteration


class ProtRelion3D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.relion3d.name, scriptname, project)
        self.Import = 'from protocol_relion3d import *'
        self.ParamsStr = ''
        if self.NumberOfMpi > 1:
            self.program = 'relion_refine_mpi'
        else:
            self.program = 'relion_refine'
        self.ParamsDict['program'] = self.program
        
        self.addParam('ORoot', self.WorkingDir + '/')        
        acquisionInfo = self.findAcquisitionInfo(self.ImgMd)

        if not acquisionInfo is None: 
            md = MetaData(acquisionInfo)
            self.addParam('SamplingRate', md.getValue(MDL_SAMPLINGRATE, md.firstObject()))
        
    def summary(self):
        md = MetaData(self.ImgMd)
        lines = ["Input images:  [%s] (<%u>)" % (self.ImgMd, md.size())]
       
        return lines
    
    def validate(self):
        errors = []
        # TODO: Check if relion is installed
        md = MetaData(self.ImgMd)
        if md.containsLabel(MDL_IMAGE):
            # Check that have same size as images:
            from protlib_xmipp import validateInputSize
            mdRef = MetaData(self.RefMd)
            references = mdRef.getColumnValues(MDL_IMAGE)
            for ref in references:
                if not xmippExists(ref):
                    errors.append("Reference: <%s> doesn't exists" % ref)
            if len(errors) == 0:
                validateInputSize(references, self.ImgMd, errors)
        else:
            errors.append("Input metadata <%s> doesn't contains image label" % self.ImgMd)
            
        return errors 
    
    def defineSteps(self):        
        restart = False
        if restart:            #Not yet implemented
            pass
#        # Restarting a previous run...
#        else:
#            # Execute protocol in the working directory
#            os.chdir(self.WorkingDir)
#            self.restart_MLrefine3D(RestartIter)
        else:
            self.ImgStar = self.workingDirPath(replaceBasenameExt(self.ImgMd, '.star'))
            self.insertStep('convertImagesMd', [self.ImgStar],
                            inputMd=self.ImgMd, outputRelion=self.ImgStar                            
                            )
            
    def insertRelionRefine(self):
        pass

        
def convertImagesMd(inputMd, outputRelion):
    """ Convert the Xmipp style MetaData to one ready for Relion.
    Main differences are: STAR labels are named different and
    in Xmipp the images rows contains a path to the CTFModel, 
    while Relion expect the explict values in the row.
    Params:
     input: input filename with the Xmipp images metadata
     output: output filename for the Relion metadata
    """
    md = MetaData(inputMd)
    # Get the values (defocus, magnification, etc) from each 
    # ctfparam files and put values in the row
    md.fillExpand(MDL_CTF_MODEL)
    # Create the mapping between relion labels and xmipp labels
    exportMdToRelion(md, outputRelion)
    
             
def renameOutput(log, WorkingDir, ProgId):
    ''' Remove ml2d prefix from:
        ml2dclasses.stk, ml2dclasses.xmd and ml2dimages.xmd'''
    prefix = '%s2d' % ProgId
    for f in ['%sclasses.stk', '%sclasses.xmd', '%simages.xmd']:
        f = join(WorkingDir, f % prefix)
        nf = f.replace(prefix, '')
        moveFile(log, f, nf)
