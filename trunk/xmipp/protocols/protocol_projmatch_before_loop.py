from distutils.errors   import DistutilsInternalError
from types              import StringTypes
import os
from xmipp import *
from protlib_utils import printLog
from protlib_filesystem import uniqueFilename
from os.path import join

#------------------------------------------------------------------------
#make ctf groups
#------------------------------------------------------------------------
def executeCtfGroups (_log, 
                                 CTFDatName, 
                                 CtfGroupDirectory,
                                 CtfGroupMaxDiff, 
                                 CtfGroupMaxResol,
                                 CtfGroupRootName,
                                 DataArePhaseFlipped,
                                 DoAutoCtfGroup,
                                 DoCtfCorrection,
                                 PaddingFactor, 
                                 SamplingRate ,
                                 SelFileName,
                                 SplitDefocusDocFile,
                                 WienerConstant, 
                                 ) :
    import glob, sys,shutil
    from protlib_utils import runJob
    if not os.path.exists(CtfGroupDirectory):
        os.makedirs(CtfGroupDirectory)
    printLog("executeCtfGroups01", _log)

    if(not DoCtfCorrection):
        MD=MetaData(SelFileName)
        block_name= 'ctfGroup'+str(1).zfill(FILENAMENUMBERLENGTH) +\
                    '@' + join(CtfGroupDirectory, CtfGroupRootName) +'_images.sel'
        MD.write(block_name)
        return 1

#    print '*********************************************************************'
#    print '* Make CTF groups'
#    remove all entries not present in sel file by
#    join between selfile and metadatafile
    MDctfdata = MetaData();
    MDctfdata.read(CTFDatName)

    MDsel = MetaData();
    MDsel.read(SelFileName)
    MDctfdata.intersection(MDsel,MDL_IMAGE)
    tmpCtfdat = uniqueFilename(CTFDatName)
    MDctfdata.write(tmpCtfdat)
    command = \
              ' --ctfdat ' + tmpCtfdat + \
              ' -o ' + CtfGroupDirectory + '/' + CtfGroupRootName + ':stk'\
              ' --wiener --wc ' + str(WienerConstant) + \
              ' --pad ' + str(PaddingFactor) + \
              ' --samplingrate ' + str (SamplingRate)

    if (DataArePhaseFlipped):
        command += ' --phase_flipped '

    if (DoAutoCtfGroup):
        command += ' --error ' + str(CtfGroupMaxDiff) + \
                   ' --resol ' + str(CtfGroupMaxResol)
    else:
        if (len(SplitDefocusDocFile) > 0):
            command += ' --split ' + SplitDefocusDocFile
        else:
            message = "Error: for non-automated ctf grouping, please provide a docfile!"
            print '* ', message
            _log.info(message)
            sys.exit()
    
    if runJob(_log,"xmipp_ctf_group",command):
        return 1
#    fn = CtfGroupDirectory + '/'+\
#                  CtfGroupRootName+\
#                 'Info.xmd'
#    MD = MetaData(fn)

def initAngularReferenceFile(_log, BlockWithAllExpImages, CtfGroupDirectory, CtfGroupRootName, DocFileName, DocFileWithOriginalAngles, SelFileName):
    '''Create Initial angular file. Either fill it with zeros or copy input'''
    printLog("initAngularReferenceFile", _log)
    import shutil
    
    MD=MetaData()

    if len(DocFileName)>1:
        #shutil.copy(DocFileName,DocFileWithOriginalAngles)
        MD.read(DocFileName)
        
        
    else:
        MD.read(SelFileName)
        MD.addLabel(MDL_ANGLEROT)
        MD.addLabel(MDL_ANGLETILT)
        MD.addLabel(MDL_ANGLEPSI)
    
    MD.write(BlockWithAllExpImages + '@'+ DocFileWithOriginalAngles)
    block_name= 'ctfGroup'+str(1).zfill(FILENAMENUMBERLENGTH) +\
                    '@' + join(CtfGroupDirectory, CtfGroupRootName) +'_images.sel'
    blocklist = getBlocksInMetaDataFile(join(CtfGroupDirectory, CtfGroupRootName) +'_images.sel')
    
    MDctf = MetaData()
    mDaux = MetaData()
    mdList = [MDL_IMAGE]
    for block in blocklist:
        #read blocks
        MDctf.read(block+'@'+join(CtfGroupDirectory, CtfGroupRootName)+'_images.sel', mdList)
        #add angles to blocks
        mDaux.join(MD,MDctf, MDL_UNDEFINED,MDL_UNDEFINED,NATURAL)
        #MDctf.intersection(MD, MDL_IMAGE)
        block_name= block+'@' +DocFileWithOriginalAngles
        MDctf.write(block_name, MD_APPEND)
        

