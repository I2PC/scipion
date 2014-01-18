from distutils.errors   import DistutilsInternalError
from types              import StringTypes
import os
from xmipp import *
from protlib_utils import printLog
from protlib_filesystem import uniqueFilename, createLink
from os.path import join, basename
#from IPython.parallel.client.client import Metadata

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
    import glob, sys, shutil
    from protlib_utils import runJob
    if not os.path.exists(CtfGroupDirectory):
        os.makedirs(CtfGroupDirectory)
    printLog("executeCtfGroups01", _log)

    if(not DoCtfCorrection):
        MD = MetaData(SelFileName)
        block_name = 'ctfGroup' + str(1).zfill(FILENAMENUMBERLENGTH) + \
                    '@' + join(CtfGroupDirectory, CtfGroupRootName) + '_images.sel'
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
    MDctfdata.intersection(MDsel, MDL_IMAGE)
    tmpCtfdat = uniqueFilename(CTFDatName)
    MDctfdata.write(tmpCtfdat)
    command = \
              ' --ctfdat ' + tmpCtfdat + \
              ' -o ' + CtfGroupDirectory + '/' + CtfGroupRootName + ':stk'\
              ' --wiener --wc ' + str(WienerConstant) + \
              ' --pad ' + str(PaddingFactor) + \
              ' --sampling_rate ' + str (SamplingRate)

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
    
    if runJob(_log, "xmipp_ctf_group", command):
        return 1
#    fn = CtfGroupDirectory + '/'+\
#                  CtfGroupRootName+\
#                 'Info.xmd'
#    MD = MetaData(fn)

def initAngularReferenceFile(_log, BlockWithAllExpImages, CtfGroupDirectory, CtfGroupRootName, DocFileName, DocFileWithOriginalAngles, SelFileName):
    '''Create Initial angular file. Either fill it with zeros or copy input'''
    printLog("initAngularReferenceFile", _log)
    import shutil
    
    MD = MetaData()

    if len(DocFileName) > 1:
        #shutil.copy(DocFileName,DocFileWithOriginalAngles)
        MD.read(DocFileName)
        
        
    else:
        MD.read(SelFileName)
        MD.addLabel(MDL_ANGLE_ROT)
        MD.addLabel(MDL_ANGLE_TILT)
        MD.addLabel(MDL_ANGLE_PSI)
    
    MD.write(BlockWithAllExpImages + '@' + DocFileWithOriginalAngles)
    block_name = 'ctfGroup' + str(1).zfill(FILENAMENUMBERLENGTH) + \
                    '@' + join(CtfGroupDirectory, CtfGroupRootName) + '_images.sel'
    blocklist = getBlocksInMetaDataFile(join(CtfGroupDirectory, CtfGroupRootName) + '_images.sel')
    
    MDctf = MetaData()
    mDaux = MetaData()
    mdList = [MDL_IMAGE]
    for block in blocklist:
        #read blocks
        MDctf.read(block + '@' + join(CtfGroupDirectory, CtfGroupRootName) + '_images.sel', mdList)
        #add angles to blocks
        mDaux.join(MD, MDctf, MDL_UNDEFINED, MDL_UNDEFINED, NATURAL)
        #MDctf.intersection(MD, MDL_IMAGE)
        block_name = block + '@' + DocFileWithOriginalAngles
        MDctf.write(block_name, MD_APPEND)

def createResults(log
                 , CTFDatName
                 , DoCtfCorrection
                 , inDocfile
                 , listWithResultVolume
                 , resultsImages
                 , resultsClasses3DRef
                 , resultsClasses3DRefDefGroup
                 , resultsVolumes
                 , selFileName
                 , workingDir
                  ):
    ''' Create standard output results_images, result_classes'''
    from os import remove
    from os.path import exists
    #create metadata file with volume names
    mdVolume = MetaData()
    for resultVolume in listWithResultVolume:
        objId = mdVolume.addObject()
        mdVolume.setValue(MDL_IMAGE,resultVolume,objId)
        #link also last iteration volumes
        createLink(log,resultVolume,join(workingDir, basename(resultVolume)))
    mdVolume.write(resultsVolumes)
    #read file with results
    allExpImagesinDocfile = FileName()
    all_exp_images="all_exp_images"
    allExpImagesinDocfile.compose(all_exp_images, inDocfile)
    md = MetaData(allExpImagesinDocfile)
    #read file with ctfs
    mdOut = MetaData()
    if DoCtfCorrection:
        #read only image and ctf_model
        mdOut = MetaData()
        mdCTF = MetaData()
        mdCTF.read(CTFDatName, [MDL_IMAGE, MDL_CTF_MODEL])
        md.addIndex(MDL_IMAGE)
        mdCTF.addIndex(MDL_IMAGE)
        mdOut.join (md, mdCTF, MDL_IMAGE, MDL_IMAGE, NATURAL_JOIN)
    else:
        mdOut = MetaData(md)#becareful with copy metadata since it only copies pointers
    mdref3D = MetaData()
    mdrefCTFgroup = MetaData()
    mdref3D.aggregate(mdOut, AGGR_COUNT, MDL_REF3D, MDL_REF3D, MDL_COUNT)
    mdrefCTFgroup.aggregate(mdOut, AGGR_COUNT, MDL_DEFGROUP, MDL_DEFGROUP, MDL_COUNT)
    #total number Image
    numberImages = mdOut.size()
    #total number CTF
    numberDefocusGroups = mdrefCTFgroup.size()
    #total number volumes 
    numberRef3D = mdref3D.size()
    fnResultImages=FileName()
    fnResultClasses=FileName()

    fnResultImages.compose("images",resultsImages)
    comment  = " numberImages=%d..................................................... "%numberImages
    comment += " numberDefocusGroups=%d................................................."%numberDefocusGroups
    comment += " numberRef3D=%d........................................................."%numberRef3D
    #
    mdAux1 = MetaData()
    mdAux2 = MetaData()
    mdAux1.read(selFileName,[MDL_IMAGE,MDL_ITEM_ID])
    mdAux2.join(mdOut, mdAux1, MDL_IMAGE, MDL_IMAGE, LEFT_JOIN)
    mdAux2.setComment(comment)
    mdAux2.write(fnResultImages)
    #
    #make query and store ref3d with all
    for id in mdref3D:
        ref3D = mdref3D.getValue(MDL_REF3D, id)
        md.importObjects(mdOut, MDValueEQ(MDL_REF3D, ref3D))
        fnResultClasses.compose(("images_ref3d%06d"%ref3D),resultsClasses3DRef)
        md.write(fnResultClasses,MD_APPEND)
    
    md2=MetaData()
    for id in mdref3D:
        ref3D = mdref3D.getValue(MDL_REF3D, id)
        #a multiquey would be better but I do not know how to implement it in python
        md.importObjects(mdOut, MDValueEQ(MDL_REF3D, ref3D))
        for id in mdrefCTFgroup:
            defocusGroup = mdrefCTFgroup.getValue(MDL_DEFGROUP, id)
            md2.importObjects(md, MDValueEQ(MDL_DEFGROUP, defocusGroup))
            fnResultClasses.compose(("images_ref3d%06d_defocusGroup%06d"%(ref3D,defocusGroup)),resultsClasses3DRefDefGroup)
            md2.write(fnResultClasses,MD_APPEND)
