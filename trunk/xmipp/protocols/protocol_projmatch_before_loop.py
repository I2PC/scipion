from distutils.errors   import DistutilsInternalError
from types              import StringTypes
import os
from xmipp import *
from protlib_utils import unique_filename, printLog

def checkVolumeProjSize(_log, ReferenceFileNames, SelFileName):
    """ check references and projection size match"""
    #5a check volumes have same size
    result = True
    message=""
    try:
        (xdim, ydim, zdim, ndim) = SingleImgSize(ReferenceFileNames[0])
        for reference in ReferenceFileNames:
            (xdim2, ydim2, zdim2, ndim2) = SingleImgSize(reference)
    except XmippError, e:
        print "\nEROR:",__name__, 'checkVolumeProjSize:', e
        exit(1)
    if (xdim2, ydim2, zdim2, ndim2) != (xdim, ydim, zdim, ndim):
        message = "Reference %s and %s have not the same size" % \
              (ReferenceFileNames[0], reference) 
        result = False
    
    if result:
        #5b check volume and projections  have same size
        (xdim2, ydim2, zdim2, ndim2) = ImgSize(SelFileName)
        if (xdim2, ydim2) != (xdim, ydim):
            message = "Volume and reference images have not the same size"
            result = False
    #exit(1)
    #print message

    printLog("checkVolumeProjSize", _log)
    if (not result):
        printLog (message,_log,False,True)
        exit(1)


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
                                 SelFileName,
                                 SplitDefocusDocFile,
                                 WienerConstant, 
                                 ) :
    import glob, sys,shutil
    from protlib_utils import runJob
    if not os.path.exists(CtfGroupDirectory):
        os.makedirs(CtfGroupDirectory)

    if(not DoCtfCorrection):
        MD=MetaData(SelFileName)
        block_name= 'ctfGroup'+str(1).zfill(FILENAMENUMBERLENTGH) +\
                    '@' + CtfGroupDirectory+"/"+ CtfGroupRootName +'_images.sel'
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
    tmpCtfdat = unique_filename(CTFDatName)
    MDctfdata.write(tmpCtfdat)
    command = \
              ' --ctfdat ' + tmpCtfdat + \
              ' -o ' + CtfGroupDirectory + '/' + CtfGroupRootName + \
              ' --wiener --wc ' + str(WienerConstant) + \
              ' --pad ' + str(PaddingFactor)

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
    fn = CtfGroupDirectory + '/'+\
                  CtfGroupRootName+\
                 'Info.xmd'
    MD = MetaData(fn)

def checkOptionsCompatibility(_log, DoAlign2D, DoCtfCorrection):
    # Never allow DoAlign2D and DoCtfCorrection together
    printLog("checkOptionsCompatibility", _log)
    if (DoAlign2D and int (DoCtfCorrection)):
        error_message = "You cannot realign classes AND perform CTF-correction. Switch either of them off!"
        printLog(error_message, _log,False,True)
        exit(1)

def initAngularReferenceFile(_log, DocFileName, DocFileWithOriginalAngles, SelFileName):
    '''Create Initial angular file. Either fill it with zeros or copy input'''
    printLog("initAngularReferenceFile", _log)
    import shutil
    if len(DocFileName)>1:
        shutil.copy(DocFileName,DocFileWithOriginalAngles)
    else:
        MD=MetaData(SelFileName)
        MD.addLabel(MDL_ANGLEROT)
        MD.addLabel(MDL_ANGLETILT)
        MD.addLabel(MDL_ANGLEPSI)
        MD.write(DocFileWithOriginalAngles)
