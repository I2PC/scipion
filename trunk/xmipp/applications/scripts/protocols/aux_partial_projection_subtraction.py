
from xmipp import *
import launch_job
import os

def joinImageCTF(_log, dict):
    tmpMD = MetaData(dict['CTFgroupName'])
    MDaux = MetaData(dict['filename_currentAngles'])
    outMD = MetaData()
    outMD.join(MDaux, tmpMD, MDL_IMAGE, NATURAL_JOIN)
    outMD.write(dict['DocFileRef'], MD_APPEND)

#reconstruct
def reconstructVolume(_log, dict):
    #xmipp_reconstruct_fourier -i ctfgroup_1@ctfgroups_Iter_13_current_angles.doc -o rec_ctfg01.vol --sym i3 --weight
    parameters  = ' -i ' +  dict['DocFileRef'] 
    parameters += ' -o ' +  dict['reconstructedVolume'] 
    parameters += ' --sym ' + dict['SymmetryGroup'] +'h'
    parameters += ' --weight'
    doParallel = dict['DoParallel']
    if (doParallel):
            parameters += ' --mpi_job_size ' + dict['MpiJobSize']
            parameters += ' --thr ' + str(dict['NumberOfThreads'])

    launch_job.launch_job('xmipp_reconstruct_fourier',
                             parameters,
                             _log,
                             doParallel,
                             dict['NumberOfMpiProcesses'],
                             dict['NumberOfThreads'],
                             dict['SystemFlavour'])


def maskVolume(_log, dict):
    parameters  = ' -i ' +  dict['reconstructedVolume'] 
    parameters += ' -o ' +  dict['maskReconstructedVolume'] 
    
    parameters += ' --mask raised_crown -%d -%d 2' %(dict['dRradiusMin'], dict['dRradiusMax'])
    launch_job.launch_job("xmipp_transform_mask",
                             parameters,
                             _log,
                             False,1,1,'')
    


    #project
def createProjections(_log, dict):
    #dirname = 'ReferenceLibrary'

    doParallel = dict['DoParallel']

    parameters  = ' -i ' +  dict['maskReconstructedVolume'] 
    parameters += ' --experimental_images ' +  dict['DocFileRef']
    
    parameters += ' -o ' +  dict['referenceStack']
    parameters += ' --sampling_rate ' + dict['AngSamplingRateDeg']
    parameters += ' --sym ' + dict['SymmetryGroup'] +'h'
    parameters += ' --compute_neighbors --near_exp_data ' 
    parameters += ' --angular_distance ' + str(dict['MaxChangeInAngles'])
    if (doParallel):
            parameters += ' --mpi_job_size ' + dict['MpiJobSize']

    launch_job.launch_job('xmipp_angular_project_library',
                             parameters,
                             _log,
                             doParallel,
                             dict['NumberOfMpiProcesses'] *dict['NumberOfThreads'],
                             1,
                             dict['SystemFlavour'])
                


    #apply transform
#
#    #resta
#def applyGeoScript(_log, dict):
##def readDocfileAndPairExperimentalAndReferenceImages(self, _szDocFile):
#
#    import launch_job
#    import os
#    import xmipp
#    #import docfiles 
#
#    #create sel file from docfile
#    _message = 'creating auxiliary sel file from doc file'
#    self.mylog.info(_message)
#    print '*************************************************'
#    print '* ', _message 
#    print '*************************************************'
#    dirname  = 'Subtraction/ShiftedImages'
#    substractedDir = 'Subtraction/Substracted'
#    if not os.path.isdir("./" + dirname + "/"):
#        os.mkdir("./" + dirname + "/")   
#    if not os.path.isdir("./" + substractedDir + "/"):
#        os.mkdir("./" + substractedDir + "/")  
#
#    #doc=docfiles.docfile(_szDocFile) 
#    doc = xmipp.MetaData(_szDocFile)
#    doc.write('checkfile.doc')
#    doc2 = xmipp.MetaData(_szDocFile)
#    doc3 = xmipp.MetaData(_szDocFile)
#
#    #size = len(doc.lineLst[0])-1 
#    selFileName  = dirname  +'/' + self.ProjOutRootName + '_v3.sel'
#    selFileName2 = dirname  +'/' + self.ProjOutRootName + '_sh_v3.sel'
#    selFileName3 = substractedDir +'/' + self.ProjOutRootName + '_v3.sel'
#    
#    #Ahora los tres selfiles tienen lo mismo, pero deberian enlazar cada uno
#    # a una carpeta (exp, shifted y subtracted)
#    doc.write(selFileName)
#    doc2.write(selFileName2)
#    doc3.write(selFileName3)
#    
#    #Bucle antiguo equivalente
#    #for i in range(len(doc.lineLst)):
#    #    newline =  doc.lineLst[i][size] + ' 1\n'
#    #    fh.write(newline)
#    #    myFileName  = os.path.basename(doc.lineLst[i][size])
#    #    newline     =  dirname +'/' + myFileName + ' 1\n'
#    #    fh2.write(newline)
#    #    newline     =  substractedDir +'/' + myFileName + ' 1\n'
#    #    fh3.write(newline)
#
#    #call to xmipp_header_apply, now xmipp_transform_geometry
#    #create image with same name but different directory 
#    xmpi_run_file='./applyGeo_v3.sh'
#    fh = open(xmpi_run_file,'w')
#
#        # Fragmento de ProjMatch para utilizar en el siguiente bucle
#        #for id in mD:
#        #    fsc = mD.getValue(MDL_RESOLUTION_FRC, id)
#        #    if(fsc < 0.5):
#        #       freq=mD.getValue(MDL_RESOLUTION_FREQ, id)
#        #       break
#
#    for id in doc:
#        inFile  = doc.getValue(MDL_IMAGE, id)
#        print inFile 
#        outFile = dirname + '/' + os.path.basename(inFile)
#        newline  = 'xmipp_transform_geometry '
#        newline += ' -i ' + inFile 
#        newline += ' -o ' + outFile 
#        newline += ' --dont_wrap ' + '\n' 
#        fh.write(newline)
#    
#    #bucle anterior para xmipp_transform_geometry
#    #for i in range(len(doc.lineLst)):
#        #    inFile  = doc.lineLst[i][size]
#        #    outFile = dirname + '/' +os.path.basename(doc.lineLst[i][size])
#        #    newline  = 'xmipp_transform_geometry '
#        #    newline += ' -i ' + inFile 
#        #    newline += ' -o ' + outFile 
#        #    newline += ' --dont_wrap ' + '\n' 
#        #    fh.write(newline)
#
#    fh.close()
#    os.chmod(xmpi_run_file,0755)
#
#    ### ! Todavia no lanzamos el script
#    #if(DoParallel):
#        #     command =' -i '+ xmpi_run_file
#        #     launch_job.launch_job("xmipp_run",
#        #                            command,
#        #                            self.mylog,
#        #                            True,
#        #                            NumberOfMpiProcesses*NumberOfThreads,
#        #                            1,
#        #                            SystemFlavour)
#    #else:
#        #     _mylog.info(xmpi_run_file )
#        #     os.system(xmpi_run_file )
#
#    
def subtractionScript(_log, dict):
    md = MetaData(dict['DocFileRef'])#experimental images
    referenceStackName = dict['referenceStack']#reference projection for a given defocus group
    imgExp =  Image()
    imgRef =  Image()

    #apply ctf to a temporary reference    
    for id in md:
        inFile  = md.getValue(MDL_IMAGE, id)
        refNum  = md.getValue(MDL_REF, id)
        refImg          =  '%06d@%s'%(refNum,referenceStackName)
        imgExp.read(inFile)
        imgRef.read(refImg)
        
        outImgName      = substractedDir + '/' + os.path.basename(inFile)
        #shiftExpImg = dirname + '/' + os.path.basename(inFile)
        a1-a2
    
