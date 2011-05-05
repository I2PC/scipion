#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for projection matching
#
# Example use:
# ./xmipp_protocol_projmatch.py
#
# Authors: Roberto Marabini, Nov 2009
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------

#Show results for iteration
""" Use data coming from iteration
"""
DisplayIterationNo=12

#{expert} Out files Root name
projOutRootName='out'

#{expert}{file} Select the input volume
"""Name of the reference volume by default Iter_X_reconstruction.vol

"""
szInputVolumeName=''
#{expert}{file} Select docfile used to compute references
"""Name of the doc file used to compute reference library, usually 
   ../Iter_(X-1)/Iter_(X-1)_current_angles.doc
"""
szDocFileRef=''

# {expert} Mask reference volume?
doMask =True

#Crown Mask radius (inner)
dRradiusMin=105

#{expert}  iCrown mask radius center (outter)
dRradiusMax=175


# {expert} Create Reference Library
doRefDirName =False

# {expert} Backup Experimental Proj. angles
doBackupProjectionAngles =True

# {expert}Angular sampling rate
"""Angular distance (in degrees) between neighboring projection  points
   usually obtained from protocol file
"""   
AngSamplingRateDeg=''

# {expert}  Angular search range 
"""Maximum change in rot & tilt  (in +/- degrees)
   usually obtained from protocol file
"""   
MaxChangeInAngles =''

# {expert} Symmetry group
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
    usually obtained from protocol file
"""
SymmetryGroup=''

# {file} Protocol Name
ProtocolName='ProjMatch/empties_sym_408_03/xmipp_protocol_projmatch_backup.py'

# {expert}{file} CTFDat file with CTF data:
""" The input selfile may be a subset of the images in the CTFDat file, but all 
    images in the input selfile must be present in the CTFDat file. This field is 
    obligatory if CTF correction is to be performed. 
    Note that this file should be positioned in the project directory, and that the
    image names and ctf parameter filenames should be in absolute paths.
    Usually optained from the protocol file
"""
CTFDatName=''
# {expert} Correct by CTF:
""" Set to True if you want to correct by CTF
"""
doCTFCorrection=False

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------

# distributed-memory parallelization (MPI)?
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
DoParallel=False

# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads=8

# Number of MPI processes to use:
NumberOfMpiProcesses=5

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. Or read: xxx
"""
SystemFlavour='HOME_MACHINEFILE'


#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
run_file1='readDocfileAndPairExperimentalAndReferenceImages.sh'

def maskInputVolume(_mylog,\
                    _szInputVolumeName,
		    _dRradiusMin,
		    _dRradiusMax):
    import launch_job
    
    #MAsked volume filename
    aux = os.path.splitext(_szInputVolumeName)
    szMask_InputVolumeName = aux[0] + '_mask' + aux[1]
    
    params = ' -i ' +  _szInputVolumeName
    params += ' -o ' +  szMask_InputVolumeName
    params += ' -mask raised_crown -%d -%d 2' %(_dRradiusMin,_dRradiusMax)
    launch_job.launch_job("xmipp_mask",
                         params,
                         _mylog,
                         False,1,1,'')
    
    return szMask_InputVolumeName
    
def CreateProjections(_mylog,\
                _szMask_InputVolumeName,\
		_cSymmetry,\
		_szDocFile,\
		_dSamplingRate,\
		_dAngularDistance,\
                _szProjOutRootName,\
		_refDirName):
    import launch_job
		
    #create directory
    #dirname = 'ReferenceLibrary'
    if not os.path.isdir("./" + _refDirName + "/"):
        os.mkdir("./" + _refDirName + "/")   

    parameters  = ' -i ' +  _szMask_InputVolumeName
    parameters += ' -experimental_images ' +  _szDocFile
    parameters += ' -o ' +  "./" + _refDirName + "/" +  _szProjOutRootName 
    parameters += ' -sampling_rate ' + str(_dSamplingRate)
    parameters += ' -sym ' + _cSymmetry +'h'
    parameters += ' -compute_neighbors -near_exp_data ' 
    parameters += ' -angular_distance ' + str(_dAngularDistance)
    if (DoParallel):
        parameters += ' -mpi_job_size 1 '
    
    launch_job.launch_job('xmipp_angular_project_library',
                         parameters,
                         _mylog,
                         DoParallel,
                         NumberOfMpiProcesses*NumberOfThreads,
                         1,
                         SystemFlavour)
   

def readDocfileAndPairExperimentalAndReferenceImages(_mylog,\
                                                     _szDocFile,\
                                                     _szProjOutRootName,\
						     _doBackupProjectionAngles,
						     _CTFDatName,\
                                                     _refDirName):

    import launch_job
    import docfiles 
    #create sel file from docfile
    _message = 'creating auxiliary sel file from doc file'
    _mylog.info(_message)
    print '*************************************************'
    print '* ', _message 
    print '*************************************************'
    dirname  = 'ShiftedImages'
    substractedDir = 'Substracted'
    if not os.path.isdir("./" + dirname + "/"):
	os.mkdir("./" + dirname + "/")   
    if not os.path.isdir("./" + substractedDir + "/"):
	os.mkdir("./" + substractedDir + "/")   
    doc=docfiles.docfile(_szDocFile) 
    size = len(doc.lineLst[0])-1 
    selFileName  = dirname  +'/' + _szProjOutRootName + '.sel'
    selFileName2 = dirname  +'/' + _szProjOutRootName + '_sh.sel'
    selFileName3 = substractedDir +'/' + _szProjOutRootName + '.sel'
    fh  = open(selFileName,'w')
    fh2 = open(selFileName2,'w')
    fh3 = open(selFileName3,'w')
    for i in range(len(doc.lineLst)):
        newline =  doc.lineLst[i][size] + ' 1\n'
        fh.write(newline)
        myFileName  = os.path.basename(doc.lineLst[i][size])
        newline     =  dirname +'/' + myFileName + ' 1\n'
        fh2.write(newline)
        newline     =  substractedDir +'/' + myFileName + ' 1\n'
        fh3.write(newline)
    fh.close()
    fh2.close()
    fh3.close()

    #may save initial header at the begining so we may restore it later
 
    #put angles in header
    if(_doBackupProjectionAngles):
         params =  ' -i ' + selFileName 
         params += ' -o saveOriginalAngles.doc '
         launch_job.launch_job("xmipp_header_extract",
                              params,
                              _mylog,
                              False,1,1,'')

    params =  ' -i ' +  _szDocFile
    params += ' -mirror '
    launch_job.launch_job("xmipp_header_assign",
                              params,
                              _mylog,
                              False,1,1,'')
    
    #call to xmipp_header_apply
    #create image with same name but different directory 
    xmpi_run_file='applyGeo.sh'
    fh = open(xmpi_run_file,'w')

    for i in range(len(doc.lineLst)):
        inFile  = doc.lineLst[i][size]
        outFile = dirname + '/' +os.path.basename(doc.lineLst[i][size])
        newline  = 'xmipp_header_apply '
        newline += ' -i ' + inFile 
        newline += ' -o ' + outFile + '\n' 
        fh.write(newline)

    fh.close()

    if(DoParallel):
         command =' -i '+ xmpi_run_file
         launch_job.launch_job("xmipp_run",
                                command,
                                _mylog,
                                True,
                                NumberOfMpiProcesses*NumberOfThreads,
                                1,
                                SystemFlavour)
    else:
         _mylog.info(xmpi_run_file )
         os.system(xmpi_run_file )

    #params  = ' -i ' + selFileName 
    #params += ' -oroot ' + dirname + '/' + _szProjOutRootName 
    #launch_job.launch_job("xmipp_header_apply",
    #                     params,
    #                     _mylog,
    #                     False,1,1,'')
    # create shift selfile
    #params  = '"' + dirname + '/' + _szProjOutRootName +'??????'+'.*"'
    #params += ">" + dirname + '/' + _szProjOutRootName +'_sh.sel'
    #launch_job.launch_job("xmipp_selfile_create",
    #                     params,
    #                     _mylog,
    #                     False,1,1,'')

    ##################################################
    # open ctfdat file
    ##################################################
    
    #apply ctf to reference and substract pairs of projections 
    #substract pairs of projections
    xmpi_run_file=run_file1
    fh = open(xmpi_run_file,'w')
    size = len(doc.lineLst[0])-1
    #apply ctf to a temporary file    
    
    import ctfdat
    if (len(_CTFDatName)>1):
        print "_CTFDatName ", _CTFDatName, len(_CTFDatName)
	myctfdat=ctfdat.ctfdat()
	myctfdat.read(_CTFDatName)
	myctfdat.produceCTFDictionary(True)
	#dirname3 = 'Substracted'
	#if not os.path.isdir("./" + dirname3 + "/"):
	#    os.mkdir("./" + dirname3 + "/")
	for i in range(len(doc.lineLst)):
	   #apply ctf to reference
	   shiftExpImg = dirname + '/' + os.path.basename((doc.lineLst[i])[10])
	   refImg      =  '%s%06d.xmp'%(_refDirName+'/'+_szProjOutRootName,\
                                	   int(float((doc.lineLst[i])[7]))) 
	   outTmpRefImg      = substractedDir + '/' + os.path.basename((doc.lineLst[i])[10]) + '.tmp'
	   expImg         = os.path.basename((doc.lineLst[i])[10])
	   command = 'xmipp_fourier_filter '	 
	   command +=  ' -i '     + refImg 
	   command +=  ' -o '     + outTmpRefImg 
	   command +=  ' -fourier_mask ctfpos ' + myctfdat.CTFDictionary[os.path.basename((doc.lineLst[i])[10])] 

	   outImg      = substractedDir + '/' + os.path.basename((doc.lineLst[i])[10])
	   #outDocFile = 'kk%06d.doc'%(i)
	   command += '; xmipp_operate '	 
	   command +=  ' -i '     + shiftExpImg 
	   command +=  ' -minus ' +  outTmpRefImg
	   command += ' -o ' + outImg 

	   command += '; rm ' + outTmpRefImg
	   fh.write(command + '\n')
    else:
	for i in range(len(doc.lineLst)):
	   refImg      =  '%s%06d.xmp'%(_refDirName+'/'+_szProjOutRootName,\
                                	   int(float((doc.lineLst[i])[7]))) 
	   outImg      = substractedDir + '/' + os.path.basename((doc.lineLst[i])[10])
	   shiftExpImg = dirname + '/' + os.path.basename((doc.lineLst[i])[10])
	   #outDocFile = 'kk%06d.doc'%(i)
	   command  = ' xmipp_operate '	 
	   command +=  ' -i '     + shiftExpImg 
	   command +=  ' -minus ' + refImg
	   command += ' -o ' + outImg 
           fh.write(command + '\n')
        	   
    fh.close()

    if(DoParallel):
         command =' -i '+ xmpi_run_file 
         launch_job.launch_job("xmipp_run",
                                command,
                                _mylog,
                                True,
                                NumberOfMpiProcesses*NumberOfThreads,
                                1,
                                SystemFlavour)
    else:
         _mylog.info(xmpi_run_file )
         os.system(xmpi_run_file )

###################################
# Main starts here
###################################

#create right enviroment
import os,sys,shutil
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path
import log,logging,arg
import visualization

# Import the projection matching protocol
# get WorkingDir and go there
pardir=os.path.abspath(os.getcwd())
_ProtocolName=ProtocolName
shutil.copy(_ProtocolName,'partial_protocol.py')
import partial_protocol
_WorkingDir=partial_protocol.WorkingDir
_LogDir=partial_protocol.LogDir
_ProjectDir=partial_protocol.ProjectDir
_iteration_number=DisplayIterationNo
_Iteration_Working_Directory='/Iter_'+  str(_iteration_number)
_SymmetryGroup=partial_protocol.SymmetryGroup
if(len(szInputVolumeName)<1):
 _szInputVolumeName = 'Iter_'+str(_iteration_number)+'_reconstruction.vol'
else:
 _szInputVolumeName = szInputVolumeName

_mylog=log.init_log_system(_ProjectDir,
                           _LogDir,
                           sys.argv[0],
                           _WorkingDir)
_mylog.setLevel(logging.DEBUG)

#change to Directory
os.chdir(_WorkingDir+'/'+_Iteration_Working_Directory)

#Mask volume
if(doMask):
    _szMask_InputVolumeName= maskInputVolume(_mylog,\
				   _szInputVolumeName,\
				  dRradiusMin,\
				  dRradiusMax)
#produce gallery of projections
#docfile from the previous iteration
if(len(szDocFileRef)<1):
   _szDocFileRef= '../Iter_'+str(_iteration_number-1)\
                 + '/Iter_'+str(_iteration_number-1)\
		 + '_current_angles.doc'
else:
   _szDocFileRef=szDocFileRef
#compute sampling rate
if(len(AngSamplingRateDeg)<1):
    _AngSamplingRateDeg=arg.getComponentFromVector(partial_protocol.AngSamplingRateDeg,\
                                                           _iteration_number-1)
else:
    _AngSamplingRateDeg=int(AngSamplingRateDeg)
if(len(MaxChangeInAngles)<1):
    _MaxChangeInAngles=arg.getComponentFromVector(partial_protocol.MaxChangeInAngles,\
                                                           _iteration_number-1)
else:
    _MaxChangeInAngles=int(MaxChangeInAngles)

if(doRefDirName):
    _ProjOutRootName=projOutRootName
else:
    _ProjOutRootName=partial_protocol.ProjectLibraryBasename
_refDirName=partial_protocol.LibraryDir

if(doRefDirName):
    CreateProjections(_mylog,\
	          _szMask_InputVolumeName,\
	          _SymmetryGroup,\
	          _szDocFileRef,\
	          _AngSamplingRateDeg,\
	          _MaxChangeInAngles,\
		  _ProjOutRootName,
		  _refDirName)
#create pair of experimentar data sinthetic data and substract
if(len(szDocFileRef)<1):
   _szDocFileRef= '../Iter_'+str(_iteration_number)\
                 + '/Iter_'+str(_iteration_number)\
		 + '_current_angles.doc'
else:
   _szDocFileRef=szDocFileRef

_CTFDatName=CTFDatName
if(len(CTFDatName)<1 and doCTFCorrection):
    if(len(partial_protocol.CTFDatName)>1):
        _CTFDatName= '../' + partial_protocol.CTFDatName

readDocfileAndPairExperimentalAndReferenceImages(_mylog,\
                                                 _szDocFileRef,\
						 _ProjOutRootName,\
						 doBackupProjectionAngles,\
						 _CTFDatName,\
                                                 _refDirName)
