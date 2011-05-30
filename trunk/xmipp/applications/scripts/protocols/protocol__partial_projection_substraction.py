#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for subtraction
#
# Example use:
# ./xmipp_partial_projection_subtraction.py
#
# Authors: 
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------

# {file} Protocol Name
ProtocolName='ProjMatch/xmipp_2.4_subtraction_crunchy/xmipp_protocol_projmatch_backup.py'

# {expert}{file} CTFDat file with CTF data:
""" The input selfile may be a subset of the images in the CTFDat file, but all 
    images in the input selfile must be present in the CTFDat file. This field is 
    obligatory if CTF correction is to be performed. 
    Note that this file should be positioned in the project directory, and that the
    image names and ctf parameter filenames should be in absolute paths.
    Usually optained from the protocol file
"""

#Show results for iteration
""" Use data coming from iteration
"""
DisplayIterationNo=2

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
doMask =False

#Crown Mask radius (inner)
dRradiusMin=39

#{expert}  iCrown mask radius center (outter)
dRradiusMax=64


# {expert} Create Reference Library
doRefDirName =True

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
DoParallel=True

# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads=1

# Number of MPI processes to use:
NumberOfMpiProcesses=7

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. Or read: xxx
"""
SystemFlavour='TORQUE-OPENMPI'

#------------------------------------------------------------------------------------------------
# {hidden} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_partial_projection_subtraction.py'


#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#

#from xmipp import *

class partial_projection_subtraction_class:

     
    def maskInputVolume(self, _szInputVolumeName):
	import launch_job

	#MAsked volume filename
	aux = os.path.splitext(_szInputVolumeName)
	szMask_InputVolumeName = aux[0] + '_mask' + aux[1]
	print 'szMask_InputVolumeName: ' + szMask_InputVolumeName

	params  = ' -i ' + _szInputVolumeName
	params += ' -o ' + szMask_InputVolumeName
	params += ' --mask raised_crown -%d -%d 2' %(self.dRradiusMin, self.dRradiusMax)
	print "xmipp_transform_mask",params
	launch_job.launch_job("xmipp_transform_mask",
                             params,
                             self.mylog,
                             False,1,1,'')

	return szMask_InputVolumeName
	

    def CreateProjections(self, _szMask_InputVolumeName, _szDocFile):
	import launch_job
	import os

	#create directory
	#dirname = 'ReferenceLibrary'
	if not os.path.isdir("./" + self.refDirName + "/"):
            os.mkdir("./" + self.refDirName + "/")   

	parameters  = ' -i ' +  _szMask_InputVolumeName
	# Ahora en _szDocFile utilizamos algo como ctfgroup_44@ctfgroups_Iter_13_current_angles_102scl.doc
	parameters += ' --experimental_images ' +  _szDocFile
	parameters += ' -o ' +  "./" + self.refDirName + "/" +  self.ProjOutRootName 
	parameters += ' --sampling_rate ' + str(self.AngSamplingRateDeg)
	parameters += ' --sym ' + self.SymmetryGroup +'h'
	parameters += ' --compute_neighbors --near_exp_data ' 
	parameters += ' --angular_distance ' + str(self.MaxChangeInAngles)
	if (DoParallel):
            parameters += ' --mpi_job_size 1 '

	launch_job.launch_job('xmipp_angular_project_library',
                             parameters,
                             self.mylog,
                             self.DoParallel,
                             self.NumberOfMpiProcesses * self.NumberOfThreads,
                             1,
                             self.SystemFlavour)
			     

    def ReconstructFourier(self, groupFile, reconstructedVolume):
	import launch_job

	#  Ya tenemos los angulos de cada grupo de foco en el fichero.
	# Reconstruir (sin ctf) cada grupo de foco para luego proyectar

	#xmipp_reconstruct_fourier -i ctfgroup_1@ctfgroups_Iter_13_current_angles.doc -o rec_ctfg01.vol --sym i3 --weight

	#parameters  = ' -i ' +  'ctfgroup_' + str(i) + '@ctfgroups_' + filename_currentAngles
	parameters  = ' -i ' +  groupFile
	parameters += ' -o ' +  reconstructedVolume 
	parameters += ' --sym ' + self.SymmetryGroup +'h'
	parameters += ' --weight' 
	if (DoParallel):
            parameters += ' --mpi_job_size 10 '

	launch_job.launch_job('xmipp_reconstruct_fourier',
                             parameters,
                             self.mylog,
                             self.DoParallel,
                             self.NumberOfMpiProcesses * self.NumberOfThreads,
                             1,
                             self.SystemFlavour)



    def readDocfileAndPairExperimentalAndReferenceImages(self, _szDocFile):

	import launch_job
	import os
	import xmipp
	#import docfiles 

	#create sel file from docfile
	_message = 'creating auxiliary sel file from doc file'
	self.mylog.info(_message)
	print '*************************************************'
	print '* ', _message 
	print '*************************************************'
	dirname  = 'Subtraction/ShiftedImages'
	substractedDir = 'Subtraction/Substracted'
	if not os.path.isdir("./" + dirname + "/"):
	    os.mkdir("./" + dirname + "/")   
	if not os.path.isdir("./" + substractedDir + "/"):
	    os.mkdir("./" + substractedDir + "/")  

	#doc=docfiles.docfile(_szDocFile) 
	doc = xmipp.MetaData(_szDocFile)
	
	doc.write('checkfile.doc')
	
        doc2 = xmipp.MetaData(_szDocFile)
        doc3 = xmipp.MetaData(_szDocFile)



	#size = len(doc.lineLst[0])-1 
	selFileName  = dirname  +'/' + self.ProjOutRootName + '_v3.sel'
	selFileName2 = dirname  +'/' + self.ProjOutRootName + '_sh_v3.sel'
	selFileName3 = substractedDir +'/' + self.ProjOutRootName + '_v3.sel'
	
	#doc.write(selFileName)
	#doc.write(selFileName2)
	#doc.write(selFileName3)
	
	#fh  = open(selFileName,'w')
	#fh2 = open(selFileName2,'w')
	#fh3 = open(selFileName3,'w')
	
	# bucle para escribir selfiles
	#for id in doc:
        #    imageName = doc.getValue(MDL_IMAGE, id)
	#    print 'imageName: ',imageName
        
	#Ahora los tres selfiles tienen lo mismo, pero deberian enlazar cada uno
	# a una carpeta (exp, shifted y subtracted)
	doc.write(selFileName)
	doc2.write(selFileName2)
	doc3.write(selFileName3)
	
	#Bucle antiguo equivalente
	#for i in range(len(doc.lineLst)):
	#    newline =  doc.lineLst[i][size] + ' 1\n'
	#    fh.write(newline)
	#    myFileName  = os.path.basename(doc.lineLst[i][size])
	#    newline     =  dirname +'/' + myFileName + ' 1\n'
	#    fh2.write(newline)
	#    newline     =  substractedDir +'/' + myFileName + ' 1\n'
	#    fh3.write(newline)

	#fh.close()
	#fh2.close()
	#fh3.close()

	#may save initial header at the begining so we may restore it later

	#put angles in header
	if(self.doBackupProjectionAngles):
             params =  ' -i ' + selFileName 
             params += ' --extract '
             params += ' -o saveOriginalAngles.doc '
             launch_job.launch_job("xmipp_image_header",
                        	  params,
                        	  self.mylog,
                        	  False,1,1,'')

	params =  ' -i ' +  _szDocFile
	params += ' --assign '
	params += ' --mirror '
	launch_job.launch_job("xmipp_image_header",
                        	  params,
                        	  self.mylog,
                        	  False,1,1,'')

	#call to xmipp_header_apply, now xmipp_transform_geometry
	#create image with same name but different directory 
	xmpi_run_file='./applyGeo_v3.sh'
	fh = open(xmpi_run_file,'w')

        # Fragmento de ProjMatch para utilizar en el siguiente bucle
        #for id in mD:
        #    fsc = mD.getValue(MDL_RESOLUTION_FRC, id)
        #    if(fsc < 0.5):
        #       freq=mD.getValue(MDL_RESOLUTION_FREQ, id)
        #       break

	for id in doc:
	    inFile  = doc.getValue(MDL_IMAGE, id)
	    print inFile 
            outFile = dirname + '/' + os.path.basename(inFile)
            newline  = 'xmipp_transform_geometry '
            newline += ' -i ' + inFile 
            newline += ' -o ' + outFile 
            newline += ' --dont_wrap ' + '\n' 
            fh.write(newline)
	
	#bucle anterior para xmipp_transform_geometry
	#for i in range(len(doc.lineLst)):
        #    inFile  = doc.lineLst[i][size]
        #    outFile = dirname + '/' +os.path.basename(doc.lineLst[i][size])
        #    newline  = 'xmipp_transform_geometry '
        #    newline += ' -i ' + inFile 
        #    newline += ' -o ' + outFile 
        #    newline += ' --dont_wrap ' + '\n' 
        #    fh.write(newline)

	fh.close()
	os.chmod(xmpi_run_file,0755)

	### ! Todavia no lanzamos el script
	#if(DoParallel):
        #     command =' -i '+ xmpi_run_file
        #     launch_job.launch_job("xmipp_run",
        #                            command,
        #                            self.mylog,
        #                            True,
        #                            NumberOfMpiProcesses*NumberOfThreads,
        #                            1,
        #                            SystemFlavour)
	#else:
        #     _mylog.info(xmpi_run_file )
        #     os.system(xmpi_run_file )

	
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
	       refImg      =  '%s%06d.xmp'%(_refDirName+'/'+self.ProjOutRootName,\
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

	       #command += '; rm ' + outTmpRefImg
	       fh.write(command + '\n')
	else:
	# solo entramos en esta opcion
	    for id in doc:
	       inFile  = doc.getValue(MDL_IMAGE, id)
	       refNum  = doc.getValue(MDL_REF, id)
	       refImg      =  '%06d@%s.xmp'%(_refDirName+'/'+self.ProjOutRootName, int(float(refNum))) 
	       outImg      = substractedDir + '/' + os.path.basename(inFile)
	       shiftExpImg = dirname + '/' + os.path.basename(inFile)
	    
	    #for i in range(len(doc.lineLst)):
	       #refImg      =  '%s%06d.xmp'%(_refDirName+'/'+self.ProjOutRootName,\
               #                 	       int(float((doc.lineLst[i])[7]))) 
	       #print 'RefImage: ',refImg
	       #outImg      = substractedDir + '/' + os.path.basename((doc.lineLst[i])[10])
	       #print 'OutImage: ',outImg
	       #shiftExpImg = dirname + '/' + os.path.basename((doc.lineLst[i])[10])
	       #print 'ShiftImage: ',shiftExpImg
	       
	       #outDocFile = 'kk%06d.doc'%(i)
	       command  = ' xmipp_image_operate '	 
	       command +=  ' -i '     + shiftExpImg 
	       command +=  ' --minus ' + refImg
	       command += ' -o ' + outImg 
               fh.write(command + '\n')

	fh.close()
	os.chmod(xmpi_run_file,0755)

	#if(DoParallel):
        #     command =' -i '+ xmpi_run_file 
        #     launch_job.launch_job("xmipp_run",
        #                            command,
        #                            _mylog,
        #                            True,
        #                            NumberOfMpiProcesses*NumberOfThreads,
        #                            1,
        #                            SystemFlavour)
	#else:
        #     _mylog.info(xmpi_run_file )
        #     os.system(xmpi_run_file )



    #init variables
    def __init__(self,
		 ProtocolName,
		 DisplayIterationNo,
		 projOutRootName,
		 szInputVolumeName,
		 szDocFileRef,
		 doMask,
		 dRradiusMin,
		 dRradiusMax,
		 doRefDirName,
		 doBackupProjectionAngles,
		 AngSamplingRateDeg,
		 MaxChangeInAngles,
		 SymmetryGroup,
		 CTFDatName,
		 doCTFCorrection,
		 DoParallel,
		 NumberOfThreads,
		 NumberOfMpiProcesses,
		 SystemFlavour
                 ):

	     self.run_file1='./readDocfileAndPairExperimentalAndReferenceImages_v3.sh'

             print '** Checkpoint 01'
	     
	     ###################################
	     # Main starts here
	     ###################################

	     #create right enviroment
	     import os,sys,shutil,time
	     #scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
	     #sys.path.append(scriptdir) # add default search path

	     #add to pythonPATH
	     scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
	     sys.path.append(scriptdir) # add default search path
	     scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
	     sys.path.append(scriptdir)

	     import log,logging,arg
	     import visualization

	     # New
	     import xmipp
	     import glob
	     #from xmipp import *
	     import launch_job
	     
	     self.DisplayIterationNo = DisplayIterationNo
	     self.projOutRootName = projOutRootName
	     self.szInputVolumeName = szInputVolumeName
	     self.szDocFileRef = szDocFileRef
	     self.doMask = doMask
	     self.dRradiusMin = dRradiusMin
	     self.dRradiusMax = dRradiusMax
	     self.doRefDirName = doRefDirName
	     self.doBackupProjectionAngles = doBackupProjectionAngles
	     self.AngSamplingRateDeg = AngSamplingRateDeg
	     self.MaxChangeInAngles = MaxChangeInAngles
	     self.CTFDatName = CTFDatName
	     self.doCTFCorrection = doCTFCorrection
	     self.DoParallel = DoParallel
	     self.NumberOfThreads = NumberOfThreads
	     self.NumberOfMpiProcesses = NumberOfMpiProcesses
	     self.SystemFlavour = SystemFlavour


	     # Import the projection matching protocol
	     # get WorkingDir and go there
	     pardir=os.path.abspath(os.getcwd())
	     self.ProtocolName=ProtocolName
	     shutil.copy(self.ProtocolName,'partial_protocol.py')
	     import partial_protocol

	     self.WorkingDir=partial_protocol.WorkingDir
	     print 'self.WorkingDir: ' + self.WorkingDir

	     self.LogDir=partial_protocol.LogDir
	     print 'self.LogDir: ' + self.LogDir

	     self.ProjectDir=partial_protocol.ProjectDir
	     print 'self.ProjectDir: ' + self.ProjectDir

	     self.Iteration_Working_Directory='/Iter_'+  str(self.DisplayIterationNo)
	     print 'self.Iteration_Working_Directory: ' + self.Iteration_Working_Directory

	     if ( len(SymmetryGroup) < 1 ):
  	         self.SymmetryGroup = partial_protocol.SymmetryGroup
	     else:
		 self.SymmetryGroup = SymmetryGroup
	     
	     print 'self.SymmetryGroup: ' + self.SymmetryGroup


	     if ( len(szInputVolumeName) < 1 ):
	      self.szInputVolumeName = 'Iter_'+str(self.DisplayIterationNo)+'_reconstruction.vol'
	     else:
	      self.szInputVolumeName = szInputVolumeName

	     print 'self.szInputVolumeName: ' + self.szInputVolumeName

	     self.mylog=log.init_log_system(self.ProjectDir,
                        		self.LogDir,
                        		sys.argv[0],
                        		self.WorkingDir)
	     self.mylog.setLevel(logging.DEBUG)

	     #change to Directory
	     os.chdir(self.WorkingDir+'/'+self.Iteration_Working_Directory)


             #parametros de el protocolo, para project_library   
             #compute sampling rate
	     if(len(AngSamplingRateDeg) < 1):
		self.AngSamplingRateDeg=arg.getComponentFromVector(partial_protocol.AngSamplingRateDeg,\
                                                        	    #self.DisplayIterationNo-1)
	                                             	            self.DisplayIterationNo)
	     else:
		self.AngSamplingRateDeg=int(AngSamplingRateDeg)
		
	     if(len(MaxChangeInAngles) < 1):
		self.MaxChangeInAngles=arg.getComponentFromVector(partial_protocol.MaxChangeInAngles,\
                                                        	    #self.DisplayIterationNo-1)
								    self.DisplayIterationNo)

	     else:
		self.MaxChangeInAngles=int(MaxChangeInAngles)

	     # Si el valor es mayor de 100, lo ponemos a -1
	     if(self.MaxChangeInAngles > 100):
		self.MaxChangeInAngles=-1

	     if(doRefDirName):
		self.ProjOutRootName=projOutRootName
	     else:
		#_ProjOutRootName=partial_protocol.ProjectLibraryBasename
		self.ProjOutRootName='ref'
	     self.refDirName=partial_protocol.LibraryDir

	     ########################

	     #filename_currentAngles = 'Iter_13_current_angles.doc'

	     filename_currentAngles = 'Iter_'+str(self.DisplayIterationNo) + '_current_angles.doc'

	     MDaux = xmipp.MetaData(filename_currentAngles)
	     print '** Checkpoint 02'
	     
	     #MDaux.write('prueba.doc')
	     # OK

	     defGroups=glob.glob('../CtfGroups/ctf_group??????.sel')
	     print '** Checkpoint 03'
	     
	     
	     if not os.path.isdir("./Subtraction/"):
                 os.mkdir("./Subtraction/") 
	     if not os.path.isdir("./Subtraction/vols/"):
                 os.mkdir("./Subtraction/vols/") 

	     #print defGroups

	     tmpMD = xmipp.MetaData()
	     outMD = xmipp.MetaData()
	     #MDaux.addLabel(MDL_CTFMODEL)

	     defGroups.sort()

	     #print defGroups
	     i = 0;

	     for group in defGroups:
	     	     
	         i = i+1 
		 
		 print '*************************************************'
	         print '* Iteration ',i,'/',len(defGroups),'started' 
                 
		 #print group
		 tmpMD.read(group)
		 #print xmipp.MDL_IMAGE
		 #print xmipp.NATURAL_JOIN
		 outMD.join(MDaux, tmpMD, xmipp.MDL_IMAGE, xmipp.NATURAL_JOIN)
		 _szDocFileRef = 'ctfgroup_' + str(i) + '@ctfgroups_' + filename_currentAngles
		 outMD.write(_szDocFileRef, xmipp.MD_APPEND)

		 # inputdocfile    = 'ctfGroup'+str(ictf).zfill(utils_xmipp.FILENAMENUMBERLENTGH) + '@' + CtfGroupName + '_images.sel'

		 reconstructedVolume = 'Subtraction/vols/rec_ctfg' + str(i) + '.vol'

		 # pasamos el SEL del grupo de foco en la variable GROUP
		 self.ReconstructFourier(group, reconstructedVolume)
		 #_mylog,group,reconstructedVolume,_SymmetryGroup)

		 # continuar el protocolo para que se haga el project library con cada sel 
		 # y cada grupo de foco por separado, pero dejando el resultado en el mismo directorio.


		 # Despues de reconstruir un grupo de foco, aplicamos la mascara (o no)
		 #Mask volume
		 if(doMask):
        	     _szMask_InputVolumeName = maskInputVolume(reconstructedVolume)
		 else:
        	     _szMask_InputVolumeName = reconstructedVolume

		 ##################
		 #produce gallery of projections
		 #docfile from the previous iteration
		 #if(len(szDocFileRef)<1):
		 # Ahora utilizamos los angulos de la misma iteracion
		 #   _szDocFileRef= '../Iter_'+str(_iteration_number)\
		 #             + '/Iter_'+str(_iteration_number)\
		 #		 + '_current_angles.doc'
		 #else:
		 #   _szDocFileRef=szDocFileRef

		 
		 # _szMask_InputVolumeName puede tener aplicada una mascara o no, dependiendo de los pasos anteriores
		 if(doRefDirName):
        	     self.CreateProjections(_szMask_InputVolumeName, _szDocFileRef)

		 # Que lo deje todo en otro directorio, por ejemplo Subtracted_v3

		 #######################				  

		 #create pair of experimentar data sinthetic data and substract
		 #if(len(szDocFileRef)<1):
		 #   _szDocFileRef= '../Iter_'+str(_iteration_number)\
		 #                 + '/Iter_'+str(_iteration_number)\
		 #		 + '_current_angles.doc'
		 #else:
		 #   _szDocFileRef=szDocFileRef

		 _CTFDatName=CTFDatName
		 if(len(CTFDatName)<1 and doCTFCorrection):
		     if(len(partial_protocol.CTFDatName)>1):
        		 _CTFDatName= '../' + partial_protocol.CTFDatName
			 
 		 # revisar si realmente hay que pasarle _szDocFileRef a la siguiente funcion	 
		 self.readDocfileAndPairExperimentalAndReferenceImages(_szDocFileRef)
		 
		 # Al finalizar cada iteracion borramos el volumen reconstruido
		 os.rm(reconstructedVolume)
                 
		 print '* Iteration ',i,'/',len(defGroups),'finished' 
                 print '*************************************************'
	         

	     print '*************************************************'
	     print '* Done! ' 
	     print '*************************************************'




# Preconditions
def preconditions(gui):
#    import os
#    retval=True
#    # Check if there is workingdir
#    if WorkingDir == "":
#        message="No working directory given"
#        if gui:
#            import tkMessageBox
#            tkMessageBox.showerror("Error", message)
#        else:
#            print message
#        retval=False
#    
#    # Check that there is a valid list of micrographs
#    if not os.path.exists(PickingDir)>0:
#        message="Cannot find "+PickingDir
#        if gui:
#            import tkMessageBox
#            tkMessageBox.showerror("Error", message)
#        else:
#            print message
#        retval=False
#    
#    # Check that all micrographs exist
#    import xmipp
#    fnPickingParameters=PickingDir+"/protocolParameters.txt"
#    isPairTilt=getParameter("IsPairList",fnPickingParameters)=="True"
#    MicrographSelfile=getParameter("MicrographSelfile",fnPickingParameters)
#    print os.path.curdir
#    mD=xmipp.MetaData();
#    xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, mD)
#    preprocessingDir,dummy=os.path.split(MicrographSelfile)
#    message="Cannot find the following micrographs:\n"
#    NnotFound=0
#    for id in mD:
#        micrograph=mD.getValue(xmipp.MDL_IMAGE)
#        if not os.path.exists(preprocessingDir+"/"+micrograph):
#            message+=preprocessingDir+"/"+micrograph+"\n"
#            NnotFound=NnotFound+1
#        if isPairTilt:
#            micrographTilted=mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1)
#            if not os.path.exists(preprocessingDir+"/"+micrographTilted):
#                message+=preprocessingDir+"/"+micrographTilted+"\n"
#                NnotFound=NnotFound+1
#    if NnotFound>0:
#        if gui:
#            import tkMessageBox
#            tkMessageBox.showerror("Error", message)
#        else:
#            print message
#        retval=False
#            
#    return retval
    import xmipp
    return True


if __name__ == '__main__':
   	# create partial_projection_subtraction_class object
	partial_projection_subtraction = partial_projection_subtraction_class(
		 ProtocolName,
		 DisplayIterationNo,
		 projOutRootName,
		 szInputVolumeName,
		 szDocFileRef,
		 doMask,
		 dRradiusMin,
		 dRradiusMax,
		 doRefDirName,
		 doBackupProjectionAngles,
		 AngSamplingRateDeg,
		 MaxChangeInAngles,
		 SymmetryGroup,
		 CTFDatName,
		 doCTFCorrection,
		 DoParallel,
		 NumberOfThreads,
		 NumberOfMpiProcesses,
		 SystemFlavour)

