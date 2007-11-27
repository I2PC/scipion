#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for projection matching
# using self-organizing maps
#
#  - delete and create workdirectory
#
# Example use:
# ./xmipp_protocol_projmatch.py
#
# Author:Roberto Marabini, March 2007
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# {file} Selfile with the input images:
SelFileName='all_images.sel'

# {file} Initial 3D reference map:
ReferenceFileName='reference.vol'
#ReferenceFileName='init_reference/LTA_rot_0.1_norm.vol'

# Working subdirectory: 
WorkDirectory='ProjMatch/Test1'

# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=True

# Number of iterations to perform
NumberofIterations=4

# Resume at iteration
""" This option may be used to finish a previously performed run.
    Set to 1 to start a new run 
    Note: Do NOT delete working directory if this option is not set to 1
"""
ContinueAtIteration=1

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir='/home2/bioinfo/rafa/Protocol/BC/PARA_Roberto'
#ProjectDir='/home/roberto2/Test/PARA_Roberto'

# {expert} Directory name for logfiles:
LogDir='Logs'

#-----------------------------------------------------------------------------
# {section} Mask
#-----------------------------------------------------------------------------
# Mask reference volume?
""" Masking the reference volume will increase the signal to noise ratio.
    Do not provide a very tight mask.
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mask for details
"""
DoMask=True

# Show masked volume
""" Masked volume will be shown. Do not set ths option to true for non iterative processing (jobs sent to queues)
"""
DisplayMask=False

# {file} Binary mask-file used to mask the reference volume
MaskFileName='mask.vol'
#MaskFileName='circular_mask.msk'
#-----------------------------------------------------------------------------
# {section} Projection Matching
#-----------------------------------------------------------------------------
# Perform projection Matching?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Projection_matching and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_projection_matching for details
"""
DoProjectionMatching=True

#  Show projection maching library and classes
""" Show average of projections. Do not set ths option to true for non iterative processing (jobs sent to queues)
"""
DisplayProjectionMatching=True

# Angular sampling rate (in degrees)
"""Angular distance between neighboring projection  points
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "8 8 2 2 " 
    specifies 4 iterations, the first two set the value to 8 
    and the last two to 2. An alternative compact notation 
    is ("2x8 2x0", i.e.,
    2 iterations with value 8, and 2 with value 2).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
AngSamplingRateDeg='3x15 3x10 4x5'

# {expert} Angular Search 
"""Maximum change in rot & tilt  (+/- degrees)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
MaxChangeInAngles='5x1000 5x20 '

# {expert} Maximum change in origin offset (+/- pixels)
"""Maximum change in shift  (+/- pixels)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
MaxChangeOffset='8x1000 2x10'

# Restrict tilt angle search?
DoRetricSearchbyTiltAngle=True

# Lower-value for restricted tilt angle search
Tilt0=40

# Higher-value for restricted tilt angle search
TiltF=140

# {file} If you want to appply symmetry supply the symmetry description file or its acronym:
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetrize
    for a description of the symmetry file and its acronym formats
    dont give anything, if no symmetry is present
"""
Symfile='C6'

# {expert} Additional options for Projection_Matching
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Projection_matching and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_projection_matching for details
    try -Ri xx -Ro yy for restricting angular search (xx and yy are
    the particle inner and outter radius)
    
"""
ProjMatchingExtra=''

#-----------------------------------------------------------------------------
# {section} 2D alignment
#-----------------------------------------------------------------------------
# Perform 2D alignment?
""" After performing a 3D projection matching iteration, each of the
    subsets of images assigned to one of the library projections is
    re-aligned using a 2D-alignment protocol.
    This may serve to remove model bias.
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d for details
    You may specify this option for each iteration. 
    This can be done by a sequence of 0 or 1 numbers (for instance, "1 1 0 0" 
    specifies 4 iterations, the first two applied alig2d while the last 2
    dont. an alternative compact notation is 
    is ("2x1 2x0", i.e.,
    2 iterations with value 1, and 2 with value 0).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    IMPORTANT: if you set this variable to 0 the output  of the projection
    muching step will be copied as output of align2d
"""
DoAlign2D='1 1 0 0'

# Display 2D alignment results
""" Show aligned classes. Do not set ths option to true for non iterative processing (jobs sent to queues)
"""
DisplayAlign2D=True

# Inner radius for rotational correlation:
""" In pixels from the image center
"""
InnerRadius=0

# Outer radius for rotational correlation
""" In pixels from the image center.
"""
OuterRadius=18
#set default OuterRadius=dim/2

#Use Class as reference
""" Use Class image as reference (True) versus use average of assigned
    images to that class (False)
"""
AlignWithClassReference=True
# {expert} Number of align2d iterations:
""" Use at least 3
"""
Align2DIterNr=4

# {expert} Maximum change in origin offset (+/- pixels)
"""Maximum change in shift  (+/- pixels)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
Align2dMaxChangeOffset='2x15 2x10'

# {expert} Maximum change in rotation (+/- degrees)
"""Maximum change in shift  (+/- pixels)
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1000 1000 10 10 " 
    specifies 4 iterations, the first two set the value to 1000 (no restriction)
    and the last two to 10degrees. An alternative compact notation 
    is ("2x1000 2x10", i.e.,
    2 iterations with value 1000, and 2 with value 10).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
Align2dMaxChangeRot='2x1000 2x10'

# {expert} Additional align2d parameters
""" For a complete description, see http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d

  Examples:
  Align2DExtraCommand=\"-trans_only  -filter 10 -sampling 2\"
  Align2DExtraCommand=\"-max_shift 2 -max_rot 3\"

  consider filtering the images with \"-filter 10 -sampling 2\"
  See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d for details
"""
Align2DExtraCommand='-filter 10 -sampling 2'

#-----------------------------------------------------------------------------
# {section} 3D Reconstruction
#-----------------------------------------------------------------------------
# Perform 3D Reconstruction?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Wbp and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_wbp and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art
        for details
"""
DoReconstruction=True

# Display reconstructed volume?
DisplayReconstruction=True

# Reconstructiom method
""" Choose between wbp or art
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "wbp wbp wbp art " 
    specifies 4 iterations, the first three set the value to wbp (no restriction)
    and the last  to art. An alternative compact notation 
    is ("3xwbp 1xart", i.e.,
    3 iterations with wbp, and 1 with art).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
ReconstructionMethod='9xwbp art'
#ReconstructionMethod='3xart wbp'

# {expert} Values of lambda for art
""" IMPORTANT: ou must specify a value of lambda for each iteration even
    if art has not been selected.
    IMPORTANT: NOte that we are using the WLS version of ART that 
    uses geater lambdas than the plain art.
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art
        for details
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, ".1 .1 .3 .3" 
    specifies 4 iterations, the first two set the value to 0.1 
    (no restriction)
    and the last  two to .3. An alternative compact notation 
    is ("2x.1 2x.3").
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
ARTLambda='0.2'

# {expert} Additional reconstruction parameters for ART
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Art
        for details
"""
ARTReconstructionExtraCommand='-k 0.5 -n 10 '

# {expert} Additional reconstruction parameters for WBP
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Wbp and
        http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Mpi_wbp and
        for details
"""
WBPReconstructionExtraCommand=' '

#-----------------------------------------------------------------------------
# {section} Compute Resolution
#-----------------------------------------------------------------------------
# Compute resolution?
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Resolution fo details
"""
DoComputeResolution=True

# Pixel size (in Ang.)
""" This will make that the X-axis in the resolution plots has units 1/Angstrom
"""
ResolSam=3.

# Display resolution?
DisplayResolution=False

#-----------------------------------------------------------------------------
# {section} Low-pass filtering
#-----------------------------------------------------------------------------
# Provide your own filtration frequecy?
"""By default the volume will be filtered at a frecuency equal to
   the  resolution computed with resolution_fsc
   plus a constant provided by the user in the next
   input box. If this option is set to true then the
   filtration will be made at the constant value provided by
   the user WITHOUT adding the resolution computed by resolution_fsc.
"""
SetResolutiontoZero=False


# Constant to by add to the estimate resolution
""" The meaning of this field depends on the previous flag.
    If set as False then  the volume will be filtered at a frecuency equal to
    the  resolution computed with resolution_fsc
    plus the value provided in this field (units pixel^-1)
    If set to False  the volume will be filtered at the resolution
    provided in this field.
    Set this value to any value larger than .5 to avoid filtration.
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, ".15 .15 .1 .1" 
    specifies 4 iterations, the first two set the constant to .15
    and the last two to 0.1. An alternative compact notation 
    is ("2x.15 2x0.1", i.e.,
    4 iterations with value 0.15, and three with value .1).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
ConstantToAddToFiltration='8x.15 2x.1'

#-----------------------------------------------------------------------------
# {section} Cleaning temporal files and Reseting original data
#-----------------------------------------------------------------------------

# Reset image header?
"""Set shifts and angles stored in the image headers to zero 
   in each iteration. You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "1 0 0 0" 
    specifies 4 iterations, the first   resets the header 
    and the last three do not reset it. An alternative compact notation 
    is ("4x1 3x0", i.e.,
    4 iterations with value 0, and three with value 1).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
    Note: Do not set this option to 1  if you want to limit the
          search in projection matching
"""
ResetImageHeader='1 3x0'

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel?
DoParallel=False

# Number of processors to use:
NumberOfCPUs=2

# {file} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
    if your queueing system using an environment variable, give it here (with the leading $, e.g. $PBS_NODEFILE
"""
MachineFile='../../../machines.dat'
#MachineFile='../../bin/machines.dat'

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_projmatch.py'

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#-----------------------------------------------------------------------------
#Do not change these variables
ReferenceVolume='reference_volume.vol'
Proj_Maching_Output_Root_Name="proj_match"
multi_align2d_sel="multi_align2d.sel"
align2d_sel="align2d.sel"
align2d_doc="align2d.doc"
docfile_with_original_angles='docfile_with_original_angles.doc'
#directory to copy images
images_dir='images'
#inpu/output of mask/fourier_filtering
Filtered_Image="filtered_reconstruction"
#inpu/output of mask/fourier_filtering
ReconstrucedVolume="reconstruction"


class projection_matching_class:

   #init variables
   
   def __init__(self,
                _NumberofIterations,
                _ContinueAtIteration,
                _DoMask, 
                _DisplayMask,
                _ReferenceFileName,
                _MaskFileName,
                _DoProjectionMatching,
                _DisplayProjectionMatching,
                _AngSamplingRateDeg,
                _DoRetricSearchbyTiltAngle,
                _Tilt0,
                _TiltF,
                _ProjMatchingExtra,
                _MaxChangeOffset,
                _MaxChangeInAngles,
                _DoAlign2D,
                _InnerRadius,
                _OuterRadius,
                _Align2DIterNr,
                _DisplayAlign2D,
                _Align2DExtraCommand,
                _Align2dMaxChangeOffset,
                _Align2dMaxChangeRot,
                _DisplayReconstruction,
                _DisplayResolution,
                _DoReconstruction,
                _ReconstructionMethod,
                _ARTLambda,
                _ARTReconstructionExtraCommand,
                _WBPReconstructionExtraCommand,
                _DoComputeResolution,
                _ResolSam,
                _SelFileName,
                _WorkDirectory,
                _ProjectDir,
                _LogDir,
                _DoParallel,
                _MyNumberOfCPUs,
                _MyMachineFile,
                _Symfile,
                _ReferenceVolume,
                _Proj_Maching_Output_Root_Name,
                _multi_align2d_sel,
                _AlignWithClassReference,
                _align2d_sel,
                _align2d_doc,
                _ResetImageHeader,
                _SetResolutiontoZero,
                _ConstantToAddToFiltration
                ):
       import os,sys
       scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
       sys.path.append(scriptdir) # add default search path
       import log,logging
       import shutil
       self._WorkDirectory=os.getcwd()+'/'+_WorkDirectory
       #self._SelFileName=_ProjectDir+'/'+str(_SelFileName)
       #self._SelFileName=os.path.abspath(str(_SelFileName))
       #self._SelFileName=os.path.abspath(_SelFileName)
       self._SelFileName=_SelFileName
       selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
       self._ReferenceFileName=os.path.abspath(_ReferenceFileName)
       self._MaskFileName=os.path.abspath(_MaskFileName)
       self._DoMask=_DoMask
       self._DoProjectionMatching=_DoProjectionMatching
       self._DisplayProjectionMatching=_DisplayProjectionMatching
       self._DoRetricSearchbyTiltAngle=_DoRetricSearchbyTiltAngle
       self._Tilt0=_Tilt0
       self._TiltF=_TiltF
       self._ProjMatchingExtra=_ProjMatchingExtra
       ##self._AngSamplingRateDeg=str(_AngSamplingRateDeg)
       ##self._MaxChangeOffset=str(_MaxChangeOffset)
       ##self._MaxChangeInAngles=str(_MaxChangeInAngles)
       self._DisplayMask=_DisplayMask
       #absolute path to main project dir
       self._ProjectDir=_ProjectDir
       #self._DoAlign2D=_DoAlign2D
       self._InnerRadius=_InnerRadius
       self._OuterRadius=_OuterRadius
       self._Align2DIterNr=_Align2DIterNr
       self._DisplayAlign2D=_DisplayAlign2D
       self._Align2DExtraCommand=_Align2DExtraCommand
       #self._Align2dMaxChangeOffset=_Align2dMaxChangeOffset
       #self._Align2dMaxChangeRot=_Align2dMaxChangeRot
       self._DisplayReconstruction=_DisplayReconstruction
       self._DisplayResolution=_DisplayResolution
       self._DoReconstruction=_DoReconstruction
       self._DoComputeResolution=_DoComputeResolution
       self._ResolSam=_ResolSam
       self._DoParallel=_DoParallel
       self._MyNumberOfCPUs=_MyNumberOfCPUs
       if (len(_Symfile)>0):
           self._Symfile=os.path.abspath(str(_Symfile))
       else:
           self._Symfile=_Symfile
       #self._iteration_number=1
       #self._ARTLambda=_ARTLambda
       self._ARTReconstructionExtraCommand=_ARTReconstructionExtraCommand
       self._WBPReconstructionExtraCommand=_WBPReconstructionExtraCommand
       #self._ReconstructionMethod=_ReconstructionMethod
       self._SetResolutiontoZero=_SetResolutiontoZero
       #self._ConstantToAddToFiltration=_ConstantToAddToFiltration
       if (_MyMachineFile[0]=="$"):
           self._MyMachineFile=_MyMachineFile
       else:
           self._MyMachineFile=os.path.abspath(_MyMachineFile)
       
       self._user_suplied_ReferenceVolume=self._ReferenceFileName
       self._Proj_Maching_Output_Root_Name=_Proj_Maching_Output_Root_Name
       self._multi_align2d_sel=_multi_align2d_sel
       self._AlignWithClassReference=_AlignWithClassReference
       self._align2d_sel=_align2d_sel
       self._align2d_doc=_align2d_doc
       #name of the masked volume
       tmp_OutPutVol=os.path.basename(self._ReferenceFileName)
       ReferenceFileName_without_ext=\
            (os.path.splitext(str(os.path.basename(tmp_OutPutVol))))[0]
       
       
       self._mylog=log.init_log_system(_ProjectDir,
                                       _LogDir,
                                       sys.argv[0],
                                       _WorkDirectory)
                                      
       #uncomment next line to get Degub level logging
       self._mylog.setLevel(logging.DEBUG)
       self._mylog.debug("Debug level logging enabled")
                                      
       #_ContinueAtIteration-=1;
       _NumberofIterations +=1;
       if _ContinueAtIteration!=1 and DoDeleteWorkingDir==True:
          print "You can not delete the working directory"
          print " and star at iteration", _ContinueAtIteration
          exit(1)
       if (DoDeleteWorkingDir): 
          delete_working_directory(self._mylog,self._WorkDirectory)
       else:
          self._mylog.info("Skipped DoDeleteWorkingDir") 
       create_working_directory(self._mylog,self._WorkDirectory)
       #made backup of this script and visualize
       log.make_backup_of_script_file(sys.argv[0],self._WorkDirectory)
       log.make_backup_of_script_file(AnalysisScript,self._WorkDirectory)
       
       #copy files to local directory
       if _ContinueAtIteration==1:
          copy_images_to_local_disk(self._mylog,
                                    self._SelFileName,
                                    self._WorkDirectory,
                                    images_dir)
       
       
       #change to working dir
       # WorkDirectory='ProjMatch/Test1'
       os.chdir(self._WorkDirectory)
       self._SelFileName=self._WorkDirectory+'/'+\
                         str(os.path.basename(self._SelFileName))
       import arg                 
       ##
       ##LOOP
       ##
       #output of reconstruction cicle
       #first value given by user
       #this names are the input of the mask program
       #in general is the output of the reconstruction plus filtration
       self._ReconstructedVolume=[]
       fill_name_vector("",
                        self._ReconstructedVolume,
                        _NumberofIterations,
                        ReconstrucedVolume)
                        
       self._ReconstructedandfilteredVolume=[]
       fill_name_vector(self._user_suplied_ReferenceVolume,
                        self._ReconstructedandfilteredVolume,
                        _NumberofIterations,
                        Filtered_Image)

       #reconstructed and filtered volume of n-1 after masking
       #called reference volume
       self._ReferenceVolume=[]
       fill_name_vector("",self._ReferenceVolume,_NumberofIterations,
                        ReferenceVolume)

       for _iteration_number in range(_ContinueAtIteration, _NumberofIterations):
          debug_string =  "ITERATION: " +  str(_iteration_number)
          print "*", debug_string
          self._mylog.info(debug_string)
          self._ResetImageHeader=arg.getComponentFromVector(_ResetImageHeader,\
                                                        _iteration_number-1)
          if (int(self._ResetImageHeader)):
              reset_image_header(self._mylog,self._SelFileName)
       
          #create working dir for this iteration
          Iteration_Working_Directory=self._WorkDirectory+'/Iter_'+\
                                      str(_iteration_number)
          create_working_directory(self._mylog,Iteration_Working_Directory)
          #change to iteration_working dir
          os.chdir(Iteration_Working_Directory)

          #save initial header
          #docfile_with_original_angles='docfile_with_original_angles.doc'
          command='xmipp_header_extract -i '
          command=command+self._SelFileName + \
                         ' -o ' + docfile_with_original_angles 
          self._mylog.debug("save original header in file: " + 
                                   docfile_with_original_angles ) 
          self._mylog.info(command)
          os.system(command)
          execute_mask(_DoMask,
                       self._mylog,
                       self._ProjectDir,
                       self._ReconstructedandfilteredVolume[_iteration_number],#in
                       self._MaskFileName,
                       self._DisplayMask,
                       _iteration_number,
                       self._ReferenceVolume[_iteration_number])#out
          #parameters for projection matching
          self._AngSamplingRateDeg=arg.getComponentFromVector(_AngSamplingRateDeg,\
                                                        _iteration_number-1)
          self._MaxChangeOffset=arg.getComponentFromVector(_MaxChangeOffset,\
                                                        _iteration_number-1)
          self._MaxChangeInAngles=arg.getComponentFromVector(_MaxChangeInAngles,\
                                                        _iteration_number-1)
          if (_DoProjectionMatching):
             #parameters for projection matching
             self._AngSamplingRateDeg=arg.getComponentFromVector(_AngSamplingRateDeg,\
                                                           _iteration_number-1)
             self._MaxChangeOffset=arg.getComponentFromVector(_MaxChangeOffset,\
                                                           _iteration_number-1)
             self._MaxChangeInAngles=arg.getComponentFromVector(_MaxChangeInAngles,\
                                                           _iteration_number-1)
             execute_projection_matching(self._mylog,
                                         self._ProjectDir,
                                         self._ReferenceVolume[_iteration_number],#in
                                         self._MaskFileName,
                                         self._SelFileName,
                                         self._Proj_Maching_Output_Root_Name, 
                                         self._AngSamplingRateDeg,
                                         self._DoRetricSearchbyTiltAngle,
                                         self._Tilt0,
                                         self._TiltF,
                                         self._InnerRadius,
                                         self._OuterRadius,
                                         self._MaxChangeOffset, 
                                         self._MaxChangeInAngles,  
                                         self._ProjMatchingExtra,
                                         self._DisplayProjectionMatching,
                                         self._DoParallel,
                                         self._MyNumberOfCPUs,
                                         self._MyMachineFile,
                                         self._WorkDirectory,
                                         self._Symfile,
                                         _iteration_number
                                         )
          else:
             self._mylog.info("Skipped ProjectionMatching") 
          #if (_DoAlign2D):
          self._Align2dMaxChangeOffset=arg.getComponentFromVector(_Align2dMaxChangeOffset,\
                                                        _iteration_number-1)
          self._Align2dMaxChangeRot=arg.getComponentFromVector(_Align2dMaxChangeRot,\
                                                        _iteration_number-1)
          self._DoAlign2D=arg.getComponentFromVector(_DoAlign2D,\
                                                           _iteration_number-1)
          execute_align2d(self._mylog,
                             self._InnerRadius,
                             self._OuterRadius,
                             self._Align2DIterNr,
                             self._Align2DExtraCommand,
                             self._Align2dMaxChangeOffset,
                             self._Align2dMaxChangeRot,
                             self._Proj_Maching_Output_Root_Name,
                             _iteration_number,
                             self._DoParallel,
                             self._MyNumberOfCPUs,
                             self._MyMachineFile,
                             self._DisplayAlign2D,
                             self._multi_align2d_sel,
                             self._AlignWithClassReference,
                             self._DoAlign2D
                             )
          #else:
             #self._mylog.info("Skipped Align2D") 
          #if (_DoAlign2D):
          command='cat ' + self._multi_align2d_sel \
                         + ' | grep  \.med\.xmp  ' \
                         + ' | grep -v \ -1 >' \
                         + self._align2d_sel
          self._mylog.info(command)
          os.system(command)
          command='xmipp_header_extract -i ' + self._align2d_sel +\
                                      ' -o ' + self._align2d_doc
          self._mylog.info(command)
          os.system(command)
          #else:
          
          self._ReconstructionMethod=arg.getComponentFromVector(_ReconstructionMethod,\
                                                        _iteration_number-1)
          self._ARTLambda=arg.getComponentFromVector(_ARTLambda,\
                                                        _iteration_number-1)
          if (_DoReconstruction):
             execute_reconstruction(self._mylog, 
                                    self._ARTReconstructionExtraCommand,
                                    self._WBPReconstructionExtraCommand,
                                    _iteration_number,
                                    self._DisplayReconstruction,
                                    self._DoParallel,
                                    self._MyNumberOfCPUs,
                                    self._MyMachineFile,
                                    self._ReconstructionMethod,
                                    self._ARTLambda,
                                    self._align2d_sel,
                                    self._Symfile,
                                    self._ReconstructedVolume[_iteration_number]
                                    )
          else:
             self._mylog.info("Skipped Reconstruction") 
          
          if (_DoComputeResolution):
              filter_frequence=execute_resolution(self._mylog,
                           self._ARTReconstructionExtraCommand,
                           self._WBPReconstructionExtraCommand,
                           self._ReconstructionMethod,
                           _iteration_number,
                           self._DisplayReconstruction,
                           self._ResolSam,
                           self._align2d_sel,
                           self._DoParallel,
                           self._MyNumberOfCPUs,
                           self._MyMachineFile,
                           self._Symfile,
                           self._DisplayResolution,
                           self._ARTLambda)

          else:
             self._mylog.info("Skipped Resolution") 
          
          self._ConstantToAddToFiltration=arg.getComponentFromVector(\
                                               ConstantToAddToFiltration,\
                                                        _iteration_number-1)
          filter_at_given_resolution(_DoComputeResolution,
                                  self._mylog, 
                                  _iteration_number,
                                  self._SetResolutiontoZero,
                                  self._ConstantToAddToFiltration,
                                  filter_frequence,
                                  self._ReconstructedVolume[_iteration_number],
                                  self._ReconstructedandfilteredVolume[1+_iteration_number]
                                 )


#------------------------------------------------------------------------
#delete_working directory
#------------------------------------------------------------------------
def delete_working_directory(_mylog,_WorkDirectory):
    import os
    import shutil
    print '*********************************************************************'
    print '* Delete working directory tree'
    _mylog.info("Delete working directory tree")

    if os.path.exists(_WorkDirectory):
       shutil.rmtree(_WorkDirectory)
       
#------------------------------------------------------------------------
#create_working directory
#------------------------------------------------------------------------
def create_working_directory(_mylog,_WorkDirectory):
    import os
    print '*********************************************************************'
    print '* Create directory ' + _WorkDirectory 
    _mylog.info("Create working directory " + _WorkDirectory )

    if not os.path.exists(_WorkDirectory):
       os.makedirs(_WorkDirectory)

#------------------------------------------------------------------------
#           copy files to local dir
#------------------------------------------------------------------------
def copy_images_to_local_disk(_mylog,_SelFileName,_WorkDirectory,images_dir):
      import os,selfile
      print '*********************************************************************'
      print '* Copying images to working directory ...'
      mysel=selfile.selfile()
      mysel.read(_SelFileName)
      newsel=mysel.copy_sel(_WorkDirectory+'/'+images_dir)
      newsel.write(_WorkDirectory+'/'+_SelFileName)
      newsel.make_abspath()
      _mylog.info("copy files to local directory")

#------------------------------------------------------------------------
#           reset_image_header(self._SelFileName)
#------------------------------------------------------------------------
def reset_image_header(_mylog,_SelFileName):
    import os
    print '**************************************************************'
    print '* Reset headers of files in selfile ' + _SelFileName 
    command = "xmipp_header_reset -i " + _SelFileName
    print '* ',command
    _mylog.info(command)
    os.system(command)

#------------------------------------------------------------------------
#execute_mask
#------------------------------------------------------------------------
def execute_mask(_DoMask,
                 _mylog,
                 _ProjectDir,
                 _ReferenceFileName,
                 _MaskFileName,
                 _DisplayMask,
                 _iteration_number,
                 _ReferenceVolume):
   import os,shutil
   _mylog.debug("execute_mask")
   if(_iteration_number==1):
      InPutVolume=_ReferenceFileName
   else:   
      InPutVolume=_ReferenceFileName+".vol"
   ##InPutVolume=InPutVolume.replace(".vol.vol",".vol")
   if (_DoMask):
       MaskVolume =_MaskFileName
       MaskedVolume=_ReferenceVolume
       print '*********************************************************************'
       print '* Mask the reference volume'
       command='xmipp_mask'+ \
               ' -i '    + InPutVolume + \
               ' -o '    + _ReferenceVolume + \
               ' -mask ' + MaskVolume 

       print '* ',command
       _mylog.info(command)
       os.system(command)
       if _DisplayMask==True:
          command='xmipp_show -vol '+ MaskedVolume +' -w 10 &'
          print '*********************************************************************'
          print '* ',command
          _mylog.info(command)
          os.system(command)
          
   else:
       shutil.copy(InPutVolume,_ReferenceVolume)
       _mylog.info("Skipped Mask")
       _mylog.info("cp" + InPutVolume +\
                   " "  + _ReferenceVolume )
       print '*********************************************************************'
       print '* Skipped Mask'

#------------------------------------------------------------------------
#execute_projection_matching
#------------------------------------------------------------------------
def execute_projection_matching(_mylog,
                                _ProjectDir,
                                _ReferenceVolume,
                                _MaskFileName,
                                _SelFileName,
                                _Proj_Maching_Output_Root_Name,    
                                _AngSamplingRateDeg,
                                _DoRetricSearchbyTiltAngle,
                                _Tilt0,
                                _TiltF,
                                _Ri,
                                _Ro,
                                _MaxChangeOffset,
                                _MaxChangeInAngles,   
                                _ProjMatchingExtra,
                                _DisplayProjectionMatching,
                                _DoParallel,
                                _MyNumberOfCPUs,
                                _MyMachineFile,
                                _WorkDirectory,
                                _Symfile,
                                _iteration_number):
                                           
   _mylog.debug("execute_projection_matching")
   import os,selfile
   Reference_Vol=_ReferenceVolume

   print '*********************************************************************'
   print '* Asign projection direction'
   parameters=' -i '           + _SelFileName         + \
              ' -vol '         + Reference_Vol + \
              ' -o '           + _Proj_Maching_Output_Root_Name      + \
              ' -sam '         + _AngSamplingRateDeg  + \
              ' -max_shift '   + _MaxChangeOffset
  
   if _DoRetricSearchbyTiltAngle:
     parameters=  parameters                           + \
              ' -tilt0 '       + str(_Tilt0)        + \
              ' -tiltF '       + str(_TiltF)
  
   parameters=  parameters                           + \
              ' -Ri '          + str(_Ri)           + \
              ' -Ro '          + str(_Ro)           + \
              ' -output_refs -output_classes ' + \
              ' -ang_search ' + str(_MaxChangeInAngles) +\
              ' '              + _ProjMatchingExtra
   if len(_Symfile)>1:
      parameters = parameters + " -sym " + _Symfile
   # -dont_modify_header          
   import launch_parallel_job
   RunInBackground=False
   launch_parallel_job.launch_job(
                       _DoParallel,
                       'xmipp_angular_projection_matching',
                       'xmipp_mpi_angular_projection_matching',
                       parameters,
                       _mylog,
                       _MyNumberOfCPUs,
                       _MyMachineFile,
                       RunInBackground)
   #make absolute path so visualization protocol can be run from
   #the same directory
   
   classes_sel_file=selfile.selfile()
   classes_sel_file.read(_Proj_Maching_Output_Root_Name+'_classes.sel')
   library_sel_file=selfile.selfile()
   library_sel_file.read(_Proj_Maching_Output_Root_Name+'_lib.sel')
   newsel=library_sel_file.intercalate_union(classes_sel_file)
   compare_sel_file=_Proj_Maching_Output_Root_Name+'_compare.sel'
   newsel=newsel.make_abspath()
   newsel.write(compare_sel_file)
   if _DisplayProjectionMatching==True:
      command='xmipp_show -sel '+ "../"+'Iter_'+\
                   str(_iteration_number) +'/'+ compare_sel_file +' -w 10 &'
      #NOTE, selection will be made in next showsel
      print '*********************************************************************'
      print '* ',command
      _mylog.info(command) 
      os.system(command)

#------------------------------------------------------------------------
#execute_align2d
#read all class*.sel files
#align then with the class*.xmp images
#save in class*_aligned.xmp
#------------------------------------------------------------------------
def execute_align2d(_mylog,
                    _InnerRadius,
                    _OuterRadius,
                    _Align2DIterNr,
                    _Align2DExtraCommand,
                    _Align2dMaxChangeOffset,
                    _Align2dMaxChangeRot,
                    _Proj_Maching_Output_Root_Name,
                    _iteration_number,
                    _DoParallel,
                    _MyNumberOfCPUs,
                    _MyMachineFile,
                    _DisplayAlign2D,
                    _multi_align2d_sel,
                    _AlignWithClassReference,
                    _DoAlign2D
                    ):
   #proj_match_class00001.corr
   #proj_match_class00001.med.xmp --> average class projection AFTER alig2d
   #                                  has weight *
   #proj_match_class00001.sig.xmp    standard deviation
   #proj_match_class00001.ref.xmp --> teporal file produed by alig2d
   #                                   with iteration reference
   #                                   this image is filtered so may lok nicer 
   #                                   than true average
   #proj_match_class00001.sel  -- > particles asigned to reference 00001
   #proj_match_class00001.xmp -- > average class projection BEFORE align *
   #proj_match_lib00001.proj -- > gallery of projections*
   #if secuential execute orden if not store it in a file
   import tempfile,shutil
   import spider_header 
   one = 1. 
   tmp_file_name = _Proj_Maching_Output_Root_Name +".tmp"
   if _DoParallel:
        fh = open(tmp_file_name,'w')
      
   _mylog.debug("execute_align2d")
   _mylog.debug("_DoAlign2D " + _DoAlign2D)
   import os,selfile, glob,shutil
   class_sel_pattern = _Proj_Maching_Output_Root_Name+'_class[0-9]*.sel'
   _mylog.debug("class_sel_pattern: " + class_sel_pattern)
   aux_sel_file=selfile.selfile()
   align2d_sel=selfile.selfile()
   #create compare sel 
   if(_DoAlign2D=='1'): 
       for class_selfile in glob.glob(class_sel_pattern):
          lib_file_name = class_selfile.replace('class','lib')
          lib_file_name = lib_file_name.replace('.sel','.proj') 
          class_file_name = class_selfile.replace('.sel','.xmp')
          if(_AlignWithClassReference):
             reference=lib_file_name
          else:  
             reference=class_file_name
          aux_sel_file.read(class_selfile)
          #first line in sel must be active

          align2d_sel.insert(lib_file_name,str(1))
          align2d_sel.insert(class_file_name,str(-1))
          aligned_file_name = class_selfile.replace('.sel','.med.xmp')
          align2d_sel.insert(aligned_file_name,str(1))
          if (aux_sel_file.length()<1):
             command="xmipp_operate -i " + \
                      class_file_name + " -mult 0 -o " + aligned_file_name +"\n"
             print '* ',command
          elif (aux_sel_file.length()==1):
             _mylog.info("cp "+aux_sel_file.sellines[0][0]+" "+aligned_file_name) 
             shutil.copy(aux_sel_file.sellines[0][0],aligned_file_name)
             # this images must  weight one
             spider_header.set_header_position(aligned_file_name,
                                               260,
                                               one
                                               )
             command=""                                  
          else:
             print '*********************************************************************'
             print '* Aligning translationally each class'
             command='xmipp_align2d'+ \
                     ' -i '  + class_selfile + \
                     ' -Ri ' +   str(_InnerRadius) + \
                     ' -Ro ' +   str(_OuterRadius) +\
                     ' -iter ' + str(_Align2DIterNr) +\
                     ' -ref ' + reference
             if  len(_Align2dMaxChangeOffset)>1 :     
               command +=' -max_shift '  + _Align2dMaxChangeOffset
             if  len(_Align2dMaxChangeRot)>1 :     
               command +=' -max_rot '  + _Align2dMaxChangeRot
             if  len(_Align2DExtraCommand)>1 :     
               command += ' '  + _Align2DExtraCommand
             print '* ',command
          if (len(command)>2):
             if _DoParallel:
                fh.writelines(command+"\n")
                _mylog.debug(command)
             else:  
                os.system(command)
                _mylog.info(command)
       align2d_sel=align2d_sel.make_abspath()         
       align2d_sel.write(_multi_align2d_sel);
       if _DoParallel:
           fh.close();
           parameters="-i " +  tmp_file_name
           _mylog.info("xmipp_mpi_run " + parameters)
           import launch_parallel_job
           RunInBackground=False
           launch_parallel_job.launch_parallel_job(
                               'xmipp_mpi_run',
                               parameters,
                               _mylog,
                               _MyNumberOfCPUs,
                               _MyMachineFile,
                               RunInBackground)
           #os.remove(tmp_file_name)
   else:
       for class_selfile in glob.glob(class_sel_pattern):
          #projection matching average
          inputfile  = class_selfile.replace('.sel','.xmp')
          #alig2n2d output
          outputfile = class_selfile.replace('.sel','.med.xmp')
          shutil.copy(inputfile,outputfile)
          aux_sel_file.read(class_selfile)
          weight=aux_sel_file.length()
          print "weight", weight
          _mylog.info("cp "+inputfile+" "+outputfile+", weight="+str(weight)) 
          # this images must  weight one
          spider_header.set_header_position(outputfile,
                                            260,
                                            weight
                                            )
          
          lib_file_name = class_selfile.replace('class','lib')
          lib_file_name = lib_file_name.replace('.sel','.proj') 
          align2d_sel.insert(lib_file_name,str(1))
          align2d_sel.insert(inputfile,str(-1))
          align2d_sel.insert(outputfile,str(1))
       align2d_sel=align2d_sel.make_abspath()         
       align2d_sel.write(_multi_align2d_sel);

#   #if a sel file is empty copy reference class
#   for class_selfile in glob.glob(class_sel_pattern):
#      aux_sel_file.read(class_selfile)
#      message= aux_sel_file.selfilename,\
#             "length ", aux_sel_file.length(), \
#             "length_t ", aux_sel_file.length_even_no_actives()
#      _mylog.debug(message)
#             
#      if (aux_sel_file.length()<1 and \
#          aux_sel_file.length_even_no_actives()>0):  
#          class_file_name = class_selfile.replace('.sel','.xmp') 
#          aligned_file_name = class_selfile.replace('.sel','.med.xmp')
#          shutil.copy(class_file_name,aligned_file_name)
#          _mylog.info("cp "+class_file_name+" "+aligned_file_name) 

   if _DisplayAlign2D==True:
      command='xmipp_show -showall -sel '+ "../"+'Iter_'+\
                   str(_iteration_number) +'/'+ _multi_align2d_sel +' -w 9'
      print '*********************************************************************'
      print '* ',command
      _mylog.info(command) 
      os.system(command)

   #Set right angle  for averages computed with align
   #we may read using headerinfo.
   _in_filename   =  _Proj_Maching_Output_Root_Name + '_classes.sel'
   _out_filename  =  _Proj_Maching_Output_Root_Name + '_classes.doc'
   command='xmipp_header_extract -i ' + _in_filename +\
                               ' -o ' + _out_filename
   _mylog.info(command)
   os.system(command)


   _input = open(_out_filename)
   _output = open(_out_filename+'.o', 'w')
   for s in _input.xreadlines():
        _output.write(s.replace("xmp", "med.xmp"))
   
   _input.close()
   _output.close()
   command='xmipp_header_assign -i ' + _out_filename+'.o'
   _mylog.info(command)
   os.system(command)
   
   

#------------------------------------------------------------------------
#execute_reconstruction
#------------------------------------------------------------------------
def execute_reconstruction(_mylog,
                           _ARTReconstructionExtraCommand,
                           _WBPReconstructionExtraCommand,
                           _iteration_number,
                           _DisplayReconstruction,
                           _DoParallel,
                           _MyNumberOfCPUs,
                           _MyMachineFile,
                           _ReconstructionMethod,
                           _ARTLambda,
                           _align2d_sel,
                           _Symfile,
                           _ReconstructedandfilteredVolume):
   import os,shutil
   _mylog.debug("execute_reconstruction")


   #the user should be able to delete images
   #this is dirty but what else can I do
   reconstruction_sel=_align2d_sel;
   Outputvolume  =_ReconstructedandfilteredVolume

   print '*********************************************************************'
   print '* Reconstruct volume using '
   if _ReconstructionMethod=='wbp':
      Outputvolume  =Outputvolume+".vol"
      program = 'xmipp_reconstruct_wbp'
      mpi_program = 'xmipp_mpi_reconstruct_wbp'
      parameters= ' -i '    + reconstruction_sel + \
                  ' -o '    + Outputvolume + \
                  ' -weight -use_each_image '
      if len(_Symfile)>1:
         parameters = parameters + ' -sym ' + _Symfile + ' '
      parameters = parameters + _WBPReconstructionExtraCommand
              
   elif _ReconstructionMethod=='art':
      program = 'xmipp_reconstruct_art'
      mpi_program = 'NULL'
      _DoParallel=False
      parameters=' -i '    + reconstruction_sel + \
                 ' -o '    + Outputvolume + ' ' + \
                 ' -WLS '
      if len(_Symfile)>1:
         parameters = parameters + ' -sym ' + _Symfile + ' '
      if len(_ARTLambda)>1:
         parameters = parameters + ' -l '   + _ARTLambda + ' '
      parameters = parameters + _ARTReconstructionExtraCommand
   else:
      _mylog.error("Reconstruction method unknown. Quiting")
      print "Reconstruction method unknown. Quiting"
      exit(1)
    
   import launch_parallel_job
   RunInBackground=False
   launch_parallel_job.launch_job(
                       _DoParallel,
                       program,
                       mpi_program,
                       parameters,
                       _mylog,
                       _MyNumberOfCPUs,
                       _MyMachineFile,
                       RunInBackground)



   #_mylog.info(command+ ' ' + parameters)
   if _DisplayReconstruction==True:
      command='xmipp_show -vol '+ Outputvolume + '&'
      print '*********************************************************************'
      print '* ',command
      _mylog.info(command)
      os.system(command)

#------------------------------------------------------------------------
#           execute_resolution(self._SelFileName)
#------------------------------------------------------------------------
def  execute_resolution(_mylog,
                        _ARTReconstructionExtraCommand,
                        _WBPReconstructionExtraCommand,
                        _ReconstructionMethod,
                        _iteration_number,
                        _DisplayReconstruction,
                        _ResolSam,
                        _align2d_sel,
                        _DoParallel,
                        _MyNumberOfCPUs,
                        _MyMachineFile,
                        _Symfile,
                        _DisplayResolution,
                        _ARTLambda):
    import os
    # Split selfiles
    #the user should be able to delete images
    #this is dirty but what else can I do
    reconstruction_sel=_align2d_sel;
    split_sel_root_name="split_sel"
    command='xmipp_selfile_split -o '+ split_sel_root_name+ \
             ' -i ' + reconstruction_sel + \
             ' -n 2'
    print '* ',command
    _mylog.info(command)
    os.system(command)

    Outputvolumes=[]
    Outputvolumes.append(split_sel_root_name+'_1')
    Outputvolumes.append(split_sel_root_name+'_2')
    
    Selfiles=[]
    Selfiles.append(split_sel_root_name+'_1'+".sel")
    Selfiles.append(split_sel_root_name+'_2'+".sel")
    for i in range(len(Outputvolumes)):
       print '*********************************************************************'
       print '* Reconstruct volume'
       if _ReconstructionMethod=='wbp':
          RunInBackground=False
          program = 'xmipp_reconstruct_wbp'
          mpi_program = 'xmipp_mpi_reconstruct_wbp'
          parameters= ' -i '    + Selfiles[i] + \
                      ' -o '    + Outputvolumes[i] +".vol"\
                  ' -weight -use_each_image '
          if len(_Symfile)>1:
             parameters = parameters + ' -sym ' + _Symfile + ' '
          parameters = parameters + _WBPReconstructionExtraCommand

       elif _ReconstructionMethod=='art':
          RunInBackground=False
          program = 'xmipp_reconstruct_art'
          mpi_program = 'NULL'
          _DoParallel=False
          parameters=' -i '    + Selfiles[i] + \
                     ' -o '    + Outputvolumes[i] +\
                     ' -WLS '
          if len(_Symfile)>1:
             parameters = parameters + ' -sym ' + _Symfile + ' '
          if len(_ARTLambda)>1:
             parameters = parameters + ' -l '   + _ARTLambda + ' '
          parameters = parameters + _ARTReconstructionExtraCommand
       else:
          _mylog.error("Reconstruction method unknown. Quiting")
          print "Reconstruction method unknown. Quiting"
          exit(1)

       import launch_parallel_job
       launch_parallel_job.launch_job(
                           _DoParallel,
                           program,
                           mpi_program,
                           parameters,
                           _mylog,
                           _MyNumberOfCPUs,
                           _MyMachineFile,
                           RunInBackground)


    
    for i in range(len(Outputvolumes)):
       Outputvolumes[i]+=".vol"
    print '**************************************************************'
    print '* Compute resolution ' 
    command = "xmipp_resolution_fsc -ref " + Outputvolumes[0] +\
              " -i " +Outputvolumes[1]  + ' -sam ' + str(_ResolSam)
    print '* ',command
    _mylog.info(command)
    os.system(command)
    import visualization
    if _DisplayResolution==True:
      plot=visualization.gnuplot()
      plot.plot_xy_file(Outputvolumes[1]+'.frc',
                          Title="Resolution",
                          X_Label="Armstrong^-1",
                          Y_Label="y",
                          X_col=1,
                          Y_col=2)
      print '*********************************************************************'
      print '* plot resolution'
      _mylog.info(" plot resolution")
    #compute resolution
    resolution_fsc_file = Outputvolumes[1]+'.frc'
    f = open(resolution_fsc_file, 'r')
    #skip first line
    fi=f.readline()
      
    filter_frequence=0. 
    for line in f:
        line = line.strip()
        if not line.startswith('#'):
            mylist = (line.split())
            if( float(mylist[1]) < 0.5):
               break
            else:
              filter_frequence=float(mylist[0])


    f.close()
    print '* maximum resolution (A^-1): ', filter_frequence
    filter_frequence *= _ResolSam
    print '* maximum resolution (px^-1): ', filter_frequence
    return filter_frequence

#------------------------------------------------------------------------
#           filter_at_given_resolution
#------------------------------------------------------------------------
def filter_at_given_resolution(_DoComputeResolution,
                               _mylog, 
                               _iteration_number,
                               _SetResolutiontoZero,
                               _ConstantToAddToFiltration,
                               _filter_frequence,
                               _ReconstructedVolume,
                               _ReconstructedandfilteredVolume
                               ):

    import os,shutil
    Inputvolume   =_ReconstructedVolume+'.vol'
    Outputvolume  =_ReconstructedandfilteredVolume+'.vol'
    if (_SetResolutiontoZero):
       filter_in_pixels_at = float(_ConstantToAddToFiltration)
    else:
       filter_in_pixels_at = float(_filter_frequence) +\
                             float(_ConstantToAddToFiltration)
    print '**************************************************************'
    print '* Filter reconstruction ' 
    if (_ConstantToAddToFiltration<0.5 or (not _DoComputeResolution)):
        shutil.copy(Inputvolume,Outputvolume) 
        command ="shutilcopy" + Inputvolume + ' ' + Outputvolume
    else:   
        command = "xmipp_fourier_filter -i " + Inputvolume +\
                  " -o " + Outputvolume + ' -low_pass ' +\
                  str (filter_in_pixels_at)
        print '* ',command
        os.system(command)
    _mylog.info(command)


def  fill_name_vector(_user_suplied_name,
                      _volume_name_list,
                      _NumberofIterations,
                      _root_name):
     _volume_name_list.append("dummy")
     if (len(_user_suplied_name)>1):
        _volume_name_list.append(_user_suplied_name)
     for _iteration_number in range(1, _NumberofIterations):
         _volume_name_list.append("../"+'Iter_'+\
                                   str(_iteration_number)+'/'+ 'Iter_'+\
                                   str(_iteration_number)+'_'+\
                                   _root_name)
                  
#
# main
#     
if __name__ == '__main__':

    # create rotational_spectra_class object
    # 
   
    #init variables
  my_projmatch=projection_matching_class(
                NumberofIterations,     
                ContinueAtIteration,  
                DoMask,                         
                DisplayMask,                    
                ReferenceFileName,              
                MaskFileName,                   
                DoProjectionMatching,           
                DisplayProjectionMatching,      
                AngSamplingRateDeg,             
                DoRetricSearchbyTiltAngle,      
                Tilt0,                          
                TiltF,                          
                ProjMatchingExtra,              
                MaxChangeOffset,
                MaxChangeInAngles,                
                DoAlign2D,                      
                InnerRadius,                    
                OuterRadius,                    
                Align2DIterNr,                  
                DisplayAlign2D,                 
                Align2DExtraCommand,
                Align2dMaxChangeOffset,
                Align2dMaxChangeRot,            
                DisplayReconstruction,
                DisplayResolution,          
                DoReconstruction,
                ReconstructionMethod,
                ARTLambda,
                ARTReconstructionExtraCommand,
                WBPReconstructionExtraCommand,
                DoComputeResolution,
                ResolSam,
                SelFileName,                    
                WorkDirectory,                  
                ProjectDir,                     
                LogDir,                         
                DoParallel,                     
                NumberOfCPUs,                   
                MachineFile,                    
                Symfile,                        
                ReferenceVolume,                
                Proj_Maching_Output_Root_Name,  
                multi_align2d_sel, 
                AlignWithClassReference,             
                align2d_sel,                    
                align2d_doc,                    
                ResetImageHeader,
                SetResolutiontoZero,
                ConstantToAddToFiltration
                )
