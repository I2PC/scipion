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
# required files: log.py, Sel_Files.py
#
#-----------------------------------------------------------------------------
# {section} Global parameters

#-----------------------------------------------------------------------------
# {expert} Root directory name for this project:
"""
All Paths should be relative to this directory
"""
ProjectDir="/home/roberto2/Test/PARA_Roberto"

# Selfile with the input images:
"""
Absolute paths are required in self file, self file localtion relative to 
ProjectDir
"""
SelFileName='../all.sel.new'

# Working directory: 
WorkDirectory='test1'

# Delete working directory if it already exists?
DoDeleteWorkingDir=True

# {expert} Directory name for logfiles:
LogDir="Logs"

# Reference file name (3D map)
ReferenceFileName="init_reference/LTA_rot_0.1_norm.vol"

# Batch submission command (use "" to launch without batch submission):
""" This will depend on your queueing system., ask your system administrator...

    Examples: LaunchJobCommand=\"bsub -q 1day\"
      or, if you do not use a queueing system: LaunchJobCommand=\"\"
"""
LaunchJobCommand="" 

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel? (see Expert options)
DoParallel=True
# Number of processors to use:
MyNumberOfCPUs=2
# {expert} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
"""
MyMachineFile="/home/roberto2/bin/machines.dat"
# {expert} Standard script to launch MPI-jobs:
""" This will also depend on your system...
    The simplest script consists of the following two lines:

    #!/usr/bin/env sh
    `which mpirun` -np MyNumberOfCPUs -machinefile MyMachineFile \

    Note that the next line with the xmipp_mpi_MLalign2D command will be
    generated automatically, and the variables MyNumberOfCPUs and MyMachineFile
    will be replaced by the corresponding values given here above.
    More scripts for different batch systems can be found at:
    [Wiki link]
"""
ParallelScript="/home/roberto2/bin/mpi.sh"

#-----------------------------------------------------------------------------
# {section} Mask
#-----------------------------------------------------------------------------
# Mask Reference volume
""" Masking the reference volume will increase the signal to noise ratio.
Do not provide a very tight mask.
"""
DoMask=True

#show masked volume
"""Masked volume will be shown.
"""
DisplayMask=True

#binary mask-file used to mask reference volume
MaskFileName='circular_mask.msk'


#-----------------------------------------------------------------------------
# {section} Projection Matching
#-----------------------------------------------------------------------------
# Projection Matching
""" 
"""
DoProjectionMatching=True

#show masked volume
"""show average of projections.
"""
DisplayProjectionMatching=True

#{expert} Root name for output projections
OutputRootName='root'

#Angular sampling rate (in degrees)
AngSamplingRateDeg=20

#Maximum change in origin offset (+/- pixels)
MaxChangeOffset=4

#Lower-value for restricted tilt angle search
Tilt0=40

#Higher-value for restricted tilt angle search
TiltF=140

#extra options for Projection_Maching
ProjMatchingExtra=""

#depretated use doparallel
#Use mpi version of projection matching
#ProjMatchingMPIon=True

#stopping rule
#projection_matching.1

#-----------------------------------------------------------------------------
# {section} Align2d
#-----------------------------------------------------------------------------
# Perform 2D alignment?
DoAlign2D=True

#display align2d results
DisplayAlign2D=True

# Inner radius for rotational correlation (also used to make spectra):
""" These values are in pixels from the image center
"""
InnerRadius=0
# Outer radius for rotational correlation
OuterRadius=18
# Number of align2d iterations to perform:
Align2DIterNr=4
# {expert} Additional align2d parameters
""" For a complete description, see http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d

  Examples:
  Align2DExtraCommand=\"-trans_only  -filter 10 -sampling 2 -max_shift 2 -max_rot 3\"
  Align2DExtraCommand=\"-max_shift 2 -max_rot 3\"

  consider filtering the images with \"-filter 10 -sampling 2\"
"""
Align2DExtraCommand="-max_shift 2"# -max_rot 10"

#-----------------------------------------------------------------------------
# {section} Reconstruction
#-----------------------------------------------------------------------------
# Perform 3D Reconstruction
DoReconstruction=False

#display reconstructed volume
DisplayReconstruction=True

# {expert} Additional reconstruction parameters
""" 
  eamples go here
"""
ReconstructionExtraCommand=""#"-max_shift 2 -max_rot 10"

#reconstructiom method
"""
wbp/Art
"""
ReconstructionMethod="xmipp_wbp"

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
class projection_matching_class:

   #init variables
   
   def __init__(self,
                _DoMask, 
                _DisplayMask,
                _ReferenceFileName,
                _MaskFileName,
                _DoProjectionMatching,
                _DisplayProjectionMatching,
                _OutputRootName,
                _AngSamplingRateDeg,
                _Tilt0,
                _TiltF,
                _ProjMatchingExtra,
                _MaxChangeOffset,
                _DoAlign2D,
                _InnerRadius,
                _OuterRadius,
                _Align2DIterNr,
                _DisplayAlign2D,
                _Align2DExtraCommand,
                _DisplayReconstruction,
                _DoReconstruction,
                _ReconstructionMethod,
                _ReconstructionExtraCommand,
                _SelFileName,
                _WorkDirectory,
                _ProjectDir,
                _LogDir,
                _DoParallel,
                _MyNumberOfCPUs,
                _MyMachineFile,
                _LaunchJobCommand,
                _ParallelScript
                ):
       import os,sys
       scriptdir=os.path.expanduser('~')+'/trunk/Xmipp/Applications/Batches/Protocols_lib'
       sys.path.append(scriptdir) # add default search path
       import log,logging
       
       self._WorkDirectory=os.getcwd()+'/'+_WorkDirectory
       self._SelFileName=os.path.abspath(str(_SelFileName))
       self._ReferenceFileName=_ReferenceFileName
       self._MaskFileName=_MaskFileName
       self._DoMask=_DoMask
       self._DoProjectionMatching=_DoProjectionMatching
       self._DisplayProjectionMatching=_DisplayProjectionMatching
       self._OutputRootName=_OutputRootName
       self._AngSamplingRateDeg=str(_AngSamplingRateDeg)
       self._Tilt0=_Tilt0
       self._TiltF=_TiltF
       self._ProjMatchingExtra=_ProjMatchingExtra
       self._MaxChangeOffset=str(_MaxChangeOffset)
       self._DisplayMask=_DisplayMask
       self._ProjectDir=_ProjectDir
       self._DoAlign2D=_DoAlign2D
       self._InnerRadius=_InnerRadius
       self._OuterRadius=_OuterRadius
       self._Align2DIterNr=_Align2DIterNr
       self._DisplayAlign2D=_DisplayAlign2D
       self._Align2DExtraCommand=_Align2DExtraCommand
       self._DisplayReconstruction=_DisplayReconstruction
       self._DoReconstruction=_DoReconstruction
       self._DoParallel=_DoParallel
       self._MyNumberOfCPUs=_MyNumberOfCPUs
       self._MyMachineFile=_MyMachineFile
       self._ParallelScript=_ParallelScript
       self._LaunchJobCommand=_LaunchJobCommand
      
       self._iteration_number=0
       #name of the masked volume
       tmp_OutPutVol=os.path.basename(self._ReferenceFileName)
       ReferenceFileName_without_ext=\
            (os.path.splitext(str(os.path.basename(tmp_OutPutVol))))[0]
       self._Masked_Vol=ReferenceFileName_without_ext+'.msk'
       
       
       self._mylog=log.init_log_system(_ProjectDir,
                                      _LogDir,
                                      sys.argv[0],
                                      _WorkDirectory)
                                      
       #uncomment next line to get Degub level logging
       self._mylog.setLevel(logging.DEBUG)
       self._mylog.debug("Debug level logging enabled")
                                      
       if (DoDeleteWorkingDir): 
          delete_working_directory(self._mylog,self._WorkDirectory)
       else:
          self._mylog.info("Skipped DoDeleteWorkingDir") 
       create_working_directory(self._mylog,self._WorkDirectory)

       #change to working dir
       os.chdir(self._WorkDirectory)
       
       #save initial header
       self.original_header_file_name =  _OutputRootName + "_init.doc"
       command='xmipp_headerinfo -extract -i '
       command=command+(self._SelFileName+' -o ' + self.original_header_file_name)
       self._mylog.debug("save original header in file: " + 
                                self.original_header_file_name ) 
       self._mylog.info(command)
       os.system(command)

       if (_DoMask):
          execute_mask(self._mylog,
                       self._ProjectDir,
                       self._ReferenceFileName,
                       self._MaskFileName,
                       self._Masked_Vol,
                       self._DisplayMask,
                       self._iteration_number)
       else:
          self._mylog.info("Skipped Mask") 
          
       if (_DoProjectionMatching):
          execute_projection_matching(self._mylog,
                                      self._ProjectDir,
                                      self._ReferenceFileName,
                                      self._MaskFileName,
                                      self._Masked_Vol,
                                      self._SelFileName,
                                      self._OutputRootName,    
                                      self._AngSamplingRateDeg,
                                      self._Tilt0,
                                      self._TiltF,
                                      self._InnerRadius,
                                      self._OuterRadius,
                                      self._MaxChangeOffset,   
                                      self._ProjMatchingExtra,
                                      self._DisplayProjectionMatching,
                                      self._iteration_number,
                                      self._DoAlign2D,
                                      self._DoParallel,
                                      self._MyNumberOfCPUs,
                                      self._MyMachineFile,
                                      self._ParallelScript,
                                      self._LaunchJobCommand,
                                      self._WorkDirectory
                                      )
       else:
          self._mylog.info("Skipped ProjectionMatching") 
#
#                       
#       from protocol_projmatch_align2D import execute_align2d
#       if (_DoAlign2D):
#          execute_align2d(self._mylog,
#                          self._InnerRadius,
#                          self._OuterRadius,
#                          self._Align2DIterNr,
#                          self._Align2DExtraCommand,
#                          self._OutputRootName,
#                          self._iteration_number,
#                          self._DisplayAlign2D)
#       else:
#          self._mylog.info("Skipped Align2D") 
#
#        
##       from protocol_projmatch_Reconstruction import execute_Reconstruction
##       if (_DoReconstruction):
##          execute_Reconstruction(self._mylog,
##                      self._ReconstructionExtraCommand,
##                      self._iteration_number,
##                      self._DisplayReconstruction,
##                      self._DoReconstruction,
##                      self._ReconstructionMethod)
##       else:
##          self._mylog.info("Skipped Reconstruction") 


       #restore initial header
       command='xmipp_headerinfo -assign -i '
       command=command+(self.original_header_file_name)
       self._mylog.debug("restore original header in file: " + 
                                self.original_header_file_name ) 
       self._mylog.info(command)
       os.system(command)

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
    print '* Create working directory'
    _mylog.info("Create working directory")

    if not os.path.exists(_WorkDirectory):
       os.makedirs(_WorkDirectory)

#------------------------------------------------------------------------
#execute_mask
#------------------------------------------------------------------------
def execute_mask(_mylog,
                 _ProjectDir,
                 _ReferenceFileName,
                 _MaskFileName,
                 _Masked_Vol,
                 _DisplayMask,
                 _iteration_number):
   import os
   _mylog.debug("execute_mask")
   InPutVolume=_ProjectDir+'/'+_ReferenceFileName
   MaskVolume =_ProjectDir+'/'+_MaskFileName
   print '*********************************************************************'
   print '* Mask the reference volume'
   command='xmipp_mask'+ \
           ' -i '    + InPutVolume + \
           ' -o '    + _Masked_Vol + \
           ' -mask ' + MaskVolume 

   print '* ',command
   _mylog.info(command)
   os.system(command)
   if _DisplayMask==True:
      command='xmipp_show -vol '+ _Masked_Vol +' &'
      print '*********************************************************************'
      print '* ',command
      _mylog.info(command)
      os.system(command)
#------------------------------------------------------------------------
#execute_projection_matching
#------------------------------------------------------------------------
def execute_projection_matching(_mylog,
                                _ProjectDir,
                                _ReferenceFileName,
                                _MaskFileName,
                                _Masked_Vol,
                                _SelFileName,
                                _OutputRootName,    
                                _AngSamplingRateDeg,
                                _Tilt0,
                                _TiltF,
                                _Ri,
                                _Ro,
                                _MaxChangeOffset,   
                                _ProjMatchingExtra,
                                _DisplayProjectionMatching,
                                _iteration_number,
                                _DoAlign2D,
                                _DoParallel,
                                _MyNumberOfCPUs,
                                _MyMachineFile,
                                _ParallelScript,
                                _LaunchJobCommand,
                                _WorkDirectory):
                                           
   _mylog.debug("execute_projection_matching")
   import os,SelFiles
   ReferenceVolume=_ProjectDir+'/'+_ReferenceFileName

   print '*********************************************************************'
   print '* Asign projection direction'
   parameters=' -i '           + _SelFileName         + \
              ' -vol '         + _WorkDirectory + '/' + _Masked_Vol + \
              ' -o '           + _OutputRootName      + \
              ' -sam '         + _AngSamplingRateDeg  + \
              ' -max_shift '   + _MaxChangeOffset     + \
              ' -tilt0 '       + str(_Tilt0)        + \
              ' -tiltF '       + str(_TiltF)        + \
              ' -Ri '          + str(_Ri)           + \
              ' -Ro '          + str(_Ro)           + \
              ' -output_refs -output_classes ' + \
              ' '              + _ProjMatchingExtra
   # -dont_modify_header          
   if  _DoParallel==False:
      command='xmipp_projection_matching ' + parameters
      # -dont_modify_header 
      print '* ',command
      _mylog.info(command)
      os.system(command)
   else:
      import launch_parallel_job  
      launch_parallel_job.launch_parallel_job('xmipp_mpi_projection_matching',
                          parameters,
                          _ParallelScript,
                          _LaunchJobCommand,
                          _mylog,
                          _MyNumberOfCPUs,
                          _MyMachineFile,
                          _WorkDirectory,
                          False)
   
   if _DisplayProjectionMatching==True:
      classes_sel_file=SelFiles.selfile()
      classes_sel_file.read(_OutputRootName+'_classes.sel')
      library_sel_file=SelFiles.selfile()
      library_sel_file.read(_OutputRootName+'_lib.sel')
      newsel=library_sel_file.intercalate_union(classes_sel_file)
      compare_sel_file=_OutputRootName+'_compare.sel'
      newsel.write(compare_sel_file)
      command='xmipp_show -sel '+ compare_sel_file +' &'
      print '*********************************************************************'
      print '* ',command
      _mylog.info(command) 
      os.system(command)






#
# main
#     
if __name__ == '__main__':

    # create rotational_spectra_class object
    # 
   
    #init variables
  my_projmatch=projection_matching_class(
                DoMask, 
                DisplayMask,
                ReferenceFileName,
                MaskFileName,
                DoProjectionMatching,
                DisplayProjectionMatching,
                OutputRootName,
                AngSamplingRateDeg,
                Tilt0,
                TiltF,
                ProjMatchingExtra,
                MaxChangeOffset,
                DoAlign2D,
                InnerRadius,
                OuterRadius,
                Align2DIterNr,
                DisplayAlign2D,
                Align2DExtraCommand,
                DisplayReconstruction,
                DoReconstruction,
                ReconstructionMethod,
                ReconstructionExtraCommand,
                SelFileName,
                WorkDirectory,
                ProjectDir,
                LogDir,
                DoParallel,
                MyNumberOfCPUs,
                MyMachineFile,
                LaunchJobCommand,
                ParallelScript
                )
