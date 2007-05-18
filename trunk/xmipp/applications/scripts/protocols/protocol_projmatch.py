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
# {file} Selfile with the input images:
SelFileName='100.sel'

# {file} Initial 3D reference map:
ReferenceFileName='init_reference/LTA_rot_0.1_norm.vol'

# Working subdirectory: 
WorkDirectory='ProjMatch/Test'

# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=False

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir='/home/roberto2/Test/PARA_Roberto'

# {expert} Directory name for logfiles:
LogDir='Logs'

#-----------------------------------------------------------------------------
# {section} Mask
#-----------------------------------------------------------------------------
# Mask Reference volume
""" Masking the reference volume will increase the signal to noise ratio.
    Do not provide a very tight mask.
"""
DoMask=False

# Show masked volume
""" Masked volume will be shown. Do not set ths option to true for non iterative processing (jobs sent to queues)
"""
DisplayMask=False

# {file} Binary mask-file used to mask the reference volume
MaskFileName=''

#-----------------------------------------------------------------------------
# {section} Projection Matching
#-----------------------------------------------------------------------------
# Projection Matching
DoProjectionMatching=False

# Show projection maching library and classes
""" Show average of projections. Do not set ths option to true for non iterative processing (jobs sent to queues)
"""
DisplayProjectionMatching=False

# Angular sampling rate (in degrees)
AngSamplingRateDeg=8

# Maximum change in origin offset (+/- pixels)
MaxChangeOffset=10

# Restrict search by tilt angle
DoRetricSearchbyTiltAngle=False

# Lower-value for restricted tilt angle search
Tilt0=40

# Higher-value for restricted tilt angle search
TiltF=140

# {file} Symmetry description file:
""" This option does not work in combination with a limited search 
    of the rot or tilt angle.

    See WIKI link for a description of the symmetry file format
    dont give anything, if no symmetry is present
"""
Symfile=''

# Extra options for Projection_Matching
""" If you want to use your own references use the -ref option here,
    references name should be proj_match_lib00001.proj
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
"""
DoAlign2D=False

# Display 2D alignment results
""" Show aligned classes. Do not set ths option to true for non iterative processing (jobs sent to queues)
"""
DisplayAlign2D=False

# Inner radius for rotational correlation:
""" These values are in pixels from the image center
"""
InnerRadius=0

# Outer radius for rotational correlation
OuterRadius=18

# Number of align2d iterations (use at least 3):
Align2DIterNr=2

# {expert} Additional align2d parameters
""" For a complete description, see http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d

  Examples:
  Align2DExtraCommand=\"-trans_only  -filter 10 -sampling 2 -max_shift 2 -max_rot 3\"
  Align2DExtraCommand=\"-max_shift 2 -max_rot 3\"

  consider filtering the images with \"-filter 10 -sampling 2\"
"""
Align2DExtraCommand='-max_shift 4'

#-----------------------------------------------------------------------------
# {section} 3D Reconstruction
#-----------------------------------------------------------------------------
# Perform 3D Reconstruction
DoReconstruction=True

# Display reconstructed volume
DisplayReconstruction=False

# {expert} Additional reconstruction parameters
""" examples go here
"""
ReconstructionExtraCommand=''#'-max_shift 2 -max_rot 10'

# Reconstruction method
""" Choose between wbp or art
"""
ReconstructionMethod='wbp'

#-----------------------------------------------------------------------------
# {section} Resolution calculation
#-----------------------------------------------------------------------------
# Perform resolution calculation
DoComputeResolution=True

# {expert} Pixel size (in Ang.)
""" This will yield resolution plots where the X-axis units are 1/Armstrong instead of digital frequencies.
"""
ResolSam=1

#-----------------------------------------------------------------------------
# {section} Cleaning temporal files and Reseting original data
#-----------------------------------------------------------------------------
## ROBERTO: PLEASE EXPLAIN MORE CLEARLY WHAT YOU MEAN HERE!
# Reset image header
""" Set shifts and angles stored in the image headers to zero
"""
ResetImageHeader=True

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel?
DoParallel=True

# Number of processors to use:
NumberOfCPUs=2

# {file} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
    if your queueing system using an environment variable, give it here (with the leading $, e.g. $PBS_NODEFILE
"""
MachineFile="/home/roberto2/bin/machines.dat"

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_projmatch.py"

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

class projection_matching_class:

   #init variables
   
   def __init__(self,
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
                _align2d_sel,
                _align2d_doc,
                _ResetImageHeader
                ):
       import os,sys
       scriptdir=os.path.expanduser('~')+'/scripts/'
       sys.path.append(scriptdir) # add default search path
       import log,logging
       
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
       self._AngSamplingRateDeg=str(_AngSamplingRateDeg)
       self._DoRetricSearchbyTiltAngle=_DoRetricSearchbyTiltAngle
       self._Tilt0=_Tilt0
       self._TiltF=_TiltF
       self._ProjMatchingExtra=_ProjMatchingExtra
       self._MaxChangeOffset=str(_MaxChangeOffset)
       self._DisplayMask=_DisplayMask
       #absolute path to main project dir
       self._ProjectDir=_ProjectDir
       self._DoAlign2D=_DoAlign2D
       self._InnerRadius=_InnerRadius
       self._OuterRadius=_OuterRadius
       self._Align2DIterNr=_Align2DIterNr
       self._DisplayAlign2D=_DisplayAlign2D
       self._Align2DExtraCommand=_Align2DExtraCommand
       self._DisplayReconstruction=_DisplayReconstruction
       self._DoReconstruction=_DoReconstruction
       self._DoComputeResolution=_DoComputeResolution
       self._ResolSam=_ResolSam
       self._DoParallel=_DoParallel
       self._MyNumberOfCPUs=_MyNumberOfCPUs
       self._MyMachineFile=os.path.abspath(_MyMachineFile)
       self._Symfile=os.path.abspath(str(_Symfile))
       self._iteration_number=1
       self._ReconstructionExtraCommand=_ReconstructionExtraCommand
       self._ReconstructionMethod=_ReconstructionMethod

       
       self._ReferenceVolume=os.path.abspath(_ReferenceVolume)
       self._Proj_Maching_Output_Root_Name=_Proj_Maching_Output_Root_Name
       self._multi_align2d_sel=_multi_align2d_sel
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
                                      
       if (DoDeleteWorkingDir): 
          delete_working_directory(self._mylog,self._WorkDirectory)
       else:
          self._mylog.info("Skipped DoDeleteWorkingDir") 
       create_working_directory(self._mylog,self._WorkDirectory)
       #made backup of this script
       log.make_backup_of_script_file(sys.argv[0],self._WorkDirectory)
       
       #copy files to local directory
       copy_images_to_local_disk(self._mylog,
                                 self._SelFileName,
                                 self._WorkDirectory)
       #change to working dir
       os.chdir(self._WorkDirectory)
       self._SelFileName=self._WorkDirectory+'/'+\
                         str(os.path.basename(self._SelFileName))
       if (_ResetImageHeader):
           reset_image_header(self._mylog,self._SelFileName)

       ##
       ##LOOP
       ##
       
       #create working dir for this iteration
       Iteration_Working_Directory=self._WorkDirectory+'/Iter_'+\
                                   str(self._iteration_number)
       create_working_directory(self._mylog,Iteration_Working_Directory)
       #change to iteration_working dir
       os.chdir(Iteration_Working_Directory)
       
       #save initial header
       docfile_with_original_angles='docfile_with_original_angles.doc'
       command='xmipp_header_extract -i '
       command=command+self._SelFileName + \
                      ' -o ' + docfile_with_original_angles 
       self._mylog.debug("save original header in file: " + 
                                docfile_with_original_angles ) 
       self._mylog.info(command)
       os.system(command)
       
       if (_DoMask):
          execute_mask(self._mylog,
                       self._ProjectDir,
                       self._ReferenceFileName,
                       self._MaskFileName,
                       self._DisplayMask,
                       self._iteration_number,
                       self._ReferenceVolume)
       else:
          self._mylog.info("Skipped Mask") 
          
       if (_DoProjectionMatching):
          execute_projection_matching(self._mylog,
                                      self._ProjectDir,
                                      self._ReferenceFileName,
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
                                      self._ProjMatchingExtra,
                                      self._DisplayProjectionMatching,
                                      self._iteration_number,
                                      self._DoAlign2D,
                                      self._DoParallel,
                                      self._MyNumberOfCPUs,
                                      self._MyMachineFile,
                                      self._WorkDirectory,
                                      self._Symfile,
                                      self._ReferenceVolume
                                      )
       else:
          self._mylog.info("Skipped ProjectionMatching") 
       if (_DoAlign2D):
          execute_align2d(self._mylog,
                          self._InnerRadius,
                          self._OuterRadius,
                          self._Align2DIterNr,
                          self._Align2DExtraCommand,
                          self._Proj_Maching_Output_Root_Name,
                          self._iteration_number,
                          self._DoParallel,
                          self._MyNumberOfCPUs,
                          self._MyMachineFile,
                          self._DisplayAlign2D,
                          self._multi_align2d_sel
                          )
       else:
          self._mylog.info("Skipped Align2D") 
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
        
       if (_DoReconstruction):
          print self._ReconstructionMetho+"\n"
          execute_reconstruction(self._mylog, 
                                 self._ReconstructionExtraCommand,
                                 self._iteration_number,
                                 self._DisplayReconstruction,
                                 self._DoParallel,
                                 self._MyNumberOfCPUs,
                                 self._MyMachineFile,
                                 self._ReconstructionMethod,
                                 self._align2d_sel,
                                 self._Symfile)
       else:
          self._mylog.info("Skipped Reconstruction") 
       if (_DoComputeResolution):
          execute_resolution(self._mylog, 
                             self._ReconstructionExtraCommand,
                             self._ReconstructionMethod,
                             self._iteration_number,
                             self._DisplayReconstruction,
                             self._ResolSam,
                             self._align2d_sel,
                             self._DoParallel,
                             self._MyNumberOfCPUs,
                             self._MyMachineFile)
       else:
          self._mylog.info("Skipped Reconstruction") 



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
def copy_images_to_local_disk(_mylog,_SelFileName,_WorkDirectory):
      import os,selfile
      print '*********************************************************************'
      print '* Copying images to working directory ...'
      mysel=selfile.selfile()
      mysel.read(_SelFileName)
      newsel=mysel.copy_sel(_WorkDirectory)
      print _WorkDirectory+'/'+_SelFileName
      newsel.write(_WorkDirectory+'/'+_SelFileName)
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
def execute_mask(_mylog,
                 _ProjectDir,
                 _ReferenceFileName,
                 _MaskFileName,
                 _DisplayMask,
                 _iteration_number,
                 _ReferenceVolume):
   import os
   _mylog.debug("execute_mask")
   MaskVolume =_MaskFileName
   MaskedVolume=_ReferenceVolume
   if _iteration_number==1:
      InPutVolume=_ReferenceFileName
   else:
      InPutVolume='../Iter_'+str(_iteration_number-1)+target
         
       
   InPutVolume=_ReferenceFileName
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
      command='xmipp_show -vol '+ MaskedVolume +' -w 10'
      print '*********************************************************************'
      print '* ',command
      _mylog.info(command)
      os.system(command)
   return MaskedVolume
#------------------------------------------------------------------------
#execute_projection_matching
#------------------------------------------------------------------------
def execute_projection_matching(_mylog,
                                _ProjectDir,
                                _ReferenceFileName,
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
                                _ProjMatchingExtra,
                                _DisplayProjectionMatching,
                                _iteration_number,
                                _DoAlign2D,
                                _DoParallel,
                                _MyNumberOfCPUs,
                                _MyMachineFile,
                                _WorkDirectory,
                                _Symfile,
                                _ReferenceVolume):
                                           
   _mylog.debug("execute_projection_matching")
   import os,selfile
   ReferenceVolume=_ReferenceFileName
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
   
   if _DisplayProjectionMatching==True:
      classes_sel_file=selfile.selfile()
      classes_sel_file.read(_Proj_Maching_Output_Root_Name+'_classes.sel')
      library_sel_file=selfile.selfile()
      library_sel_file.read(_Proj_Maching_Output_Root_Name+'_lib.sel')
      newsel=library_sel_file.intercalate_union(classes_sel_file)
      compare_sel_file=_Proj_Maching_Output_Root_Name+'_compare.sel'
      newsel.write(compare_sel_file)
      command='xmipp_show -sel '+ compare_sel_file +' -w 10 &'
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
                    _Proj_Maching_Output_Root_Name,
                    _iteration_number,
                    _DoParallel,
                    _MyNumberOfCPUs,
                    _MyMachineFile,
                    _DisplayAlign2D,
                    _multi_align2d_sel
                    ):
                    
   #if secuential execute orden if not store it in a file
   import tempfile
   tmp_file_name = _Proj_Maching_Output_Root_Name +".tmp"
   if _DoParallel:
        fh = open(tmp_file_name,'w')
      
   _mylog.debug("execute_align2d")
   import os,selfile, glob,shutil
   class_sel_pattern = _Proj_Maching_Output_Root_Name+'_class[0-9]*.sel'
   _mylog.debug("class_sel_pattern: " + class_sel_pattern)
   aux_sel_file=selfile.selfile()
   align2d_sel=selfile.selfile()
   #create compare sel 
   for class_selfile in glob.glob(class_sel_pattern):
      reference=class_selfile.replace('.sel','.xmp')
      aux_sel_file.read(class_selfile)
      lib_file_name = class_selfile.replace('class','lib')
      lib_file_name = lib_file_name.replace('.sel','.proj') 
      class_file_name = class_selfile.replace('.sel','.xmp') 
      #first line in sel must be active
      align2d_sel.insert(lib_file_name,str(1))
      align2d_sel.insert(class_file_name,str(-1))
      aligned_file_name = class_selfile.replace('.sel','.med.xmp')
      align2d_sel.insert(aligned_file_name,str(1))
      if (aux_sel_file.length()<1):
         command="xmipp_operate -i " + \
                  class_file_name + " -mult 0 -o " + aligned_file_name +"\n"
      else:
         print '*********************************************************************'
         print '* Aligning translationally each class'
         command='xmipp_align2d'+ \
                 ' -i '  + class_selfile + \
                 ' -Ri ' +   str(_InnerRadius) + \
                 ' -Ro ' +   str(_OuterRadius) +\
                 ' -iter ' + str(_Align2DIterNr) +\
                 ' -ref ' + reference +\
                 ' '  + _Align2DExtraCommand
      print '* ',command
      if _DoParallel:
         fh.writelines(command+"\n")
         _mylog.debug(command)
      else:  
         os.system(command)
         _mylog.info(command)
            
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
       os.remove(tmp_file_name)
       
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
      command='xmipp_show -showall -sel '+ _multi_align2d_sel +' -w 9'
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
                           _ReconstructionExtraCommand,
                           _iteration_number,
                           _DisplayReconstruction,
                           _DoParallel,
                           _MyNumberOfCPUs,
                           _MyMachineFile,
                           _ReconstructionMethod,
                           _align2d_sel,
                           _Symfile):
   import os,shutil
   _mylog.debug("execute_reconstruction")


   #the user should be able to delete images
   #this is dirty but what else can I do
   reconstruction_sel=_align2d_sel;

   Outputvolume  ="reconstruction_iter_"+str(_iteration_number)+".vol"
   print '*********************************************************************'
   print '* Reconstruct volume using ', _ReconstructionMethod
   if _ReconstructionMethod=='wbp':
      program = 'xmipp_reconstruct_wbp'
      mpi_program = 'xmipp_mpi_reconstruct_wbp'
      parameters= ' -i '    + reconstruction_sel + \
                  ' -o '    + Outputvolume + \
                  ' -weight -use_each_image ' + \
                  ' -sym ' + _Symfile + ' ' + \
                  _ReconstructionExtraCommand
              
   elif _ReconstructionMethod=='art':
      program = 'xmipp_reconstruct_art'
      mpi_program = 'NULL'
      _DoParallel=False
      mpi_program = 'xmipp_mpi_reconstruct_wbp'
      parameters=' -i '    + InPutVolume + \
                 ' -o '    + Outputvolume 
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
                        _ReconstructionExtraCommand,
                        _ReconstructionMethod,
                        _iteration_number,
                        _DisplayReconstruction,
                        _ResolSam,
                        _align2d_sel,
                        _DoParallel,
                        _MyNumberOfCPUs,
                        _MyMachineFile):
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
    Outputvolumes.append(split_sel_root_name+'_1_'+str(_iteration_number)+".vol")
    Outputvolumes.append(split_sel_root_name+'_2_'+str(_iteration_number)+".vol")
    for Outputvolume in Outputvolumes:
       print '*********************************************************************'
       print '* Reconstruct volume'
       if _ReconstructionMethod=='wbp':
          program = 'xmipp_reconstruct_wbp'
          mpi_program = 'xmipp_mpi_reconstruct_wbp'
          parameters= ' -i '    + reconstruction_sel + \
                      ' -o '    + Outputvolume + \
                      ' -weight -use_each_image ' + \
                      ' -sym ' + _Symfile + ' ' + \
                      _ReconstructionExtraCommand

       elif _ReconstructionMethod=='art':
          program = 'xmipp_reconstruct_art'
          mpi_program = 'NULL'
          _DoParallel=False
          mpi_program = 'xmipp_mpi_reconstruct_wbp'
          parameters=' -i '    + InPutVolume + \
                     ' -o '    + Outputvolume 
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



    #print '**************************************************************'
    #print '* Compute resolutione ' + _SelFileName 
    #command = "xmipp_resolution_fsc -ref " + _SelFileName
    #print '* ',command
    #_mylog.info(command)
    #os.system(command)

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
                AngSamplingRateDeg,
                DoRetricSearchbyTiltAngle,
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
                DoComputeResolution,
                ResolSam,
                ReconstructionMethod,
                ReconstructionExtraCommand,
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
                align2d_sel,
                align2d_doc,
                ResetImageHeader
                )
