#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for projection matching
# using self-organizing maps
#
#  - delete and create workdirectory
#  - alignment  (align2D) include statis BEFORE alignment
#  - statis 
#        - applygeo to med
#        - reverse_endian (if needed apply reverse endian to average image)
#  - reverse_endian (if needed apply reverse endian to apply geo images)
#  - findcenter
#  - make spectra
#  - kerdenSOM based classification
#
# Example use:
# ./xmipp_protocol_projmatch.py
#
# Author:Roberto Marabini, March 2007
#
# required files: log.py, s, Sel_Files.py
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# Selfile with the input images:
SelFileName='../all.sel'

# Working directory: 
WorkDirectory='test1'

# Delete working directory if it already exists?
DoDeleteWorkingDir=False

# {expert} Root directory name for this project:
ProjectDir="/home/roberto2/Test/PARA_Roberto"

# {expert} Directory name for logfiles:
LogDir="Logs"

# Reference file name (3D map)
ReferenceFileName="init_reference/LTA_rot_0.1_norm.vol"

#binary mask-file used to mask reference volume
MaskFileName='circular_mask.msk'

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel? (see Expert options)
DoParallel=False
# Number of processors to use:
MyNumberOfCPUs=10
# {expert} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
"""
MyMachineFile="/home2/bioinfo/scheres/machines.dat"
# {expert} Standard script to launch MPI-jobs:
""" This will also depend on your system...
    The simplest script consists of the following two lines:

    #!/usr/bin/env sh
    `which mpirun` -np MyNumberOfCPUs -machinefile MyMachineFile ~/machines.dat \

    Note that the next line with the xmipp_mpi_MLalign2D command will be
    generated automatically, and the variables MyNumberOfCPUs and MyMachineFile
    will be replaced by the corresponding values given here above.
    More scripts for different batch systems can be found at:
    [Wiki link]
"""

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
DisplayMask=False

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

#show masked volume
"""show average of projections.
"""

#{expert} Root name for output projections
OutputRootName='root'

#Angular sampling rate (in degrees)
AngSamplingRateDeg=20

#Maximum change in origin offset (+/- pixels)
MaxChangeOffset=4

#extra options for Projection_Maching
ProjMatchingExtra=""

#stopping rule
#projection_matching.1
#-----------------------------------------------------------------------------
# {end-of-header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
class projection_matching_class:

   #init variables
   
   def __init__(self,
                _DoMask, _DisplayMask,
                _ReferenceFileName,
                _MaskFileName,
                _DoProjectionMatching, _DisplayProjectionMatching,
                _OutputRootName,_AngSamplingRateDeg,_ProjMatchingExtra,
                _MaxChangeOffset,
                _SelFileName,
                _WorkDirectory,
                _ProjectDir,
                _LogDir
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
       self._ProjMatchingExtra=_ProjMatchingExtra
       self._MaxChangeOffset=str(_MaxChangeOffset)
       self._DisplayMask=_DisplayMask
       self._ProjectDir=_ProjectDir
       
       
       
       self.mylog=log.init_log_system(_ProjectDir,
                                      _LogDir,
                                      sys.argv[0],
                                      _WorkDirectory)
                                      
       #uncomment next line to get Degub level logging
       self.mylog.setLevel(logging.DEBUG)
       self.mylog.debug("Debug level logging enabled")
                                      
       if (DoDeleteWorkingDir): 
          self.delete_working_directory()
       self.create_working_directory()

       #change to working dir
       os.chdir(self._WorkDirectory)
       
       self.iteration_number=0
       if (DoMask):
          self.execute_mask()
       if (DoProjectionMatching):
          self.execute_projection_matching()

   #------------------------------------------------------------------------
   #delete_working directory
   #------------------------------------------------------------------------
   def delete_working_directory(self):
       import os
       import shutil
       print '*********************************************************************'
       print '* Delete working directory tree'
       self.mylog.info("Delete working directory tree")
       
       if os.path.exists(self._WorkDirectory):
          shutil.rmtree(self._WorkDirectory)
   #------------------------------------------------------------------------
   #create_working directory
   #------------------------------------------------------------------------
   def create_working_directory(self):
       import os
       print '*********************************************************************'
       print '* Create working directory'
       self.mylog.info("Create working directory")
       
       if not os.path.exists(self._WorkDirectory):
          os.makedirs(self._WorkDirectory)
        
   #------------------------------------------------------------------------
   #execute_mask
   #------------------------------------------------------------------------
   def execute_mask(self):
      import os
      self.mylog.debug("execute_mask")
      InPutVolume=self._ProjectDir+'/'+self._ReferenceFileName
      MaskVolume =self._ProjectDir+'/'+self._MaskFileName
      tmp_OutPutVol=os.path.basename(self._ReferenceFileName)
      ReferenceFileName_without_ext=\
            (os.path.splitext(str(os.path.basename(tmp_OutPutVol))))[0]
      self.Masked_Vol=ReferenceFileName_without_ext+'.msk'
      print '*********************************************************************'
      print '* Mask the reference volume'
      command='xmipp_mask'+ \
              ' -i '    + InPutVolume + \
              ' -o '    + self.Masked_Vol + \
              ' -mask ' + MaskVolume 
              
      print '* ',command
      self.mylog.info(command)
      os.system(command)
      if self._DisplayMask==True:
         command='xmipp_show -vol '+ OutPutVol +' &'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
         os.system(command)

   #------------------------------------------------------------------------
   #execute_projection_matching
   #------------------------------------------------------------------------
   def execute_projection_matching(self):
      self.mylog.debug("execute_projection_matching")
      import os,SelFiles
      ReferenceVolume=self._ProjectDir+'/'+self._ReferenceFileName

      print '*********************************************************************'
      print '* Asign projection direction'
      command='xmipp_projection_matching'                  + \
              ' -i '           + self._SelFileName         + \
              ' -vol '         + self.Masked_Vol           + \
              ' -o '           + self._OutputRootName      + \
              ' -sam '         + self._AngSamplingRateDeg  + \
              ' -max_shift '   + self._MaxChangeOffset     + \
              ' -dont_modify_header -output_refs -output_classes ' + \
              ' '              + self._ProjMatchingExtra

      print '* ',command
      self.mylog.info(command)
      ######################os.system(command)
      if self._DisplayProjectionMatching==True:
         classes_sel_file=SelFiles.selfile()
         classes_sel_file.read(self._OutputRootName+'_classes.sel')
         library_sel_file=SelFiles.selfile()
         library_sel_file.read(self._OutputRootName+'_lib.sel')
         newsel=library_sel_file.intercalate_union(classes_sel_file)
         compare_sel_file=self._OutputRootName+'_compare.sel'
         newsel.write(compare_sel_file)
         command='xmipp_show -sel '+ compare_sel_file +' &'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command) 
         os.system(command)

#
# main
#     
if __name__ == '__main__':

    # create rotational_spectra_class object
    # 
   
    #init variables
  my_projmatch=projection_matching_class(
                DoMask, DisplayMask,
                ReferenceFileName,
                MaskFileName,
                DoProjectionMatching, DisplayProjectionMatching,
                OutputRootName,AngSamplingRateDeg,ProjMatchingExtra,
                MaxChangeOffset,
                SelFileName,
                WorkDirectory,
                ProjectDir,
                LogDir
                )
