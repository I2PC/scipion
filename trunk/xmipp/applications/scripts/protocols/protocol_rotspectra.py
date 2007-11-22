#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for rotational spectra classification
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
# ./xmipp_rotational_spectra.py
#
# Author:Roberto Marabini, March 2007
#
# required files: log.py, spider_header.py, Sel_Files.py
# by defaults searchs for python files in ../python directory
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# {file} Selfile with the input images:
SelFileName='all_images.sel'
# Working subdirectory: 
WorkDirectory='RotSpectra/test1'
# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=False
# Display all intermediate results?
""" Images will be displayed after alignment, average, makespectra and kendersom...
    The display is NOT make in the background. If you want to 
    send to a queue this job set this flag to NULL 
"""
DisplayResults=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
LogDir='Logs'
#-----------------------------------------------------------------------------
# {section} 2D Pre-alignment
#-----------------------------------------------------------------------------
# Perform 2D pre-alignment?
DoAlign2D=True
# Use Selfile average as reference image
""" Is this option is set to False then the reference will be computed as a  
    piramidal combination of subset of images. Details are available at 
    http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d
"""
DoAverageAsReference=True
# Inner radius for rotational correlation:
""" These values are in pixels from the image center
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d
"""
InnerRadius=3
# Outer radius for rotational correlation (also used to make spectra):
OuterRadius=12
# Number of align2d iterations to perform:
Align2DIterNr=4
# {expert} Additional align2d parameters
""" For a complete description, 
    see http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d

  Examples:
  Align2DExtraCommand="-trans_only  -filter 10 -sampling 2 -max_shift 2 -max_rot 3"
  Align2DExtraCommand="-max_shift 2 -max_rot 3"

  consider filtering the images with "-filter 10 -sampling 2"
"""
Align2DExtraCommand=''
#-----------------------------------------------------------------------------
# {section} Find the optimal center of symmetry
#-----------------------------------------------------------------------------
# Perform search for symmetry center?
""" A search will be performed for the center of symmetry by minimizing the first harmonic
"""
DoFindCenter=True
#-----------------------------------------------------------------------------
# {section} Rotational spectra calculation
#-----------------------------------------------------------------------------
# Calculate the rotational spectra?
DoMakeSpectra=True
# Inner radius for rotational harmonics calculation:
""" These values are in pixels from the image center
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Makespectra
"""
SpectraInnerRadius=7
# Outer radius for rotational harmonics calculation:
SpectraOuterRadius=10
# Name of the output file:
""" Existing files with this name will be deleted!
"""
SpectraName='spectra.sim'
# {expert} Lowest harmonic to calculate
SpectraLowHarmonic=1
# {expert} Highest harmonic to calculate
SpectraHighHarmonic=15
#-----------------------------------------------------------------------------
# {section} Classification: classify_kerdensom 
#-----------------------------------------------------------------------------
# Perform self-organizing map calculation?
DoKerdensom=True
# Name of Output SOM:
""" Existing files with this name will be deleted!
"""
SomName='som'
# X-dimension of the self-organizing map:
SomXdim=7
# Y-dimension of the self-organizing map:
SomYdim=7
# Initial regularization factor:
""" The kerdenSOM algorithm anneals from an initial high regularization factor
    to a final lower one, in a user-defined number of steps.
    If the output map is too smooth, lower the regularization factors
    If the output map is not organized, higher the regularization factors
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
SomReg0=1000
# Final regularization factor:
SomReg1=200
# Number of steps to lower the regularization factor:
SomSteps=5
# {expert} Additional kerdenSOM parameters:
""" For a complete description 
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraCommand=''
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_rotspectra.py'
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
class rotational_spectra_class:

   #init variables
   
   def __init__(self,
                _SelFileName,
                _WorkDirectory,
                _DoDeleteWorkingDir,
                _DisplayResults,
                _ProjectDir,
                _LogDir,
                _DoAlign2D,
                _DoAverageAsReference,
                _InnerRadius,
                _OuterRadius,
                _SpectraInnerRadius,
                _SpectraOuterRadius,
                _SpectraLowHarmonic,
                _SpectraHighHarmonic,
                _Align2DIterNr,
                _Align2DExtraCommand,
                _DoFindCenter,
                _DoMakeSpectra,
                _SpectraName,
                _DoKerdensom,
                _SomName,
                _SomXdim,
                _SomYdim,
                _SomReg0,
                _SomReg1,
                _SomSteps,                
                _KerdensomExtraCommand,
                ):

       import os,sys
       scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
       sys.path.append(scriptdir) # add default search path
       import log

       self._WorkDirectory=os.getcwd()+'/'+_WorkDirectory
       self._SelFileName=_SelFileName
       #os.path.abspath(_SelFileName)
       self._DisplayResults=_DisplayResults
       self._InnerRadius=_InnerRadius
       self._OuterRadius=_OuterRadius
       self._SpectraInnerRadius=_SpectraInnerRadius
       self._SpectraOuterRadius=_SpectraOuterRadius
       self._SpectraLowHarmonic=_SpectraLowHarmonic
       self._SpectraHighHarmonic=_SpectraHighHarmonic
       self._Align2DIterNr=_Align2DIterNr
       self._Align2DExtraCommand=_Align2DExtraCommand
       self._Reverseendian=1
       self._SpectraName=_SpectraName
       self._SomName=_SomName
       self._SomXdim=_SomXdim
       self._SomYdim=_SomYdim
       self._SomReg0=_SomReg0
       self._SomReg1=_SomReg1
       self._SomSteps=_SomSteps              
       self._KerdensomExtraCommand=_KerdensomExtraCommand
       self._ProjectDir=_ProjectDir
       self._DoAverageAsReference=_DoAverageAsReference
       self.mylog=log.init_log_system(_ProjectDir,
                                      _LogDir,
                                      sys.argv[0],
                                      _WorkDirectory)
       if (_DoDeleteWorkingDir): 
          self.delete_working_directory()
       self.create_working_directory()
       
       #made backup of this script
       log.make_backup_of_script_file(sys.argv[0],self._WorkDirectory)


       if (_DoAlign2D):
          self.copy_images_to_working_dir()
       #change to working dir
       os.chdir(self._WorkDirectory)
       if (_DoAlign2D):
          self.execute_align2d()

       if (_DoFindCenter):
          self.execute_findcenter()
       
       if (_DoMakeSpectra):
           //makespectra now applies geo and understands both endians
           #if(self.true_if_file_is_NOT_in_native_endian()):
           #   self.execute_reverse_endian()
           #   self.execute_apply_geo()
           self.execute_spectra()
       if (_DoKerdensom):
           self.delete_existing_som()
           self.execute_KerDenSOM()
       self.close()
  
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
   #copy_images_to_working_dir
   #------------------------------------------------------------------------
   def copy_images_to_working_dir(self):
      import os,selfile
      print '*********************************************************************'
      print '* Copying images to working directory ...'
      mysel=selfile.selfile()
      mysel.read(self._ProjectDir+'/'+self._SelFileName)
      newsel=mysel.copy_sel(self._WorkDirectory)
      newsel.write(self._WorkDirectory+'/'+self._SelFileName)

   #------------------------------------------------------------------------
   #execute_align2d
   #------------------------------------------------------------------------
   def execute_align2d(self):
      import os,selfile
      if self._DoAverageAsReference:
         print '*********************************************************************'
         print '* Computing initial reference using average'
         command = 'xmipp_average -i ' + os.path.basename(self._SelFileName)
         print '* ',command
         os.system(command)
         self.mylog.info(command)
      
      selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
      print '*********************************************************************'
      print '* Aligning translationally and rotationally a set of 2D-images'
      command='xmipp_align2d'+ \
              ' -i '  + os.path.basename(self._SelFileName) + \
              ' -Ri ' + str(self._InnerRadius) + \
              ' -Ro ' + str(self._OuterRadius) +\
              ' -iter ' + str(self._Align2DIterNr)
      if self._DoAverageAsReference:
              command = command + ' -ref ' + selfile_without_ext + '.med.xmp'
      command = command +' '  + self._Align2DExtraCommand
      print '* ',command
      self.mylog.info(command)
      os.system(command)
      if self._DisplayResults==True:
         command='xmipp_show -sel '+ os.path.basename(self._SelFileName)
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
         os.system(command)
      if self._DisplayResults==True:
         command='xmipp_show -img '+ selfile_without_ext + '.med.xmp'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
         os.system(command)
      print "\n"

   #------------------------------------------------------------------------
   #execute_findcenter
   #------------------------------------------------------------------------
   def execute_findcenter(self):
      import os, string,spider_header
      selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
      print '*********************************************************************'
      print '* Looking for the position of the center of symmetry ...'
      filename=selfile_without_ext + '.med.xmp'
      myheader=spider_header.spiderheader(filename)
      ncolumns=myheader.nx
      nrows=myheader.ny
      nslices=myheader.nz
      command='xmipp_find_center2d'+ \
              ' -img ' + filename + \
              ' -x0 '  + str((ncolumns-1)/2) + \
              ' -y0 '  + str((nrows   -1)/2) + \
              ' -r1 '  + str(self._SpectraInnerRadius) +   ' -r2 '   + str(self._SpectraOuterRadius) +\
              ' -low ' + str(self._SpectraOuterRadius+2) + ' -high ' + str(self._SpectraOuterRadius+5)
              
      print '* ',command
      self.mylog.info(command)
      program = os.popen(command,"r")
      data = program.read()
      print data
      data=data.split()
      aux     = data[12].split(',')[0]
      self.xOffset = float (aux)
      self.yOffset = float(data[15])

   #------------------------------------------------------------------------
   #execute_apply_geo
   #------------------------------------------------------------------------
   def execute_apply_geo(self):
      import os
      print '*********************************************************************'
      print '* Applying geometrical  information in the headers. Needed for makespectra'
      command='xmipp_header_apply'+ \
              ' -i ' +os.path.basename(self._SelFileName) 
      print '* ',command
      self.mylog.info(command)
      program = os.system(command)     

   #------------------------------------------------------------------------
   #true_if_file_is_NOT_in_native_endian
   #------------------------------------------------------------------------
   def true_if_file_is_NOT_in_native_endian(self):
      import selfile,spider_header,os
      print '*********************************************************************'
      print '*  check if images are in native endian (spectra preprocesing)'
      mysel=selfile.selfile()
      #mysel=selfile.selfile()
      mysel.read(os.path.basename(self._SelFileName))
      filename,state=mysel.find_first_active_image()
      myheader=spider_header.spiderheader(filename)
      if( myheader.check_endianess()):
         print '* No endian change needed.'
      else:
         print '* Endian change needed.'
         
      if myheader.check_endianess(): return False
      else: return True

   #------------------------------------------------------------------------
   #execute_reverse_endian
   #------------------------------------------------------------------------
   def execute_reverse_endian(self):
      import os
      print '*********************************************************************'
      print '* Changing endian format. Needed for makespectra'
      command='xmipp_reverse_endian'+ \
              ' -i ' + os.path.basename(self._SelFileName) 
      print '* ',command
      self.mylog.info(command)
      program = os.system(command)     

   #------------------------------------------------------------------------
   #execute_spectra
   #------------------------------------------------------------------------
   def execute_spectra(self):
      import os, string
      if (os.path.exists(self._SpectraName)):
          os.remove(self._SpectraName)
      print '*********************************************************************'
      print '* Computing rotational power spectra'
      selFileName=os.path.basename(self._SelFileName)
      outFileName=(os.path.splitext(selFileName))[0] + '.sim'
      command='xmipp_make_spectra'+ \
              ' -i ' + selFileName + \
              ' -o ' + str(self._SpectraName) + \
              ' -x0 '  + str(self.xOffset) + \
              ' -y0 '  + str(self.yOffset) + \
              ' -r1 '  + str(self._SpectraInnerRadius) + \
              ' -r2 '   + str(self._SpectraOuterRadius) + \
              ' -low ' + str(self._SpectraLowHarmonic) + \
              ' -high ' + str(self._SpectraHighHarmonic)
              
      print '* ',command
      self.mylog.info(command)
      program = os.system(command)     
      if(self._DisplayResults==True):
          print '*********************************************************************'
          print '* Display spectra'
          selFileName=os.path.basename(self._SelFileName)
          command='xmipp_show'+ \
                  ' -spect ' + str(self._SpectraName)
          print '* ',command
          self.mylog.info(command)
          program = os.system(command)     

  
   #------------------------------------------------------------------------
   # f a SOM exists with SomName, delete all related files
   #------------------------------------------------------------------------
   def delete_existing_som(self):
       import os
       # delete existing files with this somname
       if (os.path.exists(self._SomName+'.sel')):
           print 'Deleting existing som...'
           command= 'xmipp_selfile_delete '+self._SomName+'.sel \n'
           print '* ',command
           self.log.info(command)
           os.system(command)
           command= 'rm -f '+self._SomName+'.* '+self._SomName+'_* \n'
           print '* ',command
           self.log.info(command)
           os.system(command)

   #------------------------------------------------------------------------
   #execute_KerDenSOM
   #------------------------------------------------------------------------
   def execute_KerDenSOM(self):
      import os

      selFileName=os.path.basename(self._SelFileName)
      print '*********************************************************************'
      print '* Computing kerdensom ...'
      command='xmipp_classify_kerdensom'+ ' -verb 1 -i '  + str(self._SpectraName) + \
              ' -o '    + str(self._SomName)  + \
              ' -xdim ' + str(self._SomXdim) + \
              ' -ydim ' + str(self._SomYdim) + \
              ' -reg0 ' + str(self._SomReg0) + \
              ' -reg1 ' + str(self._SomReg1) + \
              ' -steps ' + str(self._SomSteps) + \
              ' '  + str(self._KerdensomExtraCommand)
      print '* ',command
      self.mylog.info(command)
      os.system(command)
      if self._DisplayResults==True:
         command='xmipp_show -spectsom ' + \
                  str(self._SomName) + \
                  ' -din ' + str(self._SpectraName)
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
         os.system(command)
      print "\n"

   #------------------------------------------------------------------------
   #clean and close everything
   #------------------------------------------------------------------------
   def close(self):
         message="* Script finished!"
         print '*********************************************************************'
         print '* ',message
         self.mylog.info(message)
      

#
# main
#     
if __name__ == '__main__':

    # create rotational_spectra_class object
    # 
   
   #init variables
   
   myspectra=rotational_spectra_class(SelFileName,
                                      WorkDirectory,
                                      DoDeleteWorkingDir,
                                      DisplayResults,
                                      ProjectDir,
                                      LogDir,
                                      DoAlign2D,
                                      DoAverageAsReference,
                                      InnerRadius,
                                      OuterRadius,
                                      SpectraInnerRadius,
                                      SpectraOuterRadius,
                                      SpectraLowHarmonic,
                                      SpectraHighHarmonic,
                                      Align2DIterNr,
                                      Align2DExtraCommand,
                                      DoFindCenter,
                                      DoMakeSpectra,
                                      SpectraName,
                                      DoKerdensom,
                                      SomName,
                                      SomXdim,
                                      SomYdim,
                                      SomReg0,
                                      SomReg1,
                                      SomSteps,                
                                      KerdensomExtraCommand,
                                      )

    
