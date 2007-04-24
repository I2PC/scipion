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
# Selfile with the input images (relative path from ProjectDir):
""" The name of this file will be used as a seed for all new files
"""
SelFileName='100.sel'
# Working subdirectory: 
WorkDirectory='test1'
# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=True
# Display all intermediate results?
""" Images will be displayed after alignment, average, makespectra and kendersom...
"""
DisplayResults=True
# {expert} Root directory name for this project:
ProjectDir="/home/roberto2/Test/Para_Roberto/"
# {expert} Directory name for logfiles:
LogDir="Logs"
#-----------------------------------------------------------------------------
# {section} Align2d
#-----------------------------------------------------------------------------
# Perform 2D pre-alignment?
DoAlign2D=True
# Inner radius for rotational correlation (also used to make spectra):
""" These values are in pixels from the image center
"""
InnerRadius=3
# Outer radius for rotational correlation (also used to make spectra):
OuterRadius=12
# Number of align2d iterations to perform:
Align2DIterNr=4
# {expert} Additional align2d parameters
""" For a complete description, see http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d

  Examples:
  Align2DExtraCommand=\"-trans_only  -filter 10 -sampling 2 -max_shift 2 -max_rot 3\"
  Align2DExtraCommand=\"-max_shift 2 -max_rot 3\"

  consider filtering the images with \"-filter 10 -sampling 2\"
"""
Align2DExtraCommand=""
#-----------------------------------------------------------------------------
# {section} Find_center
#-----------------------------------------------------------------------------
# Perform search for symmetry center?
DoFindCenter=True
#-----------------------------------------------------------------------------
# {section} Make_spectra
#-----------------------------------------------------------------------------
# Perform Rotational Spectra calculation?
DoMakeSpectra=True
# Name of the output file:
""" Existing files with this name will be deleted!
"""
SpectraName="spectra.sim"
#-----------------------------------------------------------------------------
# {section} classify_kerdensom 
#-----------------------------------------------------------------------------
# Perform kerdensom calculation?
DoKerdensom=True
# Name of Output SOM:
""" Existing files with this name will be deleted!
"""
SomName="som"
# X-dimension of the self-organizing map:
SomXdim=7
# Y-dimension of the self-organizing map:
SomYdim=7
# Initial regularization factor:
""" The kerdenSOM algorithm anneals from an initial high regularization factor
    to a final lower one, in a user-defined number of steps.
    If the output map is too smooth, lower the regularization factors
    If the output map is not organized, higher the regularization factors
"""
SomReg0=1000
# Final regularization factor:
SomReg1=200
# Number of steps to lower the regularization factor:
SomSteps=5
# {expert} Additional kerdenSOM parameters:
""" For a complete description see http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraCommand=""
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
                _InnerRadius,
                _OuterRadius,
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
       scriptdir=os.path.expanduser('~')+'/trunk/Xmipp/Applications/Batches/Protocols_lib/'
       sys.path.append(scriptdir) # add default search path
       import log

       self._WorkDirectory=os.getcwd()+'/'+_WorkDirectory
       self._SelFileName=_ProjectDir+'/'+str(_SelFileName)
       self._DisplayResults=_DisplayResults
       self._InnerRadius=_InnerRadius
       self._OuterRadius=_OuterRadius
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

       self.mylog=log.init_log_system(_ProjectDir,
                                      _LogDir,
                                      sys.argv[0],
                                      _WorkDirectory)

       if (_DoDeleteWorkingDir): 
          self.delete_working_directory()
       self.create_working_directory()

       #change to working dir
       os.chdir(self._WorkDirectory)

       if (_DoAlign2D):
          self.execute_align2d()

       if (_DoFindCenter):
          self.execute_findcenter()
       
       if (_DoMakeSpectra):
           if(self.true_if_file_is_NOT_in_native_endian()):
              self.execute_reverse_endian()
           self.execute_apply_geo()
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
   #execute_align2d
   #------------------------------------------------------------------------
   def execute_align2d(self):
      import os,selfile
      print '*********************************************************************'
      print '* Copying images to working directory ...'
      mysel=selfile.selfile()
      mysel.read(self._SelFileName)
      newsel=mysel.copy_sel('.')
      newsel.write(os.path.basename(self._SelFileName))

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
              ' -iter ' + str(self._Align2DIterNr) +\
              ' -ref ' + selfile_without_ext + '.med.xmp' +\
              ' '  + self._Align2DExtraCommand
      print '* ',command
      self.mylog.info(command)
      os.system(command)
      if self._DisplayResults==True:
         command='xmipp_show -sel '+ os.path.basename(self._SelFileName)+' &'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info(command)
         os.system(command)
      if self._DisplayResults==True:
         command='xmipp_show -img '+ selfile_without_ext + '.med.xmp &'
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
              ' -r1 '  + str(self._InnerRadius) +   ' -r2 '   + str(self._OuterRadius) +\
              ' -low ' + str(self._OuterRadius+2) + ' -high ' + str(self._OuterRadius+5)
              
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
              ' -r1 '  + str(self._InnerRadius) +   ' -r2 '   + str(self._OuterRadius) 
              
      print '* ',command
      self.mylog.info(command)
      program = os.system(command)     
      if(self._DisplayResults==True):
          print '*********************************************************************'
          print '* Display spectra'
          selFileName=os.path.basename(self._SelFileName)
          command='xmipp_show'+ \
                  ' -spect ' + str(self._SpectraName)+' &'
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
      command='xmipp_classify_kerdensom'+ ' -verb 1 -din '  + str(self._SpectraName) + \
              ' -cout ' + str(self._SomName)  + \
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
                  ' -din ' + str(self._SpectraName) + ' &'
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
                                      InnerRadius,
                                      OuterRadius,
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

    
