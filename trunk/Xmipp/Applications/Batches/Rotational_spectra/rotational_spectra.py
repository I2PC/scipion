#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for rotational spectra classification: 
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
#
# Example use:
# ./xmipp_rotational_spectra.py
#
# Author:Roberto Marabini, March 2007
#
# required files: log.py, spider_header.py, Sel_Files.py
# by defaults searchs for python files in ../python directory
#
#

#-----------------------------------------------------------------------------
# GLOBAL VARIABLES
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#Path to the Selfile name with input images
#will be used as seed name for the
#new files
#-----------------------------------------------------------------------------
_SelFileName='../untSelect.sel'

#-----------------------------------------------------------------------------
#Work directory name
#all images will be created here except for Logs
#this directory is deleted at the begin of the script
#-----------------------------------------------------------------------------
_WorkDirectory='test1'

#-----------------------------------------------------------------------------
# If set to true images will be
# displayed after aligment, average, makespectra and kendersom... 
# valid values are True/False (note first capital letter)
#-----------------------------------------------------------------------------
_DisplayResults=True

#-----------------------------------------------------------------------------
# Only region between radii InnerRadius and OuterRadius [in pixels]
# will be considered for rotational correlation (alig2d, findcenter)
#-----------------------------------------------------------------------------
_InnerRadius=3
_OuterRadius=12


#-----------------------------------------------------------------------------
# delete working directory
#-----------------------------------------------------------------------------
# Set to True (recomended) to perform delete the working directory,
# False for skipping 
DO_DELETE_WORKING_DIRECTORY=True

#-----------------------------------------------------------------------------
# align2d
#-----------------------------------------------------------------------------
# Set to True to perform 2D alignment, False for skipping 
DO_ALIGN2D=True

#-----------------------------------------------------------------------------
# If you want to add any extra commands to align2d type them here
# examples 
# _Align2DExtraCommand="-iter 3 -trans_only  -filter 10 -sampling 2 -max_shift 2 -max_rot 3
# _Align2DExtraCommand="-iter 4 -max_shift 2 -max_rot 3"
#
# consider filtering the images with "-filter 10 -sampling 2"
# see manual at http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Align2d
_Align2DExtraCommand="-iter 4 -max_shift 2 -max_rot 3"

#-----------------------------------------------------------------------------
# findcenter
#-----------------------------------------------------------------------------
# Set to True to lok for the symmetry center, False for skipping 
DO_FINDCENTER=True

#-----------------------------------------------------------------------------
# makespectra
#-----------------------------------------------------------------------------
# Set to True to compute rotationla spectra, False for skipping 
DO_MAKESPECTRA=True

#-----------------------------------------------------------------------------
# kerdensom 
#-----------------------------------------------------------------------------
DO_KERDENSOM=True

#-----------------------------------------------------------------------------
# kerdensom has many parameter, consider playing with 
# -reg0 500 -reg1 100 -steps 5
# and
# -xdim 7 -ydim 7
# 
# -xdim: indicates the horizontal dimension of the output map
# -ydim: indicates the vertical dimension of the output map
# -reg0, reg1, steps: control the annealing of the regularization
# (smoothness) factor (from -reg0 to -reg1 in -steps). 
# If the output map looks too smooth decrease -reg0, -reg1
# If the output map looks not organized increase -reg0, -reg1
# see manual at http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
#-----------------------------------------------------------------------------
_KerdensonExtraCommand="-xdim 7 -ydim 7 -reg0 500 -reg1 100 -steps 5"

#-----------------------------------------------------------------------------
# Set to True to compute kenderson, False for skipping 

#-----------------------------------------------------------------------------
# SOM

#soy log file name

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
class rotational_spectra_class:

   #init variables
   
   def __init__(self,WorkDirectory,\
                     SelFileName,DisplayResults,
                     InnerRadius,OuterRadius,
                     Align2DExtraCommand,
                     KendersonExtraCommand
                     ):
       import os
       import sys
       sys.path.append(os.getcwd()+'/'+'../Python')#add default search path
       #import log
       import logging
       self._Base_WorkDirectory=WorkDirectory
       self._WorkDirectory=os.getcwd()+'/'+WorkDirectory
       self._SelFileName=os.path.abspath(str(SelFileName))
       self._DisplayResults=DisplayResults
       self._InnerRadius=InnerRadius
       self._OuterRadius=OuterRadius
       self._Align2DExtraCommand=Align2DExtraCommand
       self._Reverseendian=1
       self._KerdensonExtraCommand=KendersonExtraCommand
       self.myBase_FileName = (os.path.splitext(str(SelFileName)))[0] 
       self.myAbsolute_Root_Path= os.path.dirname(self._SelFileName)+'/'

       self.init_log_system()
       if (DO_DELETE_WORKING_DIRECTORY): 
          self.delete_working_directory()
       self.create_working_directory()

       #change to working dir
       os.chdir(self._WorkDirectory)

       if (DO_ALIGN2D): 
          self.execute_align2d()

       if (DO_FINDCENTER): 
          self.execute_findcenter()
       
       if (DO_MAKESPECTRA):
           if(self.true_if_file_is_NOT_in_native_endian()):
              self.execute_reverse_endian()
           self.execute_apply_geo()
           self.execute_spectra()
       if (DO_KERDENSOM):
          self.execute_KerDenSOM()
       self.close()
  
   #------------------------------------------------------------------------
   #create log file and initialize log system
   #------------------------------------------------------------------------
   def init_log_system(self):
       import os, sys
       import logging
       import socket
       """ Set up logging information
           more info at:
          http://webmaster.iu.edu/tool_guide_info/python/lib/module-logging.html
       """
       LogName = self.myAbsolute_Root_Path + 'Logs/'
       if ( not os.path.exists(LogName ) ):
                os.makedirs(LogName)
       LogName += (os.path.splitext(str(sys.argv[0])))[0]
       LogName += '_'
       LogName += self._Base_WorkDirectory
       LogName += '.log'

       self.mylog = logging.getLogger(str(sys.argv[0]))
       hdlr = logging.FileHandler(LogName)
       formatter = logging.Formatter('%(levelname)s (%(lineno)d) %(message)s (%(asctime)s)')
       hdlr.setFormatter(formatter)
       self.mylog.addHandler(hdlr) 
       self.mylog.setLevel(logging.INFO)

       # append a line with user, machine and date
       myusername = str(os.environ.get('USERNAME'))
       myhost = str(socket.gethostname())
       mypwd = str(os.environ.get('PWD'))
       event = "\n"
       event += "NEW LOG SESSION\n" 
       event += "===============\n" 
       event += myusername + '@' 
       event += myhost + ':' 
       event += mypwd 
       self.mylog.info(event)

   #------------------------------------------------------------------------
   #delete_working directory
   #------------------------------------------------------------------------
   def delete_working_directory(self):
       import os
       import shutil
       print '*********************************************************************'
       print '* Delete working directory tree'
       self.mylog.info("command:\tDelete working directory tree")
       
       if os.path.exists(self._WorkDirectory):
          shutil.rmtree(self._WorkDirectory)
   #------------------------------------------------------------------------
   #create_working directory
   #------------------------------------------------------------------------
   def create_working_directory(self):
       import os
       print '*********************************************************************'
       print '* Create working directory'
       self.mylog.info("command:\tCreate working directory")
       
       if not os.path.exists(_WorkDirectory):
          os.makedirs(_WorkDirectory)
        
   #------------------------------------------------------------------------
   #execute_align2d
   #------------------------------------------------------------------------
   def execute_align2d(self):
      import os
      print '*********************************************************************'
      print '*  cp images to working directory'
      command = 'xmipp_cpsel ' + self._SelFileName + ' ./'
      print '* ',command
      os.system(command)
      self.mylog.info("command:\t"+command)
      print '*********************************************************************'
      print '*  create selfile'
      command = 'rm *sel; xmipp_do_selfile  \"*\" | grep -v \".sel \"> ' + os.path.basename(self._SelFileName)
      print '* ',command
      os.system(command)
      self.mylog.info("command:\t"+command)
      print '*********************************************************************'
      print '*  computing initial reference using statis'
      command = 'xmipp_statis -i ' + os.path.basename(self._SelFileName)
      print '* ',command
      os.system(command)
      self.mylog.info("command:\t"+command)
      
      selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
      print '*********************************************************************'
      print '*  align translationally and rotationally a set of 2D-images'
      command='xmipp_align2D'+ \
              ' -i '  + os.path.basename(self._SelFileName) + \
              ' -Ri ' + str(self._InnerRadius) + \
              ' -Ro ' + str(self._OuterRadius) +\
              ' -ref ' + selfile_without_ext + '.med.xmp' +\
              ' '  + self._Align2DExtraCommand
      print '* ',command
      self.mylog.info("command:\t"+command)
      os.system(command)
      if self._DisplayResults==True:
         command='xmipp_show -sel '+ os.path.basename(self._SelFileName)
         print '*********************************************************************'
         print '* ',command
         print '* You may consider removing bad images'
         print ' double click them and sale sel file as DISCARTED'
         self.mylog.info("command:\t"+command)
         os.system(command)
      if self._DisplayResults==True:
         command='xmipp_show -img '+ selfile_without_ext + '.med.xmp'
         print '*********************************************************************'
         print '* ',command
         self.mylog.info("command:\t"+command)
         os.system(command)
      print "\n"

   #------------------------------------------------------------------------
   #execute_findcenter
   #------------------------------------------------------------------------
   def execute_findcenter(self):
      import os, string,spider_header
      selfile_without_ext=(os.path.splitext(str(os.path.basename(self._SelFileName))))[0]
      print '*********************************************************************'
      print '*  Look for the position of the center of symmetry of an image.'
      filename=selfile_without_ext + '.med.xmp'
      myheader=spider_header.spiderheader(filename)
      ncolumns=myheader.nx
      nrows=myheader.ny
      nslices=myheader.nz
      command='xmipp_findcenter'+ \
              ' -img ' + filename + \
              ' -x0 '  + str((ncolumns-1)/2) + \
              ' -y0 '  + str((nrows   -1)/2) + \
              ' -r1 '  + str(self._InnerRadius) +   ' -r2 '   + str(self._OuterRadius) +\
              ' -low ' + str(self._OuterRadius+2) + ' -high ' + str(self._OuterRadius+5)
              
      print '* ',command
      self.mylog.info("command:\t"+command)
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
      print '*  Applying geometrical  information in the headers. Needed for makespectra'
      command='xmipp_applygeo'+ \
              ' -i ' +os.path.basename(self._SelFileName) 
      print '* ',command
      self.mylog.info("command:\t"+command)
      program = os.system(command)     

   #------------------------------------------------------------------------
   #true_if_file_is_NOT_in_native_endian
   #------------------------------------------------------------------------
   def true_if_file_is_NOT_in_native_endian(self):
      import SelFiles,spider_header,os
      print '*********************************************************************'
      print '*  check if images are in native endian (spectra preprocesing)'
      mysel=SelFiles.selfile(os.path.basename(self._SelFileName))
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
      print '*  Change endian format. Needed for makespectra'
      command='xmipp_reverse_endian'+ \
              ' -i ' + os.path.basename(self._SelFileName) 
      print '* ',command
      self.mylog.info("command:\t"+command)
      program = os.system(command)     

   #------------------------------------------------------------------------
   #execute_spectra
   #------------------------------------------------------------------------
   def execute_spectra(self):
      import os, string
      print '*********************************************************************'
      print '*  Compute rotational power spectra'
      selFileName=os.path.basename(self._SelFileName)
      outFileName=(os.path.splitext(selFileName))[0] + '.sim'
      command='xmipp_makespectra'+ \
              ' -sel ' + selFileName + \
              ' -out ' + outFileName + \
              ' -x0 '  + str(self.xOffset) + \
              ' -y0 '  + str(self.yOffset) + \
              ' -r1 '  + str(self._InnerRadius) +   ' -r2 '   + str(self._OuterRadius) 
              
      print '* ',command
      self.mylog.info("command:\t"+command)
      program = os.system(command)     
      if(self._DisplayResults==True):
          print '*********************************************************************'
          print '*  Display spectra'
          selFileName=os.path.basename(self._SelFileName)
          spectraFileName=(os.path.splitext(selFileName))[0] + '.sim'
          command='xmipp_show'+ \
                  ' -spect ' + spectraFileName 
          print '* ',command
          self.mylog.info("command:\t"+command)
          program = os.system(command)     

  
   #------------------------------------------------------------------------
   #execute_KerDenSOM
   #------------------------------------------------------------------------
   def execute_KerDenSOM(self):
      import os

      selFileName=os.path.basename(self._SelFileName)
      spectraFileName=(os.path.splitext(selFileName))[0] + '.sim'
      print '*********************************************************************'
      print '*  compute kerdensom'
      kerdensom_out='kerd'
      command='xmipp_kerdensom'+ ' -din '  + spectraFileName + \
              ' -cout ' + kerdensom_out  + \
              ' '  + self._KerdensonExtraCommand
      print '* ',command
      self.mylog.info("command:\t"+command)
      os.system(command)
      if self._DisplayResults==True:
         command='xmipp_show -spectsom ' + \
                  kerdensom_out + \
                  ' -din ' + spectraFileName
         print '*********************************************************************'
         print '* ',command
         self.mylog.info("command:\t"+command)
         os.system(command)
      print "\n"

   #------------------------------------------------------------------------
   #clean and close everything
   #------------------------------------------------------------------------
   def close(self):
      print "Script finished"
      

#
# main
#     
if __name__ == '__main__':

    # create rotational_spectra_class object
    # 
   
   #init variables
   
   myspectra=rotational_spectra_class(_WorkDirectory,\
                                      _SelFileName,_DisplayResults,\
                                      _InnerRadius,_OuterRadius,\
                                      _Align2DExtraCommand,\
                                      _KerdensonExtraCommand
                                      )

    
