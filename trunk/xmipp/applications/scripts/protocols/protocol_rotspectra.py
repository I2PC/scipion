#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for rotational spectra classification
# using self-organizing maps
#
# Example use:
# ./xmipp_rotational_spectra.py
#
# Author:Roberto Marabini, March 2007
#        Carlos Oscar Sorzano, January 2011
#
# required files: log.py, spider_header.py, Sel_Files.py
# by defaults searchs for python files in ../python directory
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# {file} Selfile or stack with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName='all_images.sel'
# Working subdirectory: 
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='RotSpectra/test1'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {hidden} Directory name for logfiles:
LogDir='Logs'
#-----------------------------------------------------------------------------
# {section} Rotational spectra calculation
#-----------------------------------------------------------------------------
# Inner radius for rotational harmonics calculation:
""" These values are in pixels from the image center
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Makespectra
"""
SpectraInnerRadius=7
# Outer radius for rotational harmonics calculation:
SpectraOuterRadius=10
# {expert} Lowest harmonic to calculate
SpectraLowHarmonic=1
# {expert} Highest harmonic to calculate
SpectraHighHarmonic=15
#-----------------------------------------------------------------------------
# {section} Classification: classify_kerdensom 
#-----------------------------------------------------------------------------
# X-dimension of the self-organizing map:
SomXdim=7
# Y-dimension of the self-organizing map:
SomYdim=7
# {expert} Initial regularization factor:
""" The kerdenSOM algorithm anneals from an initial high regularization factor
    to a final lower one, in a user-defined number of steps.
    If the output map is too smooth, lower the regularization factors
    If the output map is not organized, higher the regularization factors
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
SomReg0=1000
# {expert} Final regularization factor:
SomReg1=200
# {expert} Number of steps to lower the regularization factor:
SomSteps=5
# {expert} Additional kerdenSOM parameters:
""" For a complete description 
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraCommand=''
#------------------------------------------------------------------------------------------------
# {hidden} Analysis of results
AnalysisScript='visualize_rotspectra.py'
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end-of-header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
import os,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path

class rotational_spectra_class:
   #init variables
   def __init__(self,
                SelFileName,
                WorkingDir,
                ProjectDir,
                SpectraInnerRadius,
                SpectraOuterRadius,
                SpectraLowHarmonic,
                SpectraHighHarmonic,
                SomXdim,
                SomYdim,
                SomReg0,
                SomReg1,
                SomSteps,                
                KerdensomExtraCommand):
       import log

       self.SelFileName=SelFileName
       self.WorkingDir=WorkingDir
       self.ProjectDir=ProjectDir
       self.SpectraInnerRadius=SpectraInnerRadius
       self.SpectraOuterRadius=SpectraOuterRadius
       self.SpectraLowHarmonic=SpectraLowHarmonic
       self.SpectraHighHarmonic=SpectraHighHarmonic
       self.SomXdim=_SomXdim
       self.SomYdim=_SomYdim
       self.SomReg0=_SomReg0
       self.SomReg1=_SomReg1
       self.SomSteps=_SomSteps              
       self.KerdensomExtraCommand=_KerdensomExtraCommand
       self.mylog=log.init_log_system(ProjectDir,
                                      LogDir,
                                      sys.argv[0],
                                      WorkingDir)
       
       #made backup of this script
       log.make_backup_of_script_file(sys.argv[0],self.WorkingDir)

       # Effectively run
       self.execute_findcenter()
       self.execute_spectra()
       self.execute_KerDenSOM()
  
   #------------------------------------------------------------------------
   #execute_findcenter
   #------------------------------------------------------------------------
   def execute_findcenter(self):
      import os, string
      print '*********************************************************************'
      print '* Looking for the position of the center of symmetry ...'
      command='xmipp_find_center2d'+ \
              ' -img ' + self.SelFileName + \
              ' -r1 '  + str(self.SpectraInnerRadius) +   ' -r2 '   + str(self.SpectraOuterRadius) +\
              ' -low ' + str(self.SpectraOuterRadius+2) + ' -high ' + str(self.SpectraOuterRadius+5)
              
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
   #execute_spectra
   #------------------------------------------------------------------------
   def execute_spectra(self):
      import os, string
      import launch_job
      if (os.path.exists(self._SpectraName)):
          os.remove(self._SpectraName)
      print '*********************************************************************'
      print '* Computing rotational power spectra'
      selFileName=os.path.basename(self._SelFileName)
      outFileName=(os.path.splitext(selFileName))[0] + '.sim'
      command=' -i ' + selFileName + \
              ' -o ' + str(self._SpectraName) + \
              ' -x0 '  + str(self.xOffset) + \
              ' -y0 '  + str(self.yOffset) + \
              ' -r1 '  + str(self._SpectraInnerRadius) + \
              ' -r2 '   + str(self._SpectraOuterRadius) + \
              ' -low ' + str(self._SpectraLowHarmonic) + \
              ' -high ' + str(self._SpectraHighHarmonic)
      launch_job.launch_job("xmipp_make_spectra",
                            command,
                            self.mylog,
                            False,1,1,'')
  
   #------------------------------------------------------------------------
   #execute_KerDenSOM
   #------------------------------------------------------------------------
   def execute_KerDenSOM(self):
      import launch_job

      selFileName=os.path.basename(self._SelFileName)
      print '*********************************************************************'
      print '* Computing kerdensom ...'
      command=' -verb 1 -i '  + str(self._SpectraName) + \
              ' -o '    + str(self._SomName)  + \
              ' -xdim ' + str(self._SomXdim) + \
              ' -ydim ' + str(self._SomYdim) + \
              ' -reg0 ' + str(self._SomReg0) + \
              ' -reg1 ' + str(self._SomReg1) + \
              ' -steps ' + str(self._SomSteps) + \
              ' '  + str(self._KerdensomExtraCommand)
      launch_job.launch_job("xmipp_classify_kerdensom",
                            command,
                            self.mylog,
                            False,1,1,'')

#
# main
#     
if __name__ == '__main__':
   myspectra=rotational_spectra_class(SelFileName,
                                      WorkingDir,
                                      ProjectDir,
                                      SpectraInnerRadius,
                                      SpectraOuterRadius,
                                      SpectraLowHarmonic,
                                      SpectraHighHarmonic,
                                      SomXdim,
                                      SomYdim,
                                      SomReg0,
                                      SomReg1,
                                      SomSteps,                
                                      KerdensomExtraCommand)