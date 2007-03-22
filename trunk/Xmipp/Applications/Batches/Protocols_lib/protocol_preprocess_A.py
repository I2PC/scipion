#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script (part A) for Xmipp-based pre-processing of single-particles: 
#  - tif2raw conversion, 
#  - downsampling
#  - power spectral density (PSD) estimation on the micrograph
#
# For each file <micrograph>.tif given on the command line,
#  this script will perform the requested operations below.
# For each micrograph a subdirectory will be created
#
# Example use:
# ./xmipp_preprocess_A.py *.tif
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# All filenames.tif on which to perform preprocessing
ALLMICROGRAPHS="./*.tif"
#------------------------------------------------------------------------------------------------
# {section} Tiff2Raw
#------------------------------------------------------------------------------------------------
# Perform tiff2raw conversion?
DO_TIF2RAW=False
#------------------------------------------------------------------------------------------------
# {section} Downsampling
#------------------------------------------------------------------------------------------------
# Perform downsampling?
DO_DOWNSAMPLE=True
# Downsampling factor 
DOWN=3
#------------------------------------------------------------------------------------------------
# {section} Power Spectral Density estimation
#------------------------------------------------------------------------------------------------
# Perform power spectral density estimation?
DO_PSDESTIMATE=False
# Visualize estimated power spectral densities?
DO_SHOWPSD=True
#------------------------------------------------------------------------------------------------
# Particle picking
#------------------------------------------------------------------------------------------------
# This step needs to be performed manually.
# In each of the subdirectories created by this script, type:
# xmipp_mark -i downX_YYYY.raw &
# Save the coordinates of the manually picked particles in these sub-directories
# (Don't forget to save regularly, and always before closing the program!)
#
# Afterwards, you may use part B of this script to continue your pre-processing
#
#
# {expert} Log File Directory
LogPath="Logs"
#
# {expert} Log File Name (.log will be added)
LogFileName="G40P"
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class preprocess_A_class:

	#init variables
	def __init__(self,ALLMICROGRAPHS,
		     DO_TIF2RAW,
		     DO_DOWNSAMPLE,
		     DOWN,
		     DO_PSDESTIMATE,
		     DO_SHOWPSD,
		     LogPath='Logs',
		     LogName='logfile',
		     debug='off',
		     theSmtpLocalServer='mail.cnb.uam.es'
		     ):
	     
		import os
		import sys
		import string
		import log

		self.ALLMICROGRAPHS=ALLMICROGRAPHS
		self.DO_TIF2RAW=DO_TIF2RAW
		self.DO_DOWNSAMPLE=DO_DOWNSAMPLE
		self.DOWN=DOWN
		self.DO_PSDESTIMATE=DO_PSDESTIMATE
		self.DO_SHOWPSD=DO_SHOWPSD
		self.log = log.logpyl(LogPath,LogName,'off')                  

		# Execute program
		self.process_all_micrographs()

	def perform_tif2raw(self):
		import os
		oname=self.shortname+'/'+self.shortname+'.raw'
		print '*********************************************************************'
		print '*  Generating RAW for micrograph: '+self.name
		command='xmipp_tiff2raw '+self.name+' '+oname
		print '* ',command
		self.log.add("command\t:",command)
		os.system(command)
    
	def perform_downsample(self):
		import os
		iname=self.shortname+'/'+self.shortname+'.raw'
		oname=self.shortname+'/'+self.downname+'.raw'
		print '*********************************************************************'
		print '*  Downsampling micrograph: '+iname
		command='xmipp_downsample -i '+iname+' -o '+oname+' -output_bits 32 -Xstep '+str(self.DOWN)+' -kernel rectangle '+str(self.DOWN)+' '+str(self.DOWN)
		print '* ',command
		self.log.add("command\t:",command)
		os.system(command )

	def perform_psdestimate(self):
		import os
		iname=self.shortname+'/'+self.downname+'.raw'
		pname=self.shortname+'/'+self.shortname+'_psd.param'
		print '*********************************************************************'
		print '*  Estimate PSD for micrograph: '+iname
		paramlist = []
		paramlist.append('image= '+iname+'\n')
		paramlist.append('micrograph_averaging= yes \n')
		paramlist.append('N_horizontal= 512 \n')
		paramlist.append('periodogram= yes \n')
		paramlist.append('dont_adjust_CTF= yes \n')
		
		FILE = open(pname,"w")
		FILE.writelines(paramlist)
		FILE.close()
		command='xmipp_assign_CTF -i '+pname
		print '* ',command
		self.log.add("command\t:",command)
		os.system(command )
		oname=self.shortname+'/'+self.downname+'_Periodogramavg.psd'
		self.psdselfile.append(oname+' 1 \n')
		
	def visualize_psds(self):
		import os
		print '*********************************************************************'
		print '*  Visualizing all PSDs: '
		FILE = open("all_psds.sel","w")
		FILE.writelines(self.psdselfile)
		command='xmipp_show -psdsel all_psds.sel &'
		print '* ',command
		self.log.add("command\t:",command)
		os.system(command )
    
	def process_all_micrographs(self):
		import os
		import glob
		print '*********************************************************************'
		print '*  Pre-processing the following micrographs: '
		for self.filename in glob.glob(self.ALLMICROGRAPHS):
			(self.filepath, self.name) = os.path.split(self.filename)
			print '*  '+self.name

		self.psdselfile = []
		for self.filename in glob.glob(self.ALLMICROGRAPHS):
			(self.filepath, self.name) = os.path.split(self.filename)

			self.shortname=self.name.replace ( '.tif', '' )
			self.downname='down'+str(self.DOWN)+'_'+self.shortname

			if (os.path.exists(self.shortname)==False):
				self.log.add("command\t:",command)
				os.system('mkdir '+self.shortname)

			if (self.DO_TIF2RAW):
				self.perform_tif2raw()

			if (self.DO_DOWNSAMPLE):
				self.perform_downsample()

			if(os.path.exists(self.shortname+'/'+self.downname+'.raw')==False):
				self.downname=self.shortname

			if (self.DO_PSDESTIMATE):
				self.perform_psdestimate()

		if (self.DO_SHOWPSD):
		    self.visualize_psds()

		print '*  Done pre-processing all'
		print '*********************************************************************'

	
	def close(self):
		self.log.close()
#		
# Main
#     
if __name__ == '__main__':

   	# create preprocess_A_class object

	preprocessA=preprocess_A_class(ALLMICROGRAPHS,\
				       DO_TIF2RAW,\
				       DO_DOWNSAMPLE,\
				       DOWN,\
				       DO_PSDESTIMATE,\
				       DO_SHOWPSD,\
				       LogPath,\
				       LogFileName,\
				       'off',\
				       'mail.cnb.uam.es'\
				       )

	# close 
	preprocessA.close()

