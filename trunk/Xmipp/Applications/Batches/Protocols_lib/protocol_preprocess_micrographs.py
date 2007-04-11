#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs: 
#  - tif2raw conversion, 
#  - downsampling
#  - power spectral density (PSD) estimation on the micrograph
#
# For each file <micrograph>.tif given, this script will perform 
# the requested operations below.
# For each micrograph a subdirectory will be created
#
# Example use:
# ./xmipp_preprocess_micrographs.py *.tif
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# All filenames.tif on which to perform preprocessing
AllMicrographs="./*.tif"
# Output selfile with preprocessed micrographs
MicrographSelfile="all_micrographs.sel"
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Tiff2Raw
#------------------------------------------------------------------------------------------------
# Perform tiff2raw conversion?
DoTif2Raw=True
#------------------------------------------------------------------------------------------------
# {section} Downsampling
#------------------------------------------------------------------------------------------------
# Perform downsampling?
DoDownSample=True
# Downsampling factor 
Down=3
#------------------------------------------------------------------------------------------------
# {section} Power Spectral Density estimation
#------------------------------------------------------------------------------------------------
# Perform power spectral density estimation?
DoEstimatePSD=True
# Visualize estimated power spectral densities?
DoVisualizePSD=True
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
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class preprocess_A_class:

	#init variables
	def __init__(self,
                     AllMicrographs,
                     MicrographSelfile,
                     ProjectDir,
                     LogDir,
		     DoTif2Raw,
		     DoDownSample,
		     Down,
		     DoEstimatePSD,
		     DoVisualizePSD):
	     
		import os,sys
                scriptdir=os.path.expanduser('~')+'/scripts/'
                sys.path.append(scriptdir) # add default search path
                import log

                self.AllMicrographs=AllMicrographs
                self.MicrographSelfile=MicrographSelfile
                self.ProjectDir=ProjectDir
                self.LogDir=LogDir
		self.DoTif2Raw=DoTif2Raw
		self.DoDownSample=DoDownSample
		self.Down=Down
		self.DoEstimatePSD=DoEstimatePSD
		self.DoVisualizePSD=DoVisualizePSD
		self.log=log.init_log_system(self.ProjectDir,
                                             self.LogDir,
                                             sys.argv[0],
                                             '.')
                
		# Execute program
		self.process_all_micrographs()

	def perform_tif2raw(self):
		import os
		oname=self.shortname+'/'+self.shortname+'.raw'
		print '*********************************************************************'
		print '*  Generating RAW for micrograph: '+self.name
		command='xmipp_tiff2raw '+self.filename+' '+oname
		print '* ',command
		self.log.info(command)
		os.system(command)
    
	def perform_downsample(self):
		import os
		iname=self.shortname+'/'+self.shortname+'.raw'
		oname=self.shortname+'/'+self.downname+'.raw'
		print '*********************************************************************'
		print '*  Downsampling micrograph: '+iname
		command='xmipp_downsample -i '+iname+' -o '+oname+' -output_bits 32 -Xstep '+str(self.Down)+' -kernel rectangle '+str(self.Down)+' '+str(self.Down)
		print '* ',command
		self.log.info(command)
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
		self.log.info(command)
		os.system(command )
		oname=self.shortname+'/'+self.downname+'_Periodogramavg.psd'
		self.psdselfile.append(oname+' 1 \n')
                FILE = open("all_psds.sel","w")
		FILE.writelines(self.psdselfile)
		FILE.close()

	def visualize_psds(self):
		import os
		print '*********************************************************************'
		print '*  Visualizing all PSDs: '
		command='xmipp_show -psdsel all_psds.sel'
		print '* ',command
                print '* You may select bad micrographs and save them as discarded.'
                print '* Make sure to use save this file again as "all_psds.sel'
                print '* This way, the discarded micrographs will no longer be processed!'
		self.log.info(command)
		os.system(command )
    
        def append_micrograph_selfile(self):
            name=self.shortname+'/'+self.downname+'.raw'
            self.micselfile.append(name+' 1 \n')
            FILE = open(self.MicrographSelfile,'w')
            FILE.writelines(self.micselfile)
            FILE.close()

        def update_micrograph_selfile(self):
            fh=open(self.MicrographSelfile,'r')
            miclines=fh.readlines()
            fh.close()
            fh=open('all_psds.sel','r')
            psdlines=fh.readlines()
            fh.close()
            if (len(miclines)!=len(psdlines)):
                message='Error: all_psds.sel and '+str(self.MicrographSelfile)+' are of unequal length'
                self.log.error(message)
                print '* ',message
                sys.exit()
            else:
                newlines=[]
                for i in range(len(miclines)):
                    args1=miclines[i].split(' ')
                    args2=psdlines[i].split(' ')
                    newlines.append(args1[0]+' '+args2[1])
                fh=open(self.MicrographSelfile,'w')
                fh.writelines(newlines)
                fh.close()
                message='Updated '+str(self.MicrographSelfile)
                print '* ',message
                self.log.info(message)

	def process_all_micrographs(self):
		import os
		import glob
		print '*********************************************************************'
		print '*  Processing the following micrographs: '
		for self.filename in glob.glob(self.AllMicrographs):
			(self.filepath, self.name) = os.path.split(self.filename)
			print '*  '+self.name

		self.psdselfile = []
                self.micselfile = []
		for self.filename in glob.glob(self.AllMicrographs):
			(self.filepath, self.name) = os.path.split(self.filename)

			self.shortname=self.name.replace ( '.tif', '' )
			self.downname='down'+str(self.Down)+'_'+self.shortname

                        if not os.path.exists(self.shortname):
                            os.makedirs(self.shortname)

			if (self.DoTif2Raw):
                            self.perform_tif2raw()

			if (self.DoDownSample):
                            self.perform_downsample()

                        if not os.path.exists(self.shortname+'/'+self.downname+'.raw'):
                            self.downname=self.shortname

			if (self.DoEstimatePSD):
                            self.perform_psdestimate()

                        self.append_micrograph_selfile()

		if (self.DoVisualizePSD):
		    self.visualize_psds()
                    self.update_micrograph_selfile()

	def close(self):
                message=" Done pre-processing all"
		print '* ',message
		print '*********************************************************************'
                self.log.info(message)
#		
# Main
#     
if __name__ == '__main__':

   	# create preprocess_A_class object

	preprocessA=preprocess_A_class(AllMicrographs,
                                       MicrographSelfile,
                                       ProjectDir,
                                       LogDir,
				       DoTif2Raw,
				       DoDownSample,
				       Down,
				       DoEstimatePSD,
				       DoVisualizePSD)

	# close 
	preprocessA.close()

