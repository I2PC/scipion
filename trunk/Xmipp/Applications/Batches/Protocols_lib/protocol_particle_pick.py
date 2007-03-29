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
AllMicrographs="./*.tif"
# Reset all "have_picked" flags
DoResetHavePicked=False
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
LogDir="Logs"
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class particle_pick_class:

	#init variables
	def __init__(self,
                     AllMicrographs,
                     DoResetHavePicked,
                     ProjectDir,
                     LogDir):
	     
		import os,sys
                scriptdir=os.path.expanduser('~')+'/scripts/'
                sys.path.append(scriptdir) # add default search path
                import log

                self.AllMicrographs=AllMicrographs
                self.ProjectDir=ProjectDir
                self.LogDir=LogDir
		self.log=log.init_log_system(self.ProjectDir,
                                             self.LogDir,
                                             sys.argv[0],
                                             '.')
                
		# Execute program
		self.process_all_micrographs()

	def process_all_micrographs(self):
		import os
		import glob
		print '*********************************************************************'
		print '*  Perform manual particle picking in the following micrographs: '
		for self.filename in glob.glob(self.AllMicrographs):
			(self.filepath, self.name) = os.path.split(self.filename)
			print '*  '+self.name

		for self.filename in glob.glob(self.AllMicrographs):
			(self.filepath, self.name) = os.path.split(self.filename)

			self.shortname=self.name.replace ( '.tif', '' )
			self.downname='down'+str(self.Down)+'_'+self.shortname
                        if not os.path.exists(self.shortname+'/'+self.downname+'.raw'):
                            self.downname=self.shortname

                        if (self.DoResetHavePicked):
                            self.reset_have_picked(self.downname)
                        else:
                            self.perform_picking(self.downname)

        def reset_have_picked(self,name):
            if os.path.exists(name+'.have_picked'):
                os.remove(name+'.have_picked')


        def perform_picking(self,name):
            directory,micrograph=os.path.split(name)
            os.chdir(directory)


            os.chdir(pardir)

	def close(self):
                message=" You have done picking all particles! :-)"
		print '* ',message
		print '*********************************************************************'
                self.log.info(message)
#		
# Main
#     
if __name__ == '__main__':

   	# create preprocess_A_class object

	particle_pick=particle_pick_class(AllMicrographs,
                                       DoResetHavePicked,
                                       ProjectDir,
                                       LogDir)

	# close 
	particle_pick.close()

