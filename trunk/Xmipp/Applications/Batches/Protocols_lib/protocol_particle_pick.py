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
MicrographSelfile="all_micrographs.sel"
# Use GUI to display list of finished micrographs?
DoUseGui=True
# {expert} Root directory name for this project:
ProjectDir="/home/scheres/work/protocols"
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
                     MicrographSelfile,
                     DoUseGui,
                     ProjectDir,
                     LogDir):
	     
		import os,sys
                scriptdir=os.path.expanduser('~')+'/scripts/'
                sys.path.append(scriptdir) # add default search path
                self.SYSTEMSCRIPTDIR=scriptdir
                import log

                self.MicrographSelfile=MicrographSelfile
                self.ProjectDir=ProjectDir
                self.LogDir=LogDir
                self.DoUseGui=DoUseGui
                self.log=log.init_log_system(self.ProjectDir,
                                             self.LogDir,
                                             sys.argv[0],
                                             '.')
                
		# Execute program
		self.process_all_micrographs()

	def process_all_micrographs(self):
            import os
            import glob
            import SelFiles
            print '*********************************************************************'
            print '*  Perform manual particle picking for micrographs in: '+str(self.MicrographSelfile)
            print '*'
            print '* DONT FORGET TO SAVE YOUR COORDINATES REGULALRLY, AND ALWAYS BEFORE CLOSING!'
            print '*'

            # Read selfile with all micrographs to process
            self.mysel=SelFiles.selfile()
            self.mysel.read(self.MicrographSelfile)
            # Prepare have_picked.py file
            self.prepare_have_picked()

            for line in self.mysel.sellines:
                micrograph,state=line[:-1].split(" ")
                if not '-1' in state:
                    self.update_have_picked()
                    if not (self.have_already_picked(micrograph)):
                        self.perform_picking(micrograph)

        def prepare_have_picked(self):
            import os
            if not os.path.exists('have_picked.py'):

                fh=open('have_picked.py','w')
                self.havepickedlines=[]
                self.havepickedlines.append('#!/usr/bin/env python \n')
                self.havepickedlines.append('# {section} Select Yes when you have finished a micrograph, and Save \n')
                for line in self.mysel.sellines:
                    micrograph,state=line[:-1].split(" ")
                    downname=os.path.basename(micrograph)
                    downname=downname.replace('.raw','')
                    if not '-1' in state:
                        self.havepickedlines.append('# Finished picking '+micrograph+'?\n')
                        self.havepickedlines.append('Done_'+downname+'=False\n')
                self.havepickedlines.append('# {end-of-header} \n')
                fh.writelines(self.havepickedlines)
                fh.close()

            if (self.DoUseGui):
                print '* Select with which micrographs you have finished and save'
                command='python '+str(self.SYSTEMSCRIPTDIR)+'/protocol_gui.py have_picked.py &'
                os.system(command)
            else:
                print '* Edit manually the have_picked.py file to tell the program '
                print '* with which micrographs you have finished '
            
        def update_have_picked(self):
            self.havepickedlines=[]
            fh=open('have_picked.py','r')
            self.havepickedlines=fh.readlines()
            fh.close()

        def have_already_picked(self,name):
            import os,sys
            downname=os.path.basename(name)
            downname=downname.replace('.raw','')
            pattern='Done_'+downname
            for line in self.havepickedlines:
                if pattern in line:
                    args=line.split('=')
                    if 'True' in (args[1]):
                        return True
                    else:
                        return False
            message='Error: pattern '+pattern+'not found in have_picked.py'
            self.log.error(message)
            print '* ',message
            sys.exit()
 
        def perform_picking(self,name):
            import os
            directory,micrograph=os.path.split(name)
            os.chdir(directory)
            command='xmipp_mark -i '+micrograph
            print '* ',command
            self.log.info(command)
            os.system(command)
            os.chdir(os.pardir)

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

	particle_pick=particle_pick_class(MicrographSelfile,
                                          DoUseGui,
                                          ProjectDir,
                                          LogDir)

	# close 
	particle_pick.close()

