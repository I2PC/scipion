#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script (part A) for Xmipp-based manual particle picking
#
# For each micrograph in the MicrographSelfile, this program will launch
# the xmipp_mark program 
# A graphical interface exists to identify micrographs that have been finished
#
# Example use:
# python xmipp_particle_pick.py &
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Selfile with all micrographs to pick particles from (in current dir):
MicrographSelfile="pairs.sel"
# Is this selfile a list of untilted-tilted pairs?
""" True for RCT-processing. In that case, provide a 3-column selfile as follows:
    untilted_pair1.raw tilted_pair1.raw 1
    untilted_pair2.raw tilted_pair2.raw 1
    etc...
    Where 1 in the third column means active pair, and -1 means inactive pair
"""
IsPairList=True
# Use GUI to display list of finished micrographs?
DoUseGui=True
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
                     MicrographSelfile,
                     IsPairList,
                     DoUseGui,
                     ProjectDir,
                     LogDir):
	     
		import os,sys
                scriptdir=os.path.expanduser('~')+'/scripts/'
                sys.path.append(scriptdir) # add default search path
                self.SYSTEMSCRIPTDIR=scriptdir
                import log

                self.MicrographSelfile=MicrographSelfile
                self.IsPairList=IsPairList
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
            import selfile
            print '*********************************************************************'
            print '*  Perform manual particle picking for micrographs in: '+str(self.MicrographSelfile)
            print '*'
            print '* DONT FORGET TO SAVE YOUR COORDINATES REGULALRLY, AND ALWAYS BEFORE CLOSING!'
            if (self.IsPairList):
                print '* AND ALSO SAVE THE ANGLES IN THE UNTILTED MICROGRAPHS!'
            print '*'

            if not self.IsPairList:
                self.mysel=selfile.selfile()
                self.mysel.read(self.MicrographSelfile)
                self.prepare_have_picked()
                
                for micrograph,state in self.mysel.sellines:
                    if (state.find('-1') == -1):
                        self.update_have_picked()
                        if not (self.have_already_picked(micrograph)):
                            self.perform_picking(micrograph)

            else:
                fh=open(self.MicrographSelfile,'r')
                self.pairlines=fh.readlines()
                fh.close()
                words=self.pairlines[0].split()
                if (len(words)<3):
                    message='Error: Selfile is not a pairlist file!'
                    print '*',message
                    self.log.error(message)
                    sys.exit()
                else:
                    self.prepare_have_picked()
                    
                for line in self.pairlines:
                    words=line.split()
                    untilted=words[0]
                    tilted=words[1]
                    state=words[2]
                    if (state.find('-1') == -1):
                        self.update_have_picked()
                        if not (self.have_already_picked(untilted)):
                            print 'untilted= ',untilted,tilted
                            self.perform_picking_pair(untilted,tilted)

        def prepare_have_picked(self):
            import os
            if not os.path.exists('have_picked.py'):

                fh=open('have_picked.py','w')
                self.havepickedlines=[]
                self.havepickedlines.append('#!/usr/bin/env python \n')
                self.havepickedlines.append('# {section} Select Yes when you have finished a micrograph, and Save \n')
                if not self.IsPairList:
                    for line in self.mysel.sellines:
                        micrograph,state=line[:-1].split(" ")
                        downname=os.path.basename(micrograph)
                        downname=downname.replace('.raw','')
                        self.append_line_have_picked(downname,state)
                else:
                    for line in self.pairlines:
                        words=line.split()
                        untilted=os.path.basename(words[0])
                        untilted=untilted.replace('.raw','')
                        tilted=words[1]
                        state=words[2]
                        self.append_line_have_picked(untilted,state)
   
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
            
        def append_line_have_picked(self,name,state):
            if (state.find('-1') == -1):
                if (self.IsPairList):
                    self.havepickedlines.append('# Finished picking pair with '+name+'?\n')
                else:
                    self.havepickedlines.append('# Finished picking micrograph '+name+'?\n')
                self.havepickedlines.append('Done_'+name+'=False\n')

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
                if (line.find(pattern) > -1):
                    args=line.split('=')
                    if (args[1].find('True') > -1):
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


        def perform_picking_pair(self,untilted,tilted):
            import os
            directory,uname=os.path.split(untilted)
            tname='../'+tilted
            os.chdir(directory)
            command='xmipp_mark -i '+uname+' -tilted '+tname
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
                                          IsPairList,
                                          DoUseGui,
                                          ProjectDir,
                                          LogDir)

	# close 
	particle_pick.close()

