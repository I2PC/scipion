#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_ML2D.py
#
# Example use:
# python visualize_ML2D.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Working subdirectory:
WorkingDir="ML3ref"
# Visualize last iterations class averages and show corresponding logfile?
DoShowLastIter=False
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_ML2D_class:

    #init variables
    def __init__(self,
                 WorkingDir,
                 DoShowLastIter):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=WorkingDir

        os.chdir(self.WorkingDir)
        if (DoShowLastIter):
            self.show_last_iter()

        # Return to parent dir
        os.chdir(os.pardir)

    def show_last_iter(self):
        # Visualize class averages:
        import os,glob
        selfiles=glob.glob('*_it?????.sel')
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            lastselfile=selfiles[-1]
            command='xmipp_show -sel '+str(lastselfile)+ ' & '
            print '* ',command
            os.system(command)

            # print logfile to screen:
            logfiles=glob.glob('*_it?????.log')
            lastlogfile=logfiles[-1]
            fh=open(lastlogfile,'r')
            loglines=fh.readlines()
            print "Logfile "+str(lastlogfile)+": "
            for line in loglines:
                print line[:-1]

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    # create ML2D_class object

    visualize_ML2D=visualize_ML2D_class(WorkingDir,
                                        DoShowLastIter)
        # close 
    visualize_ML2D.close()

