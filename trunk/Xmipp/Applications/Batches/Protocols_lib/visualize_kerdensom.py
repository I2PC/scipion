#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_kerdensom.py
#
# Example use:
# python visualize_kerdensom.py
#
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Working subdirectory:
WorkingDir="SOM_ML2ref_ref1"
# Name of the self-organizing map:
SomName="som"
# Visualize the self-organizing map?
DoShowSOM=True
# Output the som infofile to screen?
DoShowInfFile=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_kerdensom_class:

    #init variables
    def __init__(self,
                 WorkingDir,
                 SomName,
                 DoShowSOM,
                 DoShowInfFile):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=WorkingDir
        self.SomName=SomName

        os.chdir(self.WorkingDir)
        if (DoShowSOM):
            self.visualize_SOM()
        if (DoShowInfFile):
            self.output_infofile()
            
        # Return to parent dir
        os.chdir(os.pardir)

    def visualize_SOM(self):
         import os
         command='xmipp_show -som ' + str(self.SomName)+' &'
         print '* ',command
         os.system(command)

    def output_infofile(self):
        import os
        command='cat ' + str(self.SomName)+'.inf'
        print '* ',command
        os.system(command)


    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    # create kerdensom_class object
    visualize_kerdensom=visualize_kerdensom_class(WorkingDir,
                                                  SomName,
                                                  DoShowSOM,
                                                  DoShowInfFile)
    # close 
    visualize_kerdensom.close()

