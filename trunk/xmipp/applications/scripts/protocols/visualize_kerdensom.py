#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_kerdensom.py
#
# Example use:
# python visualize_kerdensom.py
#
# This script requires that protocol_kerdensom.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
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
                 DoShowSOM,
                 DoShowInfFile,
                 ProtocolName):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        # Import the corresponding protocol, get WorkingDir and go there
        pardir=os.path.abspath(os.getcwd())
        shutil.copy(ProtocolName,'protocol.py')
        import protocol
        self.WorkingDir=protocol.WorkingDir
        os.chdir(self.WorkingDir)
        self.SomName=protocol.SomName

        if (DoShowSOM):
            self.visualize_SOM()
        if (DoShowInfFile):
            self.output_infofile()
            
        # Return to parent dir and remove protocol.py(c)
        os.chdir(pardir)
        if (os.path.exists('protocol.py')):
            os.remove('protocol.py')
        if (os.path.exists('protocol.pyc')):
            os.remove('protocol.pyc')

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

    import sys
    ProtocolName=sys.argv[1]

    # create kerdensom_class object
    visualize_kerdensom=visualize_kerdensom_class(DoShowSOM,
                                                  DoShowInfFile,
                                                  ProtocolName)
    # close 
    visualize_kerdensom.close()

