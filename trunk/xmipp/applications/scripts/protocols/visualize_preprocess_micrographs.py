#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_preprocess_micrographs.py
#
# Example use:
# python visualize_preprocess_micrographs.py protocol_preprocess_micrographs.py
#
# This script requires that protocol_preprocess_micrographs.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Visualize the estimated CTFs of all micrographs?
DoShowCTFs=True
# {expert} Visualize the PSDs of all micrographs?
DoShowPSDs=False
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_micrographs_class:

    #init variables
    def __init__(self,
                 DoShowCTFs,
                 DoShowPSDs,
                 ProtocolName):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path

        # Import the corresponding protocol, get WorkingDir and go there
        import protocol
        protocol.__name__=ProtocolName.replace('py','')
        self.WorkingDir=protocol.WorkingDir
        os.chdir(self.WorkingDir)
        
        self.MicrographSelfile=protocol.MicrographSelfile

        if (DoShowPSDs):
            self.visualize_psds()
            self.update_micrograph_selfile('all_psds.sel')

        if (DoShowCTFs):
            self.visualize_ctfs()
            self.update_micrograph_selfile('all_ctfs.sel')
        
        # Return to parent dir
        os.chdir(os.pardir)
            
    def visualize_psds(self):
        import os
        print '*********************************************************************'
        print '*  Visualizing all PSDs: '
        command='xmipp_show -psdsel all_psds.sel'
        print '* ',command
        print '* '
        print '* You may select bad micrographs and save them as discarded.'
        print '* Make sure to use save this file again as "all_psds.sel'
        print '* This way, the discarded micrographs will no longer be processed!'
        os.system(command )

    def visualize_ctfs(self):
        import os
        print '*********************************************************************'
        print '*  Visualizing all CTFs: '
        command='xmipp_show -sel all_ctfs.sel'
        print '* ',command
        print '* '
        print '* You may select bad micrographs and save them as discarded.'
        print '* Make sure to use save this file again as "all_ctfs.sel'
        print '* This way, the discarded micrographs will no longer be processed!'
        os.system(command )

    def update_micrograph_selfile(self,ctfselfile):
        import sys
        import selfile
        
        sel1=selfile.selfile()
        sel2=selfile.selfile()
        sel3=selfile.selfile()
        sel1.read(self.MicrographSelfile)
        sel2.read(ctfselfile)
        
        if not (sel1.length_even_no_actives() == sel2.length_even_no_actives()):
            message='Error: '+ctfselfile+' and '+str(self.MicrographSelfile)+' are of unequal length'
            print '* ',message
            sys.exit()
        else:
            for i in range(len(sel1.sellines)):
                name1,state1=sel1.sellines[i]
                name2,state2=sel2.sellines[i]
                sel3.insert(name1,state2)
            sel3.write(str(self.MicrographSelfile))
            message='Updated '+str(self.MicrographSelfile)
            print '* ',message

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

    # create preprocess_micrographs_class object
    visualize=visualize_micrographs_class(DoShowCTFs,
                                          DoShowPSDs,
                                          ProtocolName)
    visualize.close()

