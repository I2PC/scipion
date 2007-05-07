#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_ctf_estimation.py
#
# Example use:
# python visualize_ctf_estimation.py
#
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Visualize the estimated CTF models for all micrographs?
""" You may re-estimate CTFs by selecting a micrograph, and selecting re-estimate CTF from the right-mouse button menu 
"""
DoShowCTFs=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_particles_class:

    #init variables
    def __init__(self,
                 DoShowCTFs):
	     
        import os,sys
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir)

        self.WorkingDir=WorkingDir
        os.chdir(self.WorkingDir)

        ShowSelfiles=[]
        if (DoShowCTFs):
            self.visualize_ctfs()

        # Return to parent dir
        os.chdir(os.pardir)
        
            
    def visualize_ctfs(self):
        import os
        print '*********************************************************************'
        print '*  Visualizing all CTFs: '
        command='xmipp_show -ctfsel all_ctfs.sel &'
        print '* ',command
        print '* '
        print '* You may re-estimate CTFs by selecting the corresponding micrograph (double-clicking) '
        print '* and then selecting re-estimate CTF from the right-mouse button menu '
        print '* Adjust the red circle to coincide with the first zero in the PSD '
        print '* After re-esimtation re-run the script (skipping CTF-estimation, '
        print '* but performing extraction, normalization, phase flipping etc.) '
        print '* '
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
    WorkingDir=sys.argv[1]

    # create preprocess_particles_class object
    visualize=visualize_particles_class(DoShowCTFs,
                                        WorkingDir)
    visualize.close()

