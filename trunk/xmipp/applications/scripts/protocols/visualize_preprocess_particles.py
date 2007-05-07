#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_preprocess_particles.py
#
# Example use:
# python visualize_preprocess_particles.py
#
# This script requires that protocol_preprocess_particles.py is in the current directory
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
# Visualize the individual particles, sorted by micrograph?
DoShowParticles=True
# Visualize the individual particles, sorted by statistics?
""" Most common particles will be displayed at the top, and outliers will be displayed at the bottom of the selection file"""
DoShowSortedParticles=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_particles_class:

    #init variables
    def __init__(self,
                 DoShowCTFs,
                 DoShowParticles,
                 DoShowSortedParticles,
                 WorkingDir):
	     
        import os,sys
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir)
        import visualization


        self.WorkingDir=WorkingDir
        os.chdir(self.WorkingDir)

        # import corresponding protocol
        sys.path.append(self.WorkingDir)
        import protocol_preprocess_particles_backup

        ShowSelfiles=[]
        if (DoShowCTFs):
            self.visualize_ctfs()
        if (DoShowParticles):
            selfile=protocol_preprocess_particles_backup.ProjectDir+ '/' + \
                     protocol_preprocess_particles_backup.OutSelFile
            ShowSelfiles.append(selfile)
        if (DoShowSortedParticles):
            ShowSelfiles.append('sort_junk.sel')

        visualization.visualize_images(ShowSelfiles,True)

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
                                        DoShowParticles,
                                        DoShowSortedParticles,
                                        WorkingDir)
    visualize.close()

