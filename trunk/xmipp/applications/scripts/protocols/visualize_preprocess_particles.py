#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_preprocess_particles.py
#
# Example use:
# python visualize_preprocess_particles.py protocol_preprocess_particles.py 
#
# This script requires that protocol_preprocess_particles.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Visualize the individual particles, sorted by micrograph?
DoShowParticles=True
# Visualize the individual particles, sorted by statistics?
""" Most common particles will be displayed at the top, and outliers will be displayed at the bottom of the selection file
"""
DoShowSortedParticles=False
# Visualize the Z-score used for particles sorting?
""" Particles will be sorted from low Z-scores (common particles) to high Z-scores (strange particles)
"""
DoShowZscore=False
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_particles_class:

    #init variables
    def __init__(self,
                 DoShowParticles,
                 DoShowSortedParticles,
                 DoShowZscore,
                 ProtocolName):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir)
        import visualization

        # Import the corresponding protocol, get WorkingDir and go there
        pardir=os.path.abspath(os.getcwd())
        shutil.copy(ProtocolName,'protocol.py')
        import protocol
        
        # all_images.sel is in ProjectDir!!
        ShowSelfiles=[]
        if (DoShowParticles):
            selfile=protocol.ProjectDir+ '/' + \
                     protocol.OutSelFile
            if protocol.IsPairList:
                selfile1=selfile.replace('.sel','_untilted.sel')
                selfile2=selfile.replace('.sel','_tilted.sel')
                ShowSelfiles.append(selfile1)
                ShowSelfiles.append(selfile2)
            else:
                ShowSelfiles.append(selfile)
            visualization.visualize_images(ShowSelfiles,True)
        # rest of stuff remains in WorkingDir!!
        os.chdir(protocol.WorkingDir)
        if (DoShowSortedParticles):
            ShowSelfiles.append('sort_junk.sel')
            visualization.visualize_images(ShowSelfiles,True)
        if (DoShowZscore):
            self.show_z_score()
        os.chdir(pardir)

        # Return to parent dir and remove protocol.py(c)
        if (os.path.exists('protocol.py')):
            os.remove('protocol.py')
        if (os.path.exists('protocol.pyc')):
            os.remove('protocol.pyc')

    def show_z_score(self):
        import os,sys
        import visualization
        if not os.path.exists('sort_junk.sumZ'):
            message=" Error: sort_junk.sumZ file does not exists!"
            print '* '+message
            sys.exit()
        else:
            fh=open('sort_junk.sumZ','r')
            lines=fh.readlines()
            fh.close()
            newlines=[]
            i = 0
            for line in lines:
                words=line.split()
                i = i + 1
                newlines.append(str(i)+' '+words[0]+'\n')
            fh=open('sort_junk.sumZ.dat','w')
            fh.writelines(newlines)
            fh.close()
            plot1=visualization.gnuplot()
            plot1.plot_xy_file('sort_junk.sumZ.dat',
                               'Z-score used for particle sorting (high values identify strange particles',
                               'particle number',
                               'sumZ')


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

    # create preprocess_particles_class object
    visualize=visualize_particles_class(DoShowParticles,
                                        DoShowSortedParticles,
                                        DoShowZscore,
                                        ProtocolName)
    visualize.close()

