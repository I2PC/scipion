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
# Visualize the individual particles, sorted by micrograph?
DoShowParticles=True
# Visualize the individual particles, sorted by statistics?
""" Most common particles will be displayed at the top, and outliers will be displayed at the bottom of the selection file
"""
DoShowSortedParticles=True
# Visualize the Z-score used for particles sorting?
""" Particles will be sorted from low Z-scores (common particles) to high Z-scores (strange particles)
"""
DoShowZscore=True
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
        if (DoShowParticles):
            selfile=protocol_preprocess_particles_backup.ProjectDir+ '/' + \
                     protocol_preprocess_particles_backup.OutSelFile
            ShowSelfiles.append(selfile)
        if (DoShowSortedParticles):
            ShowSelfiles.append('sort_junk.sel')
        if (DoShowZscore):
            self.show_z_score()

        visualization.visualize_images(ShowSelfiles,True)

        # Return to parent dir
        os.chdir(os.pardir)
        

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
    WorkingDir=sys.argv[1]

    # create preprocess_particles_class object
    visualize=visualize_particles_class(DoShowParticles,
                                        DoShowSortedParticles,
                                        DoShowZscore,
                                        WorkingDir)
    visualize.close()

