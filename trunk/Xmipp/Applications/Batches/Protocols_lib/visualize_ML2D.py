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
WorkingDir="ML2ref"
# Visualize last iterations class averages and output logfile to screen?
DoShowLastIter=True
# Visualize the class averages of all iterations in matrix-view?
DoMatrixAllIter=True
# Output statistics of all iterations to screen?
DoShowStatsAllIter=True

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
                 DoShowLastIter,
                 DoMatrixAllIter,
                 DoShowStatsAllIter):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=WorkingDir

        os.chdir(self.WorkingDir)
        if (DoShowLastIter):
            self.show_last_iter()
        if (DoMatrixAllIter):
            self.show_matrixview_all_iter()
        if (DoShowStatsAllIter):
            self.show_statistics_all_iter()

        # Return to parent dir
        os.chdir(os.pardir)

    def show_matrixview_all_iter(self):
        import os,glob
        import SelFiles
        
        selfiles=glob.glob('*_it?????.sel')
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            newsel=SelFiles.selfile()
            tmpsel=SelFiles.selfile()
            for selfile in selfiles:
                tmpsel.read(selfile)
                newsel.append(tmpsel.sellines)
            newsel.write('alliter.sel')
            # Get number of classes:
            fh=open(selfiles[0],'r')
            lines=fh.readlines()
            fh.close()
            nr_class=len(lines)
            command='xmipp_show -sel alliter.sel -w '+str(nr_class)+ ' & '
            print '* ',command
            os.system(command)

    def show_statistics_all_iter(self):
        import os,glob
        logfiles=glob.glob('*_it?????.log')
        if len(logfiles)==0:
            print "No logfiles yet. Visualize after job completion..."
        else:
            for logfile in logfiles:
                fh=open(logfile,'r')
                line1=fh.readline()
                line2=fh.readline()
                fh.close()
                words1=line1.split()
                words2=line2.split()
                print logfile,words1[6],words1[7],words1[8],words1[9],words2[1],words2[2],words2[3],words2[4]

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
                                        DoShowLastIter,
                                        DoMatrixAllIter,
                                        DoShowStatsAllIter)
    # close 
    visualize_ML2D.close()

