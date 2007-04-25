#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_ML2D.py
#
# Example use:
# python visualize_ml2d.py
#
# This script requires that protocol_ml2d.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Visualize the class averages of all iterations in matrix-view?
DoMatrixAllIter=False
# Visualize class averages of the last iterations separately?
DoShowLastIter=False
# Plot model (and mirror) fractions of the last iertation?
DoShowFractions=False
# Plot convergence statistics for all iterations?
DoShowStatsAllIter=False

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_ML2D_class:

    #init variables
    def __init__(self,
                 DoMatrixAllIter,
                 DoShowLastIter,
                 DoShowFractions,
                 DoShowStatsAllIter):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        # import corresponding protocol
        import protocol_ml2d
        self.WorkingDir=protocol_ml2d.WorkingDir

        os.chdir(self.WorkingDir)
        if (DoMatrixAllIter):
            self.show_matrixview_all_iter()
        if (DoShowLastIter):
            self.show_last_iter()
        if (DoShowFractions):
            self.show_model_fractions()
        if (DoShowStatsAllIter):
            self.show_LL_PmaxsumP_all_iter()

        # Return to parent dir
        os.chdir(os.pardir)

    def show_matrixview_all_iter(self):
        import os,glob
        import selfile
        
        selfiles=glob.glob('*_it?????.sel')
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            newsel=selfile.selfile()
            tmpsel=selfile.selfile()
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

    def show_LL_PmaxsumP_all_iter(self):
        import os,glob
        import visualization
        import docfiles
        logfiles=glob.glob('*_it?????.log')
        if len(logfiles)==0:
            print "No logfiles yet. Visualize after job completion..."
        else: 
            newlines1=[]
            newlines2=[]
            newlines3=[]
            newlines1.append('# iter | log-likelihood \n')
            newlines2.append('# iter | Pmax/sumP \n')
            newlines3.append('# iter | signal change \n')
            conv=0.05
            for logfile in logfiles:
                doc=docfiles.docfile(logfile)
                eps=doc.maximum_of_column(4)
                fh=open(logfile,'r')
                line1=fh.readline()
                line2=fh.readline()
                line3=fh.readline()
                fh.close()
                words1=line1.split()
                words2=line2.split()
                words3=line2.split()
                # find convergence threshold
                for i in range(len(words3)):
                    if (words2[i]=='-eps'):
                        conv=float(words3[i+1])/1000.
                # get iteration
                iter=int(words2[6])-1
                # get data
                newlines1.append(str(iter)+' '+words1[7]+'\n')
                newlines2.append(str(iter)+' '+words1[9]+'\n')
                newlines3.append(str(iter)+' '+str(eps)+'\n')
            fh=open('alliter_LL.dat','w')
            fh.writelines(newlines1)
            fh.close();
            fh=open('alliter_Pmax.dat','w')
            fh.writelines(newlines2)
            fh.close();
            fh=open('alliter_signalchange.dat','w')
            fh.writelines(newlines3)
            fh.close();
            plot1=visualization.gnuplot()
            plot1.plot_xy_file('alliter_LL.dat',
                              'Log-likelihood target function (should increase)',
                              'iterations',
                              'LL')
            plot2=visualization.gnuplot()
            plot2.plot_xy_file('alliter_Pmax.dat',
                              'The width of the probability distributions (often goes to one i.e. delta-functions)',
                              'iterations',
                              'Pmax/sumP')
            plot3=visualization.gnuplot()
            plot3.plot_xy_file('alliter_signalchange.dat',
                              'The maximum signal change over all references (green line is convergence)',
                              'iterations',
                              'signal change')
            plot3.send('set yrange [0:1]')
            plot3.send("plot 'alliter_signalchange.dat' with lines")
            plot3.send('replot '+str(conv))
                             
    def show_model_fractions(self):
        import os,glob
        import visualization
        import docfiles
        logfiles=glob.glob('*_it?????.log')
        if len(logfiles)==0:
            print "No logfiles yet. Visualize after job completion..."
        else:
            plot1=visualization.gnuplot()
            plot1.plot_xy1y2_file(logfiles[-1],
                                  'Model and mirror fractions in the last iteration',
                                  'model fraction',
                                  'mirror fraction',
                                  'reference number',
                                  'fractions',
                                  1,3,4)

    def show_last_iter(self):
        # Visualize class averages:
        import os,glob
        logfiles=glob.glob('*_it?????.log')
        if len(logfiles)==0:
            print "No logfiles yet. Visualize after job completion..."
        else:
            # print logfile to screen:
            lastlogfile=logfiles[-1]
            fh=open(lastlogfile,'r')
            loglines=fh.readlines()
            print '*********************************************************************'
            print "Logfile "+str(lastlogfile)+": "
            for line in loglines:
                print line[:-1]
            
            # Display last selfile
            lastselfile=lastlogfile.replace('.log','.sel')
            command='xmipp_show -sel '+str(lastselfile)+ ' & '
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

    # create ML2D_class object
    visualize_ML2D=visualize_ML2D_class(DoMatrixAllIter,
                                        DoShowLastIter,
                                        DoShowFractions,
                                        DoShowStatsAllIter)
    # close 
    visualize_ML2D.close()

