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
DoMatrixAllIter=True
# Separately visualize class averages of the last iteration?
DoShowLastIter=True
# Plot model (and mirror) fractions of the last iteration?
DoShowFractions=True
# Plot convergence statistics for all iterations?
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
                 DoMatrixAllIter,
                 DoShowLastIter,
                 DoShowFractions,
                 DoShowStatsAllIter,
                 WorkingDir):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        
        self.WorkingDir=WorkingDir
        os.chdir(self.WorkingDir)
        # Example how to import corresponding protocol
        #sys.path.append(self.WorkingDir)
        #import protocol_ml2d_backup

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
            newlines1.append('# iter | delta log-likelihood \n')
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
                    if (words3[i]=='-eps'):
                        conv=float(words3[i+1])/1000.
                # get iteration
                iter=int(words2[6])-1
                # get data
                newlines1.append(str(iter)+' '+words1[7]+'\n')
                newlines2.append(str(iter)+' '+words1[9]+'\n')
                newlines3.append(str(iter)+' '+str(eps)+'\n')
            fh=open('alliter_dLL.dat','w')
            fh.writelines(newlines1)
            fh.close();
            fh=open('alliter_Pmax.dat','w')
            fh.writelines(newlines2)
            fh.close();
            fh=open('alliter_signalchange.dat','w')
            fh.writelines(newlines3)
            fh.close();
            plot1=visualization.gnuplot()
            plot1.plot_xy_file('alliter_dLL.dat',
                              'Log-likelihood gain (should be positive, converges to zero)',
                              'iterations',
                              'dLL')
            plot1.send('replot 0')
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
            # Get number of classes:
            fh=open(logfiles[-1],'r')
            lines=fh.readlines()
            fh.close()
            nr_class=0.5
            for line in lines:
                if (line.find(';')<0):
                    nr_class+=1.

            plot1=visualization.gnuplot()
            plot1.prepare_empty_plot('Model and mirror fractions in the last iteration',
                                     'reference number',
                                     'fractions')
            plot1.send(" set yrange [0:1]")
            plot1.send(" set xrange [0.5:"+str(nr_class)+"]")
            plot1.send(" plot '" + logfiles[-1] + "' using 1:3 title 'model fraction' with boxes")
            plot1.send(" replot '" + logfiles[-1] + "' using 1:4 title 'mirror fraction' with boxes")

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

    import sys
    WorkingDir=sys.argv[1]

    # create ML2D_class object
    visualize_ML2D=visualize_ML2D_class(DoMatrixAllIter,
                                        DoShowLastIter,
                                        DoShowFractions,
                                        DoShowStatsAllIter,
                                        WorkingDir)
    # close 
    visualize_ML2D.close()

