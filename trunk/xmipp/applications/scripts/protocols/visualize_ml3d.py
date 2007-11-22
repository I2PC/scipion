#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_ML3D.py
#
# Example use:
# python visualize_ml3d.py protocol_ml3d.py
#
# This script requires that protocol_ml3d.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Visualize volumes in slices along Z?
VisualizeVolZ=True
# {expert} Visualize volumes in slices along X?
VisualizeVolX=False
# {expert} Visualize volumes in slices along Y?
VisualizeVolY=False
# Visualize volumes in UCSF Chimera?
""" For this to work, you need to have chimera installed!
"""
VisualizeVolChimera=False
# {expert} Width of Matrix-views (even value!):
MatrixWidth=10
#------------------------------------------------------------------------------------------------
# {section} Visualization of seeds preparation steps
#------------------------------------------------------------------------------------------------
# Visualize the library projections and the averages of the grey-scale correction?
DoVisualizeMatrixCorrectReference=False
# Visualize the grey-scale corrected reference volume?
DoVisualizeCorrectReference=False
# Visualize the low-pass filtered reference volume?
DoVisualizeFilteredReference=False
# Visualize the library projections and the averages of the seed generation runs?
DoVisualizeMatrixSeeds=False
# Visualize the generated seed volumes?
DoVisualizeGeneratedSeeds=False
#------------------------------------------------------------------------------------------------
# {section} Per-iteration ML3D Visualization
#------------------------------------------------------------------------------------------------
# Which iterations to visualize? (Separate numbers by comma's)
SelectIterations="1,2"
# Visualize the reference volumes for the given iterations?
VisualizeML3DReferences=True
# Visualize weighted 2D-averages?
VisualizeML3DAvgs=True
# Plot the angular distribution of the reference(s)?
VisualizeAngDistribution=True
# Plot data distribution over the distinct references?
VisualizeClassDistribution=True
#------------------------------------------------------------------------------------------------
# {section} Overall ML3D visualization
#------------------------------------------------------------------------------------------------
# Plot overall convergence statistics?
""" As also described in Scheres et al. (2007) Nature Methods, 4, 27-29 
"""
DoShowStatsAllIter=True
# Visualize matrixview of library projections and weighted averages of the last iteration?
VisualizeMatrixLastIter=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_ML3D_class:

    #init variables
    def __init__(self,
                 VisualizeVolZ,
                 VisualizeVolX,
                 VisualizeVolY,
                 VisualizeVolChimera,
                 MatrixWidth,
                 DoVisualizeMatrixCorrectReference,
                 DoVisualizeMatrixSeeds,
                 DoVisualizeGeneratedSeeds,
                 SelectIterations,
                 VisualizeML3DReferences,
                 VisualizeML3DAvgs,
                 VisualizeAngDistribution,
                 VisualizeClassDistribution,
                 DoShowStatsAllIter,
                 VisualizeMatrixLastIter,
                 ProtocolName
                 ):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import visualization


        self.DoVisualizeMatrixCorrectReference=DoVisualizeMatrixCorrectReference
        self.DoVisualizeCorrectReference=DoVisualizeCorrectReference
        self.DoVisualizeFilteredReference=DoVisualizeFilteredReference
        self.DoVisualizeGeneratedSeeds=DoVisualizeGeneratedSeeds
        self.DoVisualizeMatrixSeeds=DoVisualizeMatrixSeeds

        self.VisualizeML3DReferences=VisualizeML3DReferences
        self.VisualizeML3DAvgs=VisualizeML3DAvgs
        self.VisualizeAngDistribution=VisualizeAngDistribution
        self.VisualizeClassDistribution=VisualizeClassDistribution
        self.SelectIterations=SelectIterations
        self.VisualizeMatrixLastIter=VisualizeMatrixLastIter;
        self.DoShowStatsAllIter=DoShowStatsAllIter
        self.ShowVolumes=[]
        self.ShowSelfiles=[]

        # Import the corresponding protocol, get WorkingDir and go there
        pardir=os.path.abspath(os.getcwd())
        shutil.copy(ProtocolName,'protocol.py')
        import protocol
        self.WorkingDir=protocol.WorkingDir
        os.chdir(self.WorkingDir)

        self.visualize_preparation()

        self.visualize_seed_generation()

        self.visualize_ML3D()
        
        # Do the actual visualization of volumes and selfiles
        visualization.visualize_volumes(self.ShowVolumes,VisualizeVolZ,VisualizeVolX,VisualizeVolY,VisualizeVolChimera)
        visualization.visualize_images(self.ShowSelfiles,True,MatrixWidth)

        # Return to parent dir and remove protocol.py(c)
        os.chdir(pardir)
        if (os.path.exists('protocol.py')):
            os.remove('protocol.py')
        if (os.path.exists('protocol.pyc')):
            os.remove('protocol.pyc')

    def visualize_preparation(self):
        import os
        if (self.DoVisualizeMatrixCorrectReference):
            dirname='CorrectGreyscale/'
            selfile1=dirname+'corrected_reference_lib.sel'
            selfile2=dirname+'corrected_reference_classes.sel'
            if not os.path.exists(dirname):
                message='Warning: directory '+dirname+' does not exist, skipping...'
                print '* ',message
            else:
                self.join_selfiles(selfile1,selfile2,'matrixview_greyscale_correction.sel')
                self.ShowSelfiles.append('matrixview_greyscale_correction.sel')
        if (self.DoVisualizeCorrectReference):
            self.ShowVolumes.append('corrected_reference.vol')
        if (self.DoVisualizeFilteredReference):
            self.ShowVolumes.append('filtered_reference.vol')

    def visualize_seed_generation(self):
        import glob
        seeddirs=glob.glob('GenerateSeed_*/')
        seeddirs.sort()
        for i in range(len(seeddirs)):
            dirname=seeddirs[i]
            outname='seeds_split_'+str(i+1)
            if (self.DoVisualizeGeneratedSeeds):
                seedname=dirname+outname+'_it00001.vol'
                self.ShowVolumes.append(seedname)
            if (self.DoVisualizeMatrixSeeds):
                selfile1=dirname+outname+'_lib.sel'
                selfile2=dirname+outname+'_it00001.sel'
                matrixname='matrixview_'+outname+'.sel'
                self.join_selfiles(selfile1,selfile2,matrixname)
                self.ShowSelfiles.append(matrixname)


    def visualize_ML3D(self):
        import os,glob
        iters=self.SelectIterations.split(',')
        for iter in iters:
            selfiles=glob.glob('RunML3D/ml3d_it'+str(iter).zfill(5)+'_vol?????.sel')
            selfiles.sort()
            if (self.VisualizeML3DAvgs):
                for selfile in selfiles:
                    self.ShowSelfiles.append(selfile)

            if (self.VisualizeAngDistribution or self.VisualizeClassDistribution):
                self.prepare_distribution_docfiles(selfiles)
                
                if (self.VisualizeAngDistribution):
                    for i in range(len(selfiles)): 
                        self.show_ang_distribution(selfiles[i],iter,i+1)

                if (self.VisualizeClassDistribution):
                    self.show_class_distribution(selfiles,iter)

            if (self.VisualizeML3DReferences):
                volumes=glob.glob('RunML3D/ml3d_it'+str(iter).zfill(5)+'_vol?????.vol')
                volumes.sort()
                for volume in volumes:
                    self.ShowVolumes.append(volume)

        if (self.VisualizeMatrixLastIter):
            self.visualize_matrix_last_iter()

        if (self.DoShowStatsAllIter):
            self.show_convergence_stats()

    def visualize_matrix_last_iter(self):
        import os,glob
        selfiles=glob.glob('RunML3D/ml3d_it?????.sel')
        selfiles.sort()
        if len(selfiles)==0:
            print "No selfiles yet. Visualize only after job completion..."
        else:
            lastselfile=selfiles[-1]
        selfile1='RunML3D/ml3d_lib.sel'
        selfile2=lastselfile
        self.join_selfiles(selfile1,selfile2,'matrixview_ML3D_lastiter.sel')
        self.ShowSelfiles.append('matrixview_ML3D_lastiter.sel')
            
    def prepare_distribution_docfiles(self,selfiles):
        import os
        for i in range(len(selfiles)):
            docname='ang_distribution_ref'+str(i+1).zfill(5)+'.doc'
            command='xmipp_header_extract -i '+selfiles[i]+' -o '+docname
            print '* ',command
            os.system(command)
 
    def show_class_distribution(self,selfiles,iter):
        import os
        import docfiles
        import visualization

        newlines=[]
        newlines.append('reference | sum weights (# images) \n')
        for i in range(len(selfiles)):
            docname='ang_distribution_ref'+str(i+1).zfill(5)+'.doc'
            doc=docfiles.docfile(docname)
            newlines.append(str(i+1)+' '+str(doc.sum_of_column(7))+'\n')

        fh=open('class_distribution.dat','w')
        fh.writelines(newlines)
        fh.close()
        nr_class=len(newlines)
        nr_class-=0.5
        plot=visualization.gnuplot()
        plot.prepare_empty_plot('Number of images contributing to each class for iteration '+str(iter),
                                'class number',
                                'sumweight')
        plot.send(" set xrange [0.5:"+str(nr_class)+"]")
        plot.send(" set xtics 1")
        plot.send(" plot 'class_distribution.dat' using 1:2 with boxes")

    def show_ang_distribution(self,selfile,iter,ref):
        import os
        import docfiles
        import visualization

        docname='ang_distribution_ref'+str(ref).zfill(5)
        doc=docfiles.docfile(docname+'.doc')
        newname=docname+'_bin_'
        doc.check_angle_range()
        doc.write_several(newname,
                          10,
                          7,
                          doc.minimum_of_column(7),
                          doc.maximum_of_column(7)
                          )
        plot=visualization.gnuplot()
        title='Angular distribution for reference '+str(ref)+' in iteration '+str(iter)
        plot.plot_xy1y2_several_angular_doc_files(newname,
                                                  title,
                                                  'degrees',
                                                  'degrees')

    def show_convergence_stats(self):
        import os,glob
        import visualization
        logfiles=glob.glob('RunML3D/ml3d_it?????.log')
        logfiles.sort()
        if len(logfiles)==0:
            print "No logfiles yet. Visualize after job completion..."
        else: 
            newlines1=[]
            newlines2=[]
            newlines3=[]
            newlines1.append('# iter | log-likelihood \n')
            newlines2.append('# iter | Pmax/sumP \n')
            newlines3.append('# iter | moving particles \n')
            for i in range(len(logfiles)):
                fh=open(logfiles[i],'r')
                line1=fh.readline()
                line2=fh.readline()
                fh.close()
                words1=line1.split()
                words2=line2.split()
                # get iteration
                iter=int(words2[6])-1
                # get data
                newlines1.append(str(iter)+' '+words1[7]+'\n')
                newlines2.append(str(iter)+' '+words1[9]+'\n')
                # moving particles
                if (iter > 1):
                    docfile1=logfiles[i].replace('log','doc')
                    docfile0=logfiles[i-1].replace('log','doc')
                    mov=self.moving_particles(docfile0,docfile1)
                    newlines3.append(str(iter)+' '+str(mov)+'\n')
            fh=open('alliter_LL.dat','w')
            fh.writelines(newlines1)
            fh.close();
            fh=open('alliter_Pmax.dat','w')
            fh.writelines(newlines2)
            fh.close();
            fh=open('alliter_moving_particles.dat','w')
            fh.writelines(newlines3)
            fh.close();
            plot1=visualization.gnuplot()
            plot1.plot_xy_file('alliter_LL.dat',
                              'Log-likelihood (should increase)',
                              'iterations',
                              'LL')
            plot2=visualization.gnuplot()
            plot2.plot_xy_file('alliter_Pmax.dat',
                              'The width of the probability distributions (ideally goes to one i.e. delta-functions)',
                              'iterations',
                              'Pmax/sumP')
            plot3=visualization.gnuplot()
            plot3.plot_xy_file('alliter_moving_particles.dat',
                              'The number of particles that change their optimal orientation/class (ideally goes to zero)',
                              'iterations',
                              '#')
                             
    def moving_particles(self,docfile0,docfile1):
        import os, sys
        fh=open(docfile0,'r')
        lines0=fh.readlines()
        fh.close()
        fh=open(docfile1,'r')
        lines1=fh.readlines()
        fh.close()

        if not (len(lines0) == len(lines1)):
            message='Error: docfiles '+docfile0+' and '+docfile1+' are of unequal lenght'
            print '* ',message
            sys.exit()

        count=0
        for i in range(len(lines0)):
            if (lines0[i].find(';')<0):
                words0=lines0[i].split()
                words1=lines1[i].split()
                for entry in [2, 3, 4, 5, 6, 7, 8]:
                    if not (words0[entry]==words1[entry]):
                        count = count + 1
                        break
        return count

    def join_selfiles(self,selfile1,selfile2,newselfile):
        import selfile
        sel1=selfile.selfile()
        sel2=selfile.selfile()
        sel1.read(selfile1)
        sel2.read(selfile2)
        newsel=sel1.intercalate_union(sel2)
        newsel.write(newselfile)

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

    # create ML3D_class object
    visualize_ML3D=visualize_ML3D_class(VisualizeVolZ,
                                        VisualizeVolX,
                                        VisualizeVolY,
                                        VisualizeVolChimera,
                                        MatrixWidth,
                                        DoVisualizeMatrixCorrectReference,
                                        DoVisualizeMatrixSeeds,
                                        DoVisualizeGeneratedSeeds,
                                        SelectIterations,
                                        VisualizeML3DReferences,
                                        VisualizeML3DAvgs,
                                        VisualizeAngDistribution,
                                        VisualizeClassDistribution,
                                        DoShowStatsAllIter,
                                        VisualizeMatrixLastIter,
                                        ProtocolName
                                        )
    # close 
    visualize_ML3D.close()

