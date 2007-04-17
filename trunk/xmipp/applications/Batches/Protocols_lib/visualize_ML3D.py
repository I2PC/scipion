#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_ML3D.py
#
# Example use:
# python visualize_ML3D.py
#
# This script requires that protocol_ML3D.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Visualize volumes in slices along Z?
VisualizeVolZ=False
# Visualize volumes in slices along X?
VisualizeVolX=False
# Visualize volumes in slices along Y?
VisualizeVolY=False
# Visualize volumes in UCSF Chimera?
""" For this to work, you need to have chimera installed!
"""
VisualizeVolChimera=True
# {expert} Width of Matrix-views (even value!):
MatrixWidth=10
#------------------------------------------------------------------------------------------------
# {section} Visualization of seeds preparation steps
#------------------------------------------------------------------------------------------------
# Visualize the library projections and the averages of the grey-scale correction?
DoVisualizeMatrixCorrectReference=True
# Visualize the grey-scale corrected reference volume?
DoVisualizeCorrectReference=True
# Visualize the low-pass filtered reference volume?
DoVisualizeFilteredReference=True
# Visualize the library projections and the averages of the seed generation runs?
DoVisualizeMatrixSeeds=True
# Visualize the generated seeds?
DoVisualizeGeneratedSeeds=True
#------------------------------------------------------------------------------------------------
# {section} Visualization of ML3D iterations
#------------------------------------------------------------------------------------------------
# Which iterations to visualize? (Separate numbers by comma's)
SelectIterations="2"
# Visualize the reference volumes for the given iterations?
VisualizeML3DReferences=True
# Visualize weighted 2D-averages?
VisualizeML3DAvgs=True
# Output to screen the number of particles that change optimal orientation/class?
""" As described in Scheres et al. (2007) Nature Methods, 4, 27-29 
"""
DoShowMovingParticles=True
# Output overall statistics of all iterations to screen?
DoShowStatsAllIter=True
# Visualize matrixview of library projections and weighted averages of the last iteration?
VisualizeMatrixLastIter=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
import protocol_ML3D
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
                 DoShowMovingParticles,
                 DoShowStatsAllIter,
                 VisualizeMatrixLastIter
                 ):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log
        import visualization


        self.DoVisualizeMatrixCorrectReference=DoVisualizeMatrixCorrectReference
        self.DoVisualizeCorrectReference=DoVisualizeCorrectReference
        self.DoVisualizeFilteredReference=DoVisualizeFilteredReference
        self.DoVisualizeGeneratedSeeds=DoVisualizeGeneratedSeeds
        self.DoVisualizeMatrixSeeds=DoVisualizeMatrixSeeds

        self.VisualizeML3DReferences=VisualizeML3DReferences
        self.VisualizeML3DAvgs=VisualizeML3DAvgs
        self.SelectIterations=SelectIterations
        self.VisualizeMatrixLastIter=VisualizeMatrixLastIter;
        self.DoShowMovingParticles=DoShowMovingParticles
        self.DoShowStatsAllIter=DoShowStatsAllIter
        self.ShowVolumes=[]
        self.ShowSelfiles=[]

        os.chdir(protocol_ML3D.WorkingDir)

        self.visualize_preparation()

        self.visualize_seed_generation()

        self.visualize_ML3D_iterations()
        
        self.visualize_ML3D_overall()

        # Do the actual visualization of volumes and selfiles
        visualization.visualize_volumes(self.ShowVolumes,VisualizeVolZ,VisualizeVolX,VisualizeVolY,VisualizeVolChimera)
        visualization.visualize_images(self.ShowSelfiles,True,MatrixWidth)

        # Return to parent dir
        os.chdir(os.pardir)

    def visualize_preparation(self):
        if (self.DoVisualizeMatrixCorrectReference):
            dirname='CorrectGreyscale/'
            selfile1=dirname+'corrected_reference_lib.sel'
            selfile2=dirname+'corrected_reference_classes.sel'
            self.join_selfiles(selfile1,selfile2,'matrixview_greyscale_correction.sel')
            self.ShowSelfiles.append('matrixview_greyscale_correction.sel')
        if (self.DoVisualizeCorrectReference):
            self.ShowVolumes.append('corrected_reference.vol')
        if (self.DoVisualizeFilteredReference):
            self.ShowVolumes.append('filtered_reference.vol')

    def visualize_seed_generation(self):
        for i in range(protocol_ML3D.NumberOfReferences):
            dirname='GenerateSeed_'+str(i+1)+'/'
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


    def visualize_ML3D_iterations(self):
        import os,glob
        iters=self.SelectIterations.split(',')
        for iter in iters:
            if (self.VisualizeML3DAvgs):
                selfiles=glob.glob('RunML3D/ml3d_it'+str(iter).zfill(5)+'_vol?????.sel')
                for selfile in selfiles:
                    self.ShowSelfiles.append(selfile)

            if (self.VisualizeML3DReferences):
                volumes=glob.glob('RunML3D/ml3d_it'+str(iter).zfill(5)+'_vol?????.vol')
                for volume in volumes:
                    self.ShowVolumes.append(volume)

    def visualize_ML3D_overall(self):

        if (self.DoShowMovingParticles):
            self.show_moving_particles()
        
        if (self.DoShowStatsAllIter):
            self.show_statistics_alliter()

        if (self.VisualizeMatrixLastIter):
            self.visualize_matrix_last_iter()
        
    def visualize_matrix_last_iter(self):
        import os,glob
        selfiles=glob.glob('RunML3D/ml3d_it?????.sel')
        if len(selfiles)==0:
            print "No selfiles yet. Visualize only after job completion..."
        else:
            lastselfile=selfiles[-1]
        selfile1='RunML3D/ml3d_lib.sel'
        selfile2=lastselfile
        self.join_selfiles(selfile1,selfile2,'matrixview_ML3D_lastiter.sel')
        self.ShowSelfiles.append('matrixview_ML3D_lastiter.sel')
            

    def show_statistics_alliter(self):
        import os,glob
        logfiles=glob.glob('RunML3D/ml3d_it?????.log')
        if len(logfiles)==0:
            print "No logfiles yet. Visualize only after job completion..."
        else:
            for logfile in logfiles:
                fh=open(logfile,'r')
                line1=fh.readline()
                line2=fh.readline()
                fh.close()
                words1=line1.split()
                words2=line2.split()
                print logfile,words1[6],words1[7],words1[8],words1[9],\
                      words2[1],words2[2],words2[3],words2[4]

    def show_moving_particles(self):
        import os,glob
        docfiles=glob.glob('RunML3D/ml3d_it?????.doc')
        if os.path.exists('moving_particles.dat'):
            os.remove('moving_particles.dat')
        print 'The output of moving particles remains to be added to this script...'

    def join_selfiles(self,selfile1,selfile2,newselfile):
        import SelFiles
        sel1=SelFiles.selfile()
        sel2=SelFiles.selfile()
        sel1.read(selfile1)
        sel2.read(selfile2)
        newsel=sel1.intercalate_union(sel2)
        newsel.write(newselfile)

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
                                        DoShowMovingParticles,
                                        DoShowStatsAllIter,
                                        VisualizeMatrixLastIter
                                        )
    # close 
    visualize_ML3D.close()

