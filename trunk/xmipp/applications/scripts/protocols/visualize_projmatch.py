#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_projmatch.py
#
# Example use:
# python visualize_projmatch.py
#
# This script requires that protocol_projmatch.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
#show results for iteration
DisplayIterNo=1
#------------------------------------------------------------------------------------------------
# {section} Volume visualization
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
# {expert} Width of Matrix-views (multiple of 3!):
MatrixWidth=9
#------------------------------------------------------------------------------------------------
# {section}Reference volume 
#------------------------------------------------------------------------------------------------
#show reference volume 
DisplayReference=False
#Show projection maching library and aligned classes
DisplayProjectionMatching=False
#display angular distribution after projection matching
DisplayAngularDistribution=False
#display angular distribution after  align2d
DisplayAngularDistributionAlign2d=False
#display reconstructed volume
DisplayReconstruction=False



#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_projmatch_class:

    #init variables
    def __init__(self,
                _VisualizeVolZ,
                _VisualizeVolX,
                _VisualizeVolY,
                _VisualizeVolChimera,
                _DisplayIterNo,
                _DisplayReference,
                _DisplayProjectionMatching,
                _DisplayReconstruction,
                _MatrixWidth,
                _DisplayAngularDistribution,
                _DisplayAngularDistributionAlign2d,
                _WorkingDir,
                ):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log,logging
        import visualization

        # import corresponding protocol
        os.chdir(_WorkingDir)
        import protocol_projmatch_backup
        os.chdir(os.pardir)

        self._MatrixWidth=_MatrixWidth
        self._WorkDirectory=protocol_projmatch_backup.WorkDirectory
        self._LogDir=protocol_projmatch_backup.LogDir
        self._ProjectDir=protocol_projmatch_backup.ProjectDir
        self._ReferenceVolume=protocol_projmatch_backup.ReferenceVolume
        self._multi_align2d_sel=protocol_projmatch_backup.multi_align2d_sel
        self._align2d_sel=protocol_projmatch_backup.align2d_sel
        self._align2d_doc=protocol_projmatch_backup.align2d_doc
        self._SelFileName=self._ProjectDir+'/'+str(protocol_projmatch_backup.SelFileName)
        self._ReferenceVolume=protocol_projmatch_backup.ReferenceVolume
        self._Proj_Maching_Output_Root_Name=protocol_projmatch_backup.Proj_Maching_Output_Root_Name
        self._Proj_Maching_Output_Root_Name + '.doc'
        self._mylog=log.init_log_system(self._ProjectDir,
                                       self._LogDir,
                                       sys.argv[0],
                                       self._WorkDirectory)
        self._mylog.setLevel(logging.DEBUG)
        self._iteration_number=_DisplayIterNo
        self._Iteration_Working_Directory=self._WorkDirectory+'/Iter_'+\
                                      str(self._iteration_number)
        #os.chdir(Iteration_Working_Directory)

        if (_DisplayReference):
           self.ShowVolumes=[] 
           self.ShowVolumes.append(os.getcwd()+'/'+\
                                   self._Iteration_Working_Directory+'/'+\
                                   self._ReferenceVolume)
           visualization.visualize_volumes(self.ShowVolumes,
                                                _VisualizeVolZ,
                                                _VisualizeVolX,
                                                _VisualizeVolY,
                                                _VisualizeVolChimera)
        if (_DisplayAngularDistribution):
           self._ShowPlots=[]
           self._ShowPlots.append(os.getcwd()+'/'+\
                                  self._Iteration_Working_Directory+'/'+\
                                  self._Proj_Maching_Output_Root_Name+\
                                  '_classes.doc')
           title='Angular distribution after "projection matching" for iteration '+\
                    str(self._iteration_number)
           show_ang_distribution(self._ShowPlots,self._iteration_number,title)
           self._mylog.debug( self._ShowPlots[0] + " " +\
                              str(self._iteration_number) + " " +\
                              title )

        if (_DisplayProjectionMatching):
           self.ShowSelfiles=[] 
           self.ShowSelfiles.append(os.getcwd()+'/'+\
                                    self._Iteration_Working_Directory+'/'+\
                                    self._multi_align2d_sel)
           currdir=os.getcwd()
           for selfiles in self.ShowSelfiles:
              os.chdir(self._Iteration_Working_Directory)
              visualization.visualize_images(self.ShowSelfiles,
                                             True,
                                             self._MatrixWidth,
                                             self._MatrixWidth,
                                             True)
           os.chdir(currdir)

        if (_DisplayAngularDistributionAlign2d):
           self._ShowPlots=[]
           self._ShowPlots.append(os.getcwd()+'/'+\
                                  self._Iteration_Working_Directory+'/'+\
                                  self._align2d_doc)
           title='Angular distribution after "align2d" for iteration '+\
                    str(self._iteration_number)
           show_ang_distribution(self._ShowPlots,self._iteration_number,title)
           self._mylog.debug( self._ShowPlots[0] + " " +\
                              str(self._iteration_number) + " " +\
                              title )

        #if (_DisplayReconstruction):
        #    self.visualize_Reconstruction(self._SomName,self._SpectraName)
            
        # Return to parent dir
        # os.chdir(os.pardir)

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

def show_ang_distribution(_ShowPlots,_iteration_number,title):
        import os
        import docfiles
        import visualization
        for plots in _ShowPlots: 
            doc=docfiles.docfile(plots)
            doc.check_angle_range()
            mini=doc.minimum_of_column(7)
            maxi=doc.maximum_of_column(7)
            if mini==0:
               mini=1
            if maxi<mini:
               maxi=mini   
            doc.write_several(plots,
                              10,
                              7,
                              mini,
                              maxi,
                              )
            plot=visualization.gnuplot()
            title =title+' min= '+str(mini)+', max= '+str(maxi) 
            plot.plot_xy1y2_several_angular_doc_files(plots,
                                                      title,
                                                      'degrees',
                                                      'degrees',
                                                      3,
                                                      4)


#		
# Main
#     
if __name__ == '__main__':

    import sys
    WorkingDir=sys.argv[1]

    # create projmatch_class object
    visualize_projmatch=visualize_projmatch_class(VisualizeVolZ,
                                                  VisualizeVolX,
                                                  VisualizeVolY,
                                                  VisualizeVolChimera,
                                                  DisplayIterNo,
                                                  DisplayReference,
                                                  DisplayProjectionMatching,
                                                  DisplayReconstruction,
                                                  MatrixWidth,
                                                  DisplayAngularDistribution,
                                                  DisplayAngularDistributionAlign2d,
                                                  WorkingDir)
    # close 
    visualize_projmatch.close()

