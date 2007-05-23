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
#show results for iterations
""" You may specify more than one iteration here 
    This can be done by a sequence of numbers (for instance, "2 8" 
    specifies iteration 2 and 8 (but not 3, 4, 5, 6 and 7)
"""
DisplayIterationsNo="1 3"
#------------------------------------------------------------------------------------------------
# {section} Volume visualization
#------------------------------------------------------------------------------------------------
# Visualize volumes in slices along Z?
VisualizeVolZ=True
# Visualize volumes in slices along X?
VisualizeVolX=False
# Visualize volumes in slices along Y?
VisualizeVolY=False
# Visualize volumes in UCSF Chimera?
""" For this to work, you need to have chimera installed!
"""
VisualizeVolChimera=True
# {expert} Width of xmipp_show (multiple of 3!):
MatrixWidth=9
#------------------------------------------------------------------------------------------------
# {section}Reference volume 
#------------------------------------------------------------------------------------------------
#show reference volume 
DisplayReference=False
#display angular distribution after projection matching
DisplayAngularDistribution=False
#Show projection maching library and aligned classes
DisplayProjectionMatching=False
#display angular distribution after align2d
DisplayAngularDistributionAlign2d=False
#display reconstructed volume
DisplayReconstruction=False
#display resolution plots
DisplayResolutionPlots=True


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
                _DisplayIterationsNo,
                _DisplayReference,
                _DisplayProjectionMatching,
                _DisplayReconstruction,
                _DisplayResolutionPlots,
                _MatrixWidth,
                _DisplayAngularDistribution,
                _DisplayAngularDistributionAlign2d,
                _ProtocolName,
                ):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log,logging,arg
        import visualization

        # Import the corresponding protocol, get WorkingDir and go there
        import protocol
        protocol.__name__=ProtocolName.replace('py','')
        self._WorkDirectory=protocol.WorkDirectory

        self._MatrixWidth=_MatrixWidth
        self._WorkDirectory=protocol.WorkDirectory
        self._LogDir=protocol.LogDir
        self._ProjectDir=protocol.ProjectDir
        self._multi_align2d_sel=protocol.multi_align2d_sel
        self._align2d_sel=protocol.align2d_sel
        self._align2d_doc=protocol.align2d_doc
        self._SelFileName=self._ProjectDir+'/'+str(protocol.SelFileName)
        self._ReferenceVolume=protocol.ReferenceVolume
        self._Proj_Maching_Output_Root_Name=protocol.Proj_Maching_Output_Root_Name
        self._Proj_Maching_Output_Root_Name + '.doc'
        self._mylog=log.init_log_system(self._ProjectDir,
                                       self._LogDir,
                                       sys.argv[0],
                                       self._WorkDirectory)
        self._mylog.setLevel(logging.DEBUG)
        self._DisplayIterationsNo=arg.getListFromVector(_DisplayIterationsNo)
        
        for self._iteration_number in _DisplayIterationsNo:
           if self._iteration_number==' ':
              continue
           self._Iteration_Working_Directory=_WorkingDir+'/Iter_'+\
                                         str(self._iteration_number)+'/'

           self._mylog.debug ("cd " + self._Iteration_Working_Directory)
           os.chdir(self._Iteration_Working_Directory)
       
           # how many iterations should a process


           if (_DisplayReference):
              self.ShowVolumes=[] 
              self.ShowVolumes.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/'+self._ReferenceVolume)
              visualization.visualize_volumes(self.ShowVolumes,
                                                   _VisualizeVolZ,
                                                   _VisualizeVolX,
                                                   _VisualizeVolY,
                                                   _VisualizeVolChimera)
           if (_DisplayAngularDistribution):
              self._ShowPlots=[]
              self._ShowPlots.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/'+\
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
              self.ShowSelfiles.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/'+\
                                       self._multi_align2d_sel)
              for selfiles in self.ShowSelfiles:
                 visualization.visualize_images(self.ShowSelfiles,
                                                True,
                                                self._MatrixWidth,
                                                self._MatrixWidth,
                                                True)

           if (_DisplayAngularDistributionAlign2d):
              self._ShowPlots=[]
              self._ShowPlots.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/'+\
                                     self._align2d_doc)
              title='Angular distribution after "align2d" for iteration '+\
                       str(self._iteration_number)
              show_ang_distribution(self._ShowPlots,self._iteration_number,title)
              self._mylog.debug( self._ShowPlots[0] + " " +\
                                 str(self._iteration_number) + " " +\
                                 title )

           if (_DisplayReconstruction):
              self.ShowVolumes=[] 
              self.ShowVolumes.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/reconstruction.vol')
              self._mylog.debug (self.ShowVolumes[0])
              visualization.visualize_volumes(self.ShowVolumes,
                                                   _VisualizeVolZ,
                                                   _VisualizeVolX,
                                                   _VisualizeVolY,
                                                   _VisualizeVolChimera)

           if (_DisplayResolutionPlots):
              plot=visualization.gnuplot()
              plot_name='..'+'/Iter_'+\
                         str(self._iteration_number)+\
                         '/'+'split_sel_2.vol.frc'
              plot.plot_xy_file(plot_name,
                          Title="Resolution",
                          X_Label="Armstrong^-1",
                          Y_Label="y",
                          X_col=1,
                          Y_col=2)


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
    ProtocolName=sys.argv[1]

    # create projmatch_class object
    visualize_projmatch=visualize_projmatch_class(VisualizeVolZ,
                                                  VisualizeVolX,
                                                  VisualizeVolY,
                                                  VisualizeVolChimera,
                                                  DisplayIterationsNo,
                                                  DisplayReference,
                                                  DisplayProjectionMatching,
                                                  DisplayReconstruction,
                                                  DisplayResolutionPlots,
                                                  MatrixWidth,
                                                  DisplayAngularDistribution,
                                                  DisplayAngularDistributionAlign2d,
                                                  ProtocolName)
    # close 
    visualize_projmatch.close()

