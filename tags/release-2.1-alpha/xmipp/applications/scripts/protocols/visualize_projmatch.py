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
DisplayIterationsNo='1 2 3 4'
#------------------------------------------------------------------------------------------------
# {section} Display 3D volumes
#------------------------------------------------------------------------------------------------
#display reference volume
""" Volume after filtration and masking
"""
DisplayReference=False
#display reconstructed volume
""" Volume as given by the reconstruction algorithm
"""
DisplayReconstruction=True
#display reconstructed volume after filtration
""" Volume after filtration
"""
DisplayFilteredReconstruction=False
#------------------------------------------------------------------------------------------------
# {section} Display 2D projection galleries
#------------------------------------------------------------------------------------------------
#Show projection matching library and aligned classes
DisplayProjectionMatchingAlign2d=True
#------------------------------------------------------------------------------------------------
# {section} Display 2D-plots
#------------------------------------------------------------------------------------------------
#display angular distribution after projection matching
DisplayAngularDistribution=True
#display resolution plots (FSC)
DisplayResolutionPlots=True
#------------------------------------------------------------------------------------------------
# {section} Display options
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
# {expert} Width of projection galleries
""" In number of reference projections. This number will be multiplied by 3 if re-alignment of 
    classes was performed and multiplied by 2 otherwise.
"""
MatrixWidth=3
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
                _DisplayProjectionMatchingAlign2d,
                _DisplayReconstruction,
                _DisplayFilteredReconstruction,
                _DisplayResolutionPlots,
                _MatrixWidth,
                _DisplayAngularDistribution,
                _ProtocolName
                ):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,logging,arg
        import visualization

        # Import the corresponding protocol, get WorkingDir and go there
        pardir=os.path.abspath(os.getcwd())
        shutil.copy(_ProtocolName,'protocol.py')
        import protocol

        self._MatrixWidth=_MatrixWidth
        self._WorkDirectory=protocol.WorkDirectory
        self._LogDir=protocol.LogDir
        self._ProjectDir=protocol.ProjectDir
        self._multi_align2d_sel=protocol.MultiAlign2dSel
        self._SelFileName=self._ProjectDir+'/'+str(protocol.SelFileName)
        self._Filtered_Image=protocol.FilteredReconstruction
        self._mylog=log.init_log_system(self._ProjectDir,
                                       self._LogDir,
                                       sys.argv[0],
                                       self._WorkDirectory)
        self._mylog.setLevel(logging.DEBUG)
        self._DisplayIterationsNo=arg.getListFromVector(_DisplayIterationsNo)
        self._ReconstrucedVolume=protocol.ReconstructedVolume
        _user_suplied_Filtered_Image=os.path.abspath(
                                 protocol.ReferenceVolumeName)
        protocol.NumberofIterations += 1                         
        self._Reference_volume=[]
        protocol.fill_name_vector(
                        "",
                        self._Reference_volume,
                        protocol.NumberofIterations,
                        protocol.ReferenceVolumeName)
        self._mylog.debug ( "_Reference_volume "+str(self._Reference_volume))
        self._Filtered_Image=[]
        protocol.fill_name_vector(
                          "",
                          self._Filtered_Image,
                          protocol.NumberofIterations,
                          protocol.FilteredReconstruction)
        self._mylog.debug ( "_Filtered_Image "+str(self._Filtered_Image))
        
        self._DisplayReference_list=[]
        self._ShowPlotsList=[]
        self._TitleList=[]
        self._ShowSelfilesProjMachingAlign2d=[] 
        self._ShowSelfilesProjMaching=[]
        self._DisplayReconstruction_list=[]
        self._DisplayResolutionPlots_list=[]
        self._TitleResolution=[]
        self._DisplayFilteredReconstruction_list=[]
        
        self._mylog.debug ("_DisplayIterationsNo " + _DisplayIterationsNo)
        for self._iteration_number in self._DisplayIterationsNo:
           if self._iteration_number==' ':
              continue

           if (_DisplayReference):
              self._DisplayReference_list.append(
                           self._Reference_volume[int(self._iteration_number)])
           
           if (_DisplayReconstruction):
              self._DisplayReconstruction_list.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+\
                                      '/Iter_'+str(self._iteration_number)+\
                                      '_'+self._ReconstrucedVolume+'.vol')

           if (_DisplayFilteredReconstruction):
              self._DisplayFilteredReconstruction_list.append(
                          self._Filtered_Image[int(self._iteration_number)]+'.vol')


           if (_DisplayProjectionMatchingAlign2d):
               self._doAlign2d = arg.getComponentFromVector(protocol.DoAlign2D,\
                                                                int(self._iteration_number)-1)
               if (self._doAlign2d == 1):
                   self._ShowSelfilesProjMachingAlign2d.append('..'+'/Iter_'+\
                                                                   str(self._iteration_number)+\
                                                                   '/'+ self._multi_align2d_sel)
               else:
                   self._ShowSelfilesProjMaching.append('..'+'/Iter_'+\
                                                            str(self._iteration_number)+\
                                                            '/'+ self._multi_align2d_sel)
    
           if (_DisplayAngularDistribution):
              self._ShowPlotsList.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/'+\
                                     protocol.ForReconstructionDoc)
              self._TitleList.append('Angular distribution after *projection matching* for iteration '+\
                       str(self._iteration_number))

           if (_DisplayResolutionPlots):
              self._DisplayResolutionPlots_list.append('..'+'/Iter_'+\
                                         str(self._iteration_number)+\
                                         '/Iter_'+str(self._iteration_number)+\
                                         '_resolution.fsc')
              self._TitleResolution.append('Fourier shell correlation for iteration '+\
                       str(self._iteration_number))


        #NAMES ARE RELATIVE TO ITER DIRECTORY  
        self._Iteration_Working_Directory=self._WorkDirectory+'/Iter_1'
        self._mylog.debug ("cd " + self._Iteration_Working_Directory)
        os.chdir(self._Iteration_Working_Directory)

        if (_DisplayReference):
           self._mylog.debug ( "_DisplayReference_list "+str(self._DisplayReference_list))
           visualization.visualize_volumes(self._DisplayReference_list,
                                           _VisualizeVolZ,
                                           _VisualizeVolX,
                                           _VisualizeVolY,
                                           _VisualizeVolChimera)


        if (_DisplayReconstruction):
           self._mylog.debug ( "_DisplayReconstruction_list "+str(self._DisplayReconstruction_list))
           visualization.visualize_volumes(self._DisplayReconstruction_list,
                                           _VisualizeVolZ,
                                           _VisualizeVolX,
                                           _VisualizeVolY,
                                           _VisualizeVolChimera)

        if (_DisplayFilteredReconstruction):
           self._mylog.debug ( "_DisplayFilteredReconstruction_list "+\
                            str(self._DisplayFilteredReconstruction_list))
           visualization.visualize_volumes(self._DisplayFilteredReconstruction_list,
                                           _VisualizeVolZ,
                                           _VisualizeVolX,
                                           _VisualizeVolY,
                                           _VisualizeVolChimera)

        if (_DisplayProjectionMatchingAlign2d):
            if (protocol.CleanUpFiles==True):
                print '* You cannot visualize class averages after performing a cleanup! '
            else:
                self._mylog.debug ( "_ShowSelfilesProjMachingAlign2d "+\
                                        str(self._ShowSelfilesProjMachingAlign2d))
                visualization.visualize_images(self._ShowSelfilesProjMachingAlign2d,
                                               True,
                                               3*self._MatrixWidth,
                                               "",
                                               True)
                visualization.visualize_images(self._ShowSelfilesProjMaching,
                                               True,
                                               2*self._MatrixWidth,
                                               "",
                                               False)

        if (_DisplayAngularDistribution):
           self._mylog.debug ( "_ShowPlotsList "+str(self._ShowPlotsList))
           self._mylog.debug ( "_TitleList "+str(self._TitleList))
           show_ang_distribution(self._ShowPlotsList,
                                 self._iteration_number,
                                 self._TitleList,
                                 self._mylog)

        if (_DisplayResolutionPlots):
           self._mylog.debug ( "_DisplayResolutionPlots_list "+str(self._DisplayResolutionPlots_list))
           show_plots(self._DisplayResolutionPlots_list,
                      self._iteration_number,
                      self._TitleResolution,
                      self._mylog)



    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

def show_ang_distribution(_ShowPlots,_iteration_number,_title,_mylog=""):
        import os
        import docfiles
        import visualization
        for i in range(len(_ShowPlots)):
#        for plots in _ShowPlots: 
            doc=docfiles.docfile(_ShowPlots[i])
            doc.check_angle_range()
            mini=doc.minimum_of_column(7)
            maxi=doc.maximum_of_column(7)
	    if not _mylog=="":
               _mylog.debug("mini "+ str(mini) +" maxi "+ str(maxi))
            if mini==0:
               mini=1
            if maxi<mini:
               maxi=mini   
            doc.write_several(_ShowPlots[i],
                              10,
                              7,
                              mini,
                              maxi,
                              )
            plot=visualization.gnuplot()
            _title[i] =_title[i]+'\\n min= '+str(mini)+', max= '+str(maxi) 
            plot.plot_xy1y2_several_angular_doc_files(_ShowPlots[i],
                                                      _title[i],
                                                      'degrees',
                                                      'degrees',
                                                      3,
                                                      4)
def show_plots(_ShowPlots,_iteration_number,_title,_mylog):
        import os
        import docfiles
        import visualization
        plot=visualization.gnuplot()
        for i in range(len(_ShowPlots)):
#        for plots in _ShowPlots:
           plot.prepare_empty_plot(_title[i],
                                   "Armstrong^-1",
                                   "Fourier Shell Correlation")
           X_col=1
           Y_col=2
           if (i==0):
              plot.send(" plot '" + _ShowPlots[i] + "' using "+str(X_col)+":"+str(Y_col)+" with lines")   
           else:
              plot.send(" replot '" + _ShowPlots[i] + "' using "+str(X_col)+":"+str(Y_col)+" with lines")   

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
                                                  DisplayProjectionMatchingAlign2d,
                                                  DisplayReconstruction,
                                                  DisplayFilteredReconstruction,
                                                  DisplayResolutionPlots,
                                                  MatrixWidth,
                                                  DisplayAngularDistribution,
                                                  ProtocolName)
    # close 
    visualize_projmatch.close()

