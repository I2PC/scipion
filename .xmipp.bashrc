export XMIPP_HOME=/usr/local/xmipp.3.0
export PATH=$XMIPP_HOME/bin:$PATH
export LD_LIBRARY_PATH=$XMIPP_HOME/lib:$LD_LIBRARY_PATH
# Load global autocomplete file 
test -s $XMIPP_HOME/.xmipp.autocomplete && . $XMIPP_HOME/.xmipp.autocomplete || true
test -s $XMIPP_HOME/.xmipp_programs.autocomplete && . $XMIPP_HOME/.xmipp_programs.autocomplete || true
 
# Xmipp Aliases 						 
## Setup ##                        
alias xconfigure='./setup.py -j 16 configure compile ' 
alias xcompile='./setup.py -j 16 compile ' 
alias xupdate='./setup.py -j 16 update compile ' 
## Interface ##                        
alias xa='xmipp_apropos'               
alias xb='xmipp_browser'               
alias xp='xmipp_protocols'             
alias xmipp='xmipp_protocols'          
alias xs='xmipp_show'                  
alias xsj='xmipp_showj'                
alias xmipp_imagej='java -Xmx512m -jar /usr/local/xmipp.3.0/external/imagej/ij.jar'
alias xij='xmipp_imagej'               
## Image ##                            
alias xic='xmipp_image_convert'        
alias xih='xmipp_image_header'         
alias xio='xmipp_image_operate'        
alias xis='xmipp_image_statistics'     
## Metadata ##                         
alias xmu='xmipp_metadata_utilities'   
alias xmp='xmipp_metadata_plot'        
## Transformation ##                   
alias xtg='xmipp_transform_geometry'   
alias xtf='xmipp_transform_filter'     
alias xtn='xmipp_transform_normalize'  
## Other ##                            
alias xrf='xmipp_resolution_fsc'       
alias xrs='xmipp_resolution_ssnr'      
                                                                             
                                                                             
## Configuration ##                                                          
                                                                             
# This file will serve to customize some settings of you Xmipp installation  
# Each setting will be in the form o a shell variable set to some value      
                                                                             
#---------- GUI ----------                                                   
# If you set to 1 the value of this variable, by default the programs        
# will launch the gui when call without argments, default is print the help  
export XMIPP_GUI_DEFAULT=0                                                   
                                                                             
# If you set to 0 the value of this variable the script generated            
# by programs gui will not be erased and you can use the same parameters     
export XMIPP_GUI_CLEAN=1                                                     
                                                                             
#---------- Parallel ----------                                              
# This variable will point to your job submition template file               
export XMIPP_PARALLEL_LAUNCH=config_launch.py                                
                                                                             
# If you have .xmipp.cfg in your home folder it will override                
# this configurations                                                        
                                                                             
test -s ~/.xmipp.cfg && . ~/.xmipp.cfg || true                               
                                                                             
