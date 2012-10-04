export XMIPP_HOME=/home/joton/xmippGit
export PATH=$XMIPP_HOME/bin:$PATH
export LD_LIBRARY_PATH=$XMIPP_HOME/lib:$LD_LIBRARY_PATH
# Load configuration file 
test -s /home/joton/xmippGit/.xmipp.cfg && . /home/joton/xmippGit/.xmipp.cfg || true
# Load global autocomplete file 
test -s /home/joton/xmippGit/.xmipp.autocomplete && . /home/joton/xmippGit/.xmipp.autocomplete || true

 
 
# Xmipp Aliases 						 
## Setup ##                        
alias xconfigure='./setup.py -j 4 configure compile ' 
alias xcompile='./setup.py -j 4 compile ' 
alias xupdate='./setup.py -j 4 update compile ' 
## Interface ##                        
alias xa='xmipp_apropos'               
alias xb='xmipp_browser'               
alias xp='xmipp_protocols'             
alias xmipp='xmipp_protocols'          
alias xs='xmipp_show'                  
alias xsj='xmipp_showj'                
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
