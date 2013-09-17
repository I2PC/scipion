#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sanchez Sorzano, September 2013
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Image Operate
#------------------------------------------------------------------------------------------------
# {file}(images*.xmd){validate}(PathExists) 1st operand
""" This selfile points to the stack or metadata containing your images or volumes
"""
Operand1=''

# {list_combo}(plus, minus, multiply, divide, minimum, maximum, dot product, log, log10, sqrt, abs, pow, slice, column, row, radial average, reset) Operation
Operation='plus'

# {file}{condition}(Operation=='plus' or Operation=='minus' or Operation=='multiply' or Operation=='divide' or Operation=='minimum' or Operation=='maximum' or Operation=='dot product' or Operation=='column' or Operation=='slice' or Operation=='row') 2nd operand 
""" It can be a metadata, a stack of images/volumes or a value """
Operand2=''

# {eval} expandParallel(threads=0)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_image_operate import *

if __name__ == '__main__':
    protocolMain(ProtImageOperate)
