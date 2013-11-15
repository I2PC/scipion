#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Input  
#------------------------------------------------------------------------------------------------
# {file}(*.emx) EMX file name
"""EMX file name:"""
EmxFileName = ''


#------------------------------------------------------------------------------------------------
# {section} Acquisition information (fill it if not in EMX file)
#------------------------------------------------------------------------------------------------

# {validate}(IsFloatOrEmpty) Microscope voltage (in kV)
""" If left blank, the value will be read from the EMX file. """
Voltage = ""

# {validate}(IsFloatOrEmpty) Spherical aberration (in mm)
""" If left blank, the value will be read from the EMX file. """
SphericalAberration = ""

# {validate}(IsFloatOrEmpty) Sampling rate (A/pixel)
""" If left blank, the value will be read from the EMX file. """
SamplingRate = ""

# {validate}(IsFloatOrEmpty) Percentage of amplitude contrast
""" If left blank, the value will be read from the EMX file. """
AmplitudeContrast = ""

    #------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_emx_import_micrographs import *

if __name__ == '__main__':
    protocolMain(ProtEmxImportMicrographs)
