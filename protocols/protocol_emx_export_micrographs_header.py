#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Input  
#------------------------------------------------------------------------------------------------
# {file}(micrographs*.xmd) Micrographs to export:
"""Select a metadata containing your micrographs to export
as EMX data. If you select the micrographs.xmd corresponding
to a picking run, the .pos coordinates files will be also
exported.
"""
MicrographsMd = ''

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_emx_export_micrographs import *

if __name__ == '__main__':
    protocolMain(ProtEmxExportMicrographs)
