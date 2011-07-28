from xmipp import *
import os
from protlib_utils import runJob

def scaleImages(_log
                    , dimX
                    , dimY
                    , filename_currentAngles
                    , MpiJobSize
                    , NumberOfMpiProcesses
                    , NumberOfThreads
                    , scaledImages
                    , SystemFlavour
                    ):

    outFileName=scaledImages + ".stk"
    if os.path.exists(outFileName):
        os.remove(outFileName)

    if (dimY<0):
        dimY = dimX
    
    parameters  = ' -i ' +  filename_currentAngles 
    parameters += ' -o ' + outFileName 
    parameters += ' --scale fourier ' + str(dimX) + ' ' + str(dimY) + ' ' + str(NumberOfThreads)
    parameters += ' --disable_metadata' 

    runJob(_log,'xmipp_transform_geometry',
                     parameters,
                     NumberOfMpiProcesses * NumberOfThreads,
                     1, # Threads go in --scale option
                     SystemFlavour)

    #merge new scaled metadata with original metadata
    fnScaledImages = scaledImages+".xmd"
    mdScaled = MetaData(fnScaledImages)
    md = MetaData(filename_currentAngles)
    
    # copy images to original images column
    md.addLabel(MDL_IMAGE_ORIGINAL)
    for id in md:
        imageScale=mdScaled.getValue(MDL_IMAGE, id)
        imageOriginal=md.getValue(MDL_IMAGE, id)
        md.setValue(MDL_IMAGE, imageScale, id)
        md.setValue(MDL_IMAGE_ORIGINAL, imageOriginal, id)
    
    # Calculate scale factor
    (x,y,z,n) =ImgSize(filename_currentAngles)
    factorX = float(x) / dimX
    factorY = float(y) / dimY

    scaleXstr = 'shiftX=(shiftX /  ' + str(factorX) + ')'
    scaleYstr = 'shiftY=(shiftY /  ' + str(factorY) + ')'
    md.operate(scaleXstr)
    md.operate(scaleYstr)
    
    md.write(fnScaledImages)
    
