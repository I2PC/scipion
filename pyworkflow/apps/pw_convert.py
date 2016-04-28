#!/usr/bin/env python
'''
Created on Sep 2, 2015

@author: airen
'''
#!/usr/bin/env python


import argparse
from pyworkflow.em import loadSetFromDb
from pyworkflow.em.data import SetOfCoordinates 
import os       
from pyworkflow.utils import cleanPath
import pyworkflow.utils as pwutils

def main():
    parser = argparse.ArgumentParser(prog='Scipion Convert')
    parser.add_argument('--coordinates', help='Convert coordinates', action="store_true")
    parser.add_argument('--fromType', help='Convert from input type')
    parser.add_argument('--toType', help='Convert to output type')
    parser.add_argument('--input', help='Input file or folder')
    parser.add_argument('--output', help='Output file or folder')
    parser.add_argument('--extra', help='To add extra parameters')

    args = parser.parse_args()
    fromType = args.fromType
    toType = args.toType
    input = args.input
    output = args.output

    if args.coordinates:
        #print 'converting coordinates ...'
        micSet = loadSetFromDb(input)
        outputDir = output
        coordsfn = os.path.join(outputDir, 'coordinates.sqlite')
        cleanPath(coordsfn)
        coordSet = SetOfCoordinates(filename=coordsfn)
        coordSet.setMicrographs(micSet)
        if fromType == 'eman2':
            if toType == 'xmipp': 
                #print 'from eman2 to xmipp...'
                from pyworkflow.em.packages.eman2.convert import readSetOfCoordinates
                readSetOfCoordinates(outputDir, micSet, coordSet)
                from pyworkflow.em.packages.xmipp3.convert import writeSetOfCoordinates
                writeSetOfCoordinates(outputDir, coordSet, ismanual=False)
        elif fromType == 'dogpicker':
            if toType == 'xmipp': 
                #print 'from dogpicker to xmipp...'
                from pyworkflow.em.packages.appion.convert import readSetOfCoordinates
                readSetOfCoordinates(outputDir, micSet, coordSet)
                from pyworkflow.em.packages.xmipp3.convert import writeSetOfCoordinates
                writeSetOfCoordinates(outputDir, coordSet, ismanual=False)
        elif fromType == 'relion':
            if toType == 'xmipp': 
                #print 'from relion to xmipp...'
                inputCoords = args.extra
                starFiles = [os.path.join(inputCoords, pwutils.removeBaseExt(mic.getFileName()) + '_autopick.star')
                     for mic in micSet]
                from pyworkflow.em.packages.relion.convert import readSetOfCoordinates
                readSetOfCoordinates(coordSet, starFiles)
                from pyworkflow.em.packages.xmipp3.convert import writeSetOfCoordinates
                writeSetOfCoordinates(outputDir, coordSet, ismanual=False)
        elif fromType == 'gautomatch':
            if toType == 'xmipp':
                from pyworkflow.em.packages.gautomatch.convert import readSetOfCoordinates
                readSetOfCoordinates(outputDir, micSet, coordSet)
                from pyworkflow.em.packages.xmipp3.convert import writeSetOfCoordinates
                writeSetOfCoordinates(outputDir, coordSet, ismanual=False)
        elif fromType == 'gempicker':
            if toType == 'xmipp':
                from pyworkflow.em.packages.igbmc.convert import readSetOfCoordinates
                readSetOfCoordinates(outputDir, micSet, coordSet)
                from pyworkflow.em.packages.xmipp3.convert import writeSetOfCoordinates
                writeSetOfCoordinates(outputDir, coordSet, ismanual=False)
        
    
if __name__ == '__main__':
    main()
    
    

