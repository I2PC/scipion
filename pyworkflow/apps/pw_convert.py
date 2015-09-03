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

def main():
    parser = argparse.ArgumentParser(prog='Scipion Convert')
    parser.add_argument('--coordinates', help='Convert coordinates', action="store_true")
    parser.add_argument('--fromType', help='Convert from input type')
    parser.add_argument('--toType', help='Convert to output type')
    parser.add_argument('--input', help='Input file or folder')
    parser.add_argument('--output', help='Output file or folder')



    
    args = parser.parse_args()
    fromType = args.fromType
    toType = args.toType
    input = args.input
    output = args.output
    if args.coordinates:
        print 'converting coordinates ...'
        if fromType == 'dogpicker':
            if toType == 'xmipp': 
                print 'from dogpicker to xmipp...'
                micSet = loadSetFromDb(input)
                workDir = output
                coordsfn = os.path.join(workDir, 'coordinates.sqlite')
                cleanPath(coordsfn)
                coordSet = SetOfCoordinates(filename=coordsfn)
                coordSet.setMicrographs(micSet)
                from pyworkflow.em.packages.appion.convert import readSetOfCoordinates
                readSetOfCoordinates(workDir, micSet, coordSet)
                from pyworkflow.em.packages.xmipp3.convert import writeSetOfCoordinates
                writeSetOfCoordinates(workDir, coordSet, ismanual=False)
    
if __name__ == '__main__':
    main()
    
    

