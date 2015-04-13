#!/usr/bin/env python

from pyworkflow.em.viewer import ChimeraVirusClient
import os#, sys
import argparse
#import xmipp

def main():
    commonParser = argparse.ArgumentParser(add_help=False, prog='Chimera Virus Client')
    commonParser.add_argument('--input', help='Volume to visualize', required=True)
    commonParser.add_argument('--samplingRate', help='Volume sampling rate (default 1)', default=1)
    commonParser.add_argument('-h', type=int, help='Capsid parameter h (default 5)', default=5)
    commonParser.add_argument('-k', type=int, help='Capsid parameter k (default 0)', default=0)
    commonParser.add_argument('--sym', help='symmetry', default='n25')
    commonParser.add_argument('--Rhex', type=float, help='hexagomal grid radius')
    commonParser.add_argument('--Rsph', type=float, help='sphere shell radius')
    commonParser.add_argument('--rsph', type=float, help='small sphere radius')
    commonParser.add_argument('--sphere', type=float, help='ajust icosahedron to sphere. range [0-1], 0-> icosahedron')
    #parentParser = argparse.ArgumentParser(add_help=False, prog='Chimera Virus Client')
    args = commonParser.parse_args()
    #print args

    volfile   = args.input
    voxelSize = args.samplingRate if hasattr(args, 'samplingRate') else None
    #angularDistFile = args.angDistFile if hasattr(args, 'angDistFile') else None
    h = args.h
    k = args.k
    sym = args.sym
    radius = args.Rhex
    spheRadius = args.rsph
    shellRadius = args.Rsph
    sphere = args.sphere

    ChimeraVirusClient(volfile, voxelSize=voxelSize,
                       h=h, k=k, sym=sym, radius=radius,
                       color = 'red', linewidth=4,
                       spheRadius=spheRadius,sphere=sphere,
                       shellRadius=shellRadius)

    
def which(file):
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
                return path + "/" + file

    return None
    
if __name__ == '__main__':
    main()
    
    

