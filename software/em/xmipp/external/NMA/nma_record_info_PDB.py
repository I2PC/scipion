#!/usr/bin/env python

import sys

def diag(fnPDB,n_mode,rtbBlockSize):
	f=open('diagrtb.dat', 'w')
	content="""! This file can be modified and used as a command file (named diagrtb.dat) for diagrtb.
 MATRix filename            =   pdbmat.sdijf
 COORdinates filename       =     %s
 Eigenvector OUTPut filename= diagrtb.eigenfacs
 Nb of VECTors required     =         %d
 EigeNVALues chosen         =       LOWE   ! LOWEst, HIGHest.
 Type of SUBStructuring     =       NONE   ! RESIdues, SECOndary, SUBUnits, DOMAins, NONE.
 Nb of residues per BLOck   =         %d 
 Origin of MASS values      =       CONS   ! CONStant, COORdinate, PDB.
 Temporary files cleaning   =       ALL    ! ALL, NOne.
 MATRix FORMat              =       FREE   ! FREE, BINAry.
 Output PRINting level      =          2   ! =1: More detailed; =2: Debug level.
"""%(fnPDB,n_mode,rtbBlockSize)
	f.write(content)
	f.close()	

def pdbmat(fnPDB,cutoff,forceConstant):
	f = open('pdbmat.dat', 'w')
	content="""! This file can be modified and used as a command file (named pdbmat.dat) for pdbmat.
 Coordinate FILENAME        = %s
 MATRIx FILENAME            = pdbmat.sdijf
 INTERACtion DISTance CUTOF = %f ! For defining the list of interacting atoms.
 INTERACtion FORCE CONStant = %f ! For specifying frequency units.
 Origin of MASS values      =    CON ! CONstant, or from COOrdinate file.
 Output PRINTing level      =      2 ! =1: more detailled. =2: debug level.
"""%(fnPDB,cutoff,forceConstant)
	f.write(content)
	f.close()

# Main
if __name__ == '__main__':
	diag(sys.argv[3],int(sys.argv[1]),int(sys.argv[2]))
	pdbmat( sys.argv[3],float(sys.argv[4]),float(sys.argv[5]))
	sys.exit()




