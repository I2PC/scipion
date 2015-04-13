#!/usr/bin/env python

import sys

def diag(n_mode):
	f=open('diag_arpack.in', 'w')
	f.write('&inputs\n')
	f.write('nev =%3d, \n' %(n_mode))
	ncv=n_mode*2
		
	f.write('ncv =%4d,\n' %ncv)
	f.write('tol = 1e-12/')
	f.write('\n')
	f.close()	


def pdbmat(fichier,cutoff):
	f = open('pdbmat.dat', 'w')
	f.write('%10s %3d' %(fichier,cutoff))
	f.close()


# Main
if __name__ == '__main__':
	diag(int(sys.argv[1]))
	pdbmat( sys.argv[2],int(sys.argv[3]))
	sys.exit()




