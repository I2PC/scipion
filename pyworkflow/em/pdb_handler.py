import os
import fileinput
from shutil import copyfile

def fixCRYSrecordToPDBFile( inPDBFileName, tmpDir,
                           x=1, y=1, z=1, alpha=90., beta=90., gamma=90.):
    """ Check if crys record exists, if not add it before first atom line
        input: pdbfile
        output: fixed pdfile filename
    """
    def getfixedPDBFileName():
        return os.path.join(tmpDir, os.path.basename(inPDBFileName))

    # http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
    # https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
    CRYS = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n' % (x, y, z, alpha, beta, gamma)

    addCRYS = True
    with open(inPDBFileName, 'r') as inF:
        for line in inF:
            if line.startswith('CRYS'):
                addCRYS = False
                break
            elif line.startswith('ATOM'):
                break

    if not addCRYS:
        tmpFileName = inPDBFileName
    else:
        ADDCRYS = True
        tmpFileName = getfixedPDBFileName()
        with open(inPDBFileName,"r") as infile, open(tmpFileName, "w") as outfile:
            for line in infile:
                if ADDCRYS and (line.startswith('ATOM') or line.startswith('MODEL')):
                    line = CRYS + line
                    ADDCRYS = False
                outfile.write(line)

    return tmpFileName