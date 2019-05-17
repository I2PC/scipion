# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
import sys

import os
import sys

def getSubsetByDefocus(inputCTFs, inputMics, nMics):
    """ Return a subset of inputMics that covers the whole range of defocus
    from the inputCtfs set.
    This function can be used from picking wizards that wants to optimize the
    parameters for micrographs with different defocus values.
    Params:
        nMics is the number of micrographs that will be in the subset.
    """
    sortedMicIds = []

    # Sort CTFs by defocus and select only those that match with inputMics
    for ctf in inputCTFs.iterItems(orderBy='_defocusU'):
        ctfId = ctf.getObjId()
        if ctfId in inputMics:
            sortedMicIds.append(ctfId)

    # Take an equally spaced subset of micrographs
    space = len(sortedMicIds) / (nMics - 1)
    micIds = [sortedMicIds[0], sortedMicIds[-1]]
    pos = 0
    while len(micIds) < nMics:  # just add first and last
        pos += space
        micIds.insert(1, sortedMicIds[pos])

    # Return the list with selected micrographs
    return [inputMics[micId].clone() for micId in micIds]


# TODO: use biopython
def downloadPdb(pdbId, pdbFile, log=None):
    print """use AtomicStructHandler()

    aSH = AtomicStructHandler()
    pdbFileName = aSH.readFromPDBDatabase(pdbId, type='mmCif',
                                          dir=os.getcwd())
"""
    pdbGz = pdbFile + ".gz"
    result = (__downloadPdb(pdbId, pdbGz, log) and
              __unzipPdb(pdbGz, pdbFile, log))
    return result


# TODO: use biopython
def __downloadPdb(pdbId, pdbGz, log):
    import ftplib
    """Download a pdb file given its id. """
    if log:
        log.info("File to download and unzip: %s" % pdbGz)

    pdborgHostname = "ftp.wwpdb.org"
    pdborgDirectory = "/pub/pdb/data/structures/all/mmCIF/"
    prefix = ""  # use pdb for PDB and null for mmcif
    suffix = ".cif.gz"
    success = True
    # Log into serverhttp://www.rcsb.org/pdb/files/2MP1.pdb.gz
    ftp = ftplib.FTP()
    try:
        ftp.connect(pdborgHostname)
        ftp.login()
    except ftplib.error_temp:
        if log:
            log.error("ERROR! Timeout reached!")
        success = False

    if success:
        # Download  file
        _fileIn = "%s/%s%s%s" % (pdborgDirectory, prefix, pdbId.lower(), suffix)
        _fileOut = pdbGz
        try:
            ftp.retrbinary("RETR %s" % _fileIn, open(_fileOut, "wb").write)
        except ftplib.error_perm:
            os.remove(_fileOut)
            if log:
                log.error("ERROR!  %s could not be retrieved!" % _fileIn)
            success = False
        # Log out
        ftp.quit()

    return success


def __unzipPdb(pdbGz, pdbFile, log, cleanFile=True):
    """
    Unzip a pdb file.
    Params:
        pdbGz: zipped pdb file.
        pdbFile: output pdb file.
        cleanFile: remove the zipped file.
    """
    import gzip
    success = True
    try:
        f = gzip.open(pdbGz, 'r')
        g = open(pdbFile, 'w')
        g.writelines(f.readlines())
        f.close()
        g.close()
    except:
        e = sys.exc_info()[0]
        if log:
            log.error('ERROR opening gzipped file %s: %s' % (pdbGz, e))
        success = False

    try:
        if success:
            os.remove(pdbGz)
    except:
        e = sys.exc_info()[0]
        if log:
            log.error('ERROR deleting gzipped file: %s' % e)
        success = False

    return success