#   I. Foche & J. Cuenca
#   PDB downloading added by R. Marabini 
#   using code from Michael J. Harms (pdb_download.py)

from pyworkflow.em import *
from os.path import isfile
import os, ftplib, gzip

class XmippProtConvertPdb(ProtInitialVolume):
    """ Covert a PDB file to a volume.  """
    _label = 'convert a PDB'
    _pdb_file = ''
    _sampling_rate = 0.0
    _output_file = ''
    PDB_ID = 1

    HOSTNAME = "ftp.wwpdb.org"
    DIRECTORY = "/pub/pdb/data/structures/all/pdb/"
    PREFIX = "pdb"
    SUFFIX = ".ent.gz"
    # TODO unzip may go to utilities
    def unZip(self, some_file, some_output):
        """
        Unzip some_file using the gzip library and write to some_output.
        CAUTION: deletes some_file.
        """

        f = gzip.open(some_file, 'r')
        g = open(some_output, 'w')
        g.writelines(f.readlines())
        f.close()
        g.close()

        os.remove(some_file)

    def pdbDownload(self, pdbId, fileOut):
        """
        Download all pdb files in file_list and unzip them.
        """
        success = True

        # Log into server
        print "Connecting..."
        ftp = ftplib.FTP()
        ftp.connect(self.HOSTNAME)
        ftp.login()

        # Download  file
        _fileIn = "%s/%s%s%s" % (self.DIRECTORY, self.PREFIX, pdbId, self.SUFFIX) 
        _fileOut = fileOut + ".gz"
        try:
            ftp.retrbinary("RETR %s" % _fileIn, open(_fileOut, "wb").write)
            print "Unzip file %s" % _fileOut
            self.unZip(_fileOut, fileOut) 
            print "%s retrieved successfully." % fileOut
        except ftplib.error_perm:
            os.remove(_fileOut)
            print "ERROR!  %s could not be retrieved!" % _fileIn
            success = False

        # Log out
        ftp.quit()

        if success:
            return True
        else: 
            return False


    def _defineParams(self, form):
        """ Define the parameters that will be input for the Protocol.
        This definition is also used to generate automatically the GUI.
        """
        form.addSection(label='Input')
        form.addParam('inputPdbData', EnumParam, choices=['local_file', 'PDB_ID'],
                      label="Retrieve data from", default=self.PDB_ID,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Retrieve PDB data from server or use local file')
        form.addParam('pdb_file', StringParam, label="FileName or PDB ID",
                       help='type local PDB File Name or PDB ID')
        form.addParam('sampling', FloatParam, label="Sampling rate",
                      help='Sampling rate (Angstroms/pixel) ')

    def _insertAllSteps(self):
        """ In this function the steps that are going to be executed should
        be defined. Two of the most used functions are: _insertFunctionStep or _insertRunJobStep
        """
        self._pdb_file = self.pdb_file.get()
        self._sampling_rate = self.sampling.get()
        self._inputPdbData = self.inputPdbData.get()
        self._insertFunctionStep('convertPdb')
        self._insertFunctionStep('createOutput')

    def convertPdb(self):
        """ Although is not mandatory, usually is used by the protocol to
        register the resulting outputs in the database.
        """
        import xmipp
        from tempfile import NamedTemporaryFile

        if self._inputPdbData == self.PDB_ID:
            _inFile = self._getTmpPath('%s.pdb' % self._pdb_file)
            _outFile = self._getPath('%s' % self._pdb_file)
            if not self.pdbDownload(self._pdb_file, _inFile):
                return -1
        else:
            _inFile = self._pdb_file
            _outFile = self._getPath(_inFile.rsplit(".", 1)[ 0 ])

        program = "xmipp_volume_from_pdb"
        args = '-i %s --sampling %f -o %s ' % (_inFile, self._sampling_rate, _outFile)
        self.runJob(None, program, args)

    def createOutput(self):
        """ Although is not mandatory, usually is used by the protocol to
        register the resulting outputs in the database.
        """
        volume = Volume()
        if self._inputPdbData == self.PDB_ID:
            _outFile = self._getPath('%s' % self._pdb_file)
        else:
	    _outFile = self._getPath(self._pdb_file.rsplit(".", 1)[ 0 ])
        volume.setFileName(_outFile + '.vol')
        self._defineOutputs(volume=volume)
      
    def _summary(self):
        """ Even if the full set of parameters is available, this function provides
        summary information about an specific run.
        """ 
        summary = [ ] 
        # Add some lines of summary information
        if not hasattr(self, 'pdb_file'):
            summary.append("PDB file not ready yet.")
        else:
            _inFile = self.pdb_file.get()
            
            if self.inputPdbData.get() == self.PDB_ID:
                summary.append("Input PDB ID: %s" % _inFile)
            else:
                summary.append("Input PDB File: %s" % _inFile)
            
            # summary.append("Output volume: %s" % _inFile.rsplit( ".", 1 )[ 0 ]+".vol")
        return summary
      
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        errors = []
        if not (self.inputPdbData.get() == self.PDB_ID):
            if not isfile(self.pdb_file.get()):
                errors = ["File %s does not exists" % self.pdb_file.get()]

        # Add some errors if input is not valid
        return errors
    
