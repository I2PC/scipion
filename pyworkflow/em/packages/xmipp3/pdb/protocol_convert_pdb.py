#   I. Foche & J. Cuenca
#   PDB downloading added by R. Marabini 
#   using code from Michael J. Harms (pdb_download.py)

import os, ftplib, gzip
import sys

from pyworkflow.em import *
from pyworkflow.utils import *
from os.path import isfile

class XmippProtConvertPdb(ProtInitialVolume):
    """ Covert a PDB file to a volume.  """
    _label = 'convert a PDB'
    _pdb_file = ''
    _sampling_rate = 0.0
    _input_file = ''
    _output_file = ''
    _pdb_id = 1
    _size = 0

    # TODO unzip may go to utilities
    def unzipStep(self, some_file, some_output):
        """
        Unzip some_file using the gzip library and write to some_output.
        CAUTION: deletes some_file.
        """
        success = True
        try:
            f = gzip.open(some_file, 'r')
            g = open(some_output, 'w')
            g.writelines(f.readlines())
            f.close()
            g.close()
        except:
            e = sys.exc_info()[0]
            self.error('ERROR opening gzipped file %s: %s' % (some_file, e))
            success = False

        try:
            if success:
                os.remove(some_file)
        except:
            e = sys.exc_info()[0]
            self.error('ERROR deleting gzipped file: %s' % e)
            success = False

    def pdbDownloadStep(self, pdbId, fileOut):
        """
        Download all pdb files in file_list and unzip them.
        """
        pdborgHostname = "ftp.wwpdb.org"
        pdborgDirectory = "/pub/pdb/data/structures/all/pdb/"
        prefix = "pdb"
        suffix = ".ent.gz"
        success = True

        # Log into server
        print "Connecting..."
        ftp = ftplib.FTP()
        try:
            ftp.connect(pdborgHostname)
            ftp.login()
        except ftplib.error_temp:
            self.error("ERROR! Timeout reached!")
            success = False

        if success:
            # Download  file
            _fileIn = "%s/%s%s%s" % (pdborgDirectory, prefix, pdbId, suffix) 
            _fileOut = fileOut + ".gz"
            try:
                ftp.retrbinary("RETR %s" % _fileIn, open(_fileOut, "wb").write)
            except ftplib.error_perm:
                os.remove(_fileOut)
                self.error("ERROR!  %s could not be retrieved!" % _fileIn)
                success = False
            # Log out
            ftp.quit()

    def _defineParams(self, form):
        """ Define the parameters that will be input for the Protocol.
        This definition is also used to generate automatically the GUI.
        """

        form.addSection(label='Input')
        form.addParam('pdb_file', PathParam, 
                      label="Pattern",
                      help='Specify a path or an url to desired PDB structure. Also wwpdb.org IDs are accepted')
#        form.addParam('pdb_file', StringParam, label="FileName or PDB ID",
#                       help='type local PDB File Name or PDB ID')
        form.addParam('inputPdbData', EnumParam, choices=['local_file', 'pdb_id'],
                      label="Retrieve data from", default=self._pdb_id,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Retrieve PDB data from server or use local file')
        form.addParam('sampling', FloatParam, default=1.0, label="Sampling rate (A/px)",
                      help='Sampling rate (Angstroms/pixel)')
        form.addParam('setSize', BooleanParam, label='Set final size?', default=False)
        form.addParam('size', IntParam, label="Final size (px)", condition='setSize', allowsNull=True,
                      help='Final size in pixels. If no value is provided, protocol will estimate it.')
        form.addParam('centerPdb', BooleanParam, label="Center PDB", default=True, expertLevel=LEVEL_EXPERT,
                      help='Center PDB with the center of mass')

    def _insertAllSteps(self):
        """ In this function the steps that are going to be executed should
        be defined. Two of the most used functions are: _insertFunctionStep or _insertRunJobStep
        """
        # Taking all arguments
        self._pdb_file = self.pdb_file.get()
        self._sampling_rate = self.sampling.get()
        self._setSize = self.setSize.get()
        self._size = self.size.get()
        self._centerPdb = self.centerPdb.get()
        self._inputPdbData = self.inputPdbData.get()
        
        if self._inputPdbData == self._pdb_id:
            _inFile = self._getTmpPath('%s.pdb' % self._pdb_file)
            self._input_file = _inFile
            _outFile = self._getPath('%s' % self._pdb_file)
            self._pdb_file = self._pdb_file.lower()
            self.info("File to download and unzip: %s" % (_inFile+".gz"))
            self._insertFunctionStep('pdbDownloadStep', self._pdb_file, _inFile)
            self._insertFunctionStep('unzipStep',(_inFile+".gz"), (_inFile))
        self._insertFunctionStep('convertPdbStep')
        self._insertFunctionStep('createOutput')

    def convertPdbStep(self):
        """ Although is not mandatory, usually is used by the protocol to
        register the resulting outputs in the database.
        """
        import xmipp
        from tempfile import NamedTemporaryFile
        
        centerArg = _outFile = ''
        if self._centerPdb:
            centerArg = '--centerPDB'

        sizeArg = ''
        if self._setSize:
            sizeArg = '--size'
        if self.size.hasValue():
            sizeArg += ' %d' % (self._size)

        if self._inputPdbData == self._pdb_id:
            _inFile = self._input_file
        else:
            _inFile = self._pdb_file
        _outFile = self._getPath(removeBaseExt(_inFile))
        self.info("Input file: "+_inFile)
        self.info("Output file: "+_outFile)
        self._output_file = _outFile + ".vol"
        program = "xmipp_volume_from_pdb"
        args = '-i %s --sampling %f -o %s %s %s' % (_inFile, self._sampling_rate, _outFile, centerArg, sizeArg)
        self.runJob(program, args)

    def createOutput(self):
        """ Although is not mandatory, usually is used by the protocol to
        register the resulting outputs in the database.
        """
        volume = Volume()
        _outFile = self._output_file
        #    _outFile = self._getPath(self._pdb_file.rsplit(".", 1)[ 0 ])
        #_outFile = removeExt(self._output_file)
        print "createOutput output file:"+self._output_file+" "+_outFile
        volume.setFileName(self._output_file)
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
            
            if self.inputPdbData.get() == self._pdb_id:
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
        if not (self.inputPdbData.get() == self._pdb_id):
            if not isfile(self.pdb_file.get()):
                errors = ["File %s does not exists" % self.pdb_file.get()]

        # Add some errors if input is not valid
        return errors
    
