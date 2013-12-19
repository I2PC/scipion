from pyworkflow.em import *

class XmippProtConvertPdb(ProtInitialVolume):
    """ Covert a PDB file to a volume.  """
    _label = 'convert a PDB'
    _pdb_file = ''
    _sampling_rate = 0.0
    _output_file = ''
    PDB_ID=1

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

    def _defineSteps(self):
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
		import urllib
		src = "http://www.rcsb.org/pdb/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s"
		urllib.urlretrieve(src % self._pdb_file,'%s.pdb'%self._pdb_file)
		_inFile = '%s.pdb' % self._pdb_file
		_outFile= '%s' % self._pdb_file
	else:
		_inFile = self._pdb_file
		_outFile = _inFile.rsplit( ".", 1 )[ 0 ]+".vol"

        program = "xmipp_volume_from_pdb"
        args = '-i %s --sampling %f -o %s ' % (_inFile, self._sampling_rate,_outFile)
        self.runJob(None, program, args)

    def createOutput(self):
        """ Although is not mandatory, usually is used by the protocol to
        register the resulting outputs in the database.
        """
        from pyworkflow.utils import *
        volume = Volume()
        self._output_file = replaceExt(self._pdb_file,'vol') 
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
            summary.append("Input PDB: %s" % self.pdb_file.get())
            summary.append("Output volume: %s" % self._output_file)
        return summary
      
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        errors = [ ] 
        # Add some errors if input is not valid
        return errors
    
