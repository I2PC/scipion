# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from os.path import exists
import os
import pyworkflow.protocol.params as params
from base import ProtImportFiles
from pyworkflow.em.data import Sequence
from pyworkflow.em.convert.sequence import (
    SEQ_TYPE_AMINOACIDS, SEQ_TYPE_NUCLEOTIDES, IUPAC_PROTEIN_ALPHABET,
    SEQ_TYPE, EXTENDED_PROTEIN_ALPHABET, IUPAC_NUCLEOTIDE_ALPHABET,
    EXTENDED_DNA_ALPHABET, SequenceHandler, cleanSequenceScipion,
    alphabetToIndex)
from pyworkflow.em.convert.atom_struct import AtomicStructHandler


class ProtImportSequence(ProtImportFiles):
    """ Protocol to import an aminoacid/nucleotide sequence file to the
    project"""
    _label = 'import sequence'
    #SEQUENCEFILENAME = '_sequence.fasta'
    IMPORT_FROM_PLAIN_TEXT = 0
    IMPORT_FROM_STRUCTURE = 1
    IMPORT_FROM_FILES = 2
    IMPORT_FROM_UNIPROT = 3
    IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT = 0
    IMPORT_FROM_NUCLEOTIDE_STRUCTURE = 1
    IMPORT_FROM_NUCLEOTIDE_FILES = 2
    IMPORT_FROM_GENEBANK = 3
    IMPORT_STRUCTURE_FROM_ID = 0
    IMPORT_STRUCTURE_FROM_FILES = 1

    url = "http://www.uniprot.org/uniprot/"

    def __init__(self, **args):
        ProtImportFiles.__init__(self, **args)

    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputSequenceID', params.StringParam,
                      label="Sequence ID", allowsNull=True,
                      help="Write a sequence ID. Otherwise, if the "
                           "sequence derives from GeneBank/UniProt/PDB "
                           "databases, the respective database ID will be "
                           "selected as starting sequence ID; examples: if "
                           "you select GeneBank accession AJ520101, SCIPION "
                           "will assign AJ520101 as sequence ID; if "
                           "you select UniProt accession P12345, SCIPION will "
                           "assign P12345 as sequence ID; if you "
                           "select atomic structure 3lqd.cif, chain B, "
                           "SCIPION will assign 3lqd_B as sequence ID. In "
                           "the rest of cases, the Sequence name "
                           "will be selected as starting Sequence ID.")
        form.addParam('inputSequenceName', params.StringParam, important=True,
                      label="Sequence name", allowsNull=False,
                      help="Write a sequence name.")
        form.addParam('inputSequenceDescription', params.StringParam,
                      label="Sequence description",
                      allowsNull=True,
                      help="Write a description for your sequence. Otherwise, "
                           "if the "
                           "sequence derives from GeneBank/UniProt/PDB "
                           "databases, the respective database description "
                           "will be "
                           "selected as starting sequence description. In "
                           "the rest of cases, no sequence description will "
                           "be added.")
        form.addParam('inputSequence', params.EnumParam,
                      pointerClass='Sequence',
                      choices=SEQ_TYPE,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label="Import sequence of ",
                      default=SEQ_TYPE_AMINOACIDS,
                      help='Select the type of sequence to import.')
        form.addParam('inputProteinSequence', params.EnumParam,
                      choices=['plain text', 'atomic structure', 'file',
                               'UniProt ID'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='inputSequence == %d' % SEQ_TYPE_AMINOACIDS,
                      label="From ",
                      default=self.IMPORT_FROM_PLAIN_TEXT,
                      help='Select one of the four options: write the '
                           'aminoacid sequence or import it '
                           'from a previously loaded atomic structure, a local '
                           'file or an online server.')
        form.addParam('proteinIUPACalphabet', params.EnumParam,
                      choices=IUPAC_PROTEIN_ALPHABET,
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d' %
                                (SEQ_TYPE_AMINOACIDS,
                                 self.IMPORT_FROM_PLAIN_TEXT),
                      label="IUPAC Protein alphabet: ",
                      default=EXTENDED_PROTEIN_ALPHABET,
                      help='Your raw sequence will be cleaned according '
                           'a certain alphabet, i.e., only the letters '
                           'contained in the alphabet will be maintained in '
                           'the sequence. Select thus the type of protein '
                           'alphabet in order to accomplish the '
                           'cleaning:\n\nProtein alphabet: IUPAC protein '
                           'alphabet of the 20 standard amino acids; uppercase'
                           ' and single letter: *ACDEFGHIKLMNPQRSTVWY*.\n\n'
                           'Extended Protein alphabet: Extended uppercase '
                           'IUPAC '
                           'protein single letter alphabet including X etc.\n'
                           'In addition to the standard 20 single letter '
                           'protein codes, this includes:\n*B = Asx*; '
                           'Aspartic acid (R) or Asparagine (N)\n*X = Xxx*"; '
                           'Unknown or other amino acid\n*Z = Glx*; Glutamic '
                           'acid (E) or Glutamine (Q)\n*J = Xle*; Leucine ('
                           'L) or Isoleucine (I), used in mass-spec (NMR)\n'
                           '*U = Sec*; Selenocysteine\n*O = Pyl*; '
                           'Pyrrolysine\nThis alphabet is not intended to be '
                           'used with X for Selenocysteine (an ad-hoc standard'
                           ' prior to the IUPAC adoption of U instead).\n')
        form.addParam('uniProtSequence', params.StringParam,
                      condition='inputSequence == %d and '
                                'inputProteinSequence == %d' %
                                (SEQ_TYPE_AMINOACIDS,
                                 self.IMPORT_FROM_UNIPROT),
                      label="UniProt name/ID ", allowsNull=True,
                      help='Write a UniProt ID (six or ten alphanumeric '
                           'characters; examples: A2BC19, P12345, '
                           'A0A022YWF9, DGAL_ECOLI).\n You can convert other '
                           'database identifiers to UniProt accession codes '
                           'by using the "ID Mapping" tab on '
                           'https://www.uniprot.org/')
        form.addParam('inputNucleotideSequence', params.EnumParam,
                      choices=['plain text', 'atomic structure', 'file',
                               'GeneBank'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='inputSequence == %d' % SEQ_TYPE_NUCLEOTIDES,
                      label="From ",
                      default=self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT,
                      help='Select one of the four options: write the '
                           'nucleic acid sequence or import it '
                           'from a local file or an online server.')
        form.addParam('nucleotideIUPACalphabet', params.EnumParam,
                      choices=IUPAC_NUCLEOTIDE_ALPHABET,
                      display=params.EnumParam.DISPLAY_HLIST,
                      condition='inputSequence == %d and '
                                'inputNucleotideSequence == %d' %
                                (SEQ_TYPE_NUCLEOTIDES,
                                 self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT),
                      label="IUPAC Nucleic acid alphabet: ",
                      default=EXTENDED_DNA_ALPHABET,
                      help='Your raw sequence will be cleaned according '
                           'a certain alphabet, i.e., only the letters '
                           'contained in the alphabet will be maintained in '
                           'the sequence. Select thus the type of nucleic '
                           'acid alphabet in order to accomplish the '
                           'cleaning:\n\n Ambiguous DNA alphabet: Uppercase '
                           'IUPAC ambiguous DNA: *GATCRYWSMKHBVDN*.\n\n'
                           'Unambiguous DNA alphabet: Uppercase IUPAC unambiguous DNA '
                           '(letters *GATC* only).\n\nExtended DNA: Extended '
                           'IUPAC DNA alphabet.\nIn addition to the standard letter '
                           'codes GATC, this includes:\n*B* = 5-bromouridine\n'
                           '*D* = 5,6-dihydrouridine\n*S* = thiouridine\n*W* '
                           '= wyosine\n\nAmbiguous RNA: Uppercase IUPAC '
                           'ambiguous RNA; *GAUCRYWSMKHBVDN*\n\nUnambigous '
                           'RNA alphabet: Generic single letter RNA '
                           'alphabet.\n\n')
        form.addParam('inputRawSequence', params.StringParam,
                      condition='(inputSequence == %d and '
                                'inputProteinSequence == %d) or '
                                '(inputSequence == %d and '
                                'inputNucleotideSequence == %d) ' %
                                (SEQ_TYPE_AMINOACIDS,
                                self.IMPORT_FROM_PLAIN_TEXT,
                                 SEQ_TYPE_NUCLEOTIDES,
                                 self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT),
                      label="Write your sequence here:", important=True,
                      help="Write the aminoacid or nucleotide raw sequence.\n")
        form.addParam('inputStructureSequence', params.EnumParam,
                      choices=['id', 'file'],
                      condition='inputProteinSequence == %d or '
                                'inputNucleotideSequence == %d' %
                                (self.IMPORT_FROM_STRUCTURE,
                                 self.IMPORT_FROM_NUCLEOTIDE_STRUCTURE),
                      label="Atomic structure from",
                      default=self.IMPORT_STRUCTURE_FROM_ID,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Import structure data from online server or local '
                           'file',
                      pointerClass='AtomStruct',
                      allowsNull=True)
        form.addParam('pdbId', params.StringParam,
                      condition='(inputProteinSequence == %d or '
                                'inputNucleotideSequence == %d) and '
                                'inputStructureSequence == %d'
                                % (self.IMPORT_FROM_STRUCTURE,
                                   self.IMPORT_FROM_NUCLEOTIDE_STRUCTURE,
                                   self.IMPORT_STRUCTURE_FROM_ID),
                      label="Atomic structure ID ", allowsNull=True,
                      help='Type a structure ID (four alphanumeric '
                           'characters).')
        form.addParam('pdbFile', params.PathParam, label="File path",
                      condition='(inputProteinSequence == %d or '
                                'inputNucleotideSequence == %d) and '
                                'inputStructureSequence == %d'
                                % (self.IMPORT_FROM_STRUCTURE,
                                   self.IMPORT_FROM_NUCLEOTIDE_STRUCTURE,
                                   self.IMPORT_STRUCTURE_FROM_FILES),
                      allowsNull=True,
                      help='Specify a path to desired atomic structure.')
        form.addParam('inputStructureChain', params.StringParam,
                      condition='inputProteinSequence == %d or '
                                'inputNucleotideSequence == %d' %
                                (self.IMPORT_FROM_STRUCTURE,
                                 self.IMPORT_FROM_NUCLEOTIDE_STRUCTURE),
                      label="Chain ", allowsNull=True,
                      help="Select a particular chain of the atomic "
                           "structure.")
        form.addParam('fileSequence', params.PathParam,
                      label="File path",
                      condition='inputProteinSequence == %d or '
                                'inputNucleotideSequence == %d' %
                                (self.IMPORT_FROM_FILES,
                                 self.IMPORT_FROM_NUCLEOTIDE_FILES),
                      allowsNull=True,
                      help='Specify a path to desired aminoacid or '
                           'nucleic acid sequence '
                           'file.\nIf your file contains more than one '
                           'sequence, only the first one will be considered.')
        form.addParam('geneBankSequence', params.StringParam,
                      condition='inputSequence == %d and '
                                'inputNucleotideSequence == %d' %
                                (SEQ_TYPE_NUCLEOTIDES,
                                self.IMPORT_FROM_GENEBANK),
                      label="GeneBank accession ", allowsNull=True,
                      help='Write a GeneBank accession.\n')

    def _insertAllSteps(self):
        self.name = self.inputSequenceName.get()

        if self.inputSequence == SEQ_TYPE_AMINOACIDS:
            if self.inputProteinSequence == self.IMPORT_FROM_PLAIN_TEXT:
                rawSequence = self.inputRawSequence.get()
                self._insertFunctionStep('getRawSequenceStep', rawSequence)
            elif self.inputProteinSequence == self.IMPORT_FROM_STRUCTURE:
                chainId = self.inputStructureChain.get()
                self._insertFunctionStep('getSequenceOfChainStep', chainId)
            elif self.inputProteinSequence == self.IMPORT_FROM_UNIPROT:
                sequenceDB = self._getUniProtID()
                self._insertFunctionStep('sequenceDatabaseDownloadStep',
                                         sequenceDB)
            elif self.inputProteinSequence == self.IMPORT_FROM_FILES:
                self.sequenceFile = self.fileSequence.get()
                sequenceFile = self.sequenceFile
                self._insertFunctionStep('fileDownloadStep', sequenceFile)
        else:
            if self.inputNucleotideSequence == \
                    self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT:
                rawSequence = self.inputRawSequence.get()
                self._insertFunctionStep('getRawSequenceStep', rawSequence)
            elif self.inputNucleotideSequence == \
                     self.IMPORT_FROM_NUCLEOTIDE_STRUCTURE:
                chainId = self.inputStructureChain.get()
                self._insertFunctionStep('getSequenceOfChainStep', chainId)
            elif self.inputNucleotideSequence == self.IMPORT_FROM_GENEBANK:
                sequenceDB = self._getGeneBankID()
                self._insertFunctionStep('sequenceDatabaseDownloadStep',
                                         sequenceDB)
            elif self.inputNucleotideSequence == \
                self.IMPORT_FROM_NUCLEOTIDE_FILES:
                self.sequenceFile = self.fileSequence.get()
                sequenceFile = self.sequenceFile
                self._insertFunctionStep('fileDownloadStep', sequenceFile)

        self._insertFunctionStep('createOutputStep')

    def getRawSequenceStep(self, rawSequence):
        # user types sequence
        if self.inputSequenceID.get() is not None:
            self.id = self.inputSequenceID.get()
        else:
            self.id = self.name
        self.alphabet = self._getAlphabet()  # index number
        self.sequence = cleanSequenceScipion(self.inputSequence ==
                                     SEQ_TYPE_AMINOACIDS,
                                      self.alphabet, rawSequence)

    def getSequenceOfChainStep(self, chainId):
        # sequece is obtained from PDB file

        # form has a wizard that creates label with the format
        # [model: x, chain: x, xxx residues]
        import json
        chainIdDict = json.loads(self.inputStructureChain.get())

        selectedModel = chainIdDict['model']
        selectedChain = chainIdDict['chain']
        self.structureHandler = AtomicStructHandler()

        if self.pdbId.get() is not None:
        # PDB from remote database
            pdbID = self.pdbId.get()
            tmpFilePath = os.path.join("/tmp", pdbID + ".cif")
            if exists(tmpFilePath):
                # wizard already downloaded the file
                self.structureHandler.read(tmpFilePath)
            else:
                # wizard has not used and the file has not been downloaded yet
                self.structureHandler.readFromPDBDatabase(pdbID, dir="/tmp")
        else:
        # PDB from file
            self.structureHandler.read(self.pdbFile.get())

        _sequence = self.structureHandler.getSequenceFromChain(
            selectedModel, selectedChain)
        self.sequence = str(_sequence)
        self.alphabet = alphabetToIndex(self.inputSequence ==
                                        SEQ_TYPE_AMINOACIDS,
                                        _sequence.alphabet)

        # Assignation of sequence ID: if the user has provided a specific
        #  ID, this will be adopted by default; otherwise, a sequence ID
        # related with the starting structure will be selected.
        if self.inputSequenceID.get() is not None:
            self.id = self.inputSequenceID.get()
        else:
            self.id = self.structureHandler.getFullID(
                    selectedModel, selectedChain)

        print "Selected chain: %s from model: %s from structure: %s" \
              % (selectedChain, selectedModel,
                 self.structureHandler.structure.get_id())

    def sequenceDatabaseDownloadStep(self, sequenceDB):
        """Download UniProt/GeneBank sequence from its respective database
        """
        #sequenceDB = str(sequenceDB)
        if self.uniProtSequence.get() is not None:
            seqHandler = SequenceHandler()

        elif self._getGeneBankID() is not None:
            seqHandler = SequenceHandler(isAminoacid=False)

        record, error = seqHandler.downloadSeqFromDatabase(sequenceDB)
        if record is None:
            print("Error: ", error)
            self.setAborted()
            exit(0)

        if self.inputSequenceID.get() is not None:
            self.id = self.inputSequenceID.get()
        elif sequenceDB != '':
            self.id = sequenceDB
        else:
            self.id = self.name
        if record.description != '':
            self.description = record.description

        self.sequence = str(record.seq)
        self.alphabet = alphabetToIndex(self.inputSequence ==
                                     SEQ_TYPE_AMINOACIDS, record.seq.alphabet)

    def fileDownloadStep(self, sequenceFile):
        # If sequencePath contains more than one sequence, only
        # the first one will be considered
        seqHandler = SequenceHandler()
        record = seqHandler.downloadSeqFromFile(sequenceFile,
                                                        type="fasta")
        if self.inputSequenceID.get() is not None:
            self.id = self.inputSequenceID.get()
        elif record.id != '':
            self.id = record.id
        else:
            self.id = self.name
        if record.description != '':
            self.description = record.description

        self.sequence = str(record.seq)
        self.alphabet = alphabetToIndex(self.inputSequence ==
                                     SEQ_TYPE_AMINOACIDS, record.seq.alphabet)

    def createOutputStep(self):
        """ Register the output object. """

        if self.inputSequenceDescription.get() is not None:
            self.description = self.inputSequenceDescription.get()
        elif hasattr(self,'description'):
            pass
        else:
            self.description = ''

        seq = Sequence(name=self.name,
                       sequence=self.sequence,
                       alphabet = self.alphabet,
                       isAminoacids=(self.inputSequence ==
                                     SEQ_TYPE_AMINOACIDS),
                       id=self.id, description=self.description)
        outputs = {'outputSequence': seq}
        self._defineOutputs(**outputs)

    def _summary(self):
        summary = []
        self.name = self.inputSequenceName.get()
        uniProtId = self._getUniProtID()
        geneBankID = self._getGeneBankID()
        if self.inputSequence == SEQ_TYPE_AMINOACIDS:
            summary.append('Sequence of aminoacids:\n')
            if self.inputProteinSequence == self.IMPORT_FROM_PLAIN_TEXT:
                summary.append("Sequence *%s* imported from plain text\n"
                               % self.name)
            elif self.inputProteinSequence == self.IMPORT_FROM_STRUCTURE:
                if self.inputStructureSequence == \
                    self.IMPORT_STRUCTURE_FROM_ID:
                    summary.append("Sequence *%s* imported from atomic "
                                   "structure *%s.cif*\n"
                                   % (self.name, self.pdbId.get()))
                elif self.inputStructureSequence == \
                    self.IMPORT_STRUCTURE_FROM_FILES:
                    summary.append("Sequence *%s* imported from file *%s*\n"
                                   % (self.name, self.pdbFile.get()))
            elif self.inputProteinSequence == self.IMPORT_FROM_UNIPROT:
                summary.append("Sequence *%s* imported from UniProt ID "
                               "*%s*\n"
                               % (self.name, uniProtId))
            elif self.inputProteinSequence == self.IMPORT_FROM_FILES:
                summary.append("Sequence *%s* imported from file name: "
                               "*%s*\n"
                               % (self.name, self.fileSequence.get()))
        else:
            summary.append('Sequence of nucleotides:\n')
            if self.inputNucleotideSequence == \
                    self.IMPORT_FROM_NUCLEOTIDE_PLAIN_TEXT:
                summary.append("Sequence *%s* imported from plain text\n"
                               % self.name)
            elif self.inputNucleotideSequence == \
                    self.IMPORT_FROM_NUCLEOTIDE_STRUCTURE:
                if self.inputStructureSequence == \
                    self.IMPORT_STRUCTURE_FROM_ID:
                    summary.append("Sequence *%s* imported from atomic "
                                   "structure *%s.cif*\n"
                                   % (self.name, self.pdbId.get()))
                elif self.inputStructureSequence == \
                    self.IMPORT_STRUCTURE_FROM_FILES:
                    summary.append("Sequence *%s* imported from file *%s*\n"
                                   % (self.name, self.pdbFile.get()))
            elif self.inputNucleotideSequence == self.IMPORT_FROM_GENEBANK:
                summary.append("Sequence *%s* imported from geneBank ID "
                               "*%s*\n"
                               % (self.name, geneBankID))
            elif self.inputNucleotideSequence == \
                    self.IMPORT_FROM_NUCLEOTIDE_FILES:
                summary.append("Sequence *%s* imported from file name: "
                               "*%s*\n"
                               % (self.name, self.fileSequence.get()))
        return summary

    def _validate(self):
        errors = []
        return errors

    def _getSequenceName(self):
        pass

    def _getUniProtID(self):
        return self.uniProtSequence.get()

    def _getGeneBankID(self):
        return self.geneBankSequence

    def _getAlphabet(self):
        if self.inputSequence == SEQ_TYPE_AMINOACIDS:
            return self.proteinIUPACalphabet.get()
        else:
            return self.nucleotideIUPACalphabet.get()

