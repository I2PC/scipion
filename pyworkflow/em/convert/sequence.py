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


# sequence related stuff
SEQ_TYPE=['aminoacids', 'nucleotides']
SEQ_TYPE_AMINOACIDS = 0
SEQ_TYPE_NUCLEOTIDES = 1
IUPAC_PROTEIN_ALPHABET = ['Extended Protein', 'Protein']
EXTENDED_PROTEIN_ALPHABET = 0
PROTEIN_ALPHABET = 1
IUPAC_NUCLEOTIDE_ALPHABET = ['Ambiguous DNA', 'Unambiguous DNA',
                             'Extended DNA', 'Ambiguous RNA',
                             'Unambiguous RNA']
EXTENDED_DNA_ALPHABET = 0
AMBIGOUS_DNA_ALPHABET = 1
UNAMBIGOUS_DNA_ALPHABET = 2
AMBIGOUS_RNA_ALPHABET = 3
UNAMBIGOUS_RNA_ALPHABET = 4


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Entrez, SeqIO
import urllib, urllib2, sys
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline
from Bio import pairwise2


class SequenceHandler:
    def __init__(self, sequence=None,
                 iUPACAlphabet=0,
                 isAminoacid=True):

        self.isAminoacid = isAminoacid
        self.alphabet = indexToAlphabet(isAminoacid, iUPACAlphabet)

        if sequence is not None:
            #sequence = cleanSequence(self.alphabet, sequence)
            self._sequence = Seq(sequence, alphabet=self.alphabet)
        else:
            self._sequence = None
            # type(self._sequence):  <class 'Bio.Seq.Seq'>

    def saveFile(self, fileName, seqID, sequence=None, name=None,
                 seqDescription=None, type="fasta"):
        if sequence is not None:
            self._sequence = sequence
        record = SeqRecord(self._sequence, id=seqID, name=name,
                           description=seqDescription)
        # type(record): < class 'Bio.SeqRecord.SeqRecord'>
        with open(fileName, "w") as output_handle:
            SeqIO.write(record, output_handle, type)

    def downloadSeqFromFile(self, fileName, type="fasta"):
        record = next(SeqIO.parse(fileName, type))
        return record

    def downloadSeqFromDatabase(self, seqID):
        # see http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
        # for format/databases
        print "Connecting to dabase..."
        seqID = str(seqID)
        sys.stdout.flush()
        counter=1
        retries = 5
        record = None
        error = ""
        while counter <= retries:  # retry up to 5 times if server busy
            try:
                if self.isAminoacid:
                    dataBase = 'UnitProt'
                    url = "http://www.uniprot.org/uniprot/%s.xml"
                    format = "uniprot-xml"
                    handle = urllib2.urlopen(url % seqID)
                    print "URL", url % seqID
                else:
                    dataBase = 'GeneBank'
                    Entrez.email = "adam.richards@stat.duke.edu"
                    format = "fasta"
                    handle = Entrez.efetch(db="nucleotide", id=seqID,
                                           rettype=format, retmode="text")

                record = SeqIO.read(handle, format)
                break
            except urllib2.HTTPError, e:
                error = "%s is a wrong sequence ID" % seqID
                print e.code
            except urllib2.URLError, e:
                error = "Cannot connect to %s" % dataBase
                print e.args
            except Exception as ex:
                template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                error = message
            if counter == retries:
                break
            counter += 1
        return record, error

    def alignSeq(self, referenceSeq):
        if self._sequence is not None:
            alignments = pairwise2.align.globalds(self._sequence.seq,
                                                  referenceSeq.seq)
            return alignments
        else:
            print "read the sequence first"
            exit(0)

def cleanSequenceScipion(isAminoacid, iUPACAlphabet, sequence):
    return cleanSequence(indexToAlphabet(isAminoacid, iUPACAlphabet), sequence)

def cleanSequence(alphabet, sequence):
    str_list = []
    for item in sequence.upper():
        if item in alphabet.letters:
            str_list.append(item)
    value =  ''.join(str_list)
    return ''.join(str_list)

def indexToAlphabet(isAminoacid, iUPACAlphabet):
    if isAminoacid:
        if iUPACAlphabet == EXTENDED_PROTEIN_ALPHABET:
            alphabet = IUPAC.ExtendedIUPACProtein
        else:
            alphabet = IUPAC.IUPACProtein
    else:
        if iUPACAlphabet == EXTENDED_DNA_ALPHABET:
            alphabet = IUPAC.ExtendedIUPACDNA
        elif iUPACAlphabet == AMBIGOUS_DNA_ALPHABET:
            alphabet = IUPAC.IUPACAmbiguousDNA
        elif iUPACAlphabet == UNAMBIGOUS_DNA_ALPHABET:
            alphabet = IUPAC.IUPACUnambiguousDNA
        elif iUPACAlphabet == AMBIGOUS_RNA_ALPHABET:
            alphabet = IUPAC.IUPACAmbiguousRNA
        else:
            alphabet = IUPAC.IUPACUnambiguousRNA
    return alphabet

def alphabetToIndex(isAminoacid, alphabet):
    if isAminoacid:
        if type(alphabet) is type(IUPAC.ExtendedIUPACProtein):
            return EXTENDED_PROTEIN_ALPHABET
        else:
            return PROTEIN_ALPHABET
    else:
        if type(alphabet) is type(IUPAC.ExtendedIUPACDNA):
            return EXTENDED_DNA_ALPHABET
        elif type(alphabet) is type(IUPAC.IUPACAmbiguousDNA):
            return AMBIGOUS_DNA_ALPHABET
        elif type(alphabet) is type(IUPAC.IUPACUnambiguousDNA):
            return UNAMBIGOUS_DNA_ALPHABET
        elif type(alphabet) is type(IUPAC.IUPACAmbiguousRNA):
            return AMBIGOUS_RNA_ALPHABET
        else:
            return UNAMBIGOUS_RNA_ALPHABET

def saveFileSequencesToAlign(SeqDic, inFile, type="fasta"):
    # Write my sequences to a fasta file
    with open(inFile, "w") as output_handle:
        for index, seq in SeqDic.iteritems():
            record = SeqRecord(seq, id=str(index),
                           name="", description="")
            SeqIO.write(record, output_handle, type)

def alignClustalSequences(inFile, outFile):
    # Alignment of sequences with Clustal Omega program
    clustalomega_cline = ClustalOmegaCommandline(
            infile=inFile,
            outfile=outFile,
            verbose=True, auto=True)
    return clustalomega_cline

def alignMuscleSequences(inFile, outFile):
    # Alignment of sequences with Muscle program
    muscle_cline = MuscleCommandline(input=inFile, out=outFile)
    return muscle_cline

def alignBioPairwise2Sequences(structureSequenceId, structureSequence,
              userSequenceId, userSequence,
              outFileName):
    "aligns two sequences and saves them to disk using fasta format"
    # see alignment_function for globalms parameters
    alignments = pairwise2.align.globalms(structureSequence,
                                           userSequence,  3, -1, -3, -2)
    align1, align2, score, begin, end = alignments[0]
    with open(outFileName, "w") as handle:
        handle.write(">%s\n%s\n>%s\n%s\n" % (structureSequenceId,
                                             align1,
                                             userSequenceId,
                                             align2))




