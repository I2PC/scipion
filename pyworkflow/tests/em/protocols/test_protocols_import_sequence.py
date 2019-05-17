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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtImportSequence
from pyworkflow.em.convert.sequence import (
    SEQ_TYPE_NUCLEOTIDES, PROTEIN_ALPHABET, AMBIGOUS_RNA_ALPHABET,
    EXTENDED_DNA_ALPHABET, indexToAlphabet)


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportSequence(TestImportBase):
    USERID = "UserID"
    NAME = 'USER_SEQ'
    DESCRIPTION = 'User description'
    NUCLEOTIDESEQ1 = 'AATGCGGTTGGGBDSW********GGCACACG'
    AMINOACIDSSEQ1 = 'LARKJLAKPABXZJUO********VAVAVALK'
    CHAIN1 = '{"model": 0, "chain": "B", "residues": 148}' # Protein
    CHAIN2 = '{"model": 0, "chain": "A", "residues": 12}' # RNA
    CHAIN3 = '{"model": 0, "chain": "I", "residues": 146}' # DNA
    pdbID1 = "3lqd" # Protein
    pdbID2 = "205d" # RNA
    pdbID3 = "1aoi" # DNA and protein
    GENEBANKID = 'AJ520101.1'
    UNIPROTID = 'P12345'

    def testImportUserNucleotideSequence1(self):
        """
        Import a single nucleotide sequence provided by the user (nucleotide
        alphabet by default)
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'nucleotideIUPACalphabet': EXTENDED_DNA_ALPHABET,
                'inputRawSequence': self.NUCLEOTIDESEQ1
               }
        prot1 = self.newProtocol(ProtImportSequence, **args)
        prot1.setObjLabel('1_import DNA,\nseq from user\nExtended DNA '
                           'alphabet')
        self.launchProtocol(prot1)
        sequence = prot1.outputSequence
        self.assertEqual("USER_SEQ", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("AATGCGGTTG", sequence.getSequence()[:10])
        self.assertEqual("GATCBDSW",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportUserNucleotideSequence2(self):
        """
        Import a single nucleotide sequence provided by the user (nucleotide
        alphabet AMBIGOUS_RNA_ALPHABET)
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'nucleotideIUPACalphabet': AMBIGOUS_RNA_ALPHABET,
                'inputRawSequence': self.NUCLEOTIDESEQ1
               }
        prot2 = self.newProtocol(ProtImportSequence, **args)
        prot2.setObjLabel('2_import RNA,\nseq from user\nAmbigous RNA '
                           'alphabet')
        self.launchProtocol(prot2)
        sequence = prot2.outputSequence
        self.assertEqual("USER_SEQ", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("AAGCGGGGGB", sequence.getSequence()[:10])
        self.assertEqual("GAUCRYWSMKHBVDN",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportStructureNucleotideSequence1(self):
        """
        Import the sequence of chain A of atomic RNA structure 205d.cif
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_ID,
                'pdbId': self.pdbID2,
                'inputStructureChain': self.CHAIN2
                }
        prot3 = self.newProtocol(ProtImportSequence, **args)
        prot3.setObjLabel('3_import RNA seq,\n from atomic '
                           'structure')
        self.launchProtocol(prot3)
        sequence = prot3.outputSequence
        self.assertEqual("205d_A", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("GGACUUUGGU", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportStructureNucleotideSequence2(self):
        """
        Import the sequence of chain A of atomic RNA structure 205d.cif
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/205d.cif'),
                'inputStructureChain': self.CHAIN2
                }
        prot4 = self.newProtocol(ProtImportSequence, **args)
        prot4.setObjLabel('4_import RNA seq,\n from atomic '
                           'structure')
        self.launchProtocol(prot4)
        sequence = prot4.outputSequence
        self.assertEqual("UserID", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("GGACUUUGGU", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportStructureNucleotideSequence3(self):
        """
        Import the sequence of chain I of atomic DNA structure 1aoi.cif
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_ID,
                'pdbId': self.pdbID3,
                'inputStructureChain': self.CHAIN3
                }
        prot5 = self.newProtocol(ProtImportSequence, **args)
        prot5.setObjLabel('5_import DNA seq,\n from atomic '
                           'structure')
        self.launchProtocol(prot5)
        sequence = prot5.outputSequence
        self.assertEqual("1aoi_I", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("ATCAATATCC", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportStructureNucleotideSequence4(self):
        """
        Import the sequence of chain I of atomic DNA structure 1aoi.cif
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/1aoi.cif'),
                'inputStructureChain': self.CHAIN3
                }
        prot6 = self.newProtocol(ProtImportSequence, **args)
        prot6.setObjLabel('6_import DNA seq,\n from atomic '
                           'structure')
        self.launchProtocol(prot6)
        sequence = prot6.outputSequence
        self.assertEqual("UserID", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("ATCAATATCC", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportFileNucleotideSequence1(self):
        """
        Import a single nucleotide sequence from a text file
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/AJ520101.fasta')
                }
        prot7 = self.newProtocol(ProtImportSequence, **args)
        prot7.setObjLabel('7_import nucleotide,\nseq from file')
        self.launchProtocol(prot7)
        sequence = prot7.outputSequence
        self.assertEqual("AJ520101.1", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual("AJ520101.1 Rhizobium leguminosarum bv. viciae "
                         "plasmid "
                         "partial fixC gene, fixX gene, nifA gene "
                         "and nifB gene",
                        sequence.getDescription())
        self.assertEqual("GGATCCGAGA", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportFileNucleotideSequence2(self):
        """
        Import a single nucleotide sequence from a text file
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/AJ520101.fasta')
                }
        prot8 = self.newProtocol(ProtImportSequence, **args)
        prot8.setObjLabel('8_import nucleotide,\nseq from file')
        self.launchProtocol(prot8)
        sequence = prot8.outputSequence
        self.assertEqual("UserID", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("GGATCCGAGA", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportFileNucleotideSequence3(self):
        """
        Import a single nucleotide sequence from a text file
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/Several_sequences.fasta')
                }
        prot9 = self.newProtocol(ProtImportSequence, **args)
        prot9.setObjLabel('9_import nucleotide,\nseq from file')
        self.launchProtocol(prot9)
        sequence = prot9.outputSequence
        self.assertEqual("Seq1", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual("Seq1 This is the sequence number one",
                         sequence.getDescription())
        self.assertEqual("TGGCTAAATA", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportFileNucleotideSequence4(self):
        """
        Import a single nucleotide sequence from a text file that contains
        several sequences
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_NUCLEOTIDE_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/Several_sequences.fasta')
                }
        prot10 = self.newProtocol(ProtImportSequence, **args)
        prot10.setObjLabel('10_import nucleotide,\nseq from file')
        self.launchProtocol(prot10)
        sequence = prot10.outputSequence
        self.assertEqual("UserID", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("TGGCTAAATA", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportGeneBankNucleotideSequence1(self):
        """
        Import a single nucleotide sequence from GeneBank
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_GENEBANK,
                'geneBankSequence': self.GENEBANKID
                }
        prot11 = self.newProtocol(ProtImportSequence, **args)
        prot11.setObjLabel('11_import nucleotide,\nseq from '
                           'GeneBank')
        self.launchProtocol(prot11)
        sequence = prot11.outputSequence
        self.assertEqual("AJ520101.1", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual("AJ520101.1 Rhizobium leguminosarum bv. viciae "
                         "plasmid partial fixC gene, fixX gene, nifA gene "
                         "and nifB gene",
                         sequence.getDescription())
        self.assertEqual("GGATCCGAGA", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportGeneBankNucleotideSequence2(self):
        """
        Import a single nucleotide sequence from GeneBank
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputSequence': SEQ_TYPE_NUCLEOTIDES,
                'inputNucleotideSequence':
                    ProtImportSequence.IMPORT_FROM_GENEBANK,
                'geneBankSequence': self.GENEBANKID
                }
        prot12 = self.newProtocol(ProtImportSequence, **args)
        prot12.setObjLabel('12_import nucleotide,\nseq from '
                           'GeneBank')
        self.launchProtocol(prot12)
        sequence = prot12.outputSequence
        self.assertEqual("UserID", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual("User description",
                         sequence.getDescription())
        self.assertEqual("GGATCCGAGA", sequence.getSequence()[:10])
        self.assertEqual("GAUC",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportUserAminoacidSequence1(self):
        """
        Import a single aminoacid sequence provided by the user (protein
        alphabet by default)
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputRawSequence': self.AMINOACIDSSEQ1
                }
        prot13 = self.newProtocol(ProtImportSequence, **args)
        prot13.setObjLabel('13_import aminoacid seq,\n from user\nExtended '
                           'protein alphabet')
        self.launchProtocol(prot13)

        sequence = prot13.outputSequence
        self.assertEqual("USER_SEQ", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description", sequence.getDescription())
        self.assertEqual("LARKJLAKPA", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWYBXZJUO",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                        )

    def testImportUserAminoacidSequence2(self):
        """
        Import a single aminoacid sequence provided by the user (protein
        alphabet PROTEIN_ALPHABET
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'proteinIUPACalphabet': PROTEIN_ALPHABET,
                'inputRawSequence': self.AMINOACIDSSEQ1
                }
        prot14 = self.newProtocol(ProtImportSequence, **args)
        prot14.setObjLabel('14_import aminoacid seq,\n from user\nProtein '
                           'alphabet')
        self.launchProtocol(prot14)
        sequence = prot14.outputSequence
        self.assertEqual("USER_SEQ", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description", sequence.getDescription())
        self.assertEqual("LARKLAKPAV", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                        )

    def testImportStructureAminoacidSequence1(self):
        """
        Import the sequence of chain B of atomic structure 3lqd.cif
        """
        args = {'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_ID,
                'pdbId': self.pdbID1,
                'inputStructureChain': self.CHAIN1
                }
        prot15 = self.newProtocol(ProtImportSequence, **args)
        prot15.setObjLabel('15_import aminoacid seq,\n from atomic '
                           'structure')
        self.launchProtocol(prot15)
        sequence = prot15.outputSequence

        self.assertEqual("3lqd_B", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description", sequence.getDescription())
        self.assertEqual("VHLSGEEKSA", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportStructureAminoacidSequence2(self):
        """
        Import the sequence of chain B of atomic structure 3lqd.cif
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile('PDBx_mmCIF/3lqd.cif'),
                'inputStructureChain': self.CHAIN1
                }
        prot16 = self.newProtocol(ProtImportSequence, **args)
        prot16.setObjLabel('16_import aminoacid seq,\n from atomic '
                           'structure')
        self.launchProtocol(prot16)
        sequence = prot16.outputSequence
        self.assertEqual("UserID", sequence.getId())
        self.assertEqual("USER_SEQ", sequence.getSeqName())
        self.assertEqual("User description", sequence.getDescription())
        self.assertEqual("VHLSGEEKSA", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportFileAminoacidSequence1(self):
        """
        Import a single aminoacid sequence from a text file
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/COX1_human.fasta')
                }
        prot17 = self.newProtocol(ProtImportSequence, **args)
        prot17.setObjLabel('17_import aminoacid,\nseq from file')
        self.launchProtocol(prot17)
        sequence = prot17.outputSequence

        self.assertEqual("YP_003024028.1", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual('YP_003024028.1 cytochrome c oxidase '
                         'subunit I (mitochondrion) [Homo sapiens]',
                        sequence.getDescription())
        self.assertEqual("MFADRWLFST", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportFileAminoacidSequence2(self):
        """
        Import a single aminoacid sequence from a text file
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/COX1_human.fasta')
                }
        prot18 = self.newProtocol(ProtImportSequence, **args)
        prot18.setObjLabel('18_import aminoacid,\nseq from file')
        self.launchProtocol(prot18)
        sequence = prot18.outputSequence
        self.assertEqual("UserID", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual('User description',
                         sequence.getDescription())
        self.assertEqual("MFADRWLFST", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportUniprotAminoacidSequence1(self):
        """
        Import a single aminoacid sequence from UniProt
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID
                }
        prot19 = self.newProtocol(ProtImportSequence, **args)
        prot19.setObjLabel('19_import aminoacids,\nseq from '
                           'UniProt')
        self.launchProtocol(prot19)
        sequence = prot19.outputSequence
        self.assertEqual("P12345", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual('Aspartate aminotransferase, mitochondrial',
                         sequence.getDescription())
        self.assertEqual("MALLHSARVL", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)

    def testImportUniprotAminoacidSequence2(self):
        """
        Import a single aminoacid sequence from UniProt
        """
        args = {'inputSequenceID': self.USERID,
                'inputSequenceName': self.NAME,
                'inputSequenceDescription': self.DESCRIPTION,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID
                    }
        prot20 = self.newProtocol(ProtImportSequence, **args)
        prot20.setObjLabel('20_import aminoacids,\nseq from '
                            'UniProt')
        self.launchProtocol(prot20)
        sequence = prot20.outputSequence
        self.assertEqual("UserID", sequence.getId())
        self.assertEqual('USER_SEQ', sequence.getSeqName())
        self.assertEqual('User description',
                         sequence.getDescription())
        self.assertEqual("MALLHSARVL", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)



