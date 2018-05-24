# **************************************************************************
# *
# * Authors:     roberto Marabini (roberto@cnb.csic.es)
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
#
from pyworkflow.tests import *
from tempfile import NamedTemporaryFile
from pyworkflow.em.convert_atom_struct import AtomicStructHandler
from pyworkflow.em.transformations import euler_matrix, \
    translation_matrix, concatenate_matrices
from copy import deepcopy
import numpy


class TestAtomicStructHandler(unittest.TestCase):

    def testIntToChain(self):
        aSH = AtomicStructHandler(self.PDBFileName)
        solString = ['A',  'B',  'C',  'D',  'E',  'F',  'G',
                     'H',  'I',  'J',  'K',  'L',  'M',  'N',
                     'O',  'P',  'Q',  'R',  'S',  'T',  'U',
                     'V',  'W',  'X',  'Y',  'Z',  '0',  '1',
                     '2',  '3',  '4',  '5',  '6',  '7',  '8',
                     '9',  'a',  'b',  'c',  'd',  'e',  'f',
                     'g',  'h',  'i',  'j',  'k',  'l',  'm',
                     'n',  'o',  'p',  'q',  'r',  's',  't',
                     'u',  'v',  'w',  'x',  'y',  'z', "AA"]
        for i in range(63):
            self.assertEqual(solString[i], aSH._intToChain(i))

    def testReadPDB(self):
        aSH = AtomicStructHandler(self.PDBFileName)
        structure = aSH.getStructure()

        solDict = {}
        solDict['structure_method'] = 'x-ray diffraction'
        solDict['head'] = 'extracellular matrix'
        solDict['name'] = 'x-ray crystallographic determination ' \
                          'of a collagen-like peptide with the ' \
                          'repeating sequence (pro-pro-gly)'
        solDict['author'] = 'R.Z.Kramer,L.Vitagliano,J.Bella,R.Berisio,' \
                            'L.Mazzarella,B.Brodsky,A.Zagari,H.M.Berman'
        solDict['deposition_date'] = '1998-01-22'
        for k, v in solDict.iteritems():
            self.assertEqual(structure.header[k].strip(), v)

        solList = ['N', 'CA', 'C', 'O', 'CB']
        counter = 0
        for atom in structure.get_atoms():
            self.assertEqual(solList[counter], atom.get_name())
            counter += 1

    def testReadCIF(self):
        aSH = AtomicStructHandler(self.CIFFileName)
        structure = aSH.getStructure()

        solList = ['N', 'CA', 'C', 'O', 'CB']
        counter = 0
        for atom in structure.get_atoms():
            self.assertEqual(solList[counter], atom.get_name())
            counter += 1

        solDict = {}
        solDict['_exptl.method'] = 'x-ray diffraction'
        solDict['_struct_keywords.pdbx_keywords'] = 'extracellular matrix'
        solDict['_struct.title'] = 'x-ray crystallographic determination ' \
                                   'of a collagen-like peptide with the ' \
                                   'repeating sequence (pro-pro-gly)'
        _dict = aSH.readLowLevel(self.CIFFileName)

        for k, v in solDict.iteritems():
            self.assertEqual(_dict[k].strip().lower(), v.lower())

    def testRenameToChains(self):
        aSH = AtomicStructHandler(self.PDBFileName)
        structure = aSH.getStructure()

        model = structure[0]
        chain = model['B']
        chain.id = 'CC'
        aSH.renameChains(structure)
        for chain in structure.get_chains():
            self.assertEqual(chain.id, 'C')

    def testWritePDB(self):
        aSH = AtomicStructHandler(self.PDBFileName)
        PDBFileName2 = self.PDBFileName.replace(".pdb", "_2.pdb")
        aSH._write(PDBFileName2)
        aSH2 = AtomicStructHandler(PDBFileName2)
        structure1 = aSH.getStructure()
        structure2 = aSH2.getStructure()

        for atom1, atom2 in zip(structure1.get_atoms(),
                                structure2.get_atoms()):
            self.assertEqual(atom1.get_name(), atom2.get_name())
        os.unlink(PDBFileName2)

    def testWriteCIF(self):
        aSH = AtomicStructHandler(self.PDBFileName)
        CIFFileName2 = self.CIFFileName.replace(".cif", "_2.cif")
        aSH._write(CIFFileName2)
        aSH2 = AtomicStructHandler(CIFFileName2)
        structure1 = aSH.getStructure()
        structure2 = aSH2.getStructure()

        for atom1, atom2 in zip(structure1.get_atoms(),
                                structure2.get_atoms()):
            self.assertEqual(atom1.get_name(), atom2.get_name())
        os.unlink(CIFFileName2)

    def testCIFToPDB(self):
        aSH = AtomicStructHandler(self.CIFFileName)
        structure1 = aSH.getStructure()
        PDBFileName2 = self.PDBFileName.replace(".pdb", "_2.pdb")
        aSH.write(PDBFileName2)
        aSH.read(PDBFileName2)
        structure2 = aSH.getStructure()
        for atom1, atom2 in zip(structure1.get_atoms(),
                                structure2.get_atoms()):
            self.assertEqual(atom1.get_name(), atom2.get_name())
        os.unlink(PDBFileName2)

    def testPDBToCIF(self):
        aSH1 = AtomicStructHandler(self.PDBFileName)
        CIFFileName2 = self.CIFFileName.replace(".cif", "_2.cif")
        aSH1.write(CIFFileName2)
        aSH2 = AtomicStructHandler(CIFFileName2)

        structure1 = aSH1.getStructure()
        structure2 = aSH2.getStructure()
        for atom1, atom2 in zip(structure1.get_atoms(),
                                structure2.get_atoms()):
            self.assertEqual(atom1.get_name(), atom2.get_name())
        os.unlink(CIFFileName2)

    def testCenterOfMass(self):
        aSH = AtomicStructHandler()
        structure = aSH.read(self.CIFFileName)
        x, y, z = aSH.centerOfMass(structure)
        self.assertAlmostEqual(x,  8.1593891597, 2)
        self.assertAlmostEqual(y, 21.1833304818, 2)
        self.assertAlmostEqual(z, 20.0177924407, 2)

    def testTransformTranslation(self):
        aSHSource = AtomicStructHandler(self.PDBFileName)
        structure = aSHSource.getStructure()
        structure_copy = deepcopy(aSHSource.getStructure())
        shift = [100., 50., 25.]
        #        rotation_matrix = euler_matrix(deg2rad(45.), 0., 0., 'szyz')
        rotation_matrix = euler_matrix(0., 0., 0., 'szyz')
        translation = translation_matrix(shift)
        M = concatenate_matrices(rotation_matrix, translation)
        aSHSource.transform(M)
        for atom1, atom2 in zip(structure.get_atoms(),
                                structure_copy.get_atoms()):
            coord1 = atom1.get_coord()
            coord2 = [sum(x) for x in zip(atom2.get_coord(), shift)]
            for i in range(3):
                self.assertAlmostEqual(coord1[i], coord2[i], 2)

    def testTransformRotation(self):
        aSHSource = AtomicStructHandler(self.PDBFileName)
        structure = aSHSource.getStructure()
        structure_copy = deepcopy(structure)
        rot = numpy.deg2rad(10)
        theta = numpy.deg2rad(20.)
        psi = numpy.deg2rad(30.)
        rotation_matrix = euler_matrix(rot, theta, psi, 'szyz')
        translation = translation_matrix([0., 0., 0.])
        M = concatenate_matrices(rotation_matrix, translation)
        aSHSource.transform(M)
        m = M[:3, :3]
        for atom1, atom2 in zip(structure.get_atoms(),
                                structure_copy.get_atoms()):
            coord1 = atom1.get_coord()
            coord2 = m.dot(atom2.get_coord())
            for i in range(3):
                self.assertAlmostEqual(coord1[i], coord2[i], 2)

    def testTransformRotationAndTranslation(self):
        aSHSource = AtomicStructHandler(self.PDBFileName)
        structure = aSHSource.getStructure()
        structure_copy = deepcopy(structure)
        rot = numpy.deg2rad(10)
        theta = numpy.deg2rad(20.)
        psi = numpy.deg2rad(30.)
        rotation_matrix = euler_matrix(rot, theta, psi, 'szyz')
        shift = [100., 50., 25.]
        translation = translation_matrix(shift)
        M = concatenate_matrices(rotation_matrix, translation)
        aSHSource.transform(M)
        m = M[:3, :3]
        for atom1, atom2 in zip(structure.get_atoms(),
                                structure_copy.get_atoms()):
            coord1 = atom1.get_coord()
            coord2 = atom2.get_coord()
            coord2 = [sum(x) for x in zip(coord2, shift)]
            coord2 = m.dot(coord2)
            for i in range(3):
                self.assertAlmostEqual(coord1[i], coord2[i], 2)

    def testTransformTranslationCoherence(self):
        """
        Question: If I transform the PDB and the 3D map with matrix T
        do they move in the same direction?
        I do not know how to make an automatic test to check this
        The following code perform all operations but check that
        the PDB and 3D map match. This should be check using
        your eye.

        change False to True (in the following if) to
        perform all or some of the checks
        """

        # retrieve "Structure of the human TRPC3
        # both 3Dmap and PDB


        doTest = False

        if not doTest:

            print "This test is to be tested manually since it opens chimera afterwards"
            print "For testing this, edit this file and set doTest = True"
            return


        PDBID = '6CUD'
        EMDBID = '7620'

        doAll = False

        if False or doAll:  # set to False if you aready have the 3dmap file
            url = 'ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-%s/map/emd_%s.map.gz' % \
                  (EMDBID, EMDBID)
            import urllib
            urllib.urlretrieve(url, 'emd_%s.map.gz' % EMDBID)
            os.system("gunzip emd_%s.map.gz" % EMDBID)  # file is gzipped
        if False or doAll:  # set to False if you aready have the PDB file
            aSH = AtomicStructHandler()
            pdbFileName = aSH.readFromPDBDatabase(PDBID, type='pdb',
                                                  dir=os.getcwd())
        else:
            pdbFileName = 'pdb%s.ent' % PDBID.lower()

        # get 3D map sampling
        from pyworkflow.em.convert_header.CCP4.convert import Ccp4Header
        header = Ccp4Header("emd_%s.map" % EMDBID, readHeader=True)
        sampling, y, z = header.getSampling()

        def __runXmippProgram(program, args):
            """ Internal function to launch a Xmipp program. """
            import pyworkflow.em.packages.xmipp3 as xmipp3
            xmipp3.runXmippProgram(program, args)

        def __getXmippEulerAngles(matrix):
            """ Internal fuction to convert scipion to xmipp angles"""
            from pyworkflow.em.packages.xmipp3.convert \
                import geometryFromMatrix

            return geometryFromMatrix(matrix, False)

        def __applyTransform(suffix, pdbFileName, shift, angles, sampling):
            """ auxiliary function, transform PDB and 3dmap files"""
            # create a Scipion transformation matrix
            from numpy import deg2rad
            rotation_matrix = euler_matrix(deg2rad(angles[0]),
                                           deg2rad(angles[1]),
                                           deg2rad(angles[2]), 'szyz')
            translation = translation_matrix(shift)
            M = concatenate_matrices(rotation_matrix, translation)

            # apply it to the pdb file
            # if rotation move to center
            aSH = AtomicStructHandler(pdbFileName)
            if (angles[0] != 0. or angles[1] != 0. or angles[2] != 0.):
                from pyworkflow.em.convert import ImageHandler
                ih = ImageHandler()
                x, y, z, n = ih.getDimensions("emd_%s.map" % EMDBID)
                x /= 2.
                y /= 2.
                z /= 2.
                localShift = [-x, -y, -z]
                rotation_matrix = euler_matrix(0., 0., 0., 'szyz')
                translation = translation_matrix(localShift)
                localM = concatenate_matrices(rotation_matrix, translation)
                aSH.transform(localM, sampling=sampling)

            aSH.transform(M, sampling=sampling)

            if (angles[0] != 0. or angles[1] != 0. or angles[2] != 0.):
                localShift = [x, y, z]
                rotation_matrix = euler_matrix(0., 0., 0., 'szyz')
                translation = translation_matrix(localShift)
                localM = concatenate_matrices(rotation_matrix, translation)
                aSH.transform(localM, sampling=sampling)

            aSH.write("%s_%s_transformed.ent" % (suffix, PDBID.lower()))

            # get equivalent xmipp transformation
            shift, angles = __getXmippEulerAngles(M)
            # shift 3map and set sampling
            __runXmippProgram("xmipp_transform_geometry",
                              '-i emd_%s.map '
                              '-o %s_emd_%s_transform.map '
                              '--interp linear '
                              '--shift %f %f %f '
                              '--rotate_volume euler %f %f %f ' % (
                                  EMDBID,
                                  suffix,
                                  EMDBID,
                                  shift[0], shift[1], shift[2],
                                  angles[0], angles[1], angles[2]
                              )
                              )
            header = Ccp4Header("%s_emd_%s_transform.map" % (suffix, EMDBID),
                                readHeader=True)
            header.setSampling(sampling)
            # put the sampling back, xmipp_transform_geometry erased it
            header.writeHeader()

            # view the results with chimera
            from pyworkflow.em.viewers.chimera_utils \
                import runChimeraProgram
            from pyworkflow.em.viewers.chimera_utils \
                import getProgram as chimera_get_program
            args = "%s %s %s %s" % (
                   pdbFileName,
                   "emd_%s.map" % EMDBID,
                   "%s_%s_transformed.ent" % (suffix, PDBID.lower()),
                   "%s_emd_%s_transform.map" % (suffix, EMDBID)
            )
            runChimeraProgram(chimera_get_program(), args)

        # shift atomic structure
        doAll = True
        if False or doAll:
            shift = [20., 0., 0.]
            angles = [0., 0., 0.]
            __applyTransform("Xshift", pdbFileName, shift, angles, sampling)

        # repeat test this time  rotation one angle
        # problem, xmipp rotates with respect the volume center
        # pdb with respect the origin of coordinates (much better convention)
        # in order to compare both I need to
        # move pdb to origin, rotate it, put it back in the possition
        if False or doAll:
            shift = [0., 0., 0.]
            angles = [30., 0., 0.]
            __applyTransform("Rot2D", pdbFileName, shift, angles, sampling)

        # repeat test this time  rotation in 3 angles
        # problem, xmipp rotates with respect the volume center
        # pdb with respect the origin of coordinates (much better convention)
        if False or doAll:
            shift = [0., 0., 0.]
            angles = [10., 20., 30.]
            __applyTransform("Rot3D", pdbFileName, shift, angles, sampling)

        # repeat test this time  rotation in 3 angles and shift
        # problem, xmipp rotates with respect the volume center
        # pdb with respect the origin of coordinates (much better convention)
        if False or doAll:
            shift = [5., 10., 15.]
            angles = [10., 20., 30.]
            __applyTransform("Rot3DShift", pdbFileName,
                             shift, angles, sampling)

    def testReadFromPDBDatabase(self):
        PDBID = '6CUD'
        aSH = AtomicStructHandler()
        # EMD-7620
        fileName = aSH.readFromPDBDatabase(PDBID, type='pdb', dir='/tmp')
        self.assertTrue(os.path.exists(fileName))

        os.unlink(fileName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.PDBFileName):
            os.unlink(cls.PDBFileName)
        if os.path.exists(cls.CIFFileName):
            os.unlink(cls.CIFFileName)

    @classmethod
    def setUpClass(cls):
        PDBString = """HEADER    EXTRACELLULAR MATRIX                    22-JAN-98   1A3I
TITLE     X-RAY CRYSTALLOGRAPHIC DETERMINATION OF A COLLAGEN-LIKE
TITLE    2 PEPTIDE WITH THE REPEATING SEQUENCE (PRO-PRO-GLY)
EXPDTA    X-RAY DIFFRACTION
AUTHOR    R.Z.KRAMER,L.VITAGLIANO,J.BELLA,R.BERISIO,L.MAZZARELLA,
AUTHOR   2 B.BRODSKY,A.ZAGARI,H.M.BERMAN
REMARK 350 BIOMOLECULE: 1
REMARK 350 APPLY THE FOLLOWING TO CHAINS: B
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
ATOM      1  N   PRO B   1       8.316  21.206  21.530  1.00 17.44           N
ATOM      2  CA  PRO B   1       7.608  20.729  20.336  1.00 17.44           C
ATOM      3  C   PRO B   1       8.487  20.707  19.092  1.00 17.44           C
ATOM      4  O   PRO B   1       9.466  21.457  19.005  1.00 17.44           O
ATOM      5  CB  PRO B   1       6.460  21.723  20.211  1.00 22.26           C
END"""
        f = NamedTemporaryFile(delete=False, suffix=".pdb")
        f.write(PDBString)
        f.close()
        cls.PDBFileName = f.name
        CIFString = """data_1A3I
#
loop_
_audit_author.name
_audit_author.pdbx_ordinal
'Kramer, R.Z.'   1
'Vitagliano, L.' 2
'Bella, J.'      3
'Berisio, R.'    4
'Mazzarella, L.' 5
'Brodsky, B.'    6
'Zagari, A.'     7
'Berman, H.M.'   8
#
_database.entry_id        1A3I
_database.pdbx_code_NDB   1A3I
_database.pdbx_code_PDB   1A3I
_database.code_CSD        ?
#
loop_
_database_2.database_id
_database_2.database_code
PDB  1A3I
NDB  1A3I
RCSB 1A3I
#
_struct_keywords.entry_id        1A3I
_struct_keywords.pdbx_keywords   'EXTRACELLULAR MATRIX'
_struct_keywords.text            ?
#data_1A3I
#
loop_
_audit_author.name
_audit_author.pdbx_ordinal
'Kramer, R.Z.'   1
'Vitagliano, L.' 2
'Bella, J.'      3
'Berisio, R.'    4
'Mazzarella, L.' 5
'Brodsky, B.'    6
'Zagari, A.'     7
'Berman, H.M.'   8
#
_database.entry_id        1A3I
_database.pdbx_code_NDB   1A3I
_database.pdbx_code_PDB   1A3I
_database.code_CSD        ?
#
loop_
_database_2.database_id
_database_2.database_code
PDB  1A3I
NDB  1A3I
RCSB 1A3I
#
_struct_keywords.entry_id        1A3I
_struct_keywords.pdbx_keywords   'EXTRACELLULAR MATRIX'
_struct_keywords.text            ?
#
_exptl.entry_id          1A3I
_exptl.method            'X-RAY DIFFRACTION'
_exptl.crystals_number   ?
#
_struct.entry_id                  1A3I
_struct.title                     'X-RAY CRYSTALLOGRAPHIC DETERMINATION OF A COLLAGEN-LIKE PEPTIDE WITH THE REPEATING SEQUENCE (PRO-PRO-GLY)'
_struct.pdbx_descriptor           'EXTRACELLULAR MATRIX'
_struct.pdbx_model_details        ?
_struct.pdbx_details              ?
_struct.pdbx_CASP_flag            ?
_struct.pdbx_model_type_details   ?
#
_struct_biol_gen.biol_id                        1
_struct_biol_gen.asym_id                        A
_struct_biol_gen.pdbx_new_asym_id               A
_struct_biol_gen.pdbx_new_pdb_asym_id           B
_struct_biol_gen.symmetry                       1_555
_struct_biol_gen.pdbx_before_begin_residue_no   1
_struct_biol_gen.pdbx_before_end_residue_no     1
_struct_biol_gen.pdbx_after_begin_residue_no    1
_struct_biol_gen.pdbx_after_end_residue_no      1
_struct_biol_gen.details                        ?
_struct_biol_gen.pdbx_full_symmetry_operation   x,y,z
_struct_biol_gen.pdbx_PDB_order                 1
#
_struct_biol_view.biol_id            ?
_struct_biol_view.id                 ?
_struct_biol_view.rot_matrix[1][1]   ?
_struct_biol_view.rot_matrix[1][2]   ?
_struct_biol_view.rot_matrix[1][3]   ?
_struct_biol_view.rot_matrix[2][1]   ?
_struct_biol_view.rot_matrix[2][2]   ?
_struct_biol_view.rot_matrix[2][3]   ?
_struct_biol_view.rot_matrix[3][1]   ?
_struct_biol_view.rot_matrix[3][2]   ?
_struct_biol_view.rot_matrix[3][3]   ?
_struct_biol_view.pdbx_vector[1]     ?
_struct_biol_view.pdbx_vector[2]     ?
_struct_biol_view.pdbx_vector[3]     ?
#
_struct_asym.id                                              A
_struct_asym.pdbx_PDB_id                                     B
_struct_asym.pdbx_alt_id                                     B
_struct_asym.pdbx_blank_PDB_chainid_flag                     N
_struct_asym.pdbx_type                                       HETAIN
_struct_asym.pdbx_order                                      1
_struct_asym.pdbx_modified                                   N
_struct_asym.pdbx_fraction_per_asym_unit                     ?
_struct_asym.entity_id                                       1
_struct_asym.pdbx_missing_num_begin_of_chain_not_in_seqres   ?
_struct_asym.pdbx_missing_num_end_of_chain_not_in_seqres     ?
_struct_asym.pdbx_missing_num_begin_of_chain_in_seqres       ?
_struct_asym.details                                         ?
#
_pdbx_inhibitor_info.id                  1
_pdbx_inhibitor_info.name                PROLINE
_pdbx_inhibitor_info.num_per_asym_unit   1
#
_pdbx_nonpoly_scheme.asym_id         A
_pdbx_nonpoly_scheme.entity_id       1
_pdbx_nonpoly_scheme.mon_id          PRO
_pdbx_nonpoly_scheme.ndb_seq_num     1
_pdbx_nonpoly_scheme.pdb_seq_num     1
_pdbx_nonpoly_scheme.auth_seq_num    1
_pdbx_nonpoly_scheme.pdb_mon_id      PRO
_pdbx_nonpoly_scheme.auth_mon_id     PRO
_pdbx_nonpoly_scheme.pdb_strand_id   B
_pdbx_nonpoly_scheme.pdb_ins_code    .
#
_entity.id                         1
_entity.type                       non-polymer
_entity.src_method                 syn
_entity.pdbx_description           PROLINE
_entity.formula_weight             115.130
_entity.pdbx_number_of_molecules   1
_entity.details                    ?
#
_entity_keywords.entity_id               1
_entity_keywords.text                    ?
_entity_keywords.pdbx_mutation           ?
_entity_keywords.pdbx_fragment           ?
_entity_keywords.pdbx_antibody_isotype   ?
_entity_keywords.pdbx_ec                 ?
#
_entity_name_com.entity_id   1
_entity_name_com.name        ?
#
_entity_name_sys.entity_id   1
_entity_name_sys.name        ?
#
_chem_comp.id               PRO
_chem_comp.type             'L-peptide linking'
_chem_comp.mon_nstd_flag    y
_chem_comp.name             PROLINE
_chem_comp.pdbx_synonyms    ?
_chem_comp.formula          'C5 H9 N O2'
_chem_comp.formula_weight   115.130
#
_pdbx_database_status.status_code                        ?
_pdbx_database_status.entry_id                           1A3I
_pdbx_database_status.ndb_tid                            ?
_pdbx_database_status.status_coordinates_in_NDB          ?
_pdbx_database_status.recvd_deposit_form                 ?
_pdbx_database_status.date_deposition_form               ?
_pdbx_database_status.recvd_coordinates                  ?
_pdbx_database_status.date_coordinates                   ?
_pdbx_database_status.recvd_struct_fact                  ?
_pdbx_database_status.date_struct_fact                   ?
_pdbx_database_status.recvd_internal_approval            ?
_pdbx_database_status.recvd_nmr_constraints              ?
_pdbx_database_status.date_nmr_constraints               ?
_pdbx_database_status.recvd_manuscript                   ?
_pdbx_database_status.date_manuscript                    ?
_pdbx_database_status.name_depositor                     ?
_pdbx_database_status.rcsb_annotator                     ?
_pdbx_database_status.recvd_author_approval              ?
_pdbx_database_status.date_author_approval               ?
_pdbx_database_status.recvd_initial_deposition_date      ?
_pdbx_database_status.date_submitted                     ?
_pdbx_database_status.author_approval_type               ?
_pdbx_database_status.author_release_status_code         ?
_pdbx_database_status.date_revised                       ?
_pdbx_database_status.revision_id                        ?
_pdbx_database_status.replaced_entry_id                  ?
_pdbx_database_status.revision_description               ?
_pdbx_database_status.date_of_NDB_release                ?
_pdbx_database_status.date_released_to_PDB               ?
_pdbx_database_status.date_of_PDB_release                ?
_pdbx_database_status.date_hold_coordinates              ?
_pdbx_database_status.date_hold_struct_fact              ?
_pdbx_database_status.hold_for_publication               ?
_pdbx_database_status.date_hold_nmr_constraints          ?
_pdbx_database_status.dep_release_code_coordinates       ?
_pdbx_database_status.dep_release_code_struct_fact       ?
_pdbx_database_status.dep_release_code_nmr_constraints   ?
_pdbx_database_status.dep_release_code_sequence          ?
_pdbx_database_status.pdb_date_of_author_approval        ?
_pdbx_database_status.deposit_site                       ?
_pdbx_database_status.process_site                       ?
_pdbx_database_status.skip_PDB_REMARK                    ?
#
_pdbx_database_proc.entry_id           1A3I
_pdbx_database_proc.cycle_id           ?
_pdbx_database_proc.date_begin_cycle   ?
_pdbx_database_proc.date_end_cycle     ?
_pdbx_database_proc.details            ?
#
_database_PDB_remark.id     350
_database_PDB_remark.text
;BIOMOLECULE: 1
APPLY THE FOLLOWING TO CHAINS: B
  BIOMT1   1  1.000000  0.000000  0.000000        0.00000
  BIOMT2   1  0.000000  1.000000  0.000000        0.00000
;
#
_pdbx_coord.entry_id                1A3I
_pdbx_coord.chain_atoms_Y_P         ?
_pdbx_coord.hydrogen_atoms_Y_N      ?
_pdbx_coord.solvent_atoms_Y_N       ?
_pdbx_coord.structure_factors_Y_N   ?
#
_exptl_crystal.id                    1
_exptl_crystal.density_meas          ?
_exptl_crystal.density_Matthews      ?
_exptl_crystal.density_percent_sol   ?
_exptl_crystal.description           ?
#
_symmetry.entry_id                         1A3I
_symmetry.space_group_name_H-M             'P 1'
_symmetry.pdbx_full_space_group_name_H-M   ?
_symmetry.cell_setting                     ?
_symmetry.Int_Tables_number                ?
#
_cell.entry_id           1A3I
_cell.length_a           1.000
_cell.length_b           1.000
_cell.length_c           1.000
_cell.angle_alpha        90.00
_cell.angle_beta         90.00
_cell.angle_gamma        90.00
_cell.Z_PDB              1
_cell.pdbx_unique_axis   ?
#
_database_PDB_matrix.entry_id          1A3I
_database_PDB_matrix.origx[1][1]       1.000000
_database_PDB_matrix.origx[1][2]       0.000000
_database_PDB_matrix.origx[1][3]       0.000000
_database_PDB_matrix.origx[2][1]       0.000000
_database_PDB_matrix.origx[2][2]       1.000000
_database_PDB_matrix.origx[2][3]       0.000000
_database_PDB_matrix.origx[3][1]       0.000000
_database_PDB_matrix.origx[3][2]       0.000000
_database_PDB_matrix.origx[3][3]       1.000000
_database_PDB_matrix.origx_vector[1]   0.00000
_database_PDB_matrix.origx_vector[2]   0.00000
_database_PDB_matrix.origx_vector[3]   0.00000
#
_atom_sites.entry_id                    1A3I
_atom_sites.cartn_transform_axes        ?
_atom_sites.fract_transf_matrix[1][1]   1.000000
_atom_sites.fract_transf_matrix[1][2]   0.000000
_atom_sites.fract_transf_matrix[1][3]   0.000000
_atom_sites.fract_transf_matrix[2][1]   0.000000
_atom_sites.fract_transf_matrix[2][2]   1.000000
_atom_sites.fract_transf_matrix[2][3]   0.000000
_atom_sites.fract_transf_matrix[3][1]   0.000000
_atom_sites.fract_transf_matrix[3][2]   0.000000
_atom_sites.fract_transf_matrix[3][3]   1.000000
_atom_sites.fract_transf_vector[1]      0.00000
_atom_sites.fract_transf_vector[2]      0.00000
_atom_sites.fract_transf_vector[3]      0.00000
#
_struct_ncs_oper.id             ?
_struct_ncs_oper.code           ?
_struct_ncs_oper.details        ?
_struct_ncs_oper.matrix[1][1]   ?
_struct_ncs_oper.matrix[1][2]   ?
_struct_ncs_oper.matrix[1][3]   ?
_struct_ncs_oper.matrix[2][1]   ?
_struct_ncs_oper.matrix[2][2]   ?
_struct_ncs_oper.matrix[2][3]   ?
_struct_ncs_oper.matrix[3][1]   ?
_struct_ncs_oper.matrix[3][2]   ?
_struct_ncs_oper.matrix[3][3]   ?
_struct_ncs_oper.vector[1]      ?
_struct_ncs_oper.vector[2]      ?
_struct_ncs_oper.vector[3]      ?
#
_pdbx_post_process_details.entry_id      1A3I
_pdbx_post_process_details.text          ?
_pdbx_post_process_details.seq_details   ?
#
_pdbx_post_process_status.entry_id     1A3I
_pdbx_post_process_status.cycle_id     ?
_pdbx_post_process_status.date_begin   ?
_pdbx_post_process_status.date_end     ?
_pdbx_post_process_status.details      ?
_pdbx_post_process_status.annotator    ?
#
_refine_hist.pdbx_refine_id                   'X-RAY DIFFRACTION'
_refine_hist.cycle_id                         LAST
_refine_hist.pdbx_number_atoms_protein        0
_refine_hist.pdbx_number_atoms_nucleic_acid   0
_refine_hist.pdbx_number_atoms_ligand         5
_refine_hist.number_atoms_solvent             0
_refine_hist.number_atoms_total               5
_refine_hist.d_res_high                       .
_refine_hist.d_res_low                        .
#
_em_2d_projection_selection.entry_id        1A3I
_em_2d_projection_selection.software_name   ?
_em_2d_projection_selection.num_particles   ?
#
_em_3d_fitting.id                ?
_em_3d_fitting.entry_id          1A3I
_em_3d_fitting.ref_space         ?
_em_3d_fitting.ref_protocol      ?
_em_3d_fitting.target_criteria   ?
_em_3d_fitting.overall_b_value   ?
_em_3d_fitting.method            ?
#
_em_3d_reconstruction.entry_id                    1A3I
_em_3d_reconstruction.id                          ?
_em_3d_reconstruction.method                      ?
_em_3d_reconstruction.nominal_pixel_size          ?
_em_3d_reconstruction.actual_pixel_size           ?
_em_3d_reconstruction.resolution                  ?
_em_3d_reconstruction.ctf_correction_method       ?
_em_3d_reconstruction.magnification_calibration   ?
_em_3d_reconstruction.details                     ?
#
_em_assembly.id                  ?
_em_assembly.entry_id            1A3I
_em_assembly.aggregation_state   ?
_em_assembly.name                ?
_em_assembly.details             ?
#
_em_sample_preparation.entry_id               1A3I
_em_sample_preparation.id                     ?
_em_sample_preparation.buffer_id              ?
_em_sample_preparation.support_id             ?
_em_sample_preparation.sample_concentration   ?
_em_sample_preparation.ph                     ?
#
_em_vitrification.entry_id   1A3I
_em_vitrification.id         ?
_em_vitrification.details    ?
#
_em_experiment.entry_id                1A3I
_em_experiment.reconstruction_method   ?
_em_experiment.specimen_type           ?
#
_em_single_particle_entity.entry_id        1A3I
_em_single_particle_entity.symmetry_type   ?
#
_entry.id   1A3I
#
loop_
_atom_type.symbol
C
N
O
#
loop_
_pdbx_unobs_or_zero_occ_atoms.id
_pdbx_unobs_or_zero_occ_atoms.polymer_flag
_pdbx_unobs_or_zero_occ_atoms.occupancy_flag
_pdbx_unobs_or_zero_occ_atoms.PDB_model_num
_pdbx_unobs_or_zero_occ_atoms.auth_asym_id
_pdbx_unobs_or_zero_occ_atoms.auth_comp_id
_pdbx_unobs_or_zero_occ_atoms.auth_seq_id
_pdbx_unobs_or_zero_occ_atoms.PDB_ins_code
_pdbx_unobs_or_zero_occ_atoms.auth_atom_id
_pdbx_unobs_or_zero_occ_atoms.label_alt_id
1 N 1 1 B PRO 1 ? CG  ?
2 N 1 1 B PRO 1 ? CD  ?
3 N 1 1 B PRO 1 ? OXT ?
#
_pdbx_entity_nonpoly.entity_id   1
_pdbx_entity_nonpoly.name        PROLINE
_pdbx_entity_nonpoly.comp_id     PRO
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_label_seq_num
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x_esd
_atom_site.Cartn_y_esd
_atom_site.Cartn_z_esd
_atom_site.occupancy_esd
_atom_site.B_iso_or_equiv_esd
_atom_site.footnote_id
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_auth_seq_id
_atom_site.pdbx_auth_comp_id
_atom_site.pdbx_auth_asym_id
_atom_site.pdbx_auth_atom_name
_atom_site.pdbx_PDB_model_num
HETATM 1 N N  . PRO B 1 . 1 ? 8.316 21.206 21.530 1.00 17.44 ? ? ? ? ? ? ? 1 PRO B N  1 PRO B N  1
HETATM 2 C CA . PRO B 1 . 1 ? 7.608 20.729 20.336 1.00 17.44 ? ? ? ? ? ? ? 1 PRO B CA 1 PRO B CA 1
HETATM 3 C C  . PRO B 1 . 1 ? 8.487 20.707 19.092 1.00 17.44 ? ? ? ? ? ? ? 1 PRO B C  1 PRO B C  1
HETATM 4 O O  . PRO B 1 . 1 ? 9.466 21.457 19.005 1.00 17.44 ? ? ? ? ? ? ? 1 PRO B O  1 PRO B O  1
HETATM 5 C CB . PRO B 1 . 1 ? 6.460 21.723 20.211 1.00 22.26 ? ? ? ? ? ? ? 1 PRO B CB 1 PRO B CB 1
#"""
        f = NamedTemporaryFile(delete=False, suffix=".cif")
        f.write(CIFString)
        f.close()
        cls.CIFFileName = f.name
