# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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
# cif to pdb conversion is based on  cif2pdb.py by Spencer Bliven
# <spencer.bliven@gmail.com>
#
# see http://biopython.org/DIST/docs/tutorial/Tutorial.html for a description
# on the object "structure" and other Bio.xxxx modules

from __future__ import print_function
import os
import numpy
from collections import defaultdict

from Bio.PDB.Dice import ChainSelector
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO, MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import Entity
from Bio.PDB import PDBList
from collections import OrderedDict
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Polypeptide import three_to_one
from Bio.Seq import Seq
from .transformations import translation_from_matrix


class OutOfChainsError(Exception):
    pass


class scipionMMCIFIO(MMCIFIO):
    """ Class that redefines the name of the chains.
     The current biopython mmCIF parser fills label_asym_id
      with unique values and auth_asym_id with the chan ids
      of the input atomic structures. I think auth_asym_id
      should be equal to auth_asym_id if they are unique.
      Chimera uses label_asym_id for chain id while
      coot uses auth_asym_id. Furthermore, auth_asym_id field
      is a direct match to the PDB chain identifier (when it exists).
      I hope new vesion of biopython will not conflict with this change.

      Modifications of original code are marked with 'ROB'"""
    def _noRepeated(self, structure):
        chains = [c for c in structure.get_chains()]
        uniqueChains = set(chains)
        len(chains)
        len(uniqueChains)
        if len(chains) == len(uniqueChains):
            return True
        else:
            return False

    def _save_structure(self, out_file, select, preserve_atom_numbering):
        if self._noRepeated(self.structure):  # if the chain ids are unique
            atom_dict = defaultdict(list)

            for model in self.structure.get_list():
                if not select.accept_model(model):
                    continue
                # mmCIF files with a single model have it specified as model 1
                if model.serial_num == 0:
                    model_n = "1"
                else:
                    model_n = str(model.serial_num)
                # This is used to write label_entity_id and label_asym_id and
                # increments from 1, changing with each molecule
                entity_id = 0
                if not preserve_atom_numbering:
                    atom_number = 1
                for chain in model.get_list():
                    if not select.accept_chain(chain):
                        continue
                    chain_id = chain.get_id()
                    if chain_id == " ":
                        chain_id = "."
                    # This is used to write label_seq_id and increments from 1,
                    # remaining blank for hetero residues
                    #####residue_number = 1
                    # prev_residue_type = ""
                    # prev_resname = ""

                    for residue in chain.get_unpacked_list():
                        if not select.accept_residue(residue):
                            continue
                        hetfield, resseq, icode = residue.get_id()
                        label_seq_id = str(resseq)
                        if hetfield == " ":
                            residue_type = "ATOM"
                            ####residue_number += 1
                        else:
                            residue_type = "HETATM"
                            ###label_seq_id = "."
                        resseq = str(resseq)
                        if icode == " ":
                            icode = "?"
                        resname = residue.get_resname()
                        # Check if the molecule changes within the chain
                        # This will always increment for the first residue in a
                        # chain due to the starting values above
                        ## ROB if residue_type != prev_residue_type or \
                        ## ROB        (residue_type == "HETATM" and resname != prev_resname):
                        ## ROB    entity_id += 1
                        # prev_residue_type = residue_type
                        # prev_resname = resname
                        ## ROB label_asym_id = self._get_label_asym_id(entity_id)
                        for atom in residue.get_unpacked_list():
                            if select.accept_atom(atom):
                                atom_dict["_atom_site.group_PDB"].append(residue_type)
                                if preserve_atom_numbering:
                                    atom_number = atom.get_serial_number()
                                atom_dict["_atom_site.id"].append(str(atom_number))
                                if not preserve_atom_numbering:
                                    atom_number += 1
                                element = atom.element.strip()
                                if element == "":
                                    element = "?"
                                atom_dict["_atom_site.type_symbol"].append(element)
                                atom_dict["_atom_site.label_atom_id"].append(atom.get_name().strip())
                                altloc = atom.get_altloc()
                                if altloc == " ":
                                    altloc = "."
                                atom_dict["_atom_site.label_alt_id"].append(altloc)
                                atom_dict["_atom_site.label_comp_id"].append(resname.strip())
                                # modified by ROB BEGIN
                                atom_dict["_atom_site.label_asym_id"].append(chain_id)
                                # modified by ROB END
                                # The entity ID should be the same for similar chains
                                # However this is non-trivial to calculate so we write "?"
                                atom_dict["_atom_site.label_entity_id"].append("?")
                                atom_dict["_atom_site.label_seq_id"].append(label_seq_id)
                                atom_dict["_atom_site.pdbx_PDB_ins_code"].append(icode)
                                coord = atom.get_coord()
                                atom_dict["_atom_site.Cartn_x"].append("%.3f" % coord[0])
                                atom_dict["_atom_site.Cartn_y"].append("%.3f" % coord[1])
                                atom_dict["_atom_site.Cartn_z"].append("%.3f" % coord[2])
                                atom_dict["_atom_site.occupancy"].append(str(atom.get_occupancy()))
                                atom_dict["_atom_site.B_iso_or_equiv"].append(str(atom.get_bfactor()))
                                atom_dict["_atom_site.auth_seq_id"].append(resseq)
                                atom_dict["_atom_site.auth_asym_id"].append(chain_id)
                                atom_dict["_atom_site.pdbx_PDB_model_num"].append(model_n)

            # Data block name is the structure ID with special characters removed
            structure_id = self.structure.id
            for c in ["#", "$", "'", "\"", "[", "]", " ", "\t", "\n"]:
                structure_id = structure_id.replace(c, "")
            atom_dict["data_"] = structure_id

            # Set the dictionary and write out using the generic dictionary method
            self.dic = atom_dict
            self._save_dict(out_file)

        else:
            super(scipionMMCIFIO, self)._save_structure(out_file,
                                                        select,
                                                        preserve_atom_numbering)



class AtomicStructHandler:
    """ Class that contain utilities to handle pdb/cif files"""
    PDB = 0
    CIF = 1

    def __init__(self, fileName=None, permissive=1):
        # The PERMISSIVE flag indicates that a number of common problems
        # associated with PDB files will be ignored
        self.permissive = permissive
        self.pdbParser = None
        self.cifParser = None
        self.ioPDB = None
        self.ioCIF = None
        self.structure = None
        self._readDone = False
        if fileName is not None:
            self.read(fileName)

    def readFromPDBDatabase(self, pdbID, dir=None, type='mmCif'):
        """
        Retrieve structure from PDB
        :param pdbID:
        :param dir: save structure in this directory
        :param type:  mmCif or pdb
        :return: filename with pdb file
        """
        if dir is None:
            dir = os.getcwd()
        pdbl = PDBList()
        fileName = pdbl.retrieve_pdb_file(pdbID, pdir=dir, file_format=type)
        self.read(fileName)
        self.write(fileName)  # many  atom structs need to be fixed
                              # by reading and writing it we achive this
                              # goal
        return os.path.abspath(fileName)

    def getStructure(self):
        """return strcture information, model, chain,
           residues, atoms..."""
        return self.structure

    def read(self, fileName):
        """ Read and parse file."""
        # biopython assigns an ID to any read structure
        structure_id = os.path.basename(fileName)
        structure_id = structure_id[:4] if len(structure_id) > 4 else "1xxx"

        if fileName.endswith(".pdb") or fileName.endswith(".ent"):
            if self.pdbParser is None:
                self.pdbParser = PDBParser(PERMISSIVE=self.permissive)
            parser = self.pdbParser
            self.type = self.PDB
        else:
            if self.cifParser is None:
                self.cifParser = MMCIFParser()
            parser = self.cifParser
            self.type = self.CIF

        self.structure = parser.get_structure(structure_id, fileName)
        self._readDone = True

    def checkRead(self):
        if self._readDone:
            return True
        else:
            print("you must read the pdb file first")
            exit(0)

    def getModelsChains(self):
        """
        return a dic of all models and respective chains (chainID and
        length of residues) from a pdb file
        """
        self.checkRead()
        models = OrderedDict()

        for model in self.structure:
            chainDic = OrderedDict()
            for chain in model:
                if len(chain.get_unpacked_list()[0].resname) == 1: # RNA
                    seq = list()
                    for residue in chain:
                        if residue.get_resname() in ['A', 'C', 'G', 'U']:
                            seq.append(residue.get_resname())
                        else:
                            seq.append("X")
                elif len(chain.get_unpacked_list()[0].resname) == 2: # DNA
                    seq = list()
                    for residue in chain:
                        if residue.get_resname()[1] in ['A', 'C', 'G', 'T']:
                            seq.append(residue.get_resname()[1])
                        else:
                            seq.append("X")
                elif len(chain.get_unpacked_list()[0].resname) == 3: # Protein
                    seq = list()
                    counter = 0
                    for residue in chain:
                        if is_aa(residue.get_resname(), standard=True): # aminoacids
                            seq.append(three_to_one(residue.get_resname()))
                            counter += 1
                        else:
                            seq.append("X")
                    if counter == 0: # HETAM
                        for residue in chain:
                            seq.append(residue.get_resname())
                while seq[-1] == "X":
                    del seq[-1]
                while seq[0] == "X":
                    del seq[0]
                chainDic[chain.id] = len(seq)
            models[model.id] = chainDic

        return models

    def getSequenceFromChain(self, modelID, chainID):
        self.checkRead()
        seq = list()
        for model in self.structure:
            if model.id == modelID:
                for chain in model:
                    if str(chain.id) == chainID:
                        if len(chain.get_unpacked_list()[0].resname) == 1:
                            print("Your sequence is a nucleotide sequence (" \
                                  "RNA)\n")
                            # alphabet = IUPAC.IUPACAmbiguousRNA._upper()
                            for residue in chain:
                                ## Check if the residue belongs to the
                                ## standard RNA and add those residues to the
                                ## seq
                                if residue.get_resname() in ['A', 'C',
                                                                'G', 'U']:
                                    seq.append(residue.get_resname())
                                else:
                                    seq.append("X")
                        elif len(chain.get_unpacked_list()[0].resname) == 2:
                            print("Your sequence is a nucleotide sequence (" \
                                  "DNA)\n")
                            # alphabet = IUPAC.ExtendedIUPACDNA._upper()
                            for residue in chain:
                                ## Check if the residue belongs to the
                                ## standard DNA and add those residues to the
                                ## seq
                                if residue.get_resname()[1] in ['A', 'C',
                                                                'G', 'T']:
                                    seq.append(residue.get_resname()[1])
                                else:
                                    seq.append("X")
                        elif len(chain.get_unpacked_list()[0].resname) == 3:
                            counter = 0
                            for residue in chain:
                                if is_aa(residue.get_resname(), standard=True):
                                    # alphabet = IUPAC.ExtendedIUPACProtein._upper()
                                    ## The test checks if the amino acid
                                    ## is one of the 20 standard amino acids
                                    ## Some proteins have "UNK" or "XXX", or other symbols
                                    ## for missing or unknown residues
                                    seq.append(three_to_one(residue.get_resname()))
                                    counter += 1
                                else:
                                    seq.append("X")
                            if counter != 0:  # aminoacids
                                print("Your sequence is an aminoacid sequence")
                            else: # HETAM
                                print("Your sequence is a HETAM sequence")
                                for residue in chain:
                                    seq.append(residue.get_resname())
                        while seq[-1] == "X":
                            del seq[-1]
                        while seq[0] == "X":
                            del seq[0]
                        # return Seq(str(''.join(seq)), alphabet=alphabet)
                        return Seq(str(''.join(seq)))

    def getFullID(self, model_id='0', chain_id=None):
        """
        assign a label to a sequence obtained from a PDB file
        :parameter
        model_id
        chain_id
        :return: string with label
        """
        self.checkRead()  # cehck we have read the structure
        label = "%s" % self.structure.get_id()  # PDB_ID
        # check how many model are in sequence if > 1 add to id
        models = self.getModelsChains()
        if len(models) > 1:
            label = label + "__%s" % str(model_id)
        # if chain_id is None do not add it to te name
        if chain_id is not None:
            label = label + "_%s" % str(chain_id)
        return label


    def readLowLevel(self, fileName):
        """ Return a dictionary with all mmcif fields. you should parse them
            Example: get the list of the y coordinates of all atoms
              dict = readLowLevel("kk.pdb")
              y_list = dict['_atom_site.Cartn_y']
        """

        if fileName.endswith(".pdb"):
            print ("Low level access to PDB is not implemented")
        else:
            dict = MMCIF2Dict(fileName)
        return dict

    def _write(self, fileName):
        """ Do not use this function use toPDB or toCIF, they take care of some
        compatibiity issues"""
        if fileName.endswith(".pdb") or fileName.endswith(".ent"):
            if self.ioPDB is None:
                self.ioPDB = PDBIO()
            io = self.ioPDB
        else:
            if self.ioCIF is None:
                self.ioCIF = scipionMMCIFIO()
                # self.ioCIF = MMCIFIO()
            io = self.ioCIF
        io.set_structure(self.structure)
        io.save(fileName)

    def _writeLowLevel(self, fileName, dict):
        """ write a dictionary as cif file
        """

        if fileName.endswith(".pdb"):
            print ("Low level access to PDB is not implemented")
        else:
            if self.ioCIF is None:
                self.ioCIF = MMCIFIO()
            io = self.ioCIF
        io.set_dict(dict)
        io.save(fileName)

    def _intToChain(self, i, base=62):
        """
        int_to_chain(int,int) -> str
        Converts a positive integer to a chain ID. Chain IDs include uppercase
        characters, numbers, and optionally lowercase letters.
        i = a positive integer to convert
        base = the alphabet size to include. Typically 36 or 62.
        """
        if i < 0:
            raise ValueError("positive integers only")
        if base < 0 or 62 < base:
            raise ValueError("Invalid base")

        quot = int(i) // base
        rem = i % base
        if rem < 26:
            letter = chr(ord("A") + rem)
        elif rem < 36:
            letter = str(rem - 26)
        else:
            letter = chr(ord("a") + rem - 36)
        if quot == 0:
            return letter
        else:
            return self._intToChain(quot - 1, base) + letter

    def renameChains(self, structure):
        """Renames chains to be one-letter chains

        Existing one-letter chains will be kept.
        Multi-letter chains will be truncated
        or renamed to the next available letter of the alphabet.

        If more than 62 chains are present in the structure,
        raises an OutOfChainsError

        Returns a map between new and old chain IDs, as well as modifying
        the input structure

        """
        next_chain = 0  #
        # single-letters stay the same
        chainmap = {c.id: c.id
                    for c in structure.get_chains() if len(c.id) == 1}
        for o in structure.get_chains():
            if len(o.id) != 1:
                if o.id[0] not in chainmap:
                    chainmap[o.id[0]] = o.id
                    o.id = o.id[0]
                else:
                    c = self._intToChain(next_chain)
                    while c in chainmap:
                        next_chain += 1
                        c = self._intToChain(next_chain)
                        if next_chain >= 62:
                            raise OutOfChainsError()
                    chainmap[c] = o.id
                    o.id = c
        return chainmap

    def write(self, fileName):
        if fileName.endswith(".pdb") or fileName.endswith(".ent"):
            self.writeAsPdb(fileName)
        else:
            self.writeAsCif(fileName)

    def writeAsPdb(self, pdbFile):
        """ Save structure as PDB. Be aware that this is not a lossless convertion
        Returns False is conversion is not possible. True otherwise
        """
        # check input is not PDB
        if self.type == self.PDB:
            pass
        else:
            # rename long chains
            try:
                chainmap = self.renameChains(self.structure)
            except OutOfChainsError:
                print("Too many chains to represent in PDB format")
                return False

            for new, old in chainmap.items():
                if new != old:
                    print("Renaming chain {0} to {1}".format(old, new))

        self._write(pdbFile)

        return True

    def writeAsCif(self, cifFile):
        """ Save structure as CIF.
            Be aware that this is not a lossless convertion
        """
        self._write(cifFile)

    def centerOfMass(self, geometric=False):
        """
        Returns gravitic [default] or geometric center of mass of an Entity
        (anything with a get_atoms function in biopython.
        Geometric assumes all masses are equal (geometric=True)
        """
        entity = self.structure
        # Structure, Model, Chain, Residue
        if isinstance(entity, Entity.Entity):
            atom_list = entity.get_atoms()
        # List of Atoms
        elif hasattr(entity, '__iter__') and \
                [x for x in entity if x.level == 'A']:
            atom_list = entity
        else:  # Some other weirdo object
            raise ValueError("Center of Mass can only be calculated "
                             "from the following objects:\n"
                             "Structure, Model, Chain, Residue, "
                             "list of Atoms.")

        masses = []
        positions = [[], [], []]

        for atom in atom_list:
            masses.append(atom.mass)

            for i, coord in enumerate(atom.coord.tolist()):
                positions[i].append(coord)

        # If there is a single atom with undefined mass complain loudly.
        if 'ukn' in set(masses) and not geometric:
            raise ValueError("Some Atoms don't have an element assigned.\n"
                             "Try adding them manually or calculate the "
                             "geometrical center of mass instead.")

        if geometric:
            return [sum(coord_list) / len(masses) for coord_list in positions]
        else:
            w_pos = [[], [], []]
            for atom_index, atom_mass in enumerate(masses):
                w_pos[0].append(positions[0][atom_index] * atom_mass)
                w_pos[1].append(positions[1][atom_index] * atom_mass)
                w_pos[2].append(positions[2][atom_index] * atom_mass)

            return [sum(coord_list) / sum(masses) for coord_list in w_pos]

    def transform(self, transformation_matrix, sampling=1.):
        """ Geometrical transformation of a PDB structure

        :param entity: PDB biopython structure
        :param transformation_matrix -> 4x4 scipion matrix
        :paramsampling: scipion transform matrix is applied to voxels so
              length must be multiplied by samplingRate

        internal variables:
             rotation matrix -> numpy.array(
                      [[ 0.5, -0.809017,  0.309017],
                       [ 0.809017,  0.30917,  -0.5     ],
                       [ 0.309017,  0.5,       0.809017]])
             translation: translation vector -> numpy.array([1., 0., 0.], 'd')
        :return: no return, new data overwrites entity
        """
        # bioPhython and Scipion conventions do not match
        rotation_matrix = numpy.transpose(transformation_matrix[:3, :3])
        # from geometry get euler angles and recreate matrix
        translation = translation_from_matrix(transformation_matrix)
        translation = [x * sampling for x in translation]
        self.structure.transform(rotation_matrix, translation)

    def _renameChainsIfNeed(self, struct2):
        """Rename chain, we assume that there is a single model per structure"""
        repeated = False
        def RepresentsInt(s):
            try:
                int(s)
                return True
            except ValueError:
                return False
        import uuid
        chainIDs1 = [chain.id for chain in self.structure.get_chains()]
        for chain in struct2.get_chains():
            if chain.id in chainIDs1:
                repeated = True
                cId = chain.id
                l = len(cId)
                if l==1:
                    chain.id = "%s002" % cId
                elif RepresentsInt(cId[1:]): # try to fit a number and increase it by one
                    chain.id = "%s%03d" % (cId[0],int(cId[1:]) + 1)
                else: # generate a 4 byte random string
                    chain.id = uuid.uuid4().hex[:4]
        if repeated:
            self._renameChainsIfNeed(struct2)

    def addStruct(self, secondPDBfileName, outPDBfileName=None, useModel=False):
        """ Join the second structure to the first one.
            If cheon numes are the same rename them.
            if outPDBfileName id provided then new
            struct is saved to a file"""
        # read new structure
        if outPDBfileName is not None:
            pdbID = (os.path.splitext(os.path.basename(outPDBfileName))[0])[:4]
        else:
            pdbID = (os.path.splitext(os.path.basename(secondPDBfileName))[0])[:4]

        if secondPDBfileName.endswith(".pdb") or secondPDBfileName.endswith(".ent"):
            parser = PDBParser(PERMISSIVE=self.permissive)
        else:
            parser = MMCIFParser()

        struct2 = parser.get_structure(pdbID, secondPDBfileName)

        if useModel:
            modelNumber = 0
            modelID = 0
            # model.id = model.serial_num = len(self.structure)?  # not sure this
            # is valid always
            for model in self.structure:
                pass
            modelNumber = model.serial_num
            modelID = model.id
            for model in struct2:
                modelNumber +=1
                modelID +=1
                model.detach_parent()
                model.serial_num = modelNumber
                model.id = modelID
                self.structure.add(model)
        else:
            self._renameChainsIfNeed(struct2)

            for model in struct2:
                for chain in model:
                    chain.detach_parent()
                    self.structure[0].add(chain)


        # create new output file
        if outPDBfileName is not None:
            self.write(outPDBfileName)

    def extractChain(self, chainID, start=0, end=-1,
                     modelID='0', filename="output.mmcif"):
        """Code for chopping  a structure.
        This module is used internally by the Bio.PDB.extract() function.
        """
        sel = ChainSelector(chain_id=chainID, start=start,
                            end=end, model_id=int(modelID))
        io = scipionMMCIFIO()
        io.set_structure(self.structure)
        io.save(filename, sel)

def cifToPdb(fnCif,fnPdb):
    h = AtomicStructHandler()
    h.read(fnCif)
    h.writeAsPdb(fnPdb)

def pdbToCif(fnPdb, fnCif):
    h = AtomicStructHandler()
    h.read(fnPdb)
    h.writeAsCif(fnCif)

def toPdb(inFileName, outPDBFile):
    if inFileName.endswith(".pdb") or inFileName.endswith(".ent"):
        return inFileName
    elif inFileName.endswith(".cif") or inFileName.endswith(".mmcif"):
        cifToPdb(inFileName, outPDBFile)
        return outPDBFile
    else:
        print("ERROR (toPdb), Unknown file type for file = %s" % inFileName)

def toCIF(inFileName, outCIFFile):
    if inFileName.endswith(".cif") or inFileName.endswith(".mmcif"):
        return inFileName
    elif inFileName.endswith(".pdb") or inFileName.endswith(".ent"):
        pdbToCif(inFileName, outCIFFile)
        return outCIFFile
    else:
        print("ERROR (toCIF), Unknown file type for file = %s" % inFileName)

