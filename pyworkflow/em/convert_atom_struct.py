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
# cif to pdb convertion is based on  cif2pdb.py by Spencer Bliven
# <spencer.bliven@gmail.com>
#
# see http://biopython.org/DIST/docs/tutorial/Tutorial.html for a description
# on the object "structure" and other Bio.xxxx modules

import os
import numpy

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO, MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import Entity
from Bio.PDB import PDBList
from pyworkflow.em.transformations import rotation_from_matrix, \
    translation_from_matrix


class OutOfChainsError(Exception):
    pass

class AtomicStructHandler():
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
        if fileName is not None:
            self.read(fileName)

    def readFromPDBDatabase(self, pdbID, dir=None, type='mmCif'):
        """
        Retrieve structure from PDB
        :param pdbID:
        :param dir: save structure in this directory
        :param type:  mmCif or pdb
        :return:
        """
        if dir is None:
            dir = os.getcwd()
        pdbl = PDBList()
        fileName = pdbl.retrieve_pdb_file(pdbID, pdir=dir, file_format=type)
        self.read(fileName)
        return fileName

    def getStructure(self):
        return self.structure

    def read(self, fileName):
        """ Read and parse file."""
        # biopython asigns an ID to any read structure
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
                self.ioCIF = MMCIFIO()
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

def cifToPdb(fnCif,fnPdb):
    h = AtomicStructHandler()
    h.read(fnCif)
    h.writeAsPdb(fnPdb)
