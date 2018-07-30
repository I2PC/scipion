# **************************************************************************
# *
# * Authors:     Roberto Marabini
# *              Marta Martinez
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

from protocol_molprobity import PhenixProtRunMolprobity
from pyworkflow.protocol.params import LabelParam
import collections
from pyworkflow.em.packages.ccp4.convert import (getProgram, runCCP4Program)
import os

from viewer_refinement_base import PhenixProtRefinementBaseViewer

class PhenixProtRunMolprobityViewer(PhenixProtRefinementBaseViewer):
    """ Viewer for Phenix program Molprobity
    """
    _label = 'MolProbity viewer'
    _targets = [PhenixProtRunMolprobity]

    COOT = 'coot'

    def __init__(self,  **kwargs):
         PhenixProtRefinementBaseViewer.__init__(self, **kwargs)
         MOLPROBITYOUTFILENAME = self.protocol._getExtraPath(
             self.protocol.MOLPROBITYOUTFILENAME)
         self._parseFile(MOLPROBITYOUTFILENAME)

    def _defineParams(self, form):
        PhenixProtRefinementBaseViewer._defineParams(self,form)

        form.addParam('showCootOutliers', LabelParam,
                      important=True,
                      label="Coot Visualization:\nRamachandran outliers\n "
                            "Rotamer outliers\nC-beta outliers\n"
                            "Severe clashes",
                      help="Interactive visualization of outliers and clashes"
                           " with Coot.")
        group = form.addGroup('Basic Geometry: Bond Length Restraints')
        group.addParam('showBLrestraints', LabelParam,
                      label="Deviations",
                      help="Check here the number of outlier atoms "
                           "according "
                           "to the bond length restraints between pairs of "
                           "linked"
                           "atoms.\nWarning!!!: Refined structures should not "
                           "have any outliers except those are obvious in high "
                           "resolution electron density maps.")
        self.outliers = self.dictBLRestraints['Number of outliers > 4sigma']
        if self.outliers > 0:
            group.addParam('showBLoutliers', LabelParam,
                      label="Outliers", help="List of outlier atoms ("
                                             "sorted by deviation) according "
                                             "to the bond length restraints")
        group = form.addGroup('Basic Geometry: Bond Angle Restraints')
        group.addParam('showBArestraints', LabelParam,
                      label="Deviations",
                      help="Check here the number of outlier residues "
                           "according "
                           "to the tau (N-alphaC-C) bond angle restraints.\n"
                           "Warning!!!: Refined structures should not "
                           "have any outliers except those are obvious in high "
                           "resolution electron density maps.")
        self.outliers = self.dictBARestraints['Number of outliers > 4sigma']
        if self.outliers > 0:
            group.addParam('showBAoutliers', LabelParam,
                          label="Outliers", help="List of outlier residues ("
                                                 "sorted by deviation) "
                                                 "according to the bond angle "
                                                 "restraints")
        group = form.addGroup('Basic Geometry: Dihedral Angle Restraints')
        group.addParam('showDArestraints', LabelParam,
                      label="Deviations",
                      help="Check here the number of outlier residues "
                           "according "
                           "to the side chain dihedral torsion (chi) angle "
                           "restraints.\n"
                           "Warning!!!: Refined structures should not "
                           "have any outliers except those are obvious in high "
                           "resolution electron density maps.")
        self.outliers = self.dictDARestraints['Number of outliers > 4sigma']
        if self.outliers > 0:
            group.addParam('showDAoutliers', LabelParam,
                          label="Outliers", help="List of outlier residues ("
                                                 "sorted by deviation) "
                                                 "according to the dihedral "
                                                 "angle restraints")
        group = form.addGroup('Basic Geometry: Chirality Restraints')
        group.addParam('showCHILrestraints', LabelParam,
                      label="Deviations",
                      help="Check here the number of outlier residues "
                           "according to the volume chirality "
                           "restraints.\n"
                           "Warning!!!: Refined structures should not "
                           "have any outliers except those are obvious in high "
                           "resolution electron density maps.")
        self.outliers = self.dictChilRestraints['Number of outliers > 4sigma']
        if self.outliers > 0:
            group.addParam('showCHILoutliers', LabelParam,
                          label="Outliers", help="List of outlier residues ("
                                                 "sorted by deviation) "
                                                 "according to the volume "
                                                 "chirality restraints")
        group = form.addGroup('Basic Geometry: Planarity Restraints')
        group.addParam('showPLANARrestraints', LabelParam,
                      label="Deviations",
                      help="Check here the number of outliers of planar "
                           "groups, such as aromatic rings,"
                           "according to the planar "
                           "restraints.\n"
                           "Warning!!!: Refined structures should not "
                           "have any outliers except those are obvious in high "
                           "resolution electron density maps.")
        self.outliers = self.dictPlanarRestraints['Number of outliers > 4sigma']
        if self.outliers > 0:
            group.addParam('showPLANARoutliers', LabelParam,
                           label="Outliers",
                           help='List of planar group outliers '
                                '(sorted by deviation)')

    def _getVisualizeDict(self):
        myDict = super(PhenixProtRefinementBaseViewer, self)._getVisualizeDict()

        myDict.update( {
            'showFinalResults': self._visualizeFinalResults,
            'showCootOutliers': self._showCootOutliers,
            'showBLrestraints': self._showBLrestraints,
            'showBLoutliers': self._showBLoutliers,
            'showBArestraints': self._showBArestraints,
            'showBAoutliers': self._showBAoutliers,
            'showDArestraints': self._showDArestraints,
            'showDAoutliers': self._showDAoutliers,
            'showCHILrestraints': self._showCHILrestraints,
            'showCHILoutliers': self._showCHILoutliers,
            'showPLANARrestraints': self._showPLANARrestraints,
            'showPLANARoutliers': self._showPLANARoutliers
            })
        return myDict

    def _showCootOutliers(self, e=None):
        MOLPROBITYCOOTFILENAME = self.protocol._getExtraPath(
            self.protocol.MOLPROBITYCOOTFILENAME)
        args = ""
        args += " --python " + MOLPROBITYCOOTFILENAME
        pdb = os.path.abspath(self.protocol.inputStructure.get().getFileName())
        args += " " + "--pdb " + pdb
        vol = self.protocol._getInputVolume()
        if vol is not None:
            MOLPROBITYFILE = os.path.abspath(self.protocol._getExtraPath(
                self.protocol.MOLPROBITYFILE))
            args += " " + "--map " + MOLPROBITYFILE

        # run coot with args
        runCCP4Program(getProgram(self.COOT), args)

    def _showBLrestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictBLRestraints
        val = 0.3
        mesg="Bond Length Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showBLoutliers(self, e=None):
        headerList = ['Atom1', 'Atom2', 'Ideal value', 'Model value',
                      'Deviation (sigmas)']
        dataList = self.blDataList
        mesg = "List of outliers (sorted by deviation)"
        title = "Bond Length Restraints"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showBArestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictBARestraints
        val = 0.3
        mesg="Bond Angle Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showBAoutliers(self, e=None):
        headerList = ['Atoms', 'Ideal value', 'Model value', 'Deviation ('
                                                            'sigmas)']
        dataList = self.baDataList
        mesg = "List of outliers (sorted by deviation)"
        title = "Bond Angle Restraints"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showDArestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictDARestraints
        val = 0.3
        mesg="Dihedral Angle Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showDAoutliers(self, e=None):
        headerList = ['Atoms', 'Ideal value', 'Model value', 'Deviation ('
                                                            'sigmas)']
        dataList = self.daDataList
        mesg = "List of outliers (sorted by deviation)"
        title = "Dihedral Angle Restraints"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showCHILrestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictChilRestraints
        val = 0.3
        mesg="Chirality Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showCHILoutliers(self, e=None):
        headerList = ['Atoms', 'Ideal value', 'Model value', 'Deviation ('
                                                            'sigmas)']
        dataList = self.chilDataList
        mesg = "List of outliers (sorted by deviation)"
        title = "Chirality Restraints"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showPLANARrestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictPlanarRestraints
        val = 0.3
        mesg="Planarity Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showPLANARoutliers(self, e=None):
        headerList = ['Atoms', 'Max. delta', 'RMS (delta)', 'Deviation ('
                                                            'sigmas)']
        dataList = self.planarDataList
        mesg = "List of outliers (sorted by deviation)"
        title = "Planarity Restraints"
        self._showOutliers(headerList, dataList, mesg, title)


    def _parseFile(self, fileName):
        self.dictSummary = collections.OrderedDict()
        self.dictBLRestraints = collections.OrderedDict()
        self.blDataList = []
        self.dictBARestraints = collections.OrderedDict()
        self.baDataList = []
        self.dictDARestraints = collections.OrderedDict()
        self.daDataList = []
        self.dictChilRestraints = collections.OrderedDict()
        self.chilDataList = []
        self.dictPlanarRestraints = collections.OrderedDict()
        self.planarDataList = []
        with open(fileName) as f:
            line = f.readline()
            while line:
                words = line.strip().split()
                if len(words) > 1:
                    if (words[0] == 'Ramachandran' and words[1] == 'outliers'):
                        self.dictSummary['Ramachandran outliers (%)'] = \
                            float(words[3])
                    elif (words[0] == 'favored' and words[1] == '='):
                        self.dictSummary['Ramachandran favored (%)'] = \
                            float(words[2])
                    elif (words[0] == 'Rotamer' and words[1] == 'outliers'):
                        self.dictSummary['Rotamer outliers (%)'] = float(
                            words[3])
                    elif (words[0] == 'C-beta' and words[1] == 'deviations'):
                        self.dictSummary['C-beta outliers'] = int(words[3])
                    elif (words[0] == 'Clashscore' and words[1] == '='):
                        self.dictSummary['Clashscore'] = float(words[2])
                    elif (words[0] == 'RMS(bonds)' and words[1] == '='):
                        self.dictSummary['RMS (bonds)'] = float(words[2])
                    elif (words[0] == 'RMS(angles)' and words[1] == '='):
                        self.dictSummary['RMS (angles)'] = float(words[2])
                    elif (words[0] == 'MolProbity' and words[1] == 'score'):
                        self.dictSummary['Overall score'] = float(words[3])
                    elif (words[0] == 'bond' and words[1] == ':'):
                        self.dictBLRestraints['Number of restraints'] = int(
                            words[4])
                        self.dictBLRestraints['RMS (deviation)'] = float(
                            words[2])
                        self.dictBLRestraints['Max deviation'] = float(
                            words[3])
                    elif (words[0] == '----------Bond' and words[1] ==
                        'lengths----------'):
                        f.readline()
                        line = f.readline()
                        words = line.strip().split()
                        if (words[0] == 'All' and words[1] == 'restrained'):
                            self.dictBLRestraints['Number of outliers ' \
                                                  '> 4sigma'] = 0
                        elif (words[0] == 'atoms'):
                            line = f.readline()
                            words = line.strip().split()
                            self._parseFileAtom1Atom2(words, f)
                            self.dictBLRestraints['Number of outliers '\
                                                  '> 4sigma'] = len(self.Atom1)
                            self.blDataList = zip(self.Atom1, self.Atom2,
                                                  self.IdealValue,
                                                  self.ModelValue,
                                                  self.Deviation)

                    elif (words[0] == 'angle' and words[1] == ':'):
                        self.dictBARestraints['Number of restraints'] = int(
                            words[4])
                        self.dictBARestraints['RMS (deviation)'] = float(
                            words[2])
                        self.dictBARestraints['Max deviation'] = float(
                            words[3])
                    elif (words[0] == '----------Bond' and words[1] ==
                        'angles----------'):
                        f.readline()
                        line = f.readline()
                        words = line.strip().split()
                        if (words[0] == 'All' and words[1] == 'restrained'):
                            self.dictBARestraints[
                                'Number of outliers > 4sigma'] = 0
                        elif (words[0] == 'atoms'):
                            line = f.readline()
                            words = line.strip().split()
                            self._parseFileAtom123(words, f)
                            self.dictBARestraints[
                                'Number of outliers > 4sigma'] = len(self.Atom1)
                            for a1, a2, a3, iv, mv, d in zip (self.Atom1,
                                                              self.Atom2,
                                                              self.Atom3,
                                                              self.IdealValue,
                                                              self.ModelValue,
                                                              self.Deviation):
                                self.baDataList.append((a1 + ", " + a2 + ", "
                                                                         "" +
                                                        a3, iv, mv, d))
                    elif (words[0] == 'dihedral' and words[1] == ':'):
                        self.dictDARestraints['Number of restraints'] = int(
                            words[4])
                        self.dictDARestraints['RMS (deviation)'] = float(
                            words[2])
                        self.dictDARestraints['Max deviation'] = float(
                            words[3])
                    elif (words[0] == '----------Dihedral' and words[1] ==
                        'angles----------'):
                        f.readline()
                        line = f.readline()
                        words = line.strip().split()
                        if (words[0] == 'All' and words[1] == 'restrained'):
                            self.dictDARestraints[
                                'Number of outliers > 4sigma'] = 0
                        elif (words[0] == 'atoms'):
                            line = f.readline()
                            words = line.strip().split()
                            self._parseFileAtom1234(words, f)
                            self.dictDARestraints[
                                'Number of outliers > 4sigma'] = len(self.Atom1)

                            for a1, a2, a3, a4, iv, mv, d in zip (self.Atom1,
                                                                  self.Atom2,
                                                                  self.Atom3,
                                                                  self.Atom4,
                                                                  self.IdealValue,
                                                                  self.ModelValue,
                                                                  self.Deviation):
                                element = a1 + ", " + a2 + ", " + a3 + ", " + a4
                                if len(element) > 46:
                                    element = element.replace(element[46:],
                                                              "...")
                                self.daDataList.append((element,
                                                       iv, mv, d))
                    elif (words[0] == 'chirality' and words[1] == ':'):
                        self.dictChilRestraints['Number of restraints'] = int(
                            words[4])
                        self.dictChilRestraints['RMS (deviation)'] = float(
                            words[2])
                        self.dictChilRestraints['Max deviation'] = float(
                            words[3])
                    elif (words[0] == '----------Chiral' and words[1] ==
                        'volumes----------'):
                        f.readline()
                        line = f.readline()
                        words = line.strip().split()
                        if (words[0] == 'All' and words[1] == 'restrained'):
                            self.dictChilRestraints[
                                'Number of outliers > 4sigma'] = 0
                        elif (words[0] == 'atoms'):
                            line = f.readline()
                            words = line.strip().split()
                            self._parseFileAtom1234(words, f)
                            self.dictChilRestraints[
                                'Number of outliers > 4sigma'] = len(self.Atom1)

                            for a1, a2, a3, a4, iv, mv, d in zip(self.Atom1,
                                                                 self.Atom2,
                                                                 self.Atom3,
                                                                 self.Atom4,
                                                                 self.IdealValue,
                                                                 self.ModelValue,
                                                                 self.Deviation):
                                element = a1 + ", " + a2 + ", " + a3 + ", " + a4
                                if len(element) > 46:
                                    element = element.replace(element[46:],
                                                              "...")
                                self.chilDataList.append((element, iv, mv, d))
                    elif (words[0] == 'planarity' and words[1] == ':'):
                        self.dictPlanarRestraints['Number of restraints'] = \
                            int(words[4])
                        self.dictPlanarRestraints['RMS (deviation)'] = float(
                            words[2])
                        self.dictPlanarRestraints['Max deviation'] = float(
                            words[3])
                    elif (words[0] == '----------Planar' and words[1] ==
                        'groups----------'):
                        f.readline()
                        line = f.readline()
                        words = line.strip().split()
                        if (words[0] == 'All' and words[1] == 'restrained'):
                            self.dictPlanarRestraints[
                                'Number of outliers > 4sigma'] = 0
                        elif (words[0] == 'atoms'):
                            line = f.readline()
                            words = line.strip().split()
                            self._parseFileGroups(words, f)
                            self.dictPlanarRestraints[
                                'Number of outliers > 4sigma'] = int(self.cnt)

                            for a, md, rd, d in zip(self.groups,
                                                    self.MaxDelta,
                                                    self.RMSDelta,
                                                    self.Deviation):
                                element = str()
                                for i in range(len(a) -1):
                                    element += a[i] + ", "
                                element += a[len(a) - 1]
                                if len(element) > 35:
                                    element = element.replace(element[35:],
                                                             "...")
                                self.planarDataList.append(
                                    (element, md, rd, d))
                line = f.readline()

    def _parseFileAtom1Atom2(self, words, f):
        self.Atom1 = []
        self.Atom2 = []
        self.IdealValue = []
        self.ModelValue = []
        self.Deviation = []
        while (len(words) > 1):
            if len(words) == 3:
                self.Atom1.append(words[0] + ' ' +
                             words[1] + ' ' + words[2])
            elif len(words) == 4:
                self.Atom1.append(words[0] + ' ' + words[1] + ' ' +
                             words[2] + ' ' + words[3])
            elif len(words) == 10:
                self.Atom2.append(words[0] + ' ' + words[1] + ' ' +
                             words[2] + ' ' + words[3])
                self.IdealValue.append((words[4]))
                self.ModelValue.append((words[5]))
                self.Deviation.append(
                    (words[9].split('*')[0]))
            elif len(words) == 9:
                self.Atom2.append(words[0] + ' ' +
                             words[1] + ' ' + words[2])
                self.IdealValue.append((words[3]))
                self.ModelValue.append((words[4]))
                self.Deviation.append(
                    (words[8].split('*')[0]))
            line = f.readline()
            words = line.strip().split()

    def _parseFileAtom123(self, words, f):
        self.Atom1 = []
        self.Atom2 = []
        self.Atom3 = []
        Atom123 = [self.Atom1, self.Atom2, self.Atom3]
        self.IdealValue = []
        self.ModelValue = []
        self.Deviation = []
        while len(words) > 1:
            for atom in Atom123:
                if len(words) == 4:
                    atom.append(words[0] + ' ' + words[1] +
                                ' ' + words[2] + ' ' + words[3])
                if (len(words) == 10):
                    atom.append(words[0] + ' ' + words[1] +
                                ' ' + words[2] + ' ' + words[3])
                    self.IdealValue.append(float(words[4]))
                    self.ModelValue.append(float(words[5]))
                    self.Deviation.append(
                        float(words[9].split('*')[0]))
                if len(words) == 3:
                    atom.append(words[0] + ' ' + words[1] +
                                ' ' + words[2])
                if (len(words) == 9):
                    atom.append(words[0] + ' ' + words[1] +
                                ' ' + words[2])

                    self.IdealValue.append(float(words[3]))
                    self.ModelValue.append(float(words[4]))
                    self.Deviation.append(
                        float(words[8].split('*')[0]))
                line = f.readline()
                words = line.strip().split()

    def _parseFileAtom1234(self, words, f):
        self.Atom1 = []
        self.Atom2 = []
        self.Atom3 = []
        self.Atom4 = []
        Atom1234 = [self.Atom1, self.Atom2, self.Atom3, self.Atom4]
        self.IdealValue = []
        self.ModelValue = []
        self.Deviation = []

        while (len(words) > 1):
            for atom in Atom1234:
                if len(words) == 4:
                    atom.append(words[0] + ' ' + words[1] +
                                ' ' + words[2] + ' ' + words[3])
                if (len(words) == 10):
                    atom.append(words[0] + ' ' + words[1] +
                                ' ' + words[2] + ' ' + words[3])
                    self.IdealValue.append(float(words[4]))
                    self.ModelValue.append(float(words[5]))
                    self.Deviation.append(
                        float(words[9].split('*')[0]))
                if len(words) == 3:
                    atom.append(words[0] + ' ' + words[1] +
                                ' ' + words[2])
                if (len(words) == 9):
                    atom.append(words[0] + ' ' + words[1] +
                                ' ' + words[2])
                    self.IdealValue.append(float(words[3]))
                    self.ModelValue.append(float(words[4]))
                    self.Deviation.append(
                        float(words[8].split('*')[0]))
                line = f.readline()
                words = line.strip().split()

    def _parseFileGroups(self, words, f):
        self.cnt = 0
        self.groups = []
        Atoms = []
        self.MaxDelta = []
        self.RMSDelta = []
        self.Deviation = []
        while (len(words) > 1):
            if len(words) == 4:
                Atoms.append(words[0] + " " + words[1] +
                             " " + words[2] + " " +
                             words[3])
            if len(words) == 8:
                self.cnt += 1
                Atoms.append(words[0] + " " + words[1] +
                             " " + words[2] + " " +
                             words[3])
                aminoacid = Atoms
                Atoms = []
                self.groups.append(aminoacid)
                self.MaxDelta.append(float(words[5]))
                self.RMSDelta.append(float(words[4]))
                self.Deviation.append(
                    float(words[7].split('*')[0]))
            if len(words) == 3:
                Atoms.append(words[0] + " " + words[1] +
                             " " + words[2])
            if len(words) == 7:
                self.cnt += 1
                Atoms.append(words[0] + " " + words[1] +
                             " " + words[2])
                aminoacid = Atoms
                Atoms = []
                self.groups.append(aminoacid)
                self.MaxDelta.append(float(words[4]))
                self.RMSDelta.append(float(words[3]))
                self.Deviation.append(
                    float(words[6].split('*')[0]))
            line = f.readline()
            words = line.strip().split()
