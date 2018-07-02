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

from tkMessageBox import showerror
from protocol_molprobity import PhenixProtRunMolprobity
from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.em.viewer import TableView
import collections


def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except:
        print("Error:", msg)


class PhenixProtRunMolprobityViewer(ProtocolViewer):
    """ Viewer for Phenix program EMRinger
    """
    _label = 'MolProbity Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [PhenixProtRunMolprobity]

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        MOLPROBITYOUTFILENAME = self.protocol._getExtraPath(
            self.protocol.MOLPROBITYOUTFILENAME)
        self._parseFile(MOLPROBITYOUTFILENAME)

    def _defineParams(self, form):
        form.addSection(label='Visualization of MolProbity results')
        form.addParam('showFinalResults', LabelParam,
                      important=True,
                      label="Summary Table of Results",
                      help="Validation of protein geometry. Statistics "
                           "computed by the "
                           "PHENIX package using the same distributions as "
                           "the MolProbity web server."
                           "\n\nRamachandran outliers: Outlier "
                           "residues that show an unusual "
                           "combination of their phi and psi "
                           "torsion angles. Ramachandran outlier score is "
                           "computed as the percentage of Ramachandran "
                           "outliers regarding the total number of residues "
                           "in the entry showing outlier assessment."
                           "\n\nRamachandran "
                           "favored: Residues that show an normal "
                           "combination of their phi and psi "
                           "torsion angles. Ramachandran favored score is "
                           "computed as the percentage of Ramachandran "
                           "outliers regarding the total number of residues "
                           "in the entry showing outlier assessment. Between "
                           "the favored and outlier regions in the "
                           "Ramachandran plot, there is a small region of "
                           "allowed residues.\n\nRotamer outliers: Residues "
                           "that adopt unusual chi torsion angles "
                           "(non-rotameric), used to describe "
                           "the conformation of protein sidechains. The "
                           "score of rotamer outliers is computed as "
                           "the percentage of residues with non-rotameric "
                           "angles.\n\nC-beta outliers: Residues that show an "
                           "unusual deviation (higher than 0.25 A) of the "
                           "beta carbon from its ideal position. "
                           "This deviation of beta carbon "
                           "indicates incompatibilities between sidechain and "
                           "backbone.\n\nClashscore: Score associated to the "
                           "number of pairs of non-bonded atoms in the model "
                           "that are unusually close to each other. These "
                           "clashing atoms show unfavorable steric overlaps "
                           "of van der Waals shells. Steric clashes in "
                           "proteins are defined based on the Van der Waals "
                           "repulsion energy of the clashing atoms. "
                           "The clashscore is computed as the "
                           "number of serious clashes per 1000 atoms. "
                           "\n\nRMS (bonds): Root-mean-square deviation of "
                           "macromolecular bond lengths. "
                           "\n\nRMS (angles): Root-mean-square deviation of "
                           "macromolecular bond angles. "
                           "\n\nOverall score: Score that represents the "
                           "experimental resolution expected for a model of "
                           "this quality; ideally the score should be lower "
                           "than the actual resolution.\n")
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
            group.addParam('outliers', LabelParam,
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
            group.addParam('outliers', LabelParam,
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
            group.addParam('outliers', LabelParam,
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
            group.addParam('outliers', LabelParam,
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
        if self.outliers > 0:
            group.addParam('outliers', LabelParam,
                          label="Outliers", help="List of planar "
                                                 "group outliers ("
                                                 "sorted by deviation)")
        #
        # self.thresList = [str("%0.3f" % thres) for thres in self.dataDict[
        #     '_thresholds']]
        # _maxScoreIndex = self.dataDict['_maxScoreIndex']
        # form.addParam('threshold', EnumParam,
        #               choices=self.thresList,
        #               default=_maxScoreIndex,
        #               label="Density threshold",
        #               help="Choose one of the 20 map density cutoff values "
        #                    "used to compute the percentage of rotameric "
        #                    "residues. The Optimal threshold, at which the "
        #                    "maximum EMRinger score was obtained, "
        #                    "is shown by default.\n"
        #               )
        # form.addParam('showPeakCount', LabelParam,
        #               label="Peak count for the selected density threshold",
        #               help="Histograms for rotameric (blue) and non-rotameric "
        #                    "(red) residues at the selected map density.")
        # self.chainList = self.dataDict['_chains']
        # form.addParam('chain', EnumParam,
        #               choices=self.chainList,
        #               default=0,
        #               label="Chain",
        #               help="Choose one of the model chains. The first chain "
        #                    "is shown by default.")
        # form.addParam('showRollingWindows', LabelParam,
        #               label="Rolling window for the selected chain.",
        #               help="Visualize the rolling window EMRinger analysis of "
        #                    "your selected chain to distinguish regions of "
        #                    "improved model quality.This analysis was "
        #                    "performed on rolling sliding 21-residue windows "
        #                    "along the primary sequence of proteins.")
        # self.residueList = self.dataDict['_residues_format']
        # form.addParam('residue', EnumParam,
        #               choices=self.residueList,
        #               default=0,
        #               label="Residue",
        #               help="Choose one of the gamma-carbon-containing "
        #                    "residues (at least with one Chi angle) "
        #                    "located in one of the chains in the "
        #                    "position indicated. ")
        # form.addParam('showRingerResults', LabelParam,
        #               label="Ringer results for the selected residue.",
        #               help="Individual plots for each Chi angle of the "
        #                    "selected residue. Numeric values are showed in "
        #                    "the extra/*.csv file.")

    def _getVisualizeDict(self):
        return {
            'showFinalResults': self._visualizeFinalResults,
            'showBLrestraints': self._showBLrestraints,
            'showBLoutliers': self._showBLoutliers,
            'showBArestraints': self._showBArestraints,
            'showDArestraints': self._showDArestraints,
            'showCHILrestraints': self._showCHILrestraints,
            'showPLANARrestraints': self._showPLANARrestraints,
        }

    def _visualizeFinalResults(self, e=None):
        headerList = ['statistic', 'value']
        dictX = self.dictSummary
        val = 0.4
        mesg="Model Final Statistics"
        title="MolProbity: Final Results Summary"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showBLrestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictBLRestraints
        val = 0.3
        mesg="Bond Length Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showBLoutliers(selfself, e=None):
        pass

    def _showBArestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictBARestraints
        val = 0.3
        mesg="Bond Angle Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showDArestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictDARestraints
        val = 0.3
        mesg="Dihedral Angle Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showCHILrestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictChilRestraints
        val = 0.3
        mesg="Chirality Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showPLANARrestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictPlanarRestraints
        val = 0.3
        mesg="Chirality Restraints\n(Deviations from ideal values)"
        title="MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showResults(self, headerList, dictX, val, mesg, title):
        dataList = []
        for k, v in dictX.iteritems():
            if isinstance(v, int):
                dataList.append((k, v))
            elif isinstance(v, float):
                dataList.append((k, ("%" + str(val) + "f") % v))

        if not dataList:
            errorWindow(self.getTkRoot(), "No data available")
            return

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg=mesg,
                  title=title,
                  height=len(dataList), width=250, padding=40)

    def _parseFile(self, fileName):
        self.dictSummary = collections.OrderedDict()
        self.dictBLRestraints = collections.OrderedDict()
        self.dictBARestraints = collections.OrderedDict()
        self.dictDARestraints = collections.OrderedDict()
        self.dictChilRestraints = collections.OrderedDict()
        self.dictPlanarRestraints = collections.OrderedDict()
        with open(fileName) as f:
            line = f.readline()
            while line:
                words = line.strip().split()
                if len(words) > 1:
                    if (words[0] == 'Ramachandran' and words[1] == 'outliers'):
                        self.dictSummary['Ramachandran outliers'] = \
                            float(words[3])
                    elif (words[0] == 'favored' and words[1] == '='):
                        self.dictSummary['Ramachandran favored'] = \
                            float(words[2])
                    elif (words[0] == 'Rotamer' and words[1] == 'outliers'):
                        self.dictSummary['Rotamer outliers'] = float(words[3])
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
                                                  '> 4sigma'] = int(0)
                        elif (words[0] == 'atoms'):
                            cnt = 1
                            Atom1 = []
                            Atom2 = []
                            IdealValue = []
                            ModelValue = []
                            Deviation = []
                            line = f.readline()
                            words = line.strip().split()
                            Atom1.append(words[0]+ ' ' + words[1] + ' ' +
                                         words[2])
                            while (len(words) > 1 and words[1] == 'DMS'):
                                cnt += 1
                                line = f.readline()
                                words = line.strip().split()
                            self.dictBLRestraints['Number of outliers ' \
                                                  '> 4sigma'] = int(cnt / 2)
                            if words[0].startswith('C'):
                                Atom1.append(words[0] + ' ' + words[1] + ' ' +
                                             words[2])
                            if words[0].startswith('D'):
                                Atom2.append(words[0] + ' ' + words[1] + ' ' +
                                         words[2])
                                IdealValue.append(float(words[3]))
                                ModelValue.append(float(words[4]))
                                Deviation.append(float(words[8].split('*')[0]))
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
                                'Number of outliers > 4sigma'] = \
                                int(0)
                        elif (words[0] == 'atoms'):
                            cnt = 1
                            line = f.readline()
                            words = line.strip().split()
                            while (len(words) > 1 and len(words[2]) == 3):
                                cnt += 1
                                line = f.readline()
                                words = line.strip().split()
                            self.dictBARestraints[
                                'Number of outliers > 4sigma'] = \
                                int(cnt / 3)
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
                                'Number of outliers > 4sigma'] = \
                                int(0)
                        elif (words[0] == 'atoms'):
                            cnt = 1
                            line = f.readline()
                            words = line.strip().split()
                            while (len(words) > 1 and len(words[2]) == 3):
                                cnt += 1
                                line = f.readline()
                                words = line.strip().split()
                            self.dictDARestraints[
                                'Number of outliers > 4sigma'] = \
                                int(cnt / 4)
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
                                'Number of outliers > 4sigma'] = \
                                int(0)
                        elif (words[0] == 'atoms'):
                            # An example is necessary to test this part
                            cnt = 1
                            line = f.readline()
                            words = line.strip().split()
                            while (len(words) > 1 and len(words[2]) == 3):
                                cnt += 1
                                line = f.readline()
                                words = line.strip().split()
                            self.dictChilRestraints[
                                'Number of outliers > 4sigma'] = \
                                int(cnt / 4)
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
                                'Number of outliers > 4sigma'] = \
                                int(0)
                        elif (words[0] == 'atoms'):
                            # An example is necessary to test this part
                            cnt = 1
                            line = f.readline()
                            words = line.strip().split()
                            while (len(words) > 6 ):
                                cnt += 1
                                line = f.readline()
                                words = line.strip().split()
                            self.dictPlanarRestraints[
                                'Number of outliers > 4sigma'] = \
                                int(cnt / 4)
                line = f.readline()

