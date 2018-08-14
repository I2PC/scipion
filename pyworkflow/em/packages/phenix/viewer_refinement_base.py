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

from tkMessageBox import showerror
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.em.viewer import TableView
from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.em.packages.ccp4.convert import getProgram, runCCP4Program
from pyworkflow.em.viewers.chimera_utils import \
    createCoordinateAxisFile, runChimeraProgram, getProgram
import collections
import json
import os
import glob
from convert import runPhenixProgram
import matplotlib.pyplot as plt




def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except:
        print("Error:", msg)


class PhenixProtRefinementBaseViewer(ProtocolViewer):
    """ base for molprovity and real space refine programs
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    COOT = 'coot'
    ANALYSISTMPFILE = 'tmpAnalysisFile.txt'
    RAMATMPFILE = 'Rhamachandran_write_plot.py'
    ROTATMPFILE = "Rotamer_write_plot.py"
    MULTICPLOTTMPFILE = "Multi_criterion_plot.py"

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        MOLPROBITYOUTFILENAME = self.protocol._getExtraPath(
            self.protocol.MOLPROBITYOUTFILENAME)
        self._parseFile(MOLPROBITYOUTFILENAME)
        self.MOLPROBITYPKLFILENAME = self.protocol._getExtraPath(
            self.protocol.MOLPROBITYPKLFILENAME)
        self._writePickleData()
        self.dictOverall = json.loads(self.dictOverall,
                                   object_pairs_hook=collections.OrderedDict)

    def _defineParams(self, form):
        form.addSection(label="Volume and models")
        form.addParam('displayMapModel', LabelParam,
                      label="Volume and models",
                      help="Display of input volume, input pdb that has to be"
                           "refined and final refined model of the structure.")
        form.addSection(label='MolProbity results')
        group = form.addGroup('Summary MolProbity')
        group.addParam('showMolProbityResults', LabelParam,
                      important=True,
                      label="MolProbity Basic Statistics",
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
        group.addParam('showCootOutliers', LabelParam,
                      important=True,
                      label="Open in Coot",
                      help="Interactive visualization of outliers and clashes"
                           " with Coot:\n\nRamachandran outliers\n"
                            "Rotamer outliers\nC-beta outliers\n"
                            "Severe clashes ")
        if self.dictOverall['_len_missing_atoms'] > 0:
            group.addParam('showMissingAtoms', LabelParam,
                           label="Missing atoms",
                           help="For clarity, hydrogen atoms are not included")
        group = form.addGroup('Basic Geometry: Bond Length Restraints')
        group.addParam('showBLrestraints', LabelParam,
                       label="Deviations",
                       help="Check here the number of outlier atoms "
                            "according to the bond length restraints "
                            "between pairs of linked atoms.\nWarning!!!: "
                            "Refined structures should not have any outliers"
                            " except those are obvious in high "
                            "resolution electron density maps.\n")
        self.outliers = self.dictBLRestraints['Number of outliers > 4sigma']
        if self.outliers > 0:
            group.addParam('showBLoutliers', LabelParam, important=True,
                           label="Outliers",
                           help="List of outlier atoms (sorted by deviation) "
                                "according to the bond length restraints.\n")
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
            group.addParam('showBAoutliers', LabelParam, important=True,
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
            group.addParam('showDAoutliers', LabelParam, important=True,
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
            group.addParam('showCHILoutliers', LabelParam, important=True,
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
            group.addParam('showPLANARoutliers', LabelParam, important=True,
                           label="Outliers",
                           help='List of planar group outliers '
                                '(sorted by deviation)')
        if (self.dictOverall['_protein'] == True):
            group = form.addGroup('Protein')
            self.plotList = [u'Ramachandran plot', u'Chi1-Chi2 plot']
            group.addParam('plotType', EnumParam,
                           choices=self.plotList, default=0,
                           label="Select plot:",
                           help="Select a plot type. Ramachandran plot is "
                                "chosen by default. Chi1-Chi2 plot shows "
                                "rotameric angles\n")
            group.addParam('showPlotType', LabelParam, important=True,
                           label="View plot",
                           help="")
            if (self.dictOverall['_percent_rama_outliers'] > 0.):
                group.addParam('showRamaOutliersTable', LabelParam,
                                important=True,
                                label="Ramachandran outliers:",
                                help="")
            else:
                group.addParam('showMesgNoRamaOutliers', LabelParam,
                                label="No Ramachandran outliers detected",
                                help="")
            if (self.dictOverall['_percent_rota_outliers'] > 0.):
                group.addParam('showRotaOutliersTable', LabelParam,
                                important=True,
                                label="Rotamer outliers: ",
                                help="Although a residue may lie in the "
                                     "favored regions of the Chi1-Chi2 plot, "
                                     "outliers are flagged based on the "
                                     "distribution of all non-branched Chi "
                                     "angles in a residue.\nZero outliers is "
                                     "not the goal. Rotamer outliers can be  "
                                     "justified by sufficiently strong " \
                                     "electron density, van der Waals "
                                     "packing, and/or hydrogen bonds.\n")
            else:
                group.addParam('showMesgNoRotaOutliers', LabelParam,
                                label="No Rotamer outliers detected",
                                help="Although a residue may lie in the "
                                     "favored regions of the Chi1-Chi2 plot,"
                                     " outliers are flagged based on the "
                                     "distribution of all non-branched Chi "
                                     "angles in a residue.\nZero outliers "
                                     "is not the goal. Rotamer outliers can"
                                     " be  justified by sufficiently strong " \
                                     "electron density, van der Waals "
                                     "packing, and/or hydrogen bonds.\n")
            if (self.dictOverall['_n_cbeta_outliers'] > 0):
                group.addParam('showCbetaOutliersTable', LabelParam,
                                important=True,
                                label="C-beta outliers:",
                                help="C-beta position outliers (position "
                                     "deviates from ideal by more than "
                                     "0.25A).\n\nIdeal CB position is "
                                     "determined from the average of the "
                                     "ideal C-N-CA-CB and N-C-CA-CB dihedrals."
                                     " This measure is more sensitive than"
                                     " individual measures to both sidechain "
                                     "and mainchain misfittings.\n")
            else:
                group.addParam('showMesgNoCbetaOutliers', LabelParam,
                                label="No C-beta position outliers detected",
                                help="C-beta position outliers (position "
                                     "deviates from ideal by more than "
                                     "0.25A).\n\nIdeal CB position is "
                                     "determined from the average of the "
                                     "ideal C-N-CA-CB and N-C-CA-CB dihedrals."
                                     " This measure is more sensitive than"
                                     " individual measures to both sidechain "
                                     "and mainchain misfittings.\n")
            if (self.dictOverall['_n_nqh_flips_outliers'] > 0):
                group.addParam('showBackAsnGlnHisSidechains', LabelParam,
                                important=True,
                                label="Recommended Asn/Gln/His sidechain "
                                      "flips:",
                                help="Asn, Gln, and His sidechains are "
                                     "asymmetric and may require flipping to "
                                     "form favorable van der Waals contacts "
                                     "and hydrogen bonding.\n\nREDUCE "
                                     "(phenix.reduce) has been run on your "
                                     "file to add hydrogens necessary for "
                                     "identifying clashes and has also "
                                     "identified the following residues "
                                     "as needing to be flipped. Note that "
                                     "phenix.refine will often perform these "
                                     "flips by default.\n")
            else:
                group.addParam('showMesgNoSidechainFlips', LabelParam,
                                label="No sidechain flips required",
                                help="Asn, Gln, and His sidechains are "
                                     "asymmetric and may require flipping to "
                                     "form favorable van der Waals contacts "
                                     "and hydrogen bonding.\n\nREDUCE "
                                     "(phenix.reduce) has been run on your "
                                     "file to add hydrogens necessary for "
                                     "identifying clashes and has also "
                                     "identified the following residues "
                                     "as needing to be flipped. Note that "
                                     "phenix.refine will often perform these "
                                     "flips by default.\n")
            if (self.dictOverall['_n_omega_outliers'] > 0):
                group.addParam('showCisAndTwistedPeptides', LabelParam,
                                important=True,
                                label="Cis and Twisted peptides:",
                                help="Cis conformations are observed in "
                                     "about 5% of Prolines.\n\nCis "
                                     "conformations are observed in about "
                                     "0.03% of general residues.\n\nTwisted "
                                     "peptides are almost certainly "
                                     "modeling errors.\n")
            else:
                group.addParam('showMesgNoNonTransPeptides', LabelParam,
                                label="No non-trans peptides detected",
                               help="Cis conformations are observed in "
                                    "about 5% of Prolines.\n\nCis "
                                    "conformations are observed in about "
                                    "0.03% of general residues.\n\nTwisted "
                                    "peptides are almost certainly "
                                    "modeling errors.\n")
        if (self.dictOverall['_clashes'] == True):
            group = form.addGroup('Clashes')
            if (self.dictOverall['_n_clashes_outliers'] > 0):
                group.addParam('showClashes', LabelParam, important=True,
                                label="All atom-contact analysis",
                                help="This list summarizes all severe clashes "
                                     "(more than 0.4 Angstrom non-H-bond "
                                     "overlap) found by PROBE; you can view "
                                     "these graphically in Coot. If no "
                                     "hydrogens were present, REDUCE was "
                                     "used to add them prior to running "
                                     "PROBE.")
            else:
                group.addParam('showMesgNoClashes', LabelParam,
                                label="No bad contacts (> 0.4A overlap) found.",
                                help="This list summarizes all severe clashes "
                                     "(more than 0.4 Angstrom non-H-bond "
                                     "overlap) found by PROBE; you can view "
                                     "these graphically in Coot. If no "
                                     "hydrogens were present, REDUCE was "
                                     "used to add them prior to running "
                                     "PROBE.")
        if (self.dictOverall['_rna_group'] == True):
            group = form.addGroup('RNA restraint outliers')
            # TODO: Describe functions for RNA (Validation.Molprobity.py)
            if (self.dictOverall['_n_bonds_rna_outliers'] > 0):
                group.addParam('_showRNABonds', LabelParam,
                               important=True,
                               label="RNA nucleotides with excessive bond"
                                     " lengths:",
                               help="")
            else:
                group.addParam('showMesgNoBonds', LabelParam,
                               label="All bonds within expected limits",
                               help="")
            if (self.dictOverall['_n_angles_rna_outliers'] > 0):
                group.addParam('_showRNAAngles', LabelParam,
                               important=True,
                               label="RNA nucleotides with improper dihedral "
                                     "angles:",
                               help="")
            else:
                group.addParam('showMesgNoAngles', LabelParam,
                               label="All dihedral angles within expected "
                                     "limits",
                               help="")
            if (self.dictOverall['_n_puckers_rna_outliers'] > 0):
                group.addParam('_showRNASugarPuckers', LabelParam,
                               important=True,
                               label="RNA nucleotides with improper sugar "
                                     "puckers:",
                               help="")
            else:
                group.addParam('showMesgNoPuckers', LabelParam,
                               label="All sugar puckers in appropriate "
                                     "conformations",
                               help="")
            if (self.dictOverall['_n_suites_rna_outliers'] > 0):
                group.addParam('_showRNASuites', LabelParam,
                               important=True,
                               label="RNA nucleotides with bad backbone "
                                     "angles:",
                               help="")
            else:
                group.addParam('showMesgNosuites', LabelParam,
                               label="Backbone angles within expected ranges",
                               help="")
        form.addSection(label='Real-space correlation')
        group = form.addGroup('Real-space correlation to electron density')
        if (self.dictOverall['_overall_rsc'] == True):
            group.addParam('showMultiCriterionPlot', LabelParam,
                           important=True,
                           label="View Multi-criterion plot",
                           help="This plot shows simultaneously "
                                "Real-space correlation coefficients "
                                "and B-factor values for each residue.\n")
            group.addParam('showOverallRSCResults', LabelParam,
                            important=True,
                            label="Correlation Coefficients",
                            help="Real-space correlation\n\nFor a detailed "
                                 "definition of Mask CC, Volume CC and Peak "
                                 "CC, see Afonine, P. V., Klaholz, B. K., "
                                 "Moriarty, N. W., Poon, B. K., Sobolev, "
                                 "O. V., Terwilliger, T. C., Adams, P. D. "
                                 "& Urzhumtsev, A. (2018). bioRxiv. "
                                 "https://doi.org/10.1101/249607.\n\n"
                                 "Mask CC: Model-to-map correlation "
                                 "coefficient calculated in the map region "
                                 "around the model, using map values inside "
                                 "a mask calculated around the macromolecule"
                                 ".\n\nVolume CC and Peak CC compare only "
                                 "map regions with the highest density values "
                                 "and regions below a certain contouring "
                                 "threshold level are ignored.\nVolume CC: "
                                 "The map region considered is defined by "
                                 "the N highest points inside the molecular "
                                 "mask.\nPeak CC: In this case, calculations "
                                 "consider the union of regions defined by "
                                 "the N highest peaks in the model-calculated "
                                 "map and the N highest peaks in the "
                                 "experimental map.\n")
            self.residueTypeList = [u'Protein', u'Other', u'Water', u'Everything']
            group.addParam('residueType', EnumParam,
                            choices=self.residueTypeList, default=0,
                            label="Residue Type:",
                            help="Select a residue Type. Protein is chosen by "
                                 "default.\n")
            self.ccBelowList = [u'0.1', u'0.2', u'0.3', u'0.4', u'0.5', u'0.6',
                                u'0.7', u'0.8', u'0.9', u'1.0']
            group.addParam('ccIndex', EnumParam,
                            choices=self.ccBelowList, default=7,
                            label="Show CC below:",
                            help="Select a decimal. CC 0.8 is chosen by "
                                 "default.")
            group.addParam('showCCTable', LabelParam,
                            label="Table of Real-space correlation "
                                  "coefficients",
                            help="Residue: Protein chain, "
                                 "aminoacid index, and aminoacid name.\n"
                                 "B_iso: Isotropic B-factor value for "
                                 "residue (Disorder measure that quantitates "
                                 "the level of thermal motion).\nOccupancy: "
                                 "Disorder measure to indicate alternate "
                                 "location of multiple side chain atoms. "
                                 "Occupancies normally add up to unity. If "
                                 "no alternates have been observed, the "
                                 "occupancy equals to 1.00.\n2Fo-Fc: Electron "
                                 "density computed as 2Fo (Fo: structure "
                                 "factors from the diffraction patterns; "
                                 "experimental data) minus Fc (structure "
                                 "factor amplitudes calculated from model "
                                 "phases).\nFmodel: Model structure factors "
                                 "including all scales.\nCC: Correlation "
                                 "coefficent between the observed map and "
                                 "model-derived map in the real space.\n")
        group = form.addGroup('Fourier shell correlations')
        if (self.dictOverall['_fsc'] == True):
            _atomRadius = "%0.3f" % self.dictOverall['_atom_radius']
            group.addParam('showAtomRadius', LabelParam,
                            important=True,
                            label="Atom Mask Radius (Angstroms): " +
                                str(_atomRadius),
                            help='Radius of the "Fourier Shell", a spherical '
                                 'volume mask in Fourier space.\n')
            group.addParam('displayFSCplot', LabelParam,
                            label="Fourier Shell Correlation plot",
                            help="FSC regarding spatial frequency "
                                 "(1/Angstroms)")
        form.addSection(label='Atomic properties')
        group = form.addGroup('Occupancies')
        self._computeOccBFactor(self.dictOverall['_occ_bf_outliers'])
        if (len(self.occupancyList) == 0):
            group.addParam('showMesgNoOccupancies', LabelParam,
                            important=True, label="All occupancies okay")
        else:
            group.addParam('showOccupancies', LabelParam,
                                important=True,
                                label="Table of occupancies",
                                help="The table shows occupancy values lower "
                                 "than 1.0 for an atom when summed over "
                                 "alternative conformations (Type = atom) "
                                 "or for a residue when the occupancies for "
                                 "the atoms in a residue are not the same ("
                                 "Type = residue). Please check"
                                 " that they are correct.\n")
        group = form.addGroup('B-factor/ADPs')
        if ((self.dictOverall['_n_aniso_h'] is not None) and
            (self.dictOverall['_n_aniso_h'] > 0)):
           group.addParam('showWarningAnisoH', LabelParam, important= True,
                          label="WARNING: %d hydrogens have anisotropic "
                                "B-factors." % self.dictOverall['_n_aniso_h'])
        group.addParam('showIsotropicB', LabelParam,
                       important=True,
                       label="Isotropic B:",
                       help="The refinement against the map of B-factors or "
                            "Atomic Displacement Parameters (ADPs) is "
                            "performed at the last macro-cycle and using "
                            "reciprocal space.\n\nIsotropic B: Temperature "
                            "factor or displacement parameter constraining "
                            "the motion of the atoms so it is the same in "
                            "all three directions.\nMinimum, Maximum and Mean "
                            "are the statistic values of isotropic B-factor  ")
        if (len(self.suspiciousBFList) > 0):
            group.addParam('showSuspiciousBfactors', LabelParam,
                                label="Suspicious B-factors",
                                help= "The table of Suspicious B-factors (ADP"
                                " outliers) details all isotropic ADPs with "
                                "values outside a range of plus or minus 4 "
                                "sigmas around the mean value for the "
                                "considered atomic structure. We recommend "
                                "checking the atomic positions, occupancies "
                                "or element types shown in this table.\n")
            if self.dictOverall['_n_zero_b'] > 0:
                group.addParam('showWarningSuspicious', LabelParam,
                                   important=True,
                                   label="WARNING: %d atom(s) have "
                                         "B-factor(s) of zero, which indicates"
                                         " that either the occupancy or the "
                                         "element is incorrect."
                                         % self.dictOverall['_n_zero_b'])

    def _getVisualizeDict(self):
        return{
            'displayMapModel': self._displayMapModel,
            'showMolProbityResults': self._visualizeMolProbityResults,
            'showCootOutliers': self._showCootOutliers,
            'showMissingAtoms': self._showMissingAtoms,
            'showBLrestraints': self._showBLrestraints,
            'showBLoutliers': self._showBLoutliers,
            'showBArestraints': self._showBArestraints,
            'showBAoutliers': self._showBAoutliers,
            'showDArestraints': self._showDArestraints,
            'showDAoutliers': self._showDAoutliers,
            'showCHILrestraints': self._showCHILrestraints,
            'showCHILoutliers': self._showCHILoutliers,
            'showPLANARrestraints': self._showPLANARrestraints,
            'showPLANARoutliers': self._showPLANARoutliers,
            'showPlotType': self._showPlotType,
            'showRamaOutliersTable': self._showRamaOutliersTable,
            'showRotaOutliersTable': self._showRotaOutliersTable,
            'showCbetaOutliersTable': self._showCbetaOutliersTable,
            'showBackAsnGlnHisSidechains': self._showBackAsnGlnHisSidechains,
            'showCisAndTwistedPeptides': self._showCisAndTwistedPeptides,
            'showMultiCriterionPlot': self._showMultiCriterionPlot,
            'showOverallRSCResults': self._showOverallRSCResults,
            'showClashes': self._showClashes,
            'showCCTable': self._showCCTable,
            'displayFSCplot': self._displayFSCplot,
            'showOccupancies' : self._showOccupancies,
            'showIsotropicB': self._showIsotropicB,
            'showSuspiciousBfactors': self. _showSuspiciousBfactors
        }

    def _displayMapModel(self, e=None):
        bildFileName = os.path.abspath(self.protocol._getTmpPath(
            "axis_output.bild"))
        try:
            _inputVol = self.protocol.inputVolume.get()
        except:
            _inputVol = self.protocol.inputStructure.get().getVolume()

        if _inputVol is not None:
            dim = _inputVol.getDim()[0]
            sampling = _inputVol.getSamplingRate()

        else:
            # To show pdbs only
            dim = 150.
            sampling = 1.

        createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)
        counter = 0
        fnCmd = self.protocol._getTmpPath("chimera_output.cmd")
        f = open(fnCmd, 'w')
        # reference axis model = 0
        f.write("open %s\n" % bildFileName)

        # input 3D map
        counter += 1  # 1
        if _inputVol is not None:
            fnVol = self._getInputVolume()
        if fnVol is not None:
            try:
                VOLUMEFILENAME = os.path.abspath(self.protocol._getExtraPath(
                    self.protocol.MOLPROBITYFILE))
            except:
                VOLUMEFILENAME = os.path.abspath(self.protocol._getExtraPath(
                    self.protocol.REALSPACEFILE))
        f.write("open %s\n" % VOLUMEFILENAME)
        x, y, z = fnVol.getOrigin(force=True).getShifts()
        sampling = fnVol.getSamplingRate()
        f.write("volume #%d style surface voxelSize %f\nvolume #%d origin "
                "%0.2f,%0.2f,%0.2f\n" % (counter, sampling, counter, x, y, z))

        # input PDB (usually from coot)
        counter += 1  # 2
        pdbFileName = os.path.abspath(
            self.protocol.inputStructure.get().getFileName())
        f.write("open %s\n" % pdbFileName)

        # refined PDB
        if (len(os.listdir(self.protocol._getExtraPath())) > 5):
            counter += 1  # 3
            pdbFileName = os.path.abspath(self.protocol.outputPdb.getFileName())
            f.write("open %s\n" % pdbFileName)

        f.close()
        # run in the background
        runChimeraProgram(getProgram(), fnCmd + "&")
        return []

    def _visualizeMolProbityResults(self, e=None):
        headerList = ['statistic', 'value']
        dictX = self.dictSummary
        val = 0.4
        mesg = "Model Final Statistics"
        title = "MolProbity: Final Results Summary"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showCootOutliers(self, e=None):
        MOLPROBITYCOOTFILENAME = self.protocol._getExtraPath(
            self.protocol.MOLPROBITYCOOTFILENAME)
        args = ""
        args += " --script " + MOLPROBITYCOOTFILENAME
        # pdb file
        if (len(os.listdir(self.protocol._getExtraPath())) > 5):
            pdb = os.path.abspath(self.protocol.outputPdb.getFileName())
        else:
            pdb = os.path.abspath(
                self.protocol.inputStructure.get().getFileName())
        if pdb.endswith(".pdb"):
            args += " " + "--pdb " + pdb
        elif pdb.endswith(".cif"):
            args += " " + "--coords " + pdb
        # volume
        vol = self._getInputVolume()
        if vol is not None:
            try:
                VOLUMEFILENAME = os.path.abspath(self.protocol._getExtraPath(
                    self.protocol.MOLPROBITYFILE))
            except:
                VOLUMEFILENAME = os.path.abspath(self.protocol._getExtraPath(
                    self.protocol.REALSPACEFILE))
            args += " " + "--map " + VOLUMEFILENAME
        # run coot with args
        runCCP4Program(getProgram(self.COOT), args)

    def _getInputVolume(self):
        if self.protocol.inputVolume.get() is None:
            fnVol = self.protocol.inputStructure.get().getVolume()
        else:
            fnVol = self.protocol.inputVolume.get()
        return fnVol

    def _showMissingAtoms(self, e=None):
        headerList = ['Chain', 'Residue', 'AtLoc', 'Missing atoms']
        dataList = []
        chain = []
        residue = []
        atLoc = []
        missingAtoms = []
        for i in range(len(self.dictOverall['_missing_atoms'])):
            chain.append(self.dictOverall['_missing_atoms'][i][0])
            residue.append(self.dictOverall['_missing_atoms'][i][1])
            atLoc.append(self.dictOverall['_missing_atoms'][i][2])
            missingAtoms.append(self.dictOverall['_missing_atoms'][i][3])
        for c, r, aL, mA in zip(chain, residue, atLoc, missingAtoms):
            dataList.append((c, r, aL, mA))

        mesg = "Missing atoms"
        title = "Model missing atoms"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showBLrestraints(self, e=None):
        headerList = ['measure', 'value']
        dictX = self.dictBLRestraints
        val = 0.3
        mesg = "Bond Length Restraints\n(Deviations from ideal values)"
        title = "MolProbity: Basic Geometry"
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
        mesg = "Bond Angle Restraints\n(Deviations from ideal values)"
        title = "MolProbity: Basic Geometry"
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
        mesg = "Dihedral Angle Restraints\n(Deviations from ideal values)"
        title = "MolProbity: Basic Geometry"
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
        mesg = "Chirality Restraints\n(Deviations from ideal values)"
        title = "MolProbity: Basic Geometry"
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
        mesg = "Planarity Restraints\n(Deviations from ideal values)"
        title = "MolProbity: Basic Geometry"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showPLANARoutliers(self, e=None):
        headerList = ['Atoms', 'Max. delta', 'RMS (delta)', 'Deviation ('
                                                            'sigmas)']
        dataList = self.planarDataList
        mesg = "List of outliers (sorted by deviation)"
        title = "Planarity Restraints"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showPlotType(self, e=None):
        plot_index = int(self.plotType)
        self.listName = self.plotList[plot_index]
        self._writeCommand(self.listName)
        if self.listName == self.plotList[1]:
            self.TMPFILENAME = self.protocol._getTmpPath(self.ROTATMPFILE)
        else:
            self.TMPFILENAME = self.protocol._getTmpPath(self.RAMATMPFILE)
        with open(self.TMPFILENAME, "w") as f:
            f.write(self.command)
        # execute file with phenix.python
        runPhenixProgram("", self.TMPFILENAME)

    def _showRamaOutliersTable(self, e=None):
        headerList = self.dictOverall['_rama_headers']
        dataList = self.dictOverall['_rama_outliers']
        mesg = "Ramachandran outliers"
        title = "Ramachandran analysis"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showRotaOutliersTable(self, e=None):
        headerList = self.dictOverall['_rota_headers']
        dataList = self.dictOverall['_rota_outliers']
        mesg = "Rotamer outliers"
        title = "Rotamer analysis"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showCbetaOutliersTable(self, e=None):
        headerList = self.dictOverall['_cbeta_headers']
        dataList = self.dictOverall['_cbeta_outliers']
        mesg = "C-beta position outliers"
        title = "C-beta deviation analysis"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showBackAsnGlnHisSidechains(self, e=None):
        headerList = self.dictOverall['_nqh_flips_headers']
        dataList = self.dictOverall['_nqh_flips_outliers']
        mesg = "Recommended sidechain flips"
        title = "Backwards Asn/Gln/His sidechains"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showCisAndTwistedPeptides(self, e=None):
        headerList = self.dictOverall['_omega_headers']
        dataList = self.dictOverall['_omega_outliers']
        mesg = "Cis and Twisted peptides"
        title = "Cis and Twisted peptides analyis"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showClashes(self, e=None):
        headerList = self.dictOverall['_clashes_headers']
        dataList = self.dictOverall['_clashes_outliers']
        mesg = "Bad contacts from PROBE: %d overlapping atom pairs" \
               % len(dataList)
        title = "All atom-contact analyis"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showMultiCriterionPlot(self, e=None):
        self.listName = "Multi-criterion plot"
        self._writeCommand(self.listName)
        self.TMPFILENAME = self.protocol._getTmpPath(self.MULTICPLOTTMPFILE)
        with open(self.TMPFILENAME, "w") as f:
            f.write(self.command)
        # execute file with phenix.python
        runPhenixProgram("", self.TMPFILENAME)


    def _showOverallRSCResults(self, e=None):
        headerList = ['statistic', 'value']
        dictX = self.dictOverall
        val = 0.3
        mesg = "Model Final Statistics"
        title = "Real-space correlation: Final Results Summary"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showCCTable(self, e=None):
        self._computeCCTable()

        headerList = self.dictOverall['_rs_headers']
        dataList = self.RSCCList
        mesg = "Real-space correlation"
        title = "Correlation coefficients table"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showOccupancies(self, e=None):
        headerList = self.dictOverall['_occ_bf_headers']
        dataList = self.occupancyList
        mesg = "Occupancies"
        title = "Occupancies"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showIsotropicB(self, e=None):
        headerList = ['statistic', 'value']
        dictX = collections.OrderedDict()
        dictX['Minimun'] = self.dictOverall['_b_min']
        if self.dictOverall['_b_min_macromolecules'] is not None:
            dictX['Min (macromolecules)'] = \
                self.dictOverall['_b_min_macromolecules']
        if self.dictOverall['_b_min_ligands'] is not None:
            dictX['Min (ligands)'] = \
                self.dictOverall['_b_min_ligands']
        dictX['Maximun'] = self.dictOverall['_b_max']
        if self.dictOverall['_b_max_macromolecules'] is not None:
            dictX['Max (macromolecules)'] = \
                self.dictOverall['_b_max_macromolecules']
        if self.dictOverall['_b_max_ligands'] is not None:
            dictX['Max (ligands)'] = \
                self.dictOverall['_b_max_ligands']
        dictX['Mean'] = self.dictOverall['_b_mean']
        if self.dictOverall['_b_mean_macromolecules'] is not None:
            dictX['Mean (macromolecules)'] = \
                self.dictOverall['_b_mean_macromolecules']
        if self.dictOverall['_b_mean_ligands'] is not None:
            dictX['Mean (ligands)'] = \
                self.dictOverall['_b_mean_ligands']
        val = 0.2
        mesg = "Isotropic B"
        title = "B-factors/ADPs"
        self._showResults(headerList, dictX, val, mesg, title)

    def _showSuspiciousBfactors(self, e=None):
        headerList = self.dictOverall['_occ_bf_headers']
        dataList = self.suspiciousBFList
        mesg = "Suspicious B-factors"
        title = "B-factors/ADPs"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showResults(self, headerList, dictX, val, mesg, title):
        dataList = []
        for k, v in dictX.iteritems():
            if k[0] == "_":
                continue
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

    def _showOutliers(self, headerList, dataList, mesg, title):

        if not dataList:
            errorWindow(self.getTkRoot(), "No data available")
            return

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg=mesg,
                  title=title,
                  height=min(20,len(dataList)), width=250, padding=40)

    def _displayFSCplot(self, e=None):
        xList = self.dictOverall['_x_fsc']
        yList = self.dictOverall['_y_fsc']

        if not (xList or yList):
            errorWindow(self.getTkRoot(), "No data available")
            return

        title = 'Fourier Shell Correlation (Map vs. Model)'
        plt.plot(xList, yList)
        plt.axis([0, 0.8, -0.2, 1])
        plt.title(title)
        plt.xlabel('1/resolution (1/Angstrom)')
        plt.ylabel('FSC')
        plt.show()

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

    def _writePickleData(self):
        ANALYSISTMPFILENAME = self.protocol._getTmpPath(
            self.ANALYSISTMPFILE)
        command = """import pickle
import collections
import json

def pickleData(file):
    with open(file,"r") as f:
        return pickle.load(f)

# process file %s"
data = pickleData('%s')
dictOverall = collections.OrderedDict()

# missing atoms
if data.missing_atoms is not None:
    dictOverall['_missing_atoms'] = data.missing_atoms
    dictOverall['_len_missing_atoms'] = len(data.missing_atoms)
else:
    dictOverall['_len_missing_atoms'] = 0

# Protein group
if data.ramalyze is None:
    dictOverall['_protein'] = False
else:
    dictOverall['_protein'] = True

    # Ramachandran analysis
    dictOverall['_percent_rama_outliers'] = data.ramalyze.percent_outliers
    dictOverall['_rama_outliers'] = data.ramalyze.as_gui_table_data()
    dictOverall['_rama_headers'] = data.ramalyze.gui_list_headers

    # Rotamer analysis
    dictOverall['_percent_rota_outliers'] = data.rotalyze.percent_outliers
    dictOverall['_rota_outliers'] = data.rotalyze.as_gui_table_data()
    dictOverall['_rota_headers'] = data.rotalyze.gui_list_headers

    # C-beta outliers
    dictOverall['_n_cbeta_outliers'] = data.cbetadev.n_outliers
    dictOverall['_cbeta_outliers'] = data.cbetadev.as_gui_table_data()
    dictOverall['_cbeta_headers'] = data.cbetadev.gui_list_headers

    # Backwards Asn/Gln/His sidechains
    dictOverall['_n_nqh_flips_outliers'] = data.nqh_flips.n_outliers
    dictOverall['_nqh_flips_outliers'] = data.nqh_flips.as_gui_table_data()
    dictOverall['_nqh_flips_headers'] = data.nqh_flips.gui_list_headers

    # Cis and Twisted peptides
    dictOverall['_n_omega_outliers'] = data.omegalyze.n_outliers
    dictOverall['_omega_outliers'] = data.omegalyze.as_gui_table_data()
    dictOverall['_omega_headers'] = data.omegalyze.gui_list_headers

# RNA group
if data.rna is None:
    dictOverall['_rna_group'] = False
else:
    dictOverall['_rna_group'] = True
    dictOverall['_n_bonds_rna_outliers'] = data.rna.bonds.n_outliers
    dictOverall['_n_angles_rna_outliers'] = data.rna.angles.n_outliers
    dictOverall['_n_puckers_rna_outliers'] = data.puckers.angles.n_outliers
    dictOverall['_n_suites_rna_outliers'] = data.suites.angles.n_outliers

# Clashes
if data.clashes is None:
    dictOverall['_clashes'] = False
else:
    dictOverall['_clashes'] = True
    dictOverall['_n_clashes_outliers'] = data.clashes.n_outliers
    dictOverall['_clashes_outliers'] = data.clashes.as_gui_table_data()
    dictOverall['_clashes_headers'] = data.clashes.gui_list_headers

# correlation coefficients
if data.real_space is None:
    dictOverall['_overall_rsc'] = False
    dictOverall['_fsc'] = False
else:
    if data.real_space.overall_rsc is None:
        dictOverall['_overall_rsc'] = False
    else:
        dictOverall['_overall_rsc'] = True
        dictOverall['Mask CC'] = data.real_space.overall_rsc[0]
        dictOverall['Volume CC'] = data.real_space.overall_rsc[1]
        dictOverall['Peak CC'] = data.real_space.overall_rsc[2]

        # real_space_correlation_coefficients table
        dictOverall['_rs_protein'] = []
        for i in range(len(data.real_space.protein)):
            dictOverall['_rs_protein'].append(data.real_space.protein[
            i].as_table_row_phenix())
        dictOverall['_rs_other'] = [] 
        for i in range(len(data.real_space.other)):
            dictOverall['_rs_other'].append(data.real_space.other[
            i].as_table_row_phenix())
        dictOverall['_rs_water'] = []
        for i in range(len(data.real_space.water)):
            dictOverall['_rs_water'].append(data.real_space.water[
            i].as_table_row_phenix())
        dictOverall['_rs_everything'] = []
        for i in range(len(data.real_space.everything)):
            dictOverall['_rs_everything'].append(data.real_space.everything[
            i].as_table_row_phenix())
        dictOverall['_rs_headers'] = data.real_space.gui_list_headers

    # fsc
    if data.real_space.fsc is None:
        dictOverall['_fsc'] = False
    else:
        dictOverall['_fsc'] = True

        # atom mask radius
        dictOverall['_atom_radius'] = data.real_space.fsc.atom_radius

        # fsc graph data
        x_elements = []
        y_elements = []
        for x in data.real_space.fsc.d_inv:
            x_elements.append(x)
        for y in data.real_space.fsc.fsc:
            y_elements.append(y)
        dictOverall['_x_fsc'] = x_elements
        dictOverall['_y_fsc'] = y_elements

# occupancy and suspicious B-factors table
dictOverall['_occ_bf_outliers'] = data.model_stats.all.as_gui_table_data()
dictOverall['_occ_bf_headers'] = data.model_stats.all.gui_list_headers

# isotropic B (B-factors/ADPs)
dictOverall['_n_aniso_h'] = data.model_stats.all.n_aniso_h
dictOverall['_n_zero_b'] = data.model_stats.all.n_zero_b
dictOverall['_b_min'] = data.model_stats.all.b_min
dictOverall['_b_max'] = data.model_stats.all.b_max
dictOverall['_b_mean'] = data.model_stats.all.b_mean
if data.model_stats.macromolecules is not None:
    dictOverall['_b_min_macromolecules'] = data.model_stats.macromolecules.b_min
    dictOverall['_b_max_macromolecules'] = data.model_stats.macromolecules.b_max
    dictOverall['_b_mean_macromolecules'] = data.model_stats.macromolecules.b_mean
if data.model_stats.ligands is not None:
    dictOverall['_b_min_ligands'] = data.model_stats.ligands.b_min
    dictOverall['_b_max_ligands'] = data.model_stats.ligands.b_max
    dictOverall['_b_mean_ligands'] = data.model_stats.ligands.b_mean
""" % (self.MOLPROBITYPKLFILENAME, self.MOLPROBITYPKLFILENAME)

        command += """with open('%s',"w") as f:
    f.write(json.dumps(dictOverall))
""" % (ANALYSISTMPFILENAME)

        pythonFileName = ANALYSISTMPFILENAME.replace('.txt', '.py')
        # write script file
        with open(pythonFileName, "w") as f:
            f.write(command)

        # execute file with phenix.python
        runPhenixProgram("", pythonFileName)

        # read file in scipion python
        with open(ANALYSISTMPFILENAME, "r") as f:
            self.dictOverall = f.read()
            # self.dataDict = json.loads(f.read())

        self._store()

    def _writeCommand(self, listName):
        self.command ="""import pickle
        
def pickleData(file):
    with open(file,"r") as f:
        return pickle.load(f)

# process file %s"
data = pickleData('%s')
"""% (self.MOLPROBITYPKLFILENAME, self.MOLPROBITYPKLFILENAME)
        if (listName == "Multi-criterion plot"):
            self.command += """#RSC section
if data.real_space.overall_rsc is not None:
"""
        else:
            self.command += """# Protein group
if data.ramalyze is not None:
"""
        self.command +="""
    try :
        import wxtbx.app
    except ImportError, e :
        raise Sorry("wxPython not available.")
    else :
        app = wxtbx.app.CCTBXApp(0)
"""
        # TODO: Get appropriate functioning of data._multi_criterion.display_wx_plots()
        if (listName == "Multi-criterion plot"):
            self.command += """
        # data._multi_criterion.display_wx_plots() 
        data.display_wx_plots()
"""
        elif (listName == "Ramachandran plot"):
            self.command += """
        data.ramalyze.display_wx_plots()
"""
        elif (listName == "Chi1-Chi2 plot"):
            self.command += """
        data.rotalyze.display_wx_plots()
"""
        self.command +="""
        app.MainLoop()
"""

    def _computeOccBFactor(self, listL):
        self.suspiciousBFList = []
        self.occupancyList = []
        if not listL:
            pass
        else:
            for i in range(len(listL)):
                if listL[i][2] == 1.:
                    self.suspiciousBFList.append(listL[i])
                else:
                    self.occupancyList.append(listL[i])

    def _computeCCTable(self):
        decimal_index = int(self.ccIndex)
        residueType_index = int(self.residueType)
        self.RSCCList = []
        listL = []
        listName = self.residueTypeList[residueType_index]
        if listName == self.residueTypeList[0]:
            listL = self.dictOverall['_rs_protein']
        elif listName == self.residueTypeList[1]:
            listL = self.dictOverall['_rs_other']
        elif listName == self.residueTypeList[2]:
            listL = self.dictOverall['_rs_water']
        elif listName == self.residueTypeList[3]:
            listL = self.dictOverall['_rs_everything']
        if not listL:
            pass
        else:
            for i in range(len(listL)):
                if listL[i][5] <= float(self.ccBelowList[decimal_index]):
                    self.RSCCList.append((unicode(listL[i][0]),
                                          unicode(listL[i][1]),
                                          unicode(listL[i][2]),
                                          unicode(listL[i][3]),
                                          unicode(listL[i][4]),
                                          unicode(listL[i][5])))




