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

import json
import os
from tkMessageBox import showerror
from protocol_emringer import PhenixProtRunEMRinger
from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.em.viewer import TableView
import collections
import glob
from PIL import Image
from convert import runPhenixProgram


def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except:
        print("Error:", msg)


class PhenixProtRunEMRingerViewer(ProtocolViewer):
    """ Viewer for Phenix program EMRinger
    """
    _label = 'EMRinger Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [PhenixProtRunEMRinger]
    EMRINGERTOTALTHRESH = 'Total.threshold_scan.png'
    EMRINGERSUBPLOTSFILENAME = 'residue_subplots.py'

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        self.dataDict = json.loads(str(self.protocol.stringDataDict),
                                   object_pairs_hook=collections.OrderedDict)
        self.plots = glob.glob(self.protocol._getExtraPath("*_plots"))[0]
        self.EMRINGERSUBPLOTSFILENAME = self.protocol._getTmpPath(
            self.EMRINGERSUBPLOTSFILENAME)

    def _defineParams(self, form):
        form.addSection(label='Visualization of EM Ringer results')
        # group = form.addGroup('Overall results')
        form.addParam('showFinalResults', LabelParam,
                      label="Summary Table of Results",
                      help="Table of Final Results\n\nOptimal Threshold: "
                           "Electron potential map cutoff value at which the "
                           "maximum EMRinger score was obtained.\nRotamer "
                           "Ratio: Percentage of rotameric residues at the "
                           "Optimal threshold value.\nMax Zscore: Z-score "
                           "computed to determine the significance of the "
                           "distribution at the Optimal threshold value.\n"
                           "Model Length: Total of non-gamma-branched, "
                           "non-proline aminoacids with a non-H gamma atom "
                           "used in global EMRinger score computation.\n"
                           "EMRinger Score: Maximum EMRinger score calculated "
                           "at the Optimal Threshold.\n")
        form.addParam('showThresholdScan', LabelParam,
                      label="Threshold Scan",
                      help="Statistics across all thresholds\n\n"
                           "Blue line: EMRinger Score regarding"
                           "the map electron potential. "
                           "The maximum value of EMRinger Score "
                           "determines the Optimal Threshold.\nRed "
                           "line: Percentage of rotameric residues "
                           "regarding the map electron potential.")

        self.thresList = [str("%0.3f" % thres) for thres in self.dataDict[
            '_thresholds']]
        _maxScoreIndex = self.dataDict['_maxScoreIndex']
        form.addParam('threshold', EnumParam,
                      choices=self.thresList,
                      default=_maxScoreIndex,
                      label="Density threshold",
                      help="Choose one of the 20 map density cutoff values "
                           "used to compute the percentage of rotameric "
                           "residues. The Optimal threshold, at which the "
                           "maximum EMRinger score was obtained, "
                           "is shown by default.\n"
                      )
        form.addParam('showPeakCount', LabelParam,
                      label="Peak count for the selected density threshold",
                      help="Histograms for rotameric (blue) and non-rotameric "
                           "(red) residues at the selected map density.")
        self.chainList = self.dataDict['_chains']
        form.addParam('chain', EnumParam,
                      choices=self.chainList,
                      default=0,
                      label="Chain",
                      help="Choose one of the model chains. The first chain "
                           "is shown by default.")
        form.addParam('showRollingWindows', LabelParam,
                      label="Rolling window for the selected chain.",
                      help="Visualize the rolling window EMRinger analysis of "
                           "your selected chain to distinguish regions of "
                           "improved model quality.This analysis was "
                           "performed on rolling sliding 21-residue windows "
                           "along the primary sequence of proteins.")
        self.residueList = self.dataDict['_residues_format']
        form.addParam('residue', EnumParam,
                      choices=self.residueList,
                      default=0,
                      label="Residue",
                      help="Choose one of the gamma-carbon-containing "
                           "residues (at least with one Chi angle) "
                           "located in one of the chains in the "
                           "position indicated. ")
        form.addParam('showRingerResults', LabelParam,
                      label="Ringer results for the selected residue.",
                      help="Individual plots for each Chi angle of the "
                           "selected residue. Numeric values are showed in "
                           "the extra/*.csv file.")

    def _getVisualizeDict(self):
        return {
            'showFinalResults': self._visualizeFinalResults,
            'showThresholdScan': self._showThresholdScan,
            'showPeakCount': self._showPeakCount,
            'showRollingWindows': self._showRollingWindows,
            'showRingerResults': self._showRingerResults
        }

    def _visualizeFinalResults(self, e=None):

        headerList = ['statistic', 'value']
        dataList = []
        for k, v in self.dataDict.iteritems():
            if k[0] == "_":
                continue
            elif isinstance(v, int):
                dataList.append((k, v))
            elif isinstance(v, float):
                dataList.append((k, "%0.3f" % v))

        if not dataList:
            errorWindow(self.getTkRoot(), "No data available")
            return

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg="Final Statistics for Model/Map Pair\n",
                  title="EMRinger: Final Results Summary",
                  height=len(dataList), width=250, padding=40)

    def showImage(self, fileName):
        fileName = os.path.join(self.plots, fileName)
        img = Image.open(fileName)
        img.show()

    def _showThresholdScan(self, e=None):
        self.showImage(self.EMRINGERTOTALTHRESH)

    def _showPeakCount(self, e=None):
        threshold_index = int(self.threshold)
        fileName = "%s.histogram.png" % self.thresList[threshold_index]
        self.showImage(fileName)

    def _showRollingWindows(self, e=None):
        chain_index = int(self.chain)
        fileName = "%s_rolling.png" % self.chainList[chain_index]
        self.showImage(fileName)

    def _showRingerResults(self, e=None):
        i = int(self.residue)
        index = int(self.dataDict['_residues_dict'][self.residueList[i]])
        # emringer main file
        mainDataFile = glob.glob(self.protocol._getExtraPath(
            "*_emringer.pkl"))[0]
        command = """from mmtbx.ringer.em_rolling import easy_pickle
import matplotlib.pyplot as plt

# process file '%s'
file_name='%s'
ringer_results = easy_pickle.load(file_name)
figure = plt.figure() #  (figsize=(20,1000))

def show_residue (residue, show_background_boxes=True) :
    subplots = []
    for i in range(1, residue.n_chi + 1) :
        chi = residue.get_angle(i)
        if (chi is None) : continue
        if (len(subplots) > 0) :
            p = figure.add_subplot(4, 1, i, sharex=subplots[0])
        else:
            p = figure.add_subplot(4, 1, i)
            p.set_title(residue.format())
        p.set_position([0.15, 0.725 - 0.225*(i-1), 0.8, 0.225])
        x = [ k*chi.sampling for k in range(len(chi.densities)) ]
        p.plot(x, chi.densities, 'r-', linewidth=1)
        if (chi.fofc_densities is not None) :
            p.plot(x, chi.fofc_densities, linestyle='--', color=[0.5,0.0,1.0])
        p.axvline(chi.angle_current, color='b', linewidth=2, linestyle='--')
        p.axhline(0, color=(0.4,0.4,0.4), linestyle='--', linewidth=1)
        if show_background_boxes:
            p.axhspan(0.3,1,facecolor="green",alpha=0.5)
            p.axhspan(-1,0.3,facecolor="grey",alpha=0.5)
        p.set_xlim(0,360)
        p.set_ylabel("Rho")
        p.set_xlabel("Chi" + str(i))
        subplots.append(p)
        plt.subplots_adjust(left = 0.18, bottom = 0.00, right = 0.94,
        top = 0.93,
                        wspace = 0.80, hspace = 0.43)
    plt.tight_layout()
    plt.show()
index = %d
show_residue(ringer_results[index])
""" % (mainDataFile, mainDataFile, index)

        with open(self.EMRINGERSUBPLOTSFILENAME, "w") as f:
            f.write(command)
        # execute file with phenix.python
        runPhenixProgram("", self.EMRINGERSUBPLOTSFILENAME)
