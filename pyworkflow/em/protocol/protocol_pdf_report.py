# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

import os
import glob

from pyworkflow.protocol.params import PathParam
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.utils.path import moveFile, copyFile


class ProtPDFReport(EMProtocol):
    """ 
    Produce a pdf report from the files in a given directory.
    Supported file formats: *.tex, *.txt, *.jpg, *.png, *.pdf
    Files in the directory are sorted by name alphabetically,
    so if you want them to have the right order a possibility is to name them as
    0010-myText.txt
    0020-aFigure.png
    0030-anotherFigure.jpg
    0040-aPaper.pdf
    0050-anotherText.tex
    ...
    when these files are sorted, they will be sorted by the number in front.
    """    
    _label = 'pdf report'

    #--------------------------- DEFINE param functions ------------------------
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('filesPath', PathParam,
                      label="Files directory",
                      help="Directory with the input files. \n"
                           "Check protocol help for more details.")
    
    #--------------------------- INSERT steps functions ------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('createReport')        

    #--------------------------- STEPS functions -------------------------------
    def createReport(self):
        fnTex = "report.tex"
        fhTex = open(self._getExtraPath(fnTex), "w")
        template ="""
\\documentclass[12pt]{article}
\\usepackage{amsmath,amsthm,amssymb,amsfonts} 
\\usepackage{graphicx}
\\usepackage{pdfpages}
\\DeclareGraphicsExtensions{.pdf,.png,.jpg}
\\begin{document}
\\title{User Report}
\\author{Created by Scipion}
\\maketitle
"""
        fhTex.write(template)
        
        fnDir = self.filesPath.get()

        if not os.path.isdir(fnDir):
            fnDir = os.path.basename(fnDir)

        for fileName in sorted(glob.glob(os.path.join(fnDir,"*"))):
            fnDest = os.path.basename(fileName).lower()
            fnDest = fnDest.replace(" ","_")
            fnDest = fnDest.replace(":","_")
            fnDest = fnDest.replace(";","_")
            copyFile(fileName, self._getExtraPath(fnDest))

            if fnDest.endswith(".tex") or fnDest.endswith(".txt"):
                fhTex.write("\\input{%s}\n" % fnDest)
                fhTex.write("\n")
            elif fnDest.endswith(".png") or fnDest.endswith(".jpg"):
                fhTex.write("\\begin{center}\n")
                fhTex.write("\\includegraphics[width=14cm]{%s}\n" % fnDest)
                fhTex.write("\\end{center}\n")
                fhTex.write("\n")
            elif fnDest.endswith(".pdf"):
                fhTex.write("\\includepdf[pages=-]{%s}\n" % fnDest)
                fhTex.write("\\clearpage\n")
                fhTex.write("\n")
        
        template = """ 
\\end{document}
"""
        fhTex.write(template)
        fhTex.close()

        args = "-interaction=nonstopmode " + fnTex
        self.runJob("pdflatex", args, cwd=self._getExtraPath())

        fnPDF = self._getExtraPath("report.pdf")

        if os.path.exists(fnPDF):
            moveFile(fnPDF,self._getPath("report.pdf"))
        else:
            raise Exception("PDF file was not produced.")

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append("Input directory: %s" % self.filesPath.get())
        summary.append("Output report: %s" % self._getPath("report.pdf"))
        return summary
    

