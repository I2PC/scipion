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

import os

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER
from pyworkflow.em.data import Sequence
from pyworkflow.gui.text import openTextFileEditor


class SequenceViewer(Viewer):
    """ Wrapper to visualize Sequences with an editor. """
    _environments = [DESKTOP_TKINTER]
    _targets = [Sequence]

    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    # def visualize(self, obj, **kwargs):
        # fn = obj.getFileName()
        # openTextFileEditor(fn)
    def visualize(self, obj, **kwargs):
        # The sequence object is visualized in a tmp file;
        # To build this tmp file we use a method (saveFile) from the class
        # SequenceHandler that requires a Biopython sequence (Bio.Seq.Seq
        # object);

        # Step 1: transformation of the sequence of our Scipion Sequence
        # object (obj.getSequence()) in a Biopython Sequence:
        # Let keep the SequenceHandler import here to avoid a default BioPython
        # import
        from pyworkflow.em.convert.sequence import SequenceHandler
        seqHandler = SequenceHandler(obj.getSequence(),
                                     isAminoacid=obj.getIsAminoacids())
        seqBio = seqHandler._sequence  # Bio.Seq.Seq object
        # Step 2: retrieving of the other args needed in the saveFile method
        seqID = obj.getId()
        seqName = obj.getSeqName()
        seqDescription = obj.getDescription()
        seqFileName = os.path.abspath(
            self.protocol._getTmpPath(seqName + ".fasta"))
        # Step 3: Sequence saved in the tmp file
        seqHandler.saveFile(seqFileName, seqID, sequence=seqBio,
                            name=seqName, seqDescription=seqDescription,
                            type="fasta")
        # Step 4: Visualization of tmp file
        openTextFileEditor(seqFileName)