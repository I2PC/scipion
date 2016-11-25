# coding: latin-1
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
Bibtex string file for IGBMC package.
"""

_bibtexStr = """

@Article{Hoang2013,
  Title      = {gEMpicker: a highly parallel GPU-accelerated particle picking tool for cryo-electron microscopy},
  Author     = {Hoang, Thai and Cavin, Xavier and Schultz, Patrick and Ritchie, David},
  Journal    = {BMC Structural Biology},
  Year       = {2013},
  Pages      = {25},
  Volume     = {13},
  Number     = {1}
  Abstract   = {BACKGROUND:Picking images of particles in cryo-electron micrographs is an important step in solving the 3D structures of large macromolecular assemblies. However, in order to achieve sub-nanometre resolution it is often necessary to capture and process many thousands or even several millions of 2D particle images. Thus, a computational bottleneck in reaching high resolution is the accurate and automatic picking of particles from raw cryo-electron micrographs.RESULTS:We have developed "gEMpicker", a highly parallel correlation-based particle picking tool. To our knowledge, gEMpicker is the first particle picking program to use multiple graphics processor units (GPUs) to accelerate the calculation. When tested on the publicly available keyhole limpet hemocyanin dataset, we find that gEMpicker gives similar results to the FindEM program. However, compared to calculating correlations on one core of a contemporary central processor unit (CPU), running gEMpicker on a modern GPU gives a speed-up of about 27 x. To achieve even higher processing speeds, the basic correlation calculations are accelerated considerably by using a hierarchy of parallel programming techniques to distribute the calculation over multiple GPUs and CPU cores attached to multiple nodes of a computer cluster. By using a theoretically optimal reduction algorithm to collect and combine the cluster calculation results, the speed of the overall calculation scales almost linearly with the number of cluster nodes available.CONCLUSIONS:The very high picking throughput that is now possible using GPU-powered workstations or computer clusters will help experimentalists to achieve higher resolution 3D reconstructions more rapidly than before.},
  Doi        = {10.1186/1472-6807-13-25},
  Language   = {eng},
  PubMedID   = {24144335},
  ISSN       = {1472-6807},
  Url        = {http://www.biomedcentral.com/1472-6807/13/25}
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
