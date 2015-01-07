# coding: latin-1
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for Simple package.
"""

_bibtexStr = """

@Article{Chen2013,
  Title                    = {Fast and accurate reference-free alignment of subtomograms.},
  Author                   = {Chen, Y. and Pfeffer, S. and Hrabe, T. and Schuller, J. M. and Förster, F.},
  Journal                  = {J. Structural Biology},
  Year                     = {2013},
  Month                    = {Mar},
  Pages                    = {235--245},
  Volume                   = {182},
  Abstract                 = {In cryoelectron tomography alignment and averaging of subtomograms, each dnepicting the same macromolecule, improves the resolution compared to the individual subtomogram. Major challenges of subtomogram alignment are noise enhancement due to overfitting, the bias of an initial reference in the iterative alignment process, and the computational cost of processing increasingly large amounts of data. Here, we propose an efficient and accurate alignment algorithm via a generalized convolution theorem, which allows computation of a constrained correlation function using spherical harmonics. This formulation increases computational speed of rotational matching dramatically compared to rotation search in Cartesian space without sacrificing accuracy in contrast to other spherical harmonic based approaches. Using this sampling method, a reference-free alignment procedure is proposed to tackle reference bias and overfitting, which also includes contrast transfer function correction by Wiener filtering. Application of the method to simulated data allowed us to obtain resolutions near the ground truth. For two experimental datasets, ribosomes from yeast lysate and purified 20S proteasomes, we achieved reconstructions of approximately 20Å and 16Å, respectively. The software is ready-to-use and made public to the community.},
  Doi                      = {10.1016/j.jsb.2013.03.002},
  Language                 = {eng},
  Pii                      = {S1047-8477(13)00073-7},
  Pmid                     = {23523719},
  Review                   = {Very good article on subtomogram averaging with a function to compare volumes based on Spherical Harmonics},
  Url                      = {http://dx.doi.org/10.1016/j.jsb.2013.03.002}
}

@Article{Chen2014,
  Title                    = {Autofocused 3D classification of cryoelectron subtomograms.},
  Author                   = {Chen, Yuxiang and Pfeffer, Stefan and Fern{\'{a}}ndez, Jos{\'{e}} Jes{\'{u}}s and Sorzano, Carlos Oscar S. and F{\"{o}}rster, Friedrich},
  Journal                  = {Structure},
  Year                     = {2014},
  Month                    = {Oct},
  Number                   = {10},
  Pages                    = {1528--1537},
  Volume                   = {22},
  Abstract                 = {Classification of subtomograms obtained by cryoelectron tomography (cryo-ET) is a powerful approach to study the conformational landscapes of macromolecular complexes in situ. Major challenges in subtomogram classification are the low signal-to-noise ratio (SNR) of cryo-tomograms, their incomplete angular sampling, the unknown number of classes and the typically unbalanced abundances of structurally distinct complexes. Here, we propose a clustering algorithm named AC3D that is based on a similarity measure, which automatically focuses on the areas of major structural discrepancy between respective subtomogram class averages. Furthermore, we incorporate a spherical-harmonics-based fast subtomogram alignment algorithm, which provides a significant speedup. Assessment of our approach on simulated data sets indicates substantially increased classification accuracy of the presented method compared to two state-of-the-art approaches. Application to experimental subtomograms depicting endoplasmic-reticulum-associated ribosomal particles shows that AC3D is well suited to deconvolute the compositional heterogeneity of macromolecular complexes in situ.},
  Doi                      = {10.1016/j.str.2014.08.007},
  Institution              = {Department of Molecular Structural Biology, Max-Planck Institute of Biochemistry, Am Klopferspitz 18, 82152 Martinsried, Germany. Electronic address: foerster@biochem.mpg.de.},
  Language                 = {eng},
  Pii                      = {S0969-2126(14)00252-4},
  Pmid                     = {25242455},
  Url                      = {http://dx.doi.org/10.1016/j.str.2014.08.007}
}
"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)  
