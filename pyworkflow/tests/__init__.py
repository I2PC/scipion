#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
#                Laura del Cano         (ldelcano@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
import os, sys

import pyworkflow as pw
from tests import *
from pyworkflow.utils.path import makeFilePath
import model

try:
    from unittest.runner import _WritelnDecorator # Python 2.7+
except ImportError:
    from unittest import _WritelnDecorator # Python <2.6

    
DataSet(name='xmipp_tutorial', folder='xmipp_tutorial', 
        files={'micsGoldSqlite': 'gold/micrographs_gold.sqlite',
               'micsGoldSqlite2': 'gold/micrographs2_gold.sqlite',
               'micsGoldXmd': 'gold/micrographs_gold.xmd',
               'ctfGold': 'gold/xmipp_ctf.ctfparam',
               'micsSqlite': 'micrographs/micrographs.sqlite',
               'coordsGoldSqlite': 'gold/coordinates_gold.sqlite', 
               'posSupervisedDir': 'pickingXmipp/pickedSupervised',
               'posAllDir': 'pickingXmipp/pickedAll',
               'boxingDir': 'pickingEman',
               'boxingFile': 'pickingEman/info/BPV_1386_info.json',
               'allMics': 'micrographs/*.mrc',
               'rctMicsU': 'rct/micrographs/*U.mrc',
               'rctMicsT': 'rct/micrographs/*T.mrc',
               'rctCoords': 'rct/pickingXmipp',
               'mic1': 'micrographs/BPV_1386.mrc',
               'mic2': 'micrographs/BPV_1387.mrc',
               'mic3': 'micrographs/BPV_1388.mrc',
               'particles': 'particles/*.hdf',
               'particles1': 'particles/BPV_1386_ptcls.hdf',
               'volumes': 'volumes/*.mrc',
               'vol1': 'volumes/BPV_scale_filtered_windowed_64.vol',
               'vol2': 'volumes/volume_1_iter_002.mrc',
               'vol3': 'volumes/volume_1_iter_002.mrc',
               'aligned_particles': 'gold/aligned_particles.sqlite',
               'images10': 'gold/images10.xmd',
               })


DataSet(name='mda', folder='hemoglobin_mda', 
        files={'particles': 'particles/*.spi',
               'volumes': 'volumes/*.spi'})


DataSet(name='tomo', folder='xmipp_tomo_test', 
        files={'volumes': 'volumes/*.spi',
               'vol1': 'volumes/subvols_6E6.001.mrc.spi',
               'vol2': 'volumes/subvols_6E6.002.mrc.spi',
               'vol3': 'volumes/subvols_6E6.003.mrc.spi'})


DataSet(name='relion_tutorial', folder='relion_tutorial', 
        files={'posAllDir': 'pickingXmipp',
               'boxingDir': 'pickingEman',
               'allMics': 'micrographs/*.mrc',
               'volume': 'volumes/reference.mrc',
               'relion_it020_data': 'gold/relion_it020_data.star',
               'input_particles': 'gold/input_particles.star',
               'particles': 'gold/particles.sqlite'})


DataSet(name='ribo_movies', folder='ribo_movies', 
        files={'posAllDir': 'pickingXmipp',
               'movies': 'movies/1??_*.mrcs'})


DataSet(name='model',  folder='model', 
        files={'modelGoldSqlite': 'gold/model_gold.sqlite', 
               'modelGoldXml': 'gold/model_gold.xml',
               'classesSelection': 'gold/classes_selection.sqlite'})

DataSet(name='rct',  folder='rct', 
        files={'positions': 'positions',
               'untilted': 'micrographs/F_rct_u*.tif',
               'tilted': 'micrographs/F_rct_t*.tif',
               'classes': 'classes/classes2D_stable_core.sqlite'})

DataSet(name='emx',  folder='emx', 
        files={
               'coordinatesGoldT1': 'coordinates/Test1/coordinates_gold.sqlite',              
               'coordinatesT1': 'coordinates/Test1/coordinates.emx',
               'defocusParticleT2': 'defocusParticle/particles.emx',
               'micrographsGoldT2': 'defocusParticle/micrographs_gold.sqlite',              
               'particlesGoldT2': 'defocusParticle/particles_gold.sqlite',              
              })
               
               
               
               
               
               
