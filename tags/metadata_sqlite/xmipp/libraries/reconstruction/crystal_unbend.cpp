/***************************************************************************
 *
 * Authors:     Debora Gil
                Roberto Marabini (roberto@mipg.upenn.edu)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/geometry.h>
#include <interface/aph.h>

#include <stdio.h>

#include "crystal_unbend.h"

/* Main routine for transforming ------------------------------------------- */
void ROUT_Umbend(ImUmbend &ImUm)
{


    std::cout << "Reading" << std::endl;
    //Read .cor file
    ImUm.ReadMRCCord();
    //(ImUm.myUmbendfile).read(ImUm.FN_Correlation);
    std::cout << "PeaksCorrespondance" << std::endl;
    //Uniform sistematic sampling
    ImUm.PeaksCorresp();
    std::cout << "Crystal UnBending" << std::endl;
    //Displacement Interpolation
    ImUm.UnBending();
}


