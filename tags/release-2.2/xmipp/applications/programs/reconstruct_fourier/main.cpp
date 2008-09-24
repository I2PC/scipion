/***************************************************************************
 *
 * Authors:     Roberto Marabini roberto@cnb.uam.es
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
/*
read command line sel file, symmetry, output vol, pad factor image, pad factor
volume use_up_resolution_ams sampling


loop thorough all images
    pad image by factor 
    Fourier transform
    create Euler matrix
    loop over pixels in fourier transform
        multiply point position for euler matrix
        use kaiser besel interpolator
            m=2, radius=2.7,alpha=6.8
        remember simetries
        volume and projection at different scale
        you will need an extra volume for normalization
normalize by sum of (1-distance)
apply symmetry again   

what happens when no point is assigned. leave as zero
If point Z negative compute conjugate

when compute resolution do not compute extra volumes but
save data in vector
*/
#include <reconstruction/reconstruct_fourier.h>

int main(int argc, char **argv)
{

    VolumeXmipp vol;
    Prog_RecFourier_prm prm;

    try
    {

        prm.read(argc, argv);
        prm.show();
        prm.produce_Side_info();
        prm.MainLoop(vol);
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }


}

