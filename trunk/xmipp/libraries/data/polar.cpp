/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.uam.es)
 *
 *  This code is strongly based on ideas by Pawel Penczek & Zhengfan
 *  Yang as implemented in SPARX at the University of Texas - Houston 
 *  Medical School
 *
 *  see P. A. Penczek, R. Renka, and H. Schomberg,
 *      J. Opt. Soc. Am. _21_, 449 (2004)
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
#include "polar.h"

// inverse FFT of real signal

void inverseFourierTransformRings(const Polar<complex<double> > & in, 
				  Polar<double> &out, bool conjugated)
{
    Matrix1D<double> Maux;
    Matrix1D<complex<double> > Faux, Faux2;
    out.clear();

    for (int iring = 0; iring < in.rings.size(); iring++)
    { 
	Faux2 = in.rings[iring];
	if (conjugated)
	    for (int j = 0; j < XSIZE(Faux2); j++)
		Faux2(j) = conj( Faux2(j) );
	int oridim = 2 * (XSIZE(Faux2) - 1);
	InverseFourierTransformHalf(Faux,Maux,oridim);
	Maux.setStartingX(0);
	out.rings.push_back(Maux);
    }

    out.mode = in.mode;
    out.ring_radius = in.ring_radius;

}

