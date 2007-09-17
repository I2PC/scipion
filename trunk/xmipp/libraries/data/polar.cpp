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
    Matrix1D<complex<double> > Faux;
    out.clear();

    for (int iring = 0; iring < in.rings.size(); iring++)
    { 
	Faux = in.rings[iring];
	if (conjugated)
	    for (int j = 0; j < XSIZE(Faux); j++)
		Faux(j) = conj( Faux(j) );
	int oridim = 2 * (XSIZE(Faux) - 1);
	InverseFourierTransformHalf(Faux,Maux,oridim);
	Maux.setStartingX(0);
	out.rings.push_back(Maux);
    }

    out.mode = in.mode;
    out.ring_radius = in.ring_radius;

}

// conversion for complex polar
void convertPolarToSingleVector(const Polar<complex<double> > & in, 
			   vector<double> & out)
{
    out.clear();
    for (int i = 0; i < in.rings.size(); i++)
	for (int j = 0; j < XSIZE(in.rings[i]); j++)
	{
	    out.push_back((in.rings[i](j)).real());
	    out.push_back((in.rings[i](j)).imag());
	}
}

// conversion for real polar
void convertPolarToSingleVector(const Polar<double > & in, 
			   vector<double> & out)
{
    out.clear();
    for (int i = 0; i < in.rings.size(); i++)
	for (int j = 0; j < XSIZE(in.rings[i]); j++)
	{
	    out.push_back(in.rings[i](j));
	}
}

// conversion for complex polar
void convertSingleVectorToPolar(const vector<double> & in, 
				Polar<complex<double> > & out)
{
    int c = 0;
    for (int i = 0; i < out.rings.size(); i++)
    {
	for (int j = 0; j < XSIZE(out.rings[i]); j++)
	{
	    out.rings[i](j) = complex<double>(in[2*c],in[2*c+1]);
	    c++;
	}
    }

    if (2*c != in.size()) 
	REPORT_ERROR(1,"convertFromSingleVector: incorrect vector size for this template");

}

// conversion for real polar
void convertSingleVectorToPolar(const vector<double> & in, 
				Polar<double> & out)
{
    int c = 0;
    for (int i = 0; i < out.rings.size(); i++)
    {
	for (int j = 0; j < XSIZE(out.rings[i]); j++)
	{
	    out.rings[i](j) = in[c];
	    c++;
	}
    }

    if (c != in.size()) 
	REPORT_ERROR(1,"convertFromSingleVector: incorrect vector size for this template");

}

