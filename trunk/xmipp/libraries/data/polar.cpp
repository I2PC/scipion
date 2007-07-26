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

// Calculate 1D FFTs of all rings
Polar<complex<double> > Polar<double>::fourierTransformRings(bool conjugated)
{
    Polar<complex<double> > out;
    Matrix1D<complex<double> > Fring, Faux;
    out.clear();
    for (int iring = 0; iring < rings.size(); iring++)
    { 
	FourierTransform(rings[iring],Faux);
	// Store only asymmetric half of the Fourier Transform
	Whole2Half(Faux, Fring);
	if (conjugated)
	{
	    for (int i = 0; i < XSIZE(Fring); i++)
		Fring(i)=conj(Fring(i));
	}
	out.rings.push_back(Fring);
    }

    out.mode = mode;
    out.ring_radius = ring_radius;
    return out;
}

// Calculate 1D inverse FFTs of all rings
Polar<double> Polar<complex<double> >::inverseFourierTransformRings()
{
    Polar<double> out;
    Matrix1D<double> Mring;
    Matrix1D<complex<double> > Faux;
    out.clear();
    for (int iring = 0; iring < rings.size(); iring++)
    { 
	// Resize to full-length Fourier Transform
	// Original NSAM is always even!!
	int orixdim = 2 * ((int)ring_radius[iring] - 1);
	Half2Whole(rings[iring],Faux,orixdim);
	InverseFourierTransform(Faux,Mring);
	out.rings.push_back(Mring);
    }
    out.mode = mode;
    out.ring_radius = ring_radius;
    return out;
}

// Calculate real-space rotational cross-correlation via convolution theorem
void rotationalCrossCorrelation(const Polar<complex<double> > &M1,
				const Polar<complex<double> > &M2,
				Matrix1D<double> &angles, Matrix1D<double> &corr)
{
    Matrix1D<complex<double> > Fsum;
    complex<double> aux;

    int nrings = M1.getRingNo();
    if (nrings != M2.getRingNo())
	REPORT_ERROR(1,"getBestAnglePolarGridding: polar structures have unequal number of rings!");
    int nsam_last = M1.getSampleNo(nrings - 1);

    // Multiply M1 and M2 over all rings and sum
    Fsum.resize(nsam_last);
    for (int iring = 0; iring < nrings; iring++)
    { 
	int nsam = M1.getSampleNo(iring);
	for (int i = 0; i < nsam; i++)
	{
	    // Multiply M1 with M2. Assume M2 is already complex conjugated!
	    aux = M1(iring,i) * M2(iring,i);
	    // Apply weights for distinct sizes of the rings
	    Fsum(i) += (double)nsam * aux; 
	}
    }

    // Inverse FFT to get real-space correlations
    InverseFourierTransform(Fsum,corr);
    angles.resize(nsam_last);
    for (int i = 0; i < nsam_last; i++)
	angles(i)=(double)i*360./(nsam_last);

}


// Calculate Fourier-space rotational cross-correlation via convolution theorem
void rotationalCrossCorrelation(const Polar<double> &M1,
				const Polar<double> &M2,
				Matrix1D<double> &angles, Matrix1D<complex<double> > &corr)
{
    Matrix1D<double> Fsum;
    double aux;

    int nrings = M1.getRingNo();
    if (nrings != M2.getRingNo())
	REPORT_ERROR(1,"getBestAnglePolarGridding: polar structures have unequal number of rings!");
    int nsam_last = M1.getSampleNo(nrings - 1);

    // Multiply M1 and M2 over all rings and sum
    Fsum.resize(nsam_last);
    for (int iring = 0; iring < nrings; iring++)
    { 
	int nsam = M1.getSampleNo(iring);
	for (int i = 0; i < nsam; i++)
	{
	    // Multiply M1 with M2.
	    aux = M1(iring,i) * M2(iring,i);
	    // Apply weights for distinct sizes of the rings
	    Fsum(i) += (double)nsam * aux; 
	}
    }

    // Inverse FFT to get real-space correlations
    FourierTransform(Fsum,corr);
    angles.resize(nsam_last);
    for (int i = 0; i < nsam_last; i++)
	angles(i)=(double)i*360./nsam_last;

}
