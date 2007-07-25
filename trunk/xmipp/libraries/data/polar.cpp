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


vector<int> getPolarRingInfo2D(int first_ring, int last_ring, 
			       int step, int mode, bool const_sam)
{
    vector<int> out;
    out.clear();
    int nsam, totsam = 0;
    double twopi;

    if (mode == FULL_CIRCLES)
	twopi = 2. * PI;
    else if (mode == HALF_CIRCLES)
	twopi = PI;
    else
	REPORT_ERROR(1,"Incorrect mode for getPolarSampling2D!");

    for (int iring = first_ring; iring <= last_ring; iring++)
    {
	if (!const_sam)
	    // Oversample each ring twice
	    nsam = MAX(1, (int)twopi * iring * 2);
	else
	    nsam = (int)twopi * last_ring * 2;
	totsam += nsam;
	out.push_back(iring);
	out.push_back(totsam);
	out.push_back(nsam);
    }

    return out;
}

vector<double> getPolarWeights2D(vector<int> ring_info, int mode) 
{
    vector<double> out;
    out.clear();
    double w,twopi;
    
    if (mode == FULL_CIRCLES)
	twopi = 2. * PI;
    else if (mode == HALF_CIRCLES)
	twopi = PI;
    else
	REPORT_ERROR(1,"Incorrect mode for getPolarWeights2D!");

    int nring = ring_info.size() / 3;
    double maxsam = (double) ring_info[3 * nring + 2];
    for (int iring = 0; iring < nring; iring++)
    {
	// copied from sparx
	w = twopi * ring_info[3 * iring] / (double) ring_info[3 * iring + 2];
	w *= maxsam / (double) ring_info[3 * iring + 2];
	out.push_back(w);

	// OR:Weight is reciprocal of number of sampling points???;
	//w = 1. / ring_info[3 * iring + 2];
	//out.push_back(w);
    }

    return out;
}



// Calculate 1D FFTs of all rings
vector<Matrix1D<complex<double> > > fourierTransformRings(
    const vector<Matrix1D<double> > &in, bool conjugated)
{
    vector<Matrix1D<complex<double> > > out;
    Matrix1D<complex<double> > Fring;
    out.clear();

    for (int iring = 0; iring < in.size(); iring++)
    { 
	FourierTransform(in[iring],Fring);
	if (conjugated)
	{
	    for (int i = 0; i < XSIZE(Fring); i++)
		Fring(i)=conj(Fring(i));
	}
	out.push_back(Fring);
    }

    return out;
}

// Calculate 1D inverse FFTs of all rings
vector<Matrix1D<double> > inverseFourierTransformRings(
    const vector<Matrix1D<complex<double> > > &in)
{
    vector<Matrix1D<double> > out;
    Matrix1D<double > Mring;
    out.clear();

    for (int iring = 0; iring < in.size(); iring++)
    { 
	InverseFourierTransform(in[iring],Mring);
	out.push_back(Mring);
    }

    return out;
}

// Rotationally align 2 images (from FTs of their polar rings)
double getBestAnglePolarGridding(const vector<Matrix1D<complex<double> > > &M1,
				 const vector<Matrix1D<complex<double> > > &M2,
				 double &maxcorr, bool interpolate)
{

    vector<Matrix1D<complex<double> > > Fcorr;
    vector<Matrix1D<double> > Mcorr;
    Matrix1D<complex<double> > Fring;
    Matrix1D<double> corr;

    int nrings = M1.size();
    int nrings2 = M2.size();
    if (nrings != nrings2)
	REPORT_ERROR(1,"getBestAnglePolarGridding: polar structures have unequal number of rings!");

    int nsam_last = XSIZE(M1[nrings-1]);
    Fcorr.clear();
    
    // Multiply M1 and M2 over all rings
    for (int iring = 0; iring < nrings; iring++)
    { 
	int nsam = XSIZE(M1[iring]);
	Fring.resize(nsam_last);
	for (int i = 0; i < nsam; i++)
	    Fring(i) = M1[iring](i) * M2[iring](i);
	for (int i = nsam; i < nsam_last; i++)
	    Fring(i) = 0.;
	Fcorr.push_back(Fring);
    }
    
    // InverseFourierTransform each colum of Mcorr
    Mcorr=inverseFourierTransformRings(Fcorr);

    corr.resize(nsam_last);
    int iopt = 0;
    maxcorr = -99.e99;
    for (int i = 0; i < nsam_last; i++)
    {
	for (int iring = 0; iring < nrings; iring++)
	{ 
	    corr(i) += Mcorr[iring](i);
	}
	if (corr(i) > maxcorr)
	{
	    iopt = i;
	    maxcorr = corr(i);
	}
    }

    if (interpolate)
	REPORT_ERROR(1,"getBestAnglePolarGridding: interpolation not yet implemented.");

    return (double)(360. * iopt / nsam_last);
    
}
