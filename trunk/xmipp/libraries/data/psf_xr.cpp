/***************************************************************************
 *
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#include "psf_xr.h"
#include "args.h"
#include "fft.h"

/* Read -------------------------------------------------------------------- */
void XmippXRPSF::read(const FileName &fn) {
	FILE *fh_param;
	if ((fh_param = fopen(fn.c_str(), "r")) == NULL)
		REPORT_ERROR(1,
				(std::string)"XmippXROTF::read: There is a problem "
				"opening the file " + fn);

	try {
		lambda = textToFloat(getParameter(fh_param, "lambda", 0, "2.43"))
				* 1e-9;
		Flens
				= textToFloat(getParameter(fh_param, "focal_length", 0,
						"1.4742")) * 1e-3;
		Nzp = textToFloat(getParameter(fh_param, "zones_number", 0, "560"));
		Ms = textToFloat(getParameter(fh_param, "magnification", 0, "2304"));
		dxo = textToFloat(getParameter(fh_param, "sampling_rate", 0, "1")) * 1
				- 9;
		if (checkParameter(fh_param, "z_sampling_rate"))
			dzo = textToFloat(getParameter(fh_param, "z_sampling_rate", 0))
					* 1e-9;
		else
			dzo = dxo;
	} catch (Xmipp_error XE) {
		std::cout << XE << std::endl;
		REPORT_ERROR(1, (std::string)"There is an error reading " + fn);
	}
	fclose(fh_param);
}

/* Write ------------------------------------------------------------------- */
void XmippXRPSF::write(const FileName &fn) {
	std::ofstream fh_param;
	fh_param.open(fn.c_str());
	if (!fh_param)
		REPORT_ERROR(1, (std::string)"Xmipp_CTF::write: Cannot open " + fn +
				" for output");
	fh_param << *this << std::endl;
	fh_param.close();
}

/* Usage ------------------------------------------------------------------- */
void XmippXRPSF::usage() {
	std::cerr << "  [lambda=<lambda=2.43>]            : Wavelength in nm\n"
			<< "  [focal_length=<Flens=1.4742>]     : Focal length in mm\n"
			<< "  [zones_number=<Nzp=560>]          : zone plate number\n"
			<< "  [magnification=<Ms=2304>]         : Microscope magnification\n"
			<< "  [sampling_rate=<dxo=1>]           : Object XY plane resolution in nm/pixel\n"
			<< "  [z_sampling_rate=<dzo=dxo>]       : Object Z axis resolution in nm/pixel\n";
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator <<(std::ostream &out, const XmippXRPSF &psf) {

	    out << "lambda=               " << psf.lambda * 1e9 << " nm" << std::endl
			<< "focal_length=         " << psf.Flens * 1e3 << " mm" << std::endl
			<< "zones_number=         " << psf.Nzp << std::endl
			<< "magnification=        " << psf.Ms << std::endl
			<< "sampling_rate=        " << psf.dxo * 1e9 << " nm" << std::endl
			<< "z_sampling_rate=      " << psf.dzo * 1e9 << " nm" << std::endl
			<< std::endl
			<< "Lens Radius=          " << psf.Rlens * 1e6 << " um" << std::endl
			<< "Zo=                   " << psf.Zo * 1e3 << " mm" << std::endl
			<< "Zi=                   " << psf.Zi * 1e3 << " mm" << std::endl
			<< "Minimum Resolution    " << psf.dxiMax * 1e6 << " um" << std::endl

	;

	return out;
}

/* Default values ---------------------------------------------------------- */
void XmippXRPSF::clear() {
	lambda = 2.43e-9;
	Flens = 1.4742e-3;
	Nzp = 560;
	Ms = 2304;
	dzo = dxo = 1e-9;
}

/* Produce Side Information ------------------------------------------------ */
void XmippXRPSF::produceSideInfo() {
	/// Calculation of the rest of microscope parameters
	Rlens = sqrt(Nzp * lambda * Flens);
	Zo = (1 + 1 / Ms) * Flens;
	Zi = Zo * Ms;
	dxiMax = lambda * Zi / (2 * Rlens);
}

/* Apply the CTF to an image ----------------------------------------------- */
void XmippXRPSF::applyOTF(Matrix2D<double> &I) const {
	//    Matrix1D<int>    idx(2);
	//    Matrix1D<double> freq(2);
	//    FOR_ALL_ELEMENTS_IN_MATRIX2D(FFTI)
	//    {
	//        XX(idx) = j;
	//        YY(idx) = i;
	//        FFT_idx2digfreq(FFTI, idx, freq);
	//        double ctf = CTF_at(XX(freq), YY(freq));
	//        FFTI(i, j) *= ctf;
	//    }
}

/* Generate OTF Image ------------------------------------------------------ */
//#define DEBUG
void XmippXRPSF::generateOTF(int Ydim, int Xdim) const {
	//    Matrix1D<int>    idx(2);
	//    Matrix1D<double> freq(2);
	//    OTF.resize(Ydim, Xdim);
	//#ifdef DEBUG
	//    std::cout << "OTF:\n" << *this << std::endl;
	//#endif
	//    FOR_ALL_ELEMENTS_IN_MATRIX2D(OTF)
	//    {
	//        XX(idx) = j;
	//        YY(idx) = i;
	//        FFT_idx2digfreq(OTF, idx, freq);
	//        digfreq2contfreq(freq, freq, Tm);
	//        OTF(i, j) = OTF_at(XX(freq), YY(freq));
	//#ifdef DEBUG
	//        if (i == 0)
	//            std::cout << i << " " << j << " " << YY(freq) << " " << XX(freq)
	//            << " " << OTF(i, j) << std::endl;
	//#endif
	//    }
}

/* Generate the quadratic phase distribution of an ideal lens ------------- */
Matrix2D<std::complex<double> > lensPD(double Flens, double lambda, double dx,
		double Nx, double Ny) {

	Matrix2D<std::complex<double> > Im(Ny, Nx);

	double Lx0 = Nx * dx, Ly0 = Ny * dx, x, y, phase;

	Im.setXmippOrigin();

	FOR_ALL_ELEMENTS_IN_MATRIX2D(Im)
		{
			/// For indices in standard fashion
//			x = (double) j * dx + (dx - Lx0) / 2;
//			y = (double) i * dx + (dx - Ly0) / 2;

			x = (double) j * dx + (dx) / 2;
			y = (double) i * dx + (dx) / 2;

			phase = (-PI / lambda / Flens) * (x * x + y * y);

			(Im(i, j)).real() = cos(phase);
			(Im(i, j)).imag() = sin(phase);
		}

	return Im;

}

#undef DEBUG

