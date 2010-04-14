/***************************************************************************
 *
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
 *
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

#include <data/psf_xr.h>
#include <data/args.h>


void usage();

int main(int argc, char **argv) {
	FileName fnPSF, fnPSFOut;
	XmippXRPSF psf;

	try {
		if (checkParameter(argc, argv, "-psf")) {
			fnPSF = getParameter(argc, argv, "-psf");
			psf.read(fnPSF);
		} else
			psf.clear();

	} catch (Xmipp_error &XE) {
		std::cerr << XE << std::endl;
		//    	psf.usage();
		usage();
		return 1;
	}

	try {
		//    	XmippXRPSF psf;
		psf.produceSideInfo();
		std::cout << psf << std::endl;

		if (checkParameter(argc, argv, "-psfout")) {
			fnPSFOut = getParameter(argc, argv, "-psfout");
//			psf.write(fnPSFOut);


		Matrix2D < std::complex < double > > Ic(128,128);
//		Matrix2D < double > Ic(128,128);


//		Ic = (lensPD(psf.Flens, psf.lambda,psf.dxo, 128, 128));
		lensPD(Ic, psf.Flens, psf.lambda,psf.dxo,psf.dxo);

#define DEBUG

		psf.generateOTF(Ic);

#ifdef DEBUG
		ImageXmipp Im;
		Im().resize(Ic);
		FOR_ALL_ELEMENTS_IN_MATRIX2D(Ic)
		         Im(i,j) = arg(Ic(i,j));

		Im.write(fnPSFOut);
#endif
//		std::cout << Ic << std::endl;
//		std::cout << Ir << std::endl;

		}

	} catch (Xmipp_error &XE) {
		std::cerr << XE << std::endl;
		return 1;
	}
	return 0;
}

void usage() {
	std::cerr << "Usage: project_xr [options]\n"
			<< "   -psf <PSF description file>      : PSF characteristic of the microscope \n"
			<< "  [-w_dir \"[X=1,Y=0]\"             : test \n"
			<< "  [-w_step <step=0.001>]          : test\n";
}

#undef DEBUG

