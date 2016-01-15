/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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

#include "ctf_correct_wiener2d.h"

// Read arguments ==========================================================
void ProgCorrectWiener2D::readParams()
{
    XmippMetadataProgram::readParams();
    phase_flipped = checkParam("--phase_flipped");
    pad = XMIPP_MAX(1.,getDoubleParam("--pad"));
    isIsotropic = checkParam("--isIsotropic");
    wiener_constant  = getDoubleParam("--wc");
    correct_envelope = checkParam("--correct_envelope");
    sampling_rate = getDoubleParam("--sampling_rate");
}

// Define parameters ==========================================================
void ProgCorrectWiener2D::defineParams()
{
    addUsageLine("Perform CTF correction to 2D projection images with estimated ctfs using a Wiener filter.");
    each_image_produces_an_output = true;
    XmippMetadataProgram::defineParams();
    addKeywords("correct CTF by Wiener filtering");
    addParamsLine("   [--phase_flipped]       : Is the data already phase-flipped?");
    addParamsLine("   [--isIsotropic]         : Must be considered the defocus isotropic?");
    addParamsLine("   [--sampling_rate <float=1.0>]     : Sampling rate of the input particles");
    addParamsLine("   [--wc <float=-1>]       : Wiener-filter constant (if < 0: use FREALIGN default)");
    addParamsLine("   [--pad <factor=2.> ]    : Padding factor for Wiener correction");
    addParamsLine("   [--correct_envelope]     : Correct the CTF envelope");
}

// Define parameters ==========================================================
void ProgCorrectWiener2D::postProcess()
{

	MetaData &ptrMdOut=*getOutputMd();

	ptrMdOut.removeLabel(MDL_CTF_DEFOCUSA);
	ptrMdOut.removeLabel(MDL_CTF_DEFOCUSU);
	ptrMdOut.removeLabel(MDL_CTF_DEFOCUS_ANGLE);
	ptrMdOut.removeLabel(MDL_CTF_DEFOCUSV);
	ptrMdOut.removeLabel(MDL_CTF_BG_BASELINE);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_ANGLE);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_CU);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_CV);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_K);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_SIGMAU);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_SIGMAV);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_ANGLE);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_CU);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_CV);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_K);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_SIGMAU);
	ptrMdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_SIGMAV);
	ptrMdOut.removeLabel(MDL_CTF_BG_SQRT_ANGLE);
	ptrMdOut.removeLabel(MDL_CTF_BG_SQRT_K);
	ptrMdOut.removeLabel(MDL_CTF_BG_SQRT_U);
	ptrMdOut.removeLabel(MDL_CTF_BG_SQRT_V);
	ptrMdOut.removeLabel(MDL_CTF_CA);
	ptrMdOut.removeLabel(MDL_CTF_CONVERGENCE_CONE);
	ptrMdOut.removeLabel(MDL_CTF_ENERGY_LOSS);
	ptrMdOut.removeLabel(MDL_CTF_ENVELOPE);
	ptrMdOut.removeLabel(MDL_CTF_LENS_STABILITY);
	ptrMdOut.removeLabel(MDL_CTF_TRANSVERSAL_DISPLACEMENT);
	ptrMdOut.removeLabel(MDL_CTF_LONGITUDINAL_DISPLACEMENT);
	ptrMdOut.removeLabel(MDL_CTF_K);

	ptrMdOut.write(fn_out.replaceExtension("xmd"));

}

void ProgCorrectWiener2D::generateWienerFilter(MultidimArray<double> &Mwien, CTFDescription & ctf)
{
	int paddim = Ydim*pad;
	ctf.enable_CTF = true;
	ctf.enable_CTFnoise = false;
	ctf.produceSideInfo();

	MultidimArray<std::complex<double> > ctfComplex;
	MultidimArray<double> ctfIm;

	Mwien.resize(paddim,paddim);

	ctf.Tm = sampling_rate;
	//ctf.Tm /= pad;

	if (isIsotropic)
	{
		double avgdef = (ctf.DeltafU + ctf.DeltafV)/2.;
		ctf.DeltafU = avgdef;
		ctf.DeltafV = avgdef;
	}


	ctfIm.resize(1, 1, paddim, paddim,false);
	//Esto puede estar mal. Cuidado con el sampling de la ctf!!!

	if (correct_envelope)
		ctf.generateCTF(paddim, paddim, ctfComplex);
	else
		ctf.generateCTFWithoutDamping(paddim, paddim, ctfComplex);

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ctfIm)
	{
		if (phase_flipped)
			dAij(ctfIm, i, j) = fabs(dAij(ctfComplex, i, j).real());
		else
			dAij(ctfIm, i, j) = dAij(ctfComplex, i, j).real();
	}

//#define DEBUG
#ifdef DEBUG
	{
		Image<double> save;
		save()=ctfIm;
		save.write("PPPctf2.spi");
		//exit(1);
		//std::cout << "Press any key\n";
		//char c;
		//std::cin >> c;
	}
#endif
#undef DEBUG

	double result;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mwien)
	{
		result = (DIRECT_N_YX_ELEM (ctfIm, 0, 0, i, j  ));
		dAij(Mwien,i,j) = (result *result);
	}

#ifdef DEBUG

{
	Image<double> save;
	save()=Mwien;
	save.write("vienerB.spi");
}
#endif
#undef DEBUG

// Add Wiener constant
if (wiener_constant < 0.)
{

	// Use Grigorieff's default for Wiener filter constant: 10% of average over all Mwien terms
	// Grigorieff JSB 157(1) (2006), pp 117-125
	double valueW = 0.1*Mwien.computeAvg();
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mwien)
	{
		dAij(Mwien,i,j) += valueW;
		dAij(Mwien,i,j) = dAij(ctfIm, i, j)/dAij(Mwien, i, j);
	}
}
else
{
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mwien)
    {
		dAij(Mwien,i,j) += wiener_constant;
		dAij(Mwien,i,j) = dAij(ctfIm, i, j)/dAij(Mwien, i, j);
    }
}

#ifdef DEBUG

{
	Image<double> save;
	save()=Mwien;
	save.write("vienerC.spi");
}
#endif
#undef DEBUG

}

//#define DEBUG
void ProgCorrectWiener2D::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{

	rowOut = rowIn;

	img.read(fnImg);
	ctf.readFromMdRow(rowIn);
	img().setXmippOrigin();
	Ydim = YSIZE(img());
	Xdim = YSIZE(img());
	int paddim = Ydim*pad;

	generateWienerFilter(Mwien,ctf);

#ifdef DEBUG

{
	std::cout << ctf << std::endl;

	Image<double> save;
	save()=Mwien;
	save.write("PPPMwien.spi");
	std::cout << "Press any key"<< std::endl;
	char c; std::cin >> c;
}
#endif
#undef DEBUG

    if (paddim > Xdim)
    {
        // pad real-space image
        int x0 = FIRST_XMIPP_INDEX(paddim);
        int xF = LAST_XMIPP_INDEX(paddim);
        img().selfWindow(x0, x0, xF, xF);
    }

    transformer.FourierTransform(img(), Faux);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Faux)
    {
        dAij(Faux,i,j) *= dAij(Mwien,i,j);
    }

    transformer.inverseFourierTransform(Faux, img());
    if (paddim > Xdim)
    {
        // de-pad real-space image
        int x0 = FIRST_XMIPP_INDEX(Xdim);
        int xF = LAST_XMIPP_INDEX(Xdim);
        img().selfWindow(x0, x0, xF, xF);
    }

    img.write(fnImgOut);
    rowOut.setValue(MDL_IMAGE, fnImgOut);

#ifdef DEBUG
{
	Image<double> save;
	save()=img;
	save.write("imgW.spi");
	exit(0);
}
#endif
#undef DEBUG

}
