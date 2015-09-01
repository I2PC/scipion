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

ProgCorrectWiener2D::ProgCorrectWiener2D()
{
	rank=0;
	Nprocessors=1;
}


void ProgCorrectWiener2D::produceSideInfo()
{
	FileName fn_img;
	Image<double> img;
    md.read(fn_input);
 	md.getValue(MDL_IMAGE, fn_img, md.firstObject());
    img.read(fn_img);
    img.getDimensions(Xdim,Ydim,Zdim,Ndim);
}

// Read arguments ==========================================================
void ProgCorrectWiener2D::readParams()
{

    // Read command line
    fn_out = getParam("-o");
    fn_input = getParam("-i");
    phase_flipped = checkParam("--phase_flipped");
    pad = XMIPP_MAX(1.,getDoubleParam("--pad"));
    isIsotropic = checkParam("--isIsotropic");;
    wiener_constant  = getDoubleParam("--wc");
}

// Define parameters ==========================================================
void ProgCorrectWiener2D::defineParams()
{
    addUsageLine("Perform CTF correction to 2D projection images with estimated ctfs using a Wiener filter.");
    addKeywords("correct CTF by Wiener filtering");
    addParamsLine("   -i <md_file>            : Metadata file with input projections with CTF information");
    addParamsLine("   -o <md_file>            : Output metadata file with projections with corrected CTF");
    addParamsLine("   [--phase_flipped]       : Is the data already phase-flipped?");
    addParamsLine("   [--isIsotropic]         : Must be considered the defocus isotropic?");
    addParamsLine("   [--wc <float=-1>]       : Wiener-filter constant (if < 0: use FREALIGN default)");
    addParamsLine("   [--pad <factor=2.> ]    : Padding factor for Wiener correction");
}

// Define parameters ==========================================================
void ProgCorrectWiener2D::run()
{
	produceSideInfo();
	MetaData mdOut;

	mdOut = md;

	FileName fn_img;
	Image<double> img;

    if (md.size() == 0)
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE,
                     "Program should never execute this line, something went wrong");

    MultidimArray<double> Mwien;

    // Loop over all images in the metadata
    int imgNum = 1;
    FOR_ALL_OBJECTS_IN_METADATA(mdOut)
    {
    	 mdOut.getValue(MDL_IMAGE, fn_img, __iter.objId);
    	 img.read(fn_img);
    	 img.getDimensions(Xdim,Ydim,Zdim,Ndim);
    	 CTFDescription ctf;
    	 ctf.readFromMetadataRow(mdOut, __iter.objId);

    	 generateWienerFilter(Mwien,ctf);
    	 applyWienerFilter(Mwien,img());

    	 //generate Output
    	 FileName outFileName;
    	 FileName rootName = (fn_out.removeAllExtensions()).addExtension(".stk");
    	 outFileName.compose(imgNum,rootName);
    	 img.write(outFileName);
    	 imgNum++;
    	 mdOut.setValue(MDL_IMAGE, outFileName, __iter.objId);

    }

    mdOut.removeLabel(MDL_CTF_DEFOCUSA);
    mdOut.removeLabel(MDL_CTF_DEFOCUSU);
    mdOut.removeLabel(MDL_CTF_DEFOCUS_ANGLE);
    mdOut.removeLabel(MDL_CTF_DEFOCUSV);
    mdOut.removeLabel(MDL_CTF_BG_BASELINE);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_ANGLE);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_CU);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_CV);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_K);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_SIGMAU);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN2_SIGMAV);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_ANGLE);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_CU);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_CV);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_K);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_SIGMAU);
    mdOut.removeLabel(MDL_CTF_BG_GAUSSIAN_SIGMAV);
    mdOut.removeLabel(MDL_CTF_BG_SQRT_ANGLE);
    mdOut.removeLabel(MDL_CTF_BG_SQRT_K);
    mdOut.removeLabel(MDL_CTF_BG_SQRT_U);
    mdOut.removeLabel(MDL_CTF_BG_SQRT_V);
    mdOut.write("kk.xmd");

}

void ProgCorrectWiener2D::generateWienerFilter(MultidimArray<double> &Mwien,CTFDescription & ctf)
{
	int paddim = Ydim*pad;
	ctf.enable_CTF = true;
	ctf.enable_CTFnoise = false;
	ctf.produceSideInfo();

	MultidimArray<std::complex<double> > ctfComplex;
	MultidimArray<double> ctfIm;

	Mwien.resize(paddim,paddim);

	//NO entiendo esto:
	//ctf.Tm /= sqrt(2.);

	if (isIsotropic)
	{
		double avgdef = (ctf.DeltafU + ctf.DeltafV)/2.;
		ctf.DeltafU = avgdef;
		ctf.DeltafV = avgdef;
	}


	ctfIm.resize(1, 1, paddim, paddim,false);
	//Esto puede estar mal. Cuidado con el sampling de la ctf!!!
	ctf.generateCTF(paddim, paddim, ctfComplex);

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ctfIm)
	{
		if (phase_flipped)
			dAij(ctfIm, i, j) = fabs(dAij(ctfComplex, i, j).real());
		else
			dAij(ctfIm, i, j) = dAij(ctfComplex, i, j).real();
	}

#ifdef DEBUG

	{
		Image<double> save;
		save()=ctfIm;
		save.write("ctf.spi");
		std::cout << "Press any key\n";
		char c;
		std::cin >> c;
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

void ProgCorrectWiener2D::applyWienerFilter(MultidimArray<double> &Mwien, MultidimArray<double> &img)
{
	img.setXmippOrigin();
	int paddim = Ydim*pad;
	MultidimArray<std::complex<double> > Faux;


#ifdef DEBUG

{
	Image<double> save;
	save()=img;
	save.write("img.spi");
}
#endif
#undef DEBUG

    if (paddim > Xdim)
    {
        // pad real-space image
        int x0 = FIRST_XMIPP_INDEX(paddim);
        int xF = LAST_XMIPP_INDEX(paddim);
        img.selfWindow(x0, x0, xF, xF);
    }

    FourierTransform(img, Faux);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mwien)
    {
        dAij(Faux,i,j) *= dAij(Mwien,i,j);
    }

    InverseFourierTransform(Faux, img);
    if (paddim > Xdim)
    {
        // de-pad real-space image
        int x0 = FIRST_XMIPP_INDEX(Xdim);
        int xF = LAST_XMIPP_INDEX(Xdim);
        img.selfWindow(x0, x0, xF, xF);
    }

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






