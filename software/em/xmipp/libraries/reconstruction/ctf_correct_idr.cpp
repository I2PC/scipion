/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "ctf_correct_idr.h"

void ProgCtfCorrectIdr::readParams()
{
    XmippMetadataProgram::readParams();
    fn_vol = getParam("--vol");
    mu = getDoubleParam("--mu");
}

void ProgCtfCorrectIdr::preProcess()
{
    V.read(fn_vol);
    V().setXmippOrigin();
}

void ProgCtfCorrectIdr::show()
{
    if (!verbose)
        return;
    XmippMetadataProgram::show();
    std::cout << "Input volume: " << fn_vol << std::endl
    << "Relaxation factor: " << mu << std::endl;
}

void ProgCtfCorrectIdr::defineParams()
{
    each_image_produces_an_output = true;
    addUsageLine("Correct CTF by using IDR");
    addUsageLine("+This utility allows you to make reconstructions removing the effect ");
    addUsageLine("+of the CTF. The projection images are modified according to a current ");
    addUsageLine("+guess of the reconstruction and the CTF definition. The modified images ");
    addUsageLine("+should be used for reconstruction which now serves as a new current guess ");
    addUsageLine("+for a new IDR iteration.");
    addUsageLine("+The method is fully described at http://www.ncbi.nlm.nih.gov/pubmed/15005161");
    addUsageLine("+");
    addUsageLine("+The program will assume that the image phase have already been corrected");
    addUsageLine("+(otherwise the initial solution is too far from the true solution and the ");
    addUsageLine("+algorithm does not converge).");
    addSeeAlsoLine("ctf_phase_flip");
    defaultComments["-i"].clear();
    defaultComments["-i"].addComment("Metadata with images and corresponding CTFs (ctfdat)");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --vol <volume>  : Volume with the current reconstruction");
    addParamsLine("  [--mu <s=1.8>]   : Relaxation factor");
    addExampleLine("xmipp_ctf_correct_idr -i images.sel --vol currentGuess.vol -o idrCorrectedImages.stk");
}

/* IDR correction ---------------------------------------------------------- */
//#define DEBUG
void ProgCtfCorrectIdr::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    FileName fn_ctf;
    rowIn.getValue(MDL_CTF_MODEL,fn_ctf);

    // Read current input image
    Ireal.readApplyGeo(fnImg, rowIn);
    int Ydim = YSIZE(Ireal());
    int Xdim = XSIZE(Ireal());

    // Project the volume in the same direction
    projectVolume(V(), Itheo, Ydim, Xdim, Ireal.rot(), Ireal.tilt(), Ireal.psi());

    // Copy to theo_CTF and resize
    Itheo_CTF = Itheo();
    Itheo_CTF.setXmippOrigin();
    Itheo_CTF.selfWindow(FIRST_XMIPP_INDEX(2*Ydim), FIRST_XMIPP_INDEX(2*Xdim),
                         LAST_XMIPP_INDEX(2*Ydim), LAST_XMIPP_INDEX(2*Xdim));

    // Read CTF file
    if (last_fn_ctf!=fn_ctf)
    {
        ctf.FilterBand = CTF;
        ctf.ctf.read(fn_ctf);
        ctf.ctf.enable_CTFnoise = false;
        ctf.ctf.produceSideInfo();
        ctf.generateMask(Itheo_CTF);
        ctf.correctPhase();
        last_fn_ctf=fn_ctf;
    }

    // Apply CTF
    ctf.applyMaskSpace(Itheo_CTF);
    Itheo_CTF.selfWindow(FIRST_XMIPP_INDEX(Ydim), FIRST_XMIPP_INDEX(Xdim),
                         LAST_XMIPP_INDEX(Ydim), LAST_XMIPP_INDEX(Xdim));

#ifdef DEBUG

    Itheo.write("PPPtheo.xmp");
    Ireal.write("PPPreal.xmp");
#endif

    // Apply IDR process
    MultidimArray<double> &mItheo=Itheo();
    const MultidimArray<double> &mIreal=Ireal();
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mItheo)
    DIRECT_MULTIDIM_ELEM(mItheo, n) = mu * DIRECT_MULTIDIM_ELEM(mIreal, n) +
                                      (DIRECT_MULTIDIM_ELEM(mItheo, n) -
                                       mu * DIRECT_MULTIDIM_ELEM(Itheo_CTF, n));

    // Save output image
    Itheo.write(fnImgOut);

#ifdef DEBUG

    Image<double> save;
    save()=Itheo_CTF;
    save.write("PPPtheo_CTF.xmp");
    save() = Itheo() - mu * Itheo_CTF;
    save.write("PPPdiff.xmp");
    Itheo.write("PPPidr.xmp");
    std::cout << "Press any key to continue\n";
    char c;
    std::cin >> c;
#endif
}
