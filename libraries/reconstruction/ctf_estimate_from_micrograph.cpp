/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.csic.es)
 *              Carlos Oscar Sanchez Sorzano
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

#include "ctf_estimate_from_micrograph.h"
#include "ctf_enhance_psd.h"

#include <data/args.h>
#include <data/micrograph.h>
#include <data/metadata.h>
#include <data/xmipp_image.h>
#include <data/xmipp_fft.h>
#include <data/xmipp_threads.h>
#include <data/basic_pca.h>
#include <data/normalize.h>

/* Read parameters ========================================================= */
void ProgCTFEstimateFromMicrograph::readParams()
{
    fn_micrograph = getParam("--micrograph");
    fn_root = getParam("--oroot");
    if (fn_root == "")
        fn_root = fn_micrograph.withoutExtension();
    pieceDim = getIntParam("--pieceDim");
    skipBorders = getIntParam("--skipBorders");
    overlap = getDoubleParam("--overlap");
    String aux = getParam("--psd_estimator");
    if (aux == "periodogram")
        PSDEstimator_mode = Periodogram;
    else
    {
        PSDEstimator_mode = ARMA;
        ARMA_prm.readParams(this);
    }
    Nsubpiece = getIntParam("--Nsubpiece");

    String mode = getParam("--mode");
    if (mode == "micrograph")
        psd_mode = OnePerMicrograph;
    else if (mode == "regions")
    {
        psd_mode = OnePerRegion;
        fn_pos = getParam("--mode", 1);
    }
    else if (mode == "particles")
    {
        psd_mode = OnePerParticle;
        fn_pos = getParam("--mode", 1);
    }
    estimate_ctf = !checkParam("--dont_estimate_ctf");
    if (estimate_ctf)
        prmEstimateCTFFromPSD.readBasicParams(this);
    bootstrapN = getIntParam("--bootstrapFit");
}

void ProgCTFEstimateFromMicrograph::defineParams()
{
    addUsageLine("Estimate the CTF from a micrograph.");
    addUsageLine("The PSD of the micrograph is first estimated using periodogram averaging or ");
    addUsageLine("ARMA models ([[http://www.ncbi.nlm.nih.gov/pubmed/12623169][See article]]). ");
    addUsageLine("Then, the PSD is enhanced ([[http://www.ncbi.nlm.nih.gov/pubmed/16987671][See article]]). ");
    addUsageLine("And finally, the CTF is fitted to the PSD, being guided by the enhanced PSD ");
    addUsageLine("([[http://www.ncbi.nlm.nih.gov/pubmed/17911028][See article]]).");
    addParamsLine("   --micrograph <file>         : File with the micrograph");
    addParamsLine("  [--oroot <rootname=\"\">]    : Rootname for output");
    addParamsLine("                               : If not given, the micrograph without extensions is taken");
    addParamsLine("                               :++ rootname.psd or .psdstk contains the PSD or PSDs");
    addParamsLine("==+ PSD estimation");
    addParamsLine("  [--psd_estimator <method=periodogram>] : Method for estimating the PSD");
    addParamsLine("         where <method>");
    addParamsLine("                  periodogram");
    addParamsLine("                  ARMA");
    addParamsLine("  [--pieceDim <d=512>]       : Size of the piece");
    addParamsLine("  [--overlap <o=0.5>]        : Overlap (0=no overlap, 1=full overlap)");
    addParamsLine("  [--skipBorders <s=2>]      : Number of pieces around the border to skip");
    addParamsLine("  [--Nsubpiece <N=1>]        : Each piece is further subdivided into NxN subpieces.");
    addParamsLine("                              : This option is useful for small micrographs in which ");
    addParamsLine("                              : not many pieces of size pieceDim x pieceDim can be defined. ");
    addParamsLine("                              :++ Note that this is not the same as defining a smaller pieceDim. ");
    addParamsLine("                              :++ Defining a smaller pieceDim, would result in a small PSD, while ");
    addParamsLine("                              :++ subdividing the piece results in a large PSD, although smoother.");
    addParamsLine("  [--mode <mode=micrograph>]  : How many PSDs are to be estimated");
    addParamsLine("         where <mode>");
    addParamsLine("                  micrograph  : Single PSD for the whole micrograph");
    addParamsLine("                  regions <file=\"\"> : The micrograph is divided into a region grid ");
    addParamsLine("                              : and a PSD is computed for each one.");
    addParamsLine("                              : The file is metadata with the position of each particle within the micrograph");
    addParamsLine("                  particles <file> : One PSD per particle.");
    addParamsLine("                              : The file is metadata with the position of each particle within the micrograph");
    addParamsLine("==+ CTF fit");
    addParamsLine("  [--dont_estimate_ctf]       : Do not fit a CTF to PSDs");
    ARMA_parameters::defineParams(this);
    ProgCTFEstimateFromPSD::defineBasicParams(this);
    addExampleLine("Estimate PSD", false);
    addExampleLine("xmipp_ctf_estimate_from_micrograph --micrograph micrograph.mrc --dont_estimate_ctf");
    addExampleLine("Estimate a single CTF for the whole micrograph", false);
    addExampleLine("xmipp_ctf_estimate_from_micrograph --micrograph micrograph.mrc --sampling_rate 1.4 --voltage 200 --spherical_aberration 2.5");
    addExampleLine("Estimate a single CTF for the whole micrograph providing a starting point for the defocus",false);
    addExampleLine("xmipp_ctf_estimate_from_micrograph --micrograph micrograph.mrc --sampling_rate 1.4 --voltage 200 --spherical_aberration 2.5 --defocusU -15000");
    addExampleLine("Estimate a CTF per region", false);
    addExampleLine("xmipp_ctf_estimate_from_micrograph --micrograph micrograph.mrc --mode regions micrograph.pos --sampling_rate 1.4 --voltage 200 --spherical_aberration 2.5 --defocusU -15000");
    addExampleLine("Estimate a CTF per particle", false);
    addExampleLine("xmipp_ctf_estimate_from_micrograph --micrograph micrograph.mrc --mode particles micrograph.pos --sampling_rate 1.4 --voltage 200 --spherical_aberration 2.5 --defocusU -15000");
}

/* Construct piece smoother =============================================== */
void constructPieceSmoother(const MultidimArray<double> &piece,
                            MultidimArray<double> &pieceSmoother)
{
    // Attenuate borders to avoid discontinuities
    pieceSmoother.resizeNoCopy(piece);
    pieceSmoother.initConstant(1);
    pieceSmoother.setXmippOrigin();
    double iHalfsize = 2.0 / YSIZE(pieceSmoother);
    const double alpha = 0.025;
    const double alpha1 = 1 - alpha;
    const double ialpha = 1.0 / alpha;
    for (int i = STARTINGY(pieceSmoother); i <= FINISHINGY(pieceSmoother);
         i++)
    {
        double iFraction = fabs(i * iHalfsize);
        if (iFraction > alpha1)
        {
            double maskValue = 0.5
                               * (1 + cos(PI * ((iFraction - 1) * ialpha + 1)));
            for (int j = STARTINGX(pieceSmoother);
                 j <= FINISHINGX(pieceSmoother); j++)
                A2D_ELEM(pieceSmoother,i,j)*=maskValue;
        }
    }

    for (int j = STARTINGX(pieceSmoother); j <= FINISHINGX(pieceSmoother);
         j++)
    {
        double jFraction = fabs(j * iHalfsize);
        if (jFraction > alpha1)
        {
            double maskValue = 0.5
                               * (1 + cos(PI * ((jFraction - 1) * ialpha + 1)));
            for (int i = STARTINGY(pieceSmoother);
                 i <= FINISHINGY(pieceSmoother); i++)
                A2D_ELEM(pieceSmoother,i,j)*=maskValue;
        }
    }
}

/* Compute PSD by piece averaging ========================================== */
//#define DEBUG
void ProgCTFEstimateFromMicrograph::PSD_piece_by_averaging(
    MultidimArray<double> &piece, MultidimArray<double> &psd)
{
    int small_Ydim = 2 * YSIZE(piece) / Nsubpiece;
    int small_Xdim = 2 * XSIZE(piece) / Nsubpiece;
    MultidimArray<double> small_piece(small_Ydim, small_Xdim);

    int Xstep = (XSIZE(piece) - small_Xdim) / (Nsubpiece - 1);
    int Ystep = (YSIZE(piece) - small_Ydim) / (Nsubpiece - 1);
    psd.initZeros(small_piece);
#ifdef DEBUG

    Image<double> save;
    save()=piece;
    save.write("PPPpiece.xmp");
#endif

    MultidimArray<std::complex<double> > Periodogram;
    MultidimArray<double> small_psd;
    MultidimArray<int> pieceMask;
    pieceMask.resizeNoCopy(piece);
    pieceMask.initConstant(1);

    // Attenuate borders to avoid discontinuities
    MultidimArray<double> pieceSmoother;
    constructPieceSmoother(piece, pieceSmoother);

    for (int ii = 0; ii < Nsubpiece; ii++)
        for (int jj = 0; jj < Nsubpiece; jj++)
        {
            // Take the corresponding small piece from the piece
            int i0 = ii * Xstep;
            int j0 = jj * Ystep;

            int i, j, ib, jb;
            for (i = 0, ib = i0; i < small_Ydim; i++, ib++)
                for (j = 0, jb = j0; j < small_Xdim; j++, jb++)
                    DIRECT_A2D_ELEM(small_piece, i, j)=
                        DIRECT_A2D_ELEM(piece, ib, jb);
            normalize_ramp(piece, pieceMask);
            piece *= pieceSmoother;

#ifdef DEBUG

            save()=small_piece;
            save.write("PPPsmall_piece.xmp");
#endif

            // Compute the PSD of the small piece
            small_psd.initZeros(small_piece);
            if (PSDEstimator_mode == ARMA)
            {
                CausalARMA(small_piece, ARMA_prm);
                ARMAFilter(small_piece, small_psd, ARMA_prm);
            }
            else
            {
                FourierTransform(small_piece, Periodogram);
                FFT_magnitude(Periodogram, small_psd);
                small_psd *= small_psd;
                small_psd *= small_Ydim * small_Xdim;
            }

#ifdef DEBUG
            save()=small_psd;
            save.write("PPPsmall_psd.xmp");
#endif

            // Add to the average
            psd += small_psd;
        }

    // Compute the average of all the small pieces and enlarge
    psd *= 1.0 / (Nsubpiece * Nsubpiece);

#ifdef DEBUG

    save()=psd;
    save.write("PPPpsd1.xmp");
#endif

    CenterFFT(psd, true);
    selfScaleToSize(BSPLINE3, psd, YSIZE(piece), XSIZE(piece));
    CenterFFT(psd, false);
    psd.threshold("below", 0, 0);

#ifdef DEBUG

    save()=psd;
    save.write("PPPpsd2.xmp");
    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
#endif
}
#undef DEBUG

/* Main ==================================================================== */
//#define DEBUG
void ProgCTFEstimateFromMicrograph::run()
{
    // Open input files -----------------------------------------------------
    // Open coordinates
    MetaData posFile;
    if (fn_pos != "")
        posFile.read(fn_pos);
    MDIterator iterPosFile(posFile);

    // Open the micrograph --------------------------------------------------
    ImageGeneric M_in;
    size_t Zdim, Ydim, Xdim; // Micrograph dimensions
    M_in.read(fn_micrograph);
    M_in.getDimensions(Xdim, Ydim, Zdim);

    // Compute the number of divisions --------------------------------------
    int div_Number = 0;
    int div_NumberX, div_NumberY;
    if (psd_mode == OnePerParticle)
        div_Number = posFile.size();
    else if (psd_mode == OnePerMicrograph)
    {
        div_NumberX = CEIL((double)Xdim / (pieceDim *(1-overlap))) - 1;
        div_NumberY = CEIL((double)Ydim / (pieceDim *(1-overlap))) - 1;
        div_Number = div_NumberX * div_NumberY;
    }
    else if (psd_mode == OnePerRegion)
    {
        div_NumberX = CEIL((double)Xdim / pieceDim);
        div_NumberY = CEIL((double)Ydim / pieceDim);
        if (div_NumberX<=2*skipBorders || div_NumberY<=2*skipBorders)
        	REPORT_ERROR(ERR_ARG_INCORRECT,formatString("The micrograph is not big enough to skip %d pieces on each side",skipBorders));
        div_Number = div_NumberX * div_NumberY;
    }

    // Process each piece ---------------------------------------------------
    Image<double> psd_avg, psd_std, psd, psd2;
    MultidimArray<std::complex<double> > Periodogram;
    MultidimArray<double> piece(pieceDim, pieceDim);
    psd().resizeNoCopy(piece);
    MultidimArray<double> &mpsd = psd();
    MultidimArray<double> &mpsd2 = psd2();
    PCAMahalanobisAnalyzer pcaAnalyzer;
    MultidimArray<int> PCAmask;
    MultidimArray<float> PCAv;
    double pieceDim2 = pieceDim * pieceDim;
    MultidimArray<int> pieceMask;
    pieceMask.resizeNoCopy(piece);
    pieceMask.initConstant(1);

    //Multidimensional data variables to store the defocus obtained locally for plane fitting
    MultidimArray<double> defocusPlanefittingU(div_NumberX-2*skipBorders, div_NumberY-2*skipBorders);
    MultidimArray<double> defocusPlanefittingV(defocusPlanefittingU);
    MultidimArray<double> Xm(defocusPlanefittingU);
    MultidimArray<double> Ym(defocusPlanefittingU);

    // Attenuate borders to avoid discontinuities
    MultidimArray<double> pieceSmoother;
    constructPieceSmoother(piece, pieceSmoother);

    if (verbose)
        std::cerr << "Computing models of each piece ...\n";

    // Prepare these filenames in case they are needed
    FileName fn_psd;
    if (psd_mode == OnePerMicrograph)
        fn_psd = fn_root + ".psd";
    else
        fn_psd = fn_root + ".psdstk";
    if (fileExists(fn_psd))
    	fn_psd.deleteFile();
    if (fileExists(fn_root+".ctfparam"))
    	FileName(fn_root+".ctfparam").deleteFile();

    if (verbose)
        init_progress_bar(div_Number);
    int N = 1; // Index of current piece
    size_t piecei = 0, piecej = 0; // top-left corner of the current piece
    FourierTransformer transformer;
    int actualDiv_Number = 0;
    while (N <= div_Number)
    {
        bool skip = false;
        int blocki, blockj;
        // Compute the top-left corner of the piece ..........................
        if (psd_mode == OnePerParticle)
        {
            // Read position of the particle
            posFile.getValue(MDL_X, piecej, iterPosFile.objId);
            posFile.getValue(MDL_Y, piecei, iterPosFile.objId);

            // j,i are the selfWindow center, we need the top-left corner
            piecej -= (int) (pieceDim / 2);
            piecei -= (int) (pieceDim / 2);
            if (piecei < 0)
            	piecei = 0;
            if (piecej < 0)
            	piecej = 0;
        }
        else
        {
            int step = pieceDim;
            if (psd_mode == OnePerMicrograph)
                step = (int) ((1 - overlap) * step);
            blocki = (N - 1) / div_NumberX;
            blockj = (N - 1) % div_NumberX;
            if (blocki < skipBorders || blockj < skipBorders
                || blocki > (div_NumberY - skipBorders - 1)
                || blockj > (div_NumberX - skipBorders - 1))
                skip = true;
            piecei = blocki * step;
            piecej = blockj * step;
        }

        // test if the full piece is inside the micrograph
        if (piecei + pieceDim > Ydim)
        	piecei = Ydim - pieceDim;
        if (piecej + pieceDim > Xdim)
        	piecej = Xdim - pieceDim;

        if (!skip)
        {
            // Extract micrograph piece ..........................................
            M_in().window(piece, 0, 0, piecei, piecej, 0, 0, piecei + YSIZE(piece) - 1,
            		      piecej + XSIZE(piece) - 1);
            piece.statisticsAdjust(0, 1);
            normalize_ramp(piece, pieceMask);
            piece *= pieceSmoother;

            // Estimate the power spectrum .......................................
            if (Nsubpiece == 1)
                if (PSDEstimator_mode == ARMA)
                {
                    CausalARMA(piece, ARMA_prm);
                    ARMAFilter(piece, mpsd, ARMA_prm);
                }
                else
                {
                    transformer.completeFourierTransform(piece, Periodogram);
                    FFT_magnitude(Periodogram, mpsd);
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mpsd)
                    DIRECT_MULTIDIM_ELEM(mpsd,n)*=DIRECT_MULTIDIM_ELEM(mpsd,n)*pieceDim2;
                }
            else
                PSD_piece_by_averaging(piece, mpsd);
            mpsd2.resizeNoCopy(mpsd);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mpsd2)
            {
                double psdval = DIRECT_MULTIDIM_ELEM(mpsd,n);
                DIRECT_MULTIDIM_ELEM(mpsd2,n)=psdval*psdval;
            }

            // Perform averaging if applicable ...................................
            if (psd_mode == OnePerMicrograph)
            {
                actualDiv_Number += 1;
                // Compute average and standard deviation
                if (XSIZE(psd_avg()) != XSIZE(mpsd))
                {
                    psd_avg() = mpsd;
                    psd_std() = psd2();
                }
                else
                {
                    psd_avg() += mpsd;
                    psd_std() += psd2();
                }

                // Keep psd for the PCA
                if (XSIZE(PCAmask) == 0)
                {
                    PCAmask.initZeros(mpsd);
                    Matrix1D<int> idx(2);  // Indexes for Fourier plane
                    Matrix1D<double> freq(2); // Frequencies for Fourier plane
                    size_t PCAdim = 0;
                    FOR_ALL_ELEMENTS_IN_ARRAY2D(PCAmask)
                    {
                        VECTOR_R2(idx, j, i);
                        FFT_idx2digfreq(mpsd, idx, freq);
                        double w = freq.module();
                        if (w > 0.05 && w < 0.4)
                        {
                            A2D_ELEM(PCAmask,i,j)=1;
                            ++PCAdim;
                        }
                    }
                    PCAv.initZeros(PCAdim);
                }

                size_t ii = -1;
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PCAmask)
                if (DIRECT_MULTIDIM_ELEM(PCAmask,n))
                    A1D_ELEM(PCAv,++ii)=(float)DIRECT_MULTIDIM_ELEM(mpsd,n);
                pcaAnalyzer.addVector(PCAv);
            }

            // Compute the theoretical model if not averaging ....................
            if (psd_mode != OnePerMicrograph)
            {
                if (bootstrapN != -1)
                    REPORT_ERROR(ERR_VALUE_INCORRECT,
                                 "Bootstrapping is only available for micrograph averages");

                FileName fn_psd_piece;
                fn_psd_piece.compose(N, fn_psd);
                psd.write(fn_psd_piece);
                if (psd_mode == OnePerParticle)
                    posFile.setValue(MDL_PSD, fn_psd_piece, iterPosFile.objId);
                if (estimate_ctf)
                {
                    // Estimate the CTF parameters of this piece
                    prmEstimateCTFFromPSD.fn_psd = fn_psd_piece;
                    CTFDescription ctfmodel;

                    ctfmodel.isLocalCTF = true;
                    ctfmodel.x0 = piecej;
                    ctfmodel.xF = (piecej + pieceDim-1);
                    ctfmodel.y0 = piecei;
                    ctfmodel.yF = (piecei + pieceDim-1);
                    ROUT_Adjust_CTF(prmEstimateCTFFromPSD, ctfmodel, false);

                    int idxi=blocki-skipBorders;
                    int idxj=blockj-skipBorders;
                    A2D_ELEM(defocusPlanefittingU,idxi,idxj)=ctfmodel.DeltafU;
                    A2D_ELEM(defocusPlanefittingV,idxi,idxj)=ctfmodel.DeltafV;

                    A2D_ELEM(Xm,idxi,idxj)=(piecei+pieceDim/2)*ctfmodel.Tm;
                    A2D_ELEM(Ym,idxi,idxj)=(piecej+pieceDim/2)*ctfmodel.Tm;

                    if (psd_mode == OnePerParticle)
                        posFile.setValue(MDL_CTF_MODEL,
                                         fn_psd_piece.withoutExtension() + ".ctfparam",
                                         iterPosFile.objId);
                }
            }
        }
        // Increment the division counter
        ++N;
        if (verbose)
            progress_bar(N);
        if (psd_mode == OnePerParticle)
            iterPosFile.moveNext();
    }
    if (verbose)
        progress_bar(div_Number);

    // If averaging, compute the CTF model ----------------------------------
    if (psd_mode == OnePerMicrograph)
    {
        // Compute the avg and stddev of the local PSDs
        const MultidimArray<double> &mpsd_std = psd_std();
        const MultidimArray<double> &mpsd_avg = psd_avg();
        double idiv_Number = 1.0 / actualDiv_Number;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mpsd_avg)
        {
            DIRECT_MULTIDIM_ELEM(mpsd_avg,n)*=idiv_Number;
            DIRECT_MULTIDIM_ELEM(mpsd_std,n)*=idiv_Number;
            DIRECT_MULTIDIM_ELEM(mpsd_std,n)-=DIRECT_MULTIDIM_ELEM(mpsd_avg,n)*
                                              DIRECT_MULTIDIM_ELEM(mpsd_avg,n);
            if (DIRECT_MULTIDIM_ELEM(mpsd_std,n)<0)
                DIRECT_MULTIDIM_ELEM(mpsd_std,n)=0;
            else
                DIRECT_MULTIDIM_ELEM(mpsd_std,n)=sqrt(DIRECT_MULTIDIM_ELEM(mpsd_std,n));
        }
        psd_avg.write(fn_psd);

        if (estimate_ctf)
        {
            // Estimate the CTF parameters
            std::cerr << "Adjusting CTF model to the PSD ...\n";
            prmEstimateCTFFromPSD.fn_psd = fn_psd;
            CTFDescription ctfmodel;
            if (bootstrapN == -1)
            {
                try {
					// No bootstrapping
					// Compute the PCA of the local PSDs
					pcaAnalyzer.standardarizeVariables();
					// pcaAnalyzer.subtractAvg();
                    pcaAnalyzer.learnPCABasis(1, 10);
                } catch (XmippError &xe)
                {
                	if (xe.__errno==ERR_NUMERICAL)
                		REPORT_ERROR(ERR_NUMERICAL,"There is no variance in the PSD, check that the micrograph is not constant");
                	else
                		throw(xe);
                }

#ifdef DEBUG

                Image<double> save;
                save().initZeros(psd());
                int ii=-1;
                FOR_ALL_ELEMENTS_IN_ARRAY2D(PCAmask)
                if (PCAmask(i,j))
                    save(i,j)=pcaAnalyzer.PCAbasis[0](++ii);
                save.write("PPPbasis.xmp");
#endif

                Matrix2D<double> CtY;
                pcaAnalyzer.projectOnPCABasis(CtY);
                Matrix1D<double> p;
                CtY.toVector(p);
                double pavg = p.sum(true);
                double pstd = p.sum2() / VEC_XSIZE(p) - pavg * pavg;
                pstd = (pstd < 0) ? 0 : sqrt(pstd);

                std::string psign;
                FOR_ALL_ELEMENTS_IN_MATRIX1D(p)
                if (p(i) < 0)
                    psign += "-";
                else
                    psign += "+";
                double zrandomness = checkRandomness(psign);

                ctfmodel.isLocalCTF = false;
                ctfmodel.x0 = 0;
                ctfmodel.xF = (Xdim-1);
                ctfmodel.y0 = 0;
                ctfmodel.yF = (Ydim-1);
                ROUT_Adjust_CTF(prmEstimateCTFFromPSD,ctfmodel, false);

                // Evaluate PSD variance and write into the CTF
                double stdQ = 0;
                FOR_ALL_ELEMENTS_IN_ARRAY2D(mpsd_std)
                stdQ += A2D_ELEM(mpsd_std,i,j)/A2D_ELEM(mpsd_avg,i,j);
                stdQ /= MULTIDIM_SIZE(psd_std());

                MetaData MD;
                MD.read(fn_psd.withoutExtension() + ".ctfparam");
                size_t id = MD.firstObject();
                MD.setValue(MDL_CTF_CRIT_PSDVARIANCE, stdQ, id);
                MD.setValue(MDL_CTF_CRIT_PSDPCA1VARIANCE, pstd, id);
                MD.setValue(MDL_CTF_CRIT_PSDPCARUNSTEST, zrandomness, id);
                MD.write((String)"fullMicrograph@"+fn_psd.withoutExtension() + ".ctfparam");
            }
            else
            {
                // If bootstrapping
                MultidimArray<double> CTFs(bootstrapN, 32);
                prmEstimateCTFFromPSD.bootstrap = true;
                prmEstimateCTFFromPSD.show_optimization = true;
                FileName fnBase = fn_psd.withoutExtension();
                std::cerr << "Computing bootstrap ...\n";
                init_progress_bar(bootstrapN);
                for (int n = 0; n < bootstrapN; n++)
                {
                    CTFs(n, 31) = ROUT_Adjust_CTF(prmEstimateCTFFromPSD,
                                                  ctfmodel, false);
                    CTFs(n, 0) = ctfmodel.Tm;
                    CTFs(n, 1) = ctfmodel.kV;
                    CTFs(n, 2) = ctfmodel.DeltafU;
                    CTFs(n, 3) = ctfmodel.DeltafV;
                    CTFs(n, 4) = ctfmodel.azimuthal_angle;
                    CTFs(n, 5) = ctfmodel.Cs;
                    CTFs(n, 6) = ctfmodel.Ca;
                    CTFs(n, 7) = ctfmodel.espr;
                    CTFs(n, 8) = ctfmodel.ispr;
                    CTFs(n, 9) = ctfmodel.alpha;
                    CTFs(n, 10) = ctfmodel.DeltaF;
                    CTFs(n, 11) = ctfmodel.DeltaR;
                    CTFs(n, 12) = ctfmodel.Q0;
                    CTFs(n, 13) = ctfmodel.K;
                    CTFs(n, 14) = ctfmodel.gaussian_K;
                    CTFs(n, 15) = ctfmodel.sigmaU;
                    CTFs(n, 16) = ctfmodel.sigmaV;
                    CTFs(n, 17) = ctfmodel.cU;
                    CTFs(n, 18) = ctfmodel.cV;
                    CTFs(n, 19) = ctfmodel.gaussian_angle;
                    CTFs(n, 20) = ctfmodel.sqrt_K;
                    CTFs(n, 21) = ctfmodel.sqU;
                    CTFs(n, 22) = ctfmodel.sqV;
                    CTFs(n, 23) = ctfmodel.sqrt_angle;
                    CTFs(n, 24) = ctfmodel.base_line;
                    CTFs(n, 25) = ctfmodel.gaussian_K2;
                    CTFs(n, 26) = ctfmodel.sigmaU2;
                    CTFs(n, 27) = ctfmodel.sigmaV2;
                    CTFs(n, 28) = ctfmodel.cU2;
                    CTFs(n, 29) = ctfmodel.cV2;
                    CTFs(n, 30) = ctfmodel.gaussian_angle2;

                    std::string command = (std::string) "mv -i " + fnBase
                                          + ".ctfparam " + fnBase + "_bootstrap_"
                                          + integerToString(n, 4) + ".ctfparam";
                    system(command.c_str());
                    command = (std::string) "mv -i " + fnBase
                              + ".ctfmodel_quadrant " + fnBase + "_bootstrap_"
                              + integerToString(n, 4) + ".ctfmodel_quadrant";
                    system(command.c_str());
                    command = (std::string) "mv -i " + fnBase
                              + ".ctfmodel_halfplane " + fnBase + "_bootstrap_"
                              + integerToString(n, 4) + ".ctfmodel_halfplane";
                    system(command.c_str());

                    progress_bar(n);
                }
                progress_bar(bootstrapN);
            }
        }
    }

    // Assign a CTF to each particle ----------------------------------------
    if (psd_mode == OnePerRegion && estimate_ctf)
    {
        double pU0 = 0, pU1 = 0, pU2 = 0;
        planeFit(defocusPlanefittingU, Xm, Ym, pU0, pU1, pU2);
        double pV0 = 0, pV1 = 0, pV2 = 0;
        planeFit(defocusPlanefittingV, Xm, Ym, pV0, pV1, pV2);

        MetaData MDctf;
        MDctf.read(fn_root+".ctfparam");
        double Tm, downsampling;
        size_t id=MDctf.firstObject();
        MDctf.getValue(MDL_CTF_SAMPLING_RATE,Tm,id);
        MDctf.getValue(MDL_CTF_DOWNSAMPLE_PERFORMED,downsampling,id);

        MetaData MD;
        MD.setColumnFormat(false);
        id = MD.addObject();
        MD.setValue(MDL_CTF_DEFOCUS_PLANEUA, pU1, id);
        MD.setValue(MDL_CTF_DEFOCUS_PLANEUB, pU2, id);
        MD.setValue(MDL_CTF_DEFOCUS_PLANEUC, pU0, id);
        MD.setValue(MDL_CTF_DEFOCUS_PLANEVA, pV1, id);
        MD.setValue(MDL_CTF_DEFOCUS_PLANEVB, pV2, id);
        MD.setValue(MDL_CTF_DEFOCUS_PLANEVC, pV0, id);
        MD.setValue(MDL_CTF_X0, 0., id);
        MD.setValue(MDL_CTF_Y0, 0., id);
        MD.setValue(MDL_CTF_XF, (Xdim-1)*Tm*downsampling, id);
        MD.setValue(MDL_CTF_YF, (Ydim-1)*Tm*downsampling, id);
        MD.write((String)"fullMicrograph@"+fn_root+".ctfparam", MD_APPEND);

        if (fn_pos != "")
        {
            FileName fn_img, fn_psd_piece, fn_ctfparam_piece;
            int Y, X;
            FOR_ALL_OBJECTS_IN_METADATA(posFile)
            {
                posFile.getValue(MDL_IMAGE, fn_img, __iter.objId);
                posFile.getValue(MDL_X, X, __iter.objId);
                posFile.getValue(MDL_Y, Y, __iter.objId);
                int idx_X = (int)floor((double) X / pieceDim);
                int idx_Y = (int)floor((double) Y / pieceDim);
                int N = idx_Y * div_NumberX + idx_X + 1;

                fn_psd_piece.compose(N, fn_psd);
                fn_ctfparam_piece = fn_psd_piece.withoutExtension()
                                    + ".ctfparam";
                posFile.setValue(MDL_PSD, fn_psd_piece, __iter.objId);
                posFile.setValue(MDL_CTF_MODEL, fn_ctfparam_piece,
                                 __iter.objId);
            }
        }
    }
    posFile.write(fn_pos);
}

/* Fast estimate of PSD --------------------------------------------------- */
class ThreadFastEstimateEnhancedPSDParams
{
public:
    ImageGeneric *I;
    MultidimArray<double> *PSD, *pieceSmoother;
    MultidimArray<int> *pieceMask;
    Mutex *mutex;
    int Nprocessed;
};

void threadFastEstimateEnhancedPSD(ThreadArgument &thArg)
{
    ThreadFastEstimateEnhancedPSDParams *args =
        (ThreadFastEstimateEnhancedPSDParams*) thArg.workClass;
    int Nthreads = thArg.getNumberOfThreads();
    int id = thArg.thread_id;
    ImageGeneric &I = *(args->I);
    const MultidimArrayGeneric& mI = I();
    size_t IXdim, IYdim, IZdim;
    I.getDimensions(IXdim, IYdim, IZdim);
    MultidimArray<double> &pieceSmoother = *(args->pieceSmoother);
    MultidimArray<int> &pieceMask = *(args->pieceMask);
    MultidimArray<double> localPSD, piece;
    MultidimArray<std::complex<double> > Periodogram;
    piece.initZeros(pieceMask);
    localPSD.initZeros(*(args->PSD));

    FourierTransformer transformer;
    transformer.setReal(piece);

    int pieceNumber = 0;
    int Nprocessed = 0;
    double pieceDim2 = XSIZE(piece) * XSIZE(piece);
    for (size_t i = 0; i < (IYdim - YSIZE(piece)); i+=YSIZE(piece))
        for (size_t j = 0; j < (IXdim - XSIZE(piece)); j+=XSIZE(piece), pieceNumber++)
        {
            if ((pieceNumber + 1) % Nthreads != id)
                continue;
            Nprocessed++;

            // Extract micrograph piece ..........................................
            for (size_t k = 0; k < YSIZE(piece); k++)
                for (size_t l = 0; l < XSIZE(piece); l++)
                    DIRECT_A2D_ELEM(piece, k, l)= mI(i+k, j+l);
            piece.statisticsAdjust(0, 1);
            normalize_ramp(piece, pieceMask);
            piece *= pieceSmoother;

            // Estimate the power spectrum .......................................
            transformer.FourierTransform();
            transformer.getCompleteFourier(Periodogram);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(localPSD)
            {
                double *ptr = (double*) &DIRECT_MULTIDIM_ELEM(Periodogram, n);
                double re=*ptr;
                double im=*(ptr+1);
                double magnitude2=re*re+im*im;
                DIRECT_MULTIDIM_ELEM(localPSD,n)+=magnitude2*pieceDim2;
            }
        }

    // Gather results
    args->mutex->lock();
    args->Nprocessed += Nprocessed;
    *(args->PSD) += localPSD;
    args->mutex->unlock();
}

void fastEstimateEnhancedPSD(const FileName &fnMicrograph, double downsampling,
                             MultidimArray<double> &enhancedPSD, int numberOfThreads)
{
    size_t Xdim, Ydim, Zdim, Ndim;
    getImageSizeFromFilename(fnMicrograph, Xdim, Ydim, Zdim, Ndim);
    int minSize = 2 * (std::max(Xdim, Ydim) / 10);
    minSize = (int)std::min((double) std::min(Xdim, Ydim), NEXT_POWER_OF_2(minSize));
    minSize = std::min(1024, minSize);

    /*
     ProgCTFEstimateFromMicrograph prog1;
     prog1.fn_micrograph=fnMicrograph;
     prog1.fn_root=fnMicrograph.withoutExtension()+"_tmp";
     prog1.pieceDim=(int)(minSize*downsampling);
     prog1.PSDEstimator_mode=ProgCTFEstimateFromMicrograph::Periodogram;
     prog1.Nsubpiece=1;
     prog1.psd_mode=ProgCTFEstimateFromMicrograph::OnePerMicrograph;
     prog1.estimate_ctf=false;
     prog1.bootstrapN=-1;
     prog1.verbose=1;
     prog1.overlap=0;
     prog1.run();
     */
    // Prepare auxiliary variables
    ImageGeneric I;
    I.read(fnMicrograph);

    MultidimArray<double> PSD;
    PSD.initZeros(minSize, minSize);

    MultidimArray<int> pieceMask;
    pieceMask.resizeNoCopy(PSD);
    pieceMask.initConstant(1);

    MultidimArray<double> pieceSmoother;
    constructPieceSmoother(PSD, pieceSmoother);

    // Prepare thread arguments
    Mutex mutex;
    ThreadFastEstimateEnhancedPSDParams args;
    args.I = &I;
    args.PSD = &PSD;
    args.pieceMask = &pieceMask;
    args.pieceSmoother = &pieceSmoother;
    args.Nprocessed = 0;
    args.mutex = &mutex;
    ThreadManager *thMgr = new ThreadManager(numberOfThreads, &args);
    thMgr->run(threadFastEstimateEnhancedPSD);
    if (args.Nprocessed != 0)
        *(args.PSD) /= args.Nprocessed;

    ProgCTFEnhancePSD prog2;
    prog2.filter_w1 = 0.02;
    prog2.filter_w2 = 0.2;
    prog2.decay_width = 0.02;
    prog2.mask_w1 = 0.005;
    prog2.mask_w2 = 0.5;

    prog2.applyFilter(*(args.PSD));
    enhancedPSD = *(args.PSD);

    int downXdim = (int) (XSIZE(enhancedPSD) / downsampling);
    int firstIndex = FIRST_XMIPP_INDEX(downXdim);
    int lastIndex = LAST_XMIPP_INDEX(downXdim);
    enhancedPSD.setXmippOrigin();
    enhancedPSD.selfWindow(firstIndex, firstIndex, lastIndex, lastIndex);
}
