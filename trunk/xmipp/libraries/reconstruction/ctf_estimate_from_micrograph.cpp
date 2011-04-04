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

#include <data/args.h>
#include <data/micrograph.h>
#include <data/metadata.h>
#include <data/image.h>
#include <data/fft.h>
#include <data/basic_pca.h>

/* Read parameters ========================================================= */
void ProgCTFEstimateFromMicrograph::readParams()
{
    fn_micrograph=getParam("--micrograph");
    pieceDim = getIntParam("--pieceDim");
    if (getParam("--psd_estimator")=="periodogram")
        PSDEstimator_mode = Periodogram;
    else
    {
        PSDEstimator_mode = ARMA;
        ARMA_prm.readParams(this);
    }
    Nsubpiece          = getIntParam("--Nsubpiece");

    String mode=getParam("--mode");
    if (mode=="micrograph")
        psd_mode=OnePerMicrograph;
    else if (mode=="regions")
        psd_mode=OnePerRegion;
    else if (mode=="particles")
    {
        psd_mode=OnePerParticle;
        fn_pos=getParam("--mode",1);
    }
    estimate_ctf=!checkParam("--dont_estimate_ctf");
    if (estimate_ctf)
        prmEstimateCTFFromPSD.readParams();
    bootstrapN     = getIntParam("--bootstrapFit");
}

void ProgCTFEstimateFromMicrograph::defineParams()
{
    addUsageLine("Estimate the CTF from a micrograph.");
    addParamsLine("   --micrograph <file>         : File with the micrograph");
    addParamsLine("==+ PSD estimation");
    addParamsLine("  [--psd_estimator <method=periodogram>] : Method for estimating the PSD");
    addParamsLine("         where <method>");
    addParamsLine("                  periodogram");
    addParamsLine("                  ARMA");
    addParamsLine("  [--pieceDim <d=512>]       : Size of the piece");
    addParamsLine("  [--Nsubpiece <N=1>]        : Each piece is further subdivided into NxN subpieces.");
    addParamsLine("                              : This option is useful for small micrographs in which ");
    addParamsLine("                              : not many pieces of size pieceDim x pieceDim can be defined. ");
    addParamsLine("                              :++ Note that this is not the same as defining a smaller pieceDim. ");
    addParamsLine("                              :++ Defining a smaller pieceDim, would result in a small PSD, while ");
    addParamsLine("                              :++ subdividing the piece results in a large PSD, although smoother.");
    addParamsLine("  [--mode <mode=micrograph>]");
    addParamsLine("         where <mode>");
    addParamsLine("                  micrograph  : Single PSD for the whole micrograph");
    addParamsLine("                  regions     : The micrograph is divided into a region grid ");
    addParamsLine("                              : and a PSD is computed for each one.");
    addParamsLine("                  particles <file> : One PSD per particle.");
    addParamsLine("                              : The file is metadata with the position of each particle within the micrograph");
    addParamsLine("==+ CTF fit");
    addParamsLine("  [--dont_estimate_ctf]       : Do not fit a CTF to PSDs");
    ARMA_parameters::defineParams(this);
    ProgCTFEstimateFromPSD::defineBasicParams(this);
}

/* Compute PSD by piece averaging ========================================== */
//#define DEBUG
void ProgCTFEstimateFromMicrograph::PSD_piece_by_averaging(MultidimArray<double> &piece,
        MultidimArray<double> &psd)
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

    MultidimArray< std::complex<double> > Periodogram;
    MultidimArray<double> small_psd;
    for (int ii = 0; ii < Nsubpiece; ii++)
        for (int jj = 0; jj < Nsubpiece; jj++)
        {
            // Take the corresponding small piece from the piece
            int i0 = ii * Xstep;
            int j0 = jj * Ystep;

            int i, j, ib, jb;
            for (i = 0, ib = i0; i < small_Ydim; i++, ib++)
                for (j = 0, jb = j0; j < small_Xdim; j++, jb++)
                    DIRECT_A2D_ELEM(small_piece, i, j) =
                        DIRECT_A2D_ELEM(piece, ib, jb);

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
    psd *= 1.0/(Nsubpiece * Nsubpiece);

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
    int Zdim, Ydim, Xdim; // Micrograph dimensions
    M_in.readMapped(fn_micrograph);
    M_in.getDimensions(Xdim, Ydim,Zdim);

    // Compute the number of divisions --------------------------------------
    int div_Number = 0;
    int div_NumberX, div_NumberY;
    if (psd_mode==OnePerParticle)
        div_Number=posFile.size();
    else if (psd_mode==OnePerMicrograph)
    {
        div_NumberX = CEIL((double)Xdim / (pieceDim / 2)) - 1;
        div_NumberY = CEIL((double)Ydim / (pieceDim / 2)) - 1;
        div_Number = div_NumberX * div_NumberY;
    }
    else if (psd_mode==OnePerRegion)
    {
        div_NumberX = CEIL((double)Xdim / pieceDim);
        div_NumberY = CEIL((double)Ydim / pieceDim);
        div_Number = div_NumberX * div_NumberY;
    }

    // Process each piece ---------------------------------------------------
    Image<double> psd_avg, psd_std, psd, psd2;
    MultidimArray< std::complex<double> > Periodogram;
    MultidimArray<double> piece(pieceDim, pieceDim);
    psd().resizeNoCopy(piece);
    MultidimArray<double> &mpsd=psd();
    PCAMahalanobisAnalyzer pcaAnalyzer;
    MultidimArray<int> PCAmask;
    MultidimArray<float> PCAv;
    std::cerr << "Computing models of each piece ...\n";

    // Prepare these filenames in case they are needed
    FileName fn_root= fn_micrograph.withoutExtension();
    FileName fn_psd;
    if (psd_mode==OnePerMicrograph)
        fn_psd=fn_root+".psd";
    else
        fn_psd=fn_root+".psdstk";

    init_progress_bar(div_Number);
    int N = 1; // Index of current piece
    int i = 0, j = 0; // top-left corner of the current piece
    while (N <= div_Number)
    {
        // Compute the top-left corner of the piece ..........................
        if (psd_mode==OnePerParticle)
        {
            // Read position of the particle
            posFile.getValue(MDL_X,j,iterPosFile.objId);
            posFile.getValue(MDL_Y,i,iterPosFile.objId);

            // j,i are the selfWindow center, we need the top-left corner
            j -= (int)(pieceDim / 2);
            i -= (int)(pieceDim / 2);
            if (i < 0)
                i = 0;
            if (j < 0)
                j = 0;
        }
        else
        {
            int step = pieceDim;
            if (psd_mode==OnePerMicrograph)
                step /= 2;
            i = ((N - 1) / div_NumberX) * step;
            j = ((N - 1) % div_NumberX) * step;
        }

        // test if the full piece is inside the micrograph
        if (i + pieceDim > Ydim)
            i = Ydim - pieceDim;
        if (j + pieceDim > Xdim)
            j = Xdim - pieceDim;

        // Extract micrograph piece ..........................................
        const MultidimArrayGeneric& mM_in=M_in();
        for (int k = 0; k < YSIZE(piece); k++)
            for (int l = 0; l < XSIZE(piece); l++)
                DIRECT_A2D_ELEM(piece, k, l) = mM_in(i+k, j+l);
        piece.statisticsAdjust(0, 1);

        // Estimate the power spectrum .......................................
        if (Nsubpiece>1)
            if (PSDEstimator_mode == ARMA)
            {
                CausalARMA(piece, ARMA_prm);
                ARMAFilter(piece, psd(), ARMA_prm);
            }
            else
            {
                FourierTransform(piece, Periodogram);
                FFT_magnitude(Periodogram, psd());
                psd() *= psd();
                psd() *= pieceDim * pieceDim;
            }
        else
            PSD_piece_by_averaging(piece, psd());
        psd2()=psd();
        psd2()*=psd();

        // Perform averaging if applicable ...................................
        if (psd_mode==OnePerMicrograph)
        {
            // Compute average and standard deviation
            if (N == 1)
            {
                psd_avg() = psd();
                psd_std() = psd2();
            }
            else
            {
                psd_avg() += psd();
                psd_std() += psd2();
            }

            // Keep psd for the PCA
            if (XSIZE(PCAmask)==0)
            {
                PCAmask.initZeros(psd());
                Matrix1D<int>    idx(2);  // Indexes for Fourier plane
                Matrix1D<double> freq(2); // Frequencies for Fourier plane
                size_t PCAdim=0;
                FOR_ALL_ELEMENTS_IN_ARRAY2D(PCAmask)
                {
                    VECTOR_R2(idx, j, i);
                    FFT_idx2digfreq(psd(), idx, freq);
                    double w = freq.module();
                    if (w>0.05 && w<0.4)
                    {
                        A2D_ELEM(PCAmask,i,j)=1;
                        ++PCAdim;
                    }
                }
                PCAv.initZeros(PCAdim);
            }

            int ii=-1;
            FOR_ALL_ELEMENTS_IN_ARRAY2D(PCAmask)
            if (A2D_ELEM(PCAmask,i,j))
                A1D_ELEM(PCAv,++ii)=(float)A2D_ELEM(mpsd,i,j);
            pcaAnalyzer.addVector(PCAv);
        }

        // Compute the theoretical model if not averaging ....................
        if (psd_mode!=OnePerMicrograph)
        {
            if (bootstrapN!=-1)
                REPORT_ERROR(ERR_VALUE_INCORRECT,
                             "Bootstrapping is only available for micrograph averages");

            FileName fn_psd_piece;
            fn_psd_piece.compose(N,fn_psd);
            psd.write(fn_psd_piece);
            if (psd_mode==OnePerParticle)
                posFile.setValue(MDL_PSD,fn_psd_piece,iterPosFile.objId);

            if (estimate_ctf)
            {
                // Estimate the CTF parameters of this piece
                prmEstimateCTFFromPSD.fn_psd=fn_psd_piece;
                CTFDescription ctfmodel;
                double fitting_error = ROUT_Adjust_CTF(prmEstimateCTFFromPSD,
                                                       ctfmodel, false);
                if (psd_mode==OnePerParticle)
                    posFile.setValue(MDL_CTFMODEL,fn_psd_piece.withoutExtension()+".ctfparam",
                                     iterPosFile.objId);
            }
        }

        // Increment the division counter
        progress_bar(++N);
        if (psd_mode==OnePerParticle)
            iterPosFile.moveNext();
    }
    progress_bar(div_Number);

    // If averaging, compute the CTF model ----------------------------------
    if (psd_mode==OnePerMicrograph)
    {
        // Compute the avg and stddev of the local PSDs
        const MultidimArray<double> &mpsd_std=psd_std();
        const MultidimArray<double> &mpsd_avg=psd_avg();
        double idiv_Number=1.0/div_Number;
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
            if (bootstrapN==-1)
            {
            	// No bootstrapping
                // Compute the PCA of the local PSDs
                pcaAnalyzer.standardarizeVariables();
                // pcaAnalyzer.subtractAvg();
                pcaAnalyzer.learnPCABasis(1,10);

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
                double pavg=p.sum(true);
                double pstd=p.sum2()/VEC_XSIZE(p)-pavg*pavg;
                pstd=(pstd<0)?0:sqrt(pstd);

                std::string psign;
                FOR_ALL_ELEMENTS_IN_MATRIX1D(p)
                if (p(i)<0)
                    psign+="-";
                else
                    psign+="+";
                double zrandomness=checkRandomness(psign);

                double fitting_error = ROUT_Adjust_CTF(prmEstimateCTFFromPSD,
                                                       ctfmodel, false);

                // Evaluate PSD variance and write into the CTF
                double stdQ=0;
                FOR_ALL_ELEMENTS_IN_ARRAY2D(mpsd_std)
                stdQ+=A2D_ELEM(mpsd_std,i,j)/A2D_ELEM(mpsd_avg,i,j);
                stdQ/=MULTIDIM_SIZE(psd_std());

                MetaData MD;
                MD.read(fn_psd.withoutExtension() + ".ctfparam");
                size_t id = MD.firstObject();
                MD.setValue(MDL_CTF_CRITERION_PSDVARIANCE,stdQ,id);
                MD.setValue(MDL_CTF_CRITERION_PSDPCA1VARIANCE,pstd,id);
                MD.setValue(MDL_CTF_CRITERION_PSDPCARUNSTEST,zrandomness,id);
                MD.write(fn_psd.withoutExtension() + ".ctfparam");
            }
            else
            {
            	// If bootstrapping
                MultidimArray<double> CTFs(bootstrapN,32);
                prmEstimateCTFFromPSD.bootstrap=true;
                prmEstimateCTFFromPSD.show_optimization=true;
                FileName fnBase=fn_psd.withoutExtension();
                std::cerr << "Computing bootstrap ...\n";
                init_progress_bar(bootstrapN);
                for (int n=0; n<bootstrapN; n++)
                {
                    CTFs(n,31) = ROUT_Adjust_CTF(prmEstimateCTFFromPSD,
                                                 ctfmodel, false);
                    CTFs(n, 0)=ctfmodel.Tm;
                    CTFs(n, 1)=ctfmodel.kV;
                    CTFs(n, 2)=ctfmodel.DeltafU;
                    CTFs(n, 3)=ctfmodel.DeltafV;
                    CTFs(n, 4)=ctfmodel.azimuthal_angle;
                    CTFs(n, 5)=ctfmodel.Cs;
                    CTFs(n, 6)=ctfmodel.Ca;
                    CTFs(n, 7)=ctfmodel.espr;
                    CTFs(n, 8)=ctfmodel.ispr;
                    CTFs(n, 9)=ctfmodel.alpha;
                    CTFs(n,10)=ctfmodel.DeltaF;
                    CTFs(n,11)=ctfmodel.DeltaR;
                    CTFs(n,12)=ctfmodel.Q0;
                    CTFs(n,13)=ctfmodel.K;
                    CTFs(n,14)=ctfmodel.gaussian_K;
                    CTFs(n,15)=ctfmodel.sigmaU;
                    CTFs(n,16)=ctfmodel.sigmaV;
                    CTFs(n,17)=ctfmodel.cU;
                    CTFs(n,18)=ctfmodel.cV;
                    CTFs(n,19)=ctfmodel.gaussian_angle;
                    CTFs(n,20)=ctfmodel.sqrt_K;
                    CTFs(n,21)=ctfmodel.sqU;
                    CTFs(n,22)=ctfmodel.sqV;
                    CTFs(n,23)=ctfmodel.sqrt_angle;
                    CTFs(n,24)=ctfmodel.base_line;
                    CTFs(n,25)=ctfmodel.gaussian_K2;
                    CTFs(n,26)=ctfmodel.sigmaU2;
                    CTFs(n,27)=ctfmodel.sigmaV2;
                    CTFs(n,28)=ctfmodel.cU2;
                    CTFs(n,29)=ctfmodel.cV2;
                    CTFs(n,30)=ctfmodel.gaussian_angle2;

                    std::string command=(std::string)"mv -i "+fnBase+
                                        ".ctfparam "+fnBase+"_bootstrap_"+
                                        integerToString(n,4)+".ctfparam";
                    system(command.c_str());
                    command=(std::string)"mv -i "+fnBase+
                            ".ctfmodel_quadrant "+fnBase+"_bootstrap_"+
                            integerToString(n,4)+".ctfmodel_quadrant";
                    system(command.c_str());
                    command=(std::string)"mv -i "+fnBase+
                            ".ctfmodel_halfplane "+fnBase+"_bootstrap_"+
                            integerToString(n,4)+".ctfmodel_halfplane";
                    system(command.c_str());

                    progress_bar(n);
                }
                progress_bar(bootstrapN);
            }
        }
    }

    // Assign a CTF to each particle ----------------------------------------
    if (psd_mode==OnePerRegion && fn_pos != "" && estimate_ctf)
    {
        FileName fn_img, fn_psd_piece, fn_ctfparam_piece;
        int Y, X;
        FOR_ALL_OBJECTS_IN_METADATA(posFile)
        {
        	posFile.getValue(MDL_IMAGE, fn_img, __iter.objId);
        	posFile.getValue(MDL_X,X,__iter.objId);
        	posFile.getValue(MDL_Y,Y,__iter.objId);
            int idx_X = floor((double)X / pieceDim);
            int idx_Y = floor((double)Y / pieceDim);
            int N = idx_Y * div_NumberX + idx_X + 1;

            fn_psd_piece.compose(N,fn_psd);
            fn_ctfparam_piece=fn_psd_piece.withoutExtension()+".ctfparam";
            posFile.setValue(MDL_PSD,fn_psd_piece,__iter.objId);
            posFile.setValue(MDL_CTFMODEL,fn_ctfparam_piece,__iter.objId);
        }
    }
    posFile.write(fn_pos);
}
