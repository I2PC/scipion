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

#include "recons_misc.h"
#include "symmetrize.h"

#include <data/histogram.h>
#include <data/mask.h>
#include <data/wavelet.h>

/* Fill Reconstruction info structure -------------------------------------- */
void build_recons_info(SelFile &selfile, SelFile &selctf,
                       const FileName &fn_ctf, const SymList &SL,
                       Recons_info * &IMG_Inf, bool do_not_use_symproj)
{
    Matrix2D<double>  L(4, 4), R(4, 4);  // A matrix from the list
    FileName          fn_proj;
    FileName          fn_ctf1;
    Projection        read_proj;
    bool              is_there_ctf = false;
    bool              is_ctf_unique = false;

    int trueIMG = selfile.ImgNo();
    selfile.go_first_ACTIVE();
    int numIMG;
    if (!do_not_use_symproj) numIMG = trueIMG * (SL.SymsNo() + 1);
    else numIMG = trueIMG;

    // The next two ifs check whether there is a CTF file, and
    // whether it is unique
    if (fn_ctf != "")
    {
        is_there_ctf = true;
        is_ctf_unique = true;
        if (fn_ctf.get_extension() == "sel")
        {
            is_ctf_unique = false;
            int truectf = selctf.ImgNo();
            selctf.go_first_ACTIVE();
            int numctf = truectf * (SL.SymsNo() + 1);
            if (numctf != numIMG)
                REPORT_ERROR(9696, "Number of CTF and image files differ");
        }
    }

    if (IMG_Inf != NULL) delete [] IMG_Inf;
    if ((IMG_Inf = new Recons_info[numIMG]) == NULL)
        REPORT_ERROR(3008, "Build_Recons_Info: No memory for the sorting");

    int i = 0; // It will account for the number of valid projections processed
    std::cerr << "Reading angle information ...\n";
    init_progress_bar(trueIMG);
    while (!selfile.eof())
    {
        fn_proj = selfile.NextImg();
        if (fn_proj=="") break;
        if (is_there_ctf && !is_ctf_unique) fn_ctf1 = selctf.NextImg();
        if (fn_proj != "")
        {
            read_proj.read(fn_proj);

            // Filling structure
            IMG_Inf[i].fn_proj = fn_proj;
            if (is_ctf_unique) IMG_Inf[i].fn_ctf = fn_ctf;
            else if (is_there_ctf)  IMG_Inf[i].fn_ctf = fn_ctf1;
            IMG_Inf[i].sym     = -1;
            IMG_Inf[i].seed    = ROUND(65535 * rnd_unif());
            read_proj.get_eulerAngles(IMG_Inf[i].rot, IMG_Inf[i].tilt, IMG_Inf[i].psi);
            EULER_CLIPPING(IMG_Inf[i].rot, IMG_Inf[i].tilt, IMG_Inf[i].psi);

            // Any symmetry?
            if (SL.SymsNo() > 0 && !do_not_use_symproj)
            {
                for (int j = 0; j < SL.SymsNo(); j++)
                {
                    int sym_index = SYMINDEX(SL, j, i, trueIMG);
                    IMG_Inf[sym_index].fn_proj = IMG_Inf[i].fn_proj;
                    IMG_Inf[sym_index].seed   = IMG_Inf[i].seed;
                    if (is_ctf_unique)     IMG_Inf[sym_index].fn_ctf = fn_ctf;
                    else if (is_there_ctf) IMG_Inf[sym_index].fn_ctf = IMG_Inf[i].fn_ctf;
                    IMG_Inf[sym_index].sym = j;
                    SL.get_matrices(j, L, R);
                    L.resize(3, 3); // Erase last row and column
                    R.resize(3, 3); // as only the relative orientation
                    // is useful and not the translation
                    double drot, dtilt, dpsi;
                    Euler_apply_transf(L, R,
                                       IMG_Inf[i].rot, IMG_Inf[i].tilt, IMG_Inf[i].psi,
                                       drot, dtilt, dpsi);
                    IMG_Inf[sym_index].rot = (float)drot;
                    IMG_Inf[sym_index].tilt = (float)dtilt;
                    IMG_Inf[sym_index].psi = (float)dpsi;
                }
            }
        }

        i ++; // I have processed one more image
        if (i % 25 == 0) progress_bar(i);
    }
    progress_bar(trueIMG);
}

/* ------------------------------------------------------------------------- */
VariabilityClass::VariabilityClass(Basic_ART_Parameters *_prm,
                                   int _Zoutput_volume_size, int _Youtput_volume_size,
                                   int _Xoutput_volume_size)
{
    prm = _prm;
    Zoutput_volume_size = _Zoutput_volume_size;
    Youtput_volume_size = _Youtput_volume_size;
    Xoutput_volume_size = _Xoutput_volume_size;
    N = 0;
    VAR_state = VAR_none;
}

void VariabilityClass::newIteration()
{}

//#define MODE7
#define MODE8
//#define MODE15
//#define DEBUG
void VariabilityClass::newUpdateVolume(GridVolume *ptr_vol_out,
                                       Projection &read_proj)
{
    VolumeXmipp vol_voxels;

    // Convert from basis to voxels
    prm->basis.changeToVoxels(*ptr_vol_out, &(vol_voxels()),
                              Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
    (*ptr_vol_out).initZeros();
    N++;

    // Make the DWT
    Matrix3D<double> DWTV;
#ifdef MODE7
    int DWT_iterations = 2;
    int keep_from_iteration = 0;
#endif
#ifdef MODE8
    int DWT_iterations = 2;
    int keep_from_iteration = 0;
#endif
#ifdef MODE15
    int DWT_iterations = 2;
    int keep_from_iteration = 0;
#endif
    Bilib_DWT(vol_voxels(), DWTV, DWT_iterations);
#ifdef DEBUG
    vol_voxels.write("PPPVariability.vol");
    char c;
    std::cout << "Press any key\n";
    std::cin >> c;
#endif

    // Select the LLL block and keep it
    int x1, x2, y1, y2, z1, z2;
    SelectDWTBlock(keep_from_iteration, DWTV, "000", x1, x2, y1, y2, z1, z2);
    DWTV.window(z1, y1, x1, z2, y2, x2);
    STARTINGZ(DWTV) = STARTINGY(DWTV) = STARTINGX(DWTV) = 0;
    VA.push_back(DWTV);
}
#undef DEBUG

#define DEBUG
void VariabilityClass::finishAnalysis()
{
    if (VA.size() == 0) return;

    // Coocurrence matrix
    int nmax = VA.size();
    Matrix2D<int> coocurrence(nmax, nmax);

    // Study each voxel
    VolumeXmipp SignificantT2, SignificantMaxRatio, SignificantMinRatio;
    int zsize = ZSIZE(VA[0]) / 2;
    int ysize = YSIZE(VA[0]) / 2;
    int xsize = XSIZE(VA[0]) / 2;
    int zsize2 = ZSIZE(VA[0]) / 4;
    int ysize2 = YSIZE(VA[0]) / 4;
    int xsize2 = XSIZE(VA[0]) / 4;
    SignificantT2().initZeros(zsize, ysize, xsize);
    SignificantMaxRatio().initZeros(zsize, ysize, xsize);
    SignificantMinRatio().initZeros(zsize, ysize, xsize);
    std::cerr << "Classifying voxels ...\n";
    init_progress_bar(MULTIDIM_SIZE(SignificantT2()));
    int counter = 0, nxny = ysize * zsize;
#ifdef MODE7
#define NFEATURES 7
#endif
#ifdef MODE8
#define NFEATURES 8
#endif
#ifdef MODE15
#define NFEATURES 15
#endif
    FOR_ALL_ELEMENTS_IN_MATRIX3D(SignificantT2())
    {
        // Save the data for this voxel
        std::ofstream fh_dat;
        fh_dat.open("PPP.dat");
        if (!fh_dat)
            REPORT_ERROR(1, "VariabilityClass::finishAnalysis: "
                         "Cannot open PPP.dat for output");
        fh_dat << NFEATURES << " " << nmax << std::endl;
        std::vector< Matrix1D<double> > v;
        v.clear();
        for (int n = 0; n < nmax; n++)
        {
            Matrix1D<double> v_aux(NFEATURES);

#ifdef MODE7
            v_aux(0) = VA[n](k,      i, j + xsize);
            v_aux(1) = VA[n](k, i + ysize,      j);
            v_aux(2) = VA[n](k, i + ysize, j + xsize);
            v_aux(3) = VA[n](k + zsize,      i,      j);
            v_aux(4) = VA[n](k + zsize,      i, j + xsize);
            v_aux(5) = VA[n](k + zsize, i + ysize,      j);
            v_aux(6) = VA[n](k + zsize, i + ysize, j + xsize);
#endif
#ifdef MODE8
            v_aux(0) = VA[n](k,      i,      j);
            v_aux(1) = VA[n](k,      i, j + xsize);
            v_aux(2) = VA[n](k, i + ysize,      j);
            v_aux(3) = VA[n](k, i + ysize, j + xsize);
            v_aux(4) = VA[n](k + zsize,      i,      j);
            v_aux(5) = VA[n](k + zsize,      i, j + xsize);
            v_aux(6) = VA[n](k + zsize, i + ysize,      j);
            v_aux(7) = VA[n](k + zsize, i + ysize, j + xsize);
#endif
#ifdef MODE15
            v_aux(0) = VA[n](k / 2, i / 2,    j / 2);
            v_aux(1) = VA[n](k / 2, i / 2, j / 2 + xsize2);
            v_aux(2) = VA[n](k / 2, i / 2 + ysize2,    j / 2);
            v_aux(3) = VA[n](k / 2, i / 2 + ysize2, j / 2 + xsize2);
            v_aux(4) = VA[n](k / 2 + zsize2,    i / 2,       j / 2);
            v_aux(5) = VA[n](k / 2 + zsize2,    i / 2, j / 2 + xsize2);
            v_aux(6) = VA[n](k / 2 + zsize2, i / 2 + ysize2,       j / 2);
            v_aux(7) = VA[n](k / 2 + zsize2, i / 2 + ysize2, j / 2 + xsize2);
            v_aux(8) = VA[n](k,      i, j + xsize);
            v_aux(9) = VA[n](k, i + ysize,      j);
            v_aux(10) = VA[n](k, i + ysize, j + xsize);
            v_aux(11) = VA[n](k + zsize,      i,      j);
            v_aux(12) = VA[n](k + zsize,      i, j + xsize);
            v_aux(13) = VA[n](k + zsize, i + ysize,      j);
            v_aux(14) = VA[n](k + zsize, i + ysize, j + xsize);
#endif
            v_aux = v_aux * v_aux;
            // COSS: Doesn't work: v_aux=v_aux.sort();

            fh_dat << v_aux.transpose() << std::endl;
            v.push_back(v_aux);
        }
        fh_dat.close();

        // Classify
        system("xmipp_fcmeans -din PPP.dat -std::cout PPP -c 2 -saveclusters > /dev/null");

        // Pick up results
        Matrix2D<double> aux;
        int n_previous;
#define GET_RESULTS(fh,fn,avg,cov,N,idx) \
    fh.open(fn);\
    if (!fh) \
        REPORT_ERROR(1,(std::string)"VariabilityClass::finishAnalysis: " \
                     "Cannot open "+fn+" for input"); \
    n_previous=-1; \
    while (!fh.eof()) { \
        int n; fh >> n; \
        if (n!=n_previous) { \
            n_previous=n; N++; \
            idx(n)=1; \
            avg+=v[n]; \
            aux.fromVector(v[n]); \
            cov+=aux*aux.transpose(); \
        }\
    } \
    avg/=N; \
    cov/=N; \
    aux.fromVector(avg); \
    cov-=aux*aux.transpose();

        std::ifstream fh_0;
        Matrix1D<double> avg0(NFEATURES);
        Matrix1D<int>    idx0(nmax);
        Matrix2D<double> covariance0(NFEATURES, NFEATURES);
        int N0 = 0;
        GET_RESULTS(fh_0, "PPP.0", avg0, covariance0, N0, idx0);
#ifdef DEBUG
        std::cout << "Class 0 is:\n";
        for (int n = 0; n < XSIZE(idx0); n++)
        {
            if (idx0(n))
            {
                int iact_proj = prm->ordered_list(n);
                std::cout << prm->IMG_Inf[iact_proj].fn_proj << std::endl;
            }
        }
#endif

        std::ifstream fh_1;
        Matrix1D<double> avg1(NFEATURES);
        Matrix1D<int>    idx1(nmax);
        Matrix2D<double> covariance1(NFEATURES, NFEATURES);
        int N1 = 0;
        GET_RESULTS(fh_1, "PPP.1", avg1, covariance1, N1, idx1);
#ifdef DEBUG
        std::cout << "Class 1 is:\n";
        for (int n = 0; n < XSIZE(idx1); n++)
        {
            if (idx1(n))
            {
                int iact_proj = prm->ordered_list(n);
                std::cout << prm->IMG_Inf[iact_proj].fn_proj << std::endl;
            }
        }
#endif

        Matrix2D<double> T2, covariance;
        if (NFEATURES > 1)
        {
            // Perform T2-Hotelling test
            Matrix1D<double> avg_diff = avg1 - avg0;
            covariance = 1.0 / (N0 + N1 - 2) *
                         ((N0 - 1) * covariance0 + (N1 - 1) * covariance1);
            covariance *= (1.0 / N0 + 1.0 / N1);
            aux.fromVector(avg_diff);
            T2 = (double)(N0 + N1 - XSIZE(avg_diff) - 1) /
                 ((N0 + N1 - 2) * XSIZE(avg_diff)) *
                 aux.transpose() * covariance.inv() * aux;
        }
        else
        {
            // Perform t-test
            double variance = ((N0 - 1) * covariance0(0, 0) + (N1 - 1) * covariance1(0, 0)) /
                              (N0 + N1 - 2);
            double t = (avg1(0) - avg0(0)) / sqrt(variance * (1.0 / N0 + 1.0 / N1));
            T2.initZeros(1, 1);
            T2(0, 0) = t;
        }

        // Analysis of the covariance structure
        Matrix1D<double> eigenvalues;
        if (NFEATURES > 1)
        {
            Matrix2D<double> U, V;
            svdcmp(covariance, U, eigenvalues, V);
        }

        // Analysis of the coocurrences
        for (int n = 0; n < XSIZE(idx0); n++)
            for (int np = n + 1; np < XSIZE(idx0); np++)
                if (idx0(n) && idx0(np)) coocurrence(n, np)++;

        for (int n = 0; n < XSIZE(idx1); n++)
            for (int np = n + 1; np < XSIZE(idx1); np++)
                if (idx1(n) && idx1(np)) coocurrence(n, np)++;

        // Keep results
        SignificantT2(k, i, j) = T2(0, 0);
        if (NFEATURES > 1)
        {
            SignificantMinRatio(k, i, j) = eigenvalues(1) / eigenvalues(0);
            SignificantMaxRatio(k, i, j) = eigenvalues(NFEATURES - 1) / eigenvalues(0);
        }
#ifdef DEBUG
        std::cout << "T2 for this classification is " << T2(0, 0) << std::endl;
        std::cout << "Eigenvalues are " << eigenvalues.transpose() << std::endl;
#endif
        if (++counter % nxny == 0) progress_bar(counter);
    }
    progress_bar(MULTIDIM_SIZE(SignificantT2()));
    SignificantT2.write(prm->fn_root + "_variability.vol");
    SignificantMinRatio.write(prm->fn_root + "_minratio.vol");
    SignificantMaxRatio.write(prm->fn_root + "_maxratio.vol");
    system("rm PPP.dat PPP.cod PPP.vs PPP.err PPP.his PPP.inf PPP.0 PPP.1");

    for (int n = 0; n < nmax; n++)
        for (int np = n + 1; np < nmax; np++)
            coocurrence(np, n) = coocurrence(n, np);
    ImageXmipp save;
    typeCast(coocurrence, save());
    save.write(prm->fn_root + "_coocurrence.xmp");
}
#undef DEBUG

/* ------------------------------------------------------------------------- */
#define POCS_N_measure  8
#define POCS_N_use      2
/* Constructor ............................................................. */
POCSClass::POCSClass(Basic_ART_Parameters *_prm,
                     int _Zoutput_volume_size, int _Youtput_volume_size,
                     int _Xoutput_volume_size)
{
    prm = _prm;
    POCS_state = POCS_measuring;
    POCS_freq = prm->POCS_freq;
    POCS_i = 0;
    POCS_vec_i = 0;
    POCS_used = 0;
    POCS_N = 0;
    POCS_errors.initZeros(POCS_N_measure);
    Zoutput_volume_size = _Zoutput_volume_size;
    Youtput_volume_size = _Youtput_volume_size;
    Xoutput_volume_size = _Xoutput_volume_size;
    apply_POCS = (prm->surface_mask != NULL ||
                  prm->positivity || (prm->force_sym != 0 && !prm->is_crystal) ||
                  prm->known_volume != -1);
}

void POCSClass::newIteration()
{
    POCS_global_mean_error = 0;
}

void POCSClass::newProjection()
{
    POCS_N = 0;
}

/* Apply ................................................................... */
void POCSClass::apply(GridVolume &vol_basis, int it, int images)
{
    VolumeXmipp vol_POCS, theo_POCS_vol, corr_POCS_vol, vol_voxels;

    if (apply_POCS && POCS_i % POCS_freq == 0)
    {
        VolumeXmipp vol_aux;
        VolumeXmipp *desired_volume = NULL;

        // Compute the corresponding voxel volume
        prm->basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                  Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
        if (prm->tell&TELL_SAVE_AT_EACH_STEP)
        {
            vol_voxels.write("PPPvolPOCS0.vol");
            std::cout << "Stats PPPvolPOCS0.vol: ";
            vol_voxels().printStats();
            std::cout << std::endl;
        }
        // Apply surface restriction
        if (prm->surface_mask != NULL)
        {
            vol_POCS() = (*(prm->surface_mask))();
        }
        else
        {
            vol_POCS().resize(vol_voxels());
            vol_POCS().initZeros();
        }
        if (prm->tell&TELL_SAVE_AT_EACH_STEP)
        {
            vol_POCS.write("PPPvolPOCS1.vol");
            std::cout << "Stats PPPvolPOCS1.vol: ";
            vol_POCS().printStats();
            std::cout << std::endl;
        }
        // Force symmetry
        if (prm->force_sym != 0)
        {
            symmetrize(prm->SL, vol_voxels, vol_aux, false);
            desired_volume = &vol_aux;
            if (prm->tell&TELL_SAVE_AT_EACH_STEP)
            {
                vol_aux.write("PPPvolPOCS2.vol");
                std::cout << "Stats PPPvolPOCS2.vol: ";
                vol_aux().printStats();
                std::cout << std::endl;
            }
        }
        // Apply volume constraint
        if (prm->known_volume != -1)
        {
            histogram1D hist;
            Matrix3D<int> aux_mask;
            aux_mask.resize(vol_POCS());
            FOR_ALL_ELEMENTS_IN_MATRIX3D(aux_mask)
            aux_mask(k, i, j) = 1 - (int)vol_POCS(k, i, j);
            long mask_voxels = vol_POCS().countThreshold("below", 0.5, 0);
            compute_hist_within_binary_mask(
                aux_mask, vol_voxels(), hist, 300);
            double known_percentage;
            known_percentage = XMIPP_MIN(100, 100 * prm->known_volume / mask_voxels);
            double threshold;
            threshold = hist.percentil(100 - known_percentage);
            FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_voxels())
            if (vol_voxels(k, i, j) < threshold) vol_POCS(k, i, j) = 1;
        }
        if (prm->tell&TELL_SAVE_AT_EACH_STEP)
        {
            vol_POCS.write("PPPvolPOCS3.vol");
            std::cout << "Stats PPPvolPOCS3.vol: ";
            vol_POCS().printStats();
            std::cout << std::endl;
        }

        // Do not allow positivity outside interest region
        // and do not allow negativity inside the interest region
        // if positivity restrictions are to be applied
        int bg = (int) vol_POCS().sum();
        int fg = MULTIDIM_SIZE(vol_POCS()) - bg;
        int relax = 0, posi = 0;
        FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_voxels())
        if (vol_POCS(k, i, j) == 1 && vol_voxels(k, i, j) < 0)
        {
            vol_POCS(k, i, j) = 0;
            relax++;
        }
        else if (vol_POCS(k, i, j) == 0 && vol_voxels(k, i, j) < 0 &&
                 prm->positivity)
        {
            vol_POCS(k, i, j) = 1;
            posi++;
        }
        // Debugging messages
        //std::cerr << "Relaxation/Positivity " << (double)relax/(double)bg << " "
        //     << (double)posi/(double)fg << " " << std::endl;

        // Solve volumetric equations
        switch (prm->basis.type)
        {
        case Basis::blobs:
            if (desired_volume == NULL)
                ART_voxels2blobs_single_step(vol_basis, &vol_basis,
                                             prm->basis.blob, prm->D, prm->lambda(it),
                                             &(theo_POCS_vol()), NULL,
                                             &(corr_POCS_vol()),
                                             &(vol_POCS()),
                                             POCS_mean_error, POCS_max_error, VARTK);
            else
            {
                FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_POCS())
                if (vol_POCS(k, i, j) == 1)(*desired_volume)(k, i, j) = 0;
                for (int i = 0; i < prm->force_sym; i++)
                {
                    ART_voxels2blobs_single_step(vol_basis, &vol_basis,
                                                 prm->basis.blob, prm->D, prm->lambda(it),
                                                 &(theo_POCS_vol()), &((*desired_volume)()),
                                                 &(corr_POCS_vol()),
                                                 NULL,
                                                 POCS_mean_error, POCS_max_error, VARTK);
                    if (prm->tell&TELL_SAVE_AT_EACH_STEP)
                        std::cout << "    POCS Iteration " << i
                                  << " POCS Error=" <<  POCS_mean_error << std::endl;
                }
            }
            break;
        case Basis::voxels:
            if (desired_volume == NULL)
            {
                FOR_ALL_ELEMENTS_IN_MATRIX3D(vol_POCS())
                if (vol_POCS(k, i, j)) vol_basis(0)(k, i, j) = 0;
            }
            else
            {
                vol_basis(0)().initZeros();
                FOR_ALL_ELEMENTS_IN_MATRIX3D((*desired_volume)())
                vol_basis(0)(k, i, j) = (*desired_volume)(k, i, j);
            }
            POCS_mean_error = -1;
            break;
        }
        POCS_i = 1;
        POCS_global_mean_error += POCS_mean_error;
        POCS_N++;

        // Now some control logic
        if (prm->numIMG - images < 100 || images % 100 == 0 ||
            desired_volume != NULL)
        {
            POCS_freq = 1;
            POCS_state = POCS_measuring;
            POCS_vec_i = 0;
        }
        else
        {
            double dummy;
            switch (POCS_state)
            {
            case POCS_measuring:
#ifdef DEBUG_POCS
                std::cout << "M:" << POCS_vec_i << " " << POCS_mean_error << std::endl;
#endif
                POCS_errors(POCS_vec_i++) = POCS_mean_error;
                if (POCS_vec_i == POCS_N_measure)
                {
                    POCS_vec_i = 0;
                    // Change to use state
                    POCS_used = 0;
                    POCS_freq++;
                    POCS_state = POCS_use;
#ifdef DEBUG_POCS
                    std::cerr << "1: Changing to " << POCS_freq << std::endl;
#endif
                }
                break;
            case POCS_use:
                POCS_used++;
                POCS_errors.computeStats(POCS_avg,
                                         POCS_stddev, dummy, POCS_min);
#ifdef DEBUG_POCS
                std::cout << "Reference errors: " << POCS_errors.transpose() << std::endl;
                std::cout << "Checking " << ABS(POCS_mean_error - POCS_avg) << " " << 1.2*1.96*POCS_stddev << std::endl;
#endif
                if (ABS(POCS_mean_error - POCS_avg) < 1.2*1.96*POCS_stddev)
                {
                    if (POCS_mean_error < POCS_avg)
                    {
                        double max_error = POCS_errors(0);
                        POCS_vec_i = 0;
                        for (int i = 1; i < POCS_N_measure; i++)
                            if (POCS_errors(i) > max_error)
                            {
                                max_error = POCS_errors(i);
                                POCS_vec_i = i;
                            }
                        POCS_errors(POCS_vec_i) = POCS_mean_error;
                    }
                    if (POCS_used < POCS_N_use)
                    {
                        // While not enough uses
                    }
                    else if (POCS_freq < 3)
                    {
                        // increase frequency
                        POCS_freq++;
#ifdef DEBUG_POCS
                        std::cerr << "2: Changing to " << POCS_freq << std::endl;
#endif
                        POCS_used = 0;
                    }
                }
                else
                {
                    // It is behaving worse
                    if (POCS_freq > prm->POCS_freq + 1)
                    {
                        POCS_freq = prm->POCS_freq + 1;
                        POCS_used = 0;
#ifdef DEBUG_POCS
                        std::cerr << "3: Changing to " << POCS_freq << std::endl;
#endif
                    }
                    else if (POCS_used > 2)
                    {
                        POCS_freq = prm->POCS_freq;
                        // Change status
                        POCS_used = 0;
                        POCS_state = POCS_lowering;
#ifdef DEBUG_POCS
                        std::cerr << "Lowering\n";
#endif
                    }
                }
                break;
            case POCS_lowering:
                // Lower the POCS error before measuring again
                POCS_errors.computeStats(POCS_avg,
                                         POCS_stddev, POCS_max, dummy);
                POCS_used++;
                if (POCS_mean_error < POCS_max || POCS_used > 2*POCS_N_measure)
                {
                    // Change status
                    POCS_vec_i = 0;
                    POCS_state = POCS_measuring;
                }
                break;
            }
        }
    }
    else
    {
        POCS_i++;
        POCS_mean_error = -1;
    }
}

