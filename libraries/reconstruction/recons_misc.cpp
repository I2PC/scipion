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


/* Fill Reconstruction info structure -------------------------------------- */
void buildReconsInfo(MetaData &selfile,
                     const FileName &fn_ctf, const SymList &SL,
                     ReconsInfo * &IMG_Inf, bool do_not_use_symproj)
{
    Matrix2D<double>  L(4, 4), R(4, 4);  // A matrix from the list
    FileName          fn_proj;
    FileName          fn_ctf1;
    Projection        read_proj;
    bool              is_there_ctf = false;
    bool              is_ctf_unique = false;

    int trueIMG = selfile.size();
    selfile.firstObject();
    int numIMG;
    if (!do_not_use_symproj)
        numIMG = trueIMG * (SL.symsNo() + 1);
    else
        numIMG = trueIMG;

    // The next two ifs check whether there is a CTF file, and
    // whether it is unique
    if (fn_ctf != "")
    {
        is_there_ctf = true;
        is_ctf_unique = true;
    }
    else if (selfile.containsLabel(MDL_CTF_MODEL))
    {
        is_there_ctf = true;
        is_ctf_unique = false;
    }

    if (IMG_Inf != NULL)
        delete [] IMG_Inf;
    if ((IMG_Inf = new ReconsInfo[numIMG]) == NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Build_Recons_Info: No memory for the sorting");

    int i = 0; // It will account for the number of valid projections processed
    std::cout << "Reading angle information ...\n";
    init_progress_bar(trueIMG);
    FOR_ALL_OBJECTS_IN_METADATA(selfile)
    {
        ReconsInfo &imgInfo = IMG_Inf[i];
        selfile.getValue(MDL_IMAGE,fn_proj,__iter.objId);
        if (is_there_ctf && !is_ctf_unique)
            selfile.getValue(MDL_CTF_MODEL,fn_ctf1,__iter.objId);
        if (fn_proj != "")
        {
            //            read_proj.read(fn_proj, false, HEADER);
            // Filling structure
            imgInfo.fn_proj = fn_proj;
            selfile.getRow(imgInfo.row, __iter.objId);
            if (is_ctf_unique)
                imgInfo.fn_ctf = fn_ctf;
            else if (is_there_ctf)
                imgInfo.fn_ctf = fn_ctf1;
            imgInfo.sym     = -1;
            imgInfo.seed    = ROUND(65535 * rnd_unif());

            imgInfo.rot = imgInfo.tilt = imgInfo.psi = 0;

            selfile.getValue(MDL_ANGLE_ROT, imgInfo.rot,__iter.objId);
            selfile.getValue(MDL_ANGLE_TILT, imgInfo.tilt,__iter.objId);
            selfile.getValue(MDL_ANGLE_PSI, imgInfo.psi,__iter.objId);
            //            read_proj.getEulerAngles(imgInfo.rot, imgInfo.tilt, imgInfo.psi);
            EULER_CLIPPING(imgInfo.rot, imgInfo.tilt, imgInfo.psi);

            // Any symmetry?
            if (SL.symsNo() > 0 && !do_not_use_symproj)
            {
                for (int j = 0; j < SL.symsNo(); j++)
                {
                    int sym_index = SYMINDEX(SL, j, i, trueIMG);
                    IMG_Inf[sym_index].fn_proj = imgInfo.fn_proj;
                    IMG_Inf[sym_index].seed   = imgInfo.seed;
                    if (is_ctf_unique)
                        IMG_Inf[sym_index].fn_ctf = fn_ctf;
                    else if (is_there_ctf)
                        IMG_Inf[sym_index].fn_ctf = imgInfo.fn_ctf;
                    IMG_Inf[sym_index].sym = j;
                    SL.getMatrices(j, L, R);
                    L.resize(3, 3); // Erase last row and column
                    R.resize(3, 3); // as only the relative orientation
                    // is useful and not the translation
                    double drot, dtilt, dpsi;
                    Euler_apply_transf(L, R,
                                       imgInfo.rot, imgInfo.tilt, imgInfo.psi,
                                       drot, dtilt, dpsi);
                    IMG_Inf[sym_index].rot = (float)drot;
                    IMG_Inf[sym_index].tilt = (float)dtilt;
                    IMG_Inf[sym_index].psi = (float)dpsi;
                }
            }
        }

        ++i; // I have processed one more image
        if (i % 25 == 0)
            progress_bar(i);
    }
    progress_bar(trueIMG);
}



/* ------------------------------------------------------------------------- */
/* Sort_perpendicular                                                        */
/* ------------------------------------------------------------------------- */
void sortPerpendicular(int numIMG, ReconsInfo *IMG_Inf,
                       MultidimArray<int> &ordered_list, int N)
{
    int   i, j;
    MultidimArray<short> chosen(numIMG);     // 1 if that image has been already
    // chosen
    double min_prod;
    int   min_prod_proj;
    Matrix2D<double> v(numIMG, 3);
    Matrix2D<double> euler;
    MultidimArray<double> product(numIMG);

    // Initialization
    ordered_list.resize(numIMG);
    for (i = 0; i < numIMG; i++)
    {
        Matrix1D<double> z;
        // Initially no image is chosen
        A1D_ELEM(chosen, i) = 0;

        // Compute the Euler matrix for each image and keep only
        // the third row of each one
        //0.f -> double 0. It should be there is the other
        // arguments are doubles because Euler_angles2matrix
        //acepts either all doubles or all doubles
        Euler_angles2matrix(IMG_Inf[i].rot, IMG_Inf[i].tilt, 0.f, euler);
        euler.getRow(2, z);
        v.setRow(i, z);
    }

    // Pick first projection as the first one to be presented
    i = 0;
    A1D_ELEM(chosen, i) = 1;
    A1D_ELEM(ordered_list, 0) = i;

    // Choose the rest of projections
    std::cout << "Sorting projections ...\n";
    init_progress_bar(numIMG - 1);
    Matrix1D<double> rowj, rowi_1, rowi_N_1;
    for (i = 1; i < numIMG; i++)
    {
        // Compute the product of not already chosen vectors with the just
        // chosen one, and select that which has minimum product
        min_prod = MAXFLOAT;
        v.getRow(A1D_ELEM(ordered_list, i - 1),rowi_1);
        if (N != -1 && i > N)
            v.getRow(A1D_ELEM(ordered_list, i - N - 1),rowi_N_1);
        for (j = 0; j < numIMG; j++)
        {
            if (!A1D_ELEM(chosen, j))
            {
                v.getRow(j,rowj);
                A1D_ELEM(product, j) += ABS(dotProduct(rowi_1,rowj));
                if (N != -1 && i > N)
                    A1D_ELEM(product, j) -= ABS(dotProduct(rowi_N_1,rowj));
                if (A1D_ELEM(product, j) < min_prod)
                {
                    min_prod = A1D_ELEM(product, j);
                    min_prod_proj = j;
                }
            }
        }

        // Store the chosen vector and mark it as chosen
        A1D_ELEM(ordered_list, i) = min_prod_proj;
        A1D_ELEM(chosen, min_prod_proj) = 1;

        // The progress bar is updated only every 10 images
        if (i % 10 == 0)
            progress_bar(i);
    }

    // A final call to progress bar to finish a possible small piece
    progress_bar(numIMG - 1);
    std::cout << std::endl;
}

/* ------------------------------------------------------------------------- */
/* No Sort                                                                   */
/* ------------------------------------------------------------------------- */
void noSort(int numIMG, MultidimArray<int> &ordered_list)
{
    ordered_list.initLinear(0, numIMG - 1);
}

/* ------------------------------------------------------------------------- */
/* Random Sort                                                               */
/* ------------------------------------------------------------------------- */
void sortRandomly(int numIMG, MultidimArray<int> &ordered_list)
{
    MultidimArray<int> chosen;

    // Initialisation
    ordered_list.resize(numIMG);
    chosen.initZeros(numIMG);

    std::cout << "Randomizing projections ...\n";
    init_progress_bar(numIMG - 1);
    int ptr = 0;
    randomize_random_generator();
    for (int i = numIMG; i > 0; i--)
    {
        // Jump a random number starting at the pointed projection
        int rnd_indx = (int) rnd_unif(0, i) + 1;
        while (rnd_indx > 0)
        {
            // Jump one not chosen image
            ptr = (ptr + 1) % numIMG;
            // Check it is not chosen, if it is, go on skipping
            while (chosen(ptr))
            {
                ptr = (ptr + 1) % numIMG;
            }
            rnd_indx--;
        }

        // Annotate this image
        A1D_ELEM(ordered_list, i - 1) = ptr;
        A1D_ELEM(chosen, ptr) = 1;

        // The progress bar is updated only every 10 images
        if (i % 10 == 0)
            progress_bar(i);
    }

    // A final call to progress bar to finish a possible small piece
    progress_bar(numIMG - 1);
    std::cout << std::endl;
}

/* ------------------------------------------------------------------------- */
/* Update residual vector for WLS                                            */
/* ------------------------------------------------------------------------- */
void updateResidualVector(BasicARTParameters &prm, GridVolume &vol_basis,
                          double &kappa, double &pow_residual_vol, double &pow_residual_imgs)
{
    GridVolume       residual_vol;
    Projection       read_proj, dummy_proj, new_proj;
    FileName         fn_resi, fn_tmp;
    double           sqrtweight, dim2;
    Matrix2D<double> *A = NULL;
    std::vector<MultidimArray<double> > newres_imgs;
    MultidimArray<int>    mask;

    residual_vol.resize(vol_basis);
    residual_vol.initZeros();

    // Calculate volume from all backprojected residual images
    std::cout << "Backprojection of residual images " << std::endl;
    if (!(prm.tell&TELL_SHOW_ERROR))
        init_progress_bar(prm.numIMG);

    for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++)
    {
        // backprojection of the weighted residual image
        sqrtweight = sqrt(prm.residual_imgs[iact_proj].weight() / prm.sum_weight);
        read_proj = prm.residual_imgs[iact_proj];
        read_proj() *= sqrtweight;
        dummy_proj().resize(read_proj());

        dummy_proj.setAngles(prm.IMG_Inf[iact_proj].rot,
                              prm.IMG_Inf[iact_proj].tilt,
                              prm.IMG_Inf[iact_proj].psi);

        project_GridVolume(residual_vol, prm.basis, dummy_proj,
                           read_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                           prm.IMG_Inf[iact_proj].rot,
                           prm.IMG_Inf[iact_proj].tilt,
                           prm.IMG_Inf[iact_proj].psi, BACKWARD, prm.eq_mode,
                           prm.GVNeq, NULL, NULL, prm.ray_length, prm.threads);

        if (!(prm.tell&TELL_SHOW_ERROR))
            if (iact_proj % XMIPP_MAX(1, prm.numIMG / 60) == 0)
                progress_bar(iact_proj);
    }
    if (!(prm.tell&TELL_SHOW_ERROR))
        progress_bar(prm.numIMG);

    // Convert to voxels: solely for output of power of residual volume
    Image<double>      residual_vox;
    int Xoutput_volume_size = (prm.Xoutput_volume_size == 0) ?
                              prm.projXdim : prm.Xoutput_volume_size;
    int Youtput_volume_size = (prm.Youtput_volume_size == 0) ?
                              prm.projYdim : prm.Youtput_volume_size;
    int Zoutput_volume_size = (prm.Zoutput_volume_size == 0) ?
                              prm.projXdim : prm.Zoutput_volume_size;
    prm.basis.changeToVoxels(residual_vol, &(residual_vox()),
                             Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
    pow_residual_vol = residual_vox().sum2() / MULTIDIM_SIZE(residual_vox());
    residual_vox.clear();

    std::cout << "Projection of residual volume; kappa = " << kappa << std::endl;
    if (!(prm.tell&TELL_SHOW_ERROR))
        init_progress_bar(prm.numIMG);

    // Now that we have the residual volume: project in all directions
    pow_residual_imgs = 0.;
    new_proj().resize(read_proj());
    mask.resize(read_proj());
    BinaryCircularMask(mask, YSIZE(read_proj()) / 2, INNER_MASK);

    dim2 = (double)YSIZE(read_proj()) * XSIZE(read_proj());
    for (int iact_proj = 0; iact_proj < prm.numIMG ; iact_proj++)
    {
        project_GridVolume(residual_vol, prm.basis, new_proj,
                           dummy_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                           prm.IMG_Inf[iact_proj].rot,
                           prm.IMG_Inf[iact_proj].tilt,
                           prm.IMG_Inf[iact_proj].psi, FORWARD, prm.eq_mode,
                           prm.GVNeq, A, NULL, prm.ray_length, prm.threads);

        sqrtweight = sqrt(prm.residual_imgs[iact_proj].weight() / prm.sum_weight);

        // Next lines like normalization in [EHL] (2.18)?
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(new_proj())
        {
            dAij(dummy_proj(), i, j) = XMIPP_MAX(1., dAij(dummy_proj(), i, j)); // to avoid division by zero
            dAij(new_proj(), i, j) /= dAij(dummy_proj(), i, j);
        }
        new_proj() *= sqrtweight * kappa;

        /*
        fn_tmp="residual_"+integerToString(iact_proj);
        dummy_proj()=1000*prm.residual_imgs[iact_proj]();
        dummy_proj.write(fn_tmp+".old");
        */

        prm.residual_imgs[iact_proj]() -= new_proj();
        pow_residual_imgs += prm.residual_imgs[iact_proj]().sum2();

        // Mask out edges of the images
        apply_binary_mask(mask, prm.residual_imgs[iact_proj](), prm.residual_imgs[iact_proj](), 0.);

        /*
        dummy_proj()=1000*new_proj();
        dummy_proj.write(fn_tmp+".change");
        dummy_proj()=1000*prm.residual_imgs[iact_proj]();
        dummy_proj.write(fn_tmp+".new");
        */

        if (!(prm.tell&TELL_SHOW_ERROR))
            if (iact_proj % XMIPP_MAX(1, prm.numIMG / 60) == 0)
                progress_bar(iact_proj);
    }

    pow_residual_imgs /= dim2;
    newres_imgs.clear();

    if (!(prm.tell&TELL_SHOW_ERROR))
        progress_bar(prm.numIMG);

}

/* ------------------------------------------------------------------------- */
VariabilityClass::VariabilityClass(BasicARTParameters *_prm,
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
    Image<double> vol_voxels;

    // Convert from basis to voxels
    prm->basis.changeToVoxels(*ptr_vol_out, &(vol_voxels()),
                              Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
    (*ptr_vol_out).initZeros();
    N++;

    // Make the DWT
    MultidimArray<double> DWTV;
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
    DWTV.selfWindow(z1, y1, x1, z2, y2, x2);
    STARTINGZ(DWTV) = STARTINGY(DWTV) = STARTINGX(DWTV) = 0;
    VA.push_back(DWTV);
}
#undef DEBUG

//#define DEBUG
void VariabilityClass::finishAnalysis()
{
    if (VA.size() == 0)
        return;

    // Coocurrence matrix
    int nmax = VA.size();
    MultidimArray<int> coocurrence(nmax, nmax);

    // Study each voxel
    Image<double> SignificantT2, SignificantMaxRatio, SignificantMinRatio;
    int zsize = ZSIZE(VA[0]) / 2;
    int ysize = YSIZE(VA[0]) / 2;
    int xsize = XSIZE(VA[0]) / 2;
    SignificantT2().initZeros(zsize, ysize, xsize);
    SignificantMaxRatio().initZeros(zsize, ysize, xsize);
    SignificantMinRatio().initZeros(zsize, ysize, xsize);
    std::cout << "Classifying voxels ...\n";
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

    FOR_ALL_ELEMENTS_IN_ARRAY3D(SignificantT2())
    {
        // Save the data for this voxel
        std::ofstream fh_dat;
        fh_dat.open("PPP.dat");
        if (!fh_dat)
            REPORT_ERROR(ERR_IO_NOTOPEN, "VariabilityClass::finishAnalysis: "
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
        REPORT_ERROR(ERR_IO_NOTOPEN,(std::string)"VariabilityClass::finishAnalysis: " \
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
        for (size_t n = 0; n < idx0.size(); n++)
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
        for (size_t n = 0; n < idx1.size(); n++)
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
                         ((N0 - 1.0) * covariance0 + (N1 - 1.0) * covariance1);
            covariance *= (1.0 / N0 + 1.0 / N1);
            aux.fromVector(avg_diff);
            T2 = (double)(N0 + N1 - avg_diff.size() - 1) /
                 ((N0 + N1 - 2) * avg_diff.size()) *
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
        for (size_t n = 0; n < idx0.size(); n++)
            for (size_t np = n + 1; np < idx0.size(); np++)
                if (idx0(n) && idx0(np))
                    coocurrence(n, np)++;

        for (size_t n = 0; n < idx1.size(); n++)
            for (size_t np = n + 1; np < idx1.size(); np++)
                if (idx1(n) && idx1(np))
                    coocurrence(n, np)++;

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

        if (++counter % nxny == 0)
            progress_bar(counter);
    }
    progress_bar(MULTIDIM_SIZE(SignificantT2()));
    SignificantT2.write(prm->fn_root + "_variability.vol");
    SignificantMinRatio.write(prm->fn_root + "_minratio.vol");
    SignificantMaxRatio.write(prm->fn_root + "_maxratio.vol");
    system("rm PPP.dat PPP.cod PPP.vs PPP.err PPP.his PPP.inf PPP.0 PPP.1");

    for (int n = 0; n < nmax; n++)
        for (int np = n + 1; np < nmax; np++)
            coocurrence(np, n) = coocurrence(n, np);
    Image<double> save;
    typeCast(coocurrence, save());
    save.write(prm->fn_root + "_coocurrence.xmp");
}
#undef DEBUG

/* ------------------------------------------------------------------------- */
#define POCS_N_measure  8
#define POCS_N_use      2
/* Constructor ............................................................. */
POCSClass::POCSClass(BasicARTParameters *_prm,
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
    Image<double> vol_POCS, theo_POCS_vol, corr_POCS_vol, vol_voxels;

    if (apply_POCS && POCS_i % POCS_freq == 0)
    {
        Image<double> vol_aux;
        Image<double> *desired_volume = NULL;

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
            symmetrizeVolume(prm->SL, vol_voxels(), vol_aux());
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
            Histogram1D hist;
            MultidimArray<int> aux_mask;
            aux_mask.resize(vol_POCS());
            FOR_ALL_ELEMENTS_IN_ARRAY3D(aux_mask)
            aux_mask(k, i, j) = 1 - (int)vol_POCS(k, i, j);
            long mask_voxels = vol_POCS().countThreshold("below", 0.5, 0);
            compute_hist_within_binary_mask(
                aux_mask, vol_voxels(), hist, 300);
            double known_percentage;
            known_percentage = XMIPP_MIN(100, 100 * prm->known_volume / mask_voxels);
            double threshold;
            threshold = hist.percentil(100 - known_percentage);
            FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_voxels())
            if (vol_voxels(k, i, j) < threshold)
                vol_POCS(k, i, j) = 1;
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
        // int bg = (int) vol_POCS().sum();
        // int fg = MULTIDIM_SIZE(vol_POCS()) - bg;
        int relax = 0, posi = 0;
        const MultidimArray<double> &mVolVoxels=vol_voxels();
        const MultidimArray<double> &mVolPOCS=vol_POCS();
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mVolVoxels)
        if (DIRECT_A3D_ELEM(mVolPOCS,k, i, j) == 1 && DIRECT_A3D_ELEM(mVolVoxels,k, i, j) < 0)
        {
        	DIRECT_A3D_ELEM(mVolPOCS,k, i, j) = 0;
            relax++;
        }
        else if (DIRECT_A3D_ELEM(mVolPOCS,k, i, j) == 0 && DIRECT_A3D_ELEM(mVolVoxels,k, i, j) < 0 && prm->positivity)
        {
        	DIRECT_A3D_ELEM(mVolPOCS,k, i, j) = 1;
            posi++;
        }
        // Debugging messages
        //std::cout << "Relaxation/Positivity " << (double)relax/(double)bg << " "
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
                FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_POCS())
                if (vol_POCS(k, i, j) == 1)
                    (*desired_volume)(k, i, j) = 0;
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
                FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_POCS())
                if (vol_POCS(k, i, j))
                    vol_basis(0)(k, i, j) = 0;
            }
            else
            {
                vol_basis(0)().initZeros();
                FOR_ALL_ELEMENTS_IN_ARRAY3D((*desired_volume)())
                vol_basis(0)(k, i, j) = (*desired_volume)(k, i, j);
            }
            POCS_mean_error = -1;
            break;
        default:
        	REPORT_ERROR(ERR_ARG_INCORRECT,"This function cannot work with this basis");
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

                    std::cout << "1: Changing to " << POCS_freq << std::endl;
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
                    { // While not enough uses
                    }
                    else if (POCS_freq < 3)
                    { // increase frequency
                        POCS_freq++;
#ifdef DEBUG_POCS

                        std::cout << "2: Changing to " << POCS_freq << std::endl;
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

                        std::cout << "3: Changing to " << POCS_freq << std::endl;
#endif

                    }
                    else if (POCS_used > 2)
                    {
                        POCS_freq = prm->POCS_freq;
                        // Change status
                        POCS_used = 0;
                        POCS_state = POCS_lowering;
#ifdef DEBUG_POCS

                        std::cout << "Lowering\n";
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
            default:
            	REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown equation mode");
            }
        }
    }
    else
    {
        POCS_i++;
        POCS_mean_error = -1;
    }
}

