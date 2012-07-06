/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2011)
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
#include <math.h>
#include "micrograph_automatic_picking.h"
#include <data/filters.h>
#include <data/rotational_spectrum.h>
#include <reconstruction/denoise.h>
#include <data/xmipp_fft.h>
#include <data/xmipp_filename.h>
#include <algorithm>

//#define DEBUG_PREPARE
//#define DEBUG_BUILDVECTOR
//#define DEBUG_IMG_BUILDVECTOR
//#define DEBUG_CLASSIFY
//#define DEBUG_AUTO
//#define DEBUG_MORE_AUTO

// ==========================================================================
// Section: Piece Preprocessing =============================================
// ==========================================================================
/* Prepare piece ----------------------------------------------------------- */
static pthread_mutex_t preparePieceMutex = PTHREAD_MUTEX_INITIALIZER;
bool AutoParticlePicking::prepare_piece(MultidimArray<double> &piece,
                                        MultidimArray<int> &ipiece, MultidimArray<double> &original_piece)
{
    original_piece = piece;

#ifdef DEBUG_PREPARE

    Image<double> save;
    save() = piece;
    save.write("PPPpiece0.xmp");
#endif

    // Denoise the piece
    if (__fast)
        selfScaleToSize(1,piece,XSIZE(piece)/2,YSIZE(piece)/2);
    else
    {
        WaveletFilter denoiser;
        denoiser.denoising_type = WaveletFilter::BAYESIAN;
        denoiser.scale = 3;
        denoiser.output_scale = 1;
        denoiser.produceSideInfo();
        denoiser.apply(piece);
        if (!(piece(0, 0) == piece(0, 0)))
            return false;
    }

#ifdef DEBUG_PREPARE

    save() = piece;
    save.write("PPPpiece1.xmp");
#endif

    // Band pass filter
    pthread_mutex_lock(&preparePieceMutex);
    if (__filter==NULL)
    {
        // Band pass filter
        __filter=new FourierFilter();
        __filter->FilterShape = RAISED_COSINE;
        __filter->FilterBand = BANDPASS;
        __filter->w1 = __highpass_cutoff;
        __filter->w2 = 1.0 / (__particle_radius / (__reduction * 4));
        __filter->raised_w = XMIPP_MIN(0.02, __highpass_cutoff);
        __filter->generateMask(piece);

        // Plus Gaussian filter
        double w1=__filter->w2/2;
        double K2=-0.5/(w1*w1);
        double K1=1/sqrt(2*PI*w1);
        for (int i=0; i<YSIZE(__filter->maskFourierd); i++)
        {
            double wy;
            FFT_IDX2DIGFREQ(i,YSIZE(piece),wy);
            double wy2=wy*wy;
            for (int j=0; j<XSIZE(__filter->maskFourierd); j++)
            {
                double wx;
                FFT_IDX2DIGFREQ(j,XSIZE(piece),wx);
                double w2=wy2+wx*wx;
                DIRECT_A2D_ELEM(__filter->maskFourierd,i,j)*=K1*exp(K2*w2);
            }
        }
    }
    __filter->applyMaskSpace(piece);
    STARTINGX(piece) = STARTINGY(piece) = 0;
#ifdef DEBUG_PREPARE

    save() = piece;
    save.write("PPPpiece2.xmp");
#endif

    pthread_mutex_unlock(&preparePieceMutex);

    // Reject 2% of the outliers
    reject_outliers(piece, 2.0);

#ifdef DEBUG_PREPARE

    save() = piece;
    save.write("PPPpiece2_5.xmp");
#endif

    // Equalize histogram
    histogram_equalization(piece, __gray_bins);

    typeCast(piece, ipiece);

#ifdef DEBUG_PREPARE

    save() = piece;
    save.write("PPPpiece3.xmp");
#endif

    if (!original_piece.sameShape(piece))
        selfScaleToSize(LINEAR, original_piece, YSIZE(piece), XSIZE(piece));

#ifdef DEBUG_PREPARE

    save() = original_piece;
    save.write("PPPpiece4.xmp");
    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
#endif

    return true;
}

// ==========================================================================
// Section: Feature vector ==================================================
// ==========================================================================
/* Build classification vector --------------------------------------------- */
bool AutoParticlePicking::build_vector(const MultidimArray<int> &piece,
                                       const MultidimArray<double> &original_piece, int _x, int _y,
                                       Matrix1D<double> &_result)
{
#ifdef DEBUG_BUILDVECTOR
    std::cout << "build_vector(" << _x << "," << _y << "," << "_result)" << std::endl;
#endif

    // Resize the output and make same aliases
    int angleBins = __sector.size();
    _result.initZeros(32 // Histogram of the original piece
                      + __radial_bins * (__gray_bins - 1) // radial histograms
                      + (angleBins - 1) + (angleBins - 1) * angleBins // sector correlations
                      + (2 * __radial_bins - 19) * angleBins // ring correlations
                     );
    const MultidimArray<int> &mask = __mask.get_binary_mask();
    const MultidimArray<int> &classif1 = (*(__mask_classification[0]));
    const MultidimArray<int> &classif2 = (*(__mask_classification[1]));

    if (STARTINGX(mask) + _x < STARTINGX(piece))
        return false;
    if (STARTINGY(mask) + _y < STARTINGY(piece))
        return false;
    if (FINISHINGX(mask) + _x > FINISHINGX(piece))
        return false;
    if (FINISHINGY(mask) + _y > FINISHINGY(piece))
        return false;

#ifdef DEBUG_IMG_BUILDVECTOR

    bool debug_go = false;
    Image<double> save, savefg, saveOrig;
    if (true)
    {
        typeCast(piece,save());
        save.write("PPP0.xmp");
        save().initZeros(YSIZE(mask), XSIZE(mask));
        STARTINGY(save()) = STARTINGY(mask);
        STARTINGX(save()) = STARTINGX(mask);
        savefg() = save();
        savefg().initConstant(-1);
        saveOrig() = savefg();
        debug_go = true;
    }
#endif

    MultidimArray<int> radial_idx(__radial_bins);

    // Copy __radial_val, __sector, and __ring into local structures
    // so that threads are applicable
    std::vector<MultidimArray<int> > radial_val;
    for (int i = 0; i < __radial_val.size(); i++)
    {
        MultidimArray<int> dummyi = *(__radial_val[i]);
        radial_val.push_back(dummyi);
    }

    std::vector<MultidimArray<double> > sector;
    MultidimArray<double> &dummyd = *(__sector[0]);
    for (int i = 0; i < __sector.size(); i++)
        sector.push_back(dummyd);

    std::vector<MultidimArray<double> > ring;
    dummyd = *(__ring[0]);
    for (int i = 0; i < __ring.size(); i++)
        ring.push_back(dummyd);

    Histogram1D histogramOriginal;
    histogramOriginal.init(0, 255, 32);

    // Put the image values into the corresponding radial bins
    FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
    {
        int idx1 = A2D_ELEM(classif1, i, j);
        if (idx1 != -1)
        {
            int val = A2D_ELEM(piece, _y + i, _x + j);
            double orig_val = A2D_ELEM(original_piece, _y + i, _x + j);
            INSERT_VALUE(histogramOriginal,orig_val);

            MultidimArray<int> &auxi = radial_val[idx1];
            int aux = A1D_ELEM(radial_idx,idx1);
            A1D_ELEM(auxi,aux) = val;
            A1D_ELEM(radial_idx,idx1)++;
            int idx2 = A2D_ELEM(classif2, i, j);
            if (idx2 != -1)
            {
                MultidimArray<double> &auxd = sector[idx2];
                A1D_ELEM(auxd,idx1) += val;
            }
        }

        // Get particle
#ifdef DEBUG_IMG_BUILDVECTOR
        if (debug_go)
        {
            int val = A2D_ELEM(piece, _y + i, _x + j);
            save(i, j) = val;
            if (A2D_ELEM(mask, i, j))
            {
                savefg(i, j) = val;
                saveOrig(i, j) = original_piece(_y+i, _x+j);
            }
        }
#endif

    }

    // Compute the sector averages and reorganize the data in rings
    for (int j = 0; j < angleBins; j++)
    {
        MultidimArray<int> &Nsector_j = *(__Nsector[j]);
        MultidimArray<double> &sector_j = sector[j];
        FOR_ALL_ELEMENTS_IN_ARRAY1D(sector_j)
        {
            if (DIRECT_A1D_ELEM(Nsector_j,i) > 0)
            {
                DIRECT_A1D_ELEM(sector_j,i) /= DIRECT_A1D_ELEM(Nsector_j,i);
                MultidimArray<double> &ring_i = ring[i];
                DIRECT_A1D_ELEM(ring_i,j) = DIRECT_A1D_ELEM(sector_j,i);
            }
        }
    }

    // Store the histogram of
    int idx_result = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(histogramOriginal)
    VEC_ELEM(_result,idx_result++) = DIRECT_A1D_ELEM(histogramOriginal,i);

    // Compute the histogram of the radial bins and store them
    Histogram1D hist;
    for (int i = 0; i < __radial_bins; i++)
    {
        compute_hist(radial_val[i], hist, 0, __gray_bins - 1, __gray_bins);
        for (int j = 0; j < __gray_bins - 1; j++)
            VEC_ELEM(_result,idx_result++) = DIRECT_A1D_ELEM(hist,j);
    }

    // Compute the correlation of the sectors
    MultidimArray<double> sectorCorr, sectorAutocorr;
    for (int step = 1; step < angleBins; step++)
    {
        sectorCorr.initZeros(angleBins);
        for (int i = 0; i < angleBins; i++)
        {
            DIRECT_A1D_ELEM(sectorCorr,i) = correlationIndex(sector[i],
                                            sector[intWRAP(i+step,0,angleBins-1)]);
        }
        VEC_ELEM(_result,idx_result++) = sectorCorr.computeAvg();
        correlation_vector_no_Fourier(sectorCorr, sectorCorr, sectorAutocorr);
        for (int j = 0; j < XSIZE(sectorAutocorr); j++)
            VEC_ELEM(_result,idx_result++) = DIRECT_A1D_ELEM(sectorAutocorr,j);
    }

    // Compute the correlation in rings
    MultidimArray<double> ringCorr;
    for (int step = 1; step <= 2; step++)
        for (int i = 8; i < __radial_bins - step; i++)
        {
            correlation_vector_no_Fourier(ring[i], ring[i + step], ringCorr);
            for (int j = 0; j < XSIZE(ringCorr); j++)
                VEC_ELEM(_result,idx_result++) = DIRECT_A1D_ELEM(ringCorr,j);
        }

#ifdef DEBUG_IMG_BUILDVECTOR
    if (debug_go)
    {
        save.write("PPP1.xmp");
        savefg.write("PPP2.xmp");
        saveOrig.write("PPP3.xmp");

        save().initZeros(savefg());
        FOR_ALL_ELEMENTS_IN_ARRAY2D(save())
        {
            int idx1 = classif1(i, j);
            if (idx1 == -1)
                continue;
            int idx2 = classif2(i, j);
            if (idx2 == -1)
                continue;
            save(i,j)=sector[idx2](idx1);
        }
        save.write("PPP4.xmp");

        save().initZeros();
        FOR_ALL_ELEMENTS_IN_ARRAY2D(save())
        {
            int idx1 = classif1(i, j);
            if (idx1 == -1)
                continue;
            int idx2 = classif2(i, j);
            if (idx2 == -1)
                continue;
            save(i,j)=ring[idx1](idx2);
        }
        save.write("PPP5.xmp");
    }
#endif

#ifdef DEBUG_BUILDVECTOR
    MultidimArray<double> A;
    A.resize(angleBins,XSIZE(sector[0]));
    FOR_ALL_ELEMENTS_IN_ARRAY2D(A)
    A(i,j)=sector[i](j);
    A.write("PPPsector.txt");
    std::cout << _result.transpose() << std::endl;
#endif

#if defined(DEBUG_BUILDVECTOR) || defined(DEBUG_IMG_BUILDVECTOR)

    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
#endif

    return true;
}

// ==========================================================================
// Section: Training ========================================================
// ==========================================================================
/* Learn particles --------------------------------------------------------- */
void AutoParticlePicking::learnParticles(int _particle_radius)
{
    __particle_radius = _particle_radius;
    if (4 * __particle_radius <= 512)
        __piece_xsize = 512;
    else
        __piece_xsize = NEXT_POWER_OF_2(6*__particle_radius);
    __scan_overlap = round(1.8*__particle_radius);
    __piece_overlap = 2*__particle_radius;
    __learn_overlap = __particle_radius;
    createMask(2 * __particle_radius);

    std::cerr << "\n------------------Learning Phase-----------------------\n";
    buildPositiveVectors();
    getAutoTruePositives();
    buildNegativeVectors(false);
    getAutoFalsePositives();
    buildNegativeVectors(true);

    __selection_model.addMicrographScanned(count_scanning_pos());

    std::cerr << "Learning process finished..." << std::endl;
}
/* Get the false-positive particles----------------------------------------- */
void AutoParticlePicking::getAutoFalsePositives()
{
    int Nrejected = 0;
    int imax = __rejected_particles.size();
    __selection_model.addFalsePositives(imax);
    for (int i = 0; i < imax; i++)
        if (__selection_model.addParticleTraining(__rejected_particles[i], 2))
            ++Nrejected;
    std::cout << Nrejected << " false positives are considered\n";
}

/* Get the false-positive particles----------------------------------------- */
void AutoParticlePicking::getAutoTruePositives()
{
    // Add the true positives given by the threshold
    int Naccepted = 0;
    int imax = __auto_candidates.size();
    for (int i = 0; i < imax; i++)
    {
    	const Particle &p=__auto_candidates[i];
    	if (p.status==1)
        	if (__selection_model.addParticleTraining(p, 0))
        		Naccepted++;
    }
    std::cout << Naccepted << " true positives are considered\n";
}

/* Build training vectors ---------------------------------------------------*/
void AutoParticlePicking::buildPositiveVectors()
{
    int num_part = __m->ParticleNo();
    std::cout << "Building " << num_part
    << " particle vectors for this image. Please wait..." << std::endl;
    __selection_model.addMicrograph();
    int width, height;
    __m->size(width, height);
    Matrix1D<double> v;
    int numParticles = 0;

    MultidimArray<char> visited(num_part);
    MultidimArray<double> piece, original_piece;
    MultidimArray<int> ipiece;
    Particle p;
    p.status = 1;
    p.cost = 1.0;
    p.micrograph=__fn_micrograph;
    while (visited.sum() < num_part)
    {
        // get the first un-visited particle in the array
        int part_i = 0;
        while (visited(part_i) == 1)
            part_i++;
        visited(part_i) = 1;

        // Get the piece containing that particle
        int x = __m->coord(part_i).X;
        int y = __m->coord(part_i).Y;
        int posx, posy;
        get_centered_piece(piece, x, y, posx, posy);

        // Denoise, reduce, reject outliers and equalize histogram
        bool success = prepare_piece(piece, ipiece, original_piece);

        if (!success)
            continue;
        posx = ROUND(posx / __reduction);
        posy = ROUND(posy / __reduction);

        //make vector from this particle
        success = build_vector(ipiece, original_piece, posx, posy, v);
        if (success)
        {
            p.x = x;
            p.y = y;
            p.vec = v;
            if (__selection_model.addParticleTraining(p, 0))
                numParticles++;
        }

        // make vector from the neighbours
        std::vector<Matrix1D<int> > neighbour;
        neighbour.reserve(num_part);

        find_neighbour(piece, part_i, x, y, posx, posy, visited, neighbour);

        for (int i = 0; i < neighbour.size(); i++)
        {
            part_i = neighbour.at(i)(0);
            posx = neighbour.at(i)(1);
            posy = neighbour.at(i)(2);
            success = build_vector(ipiece, original_piece, posx, posy, v);
            visited(part_i) = 1;
            if (success)
            {
                p.x = __m->coord(part_i).X;
                p.y = __m->coord(part_i).Y;
                p.vec = v;
                p.status = 1;
                p.cost = 1.0;
                if (__selection_model.addParticleTraining(p, 0))
                    numParticles++;
            }
        }
    }
    __selection_model.addParticlePicked(numParticles);
}

/* Build vector from non particles------------------------------------------ */
void AutoParticlePicking::buildNegativeVectors(bool checkForPalsePostives)
{
    if (checkForPalsePostives)
        std::cerr << "Building automatic false positives ..." << std::endl;
    else
        std::cerr << "Building non particles ..." << std::endl;
    const MultidimArray<int> &mask = __mask.get_binary_mask();

    // Setup a classification model with the already known data
    std::vector<MultidimArray<double> > features;
    Matrix1D<double> probs, aux1, aux2;
    if (checkForPalsePostives)
    {
        // Prepare data to be classified
        getFeatures(features);
        getClassesProbabilities(probs);

        // Initialize classifier
        __selection_model.initNaiveBayesEnsemble(features, probs, 8,
                __penalization, 10, 1, 1, "mm");
    }

    // top,left corner of the piece
    int top = 0, left = 0, next_top = 0, next_left = 0;

    // If the piece available is small then include the scanned part
    // because we need bigger image for denoising but for scanning
    // particles we skip the already scanned part
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    Matrix1D<double> v;

    int N = 1, Nnonparticles = 0, Nfalsepositives = 0;

    // We do not want any overlap for this process,since it is only for
    // counting the non particles and calculating their features. For
    // the process of automatic selecting we will want an overlap so we
    // do not miss any particle.
    MultidimArray<double> piece, original_piece;
    MultidimArray<int> ipiece;
    while (get_corner_piece(piece, top, left, skip_y, next_skip_x, next_skip_y,
                            next_top, next_left, 0, true))
    {
        // Get a piece and prepare it
        if (!prepare_piece(piece, ipiece, original_piece))
        {
            top = next_top;
            left = next_left;
            std::cerr << "bad piece...skipping" << std::endl;
            N++;
            continue;
        }

        // Express the skip values in the reduced image
        skip_x /= __reduction;
        skip_y /= __reduction;

        // Scan this piece
        int posx = 0, next_posx = 0, posy = 0, next_posy = 0;
        next_posx = posx = skip_x + XSIZE(mask) / 2;
        next_posy = posy = skip_y + YSIZE(mask) / 2;

        // We do not want any overlap for this process, since it is only for
        // counting the non particles and calculating their features. For
        // the process of automatic selecting we will want an overlap so we
        // do not miss any particle.
        Particle P;
        P.status = 1;
        P.cost = -1;
        P.micrograph = __fn_micrograph;
        while (get_next_scanning_pos(piece, next_posx, next_posy, skip_x,
                                     skip_y, 0))
        {
            // Check if there is any particle around
            if (!anyParticle(left + posx * __reduction,
                             top + posy * __reduction, XSIZE(mask) * __reduction))
            {
                if (build_vector(ipiece, original_piece, posx, posy, v))
                {
                    // Build the Particle structure
                    P.x = left + posx * __reduction;
                    P.y = top + posy * __reduction;
                    P.vec = v;
                    if (!checkForPalsePostives)
                    {
                        if (__selection_model.addParticleTraining(P, 1))
                            Nnonparticles++;
                    }
                    else
                    {
                        double cost;
                        int votes = __selection_model.isParticle(v, cost, aux1, aux2);
                        if (votes > 5)
                        {
                            if (__selection_model.addParticleTraining(P, 2))
                                Nfalsepositives++;
                        }
                    }
                }
            }
            // Go to next scanning position
            posx = next_posx;
            posy = next_posy;
        }

        // Go to next piece in the micrograph
        top = next_top;
        left = next_left;
        skip_x = next_skip_x;
        skip_y = next_skip_y;
        N++;
    }
    if (checkForPalsePostives)
        std::cout << Nfalsepositives << " false positive chosen\n";
    else
        std::cout << Nnonparticles << " non particles randomly chosen\n";
}

// ==========================================================================
// Section: Automatic Selection =============================================
// ==========================================================================
struct SAscendingParticleSort
{
    bool operator()(const Particle& rpStart, const Particle& rpEnd)
    {
        return rpStart.cost < rpEnd.cost;
    }
};

/* Automatic phase ----------------------------------------------------------*/
int AutoParticlePicking::automaticallySelectParticles()
{
    // Check that there is a valid model
    if (__selection_model.isEmpty())
    {
        std::cerr << "No model has been created." << std::endl;
        return 0;
    }

    std::cerr << "-----------------Automatic Phase--------------------------\n";

    // Initialize some variables
    __auto_candidates.resize(0);
    const MultidimArray<int> &mask = __mask.get_binary_mask();

    // Get the training features and the a priori probabilities
    std::vector<MultidimArray<double> > features;
    Matrix1D<double> probs;
    getFeatures(features);
    getClassesProbabilities(probs);

    // Initialize classifier
    __selection_model.initNaiveBayesEnsemble(features, probs, 8,
            __penalization, 10, 1, 1, "mm");
#ifdef DEBUG_AUTO

    std::cout << "Probabilities of the classes:"
    << probs.transpose() << std::endl;
#endif

    // Automatically select particles with threads
    pthread_t * th_ids = new pthread_t[__numThreads];
    AutomaticallySelectThreadParams * th_args =
        new AutomaticallySelectThreadParams[__numThreads];
    for (int nt = 0; nt < __numThreads; nt++)
    {
        th_args[nt].autoPicking = this;
        th_args[nt].idThread = nt;
        pthread_create(&th_ids[nt], NULL, automaticallySelectParticlesThread,
                       &th_args[nt]);
    }

    // Waiting for threads to finish
    for (int nt = 0; nt < __numThreads; nt++)
        pthread_join(th_ids[nt], NULL);

    // Thread structures are not needed any more
    delete[] th_ids;
    delete[] th_args;

    // Sort particles by cost
    std::sort(__auto_candidates.begin(), __auto_candidates.end(),
              SAscendingParticleSort());

#ifdef DEBUG_AUTO

    std::cout << "Number of automatically selected particles = "
    << __auto_candidates.size() << std::endl;
#endif

    // Reject the candidates that are pointing to the same particle
    int Nalive = reject_within_distance(__auto_candidates, __particle_radius,
                                        false);

#ifdef DEBUG_AUTO

    std::cout << "Number of automatically selected particles after distance rejection = "
    << Nalive << std::endl;
#endif

    // Apply a second classifier for classifying between particle
    // and false positive. For that, remove the middle class (background)
    if (features.size() == 3)
    {
    	Matrix1D<double> aux1, aux2;
        Classification_model selection_model2(2);
        int imax;
        if (Nalive > 0)
        {
            std::vector<MultidimArray<double> >::iterator featuresIterator =
                features.begin();
            featuresIterator++;
            features.erase(featuresIterator);
            probs(1) = probs(2);
            probs.resize(2);
            probs /= probs.sum();
            selection_model2.initNaiveBayesEnsemble(features, probs, 8,
                                                    __penalization, 10, 1, 1, "mm");
#ifdef DEBUG_AUTO

            std::cout << "Second classification\n";
#endif

            imax = __auto_candidates.size();
            Nalive = 0;
            for (int i = 0; i < imax; i++)
                if (__auto_candidates[i].status == 1)
                {
                    double p;
                    int votes = selection_model2.isParticle(__auto_candidates[i].vec, p, aux1, aux2);
                    if (votes < 8)
                    {
                        __auto_candidates[i].status = 0;
                        __auto_candidates[i].cost = -1;
#ifdef DEBUG_AUTO

                        std::cout << __auto_candidates[i].x << ", "
                        << __auto_candidates[i].y
                        << " is considered as a false positive\n";
#endif

                    }
                    else
                    {
                    	__auto_candidates[i].cost = p;
                        Nalive++;
                    }
                }
        }

        // Apply a third classifier to distinguish between particles and
        // very tough particles
        if (Nalive > 0)
        {
            Classification_model selection_model3(2);
            // Remove from the error class, all those errors that
            // the previous classifier was able to classify correctly
            std::vector<int> toKeep;
            const MultidimArray<double> &features_1 = features[1]; // These are the features of the difficult particles in the model
            Matrix1D<double> trialFeatures(XSIZE(features_1));
            size_t rowLength = XSIZE(features_1) * sizeof(double);
            double cost;
            for (int i = 0; i < YSIZE(features_1); i++)
            {
                memcpy(&VEC_ELEM(trialFeatures,0), &A2D_ELEM(features_1,i,0),rowLength);
                int votes = selection_model2.isParticle(trialFeatures, cost, aux1, aux2);
                if (votes < 8)
                    toKeep.push_back(i);
            }

            if (toKeep.size() > 0)
            {
                MultidimArray<double> difficultParticles;
                difficultParticles.initZeros(toKeep.size(), XSIZE(features_1));
                for (int i=0; i<YSIZE(difficultParticles); ++i)
                {
                	int correspondingInFeatures=toKeep[i];
                	memcpy(&A2D_ELEM(difficultParticles,i, 0),&A2D_ELEM(features_1,correspondingInFeatures,0),rowLength);
                }
                features.pop_back();
                features.push_back(difficultParticles);

                selection_model3.initNaiveBayesEnsemble(features, probs, 8,
                                                        __penalization, 10, 1, 1, "mm");
#ifdef DEBUG_AUTO
                int NErrors = YSIZE(features_1);
                std::cout << "Third classification: " << YSIZE(difficultParticles)
                << " difficult particles out of " << NErrors << "\n"
                << "Before filtering there were " << Nalive << " particles\n";
#endif

                imax = __auto_candidates.size();
                Nalive = 0;
                for (int i = 0; i < imax; i++)
                    if (__auto_candidates[i].status == 1)
                    {
                        double p;
                        int votes = selection_model3.isParticle(__auto_candidates[i].vec, p, aux1, aux2);
                        if (votes < 8)
                        {
                            __auto_candidates[i].status = 0;
                            __auto_candidates[i].cost = -1;
#ifdef DEBUG_AUTO

                            std::cout << __auto_candidates[i].x << ", "
                            << __auto_candidates[i].y
                            << " is considered as a false positive 2\n";
#endif

                        }
                        else
                        {
                        	__auto_candidates[i].cost=p;
                        	Nalive++;
                        }
                    }
            }
        }
    }

    if (Nalive > 0)
    {
        int imax = __auto_candidates.size();
        // Get the maximum and minimum cost
        bool first = true;
        double minCost, maxCost;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].status == 1)
            {
            	double cost=__auto_candidates[i].cost;
                if (!ISINF(cost) && !ISNAN(cost))
                {
                    if (first || cost < minCost)
                    {
                        minCost = cost;
                        first = false;
                    }
                    if (first || cost > maxCost)
                    {
                        maxCost = cost;
                        first = false;
                    }
                }
            }

        // Insert selected particles in the result
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].status == 1)
            {
                if (ISINF(__auto_candidates[i].cost) || ISNAN(__auto_candidates[i].cost))
                    __auto_candidates[i].cost = 1;
                else
                    __auto_candidates[i].cost = (__auto_candidates[i].cost
                                                 - maxCost) / (minCost - maxCost);
#ifdef DEBUG_AUTO

                std::cout << "Particle coords " << __auto_candidates[i].x << ", "
                << __auto_candidates[i].y << " cost= " << __auto_candidates[i].cost << std::endl;
#endif
            }
            else
                __auto_candidates[i].status = 0;
    }

    std::cout << "\nAutomatic process finished. Number of particles found: "
    << Nalive << std::endl;
    return Nalive;
}

static pthread_mutex_t particleAdditionMutex = PTHREAD_MUTEX_INITIALIZER;

void * automaticallySelectParticlesThread(void * args)
{
    // Pick the input parameters
    AutomaticallySelectThreadParams * prm =
        (AutomaticallySelectThreadParams *) args;
    AutoParticlePicking *autoPicking = prm->autoPicking;
    int idThread = prm->idThread;

    int thVotes = 1;
    if (autoPicking->__selection_model.__training_particles.size() == 2
        || (autoPicking->__selection_model.__training_particles.size() == 3
            && autoPicking->__selection_model.__training_particles[2].size()
            == 0))
        thVotes = 8;

    //top,left corner of the piece
    int top = 0, left = 0, next_top = 0, next_left = 0;
    // If the piece available is small then include the scanned part
    // because we need bigger image for denoising but for scanning
    // particles we skip the already scanned part
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    Matrix1D<double> v,aux1,aux2;
    int N = 1, particle_idx = 0, Nscanned = 0;
    MultidimArray<double> piece, original_piece;
    MultidimArray<int> ipiece;
    const MultidimArray<int> &mask = autoPicking->__mask.get_binary_mask();
    std::vector<Particle> threadCandidates;
    Particle p;
    p.micrograph=autoPicking->__fn_micrograph;
    p.status = 1;
    do
    {
        bool isMine = (N % autoPicking->__numThreads == idThread);
        bool pieceOK = autoPicking->get_corner_piece(piece, top, left, skip_y,
                       next_skip_x, next_skip_y, next_top, next_left,
                       autoPicking->__piece_overlap, isMine);
        if (!pieceOK)
            break;
        if (isMine)
        {
            if (idThread == 0)
                std::cerr << ".";
#ifdef DEBUG_MORE_AUTO

            std::cout << "thread " << idThread
            << " processing piece " << N << "...\n";
            std::cout << "    (top,left)=" << top << "," << left
            << " skip y,x=" << next_skip_y << "," << next_skip_x
            << " next=" << next_top << "," << next_left << std::endl;
#endif

            // Get a piece and prepare it
            if (!autoPicking->prepare_piece(piece, ipiece, original_piece))
            {
                top = next_top;
                left = next_left;
                std::cerr << "bad piece...skipping" << std::endl;
                N++;
                continue;
            }

            // Express the skip values in the reduced image
            skip_x /= autoPicking->__reduction;
            skip_y /= autoPicking->__reduction;
#ifdef DEBUG_MORE_AUTO

            std::cout << "Skip(y,x)=" << skip_y << "," << skip_x << std::endl;
#endif

            // Scan this piece
            int posx = 0, next_posx = 0, posy = 0, next_posy = 0;
            next_posx = posx = skip_x + XSIZE(mask) / 2;
            next_posy = posy = skip_y + YSIZE(mask) / 2;

            while (autoPicking->get_next_scanning_pos(piece, next_posx,
                    next_posy, skip_x, skip_y, autoPicking->__scan_overlap))
            {
                // COSS: Uncomment the next sentence for fast debugging
                // if (rnd_unif(0,1)<0.98) continue;

#ifdef DEBUG_MORE_AUTO
                std::cout << "Pos(y,x)=" << posy << "," << posx
                << " Micro(y,x)=" << posy*autoPicking->__reduction + top
                << "," << posx*autoPicking->__reduction + left
                << " Next pos(y,x)=" << next_posy << "," << next_posx
                << std::endl;
#endif

                if (autoPicking->build_vector(ipiece, original_piece, posx,
                                              posy, v))
                {
                    double cost;
                    int votes = autoPicking->__selection_model.isParticle(v,cost,aux1,aux2);
#ifdef DEBUG_MORE_AUTO
                    std::cout << "votes= " << votes << " cost= " << cost << std::endl;
#endif
                    if (votes > thVotes)
                    {
#ifdef DEBUG_MORE_AUTO
                        std::cout << "Particle Found: "
                        << left + posx * autoPicking->__reduction << ","
                        << top + posy *autoPicking->__reduction
                        << " votes=" << votes << " cost=" << cost
                        << std::endl;
                        std::cout << "Press any key to continue...\n";
                        char c;
                        std::cin >> c;
#endif

                        // Build the Particle structure
                        p.x = left + posx * autoPicking->__reduction;
                        p.y = top + posy * autoPicking->__reduction;
                        p.vec = v;
                        p.cost = cost;
                        threadCandidates.push_back(p);
                    }
                    Nscanned++;
                }

                // Go to next scanning position
                posx = next_posx;
                posy = next_posy;
            }
        }

        // Go to next piece in the micrograph
        top = next_top;
        left = next_left;
        skip_x = next_skip_x;
        skip_y = next_skip_y;
        N++;
    }
    while (true);

    int imax = threadCandidates.size();
    pthread_mutex_lock(&particleAdditionMutex);
    for (int i = 0; i < imax; i++)
        autoPicking->__auto_candidates.push_back(threadCandidates[i]);
    pthread_mutex_unlock(&particleAdditionMutex);
}

/* Filter particles --------------------------------------------------------*/
//To calculate the euclidean distance between to points
double euclidean_distance(const Particle &p1, const Particle &p2)
{
    double dx = (p1.x - p2.x);
    double dy = (p1.y - p2.y);
    return sqrt(dx * dx + dy * dy);
}

int AutoParticlePicking::reject_within_distance(std::vector<Particle> &_Input,
        double _min_dist, bool _reject_both)
{
    int imax = _Input.size();
    int n = 0;
    for (int i = 0; i < imax; i++)
    {
        if (_Input.at(i).status == 0)
            continue;
        for (int j = i + 1; j < imax; j++)
        {
            if (_Input.at(j).status == 0)
                continue;
            double dist = euclidean_distance(_Input.at(i), _Input.at(j));
            if (dist < _min_dist)
            {
                _Input.at(j).status = 0;
                _Input.at(j).cost = -1;
                if (_reject_both)
                {
                    _Input.at(i).status = 0;
                    _Input.at(i).cost = -1;
                }
            }
        }
        if (_Input.at(i).status == 1)
            n++;
    }
    return n;
}

// ==========================================================================
// Section: Mask for feature construction ===================================
// ==========================================================================
/* Create mask for learning particles ---------------------------------------*/
void AutoParticlePicking::createMask(int mask_size)
{
    if (XSIZE(__mask.get_binary_mask()) != 0)
        return;

    int xsize = mask_size / __reduction;
    int ysize = xsize;
    int radius = __particle_radius / __reduction;

    __mask.type = BINARY_CIRCULAR_MASK;
    __mask.mode = INNER_MASK;
    __mask.R1 = radius;
    __mask.resize(xsize, ysize);
    __mask.generate_mask();
    __mask.get_binary_mask().setXmippOrigin();

    classifyMask();
}

/* Classify the mask pixels -------------------------------------------------*/
void AutoParticlePicking::classifyMask()
{
    const MultidimArray<int> &mask = __mask.get_binary_mask();
    if (XSIZE(mask) == 0)
        return;
    double max_radius_particle = __particle_radius / __reduction;

    // Determine max_radius
    double max_radius = max_radius_particle;

    // Initialize some variables
    // 6 is the minimum radius to be informative
    double radial_step = (max_radius - 6) / __radial_bins;

    MultidimArray<int> *classif1 = new MultidimArray<int> ;
    classif1->resize(mask);
    classif1->initConstant(-1);

    MultidimArray<int> *classif2 = new MultidimArray<int> ;
    classif2->resize(mask);
    classif2->initConstant(-1);
    const double deltaAng = (PI / 8.0);

    MultidimArray<int> Nrad(__radial_bins);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
    {
        double radius = sqrt((double) (i * i + j * j));
        double angle = atan2((double) i, (double) j);
        if (angle < 0)
            angle += 2 * PI;

        if (radius < max_radius)
        {
            // Classif1 is the classification for the radial mass distribution
            int radius_idx;
            if (radius > 6)
                radius_idx = XMIPP_MIN(__radial_bins - 1, 1 +
                                       FLOOR((radius - 6) / radial_step));
            else
                radius_idx = 0;
            (*classif1)(i, j) = radius_idx;
            Nrad(radius_idx)++;

            // Classif2 is the classification by angles
            (*classif2)(i, j) = FLOOR(angle/deltaAng);
        }
    }
    __mask_classification.push_back(classif1);
    __mask_classification.push_back(classif2);

    // Create the holders for the radius values in classif1
    for (int i = 0; i < __radial_bins; i++)
    {
        MultidimArray<int> *aux1 = new MultidimArray<int> ;
        aux1->initZeros(Nrad(i));
        __radial_val.push_back(aux1);
    }

    // Create the holders for the radius values in classif1 and classif2
    int angleBins = classif2->computeMax() + 1;
    for (int i = 0; i < angleBins; i++)
    {
        MultidimArray<double> *aux2 = new MultidimArray<double> ;
        aux2->initZeros(__radial_bins);
        __sector.push_back(aux2);

        MultidimArray<int> *iaux2 = new MultidimArray<int> ;
        iaux2->initZeros(__radial_bins);
        __Nsector.push_back(iaux2);
    }

    // Compute the number of elements in each sector
    FOR_ALL_ELEMENTS_IN_ARRAY2D(*classif1)
    {
        int idx1 = (*classif1)(i, j);
        if (idx1 != -1)
        {
            int idx2 = (*classif2)(i, j);
            if (idx2 != -1)
                (*__Nsector[idx2])(idx1)++;
        }
    }

    for (int i = 0; i < __radial_bins; i++)
    {
        MultidimArray<double> *aux3 = new MultidimArray<double> ;
        aux3->initZeros(angleBins);
        __ring.push_back(aux3);
    }

#ifdef DEBUG_CLASSIFY
    Image<double> save;
    typeCast(*classif1, save());
    save.write("PPPmaskClassification1.xmp");

    typeCast(*classif2, save());
    save.write("PPPmaskClassification2.xmp");
#endif
}

// ==========================================================================
// Section: Scanning and pieces =============================================
// ==========================================================================
/* Count the number of points scanned in a micrograph --------------------- */
int AutoParticlePicking::count_scanning_pos() const
{
    const MultidimArray<int> &mask = __mask.get_binary_mask();
    int top = 0, left = 0, next_top = 0, next_left = 0;
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    int Nscanned = 0;
    MultidimArray<double> piece;
    while (get_corner_piece(piece, top, left, skip_y, next_skip_x, next_skip_y,
                            next_top, next_left, __piece_overlap, true))
    {
        top = next_top;
        left = next_left;
        skip_x /= __reduction;
        skip_y /= __reduction;
        int posx = 0, next_posx = 0, posy = 0, next_posy = 0;
        next_posx = posx = skip_x + XSIZE(mask) / 2;
        next_posy = posy = skip_y + YSIZE(mask) / 2;

        while (get_next_scanning_pos(piece, next_posx, next_posy, skip_x,
                                     skip_y, __scan_overlap))
        {
            Nscanned++;
            posx = next_posx;
            posy = next_posy;
        }

        top = next_top;
        left = next_left;
        skip_x = next_skip_x;
        skip_y = next_skip_y;
    }
    return Nscanned;
}

/* Get next scanning position ---------------------------------------------- */
bool AutoParticlePicking::get_next_scanning_pos(
    const MultidimArray<double> &piece, int &_x, int &_y, int _skip_x,
    int _skip_y, int overlap) const
{
    const MultidimArray<int> &mask = __mask.get_binary_mask();

    if (_x + XSIZE(mask) / 2 > XSIZE(piece) || _y + YSIZE(mask) / 2
        > YSIZE(piece))
        return false;

    if (_x == 0 && _y == 0)
    {
        _x += XSIZE(mask) - overlap / __reduction;
        _y += YSIZE(mask) - overlap / __reduction;
    }
    else
    {
        int nextx = _x + XSIZE(mask) - overlap / __reduction;
        int nexty = _y + YSIZE(mask) - overlap / __reduction;

        if (nextx + (XSIZE(mask) / 2) > XSIZE(piece))
        {
            if (nexty + (YSIZE(mask) / 2) > YSIZE(piece))
            {
                _x = nextx;
                _y = nexty;
            }
            else
            {
                _x = _skip_x + (XSIZE(mask) / 2);
                _y = nexty;
            }
        }
        else
        {
            _x = nextx;
        }
    }
    return true;
}

/* Get neighbours ---------------------------------------------------------- */
//To get the neighbours and their positions in the piece image
void AutoParticlePicking::find_neighbour(const MultidimArray<double> &piece,
        int _index, int _x, int _y, int _posx, int _posy,
        MultidimArray<char> &_visited, std::vector<Matrix1D<int> > &_neighbours)
{
    int piece_xsize = XSIZE(piece);
    int piece_ysize = YSIZE(piece);
    const MultidimArray<int> &mask = __mask.get_binary_mask();

    // If all the particles are visited
    if (_visited.sum() == XSIZE(_visited))
        return;

    int current_part = _index + 1;
    _neighbours.clear();
    _neighbours.reserve(XSIZE(_visited));

    int top = CEIL((double)_y / __reduction) - _posy;
    int left = CEIL((double)_x / __reduction) - _posx;
    int bottom = top + FLOOR((double)piece_ysize / __reduction) - 1;
    int right = left + FLOOR((double)piece_xsize / __reduction) - 1;
    int xmask2 = CEIL(XSIZE(mask) / 2);
    int ymask2 = CEIL(YSIZE(mask) / 2);

    Matrix1D<int> current_nbr(3);
    while (current_part < XSIZE(_visited))
    {
        //find the next unvisited particle
        while (A1D_ELEM(_visited,current_part))
        {
            current_part++;
            if (current_part == XSIZE(_visited))
                break;
        }
        if (current_part == XSIZE(_visited))
            break;

        //check if it is neighbour or not
        int nx = round((double) __m->coord(current_part).X / __reduction);
        int ny = round((double) __m->coord(current_part).Y / __reduction);
        if ((nx - xmask2 > left) && (ny - ymask2 > top)
            && (nx + xmask2 < right) && (ny + ymask2 < bottom))
        {
            current_nbr(0) = current_part;
            current_nbr(1) = nx - left;
            current_nbr(2) = ny - top;
            _neighbours.push_back(current_nbr);
        }
        current_part++;
    }
}

bool AutoParticlePicking::anyParticle(int posx, int posy, int rect_size)
{
    int num_part = __m->ParticleNo();
    const std::vector<Particle_coords> &selected_particles = __m->Particles();
    for (int i = 0; i < num_part; i++)
    {
        int _x = selected_particles[i].X;
        int _y = selected_particles[i].Y;

        if ((_x > posx - rect_size) && (_x < posx + rect_size))
        {
            if ((_y > posy - rect_size) && (_y < posy + rect_size))
                return true;
        }
    }
    return false;
}

/* Get piece --------------------------------------------------------------- */
// to get the piece containing (x,y) of size xsize,ysize
// return the position of x,y in the piece in posx,posy
void AutoParticlePicking::get_centered_piece(MultidimArray<double> &piece,
        int _x, int _y, int &_posx, int &_posy)
{
    piece.initZeros(__piece_xsize, __piece_xsize);
    int startx = _x - ROUND(__piece_xsize / 2);
    int endx = _x + ROUND(__piece_xsize / 2);
    int starty = _y - ROUND(__piece_xsize / 2);
    int endy = _y + ROUND(__piece_xsize / 2);
    int maxx, maxy;
    __m->size(maxx, maxy);
    _posx = ROUND(__piece_xsize / 2);
    _posy = ROUND(__piece_xsize / 2);

    // boundary adjustments
    if (startx < 0)
    {
        _posx += startx;
        startx = 0;
        endx = __piece_xsize - 1;
    }
    if (starty < 0)
    {
        _posy += starty;
        starty = 0;
        endy = __piece_xsize - 1;
    }
    if (endx > maxx - 1)
    {
        _posx += endx - (maxx - 1);
        endx = maxx - 1;
        startx = endx - __piece_xsize;
    }
    if (endy > maxy - 1)
    {
        _posy += endy - (maxy - 1);
        endy = maxy - 1;
        starty = endy - __piece_xsize;
    }


    //read the matrix from the micrograph
    if (__incore)
    {
        const MultidimArray<double> &micrograph = MULTIDIM_ARRAY(__I);
        size_t rowLength=__piece_xsize*sizeof(double);
        for (int i = 0; i < __piece_xsize; i++)
            memcpy(&DIRECT_A2D_ELEM(piece, i, 0),&DIRECT_A2D_ELEM(micrograph,starty + i, startx),rowLength);
    }
    else
    {
        const Micrograph &micrograph = *__m;
        for (int i = 0; i < __piece_xsize; i++)
            for (int j = 0; j < __piece_xsize; j++)
                DIRECT_A2D_ELEM(piece, i, j) = micrograph(starty + i, startx + j);
    }
}

// Get a piece whose top-left corner is at the desired position (if possible)
bool AutoParticlePicking::get_corner_piece(MultidimArray<double> &piece,
        int _top, int _left, int _skip_y, int &_next_skip_x, int &_next_skip_y,
        int &_next_top, int &_next_left, int overlap, bool copyPiece) const
{
    const MultidimArray<int> &mask = __mask.get_binary_mask();

    int maxx, maxy;
    __m->size(maxx, maxy);

    if (maxx < _left + __piece_xsize || maxy < _top + __piece_xsize)
        return false;

    _next_skip_x = _next_skip_y = 0;
    bool increase_Y = false;
    if (_left + __piece_xsize != maxx)
    {
        _next_left = _left + __piece_xsize - overlap;
        if (_next_left + __piece_xsize >= maxx)
            _next_left = maxx - __piece_xsize;
    }
    else
    {
        _next_left = 0;
        increase_Y = true;
    }
    if (increase_Y)
    {
        if (_top + __piece_xsize != maxy)
        {
            _next_top = _top + __piece_xsize - overlap;
            if (_next_top + __piece_xsize >= maxy)
                _next_top = maxy - __piece_xsize;
        }
        else
        {
            _next_top = maxy;
        }
    }

    //read the matrix from the micrograph
    if (copyPiece)
    {
        piece.resizeNoCopy(__piece_xsize, __piece_xsize);
        if (__incore)
        {
            const MultidimArray<double> &micrograph = MULTIDIM_ARRAY(__I);
            size_t rowLength=__piece_xsize*sizeof(double);
            for (int i = 0; i < __piece_xsize; i++)
                memcpy(&DIRECT_A2D_ELEM(piece, i, 0),&DIRECT_A2D_ELEM(micrograph,_top + i, _left),rowLength);
        }
        else
        {
            const Micrograph& micrograph = *__m;
            for (int i = 0; i < __piece_xsize; i++)
                for (int j = 0; j < __piece_xsize; j++)
                    DIRECT_A2D_ELEM(piece,i, j) = micrograph(_top + i, _left + j);
        }
    }

    return true;
}

// ==========================================================================
// Section: Constructors, Initializers, Gets ================================
// ==========================================================================
/* Initialize -------------------------------------------------------------- */
void Classification_model::initNaiveBayesEnsemble(
    const std::vector<MultidimArray<double> > &features,
    const Matrix1D<double> &probs, int discreteLevels, double penalization,
    int numberOfClassifiers, double samplingFeatures,
    double samplingIndividuals, const String &newJudgeCombination)
{
    __bayesEnsembleNet = new EnsembleNaiveBayes(features, probs,
                         discreteLevels, numberOfClassifiers, samplingFeatures,
                         samplingIndividuals, newJudgeCombination);
    int K = features.size();
    Matrix2D<double> cost(K, K);
    cost.initConstant(1);
    for (int i = 0; i < MAT_XSIZE(cost); i++)
        cost(i, i) = 0;
    cost(0, K - 1) = penalization;
    __bayesEnsembleNet->setCostMatrix(cost);
}

AutoParticlePicking::AutoParticlePicking(const FileName &fn, Micrograph *_m, bool _fast)
{
    __m = _m;
    __fn_micrograph = fn;
    __numThreads = 1;
    __piece_xsize = 512;
    __particle_radius = 0;
    __piece_overlap = 0;
    __scan_overlap = 0;
    __learn_overlap = 0;
    __fast=_fast;
    __incore=false;
    __filter=NULL;
}

AutoParticlePicking::~AutoParticlePicking()
{
    delete __filter;
}

void AutoParticlePicking::readMicrograph()
{
    __I.read(__fn_micrograph);
    __incore=true;
}

/* produceFeatures --------------------------------------------------------- */
void AutoParticlePicking::getFeatures(std::vector<MultidimArray<double> > &_features)
{
    _features.clear();
    std::vector<std::vector<Particle> > particles = __selection_model.__training_particles;
    int vec_size = 0;
    int classNo=particles.size();
    for (int i = 0; i < particles.size(); i++)
        if (particles[i].size() != 0)
        {
            vec_size = (particles[i][0].vec).size();
            break;
        }

    for (int j = 0; j < classNo; j++)
    {
        int imax = particles[j].size();
        //if we do not have false-positives, then we only use 2 classes.
        if (imax == 0)
        {
            _features.resize(classNo - 1);
            break;
        }
        _features.push_back(*(new MultidimArray<double> ));
        MultidimArray<double> &features_j=_features[j];
        const std::vector<Particle> &particles_j=particles[j];
        features_j.resizeNoCopy(imax, vec_size);
        for (int i = 0; i < imax; i++)
            features_j.setRow(i, particles_j[i].vec);
    }
}

/* produceClassesProbabilities --------------------------------------------- */
void AutoParticlePicking::getClassesProbabilities(Matrix1D<double> &probabilities)
{
    double micrographsScanned = 0.0;
    double particlesMarked = __selection_model.__training_particles[0].size();
    double falsePositives = __selection_model.__training_particles[2].size();
    probabilities.initZeros(3);

    for (int i = 0; i < __selection_model.__micrographs_number; i++)
        if (i < __selection_model.__micrographs_scanned.size())
            micrographsScanned += __selection_model.__micrographs_scanned[i];
    micrographsScanned = XMIPP_MAX(micrographsScanned,
                                   2*(particlesMarked+falsePositives));
    int NtrainingVectors = 0;
    for (int i = 0; i < __selection_model.__training_particles.size(); i++)
        NtrainingVectors += __selection_model.__training_particles[i].size();
    micrographsScanned = XMIPP_MAX(micrographsScanned, NtrainingVectors);

    probabilities(0) = particlesMarked / micrographsScanned;
    probabilities(1) = (micrographsScanned - particlesMarked - falsePositives)
                       / micrographsScanned;
    probabilities(2) = falsePositives / micrographsScanned;
}

Classification_model::Classification_model(int _classNo,
        int _maxTrainingVectors)
{
    __maxTrainingVectors = _maxTrainingVectors;
    __classNo = _classNo;
    __training_particles.resize(__classNo);
    __micrographs_number = 0;
    __falsePositives.resize(0);
}

bool Classification_model::addParticleTraining(const Particle &p, int classIdx)
{
    if (__training_particles[classIdx].size() < __maxTrainingVectors)
    {
        __training_particles[classIdx].push_back(p);
        return true;
    }
    return false;
}

// ==========================================================================
// Section: I/O =============================================================
// ==========================================================================
/* Show -------------------------------------------------------------------- */
std::ostream & operator <<(std::ostream &_out, const Particle &_p)
{
    _out << _p.micrograph << " " << _p.x << " " << _p.y << " " << _p.cost << " ";
    FOR_ALL_ELEMENTS_IN_MATRIX1D(_p.vec)
    _out << VEC_ELEM(_p.vec,i) << " ";
    _out << std::endl;
    return _out;
}

/* Read--------------------------------------------------------------------- */
void Particle::read(std::istream &_in, int _vec_size)
{
    _in >> micrograph >> x >> y >> cost;
    if (cost>=0)
    	status=1;
    else
    	status=0;
    vec.resize(_vec_size);
    _in >> vec;
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator <<(std::ostream &_out, const Classification_model &_m)
{
    _out << "#Already_processed_parameters...\n";
    _out << "#Micrographs_processed= " << _m.__micrographs_number << std::endl;
    _out << "#Micrographs_processed_points= ";

    for (int i = 0; i < _m.__micrographs_number; i++)
        _out << _m.__micrographs_scanned[i] << " ";
    _out << "\n";

    _out << "#Particles_picked_per_micrograph= ";
    for (int i = 0; i < _m.__micrographs_number; i++)
        _out << _m.__particles_picked[i] << " ";
    _out << "\n";

    _out << "#FalsePos_Picked= ";
    for (int i = 0; i < _m.__micrographs_number; i++)
        _out << _m.__falsePositives[i] << " ";
    _out << "\n";

    _out << "#Model_parameters..." << std::endl;
    _out << "#Vector_size= " << (_m.__training_particles[0][0].vec).size()
    << std::endl;
    _out << "#Class_Number= " << _m.__classNo << std::endl;

    for (int i = 0; i < _m.__classNo; i++)
    {
        int particlesNo = _m.__training_particles[i].size();
        _out << "#Particles_No= " << particlesNo << std::endl;
        for (int j = 0; j < particlesNo; j++)
            _out << _m.__training_particles[i][j];
    }

    return _out;
}

/* Read -------------------------------------------------------------------- */
std::istream & operator >>(std::istream &_in, Classification_model &_m)
{
    String dummy;
    int classNo, vec_size;
    _in >> dummy;
    _in >> dummy >> _m.__micrographs_number;
    _in >> dummy;

    _m.__micrographs_scanned.resize(_m.__micrographs_number);
    for (int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__micrographs_scanned[i];
    _in >> dummy;

    _m.__particles_picked.resize(_m.__micrographs_number);
    for (int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__particles_picked[i];
    _in >> dummy;

    _m.__falsePositives.resize(_m.__micrographs_number);
    for (int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__falsePositives[i];
    _in >> dummy;

    _in >> dummy >> vec_size;
    _in >> dummy >> _m.__classNo;
    _m.__training_particles.resize(_m.__classNo);
    Particle p;
    for (int i = 0; i < _m.__classNo; i++)
    {
        int particlesNo;
        _in >> dummy >> particlesNo;
        for (int j = 0; j < particlesNo; j++)
        {
            p.read(_in, vec_size);
            _m.addParticleTraining(p, i);
        }
    }
    return _in;
}

void AutoParticlePicking::loadModels(const FileName &fn_root)
{
    FileName fn_training = fn_root + "_training.txt";
    if (!fn_training.exists())
        return;

    // Load training vectors
    String dummy;
    std::ifstream fh_training;
    fh_training.open(fn_training.c_str());
    if (!fh_training)
        REPORT_ERROR(ERR_IO_NOTOPEN, fn_training);
    fh_training
    >> dummy >> __piece_xsize
    >> dummy >> __particle_radius;
    fh_training >> __selection_model;
    fh_training.close();
    __piece_overlap=2*__particle_radius;
    __scan_overlap = round(1.8*__particle_radius);

    // Load the mask
    __mask.type = READ_MASK;
    __mask.fn_mask = fn_root + "_mask.xmp";
    __mask.generate_mask();
    __mask.get_binary_mask().setXmippOrigin();
    classifyMask();
}

/* Save models ------------------------------------------------------------- */
void AutoParticlePicking::saveModels(const FileName &fn_root) const
{
    // Save the mask
    Image<double> save;
    typeCast(__mask.get_binary_mask(), save());
    save.write(fn_root + "_mask.xmp");

    // Save training vectors
    std::ofstream fh_training;
    fh_training.open((fn_root + "_training.txt").c_str());
    if (!fh_training)
        REPORT_ERROR(ERR_IO_NOWRITE, fn_root + ".training.txt");
    fh_training << "#piece_xsize=      " << __piece_xsize << std::endl
    << "#particle_radius=  " << __particle_radius << std::endl;
    fh_training << __selection_model << std::endl;
    fh_training.close();
}

/* Save particles ---------------------------------------------------------- */
int AutoParticlePicking::saveAutoParticles(const FileName &fn) const
{
    MetaData MD;
    size_t nmax = __auto_candidates.size();
    for (size_t n = 0; n < nmax; ++n)
    {
        const Particle &p = __auto_candidates[n];
        if (p.cost>0 && p.status==1)
        {
            size_t id = MD.addObject();
			MD.setValue(MDL_XCOOR, p.x, id);
			MD.setValue(MDL_YCOOR, p.y, id);
			MD.setValue(MDL_COST, p.cost, id);
			MD.setValue(MDL_ENABLED,1,id);
        }
    }
    MD.write(fn,MD_APPEND);
    return MD.size();
}

void AutoParticlePicking::saveAutoFeatureVectors(const FileName &fn, int Nvectors) const
{
    std::ofstream fh_auto;
    fh_auto.open(fn.c_str());
    if (!fh_auto)
        REPORT_ERROR(ERR_IO_NOWRITE,fn);
    fh_auto << Nvectors << " ";
    if (Nvectors == 0)
        fh_auto << "0\n";
    else
        fh_auto << VEC_XSIZE(__auto_candidates[0].vec) << std::endl;
    size_t Nauto = __auto_candidates.size();
    for (size_t i = 0; i < Nauto; i++)
    {
        const Particle &p = __auto_candidates[i];
        if (p.cost>0)
        	fh_auto << p;
    }
    fh_auto.close();
}

void AutoParticlePicking::loadAutoFeatureVectors(const FileName &fn)
{
    std::ifstream fh_auto;
    fh_auto.open(fn.c_str());
    if (!fh_auto)
        REPORT_ERROR(ERR_IO_NOTOPEN,fn);
    size_t Nauto;
    int vec_size;
    fh_auto >> Nauto >> vec_size;
    __auto_candidates.reserve(Nauto);
    Particle p;
    for (size_t i = 0; i < Nauto; i++)
    {
        p.read(fh_auto, vec_size);
        __auto_candidates.push_back(p);
    }
    fh_auto.close();
}

// ==========================================================================
// Section: Program interface ===============================================
// ==========================================================================
void ProgMicrographAutomaticPicking::readParams()
{
    fn_micrograph = getParam("-i");
    fn_model = getParam("--model");
    size = getIntParam("--particleSize");
    mode = getParam("--mode");
    if (mode == "train")
        fn_train = getParam("--mode", 1);
    Nthreads = getIntParam("--thr");
    fn_root = getParam("--outputRoot");
    fast = checkParam("--fast");
    incore = checkParam("--in_core");
}

void ProgMicrographAutomaticPicking::show()
{
    if (!verbose)
        return;
    std::cout << "Micrograph   : " << fn_micrograph << std::endl
    << "Model        : " << fn_model << std::endl << "Mode         : "
    << mode << std::endl << "Training file: " << fn_train << std::endl
    << "Threads      : " << Nthreads << std::endl << "Output       : "
    << fn_root << std::endl
    << "Fast         : " << fast << std::endl
    << "In core      : " << incore << std::endl
    ;
}

void ProgMicrographAutomaticPicking::defineParams()
{
    addUsageLine("Automatic particle picking for micrographs");
    addUsageLine("+The algorithm is designed to learn the particles from the user, as well as from its own errors.");
    addUsageLine("+The algorithm is fully described in [[http://www.ncbi.nlm.nih.gov/pubmed/19555764][this paper]].");
    addParamsLine("  -i <micrograph>               : Micrograph image");
    addParamsLine("  --outputRoot <rootname>       : Output rootname");
    addParamsLine("  --mode <mode>                 : Operation mode");
    addParamsLine("         where <mode>");
    addParamsLine("                    try              : Try to autoselect within the training phase.");
    addParamsLine("                    train <posfile=\"\">  : posfile contains the coordinates of manually picked particles");
    addParamsLine("                                     : <rootname>_auto_feature_vectors.txt contains the particle structure created by this program when used in automatic selection mode");
    addParamsLine("                                     : <rootname>_false_positives.xmd contains the list of false positives among the automatically picked particles");
    addParamsLine("                    autoselect  : Autoselect");
    addParamsLine("  --model <model_rootname>      : Bayesian model of the particles to pick");
    addParamsLine("  --particleSize <size>         : Particle size in pixels");
    addParamsLine("  [--thr <p=1>]                 : Number of threads for automatic picking");
    addParamsLine("  [--fast]                      : Perform a fast preprocessing of the micrograph (Fourier filter instead of Wavelet filter)");
    addParamsLine("  [--in_core]                   : Read the micrograph in memory");
    addExampleLine("Automatically select particles during training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode try ");
    addExampleLine("Training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode train manual.pos");
    addExampleLine("Automatically select particles after training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode autoselect");
}

void ProgMicrographAutomaticPicking::run()
{
    Micrograph m;
    m.open_micrograph(fn_micrograph);

    AutoParticlePicking *autoPicking = new AutoParticlePicking(fn_micrograph,&m,fast);
    autoPicking->setNumThreads(Nthreads);
    autoPicking->loadModels(fn_model);
    if (incore)
        autoPicking->readMicrograph();
    FileName familyName=fn_model.removeDirectories();
    FileName fnAutoParticles=familyName+"@"+fn_root+"_auto.pos";
    FileName fnVectors=fn_root + "_auto_feature_vectors_"+familyName+".txt";
    if (mode == "autoselect" || mode=="try")
    {
        autoPicking->automaticallySelectParticles();
        int Nparticles=autoPicking->saveAutoParticles(fnAutoParticles);
        if (mode=="try")
            autoPicking->saveAutoFeatureVectors(fnVectors,Nparticles);
    }
    else
    {
        MetaData MD;
        // Insert all true positives
        if (fn_train!="")
        {
            MD.read(fn_train);
            int x, y;
            FOR_ALL_OBJECTS_IN_METADATA(MD)
            {
                MD.getValue(MDL_XCOOR, x, __iter.objId);
                MD.getValue(MDL_YCOOR, y, __iter.objId);
                m.add_coord(x, y, 0, 1);
            }
        }

        // Insert all false positives
        if (fnAutoParticles.existsTrim())
        {
            MD.read(fnAutoParticles);
            if (MD.size() > 0)
            {
                autoPicking->loadAutoFeatureVectors(fnVectors);
                int idx=0;
                FOR_ALL_OBJECTS_IN_METADATA(MD)
                {
                    int enabled;
                    MD.getValue(MDL_ENABLED,enabled,__iter.objId);
                    if (enabled==-1)
                    {
                        autoPicking->__auto_candidates[idx].status = 0;
                        autoPicking->__auto_candidates[idx].cost = -1;
                        autoPicking->__rejected_particles.push_back(autoPicking->__auto_candidates[idx]);
                    }
                    else
                    {
                    	double cost;
                        MD.getValue(MDL_COST,cost,__iter.objId);
                        if (cost>0)
                        {
                        	autoPicking->__auto_candidates[idx].status = 1;
                        }
                        else
                            autoPicking->__auto_candidates[idx].status = 0;
                    }
                    ++idx;
                }
            }
        }
        autoPicking->learnParticles(size / 2);
        autoPicking->saveModels(fn_model);
        if (fileExists(fnVectors))
        	unlink(fnVectors.c_str());
        if (fnAutoParticles.existsTrim())
        	autoPicking->saveAutoParticles(fnAutoParticles);
    }
}
