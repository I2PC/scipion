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

#include "micrograph_automatic_picking_for_qt.h"
#include <data/filters.h>
#include <data/rotational_spectrum.h>
#include <reconstruction/denoise.h>
#include <data/xmipp_fft.h>
#include <reconstruction/fourier_filter.h>
#include <algorithm>

//#define DEBUG_AUTO
//#define DEBUG_MORE_AUTO
//#define DEBUG_PREPARE
//#define DEBUG_CLASSIFY
//#define DEBUG_BUILDVECTOR
//#define DEBUG_IMG_BUILDVECTOR

// ==========================================================================
// Section: Piece Preprocessing =============================================
// ==========================================================================
/* Prepare piece ----------------------------------------------------------- */
static pthread_mutex_t preparePieceMutex = PTHREAD_MUTEX_INITIALIZER;
bool AutoParticlePickingQt::prepare_piece(MultidimArray<double> &piece,
                                        MultidimArray<int> &ipiece,
                                        MultidimArray<double> &original_piece)
{
    original_piece=piece;

#ifdef DEBUG_PREPARE

    Image<double> save;
    save() = piece;
    save.write("PPPpiece0.xmp");
#endif

    // Denoise the piece
    WaveletFilter denoiser;
    denoiser.denoising_type = WaveletFilter::BAYESIAN;
    denoiser.scale = 3;
    denoiser.output_scale = 1;
    denoiser.produceSideInfo();
    denoiser.apply(piece);
    if (!(piece(0, 0) == piece(0, 0)))
        return false;
#ifdef DEBUG_PREPARE

    save() = piece;
    save.write("PPPpiece1.xmp");
#endif

    if (__output_scale==0)
    {
        MultidimArray<double> auxPiece=piece;
        pyramidExpand(1, piece, auxPiece);
#ifdef DEBUG_PREPARE

        save() = piece;
        save.write("PPPpiece1.5.xmp");
#endif

    }

    // Band pass filter
    pthread_mutex_lock( &preparePieceMutex );
    FourierFilter Filter;
    Filter.FilterShape = RAISED_COSINE;
    Filter.FilterBand = BANDPASS;
    Filter.w1 = __highpass_cutoff;
    Filter.w2 = 1.0/(__particle_radius/(__reduction*5.0));
    Filter.raised_w = XMIPP_MIN(0.02, __highpass_cutoff);
    Filter.generateMask(piece);
    Filter.applyMaskSpace(piece);
    STARTINGX(piece) = STARTINGY(piece) = 0;
#ifdef DEBUG_PREPARE

    save() = piece;
    save.write("PPPpiece2.xmp");
#endif

    pthread_mutex_unlock( &preparePieceMutex );

    // Reject 5% of the outliers
    reject_outliers(piece, 5.0);

#ifdef DEBUG_PREPARE

    save() = piece;
    save.write("PPPpiece2_5.xmp");
#endif

    // Equalize histogram
    histogram_equalization(piece, __gray_bins);

    typeCast(piece,ipiece);

#ifdef DEBUG_PREPARE

    save() = piece;
    save.write("PPPpiece3.xmp");
#endif

    if (!original_piece.sameShape(piece))
        selfScaleToSize(LINEAR, original_piece, YSIZE(piece), XSIZE(piece));

    // Reject 5% of the outliers
    reject_outliers(original_piece, 5.0);

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
bool AutoParticlePickingQt::build_vector(const MultidimArray<int> &piece,
                                       const MultidimArray<double> &original_piece,
                                       int _x, int _y,
                                       Matrix1D<double> &_result)
{
#ifdef DEBUG_BUILDVECTOR
    std::cout << "build_vector(" << _x << "," << _y << "," << "_result)" << std::endl;
#endif

    // Resize the output and make same aliases
    int angleBins=__sector.size();
    _result.initZeros(32 // Histogram of the original piece
                      +__radial_bins*(__gray_bins-1) // radial histograms
                      +(angleBins-1)+(angleBins-1)*angleBins // sector correlations
                      +(2*__radial_bins-19)*angleBins // ring correlations
                     );
    const MultidimArray<int> &mask =         __mask.get_binary_mask();
    const MultidimArray<int> &classif1 =  (*(__mask_classification[0]));
    const MultidimArray<int> &classif2 =  (*(__mask_classification[1]));

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
    std::vector< MultidimArray<int> > radial_val;
    for (int i=0; i<__radial_val.size(); i++)
    {
        MultidimArray<int> dummyi=*(__radial_val[i]);
        radial_val.push_back(dummyi);
    }

    std::vector< MultidimArray<double> > sector;
    MultidimArray<double> dummyd=*(__sector[0]);
    for (int i=0; i<__sector.size(); i++)
        sector.push_back(dummyd);

    std::vector< MultidimArray<double> > ring;
    dummyd=*(__ring[0]);
    for (int i=0; i<__ring.size(); i++)
        ring.push_back(dummyd);

    Histogram1D histogramOriginal;
    histogramOriginal.init(0,255,32);

    // Put the image values into the corresponding radial bins
    FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
    {
        int idx1 = A2D_ELEM(classif1, i, j);
        if (idx1 != -1)
        {
            int val = A2D_ELEM(piece, _y + i, _x + j);
            double orig_val = A2D_ELEM(original_piece, _y + i, _x + j);
            INSERT_VALUE(histogramOriginal,orig_val);

            MultidimArray<int> &auxi=radial_val[idx1];
            int aux=A1D_ELEM(radial_idx,idx1);
            A1D_ELEM(auxi,aux) = val;
            A1D_ELEM(radial_idx,idx1)++;
            int idx2 = A2D_ELEM(classif2, i, j);
            if (idx2 != -1)
            {
                MultidimArray<double> &auxd=sector[idx2];
                A1D_ELEM(auxd,idx1) += val;
            }
        }

        // Get particle
#ifdef DEBUG_IMG_BUILDVECTOR
        if (debug_go)
        {
            save(i, j) = val;
            if (A2D_ELEM(mask, i, j))
            {
                savefg(i, j) = val;
                saveOrig(i, j) =  original_piece(_y+i, _x+j);
            }
        }
#endif

    }

    // Compute the sector averages and reorganize the data in rings
    for (int j = 0; j < angleBins; j++)
    {
        MultidimArray<int> &Nsector_j=*(__Nsector[j]);
        MultidimArray<double> &sector_j=sector[j];
        FOR_ALL_ELEMENTS_IN_ARRAY1D(sector_j)
        {
            if (DIRECT_A1D_ELEM(Nsector_j,i)>0)
            {
                DIRECT_A1D_ELEM(sector_j,i)/=DIRECT_A1D_ELEM(Nsector_j,i);
                MultidimArray<double> &ring_i=ring[i];
                DIRECT_A1D_ELEM(ring_i,j)=DIRECT_A1D_ELEM(sector_j,i);
            }
        }
    }

    // Store the histogram of
    int idx_result=0;
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
    for (int step = 1; step<angleBins; step++)
    {
        sectorCorr.initZeros(angleBins);
        for (int i = 0; i<angleBins; i++)
        {
            DIRECT_A1D_ELEM(sectorCorr,i)=correlationIndex(
                                              sector[i],
                                              sector[intWRAP(i+step,0,angleBins-1)]);
        }
        VEC_ELEM(_result,idx_result++) = sectorCorr.computeAvg();
        correlation_vector_no_Fourier(sectorCorr,sectorCorr,sectorAutocorr);
        for (int j = 0; j < XSIZE(sectorAutocorr); j++)
            VEC_ELEM(_result,idx_result++) = DIRECT_A1D_ELEM(sectorAutocorr,j);
    }

    // Compute the correlation in rings
    MultidimArray<double> ringCorr;
    for (int step = 1; step<=2; step++)
        for (int i = 8; i<__radial_bins-step; i++)
        {
            correlation_vector_no_Fourier(ring[i],ring[i+step],ringCorr);
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
void AutoParticlePickingQt::learnParticles(int _ellipse_radius)
{
    if (__particle_radius==0)
    {
        __particle_radius = _ellipse_radius;
        __penalization = 10;
        __gray_bins = 8;
        __radial_bins = 16;
        if (4*__particle_radius<512)
            __piece_xsize = 512;
        else
            __piece_xsize = NEXT_POWER_OF_2(6*__particle_radius);
        __highpass_cutoff = 0.02;
        __mask_size = 2*__particle_radius;
        __min_distance_between_particles = __particle_radius/2;
        __scan_overlap = ROUND(__mask_size*0.9);
        __reduction=(int)std::pow(2.0, __output_scale);
    }
    __learn_overlap=0;

    std::cerr << "\n------------------Learning Phase-----------------------\n";
    createMask();

    // Find in the particle list the images to learn
    std::vector<int> all_idx;
    int num_part = __m->ParticleNo();
    for (int i = 0; i < num_part; i++)
        if (__m->coord(i).valid && __m->coord(i).label != __auto_label)
            all_idx.push_back(i);

    // If there is nothing to learn, return
    if (all_idx.size() == 0 && !__is_model_loaded)
    {
        std::cerr << "No valid particles marked." << std::endl;
        return;
    }

    // If we have already learned or autoselected, delete the training vectors
    if (__learn_particles_done || __autoselection_done)
    {
        __training_model.clear();
        __training_model=__training_loaded_model;
    }

    // Actually learn
    buildPositiveVectors(all_idx, __training_model);
    getAutoTruePositives(__training_model);
    buildNegativeVectors(__training_model,false);
    getAutoFalsePositives(__training_model);
    buildNegativeVectors(__training_model,true);

    __training_model.addMicrographScanned(count_scanning_pos());
    __selection_model = __training_model;

    __learn_particles_done = true;
    std::cerr << "Learning process finished..." << std::endl;
}
/* Get the false-positive particles----------------------------------------- */
void AutoParticlePickingQt::getAutoFalsePositives(Classification_modelQt &_training_model)
{
    // Add the false positives given by the threshold
    int Nrejected=0;
    int imax = __auto_candidates.size();
    for (int i = 0; i < imax; i++)
        if (__auto_candidates[i].status == 1 &&
            __auto_candidates[i].cost<__minCost)
        {
            if (_training_model.addParticleTraining(__auto_candidates[i], 2))
                Nrejected++;
        }

    // Add the manually selected false positives
    imax = __rejected_particles.size();
    _training_model.addFalsePositives(imax);
    for (int i = 0; i < imax; i++)
        if (_training_model.addParticleTraining(__rejected_particles[i], 2))
            ++Nrejected;
    std::cout << Nrejected << " false positives are considered\n";
}

/* Get the false-positive particles----------------------------------------- */
void AutoParticlePickingQt::getAutoTruePositives(
    Classification_modelQt &_training_model)
{
    // Add the true positives given by the threshold
    int Naccepted=0;
    int imax = __auto_candidates.size();
    for (int i = 0; i < imax; i++)
        if (__auto_candidates[i].status == 1 &&
            __auto_candidates[i].cost>=__minCost)
        {
            if (_training_model.addParticleTraining(__auto_candidates[i], 0))
                Naccepted++;
        }
    std::cout << Naccepted << " true positives are considered\n";
}

/* Build training vectors ---------------------------------------------------*/
void AutoParticlePickingQt::buildPositiveVectors(std::vector<int> &_idx,
        Classification_modelQt &_model)
{
    std::cerr << "Building " << _idx.size()
    << " particle vectors for this image. Please wait..." << std::endl;
    _model.addMicrographItem();
    int num_part = _idx.size();
    int width, height;
    __m->size(width, height);
    Matrix1D<double> v;
    int numParticles = 0;

    MultidimArray<char> visited(num_part);
    MultidimArray<double> piece, original_piece;
    MultidimArray<int> ipiece;
    while (visited.sum() < num_part)
    {
        int part_i = 0;

        // get the first un-visited particle in the array
        while (visited(part_i) == 1)
            part_i++;
        if (part_i >= num_part)
            break;
        visited(part_i) = 1;

        // Get the piece containing that particle
        int part_idx = _idx.at(part_i);
        int x = __m->coord(part_idx).X;
        int y = __m->coord(part_idx).Y;
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
            ParticleQt p;
            p.x = x;
            p.y = y;
            p.idx = part_idx;
            p.vec = v;
            p.status = 1;
            p.cost = 1.0;
            if (_model.addParticleTraining(p, 0))
                numParticles++;
        }

        // make vector from the neighbours
        std::vector< Matrix1D<int> > neighbour;
        neighbour.reserve(num_part);

        find_neighbour(piece, _idx, part_i, x, y, posx, posy, visited, neighbour);

        for (int i = 0; i < neighbour.size(); i++)
        {
            part_i = neighbour.at(i)(0);
            part_idx = _idx.at(part_i);
            posx = neighbour.at(i)(1);
            posy = neighbour.at(i)(2);
            success = build_vector(ipiece, original_piece, posx, posy, v);
            visited(part_i) = 1;
            if (success)
            {
                ParticleQt p;
                p.x = __m->coord(part_idx).X;
                p.y = __m->coord(part_idx).Y;
                p.idx = part_idx;
                p.vec = v;
                p.status = 1;
                p.cost = 1.0;
                if (_model.addParticleTraining(p, 0))
                    numParticles++;
            }
        }
    }
    _model.addParticlePicked(numParticles);
}

/* Build vector from non particles------------------------------------------ */
void AutoParticlePickingQt::buildNegativeVectors(Classification_modelQt &_model,
        bool checkForPalsePostives)
{
    if (checkForPalsePostives)
        std::cerr << "Building automatic false positives ..." << std::endl;
    else
        std::cerr << "Building non particles ..." << std::endl;
    const MultidimArray<int> &mask = __mask.get_binary_mask();

    // Setup a classification model with the already known data
    std::vector < MultidimArray<double> > features;
    Matrix1D<double> probs;
    if (checkForPalsePostives)
    {
        // Gather all information for classification
        __selection_model = __training_model;

        // Prepare data to be classified
        getFeatures(__selection_model,features);
        getClassesProbabilities(__selection_model,probs);

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
    Matrix1D<double> v, aux1, aux2;

    int N = 1, Nnonparticles=0, Nfalsepositives=0;

    // We do not want any overlap for this process,since it is only for
    // counting the non particles and calculating their features. For
    // the process of automatic selecting we will want an overlap so we
    // do not miss any particle.
    MultidimArray<double> piece, original_piece;
    MultidimArray<int> ipiece;
    while (get_corner_piece(piece, top, left, skip_y,
                            next_skip_x, next_skip_y, next_top, next_left, 0,
                            true))
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
        while (get_next_scanning_pos(piece, next_posx, next_posy, skip_x,
                                     skip_y, __learn_overlap))
        {
            // Check if there is any particle around
            if (!anyParticle(left + posx  * __reduction,
                             top + posy  * __reduction,
                             XSIZE(mask) * __reduction))
            {
                if (build_vector(ipiece, original_piece, posx, posy, v))
                {
                    // Build the Particle structure
                    ParticleQt P;
                    P.x = left + posx * __reduction;
                    P.y = top + posy * __reduction;
                    P.idx = -1;
                    P.status = 1;
                    P.vec = v;
                    P.cost = -1;
                    double cost;
                    if (!checkForPalsePostives)
                    {
                        if (_model.addParticleTraining(P, 1))
                            Nnonparticles++;
                    }
                    else
                    {
                        int votes=__selection_model.isParticle(v,cost,aux1,aux2);
                        if (votes>5)
                        {
                            if (_model.addParticleTraining(P, 2))
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
    bool operator()(const ParticleQt& rpStart, const ParticleQt& rpEnd)
    {
        return rpStart.cost < rpEnd.cost;
    }
};

/* Automatic phase ----------------------------------------------------------*/
int AutoParticlePickingQt::automaticallySelectParticles()
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
    std::vector < MultidimArray<double> > features;
    Matrix1D<double> probs;
    getFeatures(__selection_model,features);
    getClassesProbabilities(__selection_model,probs);

    // Initialize classifier
    __selection_model.initNaiveBayesEnsemble(features, probs, 8,
            __penalization, 10, 1, 1, "mm");
#ifdef DEBUG_AUTO

    std::cout << "Probabilities of the classes:"
    << probs.transpose() << std::endl;
#endif

    // Automatically select particles with threads
    pthread_t * th_ids = new pthread_t[__numThreads];
    AutomaticallySelectThreadParamsQt * th_args=
        new AutomaticallySelectThreadParamsQt[__numThreads];
    for( int nt = 0 ; nt < __numThreads ; nt ++ )
    {
        th_args[nt].autoPicking = this;
        th_args[nt].idThread = nt;
        pthread_create( &th_ids[nt], NULL,
                        automaticallySelectParticlesThreadQt, &th_args[nt]);
    }

    // Waiting for threads to finish
    for( int nt=0; nt<__numThreads; nt++)
        pthread_join(th_ids[nt], NULL);

    // Thread structures are not needed any more
    delete [] th_ids;
    delete [] th_args;

    // Sort particles by cost
    std::sort(__auto_candidates.begin(),__auto_candidates.end(),
              SAscendingParticleSort());

#ifdef DEBUG_AUTO

    std::cerr << "Number of automatically selected particles = "
    << __auto_candidates.size() << std::endl;
#endif

    // Reject the candidates that are pointing to the same particle
    int Nalive = reject_within_distance(__auto_candidates, __particle_radius,
                                        false);

#ifdef DEBUG_AUTO

    std::cerr << "Number of automatically selected particles after distance rejection = "
    << __auto_candidates.size() << std::endl;
#endif

    // Apply a second classifier for classifying between particle
    // and false positive. For that, remove the middle class (background)
    if (features.size()==3)
    {
        int imax;
    	Matrix1D<double> aux1, aux2;
        if (Nalive > 0)
        {
            __selection_model2.clear();
            __selection_model2.init(2);
            std::vector < MultidimArray<double> >::iterator featuresIterator=
                features.begin();
            featuresIterator++;
            features.erase(featuresIterator);
            probs(1)=probs(2);
            probs.resize(2);
            probs/=probs.sum();
            __selection_model2.initNaiveBayesEnsemble(features, probs, 8,
                    __penalization,10,1,1,"mm");
#ifdef DEBUG_AUTO

            std::cout << "Second classification\n";
#endif

            imax = __auto_candidates.size();
            Nalive = 0;
            for (int i = 0; i < imax; i++)
                if (__auto_candidates[i].status == 1)
                {
                    double p;
                    int votes=__selection_model2.isParticle(__auto_candidates[i].vec,p,aux1,aux2);
                    if (votes<8)
                    {
                        __auto_candidates[i].status=0;
                        __auto_candidates[i].cost=-1;
#ifdef DEBUG_AUTO

                        std::cout << __auto_candidates[i].x << ", "
                        << __auto_candidates[i].y
                        << " is considered as a false positive\n";
#endif

                    }
                    else
                        Nalive++;
                }
        }

        // Apply a third classifier to distinguish between particles and
        // very tough particles
        if (Nalive > 0)
        {
            __selection_model3.clear();
            __selection_model3.init(2);
            // Remove from the error class, all those errors that
            // the previous classifier was able to classify correctly
            std::vector<int> toKeep;
            const MultidimArray<double> &features_1=features[1];
            Matrix1D<double> trialFeatures(XSIZE(features_1));
            size_t rowLength=XSIZE(features_1)*sizeof(double);
            double cost;
            for (int i=0; i<YSIZE(features_1); i++)
            {
                memcpy(&VEC_ELEM(trialFeatures,0),&A2D_ELEM(features_1,i,0),rowLength);
                int votes=__selection_model2.isParticle(trialFeatures,cost, aux1, aux2);
                if (votes<8)
                    toKeep.push_back(i);
            }

            if (toKeep.size()>0)
            {
                MultidimArray<double> difficultParticles;
                difficultParticles.initZeros(toKeep.size(),XSIZE(features[1]));
                FOR_ALL_ELEMENTS_IN_ARRAY2D(difficultParticles)
                difficultParticles(i,j)=features[1](toKeep[i],j);
                int NErrors=YSIZE(features[1]);
                features.pop_back();
                features.push_back(difficultParticles);

                __selection_model3.initNaiveBayesEnsemble(features, probs, 8,
                        __penalization, 10, 1, 1, "mm");
#ifdef DEBUG_AUTO

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
                        int votes=__selection_model3.isParticle(__auto_candidates[i].vec,p, aux1, aux2);
                        if (votes<8)
                        {
                            __auto_candidates[i].status=0;
                            __auto_candidates[i].cost=-1;
#ifdef DEBUG_AUTO

                            std::cout << __auto_candidates[i].x << ", "
                            << __auto_candidates[i].y
                            << " is considered as a false positive 2\n";
#endif

                        }
                        else
                            Nalive++;
                    }
            }
        }
    }

    if (Nalive>0)
    {
        __auto_label=-1;
        for (int i = 0; i < __m->LabelNo(); i++)
            if (__m->get_label(i)=="auto")
            {
                __auto_label=i;
                break;
            }
        if (__auto_label==-1)
        {
            __m->add_label("auto");
            for (int i = 0; i < __m->LabelNo(); i++)
                if (__m->get_label(i)=="auto")
                {
                    __auto_label=i;
                    break;
                }
        }

        int imax = __auto_candidates.size();
        // Get the maximum and minimum cost
        bool first=true;
        double minCost, maxCost;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].status == 1)
            {
                if (!std::isinf(__auto_candidates[i].cost))
                {
                    if (first || __auto_candidates[i].cost<minCost)
                    {
                        minCost=__auto_candidates[i].cost;
                        first=false;
                    }
                    if (first || __auto_candidates[i].cost>maxCost)
                    {
                        maxCost=__auto_candidates[i].cost;
                        first=false;
                    }
                }
            }

        // Insert selected particles in the result
        int idxMicrograph=0;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].status == 1)
            {
                __auto_candidates[i].idx = idxMicrograph;
#ifdef DEBUG_AUTO

                std::cout << "Particle coords " << __auto_candidates[i].x << ", "
                << __auto_candidates[i].y << std::endl;
#endif

                if (std::isinf(__auto_candidates[i].cost))
                    __auto_candidates[i].cost=1;
                else
                    __auto_candidates[i].cost=
                        (__auto_candidates[i].cost-maxCost)/
                        (minCost-maxCost);
                __m->add_coord(__auto_candidates[i].x,
                               __auto_candidates[i].y,
                               __auto_label,
                               __auto_candidates[i].cost);
                idxMicrograph++;
            }
            else
                __auto_candidates[i].status = 0;
    }

    __autoselection_done = true;
    std::cout << "\nAutomatic process finished. Number of particles found: "
    << Nalive << std::endl;
    return Nalive;
}

static pthread_mutex_t particleAdditionMutex = PTHREAD_MUTEX_INITIALIZER;

void * automaticallySelectParticlesThreadQt(void * args)
{
    // Pick the input parameters
    AutomaticallySelectThreadParamsQt * prm=
        (AutomaticallySelectThreadParamsQt *)args;
    AutoParticlePickingQt *autoPicking=prm->autoPicking;
    int idThread=prm->idThread;

    int thVotes=1;
    if (autoPicking->__selection_model.__training_particles.size()==2 ||
        (autoPicking->__selection_model.__training_particles.size()==3 &&
         autoPicking->__selection_model.__training_particles[2].size()==0))
        thVotes=8;

    //top,left corner of the piece
    int top = 0, left = 0, next_top = 0, next_left = 0;
    // If the piece available is small then include the scanned part
    // because we need bigger image for denoising but for scanning
    // particles we skip the already scanned part
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    Matrix1D<double> v,aux1,aux2;
    int N = 1, particle_idx = 0, Nscanned=0;
    MultidimArray<double> piece, original_piece;
    MultidimArray<int> ipiece;
    const MultidimArray<int> &mask = autoPicking->__mask.get_binary_mask();
    std::vector< ParticleQt > threadCandidates;
    do
    {
        bool isMine=(N%autoPicking->__numThreads==idThread);
        bool pieceOK=autoPicking->get_corner_piece(piece, top, left,
                     skip_y, next_skip_x, next_skip_y, next_top, next_left,
                     autoPicking->__piece_overlap,isMine);
        if (!pieceOK)
            break;
        if (isMine)
        {
            if (idThread==0)
                std::cerr << ".";
#ifdef DEBUG_MORE_AUTO

            std::cerr << "thread " << idThread
            << " processing piece " << N << "...\n";
            std::cerr << "    (top,left)=" << top << "," << left
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

            std::cerr << "Skip(y,x)=" << skip_y << "," << skip_x << std::endl;
#endif

            // Scan this piece
            int posx = 0, next_posx = 0, posy = 0, next_posy = 0;
            next_posx = posx = skip_x + XSIZE(mask) / 2;
            next_posy = posy = skip_y + YSIZE(mask) / 2;

            while (autoPicking->get_next_scanning_pos(
                       piece, next_posx, next_posy, skip_x, skip_y,
                       autoPicking->__scan_overlap))
            {
                // COSS: Uncomment the next sentence for fast debugging
                // if (rnd_unif(0,1)<0.98) continue;

#ifdef DEBUG_MORE_AUTO
                std::cerr << "Pos(y,x)=" << posy << "," << posx
                << " Micro(y,x)=" << posy*autoPicking->__reduction + top
                << "," << posx*autoPicking->__reduction + left
                << " Next pos(y,x)=" << next_posy << "," << next_posx
                << std::endl;
#endif

                if (autoPicking->build_vector(ipiece, original_piece, posx, posy, v))
                {
                    double cost;
                    int votes=autoPicking->__selection_model.isParticle(v,cost,aux1,aux2);
                    if (votes>thVotes)
                    {
#ifdef DEBUG_MORE_AUTO
                        std::cout << "Particle Found: "
                        << left + posx * autoPicking->__reduction << ","
                        << top  + posy *autoPicking->__reduction
                        << " votes=" << votes
                        << std::endl;
                        std::cout << "Press any key to continue...\n";
                        char c;
                        std::cin >> c;
#endif

                        // Build the Particle structure
                        ParticleQt P;
                        P.x = left + posx * autoPicking->__reduction;
                        P.y = top + posy * autoPicking->__reduction;
                        P.idx = -1;
                        P.status = 1;
                        P.vec = v;
                        P.cost = cost;
                        threadCandidates.push_back(P);
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

    int imax=threadCandidates.size();
    pthread_mutex_lock( &particleAdditionMutex );
    for (int i=0; i<imax; i++)
        autoPicking->__auto_candidates.push_back(threadCandidates[i]);
    pthread_mutex_unlock( &particleAdditionMutex );
}

/* Filter particles --------------------------------------------------------*/
//To calculate the euclidean distance between to points
double euclidean_distance(const ParticleQt &p1, const ParticleQt &p2)
{
    double dx=(p1.x - p2.x);
    double dy=(p1.y - p2.y);
    return sqrt(dx*dx + dy*dy);
}

int AutoParticlePickingQt::reject_within_distance(
    std::vector<ParticleQt> &_Input, double _min_dist,
    bool _reject_both)
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
void AutoParticlePickingQt::createMask()
{
    if (XSIZE(__mask.get_binary_mask()) != 0)
        return;

    int xsize = __mask_size / __reduction;
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
void AutoParticlePickingQt::classifyMask()
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

    MultidimArray<int> *classif1 = new MultidimArray<int>;
    classif1->resize(mask);
    classif1->initConstant(-1);

    MultidimArray<int> *classif2 = new MultidimArray<int>;
    classif2->resize(mask);
    classif2->initConstant(-1);
    const double deltaAng=(PI/8.0);

    MultidimArray<int> Nrad(__radial_bins);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
    {
        double radius = sqrt((double)(i * i + j * j));
        double angle = atan2((double)i, (double)j);
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
        MultidimArray<int> *aux1 = new MultidimArray<int>;
        aux1->initZeros(Nrad(i));
        __radial_val.push_back(aux1);
    }

    // Create the holders for the radius values in classif1 and classif2
    int angleBins=classif2->computeMax()+1;
    for (int i = 0; i < angleBins; i++)
    {
        MultidimArray<double> *aux2 = new MultidimArray<double>;
        aux2->initZeros(__radial_bins);
        __sector.push_back(aux2);

        MultidimArray<int> *iaux2 = new MultidimArray<int>;
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
        MultidimArray<double> *aux3 = new MultidimArray<double>;
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
int AutoParticlePickingQt::count_scanning_pos()
{
	__piece_overlap=2*__particle_radius;
	const MultidimArray<int> &mask = __mask.get_binary_mask();
    int top = 0, left = 0, next_top = 0, next_left = 0;
    int skip_x = 0, skip_y = 0, next_skip_x = 0, next_skip_y = 0;
    int Nscanned=0;
    MultidimArray<double> piece;
    while (get_corner_piece(piece, top, left, skip_y,
                            next_skip_x, next_skip_y, next_top,
                            next_left, __piece_overlap, true))
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
bool AutoParticlePickingQt::get_next_scanning_pos(
    const MultidimArray<double> &piece,
    int &_x, int &_y, int _skip_x, int _skip_y, int overlap) const
{
    const MultidimArray<int> &mask = __mask.get_binary_mask();

    if (_x + XSIZE(mask) / 2 > XSIZE(piece) ||
        _y + YSIZE(mask) / 2 > YSIZE(piece))
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
void AutoParticlePickingQt::find_neighbour(const MultidimArray<double> &piece,
        std::vector<int> &_idx, int _index,
        int _x, int _y, int _posx, int _posy, MultidimArray<char> &_visited,
        std::vector< Matrix1D<int> > &_nbr)
{
    int piece_xsize = XSIZE(piece);
    int piece_ysize = YSIZE(piece);
    const MultidimArray<int> &mask = __mask.get_binary_mask();

    // If all the particles are visited
    if (_visited.sum() == XSIZE(_visited))
        return;

    int current_part = _index + 1;
    int current_part_idx = _idx.at(current_part);
    _nbr.clear();
    _nbr.reserve(XSIZE(_visited));

    int top = CEIL((double)_y / __reduction) - _posy;
    int left = CEIL((double)_x / __reduction) - _posx;
    int bottom = top + FLOOR((double)piece_ysize / __reduction)-1;
    int right = left + FLOOR((double)piece_xsize / __reduction)-1;
    int xmask2 = CEIL(XSIZE(mask) / 2);
    int ymask2 = CEIL(YSIZE(mask) / 2);

    while (current_part < XSIZE(_visited))
    {
        //find the next unvisited particle
        while (_visited(current_part))
        {
            current_part++;
            if (current_part==XSIZE(_visited))
                break;
        }
        if (current_part==XSIZE(_visited))
            break;
        current_part_idx = _idx.at(current_part);

        //check if it is neighbour or not
        int nx = ROUND((double)__m->coord(current_part_idx).X / __reduction);
        int ny = ROUND((double)__m->coord(current_part_idx).Y / __reduction);
        if ((nx - xmask2 > left) &&
            (ny - ymask2 > top) &&
            (nx + xmask2 < right) &&
            (ny + ymask2 < bottom))
        {
            Matrix1D<int> current_nbr;
            current_nbr.initZeros(3);
            current_nbr(0) = current_part;
            current_nbr(1) = nx - left;
            current_nbr(2) = ny - top;
            _nbr.push_back(current_nbr);
        }
        current_part++;
    }
}

bool AutoParticlePickingQt::anyParticle(int posx, int posy, int rect_size)
{
    int num_part = __m->ParticleNo();
    std::vector<Particle_coords> selected_particles = __m->Particles();
    for (int i = 0; i < num_part; i++)
    {
        if (__m->coord(i).valid && __m->coord(i).label != __auto_label)
        {
            int _x = selected_particles[i].X;
            int _y = selected_particles[i].Y;

            if((_x > posx - rect_size) && (_x < posx + rect_size))
            {
                if((_y > posy - rect_size) && (_y < posy + rect_size))
                    return true;
            }
        }
    }
    return false;
}

/* Get piece --------------------------------------------------------------- */
// to get the piece containing (x,y) of size xsize,ysize
// return the position of x,y in the piece in posx,posy
void AutoParticlePickingQt::get_centered_piece(MultidimArray<double> &piece,
        int _x, int _y, int &_posx, int &_posy)
{
    piece.initZeros(__piece_xsize, __piece_xsize);
    int startx = _x - ROUND(__piece_xsize / 2);
    int endx   = _x + ROUND(__piece_xsize / 2);
    int starty = _y - ROUND(__piece_xsize / 2);
    int endy   = _y + ROUND(__piece_xsize / 2);
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
    const Micrograph &micrograph=*__m;
    for (int i = 0; i < __piece_xsize; i++)
        for (int j = 0; j < __piece_xsize; j++)
            DIRECT_A2D_ELEM(piece, i, j) = micrograph(starty + i, startx + j);
}

// Get a piece whose top-left corner is at the desired position (if possible)
bool AutoParticlePickingQt::get_corner_piece(
    MultidimArray<double> &piece,
    int _top, int _left, int _skip_y,
    int &_next_skip_x, int &_next_skip_y, int &_next_top, int &_next_left,
    int overlap, bool copyPiece) const
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
        const Micrograph& micrograph=*__m;
        for (int i = 0; i < __piece_xsize; i++)
            for (int j = 0; j < __piece_xsize; j++)
                DIRECT_A2D_ELEM(piece,i, j) = micrograph(_top + i, _left + j);
    }

    return true;
}

// ==========================================================================
// Section: Constructors, Initializers, Gets ================================
// ==========================================================================
/* Clear ------------------------------------------------------------------- */
void Classification_modelQt::clear()
{
    __particles_picked.clear();
    __micrographs_number = 0;
    __falsePositives.clear();

    int imax = __training_particles.size();
    for (int i = 0; i < imax; i++)
        __training_particles[i].clear();
}

/* Import model ------------------------------------------------------------ */
void Classification_modelQt::import_model(const Classification_modelQt &_model)
{
    // Import particles
    int jmax = _model.__classNo;
    for (int j = 0; j < jmax; j++)
    {
        int imax = _model.__training_particles[j].size();
        for (unsigned long i = 0; i < imax; i++)
            addParticleTraining(_model.__training_particles[j][i], j);
    }

    // Import parameters
    int imax = _model.__micrographs_number;
    for (int i = 0; i < imax; i++)
    {
        addParticlePicked(_model.__particles_picked[i]);
        addMicrographScanned(_model.__micrographs_scanned[i]);
        addFalsePositives(_model.__falsePositives[i]);
    }

    // Update the number of micrographs
    __micrographs_number += _model.__micrographs_number;
}

/* Initialize -------------------------------------------------------------- */
void Classification_modelQt::initNaiveBayesEnsemble(
    const std::vector < MultidimArray<double> > &features,
    const Matrix1D<double> &probs, int discreteLevels,
    double penalization, int numberOfClassifiers,
    double samplingFeatures, double samplingIndividuals,
    const std::string &newJudgeCombination)
{
    __bayesEnsembleNet=new EnsembleNaiveBayes(features, probs, discreteLevels,
                       numberOfClassifiers, samplingFeatures, samplingIndividuals,
                       newJudgeCombination);
    int K=features.size();
    Matrix2D<double> cost(K,K);
    cost.initConstant(1);
    for (int i=0; i<MAT_XSIZE(cost); i++)
        cost(i,i)=0;
    cost(0,K-1)=penalization;
    __bayesEnsembleNet->setCostMatrix(cost);
}

AutoParticlePickingQt::AutoParticlePickingQt(Micrograph *_m)
{
    __m                              = _m;
    __learn_particles_done           = false;
    __autoselection_done             = false;
    __modelRootName                  = "";
    __numThreads                     = 1;
    __auto_label                     = -1;
    __gray_bins                      = 8;
    __radial_bins                    = 16;
    __piece_xsize                    = 512;
    __highpass_cutoff                = 0.02;
    __penalization                   = 10;
    __minCost                        = 0;
    __output_scale                   = 1;
    __reduction                      = (int)std::pow(2.0, __output_scale);
    __particle_radius                = 0;
    __min_distance_between_particles = 0.5 * __particle_radius;
    __mask_size                      = 2 * __particle_radius;
    __piece_overlap                  = 2 * __particle_radius;
    __scan_overlap                   = (int)(2 * __particle_radius * 0.9);
    __learn_overlap                  = (int)(2 * __particle_radius * 0.5);
    __classNo                        = 3;
    __is_model_loaded                = false;
}

/* produceFeatures --------------------------------------------------------- */
void AutoParticlePickingQt::getFeatures(
    const Classification_modelQt &_model,
    std::vector < MultidimArray<double> > &_features)
{
    _features.clear();
    std::vector < std::vector<ParticleQt> > particles=
        _model.__training_particles;
    int vec_size=0;
    for (int i=0; i<particles.size(); i++)
        if (particles[i].size()!=0)
        {
            vec_size=(particles[i][0].vec).size();
            break;
        }

    for(int j = 0; j < __classNo; j++)
    {
        int imax = particles[j].size();
        //if we do not have false-positives, then we only use 2 classes.
        if (imax == 0)
        {
            _features.resize(__classNo - 1);
            break;
        }
        _features.push_back(*(new MultidimArray<double>));
        _features[j].initZeros(imax, vec_size);
        for(int i = 0; i < imax; i++)
        {
            particles[j][i].vec.setRow();
            _features[j].setRow(i, particles[j][i].vec);
        }
    }
}

/* produceClassesProbabilities --------------------------------------------- */
void AutoParticlePickingQt::getClassesProbabilities(
    const Classification_modelQt &_model, Matrix1D<double> &probabilities)
{
    double micrographsScanned = 0.0;
    double particlesMarked = 0.0;
    if (_model.__training_particles.size()>0)
        particlesMarked = _model.__training_particles[0].size();
    double falsePositives = 0.0;
    if (_model.__training_particles.size()>2)
        falsePositives = _model.__training_particles[2].size();
    probabilities.initZeros(__classNo);

    for (int i = 0; i < _model.__micrographs_number; i++)
        if (i<_model.__micrographs_scanned.size())
            micrographsScanned += _model.__micrographs_scanned[i];
    micrographsScanned=XMIPP_MAX(micrographsScanned,
                                 2*(particlesMarked+falsePositives));
    int NtrainingVectors=0;
    for (int i=0; i<_model.__training_particles.size(); i++)
        NtrainingVectors+=_model.__training_particles[i].size();
    micrographsScanned=XMIPP_MAX(micrographsScanned, NtrainingVectors);

    probabilities(0) = particlesMarked / micrographsScanned;
    probabilities(1) = (micrographsScanned - particlesMarked - falsePositives)/
                       micrographsScanned;
    probabilities(2) = falsePositives / micrographsScanned;
}

Classification_modelQt::Classification_modelQt(int _classNo, int _maxTrainingVectors)
{
    __maxTrainingVectors=_maxTrainingVectors;
    init(_classNo);
}

void Classification_modelQt::init(int _classNo)
{
    __classNo = _classNo;
    __training_particles.resize(__classNo);
    __micrographs_number = 0;
    __falsePositives.resize(0);
}

bool Classification_modelQt::addParticleTraining(const ParticleQt &p, int classIdx)
{
	if (__training_particles[classIdx].size()<__maxTrainingVectors)
	{
		__training_particles[classIdx].push_back(p);
		return true;
	}
	return false;
}

// ==========================================================================
// Section: xmipp_mark ======================================================
// ==========================================================================
/* Correct particles ------------------------------------------------------- */
void AutoParticlePickingQt::move_particle(int _idx)
{
    if (__autoselection_done)
    {
        int imax = __auto_candidates.size();
        int idxAuto=-1;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].idx == _idx)
            {
                idxAuto=i;
                break;
            }
        if (idxAuto==-1)
            return;

        double diffx=__auto_candidates[idxAuto].x - __m->coord(_idx).X;
        double diffy=__auto_candidates[idxAuto].y - __m->coord(_idx).Y;
        if (sqrt(diffx*diffx+diffy*diffy)>2*__particle_radius*0.1)
            __rejected_particles.push_back(__auto_candidates[idxAuto]);
        __auto_candidates[idxAuto].x = __m->coord(_idx).X;
        __auto_candidates[idxAuto].y = __m->coord(_idx).Y;
    }
}

/* Delete particles ---------------------------------------------------------*/
void AutoParticlePickingQt::delete_particle(int _idx)
{
    if (__autoselection_done)
    {
        // Look for this particle in the autolist
        int imax = __auto_candidates.size();
        int idxAuto=-1;
        for (int i = 0; i < imax; i++)
            if (__auto_candidates[i].idx == _idx)
            {
                idxAuto=i;
                break;
            }
        if (idxAuto==-1)
            return;

        __auto_candidates[idxAuto].status = 0;
        __auto_candidates[idxAuto].cost = -1;
        __rejected_particles.push_back(__auto_candidates[idxAuto]);
    }
}

void AutoParticlePickingQt::restrictSelection(float _cost)
{
    __minCost=_cost;
    int Nremaining=0;
    int imax = __auto_candidates.size();
    for (int i = 0; i < imax; i++)
        if (__auto_candidates[i].status == 1 && __auto_candidates[i].cost>_cost)
            Nremaining++;
    std::cout << Nremaining << " particles remain" << std::endl;
}

// ==========================================================================
// Section: I/O =============================================================
// ==========================================================================
/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out, const ParticleQt &_p)
{
    _out << _p.x      << " " << _p.y << " "
    << _p.idx    << " "
    << (int)_p.status << " "
    << _p.cost << " ";
    FOR_ALL_ELEMENTS_IN_MATRIX1D(_p.vec)
    _out << VEC_ELEM(_p.vec,i) << " ";
    _out << std::endl;
    return _out;
}

/* Read--------------------------------------------------------------------- */
void ParticleQt::read(std::istream &_in, int _vec_size)
{
    _in >> x >> y
    >> idx
    >> status
    >> cost;
    status -= '0';
    vec.resize(_vec_size);
    _in >> vec;
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out, const Classification_modelQt &_m)
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
    _out << "#Vector_size= " <<(_m.__training_particles[0][0].vec).size() << std::endl;
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
std::istream & operator >> (std::istream &_in, Classification_modelQt &_m)
{
    std::string dummy;
    int classNo, vec_size;
    _in >> dummy;
    _in >> dummy >> _m.__micrographs_number;
    _in >> dummy;

    _m.__micrographs_scanned.resize(_m.__micrographs_number);
    for(int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__micrographs_scanned[i];
    _in >> dummy;

    _m.__particles_picked.resize(_m.__micrographs_number);
    for(int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__particles_picked[i];
    _in >> dummy;

    _m.__falsePositives.resize(_m.__micrographs_number);
    for(int i = 0; i < _m.__micrographs_number; i++)
        _in >> _m.__falsePositives[i];
    _in >> dummy;

    _in >> dummy >> vec_size;
    _in >> dummy >> _m.__classNo;
    _m.__training_particles.resize(_m.__classNo);
    ParticleQt p;
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

void AutoParticlePickingQt::loadModels(const FileName &fn)
{
    __modelRootName = fn;

    // Load parameters
    std::string dummy;
    std::ifstream fh_params;
    fh_params.open((__modelRootName + ".param.txt").c_str());
    if (!fh_params)
        return;
    fh_params >> dummy >> __gray_bins
    >> dummy >> __radial_bins
    >> dummy >> __piece_xsize
    >> dummy >> __highpass_cutoff
    >> dummy >> __particle_radius
    >> dummy >> __min_distance_between_particles
    >> dummy >> __scan_overlap
    >> dummy >> __penalization
    ;
    fh_params.close();
    __piece_overlap=2*__particle_radius;
    // COSS:    if (__particle_radius>100) __output_scale=1;
    // COSS:   else __output_scale=0;
    __reduction=(int)std::pow(2.0, __output_scale);

    // Load the mask
    __mask.type = READ_MASK;
    __mask.fn_mask = __modelRootName + ".mask.xmp";
    __mask.generate_mask();
    __mask.get_binary_mask().setXmippOrigin();
    __mask_size = XSIZE(__mask.get_binary_mask()) * __reduction;

    classifyMask();

    // Load training vectors
    std::ifstream fh_training;
    fh_training.open((__modelRootName + ".training.txt").c_str());
    if (!fh_training)
        REPORT_ERROR(ERR_IO_NOTOPEN, (std::string)"AutoParticlePicking::write: Cannot open file " +
                     __modelRootName + ".training.txt" + " for input");

    fh_training >> __training_loaded_model;
    fh_training.close();

    // If we have already learnt or loaded the model, clear the training
    if (__learn_particles_done || __is_model_loaded)
        __training_model.clear();

    // Particles have not been learnt but loaded from a file
    __learn_particles_done = false;
    __is_model_loaded = true;
    __training_model.import_model(__training_loaded_model);
    __selection_model = __training_model;
    std::cout << "The model has been loaded..." << std::endl;
}

/* Save models ------------------------------------------------------------- */
void AutoParticlePickingQt::saveModels(const FileName &_fn_root) const
{
    FileName fn_root=_fn_root;
    if (fn_root=="")
        fn_root=__modelRootName;

    // Save the automatically selected particles
    saveAutoParticles();

    // Save the mask
    Image<double> save;
    typeCast(__mask.get_binary_mask(), save());
    save.write(fn_root + ".mask.xmp");

    // Save parameters
    std::ofstream fh_params;
    fh_params.open((fn_root + ".param.txt").c_str());
    if (!fh_params)
        REPORT_ERROR(ERR_IO_NOTOPEN, (std::string)"AutoParticlePicking::write: Cannot open file " +
                     fn_root + ".param.txt for output");
    fh_params << "gray_bins=                      " << __gray_bins                      << std::endl
    << "radial_bins=                    " << __radial_bins                    << std::endl
    << "piece_xsize=                    " << __piece_xsize                    << std::endl
    << "highpass=                       " << __highpass_cutoff                << std::endl
    << "particle_radius=                " << __particle_radius                << std::endl
    << "min_distance_between_particles= " << __min_distance_between_particles << std::endl
    << "particle_overlap=               " << __scan_overlap                   << std::endl
    << "penalization=                   " << __penalization                   << std::endl
    ;
    fh_params.close();

    // Save training vectors
    std::ofstream fh_training;
    fh_training.open((fn_root + ".training.txt").c_str());
    if (!fh_training)
        REPORT_ERROR(ERR_IO_NOTOPEN, (std::string)"AutoParticlePicking::write: Cannot open file " +
                     fn_root + ".training.txt for output");
    fh_training << __training_model << std::endl;
    fh_training.close();
    std::cout << "The model has been saved..." << std::endl;
}

/* Save particles ---------------------------------------------------------- */
void AutoParticlePickingQt::saveAutoParticles() const
{
    if (__autoselection_done && __auto_label!=-1)
        __m->write_coordinates(__auto_label, __minCost, __outputRoot+
                               "."+__modelRootName.removeDirectories()+".pos");
}

void AutoParticlePickingQt::saveAutoFeatureVectors(const FileName &fn) const
{
    std::ofstream fh_auto;
    fh_auto.open(fn.c_str());
    if (!fh_auto)
        REPORT_ERROR(ERR_IO_NOWRITE,fn);
    size_t Nauto=__auto_candidates.size();
    fh_auto << Nauto << " ";
    if (Nauto==0)
        fh_auto << "0\n";
    else
        fh_auto << VEC_XSIZE(__auto_candidates[0].vec) << std::endl;
    for (size_t i=0; i<Nauto; i++)
        fh_auto << __auto_candidates[i];
    fh_auto.close();
}

void AutoParticlePickingQt::loadAutoFeatureVectors(const FileName &fn)
{
    std::ifstream fh_auto;
    fh_auto.open(fn.c_str());
    if (!fh_auto)
        REPORT_ERROR(ERR_IO_NOTOPEN,fn);
    size_t Nauto;
    int vec_size;
    fh_auto >> Nauto >> vec_size;
    __auto_candidates.reserve(Nauto);
    ParticleQt p;
    for (size_t i=0; i<Nauto; i++)
    {
        p.read(fh_auto, vec_size);
        __auto_candidates.push_back(p);
    }
    fh_auto.close();
}
