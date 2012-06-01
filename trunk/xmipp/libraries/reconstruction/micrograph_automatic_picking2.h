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

#ifndef MICROGRAPH_AUTOMATIC_PICKING_H
#define MICROGRAPH_AUTOMATIC_PICKING_H

#include <data/micrograph.h>
#include <data/mask.h>
#include <data/xmipp_program.h>
#include <reconstruction/fourier_filter.h>
#include <reconstruction/transform_geometry.h>
#include <classification/naive_bayes.h>
#include <classification/svm.h><

#include <data/xmipp_image.h>
#include <data/polar.h>
#include <data/normalize.h>
#include <data/basic_pca.h>

/// @defgroup AutomaticPicking Image denoising
/// @ingroup ReconsLibrary
//@{
/* Particle ---------------------------------------------------------------- */
class Particle2
{
public:
    FileName micrograph;       // Micrograph
    int x, y;                  // position in micrograph
    char status;               // rejected=0, selected=1 or moved=2
    Matrix1D<double> vec;      // vector of that particle
    double cost;               // Associated cost

    // Print
    friend std::ostream & operator << (std::ostream &_out, const Particle2 &_p);

    // Read
    void read(std::istream &_in, int _vec_size);
};

/* Classification model ---------------------------------------------------- */
class SVMClassifier
{
public:

	svm_parameter param;
    svm_problem prob;
    svm_model *model;
public:

    SVMClassifier();
    void SVMTrain(MultidimArray<double> &trainSet,MultidimArray<int> &lable);

};

/* Automatic particle picking ---------------------------------------------- */
/** Class to perform the automatic particle picking */
class AutoParticlePicking2
{
public:

    static const int NangSteps=120;

    Micrograph                *__m;
    Image<double>              microImage;
    PCAMahalanobisAnalyzer     pcaAnalyzer;
    SVMClassifier              classifier;
    FileName                   fn_micrograph;
    int                        __numThreads;
    Mask                       __mask;
    int                          piece_xsize;
    int                          particle_size;
    int                          particle_radius;
    int                          filter_num;
    int                          NPCA;
    int                          corr_num;
    double                       scaleRate;
    int                          NRsteps;
    bool                       __fast;
    bool                       __incore;
    std::vector<Particle2>      __auto_candidates;
    std::vector<Particle2>      __rejected_particles;

public:

    /// Empty constructor
    AutoParticlePicking2(const FileName &fn, Micrograph *_m,int size, int filterNum,int pcaNum,int corrNum);

    /// Destructor
    ~AutoParticlePicking2();

    /// Read the micrograph in memory
    void readMicrograph();

    //Extract the particles from the Micrograph
    void extractParticle(const int x, const int y, ImageGeneric & filter,
                         MultidimArray<double> &particleImage);
    //Convert an image to its polar form
    void convert2Polar(MultidimArray<double> &particleImage, MultidimArray<double> &polar);

    //Calculate the correlation of different polar channels
    void polarCorrelation(MultidimArray<double> &Ipolar,MultidimArray<double> &IpolarCorr);

    /// Produce the training dataset for different samples
    void buildVector(MultidimArray<double> &trainSet);

    /// Extract Invariant Features from a particle at x and y position
    void buildInvariant(MultidimArray<double> &invariantChannel,int x,int y,
                        const FileName &fnInvariantFeat);

    /// Extract different channels from particles and Non-Particles within a Micrograph
    void extractInvariant(const FileName &fnFilterBank,const FileName &fnInvariantFeat);

    /// Extract different filter channels from particles within a Micrograph
    void extractPositiveInvariant(const FileName &fnFilterBank,const FileName &fnInvariantFeat);

    /// Extract different filter channels from Non-Particles within a Micrograph
    void extractNegativeInvariant(const FileName &fnFilterBank,const FileName &fnInvariantFeat);

    /// Project a vector in PCA space
    void PCAProject(MultidimArray<double> &inputVec,int firstIndex);

    /// Train a PCA with negative and positive vectors
    void trainPCA(const FileName &fnPositiveFeat);

    /// Train SVM
    void trainSVM(const FileName &fnInvariantFeat);

    /// Save PCA model
    void savePCAModel(const FileName &fn_root) const;

    /// Load PCA model
    void loadPCAModel(const FileName &fn_root,int channelNum);


};

class ProgMicrographAutomaticPicking2: public XmippProgram
{
public:
    /// Micrograph filename
    FileName fn_micrograph;
    /// Model rootname
    FileName fn_model;
    /// Training coordinates
    FileName fn_train;
    /// Particle size
    int size;
    /// The number of filters
    int filter_num;
    /// The number of PCA components
    int NPCA;
    /// The number of correlation for a bank
    int corr_num;
    /// Mode
    String mode;
    /// Number of threads
    int Nthreads;
    /// Output rootname
    FileName fn_root;
    /// Fast
    bool fast;
    /// In core
    bool incore;
public:
    /// Read parameters
    void readParams();

    /// Show parameters
    void show();

    /// Define Parameters
    void defineParams();

    /** Run */
    void run();
};

//@}
#endif
