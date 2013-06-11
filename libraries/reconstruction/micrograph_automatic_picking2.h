/***************************************************************************
*
* Authors:    Carlos Oscar            coss@cnb.csic.es (2011)
* 			  Vahid Abrishami		  vabrishami@cnb.csic (2012)
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

#ifndef MICROGRAPH_AUTOMATIC_PICKING2_H
#define MICROGRAPH_AUTOMATIC_PICKING2_H

#include <data/micrograph.h>
#include <data/mask.h>
#include <data/xmipp_program.h>
#include <data/transform_geometry.h>
#include <reconstruction/fourier_filter.h>
#include <classification/naive_bayes.h>
#include <classification/svm_classifier.h>
#include <classification/knn_classifier.h>
#include <reconstruction/image_rotational_pca.h>

#include <data/xmipp_image.h>
#include <data/polar.h>
#include <data/normalize.h>
#include <data/basic_pca.h>
#include <data/morphology.h>

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
    MultidimArray<double> vec; // vector of that particle
    double cost;               // Associated cost

    // Print
    friend std::ostream & operator << (std::ostream &_out, const Particle2 &_p);

    // Read
    void read(std::istream &_in, int _vec_size);
};

/* Automatic particle picking ---------------------------------------------- */
/** Class to perform the automatic particle picking */
class AutoParticlePicking2
{
public:

    static const int NangSteps=120;

    Micrograph                  *__m;
    Micrograph                   m;
    Image<double>                microImage;
    PCAMahalanobisAnalyzer       pcaAnalyzer;
    ProgImageRotationalPCA       rotPcaAnalyzer;
    SVMClassifier                classifier;
    SVMClassifier                classifier2;
    FileName                     fn_micrograph;
    int                          piece_xsize;
    int                          particle_size;
    int                          particle_radius;
    int                          filter_num;
    int                          proc_prec;
    int                          NPCA;
    int                          NRPCA;
    int                          corr_num;
    int                          num_correlation;
    int                          num_features;
    int 					     Nthreads;
    int 						 fast;
    double                       scaleRate;
    int                          NRsteps;

    MultidimArray<double>        convolveRes;
    MultidimArray<double>        filterBankStack;
    MultidimArray<double>        positiveParticleStack;
    MultidimArray<double>        negativeParticleStack;
    MultidimArray<double>        positiveInvariatnStack;
    MultidimArray<double>        negativeInvariatnStack;
    MultidimArray<double>        pcaModel;
    MultidimArray<double>        pcaRotModel;
    MultidimArray<double>        particleAvg;
    MultidimArray<double>        dataSet;
    MultidimArray<double>        dataSet1;
    MultidimArray<double>        classLabel;
    MultidimArray<double>        classLabel1;
    MultidimArray<double>        labelSet;

    std::vector<Particle2>       auto_candidates;
    std::vector<Particle2>       rejected_particles;
    std::vector<Particle2>       accepted_particles;
    Image<double>                micrographStack;

    FileName                     fn_model;
    FileName                     fnPCAModel;
    FileName                     fnPCARotModel;
    FileName                     fnAvgModel;
    FileName                     fnVector;
    FileName                     fnSVMModel;
    FileName                     fnSVMModel2;
public:

    /// Constructor
    AutoParticlePicking2(int particle_size, int filter_num = 6, int corr_num = 2, int NPCA = 4, const FileName &model_name=NULL);

    AutoParticlePicking2();

    /// Destructor
    ~AutoParticlePicking2();

    void setSize(int pSize);

    /// Read micrograph from the file
    void readMic(FileName fn_micrograph);

    void filterBankGenerator();

    void batchBuildInvariant(MetaData MD);

    void buildInvariant(MetaData MD);

    void extractInvariant();

    void extractPositiveInvariant();

    void extractNegativeInvariant();

    void trainPCA();

    void add2Dataset(int flagNegPos);

    void train(MetaData MD, bool corrFlag);

    void correction(MetaData addedParticlesMD,MetaData removedParticlesMD);

    void add2Dataset(MetaData removedParticlesMD);

    void saveTrainingSet();

    int automaticallySelectParticles(FileName fnmicrograph, int proc_prec, MetaData &md);

    void saveAutoParticles(MetaData &md);

    /// Define the parameters of the main program
    static void defineParams(XmippProgram * program);

    /// Read the parmaeters of the main program
    void readParams(XmippProgram * program);

    /// Read the micrograph in memory
    void readMicrograph();

    void produceSideInfo(Micrograph *m);

    //Check the distance between a point and positive samples in micrograph
    bool checkDist(Particle2 &p);

    /*
     * This method extracts statics features such as
     * average, variance and quantiles of a particle.
     */
    void extractStatics(MultidimArray<double> &inputVec,
                        MultidimArray<double> &features);

    /// Convert an image to its polar form
    void convert2Polar(MultidimArray<double> &particleImage,
                       MultidimArray<double> &polar);

    /// Calculate the correlation of different polar channels
    void polarCorrelation(MultidimArray<double> &Ipolar,
                          MultidimArray<double> &IpolarCorr);

    /// Convolve the micrograph with the different templates
    void applyConvolution(bool fast);

    /// Project a vector on one pca basis
    double PCAProject(MultidimArray<double> &pcaBasis,
                      MultidimArray<double> &vec);

    /*
     * Extract the particle from the micrograph at x and y
     * and normalize it if the flag normal is true.
     */
    void extractParticle(const int x, const int y,
                         MultidimArray<double> &filter,
                         MultidimArray<double> &particleImage,
                         bool normal);

    /* Extract non partiles from the micrograph and put these
     * positions in a vector. Later we select some of them in a
     * random manner in order to train the classifier with negative
     * samples.
     */
    void extractNonParticle(std::vector<Particle2> &negativePosition);

    /*
     * Extract the invariants from the particles and non particles
     * in the micrograph.The invariants are the correlations between
     * different channels in polar form.
     */
    void extractInvariant(const FileName &fnInvariantFeat,
                          const FileName &fnParticles,
                          bool avgFlag);

    /*
     * This method extracts the invariants from the particles
     * in a micrograph. If the flag is true then the average
     * of the particles is also computed.The particles are
     * saved in fnParticles.
     */
    void extractPositiveInvariant(const FileName &fnInvariantFeat,
                                  const FileName &fnParticles,
                                  bool avgFlag);

    /*
     * This method extracts the invariants from the non particles
     * in a micrograph.The non particles are saved in fnParticles.
     */
    void extractNegativeInvariant(const FileName &fnInvariantFeat,
                                  const FileName &fnParticles);

    /*
     * This method is used in order to extract all the features
     * from an input (Can be particle or non particle)
     */
    void buildVector(MultidimArray<double> &inputVec,
                     MultidimArray<double> &staticVec,
                     MultidimArray<double> &featureVec,
                     MultidimArray<double> &pieceImage);

    /*Extract the invariants from just one particle at x,y
     *The invariants are the correlations between different channels
     *in polar form.
     */
    void buildInvariant(MultidimArray<double> &invariantChannel,
                        int x,int y);

    /*
     * This method does a convolution in order to find an approximation
     * about the place of the particles.
     */
    void buildSearchSpace(std::vector<Particle2> &positionArray,bool fast);

    /*
     * This method is used in order to train an support vector
     * machine. It receives the file which contains the dataset
     * and also which classifier we want to train. In this case
     * we have two classifiers.
     */
    void trainSVM(const FileName &fnModel,int numClassifier);

    /*
     * This method uses the extracted invariants from the
     * particles in order to find some pca basis to reduce
     * the number of features.
     */
    void trainPCA(const FileName &fnPositiveFeat);

    /*
    * This method is used to generate some pca basis according
    * to different rotation of the template.
    */
    void trainRotPCA(const FileName &fnAvgModel,const FileName &fnPCARotModel);

    /*
     * Extracts all the features related to particles and non
     * particles and put them in an array with the related
     * labels.
     */
    void add2Dataset(const FileName &fnInvariantFeat,
                     const FileName &fnParticles,int lable);

    /*
     * It is the same as previous one but it does not extract
     * the features. It just puts the data in an array for dataset.
     */
    void add2Dataset();

    /* Normalize the data of a dataset according to
     * a and b.
     */
    void normalizeDataset(int a,int b);

    /// Save automatically selected particles
    int saveAutoParticles(const FileName &fn) const;

    /*
     * In Semi-Automatic step, we save all the feature
     * vectors in order to have them to retrain the
     * classifier.
     */
    void saveAutoVectors(const FileName &fn);

    /// Save the extracted features for both rejected and found features
    void saveVectors(const FileName &fn);

    /// Save the PCA basis and average for each channel
    void savePCAModel(const FileName &fn_root);

    /// Save training set into memory
    void saveTrainingSet(const FileName &fn_root);

    /*
     * In Semi-Automatic step, we save all the feature
     * vectors in order to have them to retrain the
     * classifier.
     */
    void loadAutoVectors(const FileName &fn);

    /// Load training set into the related array.
    void loadTrainingSet(const FileName &fn_root);

    /// Load the features for particles and non-particles (from the supervised)
    void loadVectors(const FileName &fn);

    /// Select particles from the micrograph in an automatic way
    int automaticallySelectParticles(bool use2Classifier);

    /*
     * This method generates two different datasets. One for the
     * particles and non particles and the other one for the
     * particles and the false positives.
     */
    void generateTrainSet();
};

struct AutoPickThreadParams
{
	AutoParticlePicking2 *autoPicking;
	std::vector<Particle2> positionArray;
	bool use2Classifier;
	int idThread;
	int Nthreads;
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
    /// Mode
    String mode;
    /// Number of threads
    /// Output rootname
    FileName fn_root;
	AutoParticlePicking2 *autoPicking;
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
