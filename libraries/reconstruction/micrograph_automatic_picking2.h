/***************************************************************************
*
* Authors:    Carlos Oscar            coss@cnb.csic.es (2011)
*      Vahid Abrishami    vabrishami@cnb.csic (2012)
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

class FeaturesThread;

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
    int particle_size, particle_radius, filter_num, proc_prec, NPCA, NRPCA, corr_num;
    int num_correlation, num_features, Nthreads, fast, NRsteps;

    SVMClassifier classifier, classifier2;
    PCAMahalanobisAnalyzer pcaAnalyzer;
    ProgImageRotationalPCA rotPcaAnalyzer;
    Point p1,p2;
    FeaturesThread * thread;
//    MetaData micList;
    Micrograph m, mPrev;

    Image<double> microImage, micrographStack, micrographStackPre, microImagePrev;

    std::vector<Particle2> auto_candidates;
    std::vector<Particle2> rejected_particles;
    std::vector<Particle2> accepted_particles;
    std::vector<Particle2> negative_candidates;
    std::vector<MDRow> micList;

    FileName fn_micrograph, fn_model, fnPCAModel, fnPCARotModel, fnAvgModel;
    FileName fnVector, fnSVMModel, fnSVMModel2, fnInvariant, fnParticles;

    double scaleRate;
    MultidimArray<double> convolveRes, filterBankStack, positiveParticleStack, negativeParticleStack;
    MultidimArray<double> positiveInvariatnStack, negativeInvariatnStack, autoFeatVec;
    MultidimArray<double> pcaModel, pcaRotModel, particleAvg, dataSet, dataSet1, classLabel;
    MultidimArray<double> classLabel1, labelSet, dataSetNormal;

public:
    /// Constructor
//    AutoParticlePicking2(int particle_size, int filter_num = 6, int corr_num = 2, int NPCA = 4,
//                         const FileName &model_name=NULL, const FileName &micsFn=NULL);
    AutoParticlePicking2(int particle_size, int filter_num = 6, int corr_num = 2, int NPCA = 4,
                         const FileName &model_name=NULL, const std::vector<MDRow> &vMicList = std::vector<MDRow>());

    AutoParticlePicking2();

    /// Destructor
    ~AutoParticlePicking2();

    void setSize(int pSize);

    /// Read micrograph from the file
    void readMic(const FileName &fn_micrograph, int keepPrev);

    void filterBankGenerator();

//    void batchBuildInvariant(const MetaData &MD);

    void batchBuildInvariant(const std::vector<MDRow> &MD);

//    void buildInvariant(const MetaData &MD);
    void buildInvariant(const std::vector<MDRow> &MD);
    void extractInvariant();

    void extractPositiveInvariant();

    void extractNegativeInvariant();

    void trainPCA();

    void add2Dataset(int flagNegPos);

//    void train(const MetaData &MD, bool corrFlag, int x, int y, int width, int height);

    void train(const std::vector<MDRow> &MD, bool corrFlag, int x, int y, int width, int height);

//    void correction(const MetaData &addedParticlesMD,const MetaData &removedParticlesMD);
    void correction(const std::vector<MDRow> &addedParticlesMD,const std::vector<MDRow> &removedParticlesMD);

    void add2Dataset(const MetaData &removedParticlesMD);

    void saveTrainingSet();

//    int automaticallySelectParticles(FileName fnmicrograph, int proc_prec, MetaData &md);

    int automaticallySelectParticles(FileName fnmicrograph, int proc_prec, std::vector<MDRow> &md);

    int automaticWithouThread(FileName fnmicrograph, int proc_prec, const FileName &fn);

    void saveAutoParticles(MetaData &md);

    void saveAutoParticles(std::vector<MDRow> &md);
    /// Define the parameters of the main program
    static void defineParams(XmippProgram * program);

    /// Read the parmaeters of the main program
    void readParams(XmippProgram * program);

    /// Read the micrograph in memory
    void readMicrograph();

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
                        int x,int y, int pre);

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

    /// Save the PCA basis and average for each channel
    void savePCAModel(const FileName &fn);

    /// Save training set into memory
    void saveTrainingSet(const FileName &fn_root);


    /// Load training set into the related array.
    void loadTrainingSet(const FileName &fn_root);

    /*
     * This method generates two different datasets. One for the
     * particles and non particles and the other one for the
     * particles and the false positives.
     */
    void generateTrainSet();

    /*
     * This method generate feature vectors for the candidates
     * of a micrographs which are obtained by a cross-correlation
     */
    void generateFeatVec(const FileName &fnmicrograph, int proc_prec,  std::vector<Particle2> &positionArray);

    /*
     * Read the next micrograph from the list of the micrographs
     */
    int readNextMic(FileName &fnmicrograph);

};

/**
 *  Structure to define random generation mode
 */
enum FeatureStatus
{
    TH_WAITING, TH_WORKING, TH_FINISHED, TH_ABORT
} ;

/** This class will compute the features calculation
 * in a separate thread.
 */
class FeaturesThread: public Thread
{
private:
    /* Condtions to be used to notifiy the thread to work(condIn)
     * and to know the thread has finished (condOut)
     */
    Condition condIn, condOut;
    bool waitingForResults; // Flag to know if the main thread is waiting
    FeatureStatus status;

public:
    AutoParticlePicking2 * picker;
    std::vector<Particle2> positionArray;
    FileName fnmicrograph;
    int proc_prec;

    FeaturesThread(AutoParticlePicking2 * picker);
    ~FeaturesThread();

    void setMicrograph(const FileName &fnMic, int proc_prec);
    void run();
    /* This function should be called from outside to wait for results. */
    void waitForResults();
    /* Notify the thread to start working another micrograph */
    void workOnMicrograph(const FileName &fnMic, int proc_prec);
    /* Stop the work because the micrograph that the thread is working is not
     * the next one the user has clicked.
     */
    void cancelWork();
    /* Call the picker generateFeatVec */
    void generateFeatures();
}
; //class FeaturesThread

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
