/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
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

#ifndef _MLALIGN2D_H
#define _MLALIGN2D_H

#include <data/fftw.h>
#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/ctf.h>
#include <data/threads.h>
#include <pthread.h>
#include <vector>
#include <sys/time.h>

#define SIGNIFICANT_WEIGHT_LOW 1e-8
#define SMALLANGLE 1.75
#define DATALINELENGTH 12

class Prog_MLalign2D_prm;

#define FOR_ALL_THREAD_REFNO() \
int refno, load; \
while ((load = getThreadRefnoJob(refno)) > 0) \
    for (int i = 0; i < load; i++, refno = (refno + 1) % model.n_ref)

#ifdef TIMING
//testing time...

#define TB_TOTAL 18
typedef enum TimingBlocks { ITER, ITER_E, ITER_M, E_RR, E_PRE, E_FOR, E_RRR, E_OUT, FOR_F1, FOR_PFS, FOR_ESI, FOR_F2, ESI_E1, ESI_E2TH, ESI_E3, ESI_E4, ESI_E5, ESI_E6TH } TimingBlocks;

class JMTimer
{
public:
    ///Some timing stuff
    timeval start_times[TB_TOTAL];
    int counts[TB_TOTAL];
    long int times[TB_TOTAL];
    char * tags[];
    //timeval start_time;
    timeval end_time;

    JMTimer()
    {
        clear();
    }

    void clear()
    {
        for (int i = 0; i < 25; i++)
            counts[i] = times[i] = 0;
    }

    void tic(TimingBlocks tb)
    {
        int i = (int)tb;
        gettimeofday(start_times + i, NULL);
        counts[i]++;
    }

    int toc(TimingBlocks tb)
    {
        int i = (int)tb;
        gettimeofday(&end_time, NULL);
        times[i] += (end_time.tv_sec - start_times[i].tv_sec) * 1000000 +
                       (end_time.tv_usec - start_times[i].tv_usec);
    }

    void printTimes(bool doClear)
    {
        char * tags[] = { "ITER", "ITER_E", "ITER_M", "E_RR", "E_PRE", "E_FOR", "E_RRR", "E_OUT", "FOR_F1", "FOR_PFS", "FOR_ESI", "FOR_F2", "ESI_E1", "ESI_E2TH", "ESI_E3", "ESI_E4", "ESI_E5", "ESI_E6TH"};

        for (int i = 0; i < TB_TOTAL; i++)
        {
            std::cout << tags[i] << " took: " << times[i] / counts[i] << " microseconds" << std::endl;
        }

        if (doClear)
            clear();
    }
};
#endif

//threadTask constants
#define TH_EXIT 0
#define TH_ESI_REFNO 1
#define TH_ESI_UPDATE_REFNO 2
#define TH_RR_REFNO 3
#define TH_RRR_REFNO 4
#define TH_PFS_REFNO 5

//output types constants
#define OUT_BLOCK 0
#define OUT_ITER  1
#define OUT_FINAL 2


// This structure is needed to pass parameters to the threads
typedef struct{
    int thread_id;
    Prog_MLalign2D_prm * prm;
} structThreadTasks;

void * doThreadsTasks(void * data);

class Model_MLalign2D
{
public:
    /** Number of reference images */
    int n_ref;
    /** References images */
    std::vector < ImageXmippT<double> > Iref;
    /** Sigma value for expected pixel noise */
    double sigma_noise;
    /** sigma-value for origin offsets */
    double sigma_offset;
    /** Vector containing estimated fraction for each model */
    std::vector<double> alpha_k;
    /** Vector containing estimated fraction for mirror of each model */
    std::vector<double> mirror_fraction;
    /** Average scales for each refno from do_norm */
    std::vector<double> scale;
    /** Sums of weights */
    double sumw_allrefs, sumw_allrefs2;
    /** Average height of the probability distribution at its maximum */
    double avePmax;
    /** the Log Likelihood */
    double LL;

    /** Dimension of images */
    int dim;
    /** Algorithmic variants */
    bool do_student, do_student_sigma_trick, do_norm;

    Model_MLalign2D();
    Model_MLalign2D(int n_ref);

    void initData();
    void setSize();
    void combineModel(Model_MLalign2D model, int sign);
    void addModel(Model_MLalign2D model);
    void substractModel(Model_MLalign2D model);

    double get_sumw(int refno);
    double get_sumw_mirror(int refno);
    double get_sumwsc(int refno);
    Matrix2D<double> get_wsum_Mref(int refno);
    double get_sumw_allrefs();
    double get_wsum_sigma_offset();
    double get_wsum_sigma_noise();
    double get_sumfracweight();
    void updateSigmaOffset(double wsum_sigma_offset);
    void updateSigmaNoise(double wsum_sigma_noise);
    void updateAvePmax(double sumfracweight);
    void updateFractions(int refno, double sumw, double sumw_mirror, double sumw_allrefs);
    void updateScale(int refno, double sumwsc, double sumw);
    ///Just for debugging now
    void print();
};//close class Model_MLalign2D

/**@defgroup MLalign2D ml_align2d (Maximum likelihood in 2D)
   @ingroup ReconsLibraryPrograms */
//@{
/** MLalign2D parameters. */
class Prog_MLalign2D_prm
{
public:
    /** Filenames reference selfile/image, fraction docfile & output rootname */
    FileName fn_sel, fn_ref, fn_root, fn_frac, fn_sig, fn_doc, fn_oext, fn_scratch, fn_control;
    /** Command line */
    std::string cline;
    double sigma_noise2;
    /** Flag for checking mirror images of all references */
    bool do_mirror;
    /** Flag whether to fix estimates for model fractions */
    bool fix_fractions;
    /** Flag whether to fix estimate for sigma of origin offset */
    bool fix_sigma_offset;
    /** Flag whether to fix estimate for sigma of noise */
    bool fix_sigma_noise;
    /** Starting iteration */
    int istart;
    /** Current iteration */
    int iter;
    /** Number of iterations to be performed */
    int Niter;
    /** dimension of the images */
    int dim, dim2, hdim;
    double ddim2;
    /** Number of steps to sample in-plane rotation in 90 degrees */
    int nr_psi;
    /** Number of operations in "flip-array" (depending on do_mirror) */
    int nr_flip;
    /** Sampling rate for in-plane rotation */
    float psi_step;
    /** Total degrees in FOR_ALL_ROTATIONS */
    double psi_max;
    /** Total number of no-mirror rotations in FOR_ALL_FLIPS */
    int nr_nomirror_flips;
    /** Total number of experimental images */
    int nr_exp_images;
    /** Sum of squared amplitudes of the references */
    std::vector<double> A2;
    /** Sum of squared amplitudes of the experimental image */
    double Xi2;
    /** Verbose level:
        1: gives progress bar (=default)
        0: gives no output to screen at all */
    int verb;
    /** Stopping criterium */
    double eps;
    /** SelFiles with experimental and reference images */
    SelFile SF, SFr;
    /** vector for flipping (i.e. 90/180-degree rotations) matrices */
    std::vector<Matrix2D<double> > F;
    /** Vector for images to hold references (new & old) */
    std::vector < ImageXmippT<double> > Iold;
    /** Matrices for calculating PDF of (in-plane) translations */
    Matrix2D<double> P_phi, Mr2;
    /** Masks for rotated references */
    Matrix2D<int> mask, omask;
    /** Fast mode */
    bool fast_mode;
    /** Fast mode */
    double C_fast;
    /** Offsets for limited translations */
    std::vector<Matrix1D<double> > Vtrans;
    /** Start all optimal offsets from zero values */
    bool zero_offsets;
    /** Limit orientational searches */
    bool limit_rot;
    /** Limited search range for projection directions */
    double search_rot;
    /** Save memory options */
    bool save_mem1, save_mem2, save_mem3;
    /** Output document file with output optimal assignments*/
    DocFile DFo;
    /** Vectors to store old phi and theta for all images */
    std::vector<float> imgs_oldphi, imgs_oldtheta;
    /** Flag for using ML3D */
    bool do_ML3D;
    /** Flag for generation of initial models from random subsets */
    bool do_generate_refs;
    /** Vector to store optimal origin offsets (if not written to disc) */
    std::vector<std::vector<double> > imgs_offsets;
    /** For initial guess of mindiff */
    double trymindiff_factor;

    /// Students t-distribution
    /** Use t-student distribution instead of normal one */
    bool do_student;
    /** Degrees of freedom for the t-student distribution */
    double df, df2, dfsigma2;
    /** Do sigma-division trick in student-t*/
    bool do_student_sigma_trick;

    /// Re-normalize internally
    /** Flag to refine normalization of each experimental image */
    bool do_norm;
    /** Grey-scale correction values */
    std::vector<double> imgs_scale, imgs_bgmean, imgs_trymindiff;
    /** Optimal refno from previous iteration */
    std::vector<int> imgs_optrefno;
    /** Overall average scale (to be forced to one)*/
    double average_scale;
 
    /** Thread stuff */
    int threads, threadTask;
    barrier_t barrier, barrier2;
    pthread_t * th_ids;
    structThreadTasks * threads_d;

    /** debug flag */
    int debug;

    /** New class variables */
    double LL, sumfracweight;
    double wsum_sigma_noise, wsum_sigma_offset, sumw_allrefs;
    std::vector<double> sumw, sumw2, sumwsc, sumwsc2, sumw_mirror;
    std::vector<Matrix2D<double > > wsum_Mref;
    std::vector<double> conv;
    Matrix2D<int> Msignificant;
    Matrix2D<double> Mimg;
    std::vector<double> allref_offsets;
    std::vector<double> pdf_directions;
    std::vector<Matrix2D<double> > mref;
    std::vector<Matrix2D<std::complex<double> > > fref;
    std::vector<Matrix2D<std::complex<double > > > wsumimgs;
    int opt_refno, iopt_psi, iopt_flip;
    double trymindiff, opt_scale, bgmean, opt_psi;
    double fracweight, maxweight2, dLL;
    std::vector<Matrix1D<double> > docfiledata;

    /** From expectationSingleImage */
    std::vector<Matrix2D<std::complex<double> > > Fimg_flip, mysumimgs;
    std::vector<double> refw, refw2, refwsc2, refw_mirror, sumw_refpsi;
    double wsum_corr, sum_refw, sum_refw2, maxweight;
    double wsum_sc, wsum_sc2, wsum_offset, old_bgmean;
    double mindiff;
    int sigdim;
    int ioptx, iopty;
    std::vector<int> ioptx_ref, iopty_ref, ioptflip_ref;
    std::vector<double> maxw_ref;
    //These are for refno work assigns to threads
    int refno_index, refno_count, refno_load, refno_load_param;

    //Some incremental stuff
    /** Model */
    Model_MLalign2D model, *current_model;
    std::vector<int> img_blocks;
    //Number of blocks for IEM
    int blocks;
    //Current processing block
    int current_block;
    //Dont randomize for tests
    bool randomize;

#ifdef TIMING
    JMTimer timer;
#endif

public:
    /// Read arguments from command line
    void read(int argc, char **argv, bool ML3D = false);

    /// Show
    void show(bool ML3D = false);

    /// Usage for ML mode
    void usage();

    /// Extended Usage
    void extendedUsage(bool ML3D = false);

    /// Setup lots of stuff
    void produceSideInfo();

    /// Generate initial references from random subset averages
    void generateInitialReferences();

    /** Read reference images in memory & set offset vectors
        (This produce_side_info is Selfile-dependent!) */
    void produceSideInfo2(int nr_vols = 1);

    /// Calculate probability density distribution for in-plane transformations
    void calculatePdfInplane();

    /// Fill vector of matrices with all rotations of reference
    void rotateReference();

    /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
    void reverseRotateReference();

    /** Calculate which references have projection directions close to
        phi and theta */
    void preselectLimitedDirections(float &phi, float &theta);

    /** Pre-calculate which model and phi have significant probabilities
       without taking translations into account! */
    void preselectFastSignificant();

    /// ML-integration over all (or -fast) translations
    void expectationSingleImage(Matrix1D<double> &opt_offsets);

    /*** Threads functions */
    /// Create working threads
    void createThreads();

    /// Exit threads and free memory
    void destroyThreads();

    /// Assign refno jobs to threads
    int getThreadRefnoJob(int &refno);

    /// Awake threads for different tasks
    void awakeThreads(int task, int start_refno, int load = 1);

    /// Thread code to parallelize refno loop in rotateReference
    void doThreadRotateReferenceRefno();

    ///Thread code to parallelize refno loop in reverseRotateReference
    void doThreadReverseRotateReferenceRefno();

    /// Thread code to parallelize refno loop in preselectFastSignificant
    void doThreadPreselectFastSignificantRefno();

    /// Thread code to parallelize refno loop in expectationSingleImage
    void doThreadExpectationSingleImageRefno();

    /// Thread code to parallelize update loop in ESI
    void doThreadESIUpdateRefno();

    /// Integrate over all experimental images
    void expectation();

    /// Update all model parameters
    void maximization(Model_MLalign2D &model, int refs_per_class=1);

    /// Correct references scale
    void correctScaleAverage(int refs_per_class=1);

    /// check convergence
    bool checkConvergence();

    /// Write Docfile
    void writeDocfile(FileName fn_base);

    /// Write model parameters
    void writeOutputFiles(Model_MLalign2D model, int outputType = OUT_FINAL);

    /// Read model from file
    void readModel(Model_MLalign2D &model, FileName fn_base);

    /// Get base name based on fn_root and some number
    FileName getBaseName(std::string suffix = "", int number = -1);

};



//@}
#endif
