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

#include "ml2d.h"

///******** Some macro definitions ****************
///Useful macro for thread iteration and work over all refno
#define FOR_ALL_THREAD_REFNO() \
int refno, load; \
while ((load = getThreadRefnoJob(refno)) > 0) \
    for (int i = 0; i < load; i++, refno = (refno + 1) % model.n_ref)
///Same macro as before, but without declaring refno and load
#define FOR_ALL_THREAD_REFNO_NODECL() \
while ((load = getThreadRefnoJob(refno)) > 0) \
    for (int i = 0; i < load; i++, refno = (refno + 1) % model.n_ref)

class ProgML2D;
// This structure is needed to pass parameters to the threads
typedef struct
{
    int thread_id;
    ProgML2D * prm;
}
structThreadTasks;

#define SPECIAL_ITER 0

void * doThreadsTasks(void * data);

/**@defgroup MLalign2D ml_align2d (Maximum likelihood in 2D)
   @ingroup ReconsLibrary */
//@{
/** MLalign2D parameters. */
class ProgML2D: public ML2DBaseProgram
{

public:
  bool no_iem;

    MultidimArray<int> mask, omask;
    /** Thread stuff */
    int threadTask;
    barrier_t barrier, barrier2, barrier3;
    pthread_t * th_ids;
    structThreadTasks * threads_d;

    /** New class variables, taken from old MAIN */
    double LL, sumfracweight;
    double wsum_sigma_noise, wsum_sigma_offset, sumw_allrefs;
    std::vector<double> sumw, sumw2, sumwsc, sumwsc2, sumw_mirror;
    std::vector<MultidimArray<double > > wsum_Mref;
    MultidimArray<int> Msignificant;
    MultidimArray<double> Mimg;
    std::vector<double> allref_offsets;
    std::vector<double> pdf_directions;
    std::vector<MultidimArray<double> > mref;
    std::vector<MultidimArray<std::complex<double> > > fref;
    std::vector<MultidimArray<std::complex<double > > > wsumimgs;
    int opt_refno, iopt_psi, iopt_flip;
    double trymindiff, opt_scale, bgmean, opt_psi;
    double fracweight, maxweight2, dLL;

    ///////// Some global variables just for thread visibility ////
    /** Taken from expectationSingleImage */
    std::vector<MultidimArray<std::complex<double> > > Fimg_flip, mysumimgs;
    std::vector<double> refw, refw2, refwsc2, refw_mirror, sumw_refpsi;
    double wsum_corr, sum_refw, sum_refw2, maxweight;
    double wsum_sc, wsum_sc2, wsum_offset, old_bgmean;
    double mindiff;
    size_t sigdim;
    int ioptx, iopty;
    std::vector<int> ioptx_ref, iopty_ref, ioptflip_ref;
    /** Taken from PreselectFastSignificant. */
    double pfs_mindiff;
    MultidimArray<double> pfs_maxweight;
    MultidimArray<double> pfs_weight;
    int pfs_count;

    /** Sum of squared amplitudes of the references */
    std::vector<double> A2;
    /** Sum of squared amplitudes of the experimental image */
    double Xi2;

    //These are for refno work assigns to threads
    size_t refno_index, refno_load_param;
    int refno_load, refno_count;
    // Which group does this image belong to in iteration 0 (generation of K references)
    int mygroup;
    /// Read arguments from command line
    void readParams();
    /// Params definition
    void defineParams();

public:
    ProgML2D();
    ///Show info at starting program
    virtual void show();
    ///Try to merge produceSideInfo1 and 2
    virtual void produceSideInfo();
    ///Try to merge produceSideInfo1 and 2
    virtual void produceSideInfo2();

    /// Calculate probability density distribution for in-plane transformations
    void calculatePdfInplane();

    /// Fill vector of matrices with all rotations of reference
    void rotateReference();

    /// Apply reverse rotations to all matrices in vector and fill new matrix with their sum
    void reverseRotateReference();

    /** Calculate which references have projection directions close to
        phi and theta */
    void preselectLimitedDirections(double &phi, double &theta);

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
    void awakeThreads(ThreadTask task, int start_refno, int load = 1);

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

    /// Perform an iteration
    virtual void iteration();

    /// Integrate over all experimental images
    virtual void expectation();

    /// Update all model parameters
    void maximizeModel(ModelML2D &model);

    /// Update all model parameters, adapted for IEM blocks use
    virtual void maximization();

    /// Correct references scale
    void correctScaleAverage();

    /// Add docfiledata to docfile
    virtual void addPartialDocfileData(const MultidimArray<double> &data, size_t first, size_t last);

    /// Write model parameters
    virtual void writeOutputFiles(const ModelML2D &model, OutputType outputType = OUT_FINAL);

    /// Read model from file
    virtual void readModel(ModelML2D &model, int block);
    /// Get base name based on fn_root and some number
    FileName getBaseName(String suffix = "", int number = -1);

    virtual void printModel(const String &msg, const ModelML2D & model);

};



//@}
#endif
