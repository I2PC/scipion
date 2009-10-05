/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2009)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
#ifndef _PROG_VQ_PROJECTIONS
#define _PROG_VQ_PROJECTIONS

#include <data/selfile.h>
#include <data/polar.h>
#include <data/fftw.h>
#include <data/histogram.h>
#include <vector>
#include <mpi.h>

/**@defgroup VQforProjections
   @ingroup ReconsLibraryPrograms */
//@{
/** VQProjection class */
class VQProjection {
public:
    // Use correlation instead of correntropy
    bool useCorrelation;

    // Gaussian interpolator
    const GaussianInterpolator *gaussianInterpolator;

    // Mask of the background
    const Matrix2D<int> *mask;

    // Projection
    Matrix2D<double> P;
    
    // Update for next iteration
    Matrix2D<double> Pupdate;

    // Polar Fourier transform of the projection at full size
    Polar<std::complex <double> > polarFourierP;

    // Rotational correlation for best_rotation
    Matrix1D<double> rotationalCorr;

    // Local transformer for the best rotation
    XmippFftw local_transformer;
    
    // Plans for the best_rotation
    Polar_fftw_plans *plans;

    // List of images assigned
    std::vector<int> currentListImg;

    // List of images assigned
    std::vector<int> nextListImg;

    // Correlations of the next class members
    std::vector<double> nextClassCorr;

    // Correlations of the next non-class members
    std::vector<double> nextNonClassCorr;

    // Histogram of the correlations of the current class members
    histogram1D histClass;

    // Histogram of the correlations of the current non-class members
    histogram1D histNonClass;

    // List of neighbour indexes
    std::vector<int> neighboursIdx;
public:

    int sendMPI(int d_rank);
    
    int receiveMPI(int s_rank);

    /** Update projection. */
    void updateProjection(const Matrix2D<double> &I,
        double corrCode, int idx);

    /** Update non-projection */
    inline void updateNonProjection(double corrCode)
    {
        nextNonClassCorr.push_back(corrCode);
    }

    /** Transfer update */
    void transferUpdate();

    /** Compute different transforms */
    void computeTransforms();
    
    /** Compute the fit of the input image with this node.
        The input image is rotationally and translationally aligned
        (2 iterations), to make it fit with the node. */
    void fitBasic(Matrix2D<double> &I, double sigma, double &corrCode);

    /** Compute the fit of the input image with this node (check mirrors). */
    void fit(Matrix2D<double> &I, 
        double sigma, bool noMirror, double &corrCode, double &likelihood);

    /// Look for K-nearest neighbours
    void lookForNeighbours(const std::vector<VQProjection *> listP,
        double sigma, bool noMirror, int K);
    
    /// Show
    void show() const;
};

struct SDescendingClusterSort
{
     bool operator()(VQProjection* const& rpStart, VQProjection* const& rpEnd)
     {
          return rpStart->currentListImg.size() > rpEnd->currentListImg.size();
     }
};

/** Class for a VQ */
class VQ {
public:
    /// Mask for the background
    Matrix2D<int> mask;

    /// List of nodes
    std::vector<VQProjection *> P;
    
    /// Number of neighbours
    int Nneighbours;
    
    /// Minimum size of a node
    double PminSize;

    /// No mirror
    bool noMirror;

    /// Corr split
    bool corrSplit;

    /// Fast
    bool fast;

    /// Use correlation instead of correntropy
    bool useCorrelation;

    /// Maximum number of iterations
    int Niter;

    /// Noise in the images
    double sigma;

    // SelFile with the images
    SelFile *SF;

    std::vector< FileName > SFv;
    
    // Current assignment
    Matrix1D<int> currentAssignment;
    
    // Gaussian interpolator
    GaussianInterpolator gaussianInterpolator;
public:
    /// Initialize
    void initialize(SelFile &_SF, int _Niter,  int _Nneighbours,
        double _PminSize, std::vector< Matrix2D<double> > _codes0,
        int _Ncodes0, bool _noMirror, bool _corrSplit, 
        bool _useCorrelation, bool _fast, int rank);
    
    /// Write the nodes
    void write(const FileName &fnRoot) const;

    /** Look for a node suitable for this image.
        The image is rotationally and translationally aligned with
        the best node. */
    void lookNode(Matrix2D<double> &I, int idx, int oldnode, int &newnode,
        double &corrCode, double &likelihood);
    
    /** Update non codes. */
    void updateNonCode(Matrix2D<double> &I, int newnode);

    /** Transfer all updates */
    void transferUpdates();

    /** Quantize with the current number of codevectors */
    void run(const FileName &fnOut, int level, int rank);

    /** Clean empty nodes.
        The number of nodes removed is returned. */
    int cleanEmptyNodes();

    /** Split node */
    void splitNode(VQProjection *node,
        VQProjection *&node1, VQProjection *&node2, int rank,
	std::vector<int> &finalAssignment) const;

    /** Split the widest node */
    void splitFirstNode(int rank);
};

/** VQ parameters. */
class Prog_VQ_prm {
public:
    /// Input selfile with the images to quantify
    FileName fnSel;
    
    /// Input selfile with initial codes
    FileName fnCodes0;

    /// Output rootname
    FileName fnOut;

    /// Number of iterations
    int Niter;

    /// Initial number of code vectors
    int Ncodes0;

    /// Final number of code vectors
    int Ncodes;

    /// Number of neighbours
    int Nneighbours;

    /// Minimum size of a node
    double PminSize;
    
    /// Use Correlation instead of Correntropy
    bool useCorrelation;

    /// Fast
    bool fast;
    
    /// CorrSplit
    bool corrSplit;
    
    /// No mirror
    bool noMirror;

    /// Read
    void read(int argc, char** argv);
    
    /// Show
    void show() const;
    
    /// Usage
    void usage() const;
    
    /// Produce side info
    void produce_side_info(int rank);
    
    /// Run
    void run(int rank);
public:
    // Selfile with all the input images
    SelFile SF;
    
    // Structure for the classes
    VQ vq;
};
//@}
#endif
