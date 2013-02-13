/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2009)
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
#ifndef _PROG_VQ_VOLUMES
#define _PROG_VQ_VOLUMES

#include <parallel/xmipp_mpi.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/filters.h>
#include <data/polar.h>
#include <data/xmipp_fftw.h>
#include <data/histogram.h>
#include <data/numerical_tools.h>
#include <data/xmipp_program.h>
#include <data/symmetries.h>
#include <vector>

/**@defgroup VQforProjections Vector Quantization for Projections
   @ingroup ClassificationLibrary */
//@{
/** AssignedImage */
class CL3DAssignment
{
public:
	double score;   // Negative corrCodes indicate invalid particles
	double shiftx;
	double shifty;
	double shiftz;
	double psi;
	double tilt;
	double rot;
	size_t objId;

	/// Empty constructor
	CL3DAssignment();

	/// Read alignment parameters
	void readAlignment(const Matrix2D<double> &M);

	/// Copy alignment
	void copyAlignment(const CL3DAssignment &alignment);
};

/// Show
std::ostream & operator << (std::ostream &out, const CL3DAssignment& assigned);

/** CL3DClass class */
class CL3DClass {
public:
    // Centroid
    MultidimArray<double> P, Paux, bgMask;
    
    // Centroid in Fourier space
    MultidimArray< std::complex<double> > Pfourier;

    // Fourier transformer
    FourierTransformer transformer;

    // Auxiliary Fourier image
    MultidimArray< std::complex<double> > Ifourier;

    // Auxiliary Fourier image magnitude
    MultidimArray<double> IfourierMag, IfourierMagSorted;

    // Update for next iteration
    MultidimArray< std::complex<double> > Pupdate;

    // Update for next iteration
    MultidimArray< double > PupdateMask;

    // P in cylindrical coordinates around Z
    MultidimArray<double> PcylZ;

    // P in cylindrical coordinates around Y
    MultidimArray<double> PcylY;

    // P in cylindrical coordinates around X
    MultidimArray<double> PcylX;

    // Correlation aux
    CorrelationAux corrAux, corrAux2;

    // Volume Alignment aux
    VolumeAlignmentAux volAlignmentAux;

    // List of images assigned
    std::vector<CL3DAssignment> currentListImg;

    // List of images assigned
    std::vector<CL3DAssignment> nextListImg;

    // List of neighbour indexes
    std::vector<int> neighboursIdx;

public:
    /** Empty constructor */
    CL3DClass();

    /** Copy constructor */
    CL3DClass(const CL3DClass &other);

    /** Update projection. */
    void updateProjection(MultidimArray<double> &I, const CL3DAssignment &assigned, bool force=false);

    /** Transfer update */
    void transferUpdate();

    /** Sparse distance to centroid */
    void sparseDistanceToCentroid(MultidimArray<double> &I, double &avgK, double &stdK, double &L1distance);

    /** Compute the fit of the input image with this node.
        The input image is rotationally and traslationally aligned
        (2 iterations), to make it fit with the node. */
    void fitBasic(MultidimArray<double> &I, CL3DAssignment &result);

    /// Look for K-nearest neighbours
    void lookForNeighbours(const std::vector<CL3DClass *> listP, int K);
};

struct SDescendingClusterSort
{
     bool operator()(CL3DClass* const& rpStart, CL3DClass* const& rpEnd)
     {
          return rpStart->currentListImg.size() > rpEnd->currentListImg.size();
     }
};

/** Class for a CL3D */
class CL3D {
public:
	/// Number of images
	size_t Nimgs;

	/// Pointer to input metadata
	MetaData *SF;

    /// List of nodes
    std::vector<CL3DClass *> P;
    
public:
    /// Read Image
    void readImage(Image<double> &I, size_t objId, bool applyGeo) const;

    /// Initialize
    void initialize(MetaData &_SF,
    		        std::vector< MultidimArray<double> > &_codes0);
    
    /// Share assignments
    void shareAssignments(bool shareAssignment, bool shareUpdates);

    /// Share split assignment
    void shareSplitAssignments(Matrix1D<int> &assignment, CL3DClass *node1, CL3DClass *node2) const;

    /// Write the nodes
    void write(const FileName &fnRoot, int level) const;

    /** Look for a node suitable for this image.
        The image is rotationally and translationally aligned with
        the best node. */
    void lookNode(MultidimArray<double> &I, int oldnode,
    			  int &newnode, CL3DAssignment &bestAssignment);
    
    /** Transfer all updates */
    void transferUpdates();

    /** Quantize with the current number of codevectors */
    void run(const FileName &fnOut, int level);

    /** Clean empty nodes.
        The number of nodes removed is returned. */
    int cleanEmptyNodes();

    /** Split node */
    void splitNode(CL3DClass *node,
        CL3DClass *&node1, CL3DClass *&node2,
        std::vector<size_t> &finalAssignment) const;

    /** Split the widest node */
    void splitFirstNode();
};

/** CL3D parameters. */
class ProgClassifyCL3D: public XmippProgram {
public:
    /// Input selfile with the images to quantify
    FileName fnSel;
    
    /// Input selfile with initial codes
    FileName fnCodes0;

    /// Output rootname
    FileName fnOut;

    /// Symmetry file or code
    FileName fnSym;

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
    
    /// Sparsity factor (0<f<1; 1=drop all coefficients, 0=do not drop any coefficient)
    double sparsity;

    /// DWT Sparsity factor (0<f<1; 1=drop all coefficients, 0=do not drop any coefficient)
    double DWTsparsity;

    /// Clasify all images
    bool classifyAllImages;

    /// Maximum shift Z
    double maxShiftZ;

    /// Maximum shift Y
    double maxShiftY;

    /// Maximum shift X
    double maxShiftX;

    /// Maximum rot
    double maxRot;

    /// Maximum tilt
    double maxTilt;

    /// Maximum psi
    double maxPsi;

    /// Mask
    Mask mask;

    // Symmetry List
    SymList SL;

    /// MPI constructor
    ProgClassifyCL3D(int argc, char** argv);

    /// Destructor
    ~ProgClassifyCL3D();

    /// Read
    void readParams();
    
    /// Show
    void show() const;
    
    /// Usage
    void defineParams();
    
    /// Produce side info
    void produceSideInfo();
    
    /// Run
    void run();
public:
    // Selfile with all the input images
    MetaData SF;
    
    // Task distributor
    FileTaskDistributor *taskDistributor;

    // Object Ids
    std::vector<size_t> objId;

    // Structure for the classes
    CL3D vq;

    // Mpi node
    MpiNode *node;

    // Gaussian interpolator
    GaussianInterpolator gaussianInterpolator;

    // Image dimensions
    size_t Zdim, Ydim, Xdim;

	/// Noise in the images
    double sigma;
};
//@}
#endif
