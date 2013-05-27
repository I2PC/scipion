/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2008)
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
#ifndef _PROG_ANGULAR_ASSIGN_FOR_TILTSERIES
#define _PROG_ANGULAR_ASSIGN_FOR_TILTSERIES

#include <vector>
#include <data/matrix1d.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/xmipp_program.h>
#include <pthread.h>

/**@defgroup AngularAssignTiltSeries angular_assign_for_tilt_series
             (Simultaneous assignment)
   @ingroup ReconsLibrary */
//@{

/** Compute the affine transformation between two images. */
double computeAffineTransformation(const MultidimArray<unsigned char> &I1,
                                   const MultidimArray<unsigned char> &I2, int maxShift, int maxIterDE,
                                   const FileName &fn_affine,
                                   Matrix2D<double> &A12, Matrix2D<double> &A21, bool show,
                                   double thresholdAffine, bool globalAffine, bool isMirror,
                                   bool checkRotation);

/** Landmark class.
    A landmark is a position (x,y) and the index of the image
    from which this landmark has been taken.*/
class Landmark
{
public:
    double x;
    double y;
    int    imgIdx;
    Landmark& operator=(const Landmark &rhs)
    {
        x=rhs.x;
        y=rhs.y;
        imgIdx=rhs.imgIdx;
        return *this;
    }
    int operator<(const Landmark &rhs) const
    {
        return imgIdx<rhs.imgIdx;
    }
};

/** A landmark chain is simply a vector of landmarks. */
typedef std::vector<Landmark> LandmarkChain;

/* Forward prototype */
class Alignment;

/** This is the main class */
class ProgTomographAlignment: public XmippProgram
{
public:
    /// MetaData File with all images
    FileName fnSel;

    /// MetaData File with all images at the original scale
    FileName fnSelOrig;

    /// Output root
    FileName fnRoot;

    /// Number of threads to use for parallel computing
    int numThreads;

    /// Look for local affine transformation
    bool globalAffine;

    /// Use critical points
    bool useCriticalPoints;

    /// Number of critical points
    size_t Ncritical;

    /// Sequence Length
    size_t seqLength;

    /** Add blind landmarks.
        Set to -1 for no blind landmarks */
    size_t blindSeqLength;

    /// Grid samples
    int gridSamples;

    /// Maxshift percentage
    double maxShiftPercentage;

    /// Maximum number of iterations in Differential Evolution
    int maxIterDE;

    /// Maximum rotation
    double psiMax;

    /// Maximum Step for chain refinement
    int maxStep;

    /// Delta rot
    double deltaRot;

    /// Local refinement size
    double localSize;

    /// Optimize tilt angle of tilt axis
    bool optimizeTiltAngle;

    /// This tilt series comes from a capillar
    bool isCapillar;

    /// Don't normalize the corrected images
    bool dontNormalize;

    /// Difficult
    bool difficult;

    /// Correlation threshold for a good landmark
    double corrThreshold;

    /// Show affine transformations
    bool showAffine;

    /// Threshold affine
    double thresholdAffine;

    /// Identify outlier micrographs
    double identifyOutliersZ;

    /// Do not identify outlier micrographs
    bool doNotIdentifyOutliers;

    /// Pyramid
    int pyramidLevel;

    /// Last step to run (-1=run them all)
    int lastStep;

    // iteration counter as a progress measure
    int iteration;

    /// Destructor
    ~ProgTomographAlignment();

    /// Usage
    void defineParams();

    /// Read parameters from argument line
    void readParams();

    /// Show parameters
    void show();

    /// Produce side info
    void produceSideInfo();

    /// Compute affine transformations
    void computeAffineTransformations(bool globalAffineToUse);

    /// Identify outliers
    void identifyOutliers(bool mark);

    /// Produce information from landmarks
    void produceInformationFromLandmarks();

    /** Refine landmark.
        ii is the index of the original image, jj is the index in the
        image at which the landmark is being refined. rii and rjj are
        the corresponding landmark positions in both images.
        
        The function returns whether the landmark is accepted or not. */
    bool refineLandmark(int ii, int jj, const Matrix1D<double> &rii,
                        Matrix1D<double> &rjj, double &maxCorr, bool tryFourier) const;


    /** Refine landmark.
        The same as the previous function but an image is provided
        as pattern (ii) instead of an index and a position. */
    bool refineLandmark(const MultidimArray<double> &pieceii, int jj,
                        Matrix1D<double> &rjj, double actualCorrThreshold,
                        bool reversed, double &maxCorr) const;

    /** Refine chain. */
    bool refineChain(LandmarkChain &chain, double &corrChain);

    /// Generate landmark set using a grid
    void generateLandmarkSet();

    /// Write landmark set
    void writeLandmarkSet(const FileName &fnLandmark) const;

    /// Write affine transformations
    void writeTransformations(const FileName &fnTransformations) const;

    /// Read landmark set
    void readLandmarkSet(const FileName &fnLandmark);

    /// Remove Outliers
    void removeOutlierLandmarks(const Alignment &alignment);

    /// Align images
    void alignImages(const Alignment &alignment);

    /// Run: do the real work
    void run();

public:
    // MetaData with the input images
    MetaData SF;

    // MetaData with the input images
    MetaData SForig;

    // List of image pointers
    std::vector < MultidimArray<unsigned char> *> img;

    // List of mask pointers
    std::vector < MultidimArray<unsigned char> *> maskImg;

    // Index of the image closest to 0 degrees in tilt
    int iMinTilt;

    // List of tilt angles
    std::vector < double > tiltList;

    // List of names
    std::vector < std::string > name_list;

    // Average forward patch correlation
    Matrix1D<double> avgForwardPatchCorr;

    // Average backward patch correlation
    Matrix1D<double> avgBackwardPatchCorr;

    // Outlier
    Matrix1D<int> isOutlier;

    // Total number of images
    int Nimg;

    // Element [i][j] is the transformation of coordinates of i into
    // coordinates of j
    std::vector< std::vector< Matrix2D<double> > > affineTransformations;

    // List of affine costs
    std::vector < double > correlationList;

    // Landmark matrix (X component)
    Matrix2D<double> allLandmarksX;

    // Landmark matrix (Y component)
    Matrix2D<double> allLandmarksY;

    // Set of all landmarks seen in image i
    std::vector< std::vector<int> > Vseti;

    // Set of all images in which the landmark j is seen
    std::vector< std::vector<int> > Vsetj;

    // Average of the landmarks in a given image
    std::vector< Matrix1D<double> > barpi;

    // Number of landmarks in image i
    Matrix1D<int> ni;

    // Best exhaustive alignmnet
    Alignment* bestPreviousAlignment;

    // Show refinement
    bool showRefinement;
};

class Alignment
{
public:
    const ProgTomographAlignment *prm;
    std::vector< Matrix1D<double> > di;
    std::vector< Matrix1D<double> > rj;
    Matrix1D<double> psi;
    double rot;
    double tilt;
    Matrix1D<double> raxis;

public:
    Alignment(const ProgTomographAlignment *_prm)
    {
        prm=_prm;
        Nimg=MAT_XSIZE(_prm->allLandmarksX);
        Nlandmark=MAT_YSIZE(_prm->allLandmarksX);
        clear();
    }

    /** Clear */
    void clear()
    {
        psi.initZeros(Nimg);
        rot=90;
        tilt=90;
        raxis.initZeros(3);

        Ai.clear();
        Ait.clear();
        Aip.clear();
        Aipt.clear();
        di.clear();
        diaxis.clear();
        B1i.clear();
        B2i.clear();
        barri.clear();
        Matrix2D<double> dummy23(2,3);
        Matrix2D<double> dummy32(3,2);
        Matrix2D<double> dummy33(3,3);
        Matrix2D<double> dummy22(2,2);
        Matrix1D<double> dummy2(2);
        Matrix1D<double> dummy3(3);
        for (int i=0; i<Nimg; i++)
        {
            Ai.push_back(dummy23);
            Ait.push_back(dummy32);
            Aip.push_back(dummy23);
            Aipt.push_back(dummy32);
            di.push_back(dummy2);
            diaxis.push_back(dummy2);
            B1i.push_back(dummy33);
            B2i.push_back(dummy33);
            barri.push_back(dummy3);
        }
        allLandmarksPredictedX.initZeros(Nlandmark,Nimg);
        allLandmarksPredictedY.initZeros(Nlandmark,Nimg);
        errorLandmark.initZeros(Nlandmark);

        rj.clear();
        for (int j=0; j<Nlandmark; j++)
            rj.push_back(dummy3);
    }

    /** Assignment operator */
    Alignment & operator= (const Alignment &op)
    {
        if (this!=&op)
        {
            prm=op.prm;
            Nimg=op.Nimg;
            Nlandmark=op.Nlandmark;
            psi=op.psi;
            rot=op.rot;
            tilt=op.tilt;
            Ai=op.Ai;
            Ait=op.Ait;
            Aip=op.Aip;
            Aipt=op.Aipt;
            diaxis=op.diaxis;
            B1i=op.B1i;
            B2i=op.B2i;
            raxis=op.raxis;
            di=op.di;
            rj=op.rj;
            barri=op.barri;
            allLandmarksPredictedX=op.allLandmarksPredictedX;
            allLandmarksPredictedY=op.allLandmarksPredictedY;
            errorLandmark=op.errorLandmark;
        }
        return *this;
    }

    /** Optimize for rot */
    double optimizeGivenAxisDirection();

    /** Compute Aip, and its transpose */
    void computeGeometryDependentOfAxis();

    /** Compute Ai, and its transposes */
    void computeGeometryDependentOfRotation();

    /** Compute error */
    double computeError() const;

    /** Compute error for landmarks */
    void computeErrorForLandmarks();

    /** Update 3D model */
    void updateModel();

    /** Print an alignment */
    friend std::ostream& operator<<(std::ostream &out,
                                    Alignment &alignment);

public:
    // Number of images
    int Nimg;

    // Number of landmarks
    int Nlandmark;

    // Set of Ai matrices
    std::vector< Matrix2D<double> > Ai;

    // Set of Ait matrices
    std::vector< Matrix2D<double> > Ait;

    // Set of Ai prime matrices
    std::vector< Matrix2D<double> > Aip;

    // Set of Ai prime transpose matrices
    std::vector< Matrix2D<double> > Aipt;

    // Set of bar ri
    std::vector< Matrix1D<double> > barri;

    // Set of shifts due to raxis
    std::vector< Matrix1D<double> > diaxis;

    // Set of B1i matrices
    std::vector< Matrix2D<double> > B1i;

    // Set of B2i matrices
    std::vector< Matrix2D<double> > B2i;

    // Matrix used for updating raxis
    Matrix2D<double> Binvraxis;

    // Set of predicted landmarks (component X)
    Matrix2D<double> allLandmarksPredictedX;

    // Set of predicted landmarks (component Y)
    Matrix2D<double> allLandmarksPredictedY;

    // Set of errors associated to each landmark
    MultidimArray<double> errorLandmark;
};

/** Compute the optimal affine transformation between two images.
    The maximum shift is limited. 
    A12 is the homogeneous matrix that transforms coordinates of 1
    into coordinates of 2, A21 is just the opposite.*/
double computeAffineTransformation(const MultidimArray<double> &I1,
                                   const MultidimArray<double> &I2, int maxShift, int maxIterDE,
                                   Matrix2D<double> &A12, Matrix2D<double> &A21, bool show,
                                   double thresholdAffine, bool globalAffine, bool isMirror,
                                   int pyramidLevel);
//@}
#endif
