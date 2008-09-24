/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2008)
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
#ifndef _PROG_ANGULAR_ASSIGN_FOR_TILTSERIES
#define _PROG_ANGULAR_ASSIGN_FOR_TILTSERIES

#include <vector>
#include <data/matrix1d.h>
#include <data/selfile.h>

/**@defgroup AngularAssignTiltSeries angular_assign_for_tilt_series
             (Simultaneous assignment)
   @ingroup ReconsLibraryPrograms */
//@{
/** Landmark class.
    A landmark is a position (x,y) and the index of the image
    from which this landmark has been taken.*/
class Landmark {
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
class Prog_tomograph_alignment {
public:
    /// SelFile with all images
    FileName fnSel;
   
    /// SelFile with all images at the original scale
    FileName fnSelOrig;
   
    /// Output root
    FileName fnRoot;
   
    /// Look for local affine transformation
    bool localAffine;
   
    /// Sequence Length
    int seqLength;

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

    /// Correlation threshold for a good landmark
    double corrThreshold;
    
    /// Show affine transformations
    bool showAffine;
    
    /// Threshold affine
    double thresholdAffine;

    /// Read parameters from argument line
    void read(int argc, char **argv);

    /// Show parameters
    void show();
   
    /// Usage
    void usage() const;

    /// Produce side info
    void produceSideInfo();

    /// Produce information from landmarks
    void produceInformationFromLandmarks();
   
    /** Refine landmark.
        ii is the index of the original image, jj is the index in the
        image at which the landmark is being refined. rii and rjj are
        the corresponding landmark positions in both images.
        
        The function returns whether the landmark is accepted or not. */
    bool refineLandmark(int ii, int jj, const Matrix1D<double> &rii,
        Matrix1D<double> &rjj) const;

    /** Refine chain. */
    bool refineChain(LandmarkChain &chain);

    /// Generate landmark set
    void generateLandmarkSet();

    /// Write landmark set
    void writeLandmarkSet(const FileName &fnLandmark) const;

    /// Read landmark set
    void readLandmarkSet(const FileName &fnLandmark);

    /// Remove Outliers
    void removeOutlierLandmarks(const Alignment &alignment);

    /// Align images
    void alignImages(const Alignment &alignment);

    /// Run: do the real work
    void run();
   
public:
    // Selfile with the input images
    SelFile SF;

    // Selfile with the input images
    SelFile SForig;

    // List of image pointers
    std::vector < Matrix2D<double> *> img;

    // List of tilt angles
    std::vector < double > tiltList;

    // List of names
    std::vector < std::string > name_list;

    // Total number of images
    int Nimg;
   
    // Element [i][j] is the transformation of coordinates of i into
    // coordinates of j
    std::vector< std::vector< Matrix2D<double> > > affineTransformations;

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

class Alignment {
public:
    const Prog_tomograph_alignment *prm;
    std::vector< Matrix1D<double> > di;
    std::vector< Matrix1D<double> > rj;
    Matrix1D<double> psi;
    double rot;
    double tilt;
    bool optimizePsi;

public:
    Alignment(const Prog_tomograph_alignment *_prm)
    {
        prm=_prm;
        Nimg=XSIZE(_prm->allLandmarksX);
        Nlandmark=YSIZE(_prm->allLandmarksX);
        clear();
    }

    /** Clear */
    void clear()
    {
        psi.initZeros(Nimg);
        rot=90;
        tilt=90;
        optimizePsi=false;
        
        Ai.clear();
        Ait.clear();
        Aip.clear();
        Aipt.clear();
        di.clear();
        barri.clear();
        Matrix2D<double> dummy23(2,3);
        Matrix2D<double> dummy32(3,2);
        Matrix1D<double> dummy2(2);
        Matrix1D<double> dummy3(3);
        for (int i=0; i<Nimg; i++)
        {
            Ai.push_back(dummy23);
            Ait.push_back(dummy32);
            Aip.push_back(dummy23);
            Aipt.push_back(dummy32);
            di.push_back(dummy2);
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
            optimizePsi=op.optimizePsi;
            Ai=op.Ai;
            Ait=op.Ait;
            Aip=op.Aip;
            Aipt=op.Aipt;
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

    // Set of predicted landmarks (component X)
    Matrix2D<double> allLandmarksPredictedX;

    // Set of predicted landmarks (component Y)
    Matrix2D<double> allLandmarksPredictedY;
    
    // Set of errors associated to each landmark
    Matrix1D<double> errorLandmark;
};

/** Compute the optimal affine transformation between two images.
    The maximum shift is limited. 
    A12 is the homogeneous matrix that transforms coordinates of 1
    into coordinates of 2, A21 is just the opposite.*/
void computeAffineTransformation(const Matrix2D<double> &I1,
    const Matrix2D<double> &I2, int maxShift, int maxIterDE,
    const FileName &fn_affine,
    Matrix2D<double> &A12, Matrix2D<double> &A21);
//@}
#endif
