/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#ifndef METADATALABEL_H
#define METADATALABEL_H

#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "xmipp_funcs.h"
#include "xmipp_strings.h"

class MDLabelData;
class MDObject;
class MDLabelStaticInit;

/** @addtogroup MetaData
 * @{
 */

/** Enumerate all posible labels to use in MetaData.
 */
enum MDLabel
{
    MDL_UNDEFINED = -1,
    MDL_FIRST_LABEL, ///< The label MDL_OBJID is special and should not be used
    MDL_OBJID = MDL_FIRST_LABEL, ///< object id (int), NOTE: This label is special and shouldn't be used

    MDL_ANGLE_PSI, ///< Psi angle of an image (double,degrees)
    MDL_ANGLE_PSI2, ///< Psi angle of an image (double,degrees)
    MDL_ANGLE_PSI_DIFF, ///< difference between psi angles (double,degrees)
    MDL_ANGLE_ROT, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLE_ROT2, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLE_ROT_DIFF, ///< difference between rot angles (double,degrees)
    MDL_ANGLE_TILT, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLE_TILT2, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLE_TILT_DIFF, ///< difference between tilt angles (double,degrees)
    MDL_ANGLE_DIFF, ///< difference between two angles (double,degrees)
    MDL_ANGLE_Y,   ///< Angle between y-axis and tilt-axis (double, degrees) for untilted micrographs
    MDL_ANGLE_Y2,   ///< Angle between y-axis and tilt-axis (double, degrees) for tilted micrographs
    MDL_AVG, ///< average value (double)
    MDL_AVG_CHANGES_ORIENTATIONS, /// Average change in angular orientation (double degrees)
    MDL_AVG_CHANGES_OFFSETS, /// Average change in offset (double pixels)
    MDL_AVG_CHANGES_CLASSES, /// Average change in class assignment(double dimensionaless)
    MDL_BGMEAN, ///< Mean background value for an image
    MDL_BLOCK_NUMBER, ///< Current block number (for incremental EM)

    MDL_CL2D_CHANGES, ///< Number of changes between iterations
    MDL_CL2D_SIMILARITY, ///< Average cross-correlation for the image (double)
    MDL_CLASS_COUNT, ///< Number of images assigned to the same class as this image
    MDL_CLASSIFICATION_DATA, ///< Data vector for classification (vector double)
    MDL_CLASSIFICATION_DATA_SIZE, ///< Size of data vectors for classification (int)
    MDL_CLASSIFICATION_DPR_05, ///< Differential Phase Residual evaluated at FRC=0.5
    MDL_CLASSIFICATION_INTRACLASS_DISTANCE, ///< Average intraclass distance (double)
    MDL_CLASSIFICATION_FRC_05, ///< Digital frequency at which the FRC drops below 0.5 (double)
    MDL_COMMENT, ///< A comment for this object /*** NOTE THIS IS A SPECIAL CASE AND SO IS TREATED ***/
    MDL_COST, ///< Cost for the image (double)
    MDL_COUNT, ///< Number of elements of a type (int) [this is a genereic type do not use to transfer information to another program]
    MDL_COUNT2, ///< Number of elements of a type (int) [this is a genereic type do not use to transfer information to another program]

    MDL_CRYSTAL_CELLX, ///< Cell location for crystals
    MDL_CRYSTAL_CELLY, ///< Cell location for crystals
    MDL_CRYSTAL_LATTICE_A,   /// < Lattice vector for projection (vector double)
    MDL_CRYSTAL_LATTICE_B,   /// < Lattice vector for projection (vector double)
    MDL_CRYSTAL_DISAPPEAR_THRE,   /// < Disappearing threshold (double)
    MDL_CRYSTAL_SHFILE,   /// < Shift file for crystal projection
    MDL_CRYSTAL_ORTHO_PRJ,   /// <Orthogonal projection or not (bool)
    MDL_CRYSTAL_PROJ,   /// < Have a crystal projection (bool)
    MDL_CRYSTAL_SHIFTX, ///< Shift for the image in the X axis (double) for crystals
    MDL_CRYSTAL_SHIFTY, ///< Shift for the image in the Y axis (double) for crystals
    MDL_CRYSTAL_SHIFTZ, ///< Shift for the image in the Z axis (double) for crystals

    MDL_CTF_INPUTPARAMS, ///< Parameters file for the CTF Model (std::string)
    MDL_CTF_MODEL, ///< Name for the CTF Model (std::string)
    MDL_CTF_MODEL2, ///< Name for another CTF model (std::string)
    MDL_CTF_SAMPLING_RATE, ///< Sampling rate
    MDL_CTF_SAMPLING_RATE_Z, ///< Sampling rate in Z direction
    MDL_CTF_VOLTAGE, ///< Microscope voltage (kV)
    MDL_CTF_DEFOCUSA, ///< average defocus (Angtroms)
    MDL_CTF_DEFOCUSU, ///< Defocus U (Angstroms)
    MDL_CTF_DEFOCUSV, ///< Defocus V (Angstroms)
    MDL_CTF_X0, ///< The CTF is valid within (x0,y0) to (xF,yF) in the micrograph coordinates
    MDL_CTF_Y0, ///< The CTF is valid within (x0,y0) to (xF,yF) in the micrograph coordinates
    MDL_CTF_XF, ///< The CTF is valid within (x0,y0) to (xF,yF) in the micrograph coordinates
    MDL_CTF_YF, ///< The CTF is valid within (x0,y0) to (xF,yF) in the micrograph coordinates
    MDL_CTF_DEFOCUS_PLANEUA, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEUB, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEUC, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEVA, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEVB, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_PLANEVC, ///< Defocus = A*x+B*y+C
    MDL_CTF_DEFOCUS_ANGLE, ///< Defocus angle (degrees)
    MDL_CTF_CS, ///< Spherical aberration
    MDL_CTF_CA, ///< Chromatic aberration
    MDL_CTF_GROUP, ///< group images by defocus
    MDL_CTF_ENERGY_LOSS, ///< Energy loss
    MDL_CTF_LENS_STABILITY, ///< Lens stability
    MDL_CTF_CONVERGENCE_CONE, ///< Convergence cone
    MDL_CTF_LONGITUDINAL_DISPLACEMENT, ///< Longitudinal displacement
    MDL_CTF_TRANSVERSAL_DISPLACEMENT, ///< Transversal displacemente
    MDL_CTF_Q0, ///< Inelastic absorption
    MDL_CTF_K, ///< CTF gain
    MDL_CTF_BG_GAUSSIAN_K, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_SIGMAU, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_SIGMAV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_CU, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_CV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_ANGLE, ///< CTF Background parameter
    MDL_CTF_BG_SQRT_K, ///< CTF Background parameter
    MDL_CTF_BG_SQRT_U, ///< CTF Background parameter
    MDL_CTF_BG_SQRT_V, ///< CTF Background parameter
    MDL_CTF_BG_SQRT_ANGLE, ///< CTF Background parameter
    MDL_CTF_BG_BASELINE, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_K, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_SIGMAU, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_SIGMAV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_CU, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_CV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN2_ANGLE, ///< CTF Background parameter
    MDL_CTF_CRIT_PSDCORRELATION90, ///< PSD correlation at 90 degrees
    MDL_CTF_CRIT_FIRSTZERORATIO, ///< First zero ratio
    MDL_CTF_CRIT_FIRSTZEROAVG, ///< First zero average (in Angstroms)
    MDL_CTF_CRIT_FIRSTZERODISAGREEMENT, ///< First zero disagreement with second model (in Angstroms)
    MDL_CTF_CRIT_MAXFREQ, ///< Maximum frequency at which the envelope drops below 0.1 (in Angstroms)
    MDL_CTF_CRIT_DAMPING, ///< Minimum damping at border
    MDL_CTF_CRIT_PSDRADIALINTEGRAL, ///< Integral of the radial PSD
    MDL_CTF_CRIT_FITTINGSCORE, ///< Score of the fitting
    MDL_CTF_CRIT_FITTINGCORR13, ///< Correlation between the 1st and 3rd ring of the CTF
    MDL_CTF_CRIT_PSDVARIANCE, ///< PSD variance
    MDL_CTF_CRIT_PSDPCA1VARIANCE, ///< Variance in the first principal component of the PSDs
    MDL_CTF_CRIT_PSDPCARUNSTEST, ///< Runs test on the projection of the PSD on the first principal component
    MDL_CTF_CRIT_NORMALITY, ///< Normality test between histogram of micrography and gaussian distribution
    MDL_CTF_DOWNSAMPLE_PERFORMED, ///< Downsampling performed to estimate the CTF
    MDL_CTF_DIMENSIONS, // Size in pixels of the 3D PSF to be created (Xdim, Ydim, Zdim)
    MDL_CTF_LAMBDA, /// Wavelength (nm)
    MDL_CTF_XRAY_LENS_TYPE, ///Algorithm used to generate Xray PSF
    MDL_CTF_XRAY_OUTER_ZONE_WIDTH, /// Outermost zone width of the X-ray Fresnel lens (nm)
    MDL_CTF_XRAY_ZONES_NUMBER, // Number of zones of the X-ray Fresnel lens
    MDL_DATATYPE, ///< if read from file original image datatype, this is an struct defined in image
    MDL_DEFGROUP, ///< Defocus group
    MDL_DIMRED, ///< Projection onto a reduced manifold (vector double)
    MDL_DIRECTION, ///< Direction in 3D

    MDL_DM3_IDTAG,
    MDL_DM3_NODEID,
    MDL_DM3_NUMBER_TYPE,
    MDL_DM3_PARENTID,
    MDL_DM3_TAGCLASS,
    MDL_DM3_TAGNAME,
    MDL_DM3_SIZE,
    MDL_DM3_VALUE,

    //End of labels

    MDL_ENABLED, ///< Is this image enabled? (int [-1 or 1])

    MDL_DATE,// < timestamp (string)
    MDL_TIME,// <  time in seconds (double)

    MDL_FLIP, ///< Flip the image? (bool)
    MDL_FOM, ///< Figure of Merit in 0-1 range (double)
    MDL_IDX, ///< Index within a list (size_t)
    MDL_IMAGE, ///< Name of an image (std::string)
    MDL_IMAGE_ORIGINAL, ///< Name of an image from which MDL_IMAGE is coming from
    MDL_IMAGE_REF, ///< Name of of the class image from which MDL_IMAGE is coming from
    MDL_IMAGE_TILTED, ///< Name of the tilted images associated to MDL_IMAGE
    MDL_IMGMD, ///< Name of Metadata file for all images (string)
    MDL_IMAGE1, ///< Image associated to this object (std::string)
    MDL_IMAGE2, ///< Image associated to this object (std::string)
    MDL_IMAGE3, ///< Image associated to this object (std::string)
    MDL_IMAGE4, ///< Image associated to this object (std::string)
    MDL_IMAGE5, ///< Image associated to this object (std::string)
    MDL_INTSCALE, ///< Intensity scale for an image
    MDL_ITEM_ID, ///< Unique identifier for items inside a list or set (std::size_t)
    MDL_ITER, ///< Current iteration number (int)
    MDL_KERDENSOM_FUNCTIONAL, ///< Functional value (double)
    MDL_KERDENSOM_REGULARIZATION, ///< Regularization value (double)
    MDL_KERDENSOM_SIGMA, ///< Sigma value (double)
    MDL_KEYWORDS, ///< Keywords associated to this line, should be a single string block (do not use spaces as separators)
    MDL_KSTEST, ///<KS-test statistics
    MDL_LL, ///< contribution of an image to log-likelihood value
    MDL_MAGNIFICATION, /// Magnification of microscope
    MDL_MAPTOPOLOGY, ///< Map topology (KerDenSOM, ...)
    MDL_MASK, ///< Name of a mask associated to image
    MDL_MAXCC, ///< Maximum cross-correlation for the image (double)
    MDL_MAX, ///< Maximum value (double)
    MDL_MICROGRAPH, ///< Name of a micrograph (std::string)
    MDL_MICROGRAPH_PARTICLES, ///< Name of a position file (std::string)
    MDL_MICROGRAPH_ORIGINAL, ///< Name of the original micrograph, MDL_MICROGRAPH is probably a downsampled version of this one (std::string)
    MDL_MICROGRAPH_TILTED, ///< Name of the corresponding tilted micrograph (std::string)
    MDL_MICROGRAPH_TILTED_ORIGINAL, ///< Name of the corresponding original, tilted micrograph (std::string)
    MDL_MIN, ///< Minimum value (double)
    MDL_MIRRORFRAC, ///< Mirror fraction for a Maximum Likelihood model
    MDL_MISSINGREGION_NR, ///< Number of missing region in subtomogram
    MDL_MISSINGREGION_TYPE, ///< Type of missing region in subtomogram
    MDL_MISSINGREGION_THY0, ///< Initial tilt angle in Y for missing region in subtomogram
    MDL_MISSINGREGION_THYF, ///< Final tilt angle in Y for missing region in subtomogram
    MDL_MISSINGREGION_THX0, ///< Initial tilt angle in X for missing region in subtomogram
    MDL_MISSINGREGION_THXF, ///< Final tilt angle in X for missing region in subtomogram

    MDL_MLF_CTF,    ///< MLF CTF estimated value (double)
    MDL_MLF_WIENER, ///< MLF Wiener filter (double)
    MDL_MLF_SIGNAL, ///< MLF signal (double)
    MDL_MLF_NOISE,  ///< MLF Wiener filter (double)

    MDL_MODELFRAC, ///< Model fraction (alpha_k) for a Maximum Likelihood model
    MDL_NEIGHBORS, ///< Vector of indexes to points some "neighbors"
    MDL_NEIGHBOR, ///< Particular neighbor (pointed myNEIGHBORS)
    MDL_NEIGHBORHOOD_RADIUS, ///< Radius of the neigborhood (radians)
    MDL_NMA, ///< Normal mode displacements (vector double)
    MDL_NMA_MODEFILE, ///< File with an NMA mode
    MDL_NMA_COLLECTIVITY, ///< NMA Collectivity of a given mode
    MDL_NMA_SCORE, ///< NMA Score of a given mode
    MDL_NOISE_ANGLES, ///< Noise description for projected angles
    MDL_NOISE_PARTICLE_COORD, ///< Noise description for particle's center coordenates (when projecting)
    MDL_NOISE_COORD,  //Use instead of MDL_NOISE_PARTICLE_COORD in future
    MDL_NOISE_PIXEL_LEVEL, ///< Noise description for pixels' gray level (when projecting)
    MDL_ORDER, /// auxiliary label to be used as an index (long)
    MDL_ORIGIN_X, ///< Origin for the image in the X axis (double)
    MDL_ORIGIN_Y, ///< Origin for the image in the Y axis (double)
    MDL_ORIGIN_Z, ///< Origin for the image in the Z axis (double)

    MDL_PHANTOM_BGDENSITY, ///< Phantom background density (double)
    MDL_PHANTOM_FEATURE_CENTER, ///< Center of the feature (vector double)
    MDL_PHANTOM_FEATURE_DENSITY, ///< The density of the feature (double)
    MDL_PHANTOM_FEATURE_OPERATION, ///< Operation in case of overlapping features (+,-)
    MDL_PHANTOM_FEATURE_SPECIFIC, ///< Specific parameters for a feature (vector double)
    MDL_PHANTOM_FEATURE_TYPE, ///< Type of the feature (Sphere, Blob, ...) (std::string)
    MDL_PHANTOM_SCALE, ///< Number which will multiply all features (double)

    MDL_MACRO_CMD, //ImageJ macro command on picker
    MDL_MACRO_CMD_ARGS, //ImageJ macro args on picker
    MDL_COLOR, ///< Color for particle picking
    MDL_PICKING_TEMPLATES, ///< Number of templates
    MDL_PICKING_STATE, ///< State for particle picking
    MDL_PICKING_MICROGRAPH_STATE, ///< Micrograph state for particle picking
    MDL_PICKING_AUTOPICKPERCENT,
    MDL_PICKING_PARTICLE_SIZE, ///< Particle size for particle picking
    MDL_PMAX, ///< Maximum value of normalized probability function (now called "Pmax/sumP") (double)
    MDL_AVGPMAX, ///< Average (per class) of the maximum value of normalized probability function) (double)
    MDL_POINTSASYMETRICUNIT, /// < Number of non-redundant projections directions (size_t)

    MDL_PRJ_DIMENSIONS, // X,Y dimensions for the generated projections
    MDL_PRJ_ANGFILE,  ///< File for generated angles
    MDL_PRJ_PSI_NOISE,  /// < Psi angle dev and mean noise (vector double)
    MDL_PRJ_PSI_RANDSTR, /// < Type of randomness for Psi (std::string)
    MDL_PRJ_PSI_RANGE,  /// < Psi angle range (vector double)
    MDL_PRJ_ROT_NOISE ,  /// < Rotational angle dev and mean noise (vector double)
    MDL_PRJ_ROT_RANDSTR,  /// < Type of randomness for Rotational (std::string)
    MDL_PRJ_ROT_RANGE,
    MDL_PRJ_TILT_NOISE,  /// < Tilt angle dev and mean noise (vector double)
    MDL_PRJ_TILT_RANDSTR,  /// < Type of randomness for Tilt (std::string)
    MDL_PRJ_TILT_RANGE, // Vector with the initial and final tilt angle values, and step size
    MDL_PRJ_VOL,        // Volume file name to generate projections from

    MDL_PROGRAM,// <  program name
    MDL_USER,// <  user name

    MDL_DIMENSIONS_3D,  // X,Y,Z dimensions
    MDL_DIMENSIONS_2D,  // X,Y dimensions
    MDL_PSD, ///< A Power Spectrum Density file name (std::string)
    MDL_PSD_ENHANCED, ///< A enhanced Power Spectrum Density file name (std::string)
    MDL_RANDOMSEED, ///< Seed for random number generator
    MDL_REF3D, ///< 3D Class to which the image belongs (int)
    MDL_REF, ///< Class to which the image belongs (int)
    MDL_REF2, ///< Store a second class (int)
    MDL_REFMD, ///< Name of Metadata file for all references(string)

    MDL_RESOLUTION_DPR, ///<differential phase residual (double)
    MDL_RESOLUTION_ERRORL2, ///<Error in l2 (double)
    MDL_RESOLUTION_FRC, ///<Fourier shell correlation (double)
    MDL_RESOLUTION_FRCRANDOMNOISE, ///<Fourier shell correlation noise (double)
    MDL_RESOLUTION_FREQ, ///<Frequency in 1/A (double)
    MDL_RESOLUTION_FREQREAL, ///< Frequency in A (double)
    MDL_RESOLUTION_SSNR, ///<Fourier shell correlation (double)

    MDL_SAMPLINGRATE, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_ORIGINAL, ///< original sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_X, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_Y, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_Z, ///< sampling rate in A/pixel (double)

    MDL_SCALE, ///< scaling factor for an image or volume (double)
    MDL_SELFILE, ///< Name of an image (std::string)
    MDL_SERIE, ///< A collection of micrographs, e.g. a tilt serie (std::string)
    MDL_SHIFT_X, ///< Shift for the image in the X axis (double)
    MDL_SHIFT_X2, ///< Shift for the image in the X axis (double)
    MDL_SHIFT_X_DIFF, ///< difference in Shift along X axis (double)
    MDL_SHIFT_Y, ///< Shift for the image in the Y axis (double)
    MDL_SHIFT_Y2, ///< Shift for the image in the Y axis (double)
    MDL_SHIFT_Y_DIFF, ///< difference in Shift along  Y axis (double)
    MDL_SHIFT_Z, ///< Shift for the image in the Z axis (double)
    MDL_SHIFT_DIFF, ///< shift difference (double)
    MDL_SIGMANOISE, ///< Standard deviation of the noise in ML model
    MDL_SIGMAOFFSET, ///< Standard deviation of the offsets in ML model
    MDL_SIGNALCHANGE, ///< Signal change for an image
    MDL_STDDEV, ///<stdandard deviation value (double)
    MDL_SUM, ///< Sum of elements of a given type (double) [this is a genereic type do not use to transfer information to another program]
    MDL_SUMWEIGHT, ///< Sum of all weights in ML model
    MDL_SYMNO, ///< Symmetry number for a projection (used in ART)
    MDL_TRANSFORMATIONMTRIX, ///< transformation matrix(vector double)

    MDL_TEST_SIZE,// < number of test assigned to a program

    MDL_VOLUME_SCORE1,/// < Score 1 for volumes
    MDL_VOLUME_SCORE2,/// < Score 2 for volumes
    MDL_VOLUME_SCORE3,/// < Score 3 for volumes
    MDL_VOLUME_SCORE4,/// < Score 4 for volumes
    MDL_VOLTAGE, ///< microscope voltage (double)
    MDL_WEIGHT, ///< Weight assigned to the image (double)
    MDL_WROBUST, ///< Weight of t-student distribution in robust Maximum likelihood
    MDL_X, ///< X component (double)
    MDL_XCOOR, ///< X component (int)
    MDL_XCOOR_TILT, ///< X component in tilted micrograph (int)
    MDL_XSIZE, ///< X size (int)
    MDL_Y, ///< Y component (double)
    MDL_YCOOR, ///< Y component (int)
    MDL_YCOOR_TILT, ///< Y component in tilted micrograph (int)
    MDL_YSIZE, ///< Y size (int)
    MDL_Z, ///< Z component (double)
    MDL_ZCOOR, ///< Z component (int)
    MDL_ZSCORE, ///< Global Z Score (double)
    MDL_ZSCORE_SHAPE1, ///< Z Score (double)
    MDL_ZSCORE_SHAPE2, ///< Z Score (double)
    MDL_ZSCORE_SNR1, ///< Z Score (double)
    MDL_ZSCORE_SNR2, ///< Z Score (double)
    MDL_ZSCORE_HISTOGRAM, ///< Z Score (double)
    MDL_ZSIZE, ///< Z size (int)

    MDL_LAST_LABEL  // **** NOTE ****: Do keep this label always at the end,it is here for looping purposes
};//close enum Label

/** Macro for iterate over all labels */
#define FOR_ALL_LABELS() for (int _label = MDL_FIRST_LABEL; _label < MDL_LAST_LABEL; ++_label)

/** Possible types of the values of labels */
enum MDLabelType
{
    LABEL_NOTYPE = -1,
    LABEL_INT,
    LABEL_BOOL,
    LABEL_DOUBLE,
    LABEL_STRING,
    LABEL_VECTOR_DOUBLE,
    LABEL_SIZET,
    LABEL_VECTOR_SIZET
};

/** Possible types of the values of labels */
enum MDLabelTag
{
    TAGLABEL_NOTAG = 0,
    TAGLABEL_TEXTFILE=0x1,
    TAGLABEL_METADATA=0x3,
    TAGLABEL_CTFPARAM=0x5,
    TAGLABEL_IMAGE=0x8,
    TAGLABEL_VOLUME=0x10,
    TAGLABEL_STACK=0x20,
    TAGLABEL_MICROGRAPH=0x48,
    TAGLABEL_PSD=0x88
};

/**Just an utility function */
bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label);

//Just an struct to store type and string alias
class MDLabelData
{
public:
    MDLabelType type;
    String str;
    int tags;
    //Default constructor
    MDLabelData()
    {
        type = LABEL_NOTYPE;
        tags=TAGLABEL_NOTAG;
    }

    MDLabelData(MDLabelType t, const String &s, int tags)
    {
        type = t;
        str = s;
        this->tags=tags;
    }
}
;//close class MDLabelData

/** Union to store values */
typedef union
{
    bool boolValue;
    int intValue;
    size_t longintValue;
    double doubleValue;
    String * stringValue;
    std::vector<double> * vectorValue;
    std::vector<size_t> * vectorValueLong;
} ObjectData;

#define _SPACE ' '
#define _QUOT '\''
#define _DQUOT '"'

/** Class to hold the labels values and type
 *
 */
class MDObject
{
public:

    ObjectData data;
    char chr; //literal char for string, could be SPACE, QUOT or DQUOT

    void labelTypeCheck(MDLabelType checkingType) const;
    void copy(const MDObject &obj);

    MDLabel label;
    MDLabelType type;
    /** Copy constructor */
    MDObject(const MDObject & obj);
    /** Assign operator */
    MDObject & operator = (const MDObject &obj);
    //Just a simple constructor with the label
    //dont do any type checking as have not value yet
    MDObject(MDLabel label);
    ///Constructors for each Label supported type
    ///these constructor will do the labels type checking
    MDObject(MDLabel label, const int &intValue);
    MDObject(MDLabel label, const double &doubleValue);
    MDObject(MDLabel label, const bool &boolValue);
    MDObject(MDLabel label, const String &stringValue);
    MDObject(MDLabel label, const std::vector<double> &vectorValue);
    MDObject(MDLabel label, const std::vector<size_t> &vectorValueLong);
    MDObject(MDLabel label, const size_t &longintValue);
    MDObject(MDLabel label, const float &floatValue);
    MDObject(MDLabel label, const char * &charValue);

    /// Destructor
    ~MDObject();

    //These getValue also do a compilation type checking
    //when expanding templates functions and only
    //will allow the supported types
    //TODO: think if the type check if needed here
    void  getValue(int &iv) const;
    void  getValue(double &dv) const;
    void  getValue(bool &bv) const;
    void  getValue(String &sv) const;
    void  getValue(std::vector<double> &vv) const;
    void  getValue(std::vector<size_t> &vv) const;
    void  getValue(size_t &lv) const;
    void  getValue(float &floatvalue) const;
    void  getValue(char*  &charvalue) const;

    void  setValue(const int &iv);
    void  setValue(const double &dv);
    void  setValue(const bool &bv);
    void  setValue(const String &sv);
    void  setValue(const std::vector<double> &vv);
    void  setValue(const std::vector<size_t> &vv);
    void  setValue(const size_t &lv);
    void  setValue(const float &floatvalue);
    void  setValue(const char*  &charvalue);

#define DOUBLE2STREAM(d) \
        if (withFormat) {\
                (os) << std::setw(12); \
                (os) << (((d) != 0. && ABS(d) < 0.001) ? std::scientific : std::fixed);\
            } os << d;

#define INT2STREAM(i) \
        if (withFormat) os << std::setw(10); \
        os << i;

    void toStream(std::ostream &os, bool withFormat = false, bool isSql=false, bool escape=true) const;
    String toString(bool withFormat = false, bool isSql=false) const;
    bool fromStream(std::istream &is, bool fromString=false);
    friend std::istream& operator>> (std::istream& is, MDObject &value);
    friend std::ostream& operator<< (std::ostream& is, MDObject &value);
    bool fromString(const String &str);
    bool fromChar(const char * str);

    friend class MDSql;
}
; //close class MDValue

/** Explicit instantiation */
#ifndef __APPLE__
template class std::vector<MDObject *>
;
#endif

/** Class for holding an entire row of posible MDObject */
class MDRow
{
public:
    //Reserve space for the maximum different labels
    //this will allows constant access to each object indexing by labels
    MDObject * objects[MDL_LAST_LABEL];
    MDLabel order[MDL_LAST_LABEL];
    int _size; //Number of active labels

public:
    /** Empty constructor */
    MDRow();
    /** Copy constructor */
    MDRow(const MDRow & row);

    /** Assignment */
    MDRow& operator = (const MDRow &row);

    /** Destructor */
    ~MDRow();
    /** True if this row contains this label */
    bool containsLabel(MDLabel label) const;

    /** Add a new label */
    void addLabel(MDLabel label);

    /** Clear elements of the row */
    void clear();
    /** Return number of labels present */
    int size() const;
    /** Function to test whether is empty */
    bool empty() const;

    /** Reset the values of the labels related to
     *  geometry to their default values
     */
    void resetGeo(bool addLabels = true);

    /** Get object */
    MDObject * getObject(MDLabel label) const;

    /** Get value */
    template <typename T>
    bool getValue(MDLabel label, T &d) const
    {
        if (objects[label] == NULL)
            return false;

        objects[label]->getValue(d);
        return true;
    }
    /** Get value */
    //    template <typename T >
    //    void getValueOrAbort(MDLabel label, T &d) const
    //    {
    //        if (!getValue(label, d))
    //         REPORT_ERROR(ERR_ARG_MISSING,(String)"Cannot find label: " + MDL::label2Str(label) );
    //        //formatString("%d",label) );
    //    }
    //weite function as macro since MDL::label2Str is not availale at class compilation time
#define rowGetValueOrAbort(__row,__label, __d)\
        if (!__row.getValue(__label, __d))\
         REPORT_ERROR(ERR_ARG_MISSING,(String)"Cannot find label: " + MDL::label2Str(__label) );

    /** Get value */
    template <typename T, typename T1>
    void getValueOrDefault(MDLabel label, T &d, T1 def) const
    {
        if (!getValue(label, d))
            d = (T) def;
    }

    bool getValue(MDObject &object) const;

    /** Set value
     *
     * @param label    Metadata label to be set
     * @param d        Value of the label
     * @param addLabel Add label if is not already contained
     */
    template <typename T>
    void setValue(MDLabel label, const T &d, bool addLabel = true)
    {
        if (objects[label] != NULL)
            objects[label]->setValue(d);
        else if (addLabel)
        {
            objects[label] = new MDObject(label, d);
            order[_size] = label;
            ++_size;
        }
    }

    void setValue(const MDObject &object);

    /** Show */
    friend std::ostream& operator << (std::ostream &out, const MDRow &row);

private:
    void copy(const MDRow &row);
};

/** Static class to group some functions with labels.
 * This class holds several function to work with labels.
 * Also performs the static initialization of the labels data.
 */
class MDL
{
public:
    /** @name String conversions
     * @{
     */
    /** Converts an string to MDLabel */
    static void str2LabelVector(const String &labelsStr, std::vector<MDLabel> &labels);
    static MDLabel str2Label(const String &labelName);
    /** Converts MDLabel to string */
    static String label2Str(const MDLabel label);
    /** Converts MDLabel to string representing SQL column*/
    static String label2SqlColumn(const MDLabel label);
    /** Return the type of the label as String */
    static String labelType2Str(MDLabelType type);
    /** @} */

    /** @name Type checks
     * Several label type checks
     * @{
     */
    static bool isInt(const MDLabel label);
    static bool isLong(const MDLabel label);
    static bool isBool(const MDLabel label);
    static bool isString(const MDLabel label);
    static bool isDouble(const MDLabel label);
    static bool isVector(const MDLabel label);
    static bool isVectorLong(const MDLabel label);
    static bool isValidLabel(const MDLabel label);
    static bool isValidLabel(const String &labelName);
    static MDLabelType labelType(const MDLabel label);
    static MDLabelType labelType(const String &labelName);
    static bool hasTag(const MDLabel label, const int tags);
    static bool isTextFile(const MDLabel label);
    static bool isMetadata(const MDLabel label);
    static bool isCtfParam(const MDLabel label);
    static bool isImage(const MDLabel label);
    static bool isVolume(const MDLabel label);
    static bool isStack(const MDLabel label);
    static bool isMicrograph(const MDLabel label);
    static bool isPSD(const MDLabel label);
    static std::map<String, MDLabel>& getLabelDict();
    static MDRow emptyHeader;
    /** @} */

private:
    //Array of MDLabelData pointers
    static MDLabelData * data[MDL_LAST_LABEL];
    static std::map<std::string, MDLabel> names;
    static MDLabelStaticInit initialization; //Just for initialization

    /** Add predefined labels to be used in metadata */
    static void addLabel(MDLabel label, MDLabelType type, const String &name, int tags=TAGLABEL_NOTAG);
    /** Add an alias for an existing label */
    static void addLabelAlias(MDLabel label, const String &alias);

    friend class MDLabelStaticInit;
public:
    /** public acces to addLabelAlias, add alias in running time,
      useful to read non xmipp star files
      */
    static void addTmpLabelAlias(MDLabel label, const String &alias);
}
;//close class MLD definition


/** Just to work as static constructor for initialize labels data.
 */
class MDLabelStaticInit
{
private:
    MDLabelStaticInit()
    {
        //Just to be safe, initialize all data to null
        for (int i = (int)MDL_FIRST_LABEL; i <= (int)MDL_LAST_LABEL; ++i)
            MDL::data[i] = NULL;

        ///==== Add labels entries from here in the SAME ORDER as declared in ENUM ==========
        //The label MDL_OBJID is special and should not be used
        MDL::addLabel(MDL_OBJID, LABEL_SIZET, "objId");

        //MDL::addLabel(MDL_ANGLE_COMPARISON, LABEL_VECTOR_DOUBLE, "angle_comparison");
        //MDL::addLabelAlias(MDL_ANGLE_COMPARISON, "angleComparison"); //3.0

        MDL::addLabel(MDL_ANGLE_PSI, LABEL_DOUBLE, "anglePsi");
        MDL::addLabelAlias(MDL_ANGLE_PSI, "psi");
        MDL::addLabel(MDL_ANGLE_PSI2, LABEL_DOUBLE, "anglePsi2");
        MDL::addLabelAlias(MDL_ANGLE_PSI2, "psi2");
        MDL::addLabel(MDL_ANGLE_PSI_DIFF, LABEL_DOUBLE, "anglePsiDiff");
        MDL::addLabel(MDL_ANGLE_ROT, LABEL_DOUBLE, "angleRot");
        MDL::addLabelAlias(MDL_ANGLE_ROT, "rot");
        MDL::addLabel(MDL_ANGLE_ROT2, LABEL_DOUBLE, "angleRot2");
        MDL::addLabelAlias(MDL_ANGLE_ROT2, "rot2");
        MDL::addLabel(MDL_ANGLE_ROT_DIFF, LABEL_DOUBLE, "angleRotDiff");
        MDL::addLabel(MDL_ANGLE_TILT, LABEL_DOUBLE, "angleTilt");
        MDL::addLabelAlias(MDL_ANGLE_TILT, "tilt");
        MDL::addLabel(MDL_ANGLE_TILT2, LABEL_DOUBLE, "angleTilt2");
        MDL::addLabelAlias(MDL_ANGLE_TILT2, "tilt2");
        MDL::addLabel(MDL_ANGLE_TILT_DIFF, LABEL_DOUBLE, "angleTiltDiff");
        MDL::addLabel(MDL_ANGLE_DIFF, LABEL_DOUBLE, "angleDiff");
        MDL::addLabel(MDL_ANGLE_Y, LABEL_DOUBLE, "angleY");
        MDL::addLabel(MDL_ANGLE_Y2, LABEL_DOUBLE, "angleY2");

        MDL::addLabel(MDL_AVG, LABEL_DOUBLE, "avg");
        MDL::addLabel(MDL_AVG_CHANGES_ORIENTATIONS, LABEL_DOUBLE, "avgChanOrient");
        MDL::addLabel(MDL_AVG_CHANGES_OFFSETS, LABEL_DOUBLE, "avgChanOffset");
        MDL::addLabel(MDL_AVG_CHANGES_CLASSES, LABEL_DOUBLE, "avgChanClass");

        MDL::addLabel(MDL_BGMEAN, LABEL_DOUBLE, "bgMean");
        MDL::addLabel(MDL_BLOCK_NUMBER, LABEL_INT, "blockNumber");

        MDL::addLabel(MDL_CL2D_CHANGES, LABEL_INT, "cl2dChanges");
        MDL::addLabel(MDL_CL2D_SIMILARITY, LABEL_DOUBLE, "cl2dSimilarity");
        MDL::addLabel(MDL_CLASS_COUNT, LABEL_SIZET, "classCount");
        MDL::addLabelAlias(MDL_CLASS_COUNT, "class_count"); //3.0
        MDL::addLabel(MDL_CLASSIFICATION_DATA, LABEL_VECTOR_DOUBLE, "classificationData");
        MDL::addLabelAlias(MDL_CLASSIFICATION_DATA, "ClassificationData");
        MDL::addLabel(MDL_CLASSIFICATION_DATA_SIZE, LABEL_SIZET, "classificationDatasize");
        MDL::addLabelAlias(MDL_CLASSIFICATION_DATA_SIZE, "ClassificationDataSize");
        MDL::addLabel(MDL_CLASSIFICATION_DPR_05, LABEL_DOUBLE, "classificationDPR05");
        MDL::addLabelAlias(MDL_CLASSIFICATION_DPR_05, "ClassificationDPR05");
        MDL::addLabel(MDL_CLASSIFICATION_FRC_05, LABEL_DOUBLE, "classificationFRC05");
        MDL::addLabelAlias(MDL_CLASSIFICATION_FRC_05, "ClassificationFRC05");
        MDL::addLabel(MDL_CLASSIFICATION_INTRACLASS_DISTANCE, LABEL_DOUBLE, "classificationIntraclassDistance");
        MDL::addLabelAlias(MDL_CLASSIFICATION_INTRACLASS_DISTANCE, "ClassificationIntraclassDistance");
        MDL::addLabel(MDL_COMMENT, LABEL_STRING, "comment");
        MDL::addLabel(MDL_COST, LABEL_DOUBLE, "cost");
        MDL::addLabel(MDL_COUNT2, LABEL_SIZET, "count2");
        MDL::addLabel(MDL_COUNT, LABEL_SIZET, "count");

        MDL::addLabel(MDL_CRYSTAL_CELLX, LABEL_INT, "crystalCellx");
        MDL::addLabel(MDL_CRYSTAL_CELLY, LABEL_INT, "crystalCelly");
        MDL::addLabel(MDL_CRYSTAL_DISAPPEAR_THRE, LABEL_DOUBLE, "crystalDisthresh");
        MDL::addLabel(MDL_CRYSTAL_LATTICE_A, LABEL_VECTOR_DOUBLE, "crystalLatticeA");
        MDL::addLabel(MDL_CRYSTAL_LATTICE_B, LABEL_VECTOR_DOUBLE, "crystalLatticeB");
        MDL::addLabel(MDL_CRYSTAL_ORTHO_PRJ, LABEL_BOOL, "crystalOrthoProj");
        MDL::addLabel(MDL_CRYSTAL_PROJ, LABEL_BOOL, "crystalProj");
        MDL::addLabel(MDL_CRYSTAL_SHFILE, LABEL_STRING, "crystalShiftFile");
        MDL::addLabel(MDL_CRYSTAL_SHIFTX, LABEL_DOUBLE, "crystalShiftX");
        MDL::addLabel(MDL_CRYSTAL_SHIFTY, LABEL_DOUBLE, "crystalShiftY");
        MDL::addLabel(MDL_CRYSTAL_SHIFTZ, LABEL_DOUBLE, "crystalShiftZ");

        MDL::addLabel(MDL_CTF_BG_BASELINE, LABEL_DOUBLE, "ctfBgBaseline");
        MDL::addLabelAlias(MDL_CTF_BG_BASELINE, "CTFBG_Baseline");//3.0
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_ANGLE, LABEL_DOUBLE, "ctfBgGaussian2Angle");
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_ANGLE, "CTFBG_Gaussian2_Angle"); //3.0

        MDL::addLabel(MDL_CTF_X0, LABEL_DOUBLE, "ctfX0");
        MDL::addLabel(MDL_CTF_XF, LABEL_DOUBLE, "ctfXF");
        MDL::addLabel(MDL_CTF_Y0, LABEL_DOUBLE, "ctfY0");
        MDL::addLabel(MDL_CTF_YF, LABEL_DOUBLE, "ctfYF");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEUA, LABEL_DOUBLE, "ctfDefocusPlaneUA");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEUB, LABEL_DOUBLE, "ctfDefocusPlaneUB");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEUC, LABEL_DOUBLE, "ctfDefocusPlaneUC");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEVA, LABEL_DOUBLE, "ctfDefocusPlaneVA");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEVB, LABEL_DOUBLE, "ctfDefocusPlaneVB");
        MDL::addLabel(MDL_CTF_DEFOCUS_PLANEVC, LABEL_DOUBLE, "ctfDefocusPlaneVC");


        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_CU, LABEL_DOUBLE, "ctfBgGaussian2CU");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_CV, LABEL_DOUBLE, "ctfBgGaussian2CV");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_K, LABEL_DOUBLE, "ctfBgGaussian2K");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_SIGMAU, LABEL_DOUBLE, "ctfBgGaussian2SigmaU");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_SIGMAV, LABEL_DOUBLE, "ctfBgGaussian2SigmaV");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_ANGLE, LABEL_DOUBLE, "ctfBgGaussianAngle");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_CU, LABEL_DOUBLE, "ctfBgGaussianCU");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_CV, LABEL_DOUBLE, "ctfBgGaussianCV");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_K, LABEL_DOUBLE, "ctfBgGaussianK");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_SIGMAU, LABEL_DOUBLE, "ctfBgGaussianSigmaU");
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN_SIGMAV, LABEL_DOUBLE, "ctfBgGaussianSigmaV");

        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_CU, "CTFBG_Gaussian2_CU");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_CV, "CTFBG_Gaussian2_CV");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_K, "CTFBG_Gaussian2_K");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_SIGMAU, "CTFBG_Gaussian2_SigmaU");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_SIGMAV, "CTFBG_Gaussian2_SigmaV");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_ANGLE, "CTFBG_Gaussian_Angle");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_CU, "CTFBG_Gaussian_CU");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_CV, "CTFBG_Gaussian_CV");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_K, "CTFBG_Gaussian_K");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_SIGMAU, "CTFBG_Gaussian_SigmaU");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN_SIGMAV, "CTFBG_Gaussian_SigmaV");//3.0

        MDL::addLabel(MDL_CTF_BG_SQRT_ANGLE, LABEL_DOUBLE, "ctfBgSqrtAngle");
        MDL::addLabel(MDL_CTF_BG_SQRT_K, LABEL_DOUBLE, "ctfBgSqrtK");
        MDL::addLabel(MDL_CTF_BG_SQRT_U, LABEL_DOUBLE, "ctfBgSqrtU");
        MDL::addLabel(MDL_CTF_BG_SQRT_V, LABEL_DOUBLE, "ctfBgSqrtV");
        MDL::addLabelAlias(MDL_CTF_BG_SQRT_ANGLE, "CTFBG_Sqrt_Angle");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_SQRT_K, "CTFBG_Sqrt_K");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_SQRT_U, "CTFBG_Sqrt_U");//3.0
        MDL::addLabelAlias(MDL_CTF_BG_SQRT_V, "CTFBG_Sqrt_V");  //3.0

        MDL::addLabel(MDL_CTF_CA, LABEL_DOUBLE, "ctfChromaticAberration");
        MDL::addLabel(MDL_CTF_CONVERGENCE_CONE, LABEL_DOUBLE, "ctfConvergenceCone");
        MDL::addLabel(MDL_CTF_CRIT_DAMPING, LABEL_DOUBLE, "ctfCritDamping");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZEROAVG, LABEL_DOUBLE, "ctfCritFirstZero");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZERODISAGREEMENT, LABEL_DOUBLE, "ctfCritDisagree");
        MDL::addLabel(MDL_CTF_CRIT_MAXFREQ, LABEL_DOUBLE, "ctfCritMaxFreq");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZERORATIO, LABEL_DOUBLE, "ctfCritfirstZeroRatio");
        MDL::addLabel(MDL_CTF_CRIT_FITTINGCORR13, LABEL_DOUBLE, "ctfCritCorr13");
        MDL::addLabel(MDL_CTF_CRIT_FITTINGSCORE, LABEL_DOUBLE, "ctfCritFitting");
        MDL::addLabel(MDL_CTF_CRIT_NORMALITY, LABEL_DOUBLE, "ctfCritNormality");
        MDL::addLabel(MDL_CTF_CRIT_PSDCORRELATION90, LABEL_DOUBLE, "ctfCritPsdCorr90");
        MDL::addLabel(MDL_CTF_CRIT_PSDPCA1VARIANCE, LABEL_DOUBLE, "ctfCritPsdPCA1");
        MDL::addLabel(MDL_CTF_CRIT_PSDPCARUNSTEST, LABEL_DOUBLE, "ctfCritPsdPCARuns");
        MDL::addLabel(MDL_CTF_CRIT_PSDRADIALINTEGRAL, LABEL_DOUBLE, "ctfCritPsdInt");
        MDL::addLabel(MDL_CTF_CRIT_PSDVARIANCE, LABEL_DOUBLE, "ctfCritPsdStdQ");
        MDL::addLabel(MDL_CTF_CS, LABEL_DOUBLE, "ctfSphericalAberration");
        MDL::addLabel(MDL_CTF_DEFOCUSA, LABEL_DOUBLE, "ctfDefocusA");//average defocus
        MDL::addLabel(MDL_CTF_DEFOCUS_ANGLE, LABEL_DOUBLE, "ctfDefocusAngle");
        MDL::addLabel(MDL_CTF_DEFOCUSU, LABEL_DOUBLE, "ctfDefocusU");
        MDL::addLabel(MDL_CTF_DEFOCUSV, LABEL_DOUBLE, "ctfDefocusV");
        MDL::addLabel(MDL_CTF_DIMENSIONS, LABEL_VECTOR_DOUBLE, "ctfDimensions");
        MDL::addLabel(MDL_CTF_DOWNSAMPLE_PERFORMED, LABEL_DOUBLE, "CtfDownsampleFactor");
        MDL::addLabel(MDL_CTF_ENERGY_LOSS, LABEL_DOUBLE, "ctfEnergyLoss");
        MDL::addLabel(MDL_CTF_GROUP, LABEL_INT, "ctfGroup");
        MDL::addLabel(MDL_CTF_INPUTPARAMS, LABEL_STRING, "ctfInputParams", TAGLABEL_TEXTFILE);
        MDL::addLabel(MDL_CTF_K, LABEL_DOUBLE, "ctfK");
        MDL::addLabel(MDL_CTF_LAMBDA, LABEL_DOUBLE, "ctfLambda");
        MDL::addLabel(MDL_CTF_LENS_STABILITY, LABEL_DOUBLE, "ctfLensStability");
        MDL::addLabel(MDL_CTF_LONGITUDINAL_DISPLACEMENT, LABEL_DOUBLE, "ctfLongitudinalDisplacement");
        MDL::addLabel(MDL_CTF_MODEL2, LABEL_STRING, "ctfModel2", TAGLABEL_CTFPARAM);
        MDL::addLabel(MDL_CTF_MODEL, LABEL_STRING, "ctfModel", TAGLABEL_CTFPARAM);
        MDL::addLabel(MDL_CTF_Q0, LABEL_DOUBLE, "ctfQ0");
        MDL::addLabel(MDL_CTF_SAMPLING_RATE, LABEL_DOUBLE, "ctfSamplingRate");
        MDL::addLabel(MDL_CTF_SAMPLING_RATE_Z, LABEL_DOUBLE, "ctfSamplingRateZ");
        MDL::addLabel(MDL_CTF_TRANSVERSAL_DISPLACEMENT, LABEL_DOUBLE, "ctfTransversalDisplacement");
        MDL::addLabel(MDL_CTF_VOLTAGE, LABEL_DOUBLE, "ctfVoltage");
        MDL::addLabel(MDL_CTF_XRAY_LENS_TYPE, LABEL_STRING, "ctfXrayLensType");
        MDL::addLabel(MDL_CTF_XRAY_OUTER_ZONE_WIDTH, LABEL_DOUBLE, "ctfXrayOuterZoneWidth");
        MDL::addLabel(MDL_CTF_XRAY_ZONES_NUMBER, LABEL_DOUBLE, "ctfXrayZonesN");

        MDL::addLabelAlias(MDL_CTF_CA, "CTF_Chromatic_aberration"); //3.0
        MDL::addLabelAlias(MDL_CTF_CONVERGENCE_CONE, "CTF_Convergence_cone"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_DAMPING, "CTFCrit_damping"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FIRSTZEROAVG, "CTFCrit_firstZero"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FIRSTZERODISAGREEMENT, "CTFCrit_disagree"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FIRSTZERORATIO, "CTFCrit_firstZeroRatio"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FITTINGCORR13, "CTFCrit_Corr13"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_FITTINGSCORE, "CTFCrit_Fitting"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_NORMALITY, "CTFCrit_Normality"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDCORRELATION90, "CTFCrit_psdcorr90"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDPCA1VARIANCE, "CTFCrit_PSDPCA1"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDPCARUNSTEST, "CTFCrit_PSDPCARuns"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDRADIALINTEGRAL, "CTFCrit_psdint"); //3.0
        MDL::addLabelAlias(MDL_CTF_CRIT_PSDVARIANCE, "CTFCrit_PSDStdQ"); //3.0
        MDL::addLabelAlias(MDL_CTF_CS, "CTF_Spherical_aberration"); //3.0
        MDL::addLabelAlias(MDL_CTF_DEFOCUSA, "CTF_Defocus_A"); //3.0//average defocus
        MDL::addLabelAlias(MDL_CTF_DEFOCUS_ANGLE, "CTF_Defocus_angle"); //3.0
        MDL::addLabelAlias(MDL_CTF_DEFOCUSU, "CTF_Defocus_U"); //3.0
        MDL::addLabelAlias(MDL_CTF_DEFOCUSV, "CTF_Defocus_V"); //3.0
        MDL::addLabelAlias(MDL_CTF_DIMENSIONS, "CTF_Xray_dimensions"); //3.0
        MDL::addLabelAlias(MDL_CTF_DOWNSAMPLE_PERFORMED, "CTFDownsampleFactor"); //3.0
        MDL::addLabelAlias(MDL_CTF_ENERGY_LOSS, "CTF_Energy_loss"); //3.0
        MDL::addLabelAlias(MDL_CTF_GROUP, "CTFGroup"); //3.0
        MDL::addLabelAlias(MDL_CTF_INPUTPARAMS, "CTFInputParams"); //3.0
        MDL::addLabelAlias(MDL_CTF_K, "CTF_K"); //3.0
        MDL::addLabelAlias(MDL_CTF_LAMBDA, "CTF_Xray_lambda"); //3.0
        MDL::addLabelAlias(MDL_CTF_LENS_STABILITY, "CTF_Lens_stability"); //3.0
        MDL::addLabelAlias(MDL_CTF_LONGITUDINAL_DISPLACEMENT, "CTF_Longitudinal_displacement"); //3.0
        MDL::addLabelAlias(MDL_CTF_MODEL2, "CTFModel2"); //3.0
        MDL::addLabelAlias(MDL_CTF_MODEL, "CTFModel"); //3.0
        MDL::addLabelAlias(MDL_CTF_Q0, "CTF_Q0"); //3.0
        MDL::addLabelAlias(MDL_CTF_SAMPLING_RATE, "CTF_Sampling_rate"); //3.0
        MDL::addLabelAlias(MDL_CTF_SAMPLING_RATE_Z, "CTF_Sampling_rate_z"); //3.0
        MDL::addLabelAlias(MDL_CTF_TRANSVERSAL_DISPLACEMENT, "CTF_Transversal_displacement"); //3.0
        MDL::addLabelAlias(MDL_CTF_VOLTAGE, "CTF_Voltage"); //3.0
        MDL::addLabelAlias(MDL_CTF_XRAY_LENS_TYPE, "CTF_Xray_lens_type"); //3.0
        MDL::addLabelAlias(MDL_CTF_XRAY_OUTER_ZONE_WIDTH, "CTF_Xray_OuterZoneWidth"); //3.0
        MDL::addLabelAlias(MDL_CTF_XRAY_ZONES_NUMBER, "CTF_Xray_ZonesN"); //3.0

        MDL::addLabel(MDL_DATATYPE, LABEL_INT, "datatype");
        MDL::addLabel(MDL_DEFGROUP, LABEL_INT, "defocusGroup");
        MDL::addLabel(MDL_DIMENSIONS_2D, LABEL_VECTOR_DOUBLE, "dimensions2D");
        MDL::addLabel(MDL_DIMENSIONS_3D, LABEL_VECTOR_DOUBLE, "dimensions3D");
        MDL::addLabel(MDL_DIMRED, LABEL_VECTOR_DOUBLE, "dimredCoeffs");
        MDL::addLabel(MDL_DIRECTION, LABEL_VECTOR_DOUBLE, "direction");
        MDL::addLabel(MDL_DM3_IDTAG, LABEL_INT, "dm3IdTag");
        MDL::addLabel(MDL_DM3_NODEID, LABEL_INT, "dm3NodeId");
        MDL::addLabel(MDL_DM3_NUMBER_TYPE, LABEL_INT, "dm3NumberType");
        MDL::addLabel(MDL_DM3_PARENTID, LABEL_INT, "dm3ParentID");
        MDL::addLabel(MDL_DM3_SIZE, LABEL_INT, "dm3Size");
        MDL::addLabel(MDL_DM3_TAGCLASS, LABEL_STRING, "dm3TagClass");
        MDL::addLabel(MDL_DM3_TAGNAME, LABEL_STRING, "dm3TagName");
        MDL::addLabel(MDL_DM3_VALUE, LABEL_VECTOR_DOUBLE, "dm3Value");

        MDL::addLabel(MDL_ENABLED, LABEL_INT, "enabled");

        //MDL_EXECUTION_DATE so far an string but may change...
        MDL::addLabel(MDL_DATE, LABEL_STRING, "date");
        MDL::addLabel(MDL_TIME, LABEL_DOUBLE, "time");

        MDL::addLabel(MDL_FLIP, LABEL_BOOL, "flip");
        MDL::addLabelAlias(MDL_FLIP, "Flip");
        MDL::addLabel(MDL_FOM, LABEL_DOUBLE, "fom");
        MDL::addLabel(MDL_IDX, LABEL_SIZET, "index");
        MDL::addLabel(MDL_IMAGE1, LABEL_STRING, "image1", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE2, LABEL_STRING, "image2", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE3, LABEL_STRING, "image3", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE4, LABEL_STRING, "image4", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE5, LABEL_STRING, "image5", TAGLABEL_IMAGE);

        MDL::addLabelAlias(MDL_IMAGE1, "associatedImage1"); //3.0
        MDL::addLabelAlias(MDL_IMAGE2, "associatedImage2"); //3.0
        MDL::addLabelAlias(MDL_IMAGE3, "associatedImage3"); //3.0
        MDL::addLabelAlias(MDL_IMAGE4, "associatedImage4"); //3.0
        MDL::addLabelAlias(MDL_IMAGE5, "associatedImage5"); //3.0

        MDL::addLabel(MDL_IMAGE, LABEL_STRING, "image", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_ORIGINAL, LABEL_STRING, "imageOriginal", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_REF, LABEL_STRING, "imageRef", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_TILTED, LABEL_STRING, "imageTilted", TAGLABEL_IMAGE);

        MDL::addLabelAlias(MDL_IMAGE_ORIGINAL, "original_image"); //3.0
        MDL::addLabelAlias(MDL_IMAGE_TILTED, "tilted_image"); //3.0

        MDL::addLabel(MDL_IMGMD, LABEL_STRING, "imageMetaData", TAGLABEL_METADATA);
        MDL::addLabel(MDL_INTSCALE, LABEL_DOUBLE, "intScale");

        MDL::addLabel(MDL_ITEM_ID, LABEL_SIZET, "itemId");
        MDL::addLabel(MDL_ITER, LABEL_INT, "iterationNumber");

        MDL::addLabel(MDL_KERDENSOM_FUNCTIONAL, LABEL_DOUBLE, "kerdensomFunctional");
        MDL::addLabel(MDL_KERDENSOM_REGULARIZATION, LABEL_DOUBLE, "kerdensomRegularization");
        MDL::addLabel(MDL_KERDENSOM_SIGMA, LABEL_DOUBLE, "kerdensomSigma");

        MDL::addLabel(MDL_KEYWORDS, LABEL_STRING, "keywords");
        MDL::addLabel(MDL_KSTEST, LABEL_DOUBLE, "kstest");
        MDL::addLabel(MDL_LL, LABEL_DOUBLE, "logLikelihood");
        MDL::addLabelAlias(MDL_LL, "LL");
        MDL::addLabel(MDL_MAGNIFICATION, LABEL_DOUBLE, "magnification");
        MDL::addLabel(MDL_MAPTOPOLOGY, LABEL_STRING, "mapTopology");
        MDL::addLabel(MDL_MASK, LABEL_STRING, "mask", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_MAXCC, LABEL_DOUBLE, "maxCC");
        MDL::addLabel(MDL_MAX, LABEL_DOUBLE, "max");
        MDL::addLabel(MDL_MICROGRAPH, LABEL_STRING, "micrograph", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_PARTICLES, LABEL_STRING, "micrographParticles", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_ORIGINAL, LABEL_STRING, "micrographOriginal", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_TILTED, LABEL_STRING, "micrographTilted", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_TILTED_ORIGINAL, LABEL_STRING, "micrographTiltedOriginal", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MIN, LABEL_DOUBLE, "min");
        MDL::addLabel(MDL_MIRRORFRAC, LABEL_DOUBLE, "mirrorFraction");
        MDL::addLabel(MDL_MISSINGREGION_NR, LABEL_INT, "missingRegionNumber");
        MDL::addLabelAlias(MDL_MISSINGREGION_NR, "Wedge");
        MDL::addLabel(MDL_MISSINGREGION_THX0, LABEL_DOUBLE, "missingRegionThetaX0");
        MDL::addLabel(MDL_MISSINGREGION_THXF, LABEL_DOUBLE, "missingRegionThetaXF");

        MDL::addLabel(MDL_MLF_CTF,    LABEL_DOUBLE, "mlfCtf");
        MDL::addLabel(MDL_MLF_WIENER, LABEL_DOUBLE, "mlfWiener");
        MDL::addLabel(MDL_MLF_SIGNAL, LABEL_DOUBLE, "mlfSignal");
        MDL::addLabel(MDL_MLF_NOISE,  LABEL_DOUBLE, "mlfNoise");

        MDL::addLabel(MDL_MISSINGREGION_THY0, LABEL_DOUBLE, "missingRegionThetaY0");
        MDL::addLabel(MDL_MISSINGREGION_THYF, LABEL_DOUBLE, "missingRegionThetaYF");
        MDL::addLabel(MDL_MISSINGREGION_TYPE, LABEL_STRING, "missingRegionType");
        MDL::addLabel(MDL_MODELFRAC, LABEL_DOUBLE, "modelFraction");
        MDL::addLabel(MDL_NEIGHBORHOOD_RADIUS, LABEL_DOUBLE, "neighborhoodRadius");
        MDL::addLabel(MDL_NEIGHBOR, LABEL_SIZET, "neighbor");
        MDL::addLabel(MDL_NEIGHBORS, LABEL_VECTOR_SIZET, "neighbors");
        MDL::addLabel(MDL_NMA, LABEL_VECTOR_DOUBLE, "nmaDisplacements");
        MDL::addLabelAlias(MDL_NMA, "NMADisplacements");//3.0
        MDL::addLabel(MDL_NMA_MODEFILE, LABEL_STRING, "nmaModefile", TAGLABEL_TEXTFILE);
        MDL::addLabelAlias(MDL_NMA_MODEFILE, "NMAModefile");//3.0
        MDL::addLabel(MDL_NMA_COLLECTIVITY, LABEL_DOUBLE, "nmaCollectivity");
        MDL::addLabel(MDL_NMA_SCORE, LABEL_DOUBLE, "nmaScore");
        MDL::addLabel(MDL_NOISE_ANGLES, LABEL_VECTOR_DOUBLE, "noiseAngles");
        MDL::addLabel(MDL_NOISE_COORD, LABEL_VECTOR_DOUBLE, "noiseCoord");
        MDL::addLabel(MDL_NOISE_PARTICLE_COORD, LABEL_VECTOR_DOUBLE, "noiseParticleCoord");
        MDL::addLabel(MDL_NOISE_PIXEL_LEVEL, LABEL_VECTOR_DOUBLE, "noisePixelLevel");
        MDL::addLabel(MDL_ORDER, LABEL_SIZET, "order_");
        MDL::addLabel(MDL_ORIGIN_X, LABEL_DOUBLE, "originX");
        MDL::addLabel(MDL_ORIGIN_Y, LABEL_DOUBLE, "originY");
        MDL::addLabel(MDL_ORIGIN_Z, LABEL_DOUBLE, "originZ");
        MDL::addLabel(MDL_PHANTOM_BGDENSITY, LABEL_DOUBLE, "phantomBGDensity");
        MDL::addLabel(MDL_PHANTOM_FEATURE_CENTER, LABEL_VECTOR_DOUBLE, "featureCenter");
        MDL::addLabel(MDL_PHANTOM_FEATURE_DENSITY, LABEL_DOUBLE, "featureDensity");
        MDL::addLabel(MDL_PHANTOM_FEATURE_OPERATION, LABEL_STRING, "featureOperation");
        MDL::addLabel(MDL_PHANTOM_FEATURE_SPECIFIC, LABEL_VECTOR_DOUBLE, "featureSpecificVector");
        MDL::addLabel(MDL_PHANTOM_FEATURE_TYPE, LABEL_STRING, "featureType");
        MDL::addLabel(MDL_MACRO_CMD, LABEL_STRING, "macroCmd");
        MDL::addLabel(MDL_MACRO_CMD_ARGS, LABEL_STRING, "macroCmdArgs");
        MDL::addLabel(MDL_COLOR, LABEL_INT, "color");
        MDL::addLabel(MDL_PICKING_STATE, LABEL_STRING, "pickingState");
        MDL::addLabelAlias(MDL_PICKING_STATE, "picking_state");//3.0
        MDL::addLabel(MDL_PICKING_MICROGRAPH_STATE, LABEL_STRING, "pickingMicrographState");
        MDL::addLabelAlias(MDL_PICKING_MICROGRAPH_STATE, "micrograph_state");//3.0
        MDL::addLabel(MDL_PICKING_PARTICLE_SIZE, LABEL_INT, "particleSize");
        MDL::addLabel(MDL_PICKING_AUTOPICKPERCENT, LABEL_INT, "autopickPercent");
        MDL::addLabel(MDL_PICKING_TEMPLATES, LABEL_INT, "templatesNum");
        MDL::addLabel(MDL_PMAX, LABEL_DOUBLE, "pMax");
        MDL::addLabel(MDL_AVGPMAX, LABEL_DOUBLE, "pMax");
        MDL::addLabelAlias(MDL_PMAX, "Pmax");
        MDL::addLabelAlias(MDL_PMAX, "sumP");
        MDL::addLabel(MDL_POINTSASYMETRICUNIT, LABEL_SIZET, "pointsAsymmetricUnit");
        MDL::addLabelAlias(MDL_POINTSASYMETRICUNIT, "pointsasymmetricUnit");
        MDL::addLabel(MDL_PRJ_ANGFILE, LABEL_STRING, "projAngleFile");
        MDL::addLabelAlias(MDL_PRJ_ANGFILE, "angleFile");//3.0
        MDL::addLabel(MDL_PRJ_DIMENSIONS, LABEL_VECTOR_DOUBLE, "projDimensions");
        MDL::addLabel(MDL_PRJ_PSI_NOISE, LABEL_VECTOR_DOUBLE, "projPsiNoise");
        MDL::addLabel(MDL_PRJ_PSI_RANDSTR, LABEL_STRING, "projPsiRandomness");
        MDL::addLabel(MDL_PRJ_PSI_RANGE, LABEL_VECTOR_DOUBLE, "projPsiRange");
        MDL::addLabel(MDL_PRJ_ROT_NOISE, LABEL_VECTOR_DOUBLE, "projRotNoise");
        MDL::addLabel(MDL_PRJ_ROT_RANDSTR, LABEL_STRING, "projRotRandomness");
        MDL::addLabel(MDL_PRJ_ROT_RANGE, LABEL_VECTOR_DOUBLE, "projRotRange");
        MDL::addLabel(MDL_PRJ_TILT_NOISE, LABEL_VECTOR_DOUBLE, "projTiltNoise");
        MDL::addLabel(MDL_PRJ_TILT_RANDSTR, LABEL_STRING, "projTiltRandomness");
        MDL::addLabel(MDL_PRJ_TILT_RANGE, LABEL_VECTOR_DOUBLE, "projTiltRange");
        MDL::addLabel(MDL_PRJ_VOL, LABEL_STRING, "projVolume", TAGLABEL_VOLUME);

        MDL::addLabel(MDL_PROGRAM, LABEL_STRING, "program");
        MDL::addLabel(MDL_USER, LABEL_STRING, "user");

        MDL::addLabel(MDL_PSD_ENHANCED, LABEL_STRING, "psdEnhanced", TAGLABEL_IMAGE);
        MDL::addLabelAlias(MDL_PSD_ENHANCED, "enhancedPowerSpectrum");//3.0
        MDL::addLabel(MDL_PSD, LABEL_STRING, "psd", TAGLABEL_PSD);
        MDL::addLabelAlias(MDL_PSD, "powerSpectrum");//3.0
        MDL::addLabel(MDL_RANDOMSEED, LABEL_INT, "randomSeed");
        MDL::addLabel(MDL_REF2, LABEL_INT, "ref2");
        MDL::addLabel(MDL_REF3D, LABEL_INT, "ref3d");
        MDL::addLabel(MDL_REF, LABEL_INT, "ref");
        MDL::addLabelAlias(MDL_REF, "Ref");
        MDL::addLabel(MDL_REFMD, LABEL_STRING, "referenceMetaData", TAGLABEL_METADATA);

        MDL::addLabel(MDL_RESOLUTION_DPR, LABEL_DOUBLE, "resolutionDPR");
        MDL::addLabel(MDL_RESOLUTION_ERRORL2, LABEL_DOUBLE, "resolutionErrorL2");
        MDL::addLabel(MDL_RESOLUTION_FRC, LABEL_DOUBLE, "resolutionFRC");
        MDL::addLabel(MDL_RESOLUTION_FRCRANDOMNOISE, LABEL_DOUBLE, "resolutionFRCRandomNoise");
        MDL::addLabel(MDL_RESOLUTION_FREQ, LABEL_DOUBLE, "resolutionFreqFourier");
        MDL::addLabel(MDL_RESOLUTION_FREQREAL, LABEL_DOUBLE, "resolutionFreqReal");
        MDL::addLabel(MDL_RESOLUTION_SSNR, LABEL_DOUBLE, "resolutionSSNR");

        MDL::addLabelAlias(MDL_RESOLUTION_DPR, "DPR");
        MDL::addLabelAlias(MDL_RESOLUTION_ERRORL2, "Error_l2");
        MDL::addLabelAlias(MDL_RESOLUTION_FRC, "FRC");
        MDL::addLabelAlias(MDL_RESOLUTION_FRCRANDOMNOISE, "FRC_random_noise");
        MDL::addLabelAlias(MDL_RESOLUTION_FREQ, "Resol_Inverse_Ang");
        MDL::addLabelAlias(MDL_RESOLUTION_FREQREAL, "Resol_Ang");

        MDL::addLabel(MDL_SAMPLINGRATE, LABEL_DOUBLE, "samplingRate");
        MDL::addLabel(MDL_SAMPLINGRATE_ORIGINAL, LABEL_DOUBLE, "samplingRateOriginal");
        MDL::addLabel(MDL_SAMPLINGRATE_X, LABEL_DOUBLE, "samplingRateX");
        MDL::addLabel(MDL_SAMPLINGRATE_Y, LABEL_DOUBLE, "samplingRateY");
        MDL::addLabel(MDL_SAMPLINGRATE_Z, LABEL_DOUBLE, "samplingRateZ");

        MDL::addLabelAlias(MDL_SAMPLINGRATE, "sampling_rate"); //3.0
        MDL::addLabelAlias(MDL_SAMPLINGRATE_ORIGINAL, "sampling_rate_original"); //3.0
        MDL::addLabelAlias(MDL_SAMPLINGRATE_X, "sampling_rateX"); //3.0
        MDL::addLabelAlias(MDL_SAMPLINGRATE_Y, "sampling_rateY"); //3.0
        MDL::addLabelAlias(MDL_SAMPLINGRATE_Z, "sampling_rateZ"); //3.0

        MDL::addLabel(MDL_SCALE, LABEL_DOUBLE, "scale");
        MDL::addLabelAlias(MDL_SCALE, "Scale");
        MDL::addLabel(MDL_SELFILE, LABEL_STRING, "selfile", TAGLABEL_METADATA);
        MDL::addLabel(MDL_SERIE, LABEL_STRING, "serie");

        MDL::addLabel(MDL_SHIFT_X, LABEL_DOUBLE, "shiftX");
        MDL::addLabelAlias(MDL_SHIFT_X, "Xoff");
        MDL::addLabel(MDL_SHIFT_X2, LABEL_DOUBLE, "shiftX2");
        MDL::addLabel(MDL_SHIFT_X_DIFF, LABEL_DOUBLE, "shiftXDiff");
        MDL::addLabel(MDL_SHIFT_Y, LABEL_DOUBLE, "shiftY");
        MDL::addLabelAlias(MDL_SHIFT_Y, "Yoff");
        MDL::addLabel(MDL_SHIFT_Y2, LABEL_DOUBLE, "shiftY2");
        MDL::addLabel(MDL_SHIFT_Y_DIFF, LABEL_DOUBLE, "shiftYDiff");
        MDL::addLabel(MDL_SHIFT_Z, LABEL_DOUBLE, "shiftZ");
        MDL::addLabelAlias(MDL_SHIFT_Z, "Zoff");
        MDL::addLabel(MDL_SHIFT_DIFF, LABEL_DOUBLE, "shiftDiff");
        MDL::addLabel(MDL_SIGMANOISE, LABEL_DOUBLE, "sigmaNoise");
        MDL::addLabel(MDL_SIGMAOFFSET, LABEL_DOUBLE, "sigmaOffset");
        MDL::addLabel(MDL_SIGNALCHANGE, LABEL_DOUBLE, "signalChange");
        MDL::addLabel(MDL_STDDEV, LABEL_DOUBLE, "stddev");
        MDL::addLabel(MDL_SUM, LABEL_DOUBLE, "sum");
        MDL::addLabel(MDL_SUMWEIGHT, LABEL_DOUBLE, "sumWeight");
        MDL::addLabel(MDL_SYMNO, LABEL_INT, "symNo");
        MDL::addLabel(MDL_TRANSFORMATIONMTRIX, LABEL_VECTOR_DOUBLE, "transMat");

        MDL::addLabel(MDL_VOLTAGE, LABEL_DOUBLE, "voltage");
        MDL::addLabel(MDL_VOLUME_SCORE1, LABEL_DOUBLE, "volScore1");
        MDL::addLabel(MDL_VOLUME_SCORE2, LABEL_DOUBLE, "volScore2");
        MDL::addLabel(MDL_VOLUME_SCORE3, LABEL_DOUBLE, "volScore3");
        MDL::addLabel(MDL_VOLUME_SCORE4, LABEL_DOUBLE, "volScore4");
        MDL::addLabel(MDL_WEIGHT, LABEL_DOUBLE, "weight");
        MDL::addLabelAlias(MDL_WEIGHT, "Weight");
        MDL::addLabel(MDL_WROBUST, LABEL_DOUBLE, "wRobust");
        MDL::addLabel(MDL_XCOOR, LABEL_INT, "xcoor");
        MDL::addLabel(MDL_XCOOR_TILT, LABEL_INT, "xcoorTilt");

        MDL::addLabel(MDL_X, LABEL_DOUBLE, "x");
        MDL::addLabel(MDL_XSIZE, LABEL_SIZET, "xSize");
        MDL::addLabel(MDL_YCOOR, LABEL_INT, "ycoor");
        MDL::addLabel(MDL_YCOOR_TILT, LABEL_INT, "ycoorTilt");
        MDL::addLabel(MDL_Y, LABEL_DOUBLE, "y");
        MDL::addLabel(MDL_YSIZE, LABEL_SIZET, "ySize");
        MDL::addLabel(MDL_ZCOOR, LABEL_INT, "zcoor");
        MDL::addLabel(MDL_Z, LABEL_DOUBLE, "z");
        MDL::addLabel(MDL_ZSCORE, LABEL_DOUBLE, "zScore");
        MDL::addLabel(MDL_ZSCORE_SHAPE1, LABEL_DOUBLE, "zScoreShape1");
        MDL::addLabel(MDL_ZSCORE_SHAPE2, LABEL_DOUBLE, "zScoreShape2");
        MDL::addLabel(MDL_ZSCORE_SNR1, LABEL_DOUBLE, "zScoreSNR1");
        MDL::addLabel(MDL_ZSCORE_SNR2, LABEL_DOUBLE, "zScoreSNR2");
        MDL::addLabel(MDL_ZSCORE_HISTOGRAM, LABEL_DOUBLE, "zScoreHistogram");
        MDL::addLabel(MDL_ZSIZE, LABEL_SIZET, "zSize");

        MDL::addLabelAlias(MDL_XCOOR, "Xcoor");//3.0
        MDL::addLabelAlias(MDL_XCOOR, "<X position>");
        MDL::addLabelAlias(MDL_XCOOR_TILT, "XcoorTilt");//3.0

        MDL::addLabelAlias(MDL_X, "X"); //3.0
        MDL::addLabelAlias(MDL_XSIZE, "Xsize"); //3.0
        MDL::addLabelAlias(MDL_YCOOR, "Ycoor"); //3.0
        MDL::addLabelAlias(MDL_YCOOR, "<Y position>");
        MDL::addLabelAlias(MDL_YCOOR_TILT, "YcoorTilt"); //3.0
        MDL::addLabelAlias(MDL_Y, "Y"); //3.0
        MDL::addLabelAlias(MDL_YSIZE, "Ysize"); //3.0
        MDL::addLabelAlias(MDL_ZCOOR, "Zcoor"); //3.0
        MDL::addLabelAlias(MDL_Z, "Z"); //3.0
        MDL::addLabelAlias(MDL_ZSCORE, "Zscore"); //3.0
        MDL::addLabelAlias(MDL_ZSIZE, "Zsize"); //3.0

        //Create an static empty header for image initialization
        MDL::emptyHeader.resetGeo();
        MDL::emptyHeader.setValue(MDL_ANGLE_ROT, 0.);
        MDL::emptyHeader.setValue(MDL_ANGLE_TILT,0.);

    }

    ~MDLabelStaticInit()
    {
        //Free memory allocated for labels data
        for (int i = MDL_FIRST_LABEL; i < MDL_LAST_LABEL; ++i)
            delete MDL::data[i];
    }
    friend class MDL;
};
/** @} */

#endif
