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

    MDL_ANGLE_COMPARISON, ///< Angular comparison (see angular_distance.cpp)
    MDL_ANGLEPSI, ///< Psi angle of an image (double,degrees)
    MDL_ANGLEPSI2, ///< Psi angle of an image (double,degrees)
    MDL_ANGLEROT, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLEROT2, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLETILT, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLETILT2, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLE_Y,   ///< Angle between y-axis and tilt-axis (double, degrees) for untilted micrographs
    MDL_ANGLE_Y2,   ///< Angle between y-axis and tilt-axis (double, degrees) for tilted micrographs
    MDL_ASSOCIATED_IMAGE1, ///< Image associated to this object (std::string)
    MDL_ASSOCIATED_IMAGE2, ///< Image associated to this object (std::string)
    MDL_ASSOCIATED_IMAGE3, ///< Image associated to this object (std::string)
    MDL_ASSOCIATED_IMAGE4, ///< Image associated to this object (std::string)
    MDL_ASSOCIATED_IMAGE5, ///< Image associated to this object (std::string)
    MDL_AVG, ///< average value (double)
    MDL_AZIMUTALANGLE, ///< ctf definition azimutal angle
    MDL_BGMEAN, ///< Mean background value for an image
    MDL_BLOCK, ///< Current block number (for incremental EM)
    MDL_CELLX, ///< Cell location for crystals
    MDL_CELLY, ///< Cell location for crystals
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
    MDL_CTFINPUTPARAMS, ///< Parameters file for the CTF Model (std::string)
    MDL_CTFMODEL, ///< Name for the CTF Model (std::string)
    MDL_CTFMODEL2, ///< Name for another CTF model (std::string)
    MDL_CTF_SAMPLING_RATE, ///< Sampling rate
    MDL_CTF_SAMPLING_RATE_Z, ///< Sampling rate in Z direction
    MDL_CTF_VOLTAGE, ///< Microscope voltage (kV)
    MDL_CTF_DEFOCUSA, ///< aver (Angage defocusstroms)
    MDL_CTF_DEFOCUSU, ///< Defocus U (Angstroms)
    MDL_CTF_DEFOCUSV, ///< Defocus V (Angstroms)
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
    MDL_CTFBG_GAUSSIAN_K, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN_SIGMAU, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN_SIGMAV, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN_CU, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN_CV, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN_ANGLE, ///< CTF Background parameter
    MDL_CTFBG_SQRT_K, ///< CTF Background parameter
    MDL_CTFBG_SQRT_U, ///< CTF Background parameter
    MDL_CTFBG_SQRT_V, ///< CTF Background parameter
    MDL_CTFBG_SQRT_ANGLE, ///< CTF Background parameter
    MDL_CTFBG_BASELINE, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN2_K, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN2_SIGMAU, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN2_SIGMAV, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN2_CU, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN2_CV, ///< CTF Background parameter
    MDL_CTFBG_GAUSSIAN2_ANGLE, ///< CTF Background parameter
    MDL_CTF_CRITERION_PSDCORRELATION90, ///< PSD correlation at 90 degrees
    MDL_CTF_CRITERION_FIRSTZERORATIO, ///< First zero ratio
    MDL_CTF_CRITERION_FIRSTZEROAVG, ///< First zero average (in Angstroms)
    MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT, ///< First zero disagreement with second model (in Angstroms)
    MDL_CTF_CRITERION_DAMPING, ///< Minimum damping at border
    MDL_CTF_CRITERION_PSDRADIALINTEGRAL, ///< Integral of the radial PSD
    MDL_CTF_CRITERION_FITTINGSCORE, ///< Score of the fitting
    MDL_CTF_CRITERION_FITTINGCORR13, ///< Correlation between the 1st and 3rd ring of the CTF
    MDL_CTF_CRITERION_PSDVARIANCE, ///< PSD variance
    MDL_CTF_CRITERION_PSDPCA1VARIANCE, ///< Variance in the first principal component of the PSDs
    MDL_CTF_CRITERION_PSDPCARUNSTEST, ///< Runs test on the projection of the PSD on the first principal component
    MDL_CTF_CRITERION_NORMALITY, ///< Normality test between histogram of micrography and gaussian distribution
    MDL_CTF_DOWNSAMPLE_PERFORMED, ///< Downsampling performed to estimate the CTF
    MDL_CTF_XRAY_DIMENSIONS, // Size in pixels of the 3D PSF to be created (Xdim, Ydim, Zdim)
    MDL_CTF_XRAY_LAMBDA, /// X-ray wavelength (nm)
    MDL_CTF_XRAY_LENS_TYPE, ///Algorithm used to generate Xray PSF
    MDL_MAGNIFICATION, /// Magnification of the X-ray microscope
    MDL_CTF_XRAY_OUTER_ZONE_WIDTH, /// Outermost zone width of the X-ray Fresnel lens (nm)
    MDL_CTF_XRAY_ZONES_NUMBER, // Number of zones of the X-ray Fresnel lens
    MDL_DATATYPE, ///< if read from file original image datatype, this is an struct defined in image
    MDL_DEFGROUP, ///< Defocus group
    MDL_DIRECTION, ///< Direction in 3D

    MDL_DM3_IDTAG,
    MDL_DM3_NODEID,
    MDL_DM3_NUMBER_TYPE,
    MDL_DM3_PARENTID,
    MDL_DM3_TAGCLASS,
    MDL_DM3_TAGNAME,
    MDL_DM3_SIZE,
    MDL_DM3_VALUE,

    MDL_ENABLED, ///< Is this image enabled? (int [-1 or 1])
    MDL_FLIP, ///< Flip the image? (bool)
    MDL_FOM, ///< Figure of Merit in 0-1 range (double)
    MDL_IDX, ///< Index within a list (size_t)
    MDL_IMAGE, ///< Name of an image (std::string)
    MDL_IMAGE_ORIGINAL, ///< Name of an image from which MDL_IMAGE is coming from
    MDL_IMAGE_TILTED, ///< Name of the tilted images associated to MDL_IMAGE
    MDL_IMGMD, ///< Name of Metadata file for all images (string)
    MDL_INTSCALE, ///< Intensity scale for an image
    MDL_ITER, ///< Current iteration number (int)
    MDL_K, ///< //ctf definition K
    MDL_KERDENSOM_FUNCTIONAL, ///< Functional value (double)
    MDL_KERDENSOM_REGULARIZATION, ///< Regularization value (double)
    MDL_KERDENSOM_SIGMA, ///< Sigma value (double)
    MDL_KEYWORDS, ///< Keywords associated to this line, should be a single string block (do not use spaces as separators)
    MDL_KSTEST, ///<KS-test statistics
    MDL_LL, ///< contribution of an image to log-likelihood value
    MDL_MAPTOPOLOGY, ///< Map topology (KerDenSOM, ...)
    MDL_MASK, ///< Name of a mask associated to image
    MDL_MAXCC, ///< Maximum cross-correlation for the image (double)
    MDL_MAX, ///<maximum value (double)
    MDL_MICROGRAPH, ///< Name of a micrograph (std::string)
    MDL_MICROGRAPH_TILTED, ///< Name of the corresponding tilted micrograph (std::string)
    MDL_MIN, ///<minimum value (double)
    MDL_MIRRORFRAC, ///< Mirror fraction for a Maximum Likelihood model
    MDL_MISSINGREGION_NR, ///< Number of missing region in subtomogram
    MDL_MISSINGREGION_TYPE, ///< Type of missing region in subtomogram
    MDL_MISSINGREGION_THY0, ///< Initial tilt angle in Y for missing region in subtomogram
    MDL_MISSINGREGION_THYF, ///< Final tilt angle in Y for missing region in subtomogram
    MDL_MISSINGREGION_THX0, ///< Initial tilt angle in X for missing region in subtomogram
    MDL_MISSINGREGION_THXF, ///< Final tilt angle in X for missing region in subtomogram
    MDL_MODELFRAC, ///< Model fraction (alpha_k) for a Maximum Likelihood model
    MDL_NEIGHBORS, ///< Vector of indexes to points some "neighbors"
    MDL_NEIGHBORHOOD_RADIUS, ///< Radius of the neigborhood (radians)
    MDL_NMA, ///< Normal mode displacements (vector double)
    MDL_NMA_MODEFILE, ///< File with an NMA mode
    MDL_NOISE_ANGLES, ///< Noise description for projected angles
    MDL_NOISE_PARTICLE_COORD, ///< Noise description for particle's center coordenates (when projecting)
    MDL_NOISE_COORD,  //Use instead of MDL_NOISE_PARTICLE_COORD in future
    MDL_NOISE_PIXEL_LEVEL, ///< Noise description for pixels' gray level (when projecting)
    MDL_ORDER, /// auxiliary label to be used as an index (long)
    MDL_ORIGINX, ///< Origin for the image in the X axis (double)
    MDL_ORIGINY, ///< Origin for the image in the Y axis (double)
    MDL_ORIGINZ, ///< Origin for the image in the Z axis (double)
    MDL_PICKING_COLOR, ///< Color for particle picking
    MDL_PICKING_FAMILY, ///< Family for particle picking
    MDL_PICKING_FAMILY_STATE, ///< Family state for particle picking	
    MDL_PICKING_MICROGRAPH_FAMILY_STATE, ///< Micrograph family state for particle picking
    MDL_PICKING_PARTICLE_SIZE, ///< Particle size for particle picking
    MDL_PMAX, ///< Maximum value of normalized probability function (now called "Pmax/sumP") (double)
    MDL_POINTSASYMETRICUNIT, //number of points in asymmetric unit
    MDL_PRJ_DIMENSIONS, // X,Y dimensions for the generated projections
    MDL_PRJ_TILT_RANGE, // Vector with the initial and final tilt angle values, and step size
    MDL_PRJ_VOL,        // Volume file name to generate projections from
    MDL_DIMENSIONS_3D,  // X,Y,Z dimensions
    MDL_DIMENSIONS_2D,  // X,Y dimensions
    MDL_PSD, ///< A Power Spectrum Density file name (std::string)
    MDL_PSD_ENHANCED, ///< A enhanced Power Spectrum Density file name (std::string)
    MDL_RANDOMSEED, ///< Seed for random number generator
    MDL_REF3D, ///< 3D Class to which the image belongs (int)
    MDL_REF, ///< Class to which the image belongs (int)
    MDL_REFMD, ///< Name of Metadata file for all references(string)
    MDL_RESOLUTION_DPR, ///<differential phase residual (double)
    MDL_RESOLUTION_ERRORL2, ///<Error in l2 (double)
    MDL_RESOLUTION_FRC, ///<Fourier shell correlation (double)
    MDL_RESOLUTION_FRCRANDOMNOISE, ///<Fourier shell correlation noise (double)
    MDL_RESOLUTION_FREQ, ///<Frequency in 1/A (double)
    MDL_RESOLUTION_FREQREAL, ///< Frequency in A (double)
    MDL_SAMPLINGRATE, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATEX, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATEY, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATEZ, ///< sampling rate in A/pixel (double)
    MDL_SCALE, ///< scaling factor for an image or volume (double)
    MDL_SELFILE, ///< Name of an image (std::string)
    MDL_SERIE, ///< A collection of micrographs, e.g. a tilt serie (std::string)
    MDL_SHIFTX, ///< Shift for the image in the X axis (double)
    MDL_SHIFTY, ///< Shift for the image in the Y axis (double)
    MDL_SHIFTZ, ///< Shift for the image in the Z axis (double)
    MDL_SHIFT_CRYSTALX, ///< Shift for the image in the X axis (double) for crystals
    MDL_SHIFT_CRYSTALY, ///< Shift for the image in the Y axis (double) for crystals
    MDL_SHIFT_CRYSTALZ, ///< Shift for the image in the Z axis (double) for crystals
    MDL_SIGMANOISE, ///< Standard deviation of the noise in ML model
    MDL_SIGMAOFFSET, ///< Standard deviation of the offsets in ML model
    MDL_SIGNALCHANGE, ///< Signal change for an image
    MDL_SPHERICALABERRATION, ///<ctf definition azimutal angle
    MDL_STDDEV, ///<stdandard deviation value (double)
    MDL_SUM, ///< Sum of elements of a given type (double) [this is a genereic type do not use to transfer information to another program]
    MDL_SUMWEIGHT, ///< Sum of all weights in ML model
    MDL_SYMNO, ///< Symmetry number for a projection (used in ART)
    MDL_TRANSFORMATIONMTRIX, ///< transformation matrix(vector double)
    MDL_VOLTAGE, ///< microscope voltage (double)
    MDL_WEIGHT, ///< Weight assigned to the image (double)
    MDL_WROBUST, ///< Weight of t-student distribution in robust Maximum likelihood
    MDL_X, ///< X component (double)
    MDL_XINT, ///< X component (int)
    MDL_XINTTILT, ///< X component in tilted micrograph (int)
    MDL_XSIZE, ///< X size (int)
    MDL_Y, ///< Y component (double)
    MDL_YINT, ///< Y component (int)
    MDL_YINTTILT, ///< Y component in tilted micrograph (int)
    MDL_YSIZE, ///< Y size (int)
    MDL_Z, ///< Z component (double)
    MDL_ZINT, ///< Z component (int)
    MDL_ZSCORE, ///< Z Score (double)
    MDL_ZSIZE, ///< Z size (int)
    MDL_PHANTOM_BGDENSITY, ///< Phantom background density (double)
    MDL_PHANTOM_SCALE, ///< Number which will multiply all features (double)
    MDL_PHANTOM_FEATURE_TYPE, ///< Type of the feature (Sphere, Blob, ...) (std::string)
    MDL_PHANTOM_FEATURE_OPERATION, ///< Operation in case of overlapping features (+,-)
    MDL_PHANTOM_FEATURE_DENSITY, ///< The density of the feature (double)
    MDL_PHANTOM_FEATURE_CENTER, ///< Center of the feature (vector double)
    MDL_PHANTOM_FEATURE_SPECIFIC, ///< Specific parameters for a feature (vector double)

    MDL_PRJ_ANGFILE,  ///< File for generated angles
    MDL_PRJ_ROT_RANGE,
    MDL_PRJ_ROT_Noise ,  /// < Rotational angle dev and mean noise (vector double)
    MDL_PRJ_ROT_RANDSTR,  /// < Type of randomness for Rotational (std::string)
    MDL_PRJ_TILT_Noise,  /// < Tilt angle dev and mean noise (vector double)
    MDL_PRJ_TILT_RANDSTR,  /// < Type of randomness for Tilt (std::string)
    MDL_PRJ_PSI_RANGE,  /// < Psi angle range (vector double)
    MDL_PRJ_PSI_Noise,  /// < Psi angle dev and mean noise (vector double)
    MDL_PRJ_PSI_RANDSTR, /// < Type of randomness for Psi (std::string)

    MDL_2D_LATTICE_VECA,   /// < Lattice vector for projection (vector double)
    MDL_2D_LATTICE_VECB,   /// < Lattice vector for projection (vector double)
    MDL_CRYSTAL_DISAPPEAR_THRE,   /// < Disappearing threshold (double)
    MDL_CRYSTAL_SHFILE,   /// < Shift file for crystal projection
    MDL_ORTHOGONAL_PROJECTION,   /// <Orthogonal projection or not (bool)
    MDL_CRYSTAL_PROJ,   /// < Have a crystal projection (boo)

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
    bool addLabel(MDLabel label);

    /** Clear elements of the row */
    void clear();
    /** Return number of labels present */
    int size() const;
    /** Function to test whether is empty */
    bool empty() const;

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

    bool getValue(MDObject &object) const;

    /** Set value */
    template <typename T>
    void setValue(MDLabel label, const T &d)
    {
      if (objects[label] == NULL)
      {
        objects[label] = new MDObject(label, d);
        order[_size] = label;
        ++_size;
      }
      else
        objects[label]->setValue(d);
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


    static void addLabel(MDLabel label, MDLabelType type, const String &name, int tags=TAGLABEL_NOTAG, const String &name2 = "", const String &name3 = "");

    friend class MDLabelStaticInit;
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

        MDL::addLabel(MDL_ANGLE_COMPARISON, LABEL_VECTOR_DOUBLE, "angleComparison", TAGLABEL_NOTAG, "psi2");
        MDL::addLabel(MDL_ANGLEPSI, LABEL_DOUBLE, "anglePsi", TAGLABEL_NOTAG, "psi");
        MDL::addLabel(MDL_ANGLEPSI2, LABEL_DOUBLE, "anglePsi2", TAGLABEL_NOTAG, "psi2");
        MDL::addLabel(MDL_ANGLEROT, LABEL_DOUBLE, "angleRot", TAGLABEL_NOTAG, "rot");
        MDL::addLabel(MDL_ANGLEROT2, LABEL_DOUBLE, "angleRot2", TAGLABEL_NOTAG, "rot2");
        MDL::addLabel(MDL_ANGLETILT, LABEL_DOUBLE, "angleTilt", TAGLABEL_NOTAG, "tilt");
        MDL::addLabel(MDL_ANGLETILT2, LABEL_DOUBLE, "angleTilt2", TAGLABEL_NOTAG, "tilt2");
        MDL::addLabel(MDL_ANGLE_Y, LABEL_DOUBLE, "angleY", TAGLABEL_NOTAG);
        MDL::addLabel(MDL_ANGLE_Y2, LABEL_DOUBLE, "angleY2", TAGLABEL_NOTAG);
        MDL::addLabel(MDL_ASSOCIATED_IMAGE1, LABEL_STRING, "associatedImage1", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_ASSOCIATED_IMAGE2, LABEL_STRING, "associatedImage2", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_ASSOCIATED_IMAGE3, LABEL_STRING, "associatedImage3", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_ASSOCIATED_IMAGE4, LABEL_STRING, "associatedImage4", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_ASSOCIATED_IMAGE5, LABEL_STRING, "associatedImage5", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_AVG, LABEL_DOUBLE, "avg");
        MDL::addLabel(MDL_AZIMUTALANGLE, LABEL_DOUBLE, "azimutalAngle");
        MDL::addLabel(MDL_BGMEAN, LABEL_DOUBLE, "bgMean");
        MDL::addLabel(MDL_BLOCK, LABEL_INT, "blockNumber");
        MDL::addLabel(MDL_CELLX, LABEL_INT, "cellX");
        MDL::addLabel(MDL_CELLY, LABEL_INT, "cellY");
        MDL::addLabel(MDL_CL2D_CHANGES, LABEL_INT, "cl2d_changes");
        MDL::addLabel(MDL_CL2D_SIMILARITY, LABEL_DOUBLE, "cl2dsimilarity");
        MDL::addLabel(MDL_CLASS_COUNT, LABEL_SIZET, "class_count");
        MDL::addLabel(MDL_CLASSIFICATION_DATA, LABEL_VECTOR_DOUBLE, "ClassificationData");
        MDL::addLabel(MDL_CLASSIFICATION_DATA_SIZE, LABEL_INT, "ClassificationDataSize");
        MDL::addLabel(MDL_CLASSIFICATION_DPR_05, LABEL_DOUBLE, "ClassificationDPR05");
        MDL::addLabel(MDL_CLASSIFICATION_INTRACLASS_DISTANCE, LABEL_DOUBLE, "ClassificationIntraclassDistance");
        MDL::addLabel(MDL_CLASSIFICATION_FRC_05, LABEL_DOUBLE, "ClassificationFRC05");
        MDL::addLabel(MDL_COMMENT, LABEL_STRING, "comment");
        MDL::addLabel(MDL_COST, LABEL_DOUBLE, "cost");
        MDL::addLabel(MDL_COUNT, LABEL_SIZET, "count");
        MDL::addLabel(MDL_COUNT2, LABEL_SIZET, "count2");
        MDL::addLabel(MDL_CTFINPUTPARAMS, LABEL_STRING, "CTFInputParams", TAGLABEL_TEXTFILE);
        MDL::addLabel(MDL_CTFMODEL, LABEL_STRING, "CTFModel", TAGLABEL_CTFPARAM);
        MDL::addLabel(MDL_CTFMODEL2, LABEL_STRING, "CTFModel2", TAGLABEL_CTFPARAM);
        MDL::addLabel(MDL_CTF_SAMPLING_RATE, LABEL_DOUBLE, "CTF_Sampling_rate");
        MDL::addLabel(MDL_CTF_SAMPLING_RATE_Z, LABEL_DOUBLE, "CTF_Sampling_rate_z");
        MDL::addLabel(MDL_CTF_VOLTAGE, LABEL_DOUBLE, "CTF_Voltage");
        MDL::addLabel(MDL_CTF_DEFOCUSA, LABEL_DOUBLE, "CTF_Defocus_A");//average defocus
        MDL::addLabel(MDL_CTF_DEFOCUSU, LABEL_DOUBLE, "CTF_Defocus_U");
        MDL::addLabel(MDL_CTF_DEFOCUSV, LABEL_DOUBLE, "CTF_Defocus_V");
        MDL::addLabel(MDL_CTF_DEFOCUS_ANGLE, LABEL_DOUBLE, "CTF_Defocus_angle");
        MDL::addLabel(MDL_CTF_CS, LABEL_DOUBLE, "CTF_Spherical_aberration");
        MDL::addLabel(MDL_CTF_CA, LABEL_DOUBLE, "CTF_Chromatic_aberration");
        MDL::addLabel(MDL_CTF_GROUP, LABEL_INT, "CTFGroup");
        MDL::addLabel(MDL_CTF_ENERGY_LOSS, LABEL_DOUBLE, "CTF_Energy_loss");
        MDL::addLabel(MDL_CTF_LENS_STABILITY, LABEL_DOUBLE, "CTF_Lens_stability");
        MDL::addLabel(MDL_CTF_CONVERGENCE_CONE, LABEL_DOUBLE, "CTF_Convergence_cone");
        MDL::addLabel(MDL_CTF_LONGITUDINAL_DISPLACEMENT, LABEL_DOUBLE, "CTF_Longitudinal_displacement");
        MDL::addLabel(MDL_CTF_TRANSVERSAL_DISPLACEMENT, LABEL_DOUBLE, "CTF_Transversal_displacement");
        MDL::addLabel(MDL_CTF_Q0, LABEL_DOUBLE, "CTF_Q0");
        MDL::addLabel(MDL_CTF_K, LABEL_DOUBLE, "CTF_K");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN_K, LABEL_DOUBLE, "CTFBG_Gaussian_K");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN_SIGMAU, LABEL_DOUBLE, "CTFBG_Gaussian_SigmaU");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN_SIGMAV, LABEL_DOUBLE, "CTFBG_Gaussian_SigmaV");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN_CU, LABEL_DOUBLE, "CTFBG_Gaussian_CU");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN_CV, LABEL_DOUBLE, "CTFBG_Gaussian_CV");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN_ANGLE, LABEL_DOUBLE, "CTFBG_Gaussian_Angle");
        MDL::addLabel(MDL_CTFBG_SQRT_K, LABEL_DOUBLE, "CTFBG_Sqrt_K");
        MDL::addLabel(MDL_CTFBG_SQRT_U, LABEL_DOUBLE, "CTFBG_Sqrt_U");
        MDL::addLabel(MDL_CTFBG_SQRT_V, LABEL_DOUBLE, "CTFBG_Sqrt_V");
        MDL::addLabel(MDL_CTFBG_SQRT_ANGLE, LABEL_DOUBLE, "CTFBG_Sqrt_Angle");
        MDL::addLabel(MDL_CTFBG_BASELINE, LABEL_DOUBLE, "CTFBG_Baseline");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN2_K, LABEL_DOUBLE, "CTFBG_Gaussian2_K");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN2_SIGMAU, LABEL_DOUBLE, "CTFBG_Gaussian2_SigmaU");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN2_SIGMAV, LABEL_DOUBLE, "CTFBG_Gaussian2_SigmaV");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN2_CU, LABEL_DOUBLE, "CTFBG_Gaussian2_CU");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN2_CV, LABEL_DOUBLE, "CTFBG_Gaussian2_CV");
        MDL::addLabel(MDL_CTFBG_GAUSSIAN2_ANGLE, LABEL_DOUBLE, "CTFBG_Gaussian2_Angle");
        MDL::addLabel(MDL_CTF_CRITERION_PSDCORRELATION90, LABEL_DOUBLE, "CTFCrit_psdcorr90");
        MDL::addLabel(MDL_CTF_CRITERION_FIRSTZERORATIO, LABEL_DOUBLE, "CTFCrit_firstZeroRatio");
        MDL::addLabel(MDL_CTF_CRITERION_FIRSTZEROAVG, LABEL_DOUBLE, "CTFCrit_firstZero");
        MDL::addLabel(MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT, LABEL_DOUBLE, "CTFCrit_disagree");
        MDL::addLabel(MDL_CTF_CRITERION_DAMPING, LABEL_DOUBLE, "CTFCrit_damping");
        MDL::addLabel(MDL_CTF_CRITERION_PSDRADIALINTEGRAL, LABEL_DOUBLE, "CTFCrit_psdint");
        MDL::addLabel(MDL_CTF_CRITERION_FITTINGSCORE, LABEL_DOUBLE, "CTFCrit_Fitting");
        MDL::addLabel(MDL_CTF_CRITERION_FITTINGCORR13, LABEL_DOUBLE, "CTFCrit_Corr13");
        MDL::addLabel(MDL_CTF_CRITERION_PSDVARIANCE, LABEL_DOUBLE, "CTFCrit_PSDStdQ");
        MDL::addLabel(MDL_CTF_CRITERION_PSDPCA1VARIANCE, LABEL_DOUBLE, "CTFCrit_PSDPCA1");
        MDL::addLabel(MDL_CTF_CRITERION_PSDPCARUNSTEST, LABEL_DOUBLE, "CTFCrit_PSDPCARuns");
        MDL::addLabel(MDL_CTF_CRITERION_NORMALITY, LABEL_DOUBLE, "CTFCrit_Normality");
        MDL::addLabel(MDL_CTF_DOWNSAMPLE_PERFORMED, LABEL_DOUBLE, "CTFDownsampleFactor");

        MDL::addLabel(MDL_CTF_XRAY_DIMENSIONS, LABEL_VECTOR_DOUBLE, "CTF_Xray_dimensions");
        MDL::addLabel(MDL_CTF_XRAY_LAMBDA, LABEL_DOUBLE, "CTF_Xray_lambda");
        MDL::addLabel(MDL_CTF_XRAY_LENS_TYPE, LABEL_STRING, "CTF_Xray_lens_type");
        MDL::addLabel(MDL_CTF_XRAY_OUTER_ZONE_WIDTH, LABEL_DOUBLE, "CTF_Xray_OuterZoneWidth");
        MDL::addLabel(MDL_CTF_XRAY_ZONES_NUMBER, LABEL_DOUBLE, "CTF_Xray_ZonesN");

        MDL::addLabel(MDL_DATATYPE, LABEL_INT, "datatype");
        MDL::addLabel(MDL_DEFGROUP, LABEL_INT, "defocusGroup");
        MDL::addLabel(MDL_DIRECTION, LABEL_VECTOR_DOUBLE, "direction");

        MDL::addLabel(MDL_DM3_IDTAG, LABEL_INT, "IdTag");
        MDL::addLabel(MDL_DM3_NODEID, LABEL_INT, "NodeID");
        MDL::addLabel(MDL_DM3_NUMBER_TYPE, LABEL_INT, "Number_Type");
        MDL::addLabel(MDL_DM3_PARENTID, LABEL_INT, "ParentID");
        MDL::addLabel(MDL_DM3_TAGCLASS, LABEL_STRING, "Tag_Class");
        MDL::addLabel(MDL_DM3_TAGNAME, LABEL_STRING, "TagName");
        MDL::addLabel(MDL_DM3_SIZE, LABEL_INT, "Size");
        MDL::addLabel(MDL_DM3_VALUE, LABEL_VECTOR_DOUBLE, "Value");

        MDL::addLabel(MDL_ENABLED, LABEL_INT, "enabled");
        MDL::addLabel(MDL_FLIP, LABEL_BOOL, "flip", TAGLABEL_NOTAG, "Flip");
        MDL::addLabel(MDL_FOM, LABEL_DOUBLE, "fom");
        MDL::addLabel(MDL_IDX, LABEL_SIZET, "index");
        MDL::addLabel(MDL_IMAGE, LABEL_STRING, "image", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_ORIGINAL, LABEL_STRING, "original_image", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_TILTED, LABEL_STRING, "tilted_image", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMGMD, LABEL_STRING, "imageMetaData", TAGLABEL_METADATA);
        MDL::addLabel(MDL_INTSCALE, LABEL_DOUBLE, "intScale");
        MDL::addLabel(MDL_ITER, LABEL_INT, "iterationNumber");
        MDL::addLabel(MDL_K, LABEL_DOUBLE, "K");
        MDL::addLabel(MDL_KERDENSOM_FUNCTIONAL, LABEL_DOUBLE, "KerDenSOM_Functional");
        MDL::addLabel(MDL_KERDENSOM_REGULARIZATION, LABEL_DOUBLE, "KerDenSOM_Regularization");
        MDL::addLabel(MDL_KERDENSOM_SIGMA, LABEL_DOUBLE, "KerDenSOM_Sigma");
        MDL::addLabel(MDL_KEYWORDS, LABEL_STRING, "Keywords");
        MDL::addLabel(MDL_KSTEST, LABEL_DOUBLE, "KStest");
        MDL::addLabel(MDL_LL, LABEL_DOUBLE, "logLikelihood", TAGLABEL_NOTAG, "LL");
        MDL::addLabel(MDL_MAGNIFICATION, LABEL_DOUBLE, "Magnification");
        MDL::addLabel(MDL_MAPTOPOLOGY, LABEL_STRING, "mapTopology");
        MDL::addLabel(MDL_MASK, LABEL_STRING, "mask", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_MAXCC, LABEL_DOUBLE, "maxCC");
        MDL::addLabel(MDL_MAX, LABEL_DOUBLE, "max");
        MDL::addLabel(MDL_MICROGRAPH, LABEL_STRING, "micrograph", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_TILTED, LABEL_STRING, "micrographTilted", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MIN, LABEL_DOUBLE, "min");
        MDL::addLabel(MDL_MIRRORFRAC, LABEL_DOUBLE, "mirrorFraction");
        MDL::addLabel(MDL_MISSINGREGION_NR, LABEL_INT, "missingRegionNumber", TAGLABEL_NOTAG, "Wedge");
        MDL::addLabel(MDL_MISSINGREGION_TYPE, LABEL_STRING, "missingRegionType");
        MDL::addLabel(MDL_MISSINGREGION_THX0, LABEL_DOUBLE, "missingRegionThetaX0");
        MDL::addLabel(MDL_MISSINGREGION_THXF, LABEL_DOUBLE, "missingRegionThetaXF");
        MDL::addLabel(MDL_MISSINGREGION_THY0, LABEL_DOUBLE, "missingRegionThetaY0");
        MDL::addLabel(MDL_MISSINGREGION_THYF, LABEL_DOUBLE, "missingRegionThetaYF");
        MDL::addLabel(MDL_MODELFRAC, LABEL_DOUBLE, "modelFraction");
        MDL::addLabel(MDL_NEIGHBORS, LABEL_VECTOR_SIZET, "neighbors");
        MDL::addLabel(MDL_NEIGHBORHOOD_RADIUS, LABEL_DOUBLE, "neighborhoodRadius");
        MDL::addLabel(MDL_NMA, LABEL_VECTOR_DOUBLE, "NMADisplacements");
        MDL::addLabel(MDL_NMA_MODEFILE, LABEL_STRING, "NMAModefile", TAGLABEL_TEXTFILE);
		MDL::addLabel(MDL_NOISE_ANGLES, LABEL_VECTOR_DOUBLE, "noiseAngles");
		MDL::addLabel(MDL_NOISE_COORD, LABEL_VECTOR_DOUBLE, "noiseCoord");
		MDL::addLabel(MDL_NOISE_PARTICLE_COORD, LABEL_VECTOR_DOUBLE, "noiseParticleCoord");
		MDL::addLabel(MDL_NOISE_PIXEL_LEVEL, LABEL_VECTOR_DOUBLE, "noisePixelLevel");
		MDL::addLabel(MDL_ORDER, LABEL_SIZET, "order_");
        MDL::addLabel(MDL_ORIGINX, LABEL_DOUBLE, "originX");
        MDL::addLabel(MDL_ORIGINY, LABEL_DOUBLE, "originY");
        MDL::addLabel(MDL_ORIGINZ, LABEL_DOUBLE, "originZ");
        MDL::addLabel(MDL_PICKING_FAMILY, LABEL_STRING, "family");
        MDL::addLabel(MDL_PICKING_COLOR, LABEL_INT, "color");
        MDL::addLabel(MDL_PICKING_PARTICLE_SIZE, LABEL_INT, "particleSize");
        MDL::addLabel(MDL_PICKING_FAMILY_STATE, LABEL_STRING, "family_state");
        MDL::addLabel(MDL_PICKING_MICROGRAPH_FAMILY_STATE, LABEL_STRING, "micrograph_family_state");
        MDL::addLabel(MDL_PMAX, LABEL_DOUBLE, "pMax", TAGLABEL_NOTAG, "Pmax", "sumP");
        MDL::addLabel(MDL_POINTSASYMETRICUNIT, LABEL_SIZET, "pointsasymmetricunit");
        MDL::addLabel(MDL_PRJ_DIMENSIONS, LABEL_VECTOR_DOUBLE, "projDimensions");
        MDL::addLabel(MDL_PRJ_TILT_RANGE, LABEL_VECTOR_DOUBLE, "projTiltRange");
        MDL::addLabel(MDL_PRJ_VOL, LABEL_STRING, "projVolume", TAGLABEL_VOLUME);
        MDL::addLabel(MDL_DIMENSIONS_3D, LABEL_VECTOR_DOUBLE, "dimensions3D");
        MDL::addLabel(MDL_DIMENSIONS_2D, LABEL_VECTOR_DOUBLE, "dimensions2D");
        MDL::addLabel(MDL_PSD, LABEL_STRING, "powerSpectrum", TAGLABEL_PSD);
        MDL::addLabel(MDL_PSD_ENHANCED, LABEL_STRING, "enhancedPowerSpectrum", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_RANDOMSEED, LABEL_INT, "randomSeed");
        MDL::addLabel(MDL_REF3D, LABEL_INT, "ref3d");
        MDL::addLabel(MDL_REF, LABEL_INT, "ref", TAGLABEL_NOTAG, "Ref");
        MDL::addLabel(MDL_REFMD, LABEL_STRING, "referenceMetaData", TAGLABEL_METADATA);
        MDL::addLabel(MDL_RESOLUTION_DPR, LABEL_DOUBLE, "DPR");
        MDL::addLabel(MDL_RESOLUTION_ERRORL2, LABEL_DOUBLE, "Error_l2");
        MDL::addLabel(MDL_RESOLUTION_FRC, LABEL_DOUBLE, "FRC");
        MDL::addLabel(MDL_RESOLUTION_FRCRANDOMNOISE, LABEL_DOUBLE, "FRC_random_noise");
        MDL::addLabel(MDL_RESOLUTION_FREQ, LABEL_DOUBLE, "Resol_Inverse_Ang");
        MDL::addLabel(MDL_RESOLUTION_FREQREAL, LABEL_DOUBLE, "Resol_Ang");
        MDL::addLabel(MDL_SAMPLINGRATE, LABEL_DOUBLE, "sampling_rate");
        MDL::addLabel(MDL_SAMPLINGRATEX, LABEL_DOUBLE, "sampling_rateX");
        MDL::addLabel(MDL_SAMPLINGRATEY, LABEL_DOUBLE, "sampling_rateY");
        MDL::addLabel(MDL_SAMPLINGRATEZ, LABEL_DOUBLE, "sampling_rateZ");
        MDL::addLabel(MDL_SCALE, LABEL_DOUBLE, "scale", TAGLABEL_NOTAG, "Scale");
        MDL::addLabel(MDL_SELFILE, LABEL_STRING, "selfile", TAGLABEL_METADATA);
        MDL::addLabel(MDL_SERIE, LABEL_STRING, "serie");
        MDL::addLabel(MDL_SHIFTX, LABEL_DOUBLE, "shiftX", TAGLABEL_NOTAG, "Xoff");
        MDL::addLabel(MDL_SHIFTY, LABEL_DOUBLE, "shiftY", TAGLABEL_NOTAG, "Yoff");
        MDL::addLabel(MDL_SHIFTZ, LABEL_DOUBLE, "shiftZ", TAGLABEL_NOTAG, "Zoff");
        MDL::addLabel(MDL_SHIFT_CRYSTALX, LABEL_DOUBLE, "crystalShiftX", TAGLABEL_NOTAG, "Xoff");
        MDL::addLabel(MDL_SHIFT_CRYSTALY, LABEL_DOUBLE, "crystalShiftY", TAGLABEL_NOTAG, "Yoff");
        MDL::addLabel(MDL_SHIFT_CRYSTALZ, LABEL_DOUBLE, "crystalShiftZ", TAGLABEL_NOTAG, "Zoff");
        MDL::addLabel(MDL_SIGMANOISE, LABEL_DOUBLE, "sigmaNoise");
        MDL::addLabel(MDL_SIGMAOFFSET, LABEL_DOUBLE, "sigmaOffset");
        MDL::addLabel(MDL_SIGNALCHANGE, LABEL_DOUBLE, "signalChange");
        MDL::addLabel(MDL_SPHERICALABERRATION, LABEL_DOUBLE, "sphericalAberration");
        MDL::addLabel(MDL_STDDEV, LABEL_DOUBLE, "stddev");
        MDL::addLabel(MDL_SUM, LABEL_DOUBLE, "sum");
        MDL::addLabel(MDL_SUMWEIGHT, LABEL_DOUBLE, "sumWeight");
        MDL::addLabel(MDL_SYMNO, LABEL_INT, "symNo");
        MDL::addLabel(MDL_TRANSFORMATIONMTRIX, LABEL_VECTOR_DOUBLE, "transMat");
        MDL::addLabel(MDL_VOLTAGE, LABEL_DOUBLE, "voltage");
        MDL::addLabel(MDL_WEIGHT, LABEL_DOUBLE, "weight", TAGLABEL_NOTAG, "Weight");
        MDL::addLabel(MDL_WROBUST, LABEL_DOUBLE, "wRobust");
        MDL::addLabel(MDL_X, LABEL_DOUBLE, "X");
        MDL::addLabel(MDL_XINT, LABEL_INT, "Xcoor", TAGLABEL_NOTAG, "<X position>");
        MDL::addLabel(MDL_XINTTILT, LABEL_INT, "XcoorTilt");
        MDL::addLabel(MDL_XSIZE, LABEL_INT, "Xsize");
        MDL::addLabel(MDL_Y, LABEL_DOUBLE, "Y");
        MDL::addLabel(MDL_YINT, LABEL_INT, "Ycoor", TAGLABEL_NOTAG, "<Y position>");
        MDL::addLabel(MDL_YINTTILT, LABEL_INT, "YcoorTilt");
        MDL::addLabel(MDL_YSIZE, LABEL_INT, "Ysize");
        MDL::addLabel(MDL_Z, LABEL_DOUBLE, "Z");
        MDL::addLabel(MDL_ZINT, LABEL_INT, "Zcoor");
        MDL::addLabel(MDL_ZSCORE, LABEL_DOUBLE, "Zscore");
        MDL::addLabel(MDL_ZSIZE, LABEL_INT, "Zsize");
        MDL::addLabel(MDL_PHANTOM_BGDENSITY, LABEL_DOUBLE, "phantomBGDensity");
        MDL::addLabel(MDL_PHANTOM_FEATURE_TYPE, LABEL_STRING, "featureType");
        MDL::addLabel(MDL_PHANTOM_FEATURE_OPERATION, LABEL_STRING, "featureOperation");
        MDL::addLabel(MDL_PHANTOM_FEATURE_DENSITY, LABEL_DOUBLE, "featureDensity");
        MDL::addLabel(MDL_PHANTOM_FEATURE_CENTER, LABEL_VECTOR_DOUBLE, "featureCenter");
        MDL::addLabel(MDL_PHANTOM_FEATURE_SPECIFIC, LABEL_VECTOR_DOUBLE, "featureSpecificVector");
        MDL::addLabel(MDL_PRJ_ANGFILE, LABEL_STRING, "angleFile");
        MDL::addLabel(MDL_PRJ_ROT_RANGE, LABEL_VECTOR_DOUBLE, "projRotRange");
        MDL::addLabel(MDL_PRJ_ROT_RANDSTR, LABEL_STRING, "projRotRandomness");
        MDL::addLabel(MDL_PRJ_ROT_Noise, LABEL_VECTOR_DOUBLE, "projRotNoise");
        MDL::addLabel(MDL_PRJ_TILT_RANDSTR, LABEL_STRING, "projTiltRandomness");
        MDL::addLabel(MDL_PRJ_TILT_Noise, LABEL_VECTOR_DOUBLE, "projTiltNoise");
        MDL::addLabel(MDL_PRJ_PSI_RANGE, LABEL_VECTOR_DOUBLE, "projPsiRange");
        MDL::addLabel(MDL_PRJ_PSI_RANDSTR, LABEL_STRING, "projPsiRandomness");
        MDL::addLabel(MDL_PRJ_PSI_Noise, LABEL_VECTOR_DOUBLE, "projPsiNoise");
        MDL::addLabel(MDL_2D_LATTICE_VECA, LABEL_VECTOR_DOUBLE, "latticeVec1");
        MDL::addLabel(MDL_2D_LATTICE_VECB, LABEL_VECTOR_DOUBLE, "latticeVec2");
        MDL::addLabel(MDL_CRYSTAL_DISAPPEAR_THRE, LABEL_DOUBLE, "crystalDisThresh");
        MDL::addLabel(MDL_ORTHOGONAL_PROJECTION, LABEL_BOOL, "orthogonalProj");
        MDL::addLabel( MDL_CRYSTAL_SHFILE, LABEL_STRING, "crystalShiftFile");
        MDL::addLabel( MDL_CRYSTAL_PROJ, LABEL_BOOL, "crystalProj");

        //Create an static empty header for image initialization
        MDL::emptyHeader.setValue(MDL_ORIGINX,  0.);
        MDL::emptyHeader.setValue(MDL_ORIGINY,  0.);
        MDL::emptyHeader.setValue(MDL_ORIGINZ,  0.);
        MDL::emptyHeader.setValue(MDL_SHIFTX,   0.);
        MDL::emptyHeader.setValue(MDL_SHIFTY,   0.);
        MDL::emptyHeader.setValue(MDL_SHIFTZ,   0.);
        MDL::emptyHeader.setValue(MDL_ANGLEROT, 0.);
        MDL::emptyHeader.setValue(MDL_ANGLETILT,0.);
        MDL::emptyHeader.setValue(MDL_ANGLEPSI, 0.);
        MDL::emptyHeader.setValue(MDL_WEIGHT,   1.);
        MDL::emptyHeader.setValue(MDL_FLIP,     false);
        MDL::emptyHeader.setValue(MDL_SCALE,    1.);

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
