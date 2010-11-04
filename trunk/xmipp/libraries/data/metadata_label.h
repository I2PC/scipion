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
#include "strings.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "funcs.h"

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

    MDL_ANGLEPSI2, ///< Psi angle of an image (double)
    MDL_ANGLEPSI, ///< Psi angle of an image (double)
    MDL_ANGLEROT2, ///< Rotation angle of an image (double)
    MDL_ANGLEROT, ///< Rotation angle of an image (double)
    MDL_ANGLETILT2, ///< Tilting angle of an image (double)
    MDL_ANGLETILT, ///< Tilting angle of an image (double)
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
    MDL_COMMENT, ///< A comment for this object /*** NOTE THIS IS A SPECIAL CASE AND SO IS TREATED ***/
    MDL_COST, ///< Cost for the image (double)
    MDL_COUNT, ///< Number of elements of a type (int) [this is a genereic type do not use to transfer information to another program]
    MDL_CTFINPUTPARAMS, ///< Parameters file for the CTF Model (std::string)
    MDL_CTFMODEL, ///< Name for the CTF Model (std::string)
    MDL_CTFMODEL2, ///< Name for another CTF model (std::string)
    MDL_CTF_SAMPLING_RATE, ///< Sampling rate
    MDL_CTF_VOLTAGE, ///< Microscope voltage (kV)
    MDL_CTF_DEFOCUSU, ///< Defocus U (Angstroms)
    MDL_CTF_DEFOCUSV, ///< Defocus V (Angstroms)
    MDL_CTF_DEFOCUS_ANGLE, ///< Defocus angle (degrees)
    MDL_CTF_CS, ///< Spherical aberration
    MDL_CTF_CA, ///< Chromatic aberration
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
    MDL_CTF_CRITERION_PSDPCA1VARIACE, ///< Variance in the first principal component of the PSDs
    MDL_DATATYPE, ///< if read from file original image datatype, this is an struct defined in image
    MDL_DEFGROUP, ///< Defocus group
    MDL_DEFOCUSU, ///< microscope defocus U direction (double)
    MDL_DEFOCUSV, ///< microscope defocus V direction (double)

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
    MDL_IMAGE_CLASS_COUNT, ///< Number of images assigned to the same class as this image
    MDL_IMAGE_CLASS_GROUP, ///< Name of the class group for this image (metadata with all the images assigned to that class)
    MDL_IMAGE_CLASS, ///< Name of the class representative for this image
    MDL_IMAGE, ///< Name of an image (std::string)
    MDL_IMAGE_ORIGINAL, ///< Name of an image from which MDL_IMAGE is coming from
    MDL_IMGMD, ///< Name of Metadata file for all images (string)
    MDL_INTSCALE, ///< Intensity scale for an image
    MDL_ITER, ///< Current iteration number (int)
    MDL_K, ///< //ctf definition K
    MDL_KSTEST, ///<KS-test statistics
    MDL_LL, ///< contribution of an image to log-likelihood value
    MDL_MASK, ///< Name of a mask associated to image
    MDL_MAXCC, ///< Cross-correlation for the image (double)
    MDL_MAX, ///<maximum value (double)
    MDL_MICROGRAPH, ///< Name of a micrograph (std::string)
    MDL_MIN, ///<minimum value (double)
    MDL_MIRRORFRAC, ///< Mirror fraction for a Maximum Likelihood model
    MDL_MISSINGREGION_NR, ///< Number of missing region in subtomogram
    MDL_MISSINGREGION_TYPE, ///< Type of missing region in subtomogram
    MDL_MISSINGREGION_THY0, ///< Initial tilt angle in Y for missing region in subtomogram
    MDL_MISSINGREGION_THYF, ///< Final tilt angle in Y for missing region in subtomogram
    MDL_MISSINGREGION_THX0, ///< Initial tilt angle in X for missing region in subtomogram
    MDL_MISSINGREGION_THXF, ///< Final tilt angle in X for missing region in subtomogram
    MDL_MODELFRAC, ///< Model fraction (alpha_k) for a Maximum Likelihood model
    MDL_NMA, ///< Normal mode displacements (vector double)

    MDL_ORIGINX, ///< Origin for the image in the X axis (double)
    MDL_ORIGINY, ///< Origin for the image in the Y axis (double)
    MDL_ORIGINZ, ///< Origin for the image in the Z axis (double)
    MDL_PMAX, ///< Maximum value of normalized probability function (now called "Pmax/sumP") (double)
    MDL_PSD, ///< A Power Spectrum Density file name (std::string)
    MDL_Q0, ///< ctf definition Q0
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
    MDL_XINT, ///< X component (int)
    MDL_X, ///< X component (double)
    MDL_YINT, ///< Y component (int)
    MDL_Y, ///< Y component (double)
    MDL_ZINT, ///< Z component (int)
    MDL_Z, ///< Z component (double)
    MDL_ZSCORE, ///< Z Score (double)
    

    MDL_LAST_LABEL                       // **** NOTE ****: Do keep this label always at the end
    // it is here for looping purposes
};//close enum Label

/** Possible types of the values of labels */
enum MDLabelType
{
    LABEL_NOTYPE = -1, LABEL_INT, LABEL_BOOL, LABEL_DOUBLE, LABEL_FLOAT, LABEL_STRING, LABEL_VECTOR, LABEL_LONG
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
    static MDLabel str2Label(const String &labelName);
    /** Converts MDLabel to string */
    static String label2Str(const MDLabel label);
    /** Converts MDLabel to string representing SQL column*/
    static String label2SqlColumn(const MDLabel label);
    /** @} */

    /** @name Type checks
     * Several label type checks
     * @{
     */
    static bool isInt(const MDLabel label);
    static bool isBool(const MDLabel label);
    static bool isString(const MDLabel label);
    static bool isDouble(const MDLabel label);
    static bool isVector(const MDLabel label);
    static bool isValidLabel(const MDLabel label);
    static bool isValidLabel(const String &labelName);
    static MDLabelType labelType(const MDLabel label);
    static MDLabelType labelType(const String &labelName);
    /** @} */

private:
    static std::map<MDLabel, MDLabelData> data;
    static std::map<std::string, MDLabel> names;
    static MDLabelStaticInit initialization; //Just for initialization

    static void addLabel(MDLabel label, MDLabelType type, const String &name, const String &name2 = "", const String &name3 = "");

    friend class MDLabelStaticInit;
}
;//close class MLD definition
/** @} */


/** Just to work as static constructor for initialize labels data.
 */
class MDLabelStaticInit
{
private:
    MDLabelStaticInit()
    {
        ///==== Add labels entries from here in the SAME ORDER as declared in ENUM ==========
        //The label MDL_OBJID is special and should not be used
        MDL::addLabel(MDL_OBJID, LABEL_LONG, "objId");

        MDL::addLabel(MDL_ANGLEPSI2, LABEL_DOUBLE, "anglePsi2", "psi2");
        MDL::addLabel(MDL_ANGLEPSI, LABEL_DOUBLE, "anglePsi", "psi");
        MDL::addLabel(MDL_ANGLEROT2, LABEL_DOUBLE, "angleRot2", "rot2");
        MDL::addLabel(MDL_ANGLEROT, LABEL_DOUBLE, "angleRot", "rot");
        MDL::addLabel(MDL_ANGLETILT2, LABEL_DOUBLE, "angleTilt2", "tilt2");
        MDL::addLabel(MDL_ANGLETILT, LABEL_DOUBLE, "angleTilt", "tilt");
        MDL::addLabel(MDL_ASSOCIATED_IMAGE1, LABEL_STRING, "associatedImage1");
        MDL::addLabel(MDL_ASSOCIATED_IMAGE2, LABEL_STRING, "associatedImage2");
        MDL::addLabel(MDL_ASSOCIATED_IMAGE3, LABEL_STRING, "associatedImage3");
        MDL::addLabel(MDL_ASSOCIATED_IMAGE4, LABEL_STRING, "associatedImage4");
        MDL::addLabel(MDL_ASSOCIATED_IMAGE5, LABEL_STRING, "associatedImage5");
        MDL::addLabel(MDL_AVG, LABEL_DOUBLE, "avg");
        MDL::addLabel(MDL_AZIMUTALANGLE, LABEL_DOUBLE, "azimutalAngle");
        MDL::addLabel(MDL_BGMEAN, LABEL_DOUBLE, "bgMean");
        MDL::addLabel(MDL_BLOCK, LABEL_INT, "blockNumber");
        MDL::addLabel(MDL_CELLX, LABEL_INT, "cellX");
        MDL::addLabel(MDL_CELLY, LABEL_INT, "cellY");
        MDL::addLabel(MDL_COMMENT, LABEL_STRING, "comment");
        MDL::addLabel(MDL_COST, LABEL_DOUBLE, "cost");
        MDL::addLabel(MDL_COUNT, LABEL_INT, "count");
        MDL::addLabel(MDL_CTFINPUTPARAMS, LABEL_STRING, "CTFInputParams");
        MDL::addLabel(MDL_CTFMODEL, LABEL_STRING, "CTFModel");
        MDL::addLabel(MDL_CTFMODEL2, LABEL_STRING, "CTFModel2");
        MDL::addLabel(MDL_CTF_SAMPLING_RATE, LABEL_DOUBLE, "CTF_Sampling_rate");
        MDL::addLabel(MDL_CTF_VOLTAGE, LABEL_DOUBLE, "CTF_Voltage");
        MDL::addLabel(MDL_CTF_DEFOCUSU, LABEL_DOUBLE, "CTF_Defocus_U");
        MDL::addLabel(MDL_CTF_DEFOCUSV, LABEL_DOUBLE, "CTF_Defocus_V");
        MDL::addLabel(MDL_CTF_DEFOCUS_ANGLE, LABEL_DOUBLE, "CTF_Defocus_angle");
        MDL::addLabel(MDL_CTF_CS, LABEL_DOUBLE, "CTF_Spherical_aberration");
        MDL::addLabel(MDL_CTF_CA, LABEL_DOUBLE, "CTF_Chromatic_aberration");
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
        MDL::addLabel(MDL_CTF_CRITERION_PSDPCA1VARIACE, LABEL_DOUBLE, "CTFCrit_PSDPCA1");
        MDL::addLabel(MDL_DEFGROUP, LABEL_INT, "defocusGroup");
        MDL::addLabel(MDL_DATATYPE, LABEL_INT, "datatype");
        MDL::addLabel(MDL_DEFOCUSU, LABEL_DOUBLE, "defocusU");
        MDL::addLabel(MDL_DEFOCUSV, LABEL_DOUBLE, "defocusV");

        MDL::addLabel(MDL_DM3_IDTAG, LABEL_INT, "IdTag");
        MDL::addLabel(MDL_DM3_NODEID, LABEL_INT, "NodeID");
        MDL::addLabel(MDL_DM3_NUMBER_TYPE, LABEL_INT, "Number_Type");
        MDL::addLabel(MDL_DM3_PARENTID, LABEL_INT, "ParentID");
        MDL::addLabel(MDL_DM3_TAGCLASS, LABEL_STRING, "Tag_Class");
        MDL::addLabel(MDL_DM3_TAGNAME, LABEL_STRING, "TagName");
        MDL::addLabel(MDL_DM3_SIZE, LABEL_INT, "Size");
        MDL::addLabel(MDL_DM3_VALUE, LABEL_VECTOR, "Value");

        MDL::addLabel(MDL_ENABLED, LABEL_INT, "enabled");
        MDL::addLabel(MDL_FLIP, LABEL_BOOL, "flip", "Flip");
        MDL::addLabel(MDL_IMAGE_CLASS_COUNT, LABEL_INT, "class_count");
        MDL::addLabel(MDL_IMAGE_CLASS_GROUP, LABEL_STRING, "class_group");
        MDL::addLabel(MDL_IMAGE_CLASS, LABEL_STRING, "class_representative");
        MDL::addLabel(MDL_IMAGE, LABEL_STRING, "image");
        MDL::addLabel(MDL_IMAGE_ORIGINAL, LABEL_STRING, "original_image");
        MDL::addLabel(MDL_IMGMD, LABEL_STRING, "imageMetaData");
        MDL::addLabel(MDL_INTSCALE, LABEL_DOUBLE, "intScale");
        MDL::addLabel(MDL_ITER, LABEL_INT, "iterationNumber");
        MDL::addLabel(MDL_K, LABEL_DOUBLE, "K");
        MDL::addLabel(MDL_KSTEST, LABEL_DOUBLE, "KStest");
        MDL::addLabel(MDL_LL, LABEL_DOUBLE, "logLikelihood", "LL");
        MDL::addLabel(MDL_MASK, LABEL_STRING, "mask");
        MDL::addLabel(MDL_MAXCC, LABEL_DOUBLE, "maxCC");
        MDL::addLabel(MDL_MAX, LABEL_DOUBLE, "max");
        MDL::addLabel(MDL_MICROGRAPH, LABEL_STRING, "micrograph");
        MDL::addLabel(MDL_MIN, LABEL_DOUBLE, "min");
        MDL::addLabel(MDL_MIRRORFRAC, LABEL_DOUBLE, "mirrorFraction");
        MDL::addLabel(MDL_MISSINGREGION_NR, LABEL_INT, "missingRegionNumber");
        MDL::addLabel(MDL_MISSINGREGION_TYPE, LABEL_STRING, "missingRegionType");
        MDL::addLabel(MDL_MISSINGREGION_THX0, LABEL_DOUBLE, "missingRegionThetaX0");
        MDL::addLabel(MDL_MISSINGREGION_THXF, LABEL_DOUBLE, "missingRegionThetaXF");
        MDL::addLabel(MDL_MISSINGREGION_THY0, LABEL_DOUBLE, "missingRegionThetaY0");
        MDL::addLabel(MDL_MISSINGREGION_THYF, LABEL_DOUBLE, "missingRegionThetaYF");
        MDL::addLabel(MDL_MODELFRAC, LABEL_DOUBLE, "modelFraction");
        MDL::addLabel(MDL_NMA, LABEL_VECTOR, "NMADisplacements");
        MDL::addLabel(MDL_ORIGINX, LABEL_DOUBLE, "originX");
        MDL::addLabel(MDL_ORIGINY, LABEL_DOUBLE, "originY");
        MDL::addLabel(MDL_ORIGINZ, LABEL_DOUBLE, "originZ");
        MDL::addLabel(MDL_PMAX, LABEL_DOUBLE, "pMax", "Pmax", "sumP");
        MDL::addLabel(MDL_PSD, LABEL_STRING, "powerSpectrum");
        MDL::addLabel(MDL_Q0, LABEL_DOUBLE, "Q0");
        MDL::addLabel(MDL_RANDOMSEED, LABEL_INT, "randomSeed");
        MDL::addLabel(MDL_REF3D, LABEL_INT, "ref3d");
        MDL::addLabel(MDL_REF, LABEL_INT, "ref", "Ref");
        MDL::addLabel(MDL_REFMD, LABEL_STRING, "referenceMetaData");
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
        MDL::addLabel(MDL_SCALE, LABEL_DOUBLE, "scale", "Scale");
        MDL::addLabel(MDL_SELFILE, LABEL_STRING, "selfile");
        MDL::addLabel(MDL_SERIE, LABEL_STRING, "serie");
        MDL::addLabel(MDL_SHIFTX, LABEL_DOUBLE, "shiftX", "Xoff");
        MDL::addLabel(MDL_SHIFTY, LABEL_DOUBLE, "shiftY", "Yoff");
        MDL::addLabel(MDL_SHIFTZ, LABEL_DOUBLE, "shiftZ", "Zoff");
        MDL::addLabel(MDL_SIGMANOISE, LABEL_DOUBLE, "sigmaNoise");
        MDL::addLabel(MDL_SIGMAOFFSET, LABEL_DOUBLE, "sigmaOffset");
        MDL::addLabel(MDL_SIGNALCHANGE, LABEL_DOUBLE, "signalChange");
        MDL::addLabel(MDL_SPHERICALABERRATION, LABEL_DOUBLE, "sphericalAberration");
        MDL::addLabel(MDL_STDDEV, LABEL_DOUBLE, "stddev");
        MDL::addLabel(MDL_SUM, LABEL_DOUBLE, "sum");
        MDL::addLabel(MDL_SUMWEIGHT, LABEL_DOUBLE, "sumWeight");
        MDL::addLabel(MDL_SYMNO, LABEL_INT, "symNo");
        MDL::addLabel(MDL_TRANSFORMATIONMTRIX, LABEL_VECTOR, "transMat");
        MDL::addLabel(MDL_VOLTAGE, LABEL_DOUBLE, "voltage");
        MDL::addLabel(MDL_WEIGHT, LABEL_DOUBLE, "weight", "Weight");
        MDL::addLabel(MDL_WROBUST, LABEL_DOUBLE, "wRobust");
        MDL::addLabel(MDL_XINT, LABEL_INT, "Xcoor", "<X position>");
        MDL::addLabel(MDL_X, LABEL_DOUBLE, "X");
        MDL::addLabel(MDL_YINT, LABEL_INT, "Ycoor", "<Y position>");
        MDL::addLabel(MDL_Y, LABEL_DOUBLE, "Y");
        MDL::addLabel(MDL_ZINT, LABEL_INT, "Zcoor");
        MDL::addLabel(MDL_Z, LABEL_DOUBLE, "Z");
        MDL::addLabel(MDL_ZSCORE, LABEL_DOUBLE, "Zscore");

    }

    ~MDLabelStaticInit()
    {}
    friend class MDL;
};

/**Just an utility function */
bool vectorContainsLabel(const std::vector<MDLabel>& labelsVector, const MDLabel label);

//Just an struct to store type and string alias
class MDLabelData
{
public:
    MDLabelType type;
    String str;
    //Default constructor
    MDLabelData()
    {
        type = LABEL_NOTYPE;
    }

    MDLabelData(MDLabelType t, const String &s)
    {
        type = t;
        str = s;
    }
}
;//close class MDLabelData

/** Union to store values */
typedef union
 {
    bool boolValue;
    int intValue;
    long int longintValue;
    double doubleValue;
    String * stringValue;
    std::vector<double> * vectorValue;
} ObjectData;

/** Class to hold the labels values and type
 *
 */
class MDObject
{
private:
  ObjectData data;

    void labelTypeCheck(MDLabelType checkingType) const;
    void copy(const MDObject &obj);

public:
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
    MDObject(MDLabel label, const long int longintValue);
    MDObject(MDLabel label, const float &floatValue);
    MDObject(MDLabel label, const char * &charValue);

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
    void  getValue(long int &lv) const;
    void  getValue(float &floatvalue) const;
    void  getValue(char*  &charvalue) const;

    void  setValue(const int &iv);
    void  setValue(const double &dv);
    void  setValue(const bool &bv);
    void  setValue(const String &sv);
    void  setValue(const std::vector<double> &vv);
    void  setValue(const long int &lv);
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

    void toStream(std::ostream &os, bool withFormat = false, bool isSql=false) const;
    String toString(bool withFormat = false, bool isSql=false) const;
    bool fromStream(std::istream &is);
    friend std::istream& operator>> (std::istream& is, MDObject &value);
    bool fromString(const String &str);

    friend class MDSql;
}
; //close class MDValue

/** Explicit instantiation */
template class std::vector<MDObject *>;

/** Class for holding an entire row of posible MDObject */
class MDRow: public std::vector<MDObject*>
{
public:
    //MDObject & operator [](MDLabel label);
    bool containsLabel(MDLabel label) const;
    MDObject * getObject(MDLabel label);
    template <typename T>
    bool getValue(MDLabel label, T &d) const
    {
        for (const_iterator it = begin(); it != end(); ++it)
            if ((*it)->label == label)
            {
                (*it)->getValue(d);
                return true;
            }

        return false;
    }
    template <typename T>
    void setValue(MDLabel label, const T &d)
    {
        for (const_iterator it = begin(); it != end(); ++it)
            if ((*it)->label == label)
            {
                (*it)->setValue(d);
                return;
            }
        push_back(new MDObject(label, d));
    }

    /** Copy constructor */
    MDRow(const MDRow & row);
    MDRow();
    MDRow& operator = (const MDRow &row);
    ~MDRow();

private:
    void copy(const MDRow &row);


};

#endif
