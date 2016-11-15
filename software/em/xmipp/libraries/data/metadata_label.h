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
    MDL_GATHER_ID, /// Special label to be used when gathering MDs in MpiMetadataPrograms

    MDL_ANGLE_PSI, ///< Psi angle of an image (double,degrees)
    MDL_ANGLE_PSI2, ///< Psi angle of an image (double,degrees)
    MDL_ANGLE_PSI_DIFF, ///< difference between psi angles (double,degrees)
    MDL_ANGLE_ROT, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLE_ROT2, ///< Rotation angle of an image (double,degrees)
    MDL_ANGLE_ROT_DIFF, ///< difference between rot angles (double,degrees)
    MDL_ANGLE_TILT, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLE_TILT2, ///< Tilting angle of an image (double,degrees)
    MDL_ANGLE_TILT_DIFF, ///< difference between tilt angles (double,degrees)
    MDL_ANGLE_DIFF0, ///< difference between two angles (double,degrees)
    MDL_ANGLE_DIFF, ///< difference between two angles (double,degrees)
    MDL_ANGLE_DIFF2, ///< difference between two angles (double,degrees)
    MDL_ANGLE_Y,   ///< Angle between y-axis and tilt-axis (double, degrees) for untilted micrographs
    MDL_ANGLE_Y2,   ///< Angle between y-axis and tilt-axis (double, degrees) for tilted micrographs
    MDL_APPLY_SHIFT,///<Apply shift when project the volume ,
    MDL_AVG, ///< average value (double)
    MDL_AVG_CHANGES_ORIENTATIONS, /// Average change in angular orientation (double degrees)
    MDL_AVG_CHANGES_OFFSETS, /// Average change in offset (double pixels)
    MDL_AVG_CHANGES_CLASSES, /// Average change in class assignment(double dimensionaless)
    MDL_AVGPMAX, ///< Average (per class) of the maximum value of normalized probability function) (double)
    MDL_BGMEAN, ///< Mean background value for an image
    MDL_BLOCK_NUMBER, ///< Current block number (for incremental EM)

    MDL_CL2D_CHANGES, ///< Number of changes between iterations
    MDL_CL2D_SIMILARITY, ///< Average cross-correlation for the image (double)
    MDL_CLASS_COUNT, ///< Number of images assigned to the same class as this image
    MDL_CLASS_PERCENTAGE, ///< Percentage of images assigned to the same class as this image
    MDL_CLASSIFICATION_DATA, ///< Data vector for classification (vector double)
    MDL_CLASSIFICATION_DATA_SIZE, ///< Size of data vectors for classification (int)
    MDL_CLASSIFICATION_DPR_05, ///< Differential Phase Residual evaluated at FRC=0.5
    MDL_CLASSIFICATION_INTRACLASS_DISTANCE, ///< Average intraclass distance (double)
    MDL_CLASSIFICATION_FRC_05, ///< Digital frequency at which the FRC drops below 0.5 (double)
    MDL_COMMENT, ///< Serve to make annotations on the metadata row
    MDL_COST, ///< Cost for the image (double)
    MDL_COST_PERCENTILE, ///< Cost percentile for the image (double)
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
    MDL_CRYSTAL_NOISE_SHIFT , ///< noise if center of unit cell (vector double)

    MDL_CONTINUOUS_X, ///< X shift of continuous assignment
    MDL_CONTINUOUS_Y, ///< Y shift of continuous assignment
    MDL_CONTINUOUS_FLIP, ///< Flip of continuous assignment
    MDL_CONTINUOUS_GRAY_A, ///< a value of continuous assignment
    MDL_CONTINUOUS_GRAY_B, ///< b value of continuous assignment
    MDL_CONTINUOUS_SCALE_X, ///< a value of continuous assignment
    MDL_CONTINUOUS_SCALE_Y, ///< b value of continuous assignment

    MDL_CTF_DATA_PHASE_FLIPPED, // Is the Data Phase-Flippled?
    MDL_CTF_CORRECTED, // Is the CTF corrected?
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
    MDL_CTF_ENVELOPE, //<Envelope function
    MDL_CTF_ENVELOPE_PLOT, //<Envelope function
    MDL_CTF_LENS_STABILITY, ///< Lens stability
    MDL_CTF_CONVERGENCE_CONE, ///< Convergence cone
    MDL_CTF_LONGITUDINAL_DISPLACEMENT, ///< Longitudinal displacement
    MDL_CTF_TRANSVERSAL_DISPLACEMENT, ///< Transversal displacemente
    MDL_CTF_Q0, ///< Inelastic absorption
    MDL_CTF_K, ///< CTF gain
	MDL_CTF_ENV_R0, ///< CTF Envelope polynomial parameter
	MDL_CTF_ENV_R1, ///< CTF Envelope polynomial parameter
	MDL_CTF_ENV_R2, ///< CTF Envelope polynomial parameter
    MDL_CTF_BG_GAUSSIAN_K, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_SIGMAU, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_SIGMAV, ///< CTF Background parameter
    MDL_CTF_BG_GAUSSIAN_CU, ///<CTF_ CTF Background parameter
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
	MDL_CTF_BG_R1, ///< CTF Background polynomial parameter
	MDL_CTF_BG_R2, ///< CTF Background polynomial parameter
	MDL_CTF_BG_R3, ///< CTF Background polynomial parameter
    MDL_CTF_CRIT_NONASTIGMATICVALIDITY, ///< Maximum frequency (in Angstroms) at which non-astigmatic CTF correction is valid
    MDL_CTF_CRIT_PSDCORRELATION90, ///< PSD correlation at 90 degrees
    MDL_CTF_CRIT_FIRSTZERORATIO, ///< First zero ratio
    MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_RATIO, ///< Ratio sigma(firstMinimum)/sigma(firstZero)
    MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_DIFF_RATIO, ///< Ratio sigma(firstMinimum-firstZero)/sigma(firstZero)
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
    MDL_CUMULATIVE_SSNR, ///<Cumulative SSNR (double)
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
    MDL_FRAME_ID, ///< Unique id of frame inside a Movie
    MDL_IDX, ///< Index within a list (size_t)
    MDL_IMAGE, ///< Name of an image (std::string)
    MDL_IMAGE_COVARIANCE, ///< Name of the covariance imagee associated to this image
    MDL_IMAGE_IDX, ///< Index of an image within a list (size_t)
    MDL_IMAGE_ORIGINAL, ///< Name of an image from which MDL_IMAGE is coming from
    MDL_IMAGE_REF, ///< Name of of the class image from which MDL_IMAGE is coming from
    MDL_IMAGE_RESIDUAL, ///< Name of a residual image associated to this image
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
    MDL_MAXCC_PERCENTILE, ///< Percentile of the maximum cross-correlation for the image (double)
    MDL_MAX, ///< Maximum value (double)
    MDL_MICROGRAPH, ///< Name of a micrograph (std::string)
    MDL_MICROGRAPH_ID, ///< Micrograph unique id for reference (MDL_ITEM_ID should be used for Micrographs list)
    MDL_MICROGRAPH_MOVIE, ///< Name of a movie (std::string)
    MDL_MICROGRAPH_MOVIE_ID, ///< Unique identifier of a movie.
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
    MDL_NEIGHBORHOOD_RADIUS, ///< Radius of the neighborhood (radians)
    MDL_NMA, ///< Normal mode displacements (vector double)
    MDL_NMA_COLLECTIVITY, ///< NMA Collectivity of a given mode
    MDL_NMA_ATOMSHIFT, ///< NMA Atom shift in Angstroms
    MDL_NMA_ENERGY, ///< NMA energy contained in the NMA displacement vector
    MDL_NMA_MINRANGE, ///< Minimum value observed for a given NMA mode
    MDL_NMA_MAXRANGE, ///< Maximum value observed for a given NMA mode
    MDL_NMA_MODEFILE, ///< File with an NMA mode
    MDL_NMA_SCORE, ///< NMA Score of a given mode
    MDL_NOISE_ANGLES, ///< Noise description for projected angles
    MDL_NOISE_PARTICLE_COORD, ///< Noise description for particle's center coordenates (when projecting)
    MDL_NOISE_COORD,  //Use instead of MDL_NOISE_PARTICLE_COORD in future
    MDL_NOISE_PIXEL_LEVEL, ///< Noise description for pixels' gray level (when projecting)
    MDL_ORDER, /// auxiliary label to be used as an index (long)
    MDL_ORIGIN_X, ///< Origin for the image in the X axis (double)
    MDL_ORIGIN_Y, ///< Origin for the image in the Y axis (double)
    MDL_ORIGIN_Z, ///< Origin for the image in the Z axis (double)

    MDL_OPTICALFLOW_MEANX, ///<Mean of the movement in x direction of the motion map> (double)
    MDL_OPTICALFLOW_MEANY, ///<Mean of the movement in y direction of the motion map> (double)
    MDL_OPTICALFLOW_STDX, ///<Standatd deviation of the movement in x direction of the motion map> (double)
    MDL_OPTICALFLOW_STDY, ///<Standard deviation of the movement in y direction of the motion map> (double)

    MDL_PARTICLE_ID, ///< Particle unique identifier for reference. (The MDL_ITEM_ID should be used when particle list)
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
    MDL_PICKING_AUTOPARTICLES_SIZE, ///< Number of automatic particles picked
    MDL_PICKING_MANUALPARTICLES_SIZE, ///< Number of manual particles picked
    MDL_PMAX, ///< Maximum value of normalized probability function (now called "Pmax/sumP") (double)
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
    MDL_RESOLUTION_FREQ2, ///< Frequency in 1/A squared (double)
    MDL_RESOLUTION_FREQREAL, ///< Frequency in A (double)
    MDL_RESOLUTION_LOG_STRUCTURE_FACTOR, ///<Logarithm of the structure factor
    MDL_RESOLUTION_SSNR, ///<Fourier shell correlation (double)
    MDL_RESOLUTION_STRUCTURE_FACTOR, ///<Structure factor
    MDL_RESOLUTION_RFACTOR, /// rfactor

    MDL_SAMPLINGRATE, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_ORIGINAL, ///< original sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_X, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_Y, ///< sampling rate in A/pixel (double)
    MDL_SAMPLINGRATE_Z, ///< sampling rate in A/pixel (double)

    MDL_SCALE, ///< scaling factor for an image or volume (double)
    MDL_SCORE_BY_PCA_RESIDUAL_PROJ,
    MDL_SCORE_BY_PCA_RESIDUAL_EXP,
    MDL_SCORE_BY_PCA_RESIDUAL,
    MDL_SCORE_BY_ALIGNABILITY, ///< score by alignability (double)
    MDL_SCORE_BY_ALIGNABILITY_PRECISION, ///< score by alignability (double)
    MDL_SCORE_BY_ALIGNABILITY_ACCURACY, ///< score by alignability (double)
    MDL_SCORE_BY_ALIGNABILITY_PRECISION_EXP, ///< score by alignability experimental particles (double)
    MDL_SCORE_BY_ALIGNABILITY_PRECISION_REF, ///< score by alignability references (double)
    MDL_SCORE_BY_ALIGNABILITY_ACCURACY_EXP, ///< score by alignability experimental particles (double)
    MDL_SCORE_BY_ALIGNABILITY_ACCURACY_REF, ///< score by alignability references (double)
    MDL_SCORE_BY_ALIGNABILITY_NOISE, ///< score by alignability noise (double)
    MDL_SCORE_BY_MIRROR, ///< score by mirror (double)
    MDL_SCORE_BY_ZSCORE,
    MDL_SELFILE, ///< Name of an image (std::string)
    MDL_SERIE, ///< A collection of micrographs, e.g. a tilt serie (std::string)
    MDL_SHIFT_X, ///< Shift for the image in the X axis (double)
    MDL_SHIFT_X2, ///< Shift for the image in the X axis (double)
    MDL_SHIFT_X_DIFF, ///< difference in Shift along X axis (double)
    MDL_SHIFT_Y, ///< Shift for the image in the Y axis (double)
    MDL_SHIFT_Y2, ///< Shift for the image in the Y axis (double)
    MDL_SHIFT_Y_DIFF, ///< difference in Shift along  Y axis (double)
    MDL_SHIFT_Z, ///< Shift for the image in the Z axis (double)
    MDL_SHIFT_DIFF0, ///< shift difference (double)
    MDL_SHIFT_DIFF, ///< shift difference (double)
    MDL_SHIFT_DIFF2, ///< shift difference (double)
    MDL_SIGMANOISE, ///< Standard deviation of the noise in ML model
    MDL_SIGMAOFFSET, ///< Standard deviation of the offsets in ML model
    MDL_SIGNALCHANGE, ///< Signal change for an image
    MDL_STDDEV, ///<stdandard deviation value (double)
    MDL_STAR_COMMENT, ///< A comment for this object /*** NOTE THIS IS A SPECIAL CASE AND SO IS TREATED ***/
    MDL_SUM, ///< Sum of elements of a given type (double) [this is a genereic type do not use to transfer information to another program]
    MDL_SUMWEIGHT, ///< Sum of all weights in ML model
    MDL_SYMNO, ///< Symmetry number for a projection (used in ART)
    MDL_TOMOGRAM_VOLUME, ///< Name for the reconstructed tomogram volume (std::string)
    MDL_TOMOGRAMMD, ///< Name for a Metadata file (std::string)
    MDL_TRANSFORM_MATRIX, ///< transformation matrix in numpy string format or space separated (std::string)

    MDL_TEST_SIZE,// < number of test assigned to a program

    MDL_VOLUME_SCORE_SUM, /// < Score corresponding to the sum of cc with cc>threshold
    MDL_VOLUME_SCORE_SUM_TH, ///< Score corresponding to the sum of cc-threshold with cc>threshold
    MDL_VOLUME_SCORE_MEAN, ///< Score corresponding to the mean of cc with cc>threshold
    MDL_VOLUME_SCORE_MIN, ///< Score corresponding to the min of cc with cc>threshold
    MDL_VOLUME_SCORE1,/// < Score 1 for volumes
    MDL_VOLUME_SCORE2,/// < Score 2 for volumes
    MDL_VOLUME_SCORE3,/// < Score 3 for volumes
    MDL_VOLUME_SCORE4,/// < Score 4 for volumes
    MDL_WEIGHT, ///< Weight assigned to the image (double)
    MDL_WEIGHT_P, ///< Weight assigned to the image accordint to its clusterability with a significance with respect noise (double)
    MDL_WEIGHT_CONTINUOUS2, ///< Weight due to angular continuous assignment
    MDL_WEIGHT_JUMPER0, ///< Weight due to angular jumping
    MDL_WEIGHT_JUMPER, ///< Weight due to angular jumping
    MDL_WEIGHT_JUMPER2, ///< Weight due to angular jumping
    MDL_WEIGHT_SIGNIFICANT, ///< Weight due to Angular significance
    MDL_WEIGHT_SSNR, ///< Weight due to SSNR
    MDL_WEIGHT_PRECISION_ALIGNABILITY, ///< Weight due to Alignability Precision
    MDL_WEIGHT_ALIGNABILITY, ///< Weight due to Alignability Precision and Accuracy
    MDL_WEIGHT_ACCURACY_ALIGNABILITY, ///< Weight due to Alignability Accuracy
    MDL_WEIGHT_PRECISION_MIRROR, ///< Weight due to Mirror Precision
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
    MDL_ZSCORE_HISTOGRAM, ///< Z Score (double)
    MDL_ZSCORE_RESMEAN, ///< Z Score of the mean of the residuals (double)
    MDL_ZSCORE_RESVAR, ///< Z Score of the stddev of the residuals (double)
    MDL_ZSCORE_RESCOV, ///< Z Score of the covariance matrix of the residuals (double)
    MDL_ZSCORE_SHAPE1, ///< Z Score (double)
    MDL_ZSCORE_SHAPE2, ///< Z Score (double)
    MDL_ZSCORE_SNR1, ///< Z Score (double)
    MDL_ZSCORE_SNR2, ///< Z Score (double)
    MDL_ZSIZE, ///< Z size (int)

    /** RELION labels */
    RLN_AREA_ID, ///< ID for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
    RLN_AREA_NAME, ///< Name for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
    RLN_COMMENT, // The RLN_COMMENT is handled specially as well

    RLN_CTF_BFACTOR, ///< B-factor
    RLN_CTF_SCALEFACTOR, ///< linear scale-factor
    RLN_CTF_SAMPLING_RATE, ///< Sampling rate
    RLN_CTF_VOLTAGE, ///< Microscope voltage (kV)
    RLN_CTF_DEFOCUSU, ///< Defocus U (Angstroms)
    RLN_CTF_DEFOCUSV, ///< Defocus V (Angstroms)
    RLN_CTF_DEFOCUS_ANGLE, ///< Defocus angle (degrees)
    RLN_CTF_CS, ///< Spherical aberration
    RLN_CTF_CA, ///< Chromatic aberration
    RLN_CTF_DETECTOR_PIXEL_SIZE, ///< Pixel size for detector as used in CTF-determination
    RLN_CTF_ENERGY_LOSS, ///< Energy loss
    RLN_CTF_FOM, ///< ctffind3 FOM (CC) for quality of CTF-fit
    RLN_CTF_IMAGE, ///< name of an image describing the CTF model
    RLN_CTF_LENS_STABILITY, ///< Lens stability
    RLN_CTF_MAGNIFICATION, ///< Magnification used for CTF-determination
    RLN_CTF_CONVERGENCE_CONE, ///< Convergence cone
    RLN_CTF_LONGITUDINAL_DISPLACEMENT, ///< Longitudinal displacement
    RLN_CTF_TRANSVERSAL_DISPLACEMENT, ///< Transversal displacemente
    RLN_CTF_Q0, ///< Amplitude contrast
    RLN_CTF_K, ///< CTF gain
    RLN_CTF_VALUE, ///< CTF value

    RLN_IMAGE_NAME,
    RLN_IMAGE_RECONSTRUCT_NAME,
    RLN_IMAGE_ID,
    RLN_IMAGE_ENABLED,
    RLN_IMAGE_DATATYPE,
    RLN_IMAGE_DIMENSIONALITY,
    RLN_IMAGE_BEAMTILT_X,
    RLN_IMAGE_BEAMTILT_Y,
    RLN_IMAGE_BEAMTILT_GROUP,
    RLN_IMAGE_COORD_X,
    RLN_IMAGE_COORD_Y,
    RLN_IMAGE_COORD_Z,
    RLN_IMAGE_FRAME_NR,
    RLN_IMAGE_MAGNIFICATION_CORRECTION,
    RLN_IMAGE_NORM_CORRECTION,
    RLN_IMAGE_ORI_NAME,
    RLN_IMAGE_SAMPLINGRATE,
    RLN_IMAGE_SAMPLINGRATE_X,
    RLN_IMAGE_SAMPLINGRATE_Y,
    RLN_IMAGE_SAMPLINGRATE_Z,
    RLN_IMAGE_SIZE,
    RLN_IMAGE_SIZEX,
    RLN_IMAGE_SIZEY,
    RLN_IMAGE_SIZEZ,
    RLN_IMAGE_STATS_MIN,
    RLN_IMAGE_STATS_MAX,
    RLN_IMAGE_STATS_AVG,
    RLN_IMAGE_STATS_STDDEV,
    RLN_IMAGE_STATS_SKEW,
    RLN_IMAGE_STATS_KURT,
    RLN_IMAGE_WEIGHT,

    RLN_MATRIX_1_1,
    RLN_MATRIX_1_2,
    RLN_MATRIX_1_3,
    RLN_MATRIX_2_1,
    RLN_MATRIX_2_2,
    RLN_MATRIX_2_3,
    RLN_MATRIX_3_1,
    RLN_MATRIX_3_2,
    RLN_MATRIX_3_3,

    RLN_MICROGRAPH_ID,
    RLN_MICROGRAPH_NAME,
    RLN_MICROGRAPH_TILT_ANGLE,
    RLN_MICROGRAPH_TILT_AXIS_DIRECTION,
    RLN_MICROGRAPH_TILT_AXIS_OUTOFPLANE,

    RLN_MLMODEL_ACCURACY_ROT,
    RLN_MLMODEL_ACCURACY_TRANS,
    RLN_MLMODEL_AVE_PMAX,
    RLN_MLMODEL_CURRENT_RESOLUTION,
    RLN_MLMODEL_CURRENT_SIZE,
    RLN_MLMODEL_DATA_VS_PRIOR_REF,
    RLN_MLMODEL_DIMENSIONALITY,
    RLN_MLMODEL_DIMENSIONALITY_DATA,
    RLN_MLMODEL_DIFF2_HALVES_REF,
    RLN_MLMODEL_FSC_HALVES_REF,
    RLN_MLMODEL_GROUP_NAME,
    RLN_MLMODEL_GROUP_NO,
    RLN_MLMODEL_GROUP_NR_PARTICLES,
    RLN_MLMODEL_GROUP_SCALE_CORRECTION,
    RLN_MLMODEL_INTERPOLATOR,
    RLN_MLMODEL_LL,
    RLN_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION,
    RLN_MLMODEL_NORM_CORRECTION_AVG,
    RLN_MLMODEL_NR_CLASSES,
    RLN_MLMODEL_NR_GROUPS,
    RLN_MLMODEL_ORIGINAL_SIZE,
    RLN_MLMODEL_ORIENTABILITY_CONTRIBUTION,
    RLN_MLMODEL_PADDING_FACTOR,
    RLN_MLMODEL_PDF_CLASS,
    RLN_MLMODEL_PRIOR_OFFX_CLASS,
    RLN_MLMODEL_PRIOR_OFFY_CLASS,
    RLN_MLMODEL_PDF_ORIENT,
    RLN_MLMODEL_PIXEL_SIZE,
    RLN_MLMODEL_POWER_REF,
    RLN_MLMODEL_PRIOR_MODE,
    RLN_MLMODEL_SIGMA_OFFSET,
    RLN_MLMODEL_SIGMA_ROT,
    RLN_MLMODEL_SIGMA_TILT,
    RLN_MLMODEL_SIGMA_PSI,
    RLN_MLMODEL_REF_IMAGE,
    RLN_MLMODEL_SIGMA2_NOISE,
    RLN_MLMODEL_SIGMA2_REF,
    RLN_MLMODEL_SSNR_REF,
    RLN_MLMODEL_TAU2_FUDGE_FACTOR,
    RLN_MLMODEL_TAU2_REF,

    RLN_OPTIMISER_ACCURACY_ROT,
    RLN_OPTIMISER_ACCURACY_TRANS,
    RLN_OPTIMISER_ADAPTIVE_FRACTION,
    RLN_OPTIMISER_ADAPTIVE_OVERSAMPLING,
    RLN_OPTIMISER_AUTO_LOCAL_HP_ORDER,
    RLN_OPTIMISER_AVAILABLE_MEMORY,
    RLN_OPTIMISER_BEST_RESOL_THUS_FAR,
    RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS,
    RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS,
    RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES,
    RLN_OPTIMISER_COARSE_SIZE,
    RLN_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED,
    RLN_OPTIMISER_DATA_STARFILE,
    RLN_OPTIMISER_DO_AUTO_REFINE,
    RLN_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES,
    RLN_OPTIMISER_DO_CORRECT_CTF,
    RLN_OPTIMISER_DO_CORRECT_MAGNIFICATION,
    RLN_OPTIMISER_DO_CORRECT_NORM,
    RLN_OPTIMISER_DO_CORRECT_SCALE,
    RLN_OPTIMISER_DO_REALIGN_MOVIES,
    RLN_OPTIMISER_DO_MAP,
    RLN_OPTIMISER_DO_SOLVENT_FLATTEN,
    RLN_OPTIMISER_DO_SKIP_ALIGN,
    RLN_OPTIMISER_DO_SKIP_ROTATE,
    RLN_OPTIMISER_DO_SPLIT_RANDOM_HALVES,
    RLN_OPTIMISER_DO_ZERO_MASK,
    RLN_OPTIMISER_FIX_SIGMA_NOISE,
    RLN_OPTIMISER_FIX_SIGMA_OFFSET,
    RLN_OPTIMISER_FIX_TAU,
    RLN_OPTIMISER_HAS_CONVERGED,
    RLN_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT,
    RLN_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO,
    RLN_OPTIMISER_HIGHRES_LIMIT_EXP,
    RLN_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK,
    RLN_OPTIMISER_INCR_SIZE,
    RLN_OPTIMISER_ITERATION_NO,
    RLN_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES,
    RLN_OPTIMISER_MAGNIFICATION_RANGE,
    RLN_OPTIMISER_MAGNIFICATION_STEP,
    RLN_OPTIMISER_MAX_COARSE_SIZE,
    RLN_OPTIMISER_MAX_NR_POOL,
    RLN_OPTIMISER_MODEL_STARFILE,
    RLN_OPTIMISER_MODEL_STARFILE2,
    RLN_OPTIMISER_NR_ITERATIONS,
    RLN_OPTIMISER_NR_ITER_WO_RESOL_GAIN,
    RLN_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES,
    RLN_OPTIMISER_OUTPUT_ROOTNAME,
    RLN_OPTIMISER_PARTICLE_DIAMETER,
    RLN_OPTIMISER_RADIUS_MASK_3D_MAP,
    RLN_OPTIMISER_RADIUS_MASK_EXP_PARTICLES,
    RLN_OPTIMISER_RANDOM_SEED,
    RLN_OPTIMISER_REFS_ARE_CTF_CORRECTED,
    RLN_OPTIMISER_SAMPLING_STARFILE,
    RLN_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES,
    RLN_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS,
    RLN_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS,
    RLN_OPTIMISER_SOLVENT_MASK_NAME,
    RLN_OPTIMISER_SOLVENT_MASK2_NAME,
    RLN_OPTIMISER_TAU_SPECTRUM_NAME,
    RLN_OPTIMISER_USE_TOO_COARSE_SAMPLING,
    RLN_OPTIMISER_WIDTH_MASK_EDGE,

    RLN_ORIENT_FLIP,
    RLN_ORIENT_ID,
    RLN_ORIENT_ORIGIN_X,
    RLN_ORIENT_ORIGIN_X_PRIOR,
    RLN_ORIENT_ORIGIN_Y,
    RLN_ORIENT_ORIGIN_Y_PRIOR,
    RLN_ORIENT_ORIGIN_Z,
    RLN_ORIENT_ORIGIN_Z_PRIOR,
    RLN_ORIENT_ROT,
    RLN_ORIENT_ROT_PRIOR,
    RLN_ORIENT_TILT,
    RLN_ORIENT_TILT_PRIOR,
    RLN_ORIENT_PSI,
    RLN_ORIENT_PSI_PRIOR,

    RLN_PARTICLE_AUTOPICK_FOM,
    RLN_PARTICLE_CLASS,
    RLN_PARTICLE_DLL,
    RLN_PARTICLE_ID,
    RLN_PARTICLE_FOM,
    RLN_PARTICLE_KL_DIVERGENCE,
    RLN_PARTICLE_RANDOM_SUBSET,
    RLN_PARTICLE_NAME,
    RLN_PARTICLE_ORI_NAME,
    RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES,
    RLN_PARTICLE_NR_FRAMES,
    RLN_PARTICLE_PMAX,

    RLN_POSTPROCESS_BFACTOR,
    RLN_POSTPROCESS_FINAL_RESOLUTION,
    RLN_POSTPROCESS_FSC_TRUE,
    RLN_POSTPROCESS_FSC_MASKED,
    RLN_POSTPROCESS_FSC_UNMASKED,
    RLN_POSTPROCESS_FSC_RANDOM_MASKED,
    RLN_POSTPROCESS_GUINIER_FIT_CORRELATION,
    RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT,
    RLN_POSTPROCESS_GUINIER_FIT_SLOPE,
    RLN_POSTPROCESS_GUINIER_VALUE_IN,
    RLN_POSTPROCESS_GUINIER_VALUE_INVMTF,
    RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED,
    RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED,
    RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT,
    RLN_POSTPROCESS_GUINIER_RESOL_SQUARED,
    RLN_POSTPROCESS_MTF_VALUE, ///< Detector MTF value

    RLN_SAMPLING_IS_3D,
    RLN_SAMPLING_IS_3D_TRANS,
    RLN_SAMPLING_HEALPIX_ORDER,
    RLN_SAMPLING_LIMIT_TILT,
    RLN_SAMPLING_OFFSET_RANGE,
    RLN_SAMPLING_OFFSET_STEP,
    RLN_SAMPLING_PERTURB,
    RLN_SAMPLING_PERTURBATION_FACTOR,
    RLN_SAMPLING_PRIOR_MODE,
    RLN_SAMPLING_PSI_STEP,
    RLN_SAMPLING_SIGMA_ROT,
    RLN_SAMPLING_SIGMA_TILT,
    RLN_SAMPLING_SIGMA_PSI,
    RLN_SAMPLING_SYMMETRY,

    RLN_SELECTED,
    RLN_SELECT_PARTICLES_ZSCORE,
    RLN_SORTED_IDX,
    RLN_PERFRAME_CUMULATIVE_WEIGHT,
    RLN_PERFRAME_RELATIVE_WEIGHT,

    RLN_RESOLUTION,
    RLN_RESOLUTION_ANGSTROM,
    RLN_RESOLUTION_INVPIXEL,
    RLN_SPECTRAL_IDX,

    /** BSOFT labels */
    //BSOFT_COMMENT,
    BSOFT_ID,
    BSOFT_PROJECT,
    BSOFT_FIELD,
    BSOFT_FIELD_ID,
    BSOFT_MAP,
    BSOFT_MAP_ID,
    BSOFT_MAP_REFERENCE,
    BSOFT_MAP_RECONSTRUCTION,
    BSOFT_MAP_TRANSFORM_FILE,
    BSOFT_MAP_POWERSPEC_FILE,
    BSOFT_MAP_SIZE_X,
    BSOFT_MAP_SIZE_Y,
    BSOFT_MAP_SIZE_Z,
    BSOFT_MAP_ORIGIN_X,
    BSOFT_MAP_ORIGIN_Y,
    BSOFT_MAP_ORIGIN_Z,
    BSOFT_MAP_SCALE_X,
    BSOFT_MAP_SCALE_Y,
    BSOFT_MAP_SCALE_Z,
    BSOFT_MAP_VOXEL_SIZE,
    BSOFT_MAP_SELECT,
    BSOFT_MAP_FOM,
    BSOFT_MAP_MAGNIFICATION,
    BSOFT_MAP_VIEW_X,
    BSOFT_MAP_VIEW_Y,
    BSOFT_MAP_VIEW_Z,
    BSOFT_MAP_VIEW_ANGLE,
    BSOFT_MAP_BACK_RWEIGHT,
    BSOFT_MAP_MODEL,
    BSOFT_MAP_SYMMETRY,
    //BSOFT_MICROGRAPH,
    BSOFT_MICROGRAPH_FILE,
    BSOFT_MICROGRAPH_PARTICLE_FILE,
    BSOFT_MICROGRAPH_FILAMENT_FILE,
    BSOFT_MICROGRAPH_TRANSFORM_FILE,
    BSOFT_MICROGRAPH_POWERSPEC_FILE,
    BSOFT_MICROGRAPH_ID,
    BSOFT_MICROGRAPH_FIELD_ID,
    BSOFT_MICROGRAPH_NUMBER,
    BSOFT_MICROGRAPH_SELECT,
    BSOFT_MICROGRAPH_FOM,
    BSOFT_MICROGRAPH_MAGNIFICATION,
    BSOFT_MICROGRAPH_SAMPLING,
    BSOFT_MICROGRAPH_PIXEL,
    BSOFT_MICROGRAPH_UNITS,
    BSOFT_MICROGRAPH_DOSE,
    BSOFT_MICROGRAPH_ORIGIN_X,
    BSOFT_MICROGRAPH_ORIGIN_Y,
    BSOFT_MICROGRAPH_ORIGIN_Z,
    BSOFT_MICROGRAPH_SCALE_X,
    BSOFT_MICROGRAPH_SCALE_Y,
    BSOFT_MICROGRAPH_SCALE_Z,
    BSOFT_MICROGRAPH_TILT_AXIS,
    BSOFT_MICROGRAPH_TILT_ANGLE,
    BSOFT_MICROGRAPH_LEVEL_ANGLE,
    BSOFT_MICROGRAPH_ROT_ANGLE,
    BSOFT_MICROGRAPH_VIEW_X,
    BSOFT_MICROGRAPH_VIEW_Y,
    BSOFT_MICROGRAPH_VIEW_Z,
    BSOFT_MICROGRAPH_VIEW_ANGLE,
    BSOFT_MICROGRAPH_MATRIX_1_1,
    BSOFT_MICROGRAPH_MATRIX_1_2,
    BSOFT_MICROGRAPH_MATRIX_1_3,
    BSOFT_MICROGRAPH_MATRIX_2_1,
    BSOFT_MICROGRAPH_MATRIX_2_2,
    BSOFT_MICROGRAPH_MATRIX_2_3,
    BSOFT_MICROGRAPH_MATRIX_3_1,
    BSOFT_MICROGRAPH_MATRIX_3_2,
    BSOFT_MICROGRAPH_MATRIX_3_3,
    BSOFT_MICROGRAPH_HVEC_X,
    BSOFT_MICROGRAPH_HVEC_Y,
    BSOFT_MICROGRAPH_HVEC_Z,
    BSOFT_MICROGRAPH_KVEC_X,
    BSOFT_MICROGRAPH_KVEC_Y,
    BSOFT_MICROGRAPH_KVEC_Z,
    BSOFT_MICROGRAPH_LVEC_X,
    BSOFT_MICROGRAPH_LVEC_Y,
    BSOFT_MICROGRAPH_LVEC_Z,
    BSOFT_MICROGRAPH_HELIX_AXIS,
    BSOFT_MICROGRAPH_HELIX_RISE,
    BSOFT_MICROGRAPH_HELIX_ANGLE,
    BSOFT_MICROGRAPH_HELIX_RADIUS,
    BSOFT_MICROGRAPH_VOLTAGE,
    BSOFT_MICROGRAPH_CTF_CS,
    BSOFT_MICROGRAPH_CTF_CC,
    BSOFT_MICROGRAPH_CTF_ALPHA,
    BSOFT_MICROGRAPH_CTF_DE,
    BSOFT_MICROGRAPH_CTF_AMP_CONT,
    BSOFT_MICROGRAPH_CTF_ZERO,
    BSOFT_MICROGRAPH_CTF_DEF_AVG,
    BSOFT_MICROGRAPH_CTF_DEF_DEV,
    BSOFT_MICROGRAPH_CTF_DEF_MIN,
    BSOFT_MICROGRAPH_CTF_DEF_MAX,
    BSOFT_MICROGRAPH_CTF_AST_ANG,
    BSOFT_MICROGRAPH_CTF_BASELINE,
    BSOFT_MICROGRAPH_CTF_ENVELOPE,
    BSOFT_MICROGRAPH_BOX_RADIUS,
    BSOFT_MICROGRAPH_BOX_RADIUS_X,
    BSOFT_MICROGRAPH_BOX_RADIUS_Y,
    BSOFT_MICROGRAPH_BOX_RADIUS_Z,
    BSOFT_MICROGRAPH_BAD,
    BSOFT_MICROGRAPH_BAD_RADIUS,
    BSOFT_MICROGRAPH_BAD_X,
    BSOFT_MICROGRAPH_BAD_Y,
    BSOFT_MICROGRAPH_BAD_Z,
    BSOFT_MICROGRAPH_MARKER_RADIUS,
    BSOFT_MICROGRAPH_MARKER_ID,
    BSOFT_MICROGRAPH_MARKER_X,
    BSOFT_MICROGRAPH_MARKER_Y,
    BSOFT_MICROGRAPH_MARKER_Z,
    BSOFT_MICROGRAPH_MARKER_ERROR_X,
    BSOFT_MICROGRAPH_MARKER_ERROR_Y,
    BSOFT_MICROGRAPH_MARKER_ERROR_Z,
    BSOFT_MICROGRAPH_MARKER_FOM,
    BSOFT_MICROGRAPH_FILAMENT_WIDTH,
    BSOFT_MICROGRAPH_FILNODE_RADIUS,
    BSOFT_CTF,
    BSOFT_CTF_VOLTAGE,
    BSOFT_CTF_CS,
    BSOFT_CTF_CC,
    BSOFT_CTF_ALPHA,
    BSOFT_CTF_DE,
    BSOFT_CTF_AMP,
    BSOFT_CTF_ZERO,
    BSOFT_CTF_DEF_AVG,
    BSOFT_CTF_DEF_DEV,
    BSOFT_CTF_DEF_MIN,
    BSOFT_CTF_DEF_MAX,
    BSOFT_CTF_AST_ANG,
    BSOFT_CTF_BASELINE,
    BSOFT_CTF_ENVELOPE,
    BSOFT_PARTICLE,
    BSOFT_PARTICLE_FILE,
    BSOFT_PARTICLE_NUMBER,
    BSOFT_PARTICLE_ID,
    BSOFT_PARTICLE_GROUP,
    BSOFT_PARTICLE_MG_ID,
    BSOFT_PARTICLE_MG_X,
    BSOFT_PARTICLE_MG_Y,
    BSOFT_PARTICLE_MG_Z,
    BSOFT_PARTICLE_X,
    BSOFT_PARTICLE_Y,
    BSOFT_PARTICLE_Z,
    BSOFT_PARTICLE_X_ORIGIN,
    BSOFT_PARTICLE_Y_ORIGIN,
    BSOFT_PARTICLE_Z_ORIGIN,
    BSOFT_PARTICLE_ORIGIN_X,
    BSOFT_PARTICLE_ORIGIN_Y,
    BSOFT_PARTICLE_ORIGIN_Z,
    BSOFT_PARTICLE_PSI,
    BSOFT_PARTICLE_THETA,
    BSOFT_PARTICLE_PHI,
    BSOFT_PARTICLE_OMEGA,
    BSOFT_PARTICLE_VIEW_X,
    BSOFT_PARTICLE_VIEW_Y,
    BSOFT_PARTICLE_VIEW_Z,
    BSOFT_PARTICLE_VIEW_ANGLE,
    BSOFT_PARTICLE_MAGNIF,
    BSOFT_PARTICLE_DEFOCUS,
    BSOFT_PARTICLE_DEF_DEV,
    BSOFT_PARTICLE_AST_ANG,
    BSOFT_PARTICLE_SELECT,
    BSOFT_PARTICLE_FOM,
    BSOFT_PARTICLE_FOM_CV,
    BSOFT_PARTICLE_FOM_AVG,
    BSOFT_PARTICLE_FOM_STD,
    BSOFT_PARTICLE_HANDA_FOM,
    BSOFT_PARTICLE_HANDB_FOM,
    BSOFT_PARTICLE_CC,
    BSOFT_PARTICLE_PFT_CC,
    BSOFT_PARTICLE_PRJ_CC,
    BSOFT_PARTICLE_CMP_CC,
    BSOFT_PARTICLE_RFACTORAB,
    BSOFT_PARTICLE_COVERAGE,
    BSOFT_PARTICLE_BOX_SIZE,
    BSOFT_PARTICLE_BOX_SIZE_X,
    BSOFT_PARTICLE_BOX_SIZE_Y,
    BSOFT_PARTICLE_BOX_SIZE_Z,
    BSOFT_PARTICLE_BOX_RADIUS,
    BSOFT_PARTICLE_BOX_RADIUS_X,
    BSOFT_PARTICLE_BOX_RADIUS_Y,
    BSOFT_PARTICLE_BOX_RADIUS_Z,
    BSOFT_PARTICLE_BAD,
    BSOFT_PARTICLE_BAD_RADIUS,
    BSOFT_PARTICLE_BAD_X,
    BSOFT_PARTICLE_BAD_Y,
    BSOFT_PARTICLE_BAD_Z,
    BSOFT_FILAMENT,
    BSOFT_FILAMENT_FILE,
    BSOFT_FILAMENT_ID,
    BSOFT_FILAMENT_NODE,
    BSOFT_FILAMENT_NODE_ID,
    BSOFT_FILAMENT_NODE_X,
    BSOFT_FILAMENT_NODE_Y,
    BSOFT_FILAMENT_NODE_Z,
    BSOFT_FILAMENT_WIDTH,
    BSOFT_FILNODE_RADIUS,
    BSOFT_ORIENT_ID,
    BSOFT_ORIENT_ORIGIN_X,
    BSOFT_ORIENT_ORIGIN_Y,
    BSOFT_ORIENT_ORIGIN_Z,
    BSOFT_ORIENT_VIEW_X,
    BSOFT_ORIENT_VIEW_Y,
    BSOFT_ORIENT_VIEW_Z,
    BSOFT_ORIENT_VIEW_ANGLE,
    BSOFT_ORIENT_FOM,
    BSOFT_ORIENT_SELECT,
    BSOFT_MARKER,
    BSOFT_MARKER_RADIUS,
    BSOFT_MARKER_ID,
    BSOFT_MARKER_X,
    BSOFT_MARKER_Y,
    BSOFT_MARKER_Z,
    BSOFT_MARKER_ERROR_X,
    BSOFT_MARKER_ERROR_Y,
    BSOFT_MARKER_ERROR_Z,
    BSOFT_MARKER_IMAGE,
    BSOFT_MARKER_RESIDUAL,
    BSOFT_MARKER_FOM,
    BSOFT_MARKER_SELECT,
    BSOFT_REFLEX,
    BSOFT_REFLEX_RADIUS,
    BSOFT_REFLEX_X,
    BSOFT_REFLEX_Y,
    BSOFT_REFLEX_Z,
    BSOFT_REFLEX_H,
    BSOFT_REFLEX_K,
    BSOFT_REFLEX_L,
    BSOFT_REFLEX_AMP,
    BSOFT_REFLEX_SIGAMP,
    BSOFT_REFLEX_PHI,
    BSOFT_REFLEX_SIGPHI,
    BSOFT_REFLEX_FOM,
    BSOFT_REFLEX_STATUS,
    BSOFT_LAYERLINE,
    BSOFT_LAYERLINE_NUMBER,
    BSOFT_LAYERLINE_ORDER,
    BSOFT_LAYERLINE_DISTANCE,
    BSOFT_LAYERLINE_FREQ,
    BSOFT_LAYERLINE_AMP,
    BSOFT_LAYERLINE_FOM,
    BSOFT_LAYERLINE_SELECT,
    BSOFT_SYMMETRY_INT_TABLES_NUMBER,
    BSOFT_SYMMETRY_SPACE_GROUP_NAME_H_M,
    BSOFT_SYMMETRY_CELL_SETTING,
    BSOFT_SYMMETRY_EQUIV_ID,
    BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ,

    MDL_LAST_LABEL  // **** NOTE ****: Do keep this label always at the end,it is here for looping purposes
};//close enum Label

typedef std::vector<MDLabel> MDLabelVector;

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
    bool failed; // Set to True if the parsing from Star files fails
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
    //don't do any type checking as have not value yet
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
        if (withFormat) os << std::setw(20); \
        os << i;
        //this must have 20 since SIZE_MAX = 18446744073709551615 size

    void toStream(std::ostream &os, bool withFormat = false, bool isSql=false, bool escape=true) const;
    String toString(bool withFormat = false, bool isSql=false) const;
    bool fromStream(std::istream &is, bool fromString=false);
    friend std::istream& operator>> (std::istream& is, MDObject &value);
    friend std::ostream& operator<< (std::ostream& is, const MDObject &value);
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

    void setValueFromStr(MDLabel label, const String &value);

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
    /** Same as label2Str but escaping with '' to use in Sqlite. */
    static String label2StrSql(const MDLabel label);
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
    static MDLabelData * data[MDL_LAST_LABEL+1];
    static std::map<std::string, MDLabel> names;
    static MDLabelStaticInit initialization; //Just for initialization

    /** Add predefined labels to be used in metadata */
    static void addLabel(MDLabel label, MDLabelType type, const String &name, int tags=TAGLABEL_NOTAG);
    /** This function will read extra label alias from XMIPP_LABEL_ALIASES var environment */
    static void addExtraAliases();
    friend class MDLabelStaticInit;
public:
    /** Add an alias for an existing label,
     * If replace=true, then new alias name will be used
     * as the label string for writing back to file */
    static void addLabelAlias(MDLabel label, const String &alias, bool replace=false);
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
        //The label MDL_GATHER_ID is special and should not be used
        MDL::addLabel(MDL_GATHER_ID, LABEL_SIZET, "gatherId");

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
        MDL::addLabel(MDL_ANGLE_DIFF0, LABEL_DOUBLE, "angleDiff0");
        MDL::addLabel(MDL_ANGLE_DIFF, LABEL_DOUBLE, "angleDiff");
        MDL::addLabel(MDL_ANGLE_DIFF2, LABEL_DOUBLE, "angleDiff2");
        MDL::addLabel(MDL_ANGLE_Y, LABEL_DOUBLE, "angleY");
        MDL::addLabel(MDL_ANGLE_Y2, LABEL_DOUBLE, "angleY2");

        MDL::addLabel(MDL_APPLY_SHIFT, LABEL_BOOL, "applyShift");
        MDL::addLabel(MDL_AVG, LABEL_DOUBLE, "avg");
        MDL::addLabel(MDL_AVG_CHANGES_ORIENTATIONS, LABEL_DOUBLE, "avgChanOrient");
        MDL::addLabel(MDL_AVG_CHANGES_OFFSETS, LABEL_DOUBLE, "avgChanOffset");
        MDL::addLabel(MDL_AVG_CHANGES_CLASSES, LABEL_DOUBLE, "avgChanClass");
        MDL::addLabel(MDL_AVGPMAX, LABEL_DOUBLE, "avgPMax");

        MDL::addLabel(MDL_BGMEAN, LABEL_DOUBLE, "bgMean");
        MDL::addLabel(MDL_BLOCK_NUMBER, LABEL_INT, "blockNumber");

        MDL::addLabel(MDL_CL2D_CHANGES, LABEL_INT, "cl2dChanges");
        MDL::addLabel(MDL_CL2D_SIMILARITY, LABEL_DOUBLE, "cl2dSimilarity");
        MDL::addLabel(MDL_CLASS_COUNT, LABEL_SIZET, "classCount");
        MDL::addLabelAlias(MDL_CLASS_COUNT, "class_count"); //3.0
        MDL::addLabel(MDL_CLASS_PERCENTAGE, LABEL_DOUBLE, "classPercentage");
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
        MDL::addLabel(MDL_COLOR, LABEL_INT, "color");
        MDL::addLabel(MDL_COMMENT, LABEL_STRING, "comment");
        MDL::addLabel(MDL_COST, LABEL_DOUBLE, "cost");
        MDL::addLabel(MDL_COST_PERCENTILE, LABEL_DOUBLE, "costPerc");
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
        MDL::addLabel(MDL_CRYSTAL_NOISE_SHIFT,LABEL_VECTOR_DOUBLE, "crystalNoiseShift");
        MDL::addLabel(MDL_CTF_BG_BASELINE, LABEL_DOUBLE, "ctfBgBaseline");
        MDL::addLabelAlias(MDL_CTF_BG_BASELINE, "CTFBG_Baseline");//3.0
        MDL::addLabel(MDL_CTF_BG_GAUSSIAN2_ANGLE, LABEL_DOUBLE, "ctfBgGaussian2Angle");
        MDL::addLabelAlias(MDL_CTF_BG_GAUSSIAN2_ANGLE, "CTFBG_Gaussian2_Angle"); //3.0
        MDL::addLabel(MDL_CTF_BG_R1, LABEL_DOUBLE, "ctfBgR1");
        MDL::addLabel(MDL_CTF_BG_R2, LABEL_DOUBLE, "ctfBgR2");
        MDL::addLabel(MDL_CTF_BG_R3, LABEL_DOUBLE, "ctfBgR3");

        MDL::addLabel(MDL_CTF_DATA_PHASE_FLIPPED, LABEL_BOOL, "ctfPhaseFlipped");
        MDL::addLabel(MDL_CTF_CORRECTED, LABEL_BOOL, "ctfCorrected");
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

        MDL::addLabel(MDL_CONTINUOUS_X, LABEL_DOUBLE, "continuousX");
        MDL::addLabel(MDL_CONTINUOUS_Y, LABEL_DOUBLE, "continuousY");
        MDL::addLabel(MDL_CONTINUOUS_FLIP, LABEL_BOOL, "continuousFlip");
        MDL::addLabel(MDL_CONTINUOUS_GRAY_A, LABEL_DOUBLE, "continuousA");
        MDL::addLabel(MDL_CONTINUOUS_GRAY_B, LABEL_DOUBLE, "continuousB");
        MDL::addLabel(MDL_CONTINUOUS_SCALE_X, LABEL_DOUBLE, "continuousScaleX");
        MDL::addLabel(MDL_CONTINUOUS_SCALE_Y, LABEL_DOUBLE, "continuousScaleY");
        MDL::addLabel(MDL_CTF_CA, LABEL_DOUBLE, "ctfChromaticAberration");
        MDL::addLabel(MDL_CTF_CONVERGENCE_CONE, LABEL_DOUBLE, "ctfConvergenceCone");
        MDL::addLabel(MDL_CTF_CRIT_NONASTIGMATICVALIDITY, LABEL_DOUBLE, "ctfCritNonAstigmaticValidty");
        MDL::addLabel(MDL_CTF_CRIT_DAMPING, LABEL_DOUBLE, "ctfCritDamping");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZEROAVG, LABEL_DOUBLE, "ctfCritFirstZero");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZERODISAGREEMENT, LABEL_DOUBLE, "ctfCritDisagree");
        MDL::addLabel(MDL_CTF_CRIT_MAXFREQ, LABEL_DOUBLE, "ctfCritMaxFreq");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTZERORATIO, LABEL_DOUBLE, "ctfCritfirstZeroRatio");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_RATIO, LABEL_DOUBLE, "ctfCritFirstMinFirstZeroRatio");
        MDL::addLabel(MDL_CTF_CRIT_FIRSTMINIMUM_FIRSTZERO_DIFF_RATIO, LABEL_DOUBLE, "ctfCritCtfMargin");
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
        MDL::addLabel(MDL_CTF_ENVELOPE, LABEL_DOUBLE, "ctfEnvelope");
        MDL::addLabel(MDL_CTF_ENVELOPE_PLOT, LABEL_STRING, "ctfEnvelopePlot");
        MDL::addLabel(MDL_CTF_GROUP, LABEL_INT, "ctfGroup");
        MDL::addLabel(MDL_CTF_INPUTPARAMS, LABEL_STRING, "ctfInputParams", TAGLABEL_TEXTFILE);
        MDL::addLabel(MDL_CTF_K, LABEL_DOUBLE, "ctfK");
        MDL::addLabel(MDL_CTF_ENV_R0, LABEL_DOUBLE, "ctfEnvR0");
        MDL::addLabel(MDL_CTF_ENV_R1, LABEL_DOUBLE, "ctfEnvR1");
        MDL::addLabel(MDL_CTF_ENV_R2, LABEL_DOUBLE, "ctfEnvR2");
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
        MDL::addLabel(MDL_CUMULATIVE_SSNR, LABEL_DOUBLE, "cumulativeSSNR");

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
        MDL::addLabel(MDL_FRAME_ID, LABEL_SIZET, "frameId");

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
        MDL::addLabel(MDL_IMAGE_COVARIANCE, LABEL_STRING, "imageCovariance", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_IDX, LABEL_SIZET, "imageIndex");
        MDL::addLabel(MDL_IMAGE_ORIGINAL, LABEL_STRING, "imageOriginal", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_REF, LABEL_STRING, "imageRef", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_IMAGE_RESIDUAL, LABEL_STRING, "imageResidual", TAGLABEL_IMAGE);
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
        MDL::addLabel(MDL_MACRO_CMD, LABEL_STRING, "macroCmd");
        MDL::addLabel(MDL_MACRO_CMD_ARGS, LABEL_STRING, "macroCmdArgs");
        MDL::addLabel(MDL_MAGNIFICATION, LABEL_DOUBLE, "magnification");
        MDL::addLabel(MDL_MAPTOPOLOGY, LABEL_STRING, "mapTopology");
        MDL::addLabel(MDL_MASK, LABEL_STRING, "mask", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_MAXCC, LABEL_DOUBLE, "maxCC");
        MDL::addLabel(MDL_MAXCC_PERCENTILE, LABEL_DOUBLE, "maxCCPerc");
        MDL::addLabel(MDL_MAX, LABEL_DOUBLE, "max");
        MDL::addLabel(MDL_MICROGRAPH_ID, LABEL_SIZET, "micrographId");
        MDL::addLabel(MDL_MICROGRAPH, LABEL_STRING, "micrograph", TAGLABEL_MICROGRAPH);
        MDL::addLabel(MDL_MICROGRAPH_MOVIE_ID, LABEL_SIZET, "micrographMovieId");
        MDL::addLabel(MDL_MICROGRAPH_MOVIE, LABEL_STRING, "movie", TAGLABEL_MICROGRAPH);
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
        MDL::addLabel(MDL_NMA_ATOMSHIFT, LABEL_DOUBLE, "nmaAtomShift");
        MDL::addLabel(MDL_NMA_COLLECTIVITY, LABEL_DOUBLE, "nmaCollectivity");
        MDL::addLabel(MDL_NMA_ENERGY, LABEL_DOUBLE, "nmaEnergy");
        MDL::addLabel(MDL_NMA_MINRANGE, LABEL_DOUBLE, "nmaMin");
        MDL::addLabel(MDL_NMA_MAXRANGE, LABEL_DOUBLE, "nmaMax");
        MDL::addLabel(MDL_NMA_MODEFILE, LABEL_STRING, "nmaModefile", TAGLABEL_TEXTFILE);
        MDL::addLabelAlias(MDL_NMA_MODEFILE, "NMAModefile");//3.0
        MDL::addLabel(MDL_NMA_SCORE, LABEL_DOUBLE, "nmaScore");
        MDL::addLabel(MDL_NOISE_ANGLES, LABEL_VECTOR_DOUBLE, "noiseAngles");
        MDL::addLabel(MDL_NOISE_COORD, LABEL_VECTOR_DOUBLE, "noiseCoord");
        MDL::addLabel(MDL_NOISE_PARTICLE_COORD, LABEL_VECTOR_DOUBLE, "noiseParticleCoord");
        MDL::addLabel(MDL_NOISE_PIXEL_LEVEL, LABEL_VECTOR_DOUBLE, "noisePixelLevel");
        MDL::addLabel(MDL_ORDER, LABEL_SIZET, "order_");
        MDL::addLabel(MDL_ORIGIN_X, LABEL_DOUBLE, "originX");
        MDL::addLabel(MDL_ORIGIN_Y, LABEL_DOUBLE, "originY");
        MDL::addLabel(MDL_ORIGIN_Z, LABEL_DOUBLE, "originZ");
        MDL::addLabel(MDL_PARTICLE_ID, LABEL_SIZET, "particleId");
        MDL::addLabel(MDL_PHANTOM_BGDENSITY, LABEL_DOUBLE, "phantomBGDensity");
        MDL::addLabel(MDL_PHANTOM_FEATURE_CENTER, LABEL_VECTOR_DOUBLE, "featureCenter");
        MDL::addLabel(MDL_PHANTOM_FEATURE_DENSITY, LABEL_DOUBLE, "featureDensity");
        MDL::addLabel(MDL_PHANTOM_FEATURE_OPERATION, LABEL_STRING, "featureOperation");
        MDL::addLabel(MDL_PHANTOM_FEATURE_SPECIFIC, LABEL_VECTOR_DOUBLE, "featureSpecificVector");
        MDL::addLabel(MDL_PHANTOM_FEATURE_TYPE, LABEL_STRING, "featureType");
        MDL::addLabel(MDL_PHANTOM_SCALE, LABEL_DOUBLE, "phantomScale");

        MDL::addLabel(MDL_OPTICALFLOW_MEANX, LABEL_DOUBLE, "opticalMeanX");
        MDL::addLabel(MDL_OPTICALFLOW_MEANY, LABEL_DOUBLE, "opticalMeanY");
        MDL::addLabel(MDL_OPTICALFLOW_STDX, LABEL_DOUBLE, "opticalStdX");
        MDL::addLabel(MDL_OPTICALFLOW_STDY, LABEL_DOUBLE, "opticalStdY");

        MDL::addLabel(MDL_PICKING_STATE, LABEL_STRING, "pickingState");
        MDL::addLabelAlias(MDL_PICKING_STATE, "picking_state");//3.0
        MDL::addLabel(MDL_PICKING_MICROGRAPH_STATE, LABEL_STRING, "pickingMicrographState");
        MDL::addLabelAlias(MDL_PICKING_MICROGRAPH_STATE, "micrograph_state");//3.0
        MDL::addLabel(MDL_PICKING_PARTICLE_SIZE, LABEL_INT, "particleSize");
        MDL::addLabel(MDL_PICKING_AUTOPICKPERCENT, LABEL_INT, "autopickPercent");
        MDL::addLabel(MDL_PICKING_TEMPLATES, LABEL_INT, "templatesNum");
        MDL::addLabel(MDL_PICKING_AUTOPARTICLES_SIZE, LABEL_INT, "autoParticlesNum");
        MDL::addLabel(MDL_PICKING_MANUALPARTICLES_SIZE, LABEL_INT, "manualParticlesNum");

        MDL::addLabel(MDL_PMAX, LABEL_DOUBLE, "pMax");
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
        MDL::addLabel(MDL_RESOLUTION_FREQ2, LABEL_DOUBLE, "resolutionFreqFourier2");
        MDL::addLabel(MDL_RESOLUTION_FREQREAL, LABEL_DOUBLE, "resolutionFreqReal");
        MDL::addLabel(MDL_RESOLUTION_LOG_STRUCTURE_FACTOR, LABEL_DOUBLE, "resolutionLogStructure");
        MDL::addLabel(MDL_RESOLUTION_STRUCTURE_FACTOR, LABEL_DOUBLE, "resolutionStructure");
        MDL::addLabel(MDL_RESOLUTION_SSNR, LABEL_DOUBLE, "resolutionSSNR");
        MDL::addLabel(MDL_RESOLUTION_RFACTOR, LABEL_DOUBLE, "resolutionRfactor");

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
        MDL::addLabel(MDL_SCORE_BY_PCA_RESIDUAL_PROJ, LABEL_DOUBLE, "scoreByPcaResidualProj");
        MDL::addLabel(MDL_SCORE_BY_PCA_RESIDUAL_EXP, LABEL_DOUBLE, "scoreByPcaResidualExp");
        MDL::addLabel(MDL_SCORE_BY_PCA_RESIDUAL, LABEL_DOUBLE, "scoreByPcaResidual");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY, LABEL_DOUBLE, "scoreByAlignability");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_PRECISION, LABEL_DOUBLE, "scoreByAlignabilityPrecision");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_ACCURACY, LABEL_DOUBLE, "scoreByAlignabilityAccuracy");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_PRECISION_EXP, LABEL_DOUBLE, "scoreByAlignabilityPrecisionExp");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_PRECISION_REF, LABEL_DOUBLE, "scoreByAlignabilityPrecisionRef");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_ACCURACY_EXP, LABEL_DOUBLE, "scoreByAlignabilityAccuracyExp");
        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_ACCURACY_REF, LABEL_DOUBLE, "scoreByAlignabilityAccuracyRef");

        MDL::addLabel(MDL_SCORE_BY_ALIGNABILITY_NOISE, LABEL_DOUBLE, "scoreByAlignabilityNoise");
        MDL::addLabel(MDL_SCORE_BY_MIRROR, LABEL_DOUBLE, "scoreByMirror");
        MDL::addLabel(MDL_SCORE_BY_ZSCORE, LABEL_DOUBLE, "scoreByZScore");

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
        MDL::addLabel(MDL_SHIFT_DIFF0, LABEL_DOUBLE, "shiftDiff0");
        MDL::addLabel(MDL_SHIFT_DIFF, LABEL_DOUBLE, "shiftDiff");
        MDL::addLabel(MDL_SHIFT_DIFF2, LABEL_DOUBLE, "shiftDiff2");
        MDL::addLabel(MDL_SIGMANOISE, LABEL_DOUBLE, "sigmaNoise");
        MDL::addLabel(MDL_SIGMAOFFSET, LABEL_DOUBLE, "sigmaOffset");
        MDL::addLabel(MDL_SIGNALCHANGE, LABEL_DOUBLE, "signalChange");
        MDL::addLabel(MDL_STDDEV, LABEL_DOUBLE, "stddev");
        MDL::addLabel(MDL_STAR_COMMENT, LABEL_STRING, "starComment");
        MDL::addLabel(MDL_SUM, LABEL_DOUBLE, "sum");
        MDL::addLabel(MDL_SUMWEIGHT, LABEL_DOUBLE, "sumWeight");
        MDL::addLabel(MDL_SYMNO, LABEL_INT, "symNo");

        MDL::addLabel(MDL_TOMOGRAM_VOLUME, LABEL_STRING, "tomogramVolume", TAGLABEL_IMAGE);
        MDL::addLabel(MDL_TOMOGRAMMD, LABEL_STRING, "tomogramMetadata", TAGLABEL_METADATA);

        MDL::addLabel(MDL_TRANSFORM_MATRIX, LABEL_STRING, "transformMatrix");

        MDL::addLabel(MDL_TEST_SIZE, LABEL_INT, "testSize");

        MDL::addLabel(MDL_VOLUME_SCORE_SUM, LABEL_DOUBLE, "volScoreSum");
        MDL::addLabel(MDL_VOLUME_SCORE_MEAN, LABEL_DOUBLE, "volScoreMean");
        MDL::addLabel(MDL_VOLUME_SCORE_MIN, LABEL_DOUBLE, "volScoreMin");
        MDL::addLabel(MDL_VOLUME_SCORE_SUM_TH, LABEL_DOUBLE, "volScoreSumTh");
        MDL::addLabel(MDL_VOLUME_SCORE1, LABEL_DOUBLE, "volScore1");
        MDL::addLabel(MDL_VOLUME_SCORE2, LABEL_DOUBLE, "volScore2");
        MDL::addLabel(MDL_VOLUME_SCORE3, LABEL_DOUBLE, "volScore3");
        MDL::addLabel(MDL_VOLUME_SCORE4, LABEL_DOUBLE, "volScore4");
        MDL::addLabel(MDL_WEIGHT, LABEL_DOUBLE, "weight");
        MDL::addLabel(MDL_WEIGHT_P, LABEL_DOUBLE, "weight_clusterability");
        MDL::addLabelAlias(MDL_WEIGHT, "Weight");
        MDL::addLabel(MDL_WEIGHT_CONTINUOUS2, LABEL_DOUBLE, "weightContinuous2");
        MDL::addLabel(MDL_WEIGHT_JUMPER0, LABEL_DOUBLE, "weightJumper0");
        MDL::addLabel(MDL_WEIGHT_JUMPER, LABEL_DOUBLE, "weightJumper");
        MDL::addLabel(MDL_WEIGHT_JUMPER2, LABEL_DOUBLE, "weightJumper2");
        MDL::addLabel(MDL_WEIGHT_SIGNIFICANT, LABEL_DOUBLE, "weightSignificant");
        MDL::addLabel(MDL_WEIGHT_SSNR, LABEL_DOUBLE, "weightSSNR");

        MDL::addLabel(MDL_WEIGHT_PRECISION_ALIGNABILITY, LABEL_DOUBLE, "weightPrecisionAlignability");
        MDL::addLabel(MDL_WEIGHT_ACCURACY_ALIGNABILITY, LABEL_DOUBLE, "weightAccuracyAlignability");
        MDL::addLabel(MDL_WEIGHT_ALIGNABILITY, LABEL_DOUBLE, "weightAlignability");
        MDL::addLabel(MDL_WEIGHT_PRECISION_MIRROR, LABEL_DOUBLE, "weightPrecisionMirror");

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
        MDL::addLabel(MDL_ZSCORE_HISTOGRAM, LABEL_DOUBLE, "zScoreHistogram");
        MDL::addLabel(MDL_ZSCORE_RESMEAN, LABEL_DOUBLE, "zScoreResMean");
        MDL::addLabel(MDL_ZSCORE_RESVAR, LABEL_DOUBLE, "zScoreResVar");
        MDL::addLabel(MDL_ZSCORE_RESCOV, LABEL_DOUBLE, "zScoreResCov");
        MDL::addLabel(MDL_ZSCORE_SHAPE1, LABEL_DOUBLE, "zScoreShape1");
        MDL::addLabel(MDL_ZSCORE_SHAPE2, LABEL_DOUBLE, "zScoreShape2");
        MDL::addLabel(MDL_ZSCORE_SNR1, LABEL_DOUBLE, "zScoreSNR1");
        MDL::addLabel(MDL_ZSCORE_SNR2, LABEL_DOUBLE, "zScoreSNR2");
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


        /*Relion labels */
        MDL::addLabel(RLN_AREA_ID, LABEL_SIZET, "rlnAreaId");
        MDL::addLabel(RLN_AREA_NAME, LABEL_STRING, "rlnAreaName");

        MDL::addLabel(RLN_CTF_BFACTOR, LABEL_DOUBLE, "rlnBfactor");
        MDL::addLabel(RLN_CTF_SCALEFACTOR, LABEL_DOUBLE, "rlnCtfScalefactor");
        MDL::addLabel(RLN_CTF_VOLTAGE, LABEL_DOUBLE, "rlnVoltage");
        MDL::addLabel(RLN_CTF_DEFOCUSU, LABEL_DOUBLE, "rlnDefocusU");
        MDL::addLabel(RLN_CTF_DEFOCUSV, LABEL_DOUBLE, "rlnDefocusV");
        MDL::addLabel(RLN_CTF_DEFOCUS_ANGLE, LABEL_DOUBLE, "rlnDefocusAngle");
        MDL::addLabel(RLN_CTF_CS, LABEL_DOUBLE, "rlnSphericalAberration");
        MDL::addLabel(RLN_CTF_CA, LABEL_DOUBLE, "rlnChromaticAberration");
        MDL::addLabel(RLN_CTF_DETECTOR_PIXEL_SIZE, LABEL_DOUBLE, "rlnDetectorPixelSize");
        MDL::addLabel(RLN_CTF_ENERGY_LOSS, LABEL_DOUBLE, "rlnEnergyLoss");
        MDL::addLabel(RLN_CTF_FOM, LABEL_DOUBLE, "rlnCtfFigureOfMerit");
        MDL::addLabel(RLN_CTF_IMAGE, LABEL_STRING, "rlnCtfImage", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_CTF_LENS_STABILITY, LABEL_DOUBLE, "rlnLensStability");
        MDL::addLabel(RLN_CTF_MAGNIFICATION, LABEL_DOUBLE, "rlnMagnification");
        MDL::addLabel(RLN_CTF_CONVERGENCE_CONE, LABEL_DOUBLE, "rlnConvergenceCone");
        MDL::addLabel(RLN_CTF_LONGITUDINAL_DISPLACEMENT, LABEL_DOUBLE, "rlnLongitudinalDisplacement");
        MDL::addLabel(RLN_CTF_TRANSVERSAL_DISPLACEMENT, LABEL_DOUBLE, "rlnTransversalDisplacement");
        MDL::addLabel(RLN_CTF_Q0, LABEL_DOUBLE, "rlnAmplitudeContrast");
        MDL::addLabel(RLN_CTF_VALUE, LABEL_DOUBLE, "rlnCtfValue");

        MDL::addLabel(RLN_IMAGE_NAME, LABEL_STRING, "rlnImageName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_IMAGE_RECONSTRUCT_NAME, LABEL_STRING, "rlnReconstructImageName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_IMAGE_ID, LABEL_SIZET, "rlnImageId");
        MDL::addLabel(RLN_IMAGE_ENABLED, LABEL_BOOL, "rlnEnabled");
        MDL::addLabel(RLN_IMAGE_DATATYPE, LABEL_INT, "rlnDataType");
        // Label rlnImageDimensionality is originally rlnDataDimensionality, which is
        // duplicated for other label. A relion bug???
        MDL::addLabel(RLN_IMAGE_DIMENSIONALITY, LABEL_INT, "rlnImageDimensionality");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_X, LABEL_DOUBLE, "rlnBeamTiltX");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_Y, LABEL_DOUBLE, "rlnBeamTiltY");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_GROUP, LABEL_STRING, "rlnBeamTiltGroupName");
        MDL::addLabel(RLN_IMAGE_COORD_X, LABEL_DOUBLE, "rlnCoordinateX");
        MDL::addLabel(RLN_IMAGE_COORD_Y, LABEL_DOUBLE, "rlnCoordinateY");
        MDL::addLabel(RLN_IMAGE_COORD_Z, LABEL_DOUBLE, "rlnCoordinateZ");
        MDL::addLabel(RLN_IMAGE_FRAME_NR, LABEL_INT, "rlnMovieFrameNumber");
        MDL::addLabel(RLN_IMAGE_MAGNIFICATION_CORRECTION, LABEL_DOUBLE, "rlnMagnificationCorrection");
        MDL::addLabel(RLN_IMAGE_NORM_CORRECTION, LABEL_DOUBLE, "rlnNormCorrection");
        MDL::addLabel(RLN_IMAGE_ORI_NAME, LABEL_STRING, "rlnImageOriginalName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE, LABEL_DOUBLE, "rlnSamplingRate");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_X, LABEL_DOUBLE, "rlnSamplingRateX");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_Y, LABEL_DOUBLE, "rlnSamplingRateY");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_Z, LABEL_DOUBLE, "rlnSamplingRateZ");
        MDL::addLabel(RLN_IMAGE_SIZE, LABEL_INT, "rlnImageSize");
        MDL::addLabel(RLN_IMAGE_SIZEX, LABEL_INT, "rlnImageSizeX");
        MDL::addLabel(RLN_IMAGE_SIZEY, LABEL_INT, "rlnImageSizeY");
        MDL::addLabel(RLN_IMAGE_SIZEZ, LABEL_INT, "rlnImageSizeZ");
        MDL::addLabel(RLN_IMAGE_STATS_MIN, LABEL_DOUBLE, "rlnMinimumValue");
        MDL::addLabel(RLN_IMAGE_STATS_MAX, LABEL_DOUBLE, "rlnMaximumValue");
        MDL::addLabel(RLN_IMAGE_STATS_AVG, LABEL_DOUBLE, "rlnAverageValue");
        MDL::addLabel(RLN_IMAGE_STATS_STDDEV, LABEL_DOUBLE, "rlnStandardDeviationValue");
        MDL::addLabel(RLN_IMAGE_STATS_SKEW, LABEL_DOUBLE, "rlnSkewnessValue");
        MDL::addLabel(RLN_IMAGE_STATS_KURT, LABEL_DOUBLE, "rlnKurtosisExcessValue");
        MDL::addLabel(RLN_IMAGE_WEIGHT, LABEL_DOUBLE, "rlnImageWeight");

        MDL::addLabel(RLN_MATRIX_1_1, LABEL_DOUBLE, "rlnMatrix_1_1");
        MDL::addLabel(RLN_MATRIX_1_2, LABEL_DOUBLE, "rlnMatrix_1_2");
        MDL::addLabel(RLN_MATRIX_1_3, LABEL_DOUBLE, "rlnMatrix_1_3");
        MDL::addLabel(RLN_MATRIX_2_1, LABEL_DOUBLE, "rlnMatrix_2_1");
        MDL::addLabel(RLN_MATRIX_2_2, LABEL_DOUBLE, "rlnMatrix_2_2");
        MDL::addLabel(RLN_MATRIX_2_3, LABEL_DOUBLE, "rlnMatrix_2_3");
        MDL::addLabel(RLN_MATRIX_3_1, LABEL_DOUBLE, "rlnMatrix_3_1");
        MDL::addLabel(RLN_MATRIX_3_2, LABEL_DOUBLE, "rlnMatrix_3_2");
        MDL::addLabel(RLN_MATRIX_3_3, LABEL_DOUBLE, "rlnMatrix_3_3");

        MDL::addLabel(RLN_MICROGRAPH_ID, LABEL_SIZET, "rlnMicrographId");
        MDL::addLabel(RLN_MICROGRAPH_NAME, LABEL_STRING, "rlnMicrographName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_MICROGRAPH_TILT_ANGLE, LABEL_DOUBLE, "rlnMicrographTiltAngle");
        MDL::addLabel(RLN_MICROGRAPH_TILT_AXIS_DIRECTION, LABEL_DOUBLE, "rlnMicrographTiltAxisDirection");
        MDL::addLabel(RLN_MICROGRAPH_TILT_AXIS_OUTOFPLANE, LABEL_DOUBLE, "rlnMicrographTiltAxisOutOfPlane");

        MDL::addLabel(RLN_MLMODEL_ACCURACY_ROT, LABEL_DOUBLE, "rlnAccuracyRotations");
        MDL::addLabel(RLN_MLMODEL_ACCURACY_TRANS, LABEL_DOUBLE, "rlnAccuracyTranslations");
        MDL::addLabel(RLN_MLMODEL_AVE_PMAX, LABEL_DOUBLE, "rlnAveragePmax");
        MDL::addLabel(RLN_MLMODEL_CURRENT_RESOLUTION, LABEL_DOUBLE, "rlnCurrentResolution");
        MDL::addLabel(RLN_MLMODEL_CURRENT_SIZE, LABEL_INT, "rlnCurrentImageSize");
        MDL::addLabel(RLN_MLMODEL_DATA_VS_PRIOR_REF, LABEL_DOUBLE, "rlnSsnrMap");
        MDL::addLabel(RLN_MLMODEL_DIMENSIONALITY, LABEL_INT, "rlnReferenceDimensionality");
        MDL::addLabel(RLN_MLMODEL_DIMENSIONALITY_DATA, LABEL_INT, "rlnDataDimensionality");
        MDL::addLabel(RLN_MLMODEL_DIFF2_HALVES_REF, LABEL_DOUBLE, "rlnDiff2RandomHalves");
        MDL::addLabel(RLN_MLMODEL_FSC_HALVES_REF, LABEL_DOUBLE, "rlnGoldStandardFsc");
        MDL::addLabel(RLN_MLMODEL_GROUP_NAME, LABEL_STRING, "rlnGroupName");
        MDL::addLabel(RLN_MLMODEL_GROUP_NO, LABEL_SIZET, "rlnGroupNumber");
        MDL::addLabel(RLN_MLMODEL_GROUP_NR_PARTICLES, LABEL_SIZET, "rlnGroupNrParticles");
        MDL::addLabel(RLN_MLMODEL_GROUP_SCALE_CORRECTION, LABEL_DOUBLE, "rlnGroupScaleCorrection");
        MDL::addLabel(RLN_MLMODEL_INTERPOLATOR, LABEL_INT, "rlnFourierSpaceInterpolator");
        MDL::addLabel(RLN_MLMODEL_LL, LABEL_DOUBLE, "rlnLogLikelihood");
        MDL::addLabel(RLN_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION,  LABEL_INT, "rlnMinRadiusNnInterpolation");
        MDL::addLabel(RLN_MLMODEL_NORM_CORRECTION_AVG, LABEL_DOUBLE, "rlnNormCorrectionAverage");
        MDL::addLabel(RLN_MLMODEL_NR_CLASSES, LABEL_INT, "rlnNrClasses");
        MDL::addLabel(RLN_MLMODEL_NR_GROUPS, LABEL_INT, "rlnNrGroups");
        MDL::addLabel(RLN_MLMODEL_ORIENTABILITY_CONTRIBUTION, LABEL_DOUBLE, "rlnSpectralOrientabilityContribution");
        MDL::addLabel(RLN_MLMODEL_ORIGINAL_SIZE, LABEL_INT, "rlnOriginalImageSize");
        MDL::addLabel(RLN_MLMODEL_PADDING_FACTOR, LABEL_INT, "rlnPaddingFactor");
        MDL::addLabel(RLN_MLMODEL_PDF_CLASS, LABEL_DOUBLE, "rlnClassDistribution");
        MDL::addLabel(RLN_MLMODEL_PRIOR_OFFX_CLASS, LABEL_DOUBLE, "rlnClassPriorOffsetX");
        MDL::addLabel(RLN_MLMODEL_PRIOR_OFFY_CLASS, LABEL_DOUBLE, "rlnClassPriorOffsetY");
        MDL::addLabel(RLN_MLMODEL_PDF_ORIENT, LABEL_DOUBLE, "rlnOrientationDistribution");
        MDL::addLabel(RLN_MLMODEL_PIXEL_SIZE, LABEL_DOUBLE, "rlnPixelSize");
        MDL::addLabel(RLN_MLMODEL_POWER_REF, LABEL_DOUBLE, "rlnReferenceSpectralPower");
        MDL::addLabel(RLN_MLMODEL_PRIOR_MODE, LABEL_INT, "rlnOrientationalPriorMode");
        MDL::addLabel(RLN_MLMODEL_REF_IMAGE, LABEL_STRING, "rlnReferenceImage", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_MLMODEL_SIGMA_OFFSET, LABEL_DOUBLE, "rlnSigmaOffsets");
        MDL::addLabel(RLN_MLMODEL_SIGMA2_NOISE, LABEL_DOUBLE, "rlnSigma2Noise");
        MDL::addLabel(RLN_MLMODEL_SIGMA2_REF, LABEL_DOUBLE, "rlnReferenceSigma2");
        MDL::addLabel(RLN_MLMODEL_SIGMA_ROT, LABEL_DOUBLE, "rlnSigmaPriorRotAngle");
        MDL::addLabel(RLN_MLMODEL_SIGMA_TILT, LABEL_DOUBLE, "rlnSigmaPriorTiltAngle");
        MDL::addLabel(RLN_MLMODEL_SIGMA_PSI, LABEL_DOUBLE, "rlnSigmaPriorPsiAngle");
        MDL::addLabel(RLN_MLMODEL_SSNR_REF, LABEL_DOUBLE, "rlnSignalToNoiseRatio");
        MDL::addLabel(RLN_MLMODEL_TAU2_FUDGE_FACTOR, LABEL_DOUBLE, "rlnTau2FudgeFactor");
        MDL::addLabel(RLN_MLMODEL_TAU2_REF, LABEL_DOUBLE, "rlnReferenceTau2");

        MDL::addLabel(RLN_OPTIMISER_ACCURACY_ROT, LABEL_DOUBLE, "rlnOverallAccuracyRotations");
        MDL::addLabel(RLN_OPTIMISER_ACCURACY_TRANS, LABEL_DOUBLE, "rlnOverallAccuracyTranslations");
        MDL::addLabel(RLN_OPTIMISER_ADAPTIVE_FRACTION, LABEL_DOUBLE, "rlnAdaptiveOversampleFraction");
        MDL::addLabel(RLN_OPTIMISER_ADAPTIVE_OVERSAMPLING, LABEL_INT, "rlnAdaptiveOversampleOrder");
        MDL::addLabel(RLN_OPTIMISER_AUTO_LOCAL_HP_ORDER, LABEL_INT, "rlnAutoLocalSearchesHealpixOrder");
        MDL::addLabel(RLN_OPTIMISER_AVAILABLE_MEMORY, LABEL_DOUBLE, "rlnAvailableMemory");
        MDL::addLabel(RLN_OPTIMISER_BEST_RESOL_THUS_FAR, LABEL_DOUBLE, "rlnBestResolutionThusFar");
        MDL::addLabel(RLN_OPTIMISER_COARSE_SIZE, LABEL_INT, "rlnCoarseImageSize");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, LABEL_DOUBLE, "rlnChangesOptimalOffsets");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, LABEL_DOUBLE, "rlnChangesOptimalOrientations");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES, LABEL_DOUBLE, "rlnChangesOptimalClasses");
        MDL::addLabel(RLN_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, LABEL_BOOL, "rlnCtfDataArePhaseFlipped");
        MDL::addLabel(RLN_OPTIMISER_DATA_STARFILE, LABEL_STRING, "rlnExperimentalDataStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_CTF, LABEL_BOOL, "rlnDoCorrectCtf");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_MAGNIFICATION, LABEL_BOOL, "rlnDoCorrectMagnification");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_NORM, LABEL_BOOL, "rlnDoCorrectNorm");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_SCALE, LABEL_BOOL, "rlnDoCorrectScale");
        MDL::addLabel(RLN_OPTIMISER_DO_REALIGN_MOVIES, LABEL_BOOL, "rlnDoRealignMovies");
        MDL::addLabel(RLN_OPTIMISER_DO_MAP, LABEL_BOOL, "rlnDoMapEstimation");
        MDL::addLabel(RLN_OPTIMISER_DO_AUTO_REFINE, LABEL_BOOL, "rlnDoAutoRefine");
        MDL::addLabel(RLN_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, LABEL_BOOL, "rlnDoOnlyFlipCtfPhases");
        MDL::addLabel(RLN_OPTIMISER_DO_SOLVENT_FLATTEN, LABEL_BOOL, "rlnDoSolventFlattening");
        MDL::addLabel(RLN_OPTIMISER_DO_SKIP_ALIGN, LABEL_BOOL, "rlnDoSkipAlign");
        MDL::addLabel(RLN_OPTIMISER_DO_SKIP_ROTATE, LABEL_BOOL, "rlnDoSkipRotate");
        MDL::addLabel(RLN_OPTIMISER_DO_SPLIT_RANDOM_HALVES, LABEL_BOOL, "rlnDoSplitRandomHalves");
        MDL::addLabel(RLN_OPTIMISER_DO_ZERO_MASK, LABEL_BOOL, "rlnDoZeroMask");
        MDL::addLabel(RLN_OPTIMISER_FIX_SIGMA_NOISE, LABEL_BOOL, "rlnFixSigmaNoiseEstimates");
        MDL::addLabel(RLN_OPTIMISER_FIX_SIGMA_OFFSET ,LABEL_BOOL, "rlnFixSigmaOffsetEstimates");
        MDL::addLabel(RLN_OPTIMISER_FIX_TAU, LABEL_BOOL, "rlnFixTauEstimates");
        MDL::addLabel(RLN_OPTIMISER_HAS_CONVERGED, LABEL_BOOL, "rlnHasConverged");
        MDL::addLabel(RLN_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, LABEL_BOOL, "rlnHasHighFscAtResolLimit");
        MDL::addLabel(RLN_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, LABEL_INT, "rlnHasLargeSizeIncreaseIterationsAgo");
        MDL::addLabel(RLN_OPTIMISER_HIGHRES_LIMIT_EXP, LABEL_DOUBLE, "rlnHighresLimitExpectation");
        MDL::addLabel(RLN_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, LABEL_BOOL, "rlnDoIgnoreCtfUntilFirstPeak");
        MDL::addLabel(RLN_OPTIMISER_INCR_SIZE, LABEL_INT, "rlnIncrementImageSize");
        MDL::addLabel(RLN_OPTIMISER_ITERATION_NO, LABEL_INT, "rlnCurrentIteration");
        MDL::addLabel(RLN_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, LABEL_DOUBLE, "rlnJoinHalvesUntilThisResolution");
        MDL::addLabel(RLN_OPTIMISER_MAGNIFICATION_RANGE, LABEL_DOUBLE, "rlnMagnificationSearchRange");
        MDL::addLabel(RLN_OPTIMISER_MAGNIFICATION_STEP, LABEL_DOUBLE, "rlnMagnificationSearchStep");
        MDL::addLabel(RLN_OPTIMISER_MAX_COARSE_SIZE, LABEL_INT, "rlnMaximumCoarseImageSize");
        MDL::addLabel(RLN_OPTIMISER_MAX_NR_POOL, LABEL_INT, "rlnMaxNumberOfPooledParticles");
        MDL::addLabel(RLN_OPTIMISER_MODEL_STARFILE, LABEL_STRING, "rlnModelStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_MODEL_STARFILE2, LABEL_STRING, "rlnModelStarFile2", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_NR_ITERATIONS, LABEL_INT, "rlnNumberOfIterations");
        MDL::addLabel(RLN_OPTIMISER_NR_ITER_WO_RESOL_GAIN, LABEL_INT, "rlnNumberOfIterWithoutResolutionGain");
        MDL::addLabel(RLN_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, LABEL_INT, "rlnNumberOfIterWithoutChangingAssignments");
        MDL::addLabel(RLN_OPTIMISER_OUTPUT_ROOTNAME, LABEL_STRING, "rlnOutputRootName");
        MDL::addLabel(RLN_OPTIMISER_PARTICLE_DIAMETER, LABEL_DOUBLE, "rlnParticleDiameter");
        MDL::addLabel(RLN_OPTIMISER_RADIUS_MASK_3D_MAP, LABEL_INT, "rlnRadiusMaskMap");
        MDL::addLabel(RLN_OPTIMISER_RADIUS_MASK_EXP_PARTICLES, LABEL_INT, "rlnRadiusMaskExpImages");
        MDL::addLabel(RLN_OPTIMISER_RANDOM_SEED, LABEL_INT, "rlnRandomSeed");
        MDL::addLabel(RLN_OPTIMISER_REFS_ARE_CTF_CORRECTED, LABEL_BOOL, "rlnRefsAreCtfCorrected");
        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, LABEL_INT, "rlnSmallestChangesClasses");
        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, LABEL_DOUBLE, "rlnSmallestChangesOffsets");
        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, LABEL_DOUBLE, "rlnSmallestChangesOrientations");
        MDL::addLabel(RLN_OPTIMISER_SAMPLING_STARFILE, LABEL_STRING, "rlnOrientSamplingStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_SOLVENT_MASK_NAME, LABEL_STRING, "rlnSolventMaskName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_OPTIMISER_SOLVENT_MASK2_NAME, LABEL_STRING, "rlnSolventMask2Name", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_OPTIMISER_TAU_SPECTRUM_NAME, LABEL_STRING, "rlnTauSpectrumName");
        MDL::addLabel(RLN_OPTIMISER_USE_TOO_COARSE_SAMPLING, LABEL_BOOL, "rlnUseTooCoarseSampling");
        MDL::addLabel(RLN_OPTIMISER_WIDTH_MASK_EDGE, LABEL_INT, "rlnWidthMaskEdge");

        MDL::addLabel(RLN_ORIENT_FLIP, LABEL_BOOL, "rlnIsFlip");
        MDL::addLabel(RLN_ORIENT_ID, LABEL_SIZET, "rlnOrientationsID");
        MDL::addLabel(RLN_ORIENT_ORIGIN_X, LABEL_DOUBLE, "rlnOriginX");
        MDL::addLabel(RLN_ORIENT_ORIGIN_X_PRIOR, LABEL_DOUBLE, "rlnOriginXPrior");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Y, LABEL_DOUBLE, "rlnOriginY");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Y_PRIOR, LABEL_DOUBLE, "rlnOriginYPrior");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Z, LABEL_DOUBLE, "rlnOriginZ");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Z_PRIOR, LABEL_DOUBLE, "rlnOriginZPrior");
        MDL::addLabel(RLN_ORIENT_ROT, LABEL_DOUBLE, "rlnAngleRot");
        MDL::addLabel(RLN_ORIENT_ROT_PRIOR, LABEL_DOUBLE, "rlnAngleRotPrior");
        MDL::addLabel(RLN_ORIENT_TILT, LABEL_DOUBLE, "rlnAngleTilt");
        MDL::addLabel(RLN_ORIENT_TILT_PRIOR, LABEL_DOUBLE, "rlnAngleTiltPrior");
        MDL::addLabel(RLN_ORIENT_PSI, LABEL_DOUBLE, "rlnAnglePsi");
        MDL::addLabel(RLN_ORIENT_PSI_PRIOR, LABEL_DOUBLE, "rlnAnglePsiPrior");

        MDL::addLabel(RLN_PARTICLE_AUTOPICK_FOM, LABEL_DOUBLE, "rlnAutopickFigureOfMerit");
        MDL::addLabel(RLN_PARTICLE_CLASS, LABEL_INT, "rlnClassNumber");
        MDL::addLabel(RLN_PARTICLE_DLL, LABEL_DOUBLE, "rlnLogLikeliContribution");
        MDL::addLabel(RLN_PARTICLE_ID, LABEL_SIZET, "rlnParticleId");
        MDL::addLabel(RLN_PARTICLE_FOM, LABEL_DOUBLE, "rlnParticleFigureOfMerit");
        MDL::addLabel(RLN_PARTICLE_KL_DIVERGENCE, LABEL_DOUBLE, "rlnKullbackLeibnerDivergence");
        MDL::addLabel(RLN_PARTICLE_RANDOM_SUBSET, LABEL_INT, "rlnRandomSubset");
        MDL::addLabel(RLN_PARTICLE_NAME, LABEL_STRING, "rlnParticleName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_PARTICLE_ORI_NAME, LABEL_STRING, "rlnOriginalParticleName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES, LABEL_INT, "rlnNrOfSignificantSamples"); /**< particle, Number of orientations contributing to weights*/
        MDL::addLabel(RLN_PARTICLE_NR_FRAMES, LABEL_INT, "rlnNrOfFrames");
        MDL::addLabel(RLN_PARTICLE_PMAX, LABEL_DOUBLE, "rlnMaxValueProbDistribution"); /**< particle, Maximum value of probability distribution */

        MDL::addLabel(RLN_POSTPROCESS_FINAL_RESOLUTION, LABEL_DOUBLE, "rlnFinalResolution");
        MDL::addLabel(RLN_POSTPROCESS_BFACTOR, LABEL_DOUBLE, "rlnBfactorUsedForSharpening");
        MDL::addLabel(RLN_POSTPROCESS_FSC_TRUE, LABEL_DOUBLE, "rlnFourierShellCorrelationCorrected");
        MDL::addLabel(RLN_POSTPROCESS_FSC_MASKED, LABEL_DOUBLE, "rlnFourierShellCorrelationMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_FSC_UNMASKED, LABEL_DOUBLE, "rlnFourierShellCorrelationUnmaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_FSC_RANDOM_MASKED, LABEL_DOUBLE, "rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT, LABEL_DOUBLE, "rlnFittedInterceptGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_SLOPE, LABEL_DOUBLE, "rlnFittedSlopeGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_CORRELATION, LABEL_DOUBLE, "rlnCorrelationFitGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_IN, LABEL_DOUBLE, "rlnLogAmplitudesOriginal");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_INVMTF, LABEL_DOUBLE, "rlnLogAmplitudesMTFCorrected");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED, LABEL_DOUBLE, "rlnLogAmplitudesWeighted");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED, LABEL_DOUBLE, "rlnLogAmplitudesSharpened");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT, LABEL_DOUBLE, "rlnLogAmplitudesIntercept");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_RESOL_SQUARED, LABEL_DOUBLE, "rlnResolutionSquared");
        MDL::addLabel(RLN_POSTPROCESS_MTF_VALUE, LABEL_DOUBLE, "rlnMtfValue");

        MDL::addLabel(RLN_SAMPLING_IS_3D, LABEL_BOOL, "rlnIs3DSampling");
        MDL::addLabel(RLN_SAMPLING_IS_3D_TRANS, LABEL_BOOL, "rlnIs3DTranslationalSampling");
        MDL::addLabel(RLN_SAMPLING_HEALPIX_ORDER, LABEL_INT, "rlnHealpixOrder");
        MDL::addLabel(RLN_SAMPLING_LIMIT_TILT, LABEL_DOUBLE, "rlnTiltAngleLimit");
        MDL::addLabel(RLN_SAMPLING_OFFSET_RANGE, LABEL_DOUBLE, "rlnOffsetRange");
        MDL::addLabel(RLN_SAMPLING_OFFSET_STEP, LABEL_DOUBLE, "rlnOffsetStep");
        MDL::addLabel(RLN_SAMPLING_PERTURB, LABEL_DOUBLE, "rlnSamplingPerturbInstance");
        MDL::addLabel(RLN_SAMPLING_PERTURBATION_FACTOR, LABEL_DOUBLE, "rlnSamplingPerturbFactor");
        MDL::addLabel(RLN_SAMPLING_PSI_STEP, LABEL_DOUBLE, "rlnPsiStep");
        MDL::addLabel(RLN_SAMPLING_SYMMETRY, LABEL_STRING, "rlnSymmetryGroup");

        MDL::addLabel(RLN_SELECTED, LABEL_BOOL, "rlnSelected");
        MDL::addLabel(RLN_SELECT_PARTICLES_ZSCORE, LABEL_DOUBLE, "rlnParticleSelectZScore");
        MDL::addLabel(RLN_SORTED_IDX, LABEL_SIZET, "rlnSortedIndex");
        MDL::addLabel(RLN_PERFRAME_CUMULATIVE_WEIGHT, LABEL_DOUBLE, "rlnPerFrameCumulativeWeight");
        MDL::addLabel(RLN_PERFRAME_RELATIVE_WEIGHT, LABEL_DOUBLE, "rlnPerFrameRelativeWeight");

        MDL::addLabel(RLN_RESOLUTION, LABEL_DOUBLE, "rlnResolution");
        MDL::addLabel(RLN_RESOLUTION_ANGSTROM, LABEL_DOUBLE, "rlnAngstromResolution");
        MDL::addLabel(RLN_RESOLUTION_INVPIXEL, LABEL_DOUBLE, "rlnResolutionInversePixel");
        MDL::addLabel(RLN_SPECTRAL_IDX, LABEL_INT, "rlnSpectralIndex");

        // "comment" already exists in xmipp
        // MDL::addLabel(BSOFT_COMMENT, LABEL_STRING, "comment");
        MDL::addLabel(BSOFT_ID, LABEL_STRING, "id");
        MDL::addLabel(BSOFT_PROJECT, LABEL_STRING, "project");
        MDL::addLabel(BSOFT_FIELD, LABEL_STRING, "field");
        MDL::addLabel(BSOFT_FIELD_ID, LABEL_STRING, "field.id");
        MDL::addLabel(BSOFT_MAP, LABEL_STRING, "map.3D_reconstruction");
        MDL::addLabel(BSOFT_MAP_ID, LABEL_STRING, "map.3D_reconstruction.id");
        MDL::addLabel(BSOFT_MAP_REFERENCE, LABEL_STRING, "map.reference.file_name");
        MDL::addLabel(BSOFT_MAP_RECONSTRUCTION, LABEL_STRING, "map.3D_reconstruction.file_name");
        MDL::addLabel(BSOFT_MAP_TRANSFORM_FILE, LABEL_STRING, "map.3D_reconstruction_fourier_transform.file_name");
        MDL::addLabel(BSOFT_MAP_POWERSPEC_FILE, LABEL_STRING, "map.3D_reconstruction_powerspectrum.file_name");
        MDL::addLabel(BSOFT_MAP_SIZE_X, LABEL_STRING, "map.3D_reconstruction.size_x");
        MDL::addLabel(BSOFT_MAP_SIZE_Y, LABEL_STRING, "map.3D_reconstruction.size_y");
        MDL::addLabel(BSOFT_MAP_SIZE_Z, LABEL_STRING, "map.3D_reconstruction.size_z");
        MDL::addLabel(BSOFT_MAP_ORIGIN_X, LABEL_STRING, "map.3D_reconstruction.origin_x");
        MDL::addLabel(BSOFT_MAP_ORIGIN_Y, LABEL_STRING, "map.3D_reconstruction.origin_y");
        MDL::addLabel(BSOFT_MAP_ORIGIN_Z, LABEL_STRING, "map.3D_reconstruction.origin_z");
        MDL::addLabel(BSOFT_MAP_SCALE_X, LABEL_STRING, "map.3D_reconstruction.scale_x");
        MDL::addLabel(BSOFT_MAP_SCALE_Y, LABEL_STRING, "map.3D_reconstruction.scale_y");
        MDL::addLabel(BSOFT_MAP_SCALE_Z, LABEL_STRING, "map.3D_reconstruction.scale_z");
        MDL::addLabel(BSOFT_MAP_VOXEL_SIZE, LABEL_STRING, "map.3D_reconstruction.voxel_size");
        MDL::addLabel(BSOFT_MAP_SELECT, LABEL_STRING, "map.3D_reconstruction.select");
        MDL::addLabel(BSOFT_MAP_FOM, LABEL_STRING, "map.3D_reconstruction.fom");
        MDL::addLabel(BSOFT_MAP_MAGNIFICATION, LABEL_STRING, "map.magnification");
        MDL::addLabel(BSOFT_MAP_VIEW_X, LABEL_STRING, "map.view_x");
        MDL::addLabel(BSOFT_MAP_VIEW_Y, LABEL_STRING, "map.view_y");
        MDL::addLabel(BSOFT_MAP_VIEW_Z, LABEL_STRING, "map.view_z");
        MDL::addLabel(BSOFT_MAP_VIEW_ANGLE, LABEL_STRING, "map.view_angle");
        MDL::addLabel(BSOFT_MAP_BACK_RWEIGHT, LABEL_STRING, "map.back_projection.rweight");
        MDL::addLabel(BSOFT_MAP_MODEL, LABEL_STRING, "map.model");
        MDL::addLabel(BSOFT_MAP_SYMMETRY, LABEL_STRING, "map.symmetry");
        // "micrograph" label is already in xmipp.
        //MDL::addLabel(BSOFT_MICROGRAPH, LABEL_STRING, "micrograph");
        MDL::addLabel(BSOFT_MICROGRAPH_FILE, LABEL_STRING, "micrograph.file_name");
        MDL::addLabel(BSOFT_MICROGRAPH_PARTICLE_FILE, LABEL_STRING, "micrograph_particle.file_name");
        MDL::addLabel(BSOFT_MICROGRAPH_FILAMENT_FILE, LABEL_STRING, "micrograph_filament.file_name");
        MDL::addLabel(BSOFT_MICROGRAPH_TRANSFORM_FILE, LABEL_STRING, "micrograph_fourier_transform.file_name");
        MDL::addLabel(BSOFT_MICROGRAPH_POWERSPEC_FILE, LABEL_STRING, "micrograph_powerspectrum.file_name");
        MDL::addLabel(BSOFT_MICROGRAPH_ID, LABEL_STRING, "micrograph.id");
        MDL::addLabel(BSOFT_MICROGRAPH_FIELD_ID, LABEL_STRING, "micrograph.field_id");
        MDL::addLabel(BSOFT_MICROGRAPH_NUMBER, LABEL_STRING, "micrograph.number");
        MDL::addLabel(BSOFT_MICROGRAPH_SELECT, LABEL_STRING, "micrograph.select");
        MDL::addLabel(BSOFT_MICROGRAPH_FOM, LABEL_STRING, "micrograph.fom");
        MDL::addLabel(BSOFT_MICROGRAPH_MAGNIFICATION, LABEL_STRING, "micrograph.magnification");
        MDL::addLabel(BSOFT_MICROGRAPH_SAMPLING, LABEL_STRING, "micrograph.sampling");
        MDL::addLabel(BSOFT_MICROGRAPH_PIXEL, LABEL_STRING, "micrograph.pixel_size");
        MDL::addLabel(BSOFT_MICROGRAPH_UNITS, LABEL_STRING, "micrograph.units");
        MDL::addLabel(BSOFT_MICROGRAPH_DOSE, LABEL_STRING, "micrograph.electron_dose");
        MDL::addLabel(BSOFT_MICROGRAPH_ORIGIN_X, LABEL_STRING, "micrograph.origin_x");
        MDL::addLabel(BSOFT_MICROGRAPH_ORIGIN_Y, LABEL_STRING, "micrograph.origin_y");
        MDL::addLabel(BSOFT_MICROGRAPH_ORIGIN_Z, LABEL_STRING, "micrograph.origin_z");
        MDL::addLabel(BSOFT_MICROGRAPH_SCALE_X, LABEL_STRING, "micrograph.scale_x");
        MDL::addLabel(BSOFT_MICROGRAPH_SCALE_Y, LABEL_STRING, "micrograph.scale_y");
        MDL::addLabel(BSOFT_MICROGRAPH_SCALE_Z, LABEL_STRING, "micrograph.scale_z");
        MDL::addLabel(BSOFT_MICROGRAPH_TILT_AXIS, LABEL_STRING, "micrograph.tilt_axis");
        MDL::addLabel(BSOFT_MICROGRAPH_TILT_ANGLE, LABEL_STRING, "micrograph.tilt_angle");
        MDL::addLabel(BSOFT_MICROGRAPH_LEVEL_ANGLE, LABEL_STRING, "micrograph.level_angle");
        MDL::addLabel(BSOFT_MICROGRAPH_ROT_ANGLE, LABEL_STRING, "micrograph.rotation_angle");
        MDL::addLabel(BSOFT_MICROGRAPH_VIEW_X, LABEL_STRING, "micrograph.view_x");
        MDL::addLabel(BSOFT_MICROGRAPH_VIEW_Y, LABEL_STRING, "micrograph.view_y");
        MDL::addLabel(BSOFT_MICROGRAPH_VIEW_Z, LABEL_STRING, "micrograph.view_z");
        MDL::addLabel(BSOFT_MICROGRAPH_VIEW_ANGLE, LABEL_STRING, "micrograph.view_angle");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_1_1, LABEL_STRING, "micrograph.matrix_1_1");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_1_2, LABEL_STRING, "micrograph.matrix_1_2");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_1_3, LABEL_STRING, "micrograph.matrix_1_3");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_2_1, LABEL_STRING, "micrograph.matrix_2_1");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_2_2, LABEL_STRING, "micrograph.matrix_2_2");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_2_3, LABEL_STRING, "micrograph.matrix_2_3");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_3_1, LABEL_STRING, "micrograph.matrix_3_1");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_3_2, LABEL_STRING, "micrograph.matrix_3_2");
        MDL::addLabel(BSOFT_MICROGRAPH_MATRIX_3_3, LABEL_STRING, "micrograph.matrix_3_3");
        MDL::addLabel(BSOFT_MICROGRAPH_HVEC_X, LABEL_STRING, "micrograph.h_x");
        MDL::addLabel(BSOFT_MICROGRAPH_HVEC_Y, LABEL_STRING, "micrograph.h_y");
        MDL::addLabel(BSOFT_MICROGRAPH_HVEC_Z, LABEL_STRING, "micrograph.h_z");
        MDL::addLabel(BSOFT_MICROGRAPH_KVEC_X, LABEL_STRING, "micrograph.k_x");
        MDL::addLabel(BSOFT_MICROGRAPH_KVEC_Y, LABEL_STRING, "micrograph.k_y");
        MDL::addLabel(BSOFT_MICROGRAPH_KVEC_Z, LABEL_STRING, "micrograph.k_z");
        MDL::addLabel(BSOFT_MICROGRAPH_LVEC_X, LABEL_STRING, "micrograph.l_x");
        MDL::addLabel(BSOFT_MICROGRAPH_LVEC_Y, LABEL_STRING, "micrograph.l_y");
        MDL::addLabel(BSOFT_MICROGRAPH_LVEC_Z, LABEL_STRING, "micrograph.l_z");
        MDL::addLabel(BSOFT_MICROGRAPH_HELIX_AXIS, LABEL_STRING, "micrograph.helix_axis_angle");
        MDL::addLabel(BSOFT_MICROGRAPH_HELIX_RISE, LABEL_STRING, "micrograph.helix_subunit_rise");
        MDL::addLabel(BSOFT_MICROGRAPH_HELIX_ANGLE, LABEL_STRING, "micrograph.helix_subunit_angle");
        MDL::addLabel(BSOFT_MICROGRAPH_HELIX_RADIUS, LABEL_STRING, "micrograph.helix_radius");
        MDL::addLabel(BSOFT_MICROGRAPH_VOLTAGE, LABEL_STRING, "micrograph.voltage");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_CS, LABEL_STRING, "micrograph.ctf.Cs");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_CC, LABEL_STRING, "micrograph.ctf.Cc");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_ALPHA, LABEL_STRING, "micrograph.ctf.alpha");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_DE, LABEL_STRING, "micrograph.ctf.energy_spread");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_AMP_CONT, LABEL_STRING, "micrograph.ctf.amp_contrast");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_ZERO, LABEL_STRING, "micrograph.ctf.first_zero");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_DEF_AVG, LABEL_STRING, "micrograph.ctf.defocus_average");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_DEF_DEV, LABEL_STRING, "micrograph.ctf.defocus_deviation");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_DEF_MIN, LABEL_STRING, "micrograph.ctf.defocus_min");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_DEF_MAX, LABEL_STRING, "micrograph.ctf.defocus_max");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_AST_ANG, LABEL_STRING, "micrograph.ctf.astigmatism_angle");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_BASELINE, LABEL_STRING, "micrograph.ctf.baseline");
        MDL::addLabel(BSOFT_MICROGRAPH_CTF_ENVELOPE, LABEL_STRING, "micrograph.ctf.envelope");
        MDL::addLabel(BSOFT_MICROGRAPH_BOX_RADIUS, LABEL_STRING, "micrograph.box_radius");
        MDL::addLabel(BSOFT_MICROGRAPH_BOX_RADIUS_X, LABEL_STRING, "micrograph.box_radius_x");
        MDL::addLabel(BSOFT_MICROGRAPH_BOX_RADIUS_Y, LABEL_STRING, "micrograph.box_radius_y");
        MDL::addLabel(BSOFT_MICROGRAPH_BOX_RADIUS_Z, LABEL_STRING, "micrograph.box_radius_z");
        MDL::addLabel(BSOFT_MICROGRAPH_BAD, LABEL_STRING, "micrograph.bad");
        MDL::addLabel(BSOFT_MICROGRAPH_BAD_RADIUS, LABEL_STRING, "micrograph.bad_radius");
        MDL::addLabel(BSOFT_MICROGRAPH_BAD_X, LABEL_STRING, "micrograph.bad_x");
        MDL::addLabel(BSOFT_MICROGRAPH_BAD_Y, LABEL_STRING, "micrograph.bad_y");
        MDL::addLabel(BSOFT_MICROGRAPH_BAD_Z, LABEL_STRING, "micrograph.bad_z");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_RADIUS, LABEL_STRING, "micrograph.marker_radius");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_ID, LABEL_STRING, "micrograph.marker_id");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_X, LABEL_STRING, "micrograph.marker_x");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_Y, LABEL_STRING, "micrograph.marker_y");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_Z, LABEL_STRING, "micrograph.marker_z");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_ERROR_X, LABEL_STRING, "micrograph.marker_error_x");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_ERROR_Y, LABEL_STRING, "micrograph.marker_error_y");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_ERROR_Z, LABEL_STRING, "micrograph.marker_error_z");
        MDL::addLabel(BSOFT_MICROGRAPH_MARKER_FOM, LABEL_STRING, "micrograph.marker_fom");
        MDL::addLabel(BSOFT_MICROGRAPH_FILAMENT_WIDTH, LABEL_STRING, "micrograph.filament_width");
        MDL::addLabel(BSOFT_MICROGRAPH_FILNODE_RADIUS, LABEL_STRING, "micrograph.filament_node_radius");
        MDL::addLabel(BSOFT_CTF, LABEL_STRING, "ctf");
        MDL::addLabel(BSOFT_CTF_VOLTAGE, LABEL_STRING, "ctf.voltage");
        MDL::addLabel(BSOFT_CTF_CS, LABEL_STRING, "ctf.Cs");
        MDL::addLabel(BSOFT_CTF_CC, LABEL_STRING, "ctf.Cc");
        MDL::addLabel(BSOFT_CTF_ALPHA, LABEL_STRING, "ctf.alpha");
        MDL::addLabel(BSOFT_CTF_DE, LABEL_STRING, "ctf.energy_spread");
        MDL::addLabel(BSOFT_CTF_AMP, LABEL_STRING, "ctf.amp_contrast");
        MDL::addLabel(BSOFT_CTF_ZERO, LABEL_STRING, "ctf.first_zero");
        MDL::addLabel(BSOFT_CTF_DEF_AVG, LABEL_STRING, "ctf.defocus_average");
        MDL::addLabel(BSOFT_CTF_DEF_DEV, LABEL_STRING, "ctf.defocus_deviation");
        MDL::addLabel(BSOFT_CTF_DEF_MIN, LABEL_STRING, "ctf.defocus_min");
        MDL::addLabel(BSOFT_CTF_DEF_MAX, LABEL_STRING, "ctf.defocus_max");
        MDL::addLabel(BSOFT_CTF_AST_ANG, LABEL_STRING, "ctf.astigmatism_angle");
        MDL::addLabel(BSOFT_CTF_BASELINE, LABEL_STRING, "ctf.baseline");
        MDL::addLabel(BSOFT_CTF_ENVELOPE, LABEL_STRING, "ctf.envelope");
        MDL::addLabel(BSOFT_PARTICLE, LABEL_STRING, "particle");
        MDL::addLabel(BSOFT_PARTICLE_FILE, LABEL_STRING, "particle.file_name");
        MDL::addLabel(BSOFT_PARTICLE_NUMBER, LABEL_STRING, "particle.number");
        MDL::addLabel(BSOFT_PARTICLE_ID, LABEL_STRING, "particle.id");
        MDL::addLabel(BSOFT_PARTICLE_GROUP, LABEL_STRING, "particle.group_id");
        MDL::addLabel(BSOFT_PARTICLE_MG_ID, LABEL_STRING, "particle.micrograph_id");
        MDL::addLabel(BSOFT_PARTICLE_MG_X, LABEL_STRING, "particle.micrograph_x");
        MDL::addLabel(BSOFT_PARTICLE_MG_Y, LABEL_STRING, "particle.micrograph_y");
        MDL::addLabel(BSOFT_PARTICLE_MG_Z, LABEL_STRING, "particle.micrograph_z");
        MDL::addLabel(BSOFT_PARTICLE_X, LABEL_STRING, "particle.x");
        MDL::addLabel(BSOFT_PARTICLE_Y, LABEL_STRING, "particle.y");
        MDL::addLabel(BSOFT_PARTICLE_Z, LABEL_STRING, "particle.z");
        MDL::addLabel(BSOFT_PARTICLE_X_ORIGIN, LABEL_STRING, "particle.x_origin");
        MDL::addLabel(BSOFT_PARTICLE_Y_ORIGIN, LABEL_STRING, "particle.y_origin");
        MDL::addLabel(BSOFT_PARTICLE_Z_ORIGIN, LABEL_STRING, "particle.z_origin");
        MDL::addLabel(BSOFT_PARTICLE_ORIGIN_X, LABEL_STRING, "particle.origin_x");
        MDL::addLabel(BSOFT_PARTICLE_ORIGIN_Y, LABEL_STRING, "particle.origin_y");
        MDL::addLabel(BSOFT_PARTICLE_ORIGIN_Z, LABEL_STRING, "particle.origin_z");
        MDL::addLabel(BSOFT_PARTICLE_PSI, LABEL_STRING, "particle.psi");
        MDL::addLabel(BSOFT_PARTICLE_THETA, LABEL_STRING, "particle.theta");
        MDL::addLabel(BSOFT_PARTICLE_PHI, LABEL_STRING, "particle.phi");
        MDL::addLabel(BSOFT_PARTICLE_OMEGA, LABEL_STRING, "particle.omega");
        MDL::addLabel(BSOFT_PARTICLE_VIEW_X, LABEL_STRING, "particle.view_x");
        MDL::addLabel(BSOFT_PARTICLE_VIEW_Y, LABEL_STRING, "particle.view_y");
        MDL::addLabel(BSOFT_PARTICLE_VIEW_Z, LABEL_STRING, "particle.view_z");
        MDL::addLabel(BSOFT_PARTICLE_VIEW_ANGLE, LABEL_STRING, "particle.view_angle");
        MDL::addLabel(BSOFT_PARTICLE_MAGNIF, LABEL_STRING, "particle.magnification");
        MDL::addLabel(BSOFT_PARTICLE_DEFOCUS, LABEL_STRING, "particle.defocus");
        MDL::addLabel(BSOFT_PARTICLE_DEF_DEV, LABEL_STRING, "particle.defocus_deviation");
        MDL::addLabel(BSOFT_PARTICLE_AST_ANG, LABEL_STRING, "particle.astigmatism_angle");
        MDL::addLabel(BSOFT_PARTICLE_SELECT, LABEL_STRING, "particle.select");
        MDL::addLabel(BSOFT_PARTICLE_FOM, LABEL_STRING, "particle.fom");
        MDL::addLabel(BSOFT_PARTICLE_FOM_CV, LABEL_STRING, "particle.fom_crossvalidation");
        MDL::addLabel(BSOFT_PARTICLE_FOM_AVG, LABEL_STRING, "particle.fom_average");
        MDL::addLabel(BSOFT_PARTICLE_FOM_STD, LABEL_STRING, "particle.fom_stdev");
        MDL::addLabel(BSOFT_PARTICLE_HANDA_FOM, LABEL_STRING, "particle.handa_fom");
        MDL::addLabel(BSOFT_PARTICLE_HANDB_FOM, LABEL_STRING, "particle.handb_fom");
        MDL::addLabel(BSOFT_PARTICLE_CC, LABEL_STRING, "particle.cc");
        MDL::addLabel(BSOFT_PARTICLE_PFT_CC, LABEL_STRING, "particle.pft_cc");
        MDL::addLabel(BSOFT_PARTICLE_PRJ_CC, LABEL_STRING, "particle.prj_cc");
        MDL::addLabel(BSOFT_PARTICLE_CMP_CC, LABEL_STRING, "particle.cmp_cc");
        MDL::addLabel(BSOFT_PARTICLE_RFACTORAB, LABEL_STRING, "particle.rfactorab");
        MDL::addLabel(BSOFT_PARTICLE_COVERAGE, LABEL_STRING, "particle.coverage");
        MDL::addLabel(BSOFT_PARTICLE_BOX_SIZE, LABEL_DOUBLE, "particle.box_size");
        MDL::addLabel(BSOFT_PARTICLE_BOX_SIZE_X, LABEL_DOUBLE, "particle.box_size_x");
        MDL::addLabel(BSOFT_PARTICLE_BOX_SIZE_Y, LABEL_DOUBLE, "particle.box_size_y");
        MDL::addLabel(BSOFT_PARTICLE_BOX_SIZE_Z, LABEL_DOUBLE, "particle.box_size_z");
        MDL::addLabel(BSOFT_PARTICLE_BOX_RADIUS, LABEL_STRING, "particle.box_radius");
        MDL::addLabel(BSOFT_PARTICLE_BOX_RADIUS_X, LABEL_STRING, "particle.box_radius_x");
        MDL::addLabel(BSOFT_PARTICLE_BOX_RADIUS_Y, LABEL_STRING, "particle.box_radius_y");
        MDL::addLabel(BSOFT_PARTICLE_BOX_RADIUS_Z, LABEL_STRING, "particle.box_radius_z");
        MDL::addLabel(BSOFT_PARTICLE_BAD, LABEL_STRING, "particle.bad");
        MDL::addLabel(BSOFT_PARTICLE_BAD_RADIUS, LABEL_STRING, "particle.bad_radius");
        MDL::addLabel(BSOFT_PARTICLE_BAD_X, LABEL_STRING, "particle.bad_x");
        MDL::addLabel(BSOFT_PARTICLE_BAD_Y, LABEL_STRING, "particle.bad_y");
        MDL::addLabel(BSOFT_PARTICLE_BAD_Z, LABEL_STRING, "particle.bad_z");
        MDL::addLabel(BSOFT_FILAMENT, LABEL_STRING, "filament");
        MDL::addLabel(BSOFT_FILAMENT_FILE, LABEL_STRING, "filament.file_name");
        MDL::addLabel(BSOFT_FILAMENT_ID, LABEL_STRING, "filament.id");
        MDL::addLabel(BSOFT_FILAMENT_NODE, LABEL_STRING, "filament.node");
        MDL::addLabel(BSOFT_FILAMENT_NODE_ID, LABEL_STRING, "filament.node_id");
        MDL::addLabel(BSOFT_FILAMENT_NODE_X, LABEL_STRING, "filament.x");
        MDL::addLabel(BSOFT_FILAMENT_NODE_Y, LABEL_STRING, "filament.y");
        MDL::addLabel(BSOFT_FILAMENT_NODE_Z, LABEL_STRING, "filament.z");
        MDL::addLabel(BSOFT_FILAMENT_WIDTH, LABEL_STRING, "filament.width");
        MDL::addLabel(BSOFT_FILNODE_RADIUS, LABEL_STRING, "filament.node_radius");
        MDL::addLabel(BSOFT_ORIENT_ID, LABEL_STRING, "orient.id");
        MDL::addLabel(BSOFT_ORIENT_ORIGIN_X, LABEL_STRING, "orient.origin_x");
        MDL::addLabel(BSOFT_ORIENT_ORIGIN_Y, LABEL_STRING, "orient.origin_y");
        MDL::addLabel(BSOFT_ORIENT_ORIGIN_Z, LABEL_STRING, "orient.origin_z");
        MDL::addLabel(BSOFT_ORIENT_VIEW_X, LABEL_STRING, "orient.view_x");
        MDL::addLabel(BSOFT_ORIENT_VIEW_Y, LABEL_STRING, "orient.view_y");
        MDL::addLabel(BSOFT_ORIENT_VIEW_Z, LABEL_STRING, "orient.view_z");
        MDL::addLabel(BSOFT_ORIENT_VIEW_ANGLE, LABEL_STRING, "orient.view_angle");
        MDL::addLabel(BSOFT_ORIENT_FOM, LABEL_STRING, "orient.fom");
        MDL::addLabel(BSOFT_ORIENT_SELECT, LABEL_STRING, "orient.select");
        MDL::addLabel(BSOFT_MARKER, LABEL_STRING, "marker");
        MDL::addLabel(BSOFT_MARKER_RADIUS, LABEL_STRING, "marker.radius");
        MDL::addLabel(BSOFT_MARKER_ID, LABEL_STRING, "marker.id");
        MDL::addLabel(BSOFT_MARKER_X, LABEL_STRING, "marker.x");
        MDL::addLabel(BSOFT_MARKER_Y, LABEL_STRING, "marker.y");
        MDL::addLabel(BSOFT_MARKER_Z, LABEL_STRING, "marker.z");
        MDL::addLabel(BSOFT_MARKER_ERROR_X, LABEL_STRING, "marker.error_x");
        MDL::addLabel(BSOFT_MARKER_ERROR_Y, LABEL_STRING, "marker.error_y");
        MDL::addLabel(BSOFT_MARKER_ERROR_Z, LABEL_STRING, "marker.error_z");
        MDL::addLabel(BSOFT_MARKER_IMAGE, LABEL_STRING, "marker.image_number");
        MDL::addLabel(BSOFT_MARKER_RESIDUAL, LABEL_STRING, "marker.residual");
        MDL::addLabel(BSOFT_MARKER_FOM, LABEL_STRING, "marker.fom");
        MDL::addLabel(BSOFT_MARKER_SELECT, LABEL_STRING, "marker.select");
        MDL::addLabel(BSOFT_REFLEX, LABEL_STRING, "refln");
        MDL::addLabel(BSOFT_REFLEX_RADIUS, LABEL_STRING, "refln.radius");
        MDL::addLabel(BSOFT_REFLEX_X, LABEL_STRING, "refln.location_x");
        MDL::addLabel(BSOFT_REFLEX_Y, LABEL_STRING, "refln.location_y");
        MDL::addLabel(BSOFT_REFLEX_Z, LABEL_STRING, "refln.location_z");
        MDL::addLabel(BSOFT_REFLEX_H, LABEL_STRING, "refln.index_h");
        MDL::addLabel(BSOFT_REFLEX_K, LABEL_STRING, "refln.index_k");
        MDL::addLabel(BSOFT_REFLEX_L, LABEL_STRING, "refln.index_l");
        MDL::addLabel(BSOFT_REFLEX_AMP, LABEL_STRING, "refln.F_meas_au");
        MDL::addLabel(BSOFT_REFLEX_SIGAMP, LABEL_STRING, "refln.F_meas_sigma_au");
        MDL::addLabel(BSOFT_REFLEX_PHI, LABEL_STRING, "refln.phase_meas");
        MDL::addLabel(BSOFT_REFLEX_SIGPHI, LABEL_STRING, "refln.phase_sigma");
        MDL::addLabel(BSOFT_REFLEX_FOM, LABEL_STRING, "refln.weight");
        MDL::addLabel(BSOFT_REFLEX_STATUS, LABEL_STRING, "refln.status");
        MDL::addLabel(BSOFT_LAYERLINE, LABEL_STRING, "layer_line");
        MDL::addLabel(BSOFT_LAYERLINE_NUMBER, LABEL_STRING, "layer_line.number");
        MDL::addLabel(BSOFT_LAYERLINE_ORDER, LABEL_STRING, "layer_line.bessel_order");
        MDL::addLabel(BSOFT_LAYERLINE_DISTANCE, LABEL_STRING, "layer_line.distance");
        MDL::addLabel(BSOFT_LAYERLINE_FREQ, LABEL_STRING, "layer_line.frequency");
        MDL::addLabel(BSOFT_LAYERLINE_AMP, LABEL_STRING, "layer_line.amplitude");
        MDL::addLabel(BSOFT_LAYERLINE_FOM, LABEL_STRING, "layer_line.fom");
        MDL::addLabel(BSOFT_LAYERLINE_SELECT, LABEL_STRING, "layer_line.select");

        MDL::addLabel(BSOFT_SYMMETRY_INT_TABLES_NUMBER, LABEL_INT, "symmetry.Int_Tables_number");
        MDL::addLabel(BSOFT_SYMMETRY_SPACE_GROUP_NAME_H_M, LABEL_STRING, "symmetry.space_group_name_H-M");
        MDL::addLabel(BSOFT_SYMMETRY_CELL_SETTING, LABEL_STRING, "symmetry.cell_setting");
        MDL::addLabel(BSOFT_SYMMETRY_EQUIV_ID, LABEL_INT, "symmetry_equiv.id");
        MDL::addLabel(BSOFT_SYMMETRY_EQUIV_POS_AS_XYZ, LABEL_STRING, "symmetry_equiv.pos_as_xyz");

        //Create an static empty header for image initialization
        MDL::emptyHeader.resetGeo();
        MDL::emptyHeader.setValue(MDL_ANGLE_ROT, 0.);
        MDL::emptyHeader.setValue(MDL_ANGLE_TILT,0.);

        // Add user-defined label aliases for allow use Xmipp-MetaData
        // with other labels
        MDL::addExtraAliases();

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
