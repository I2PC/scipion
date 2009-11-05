/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
#ifndef _PROG_EVALUATE_HH
#  define _PROG_EVALUATE_HH

#include <data/funcs.h>
#include <data/volume.h>

#include <data/phantom.h>

/**@defgroup EvaluateFOMS foms_evaluate (FOMs evaluation program)
   @ingroup ReconsLibraryPrograms */
//@{
/* Evaluation Program Parameters ------------------------------------------- */
/** Parameter class for the evaluate program */
class Prog_Evaluate_Parameters
{
public:
    /// SelFile with all reconstructions
    FileName fn_sel;

    /// Phantom filename: either Xmipp volume or phantom description
    FileName fnPhantom;
    /// Reconstruction filename (Xmipp volume)
    FileName fn_recons;
    /** Percent of mass to be left out in slice histograms.
        The slice histograms are usually dominated by the background value,
        with this percentage mass you can leave out most of the background
        and what you are presented is the histogram of all values above
        a calculated threshold (adapted to the percentage mass) slice by
        slice. This value is usually around 99% */
    double    percent_mass;
    /// Angles for the main oversampled direction
    double    RSrot, RStilt;
    /** Global sphere mask radius.
        If you give this value (i.e., different from 0) then when computing
        structural consistency and directional FOMs those voxels outside
        the sphere with center at the volume center and this radius are
        not taken into the computations */
    double    global_radius;
    /** Surface mask applied.
        It is used for not taking into account those voxels outside
        the mask for a fair comparison. */
    FileName fn_mask;
    /** Apply gray fitting between volumes before evaluating. */
    bool fit_gray_scales;

    /** Background definition.
       The background of a feature is defined either by a sphere of radius
       back_radius or a scaled version of that same feature. The way of
       knowing which one it is, is using the flag back_mode which can only
       take two possible values ENLARGE_MODE or SPHERE_MODE. */
    //@{
    /// Enlarged or sphere background
    int      back_mode;
    /// Background sphere radius
    double    back_radius;
    /// Enlarging factor
    double    back_factor;
    //@}

#define  SAVE_MAPS 0x1
#define  SHOW_VALUES 0x2
#define  SHOW_PROCESS 0X4
#define  SAVE_HISTOGRAMS 0x8
#define  ONLY_STRUCTURAL 0x10
    /** Debugging level.
        The debugging level is a bitwise field where you can set
        any of the following bits:
        \\ SAVE_MAPS: save generated phantom, label map, distance map,
           difference, absolute difference and quadratic difference maps
           (volumes) with names:
           @code
           w0001_eval_phantom.vol
           w0001_eval_label.vol
           w0001_evaluateQuadratic_map.vol
           w0001_eval_difference_map.vol
           w0001_eval_absolute_map.vol
           w0001_eval_distance_map.vol
           @endcode
        \\ SHOW_VALUES: you can ask interactively to show the voxel
           values for this feature.
        \\ SHOW_PROCESS: routines show information about the main
           variables used in the FOM computation for each feature
        \\ SAVE_HISTOGRAMS: in the histogram based FOMs this flag
           forces to dump the histogram into a file (the filename is
           given in the function call) or the screen.
        \\ ONLY_STRUCTURAL: This is used to train the algorithms.
           Only the structural consistency FOMs are computed.
    */
    int      tell;
public:
    /** Set default values for program parameters.
        Backmode=Enlarge_mode by 1.5, percentage_mass=99, directional
        FOMs from Z axis, and global_radius mask=0. */
    void default_values();

    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introdustd::cing this parameters. */
    void usage();

    /** std::cout << prm; */
    friend std::ostream & operator << (std::ostream &out,
                                  const Prog_Evaluate_Parameters &prm);
};

/** Evaluate program Side information.
    This class contains side necessary information for the Evaluate program.
    This information can be obtained from the parameters and is basically
    the Xmipp phantom and reconstruction, the labelling and mask (if
    global_radius!=0). */
class EVALUATE_Side_Info
{
public:
    /// Root of the reconstruction filename
    FileName fn_root;
    /// Phantom description file (if necessary)
    Phantom           phantom_descr;
    /// Phantom Xmipp volume (internally generated or externally given)
    VolumeXmipp       vol_phantom;
    /// Reconstructed Xmipp volume
    VolumeXmipp       vol_recons;
    /** Labelled volume.
        For instance, vol_label[i][j][k]=3 means that the pixel in
        position (i,j,k) belongs to the feature, object, number 3.
        ...=0 means that it belongs to the background */
    VolumeXmipp       vol_label;
    /** Global mask.
        If you give a global_radius then all voxels outside this
        mask are not taken into the computations */
    VolumeXmipp       vol_mask;
    /** Distance volume.
        For a voxel in the background the distance is defined as the
        minimum distance to a feature, and for a feature voxel the
        distance is the minimum distance to the background. */
    VolumeXmipp       vol_distance;
    /// Number of features present
    int               num_feat;
    /** Paremeter for generating background.
        It is either a radius or a enlarging factor */
    double             back_param;
#define XMIPP_PHANTOM 1
#define MATH_PHANTOM  2
    /** Phantom description.
        Flag to indicate if the phantom has been given with a mathematical
        description or a Xmipp volume. It can only take two values:
        XMIPP_PHANTOM or MATH_PHANTOM */
    int               descr_mode;
public:
    /** Produce Evaluate Side information.
        This function produce the side information from the evaluate
        program parameters. Basically it loads the phantom (generates a
        Xmipp volume with the phantom if it was given as a mathematical
        description), generates a global sphere mask (if necessary)
        and label the phantom volume according to the different
        features. Labelling Xmipp volumes is not implemented yet. */
    void produce_Side_Info(const Prog_Evaluate_Parameters &prm);
};

/** Evaluation result class.
    This class contains the results after the evaluation of the volume
    indicated in the program parameters. When the result is a matrix
    index 0 stands for the background, index 1 for feature 1, index 2
    for feature 2 ... */
class EVALUATE_results
{
public:
    //@{
    /// for features: L2 error
    Matrix1D<double> scL2_FOMs;
    /// for features: L1 error
    Matrix1D<double> scL1_FOMs;
    /// for features: error in the averages
    Matrix1D<double> scmu_FOMs;
    /// for features: error in the stddevs
    Matrix1D<double> scdev_FOMs;
    /// for features: error in the ranges
    Matrix1D<double> scrange_FOMs;
    /// for features: correlation between reconstruction and phantom
    Matrix1D<double> sccorr_FOMs;
    /// for features: mutual information
    Matrix1D<double> scinf_FOMs;
    /// global: L2 error
    double           scL2_FOM;
    /// global: L1 error
    double           scL1_FOM;
    /// global: weighted L2 error
    double           scL2w_FOM;
    /// global: weighted L1 error
    double           scL1w_FOM;
    /// global: error in the averages
    double           scmu_FOM;
    /// global: error in the stddevs
    double           scdev_FOM;
    /// global: error in the ranges
    double           scrange_FOM;
    /// global: correlation
    double           sccorr_FOM;
    /// global: mutual information
    double           scinf_FOM;
    /// global: resolution
    double           resol_FOM;
    //@}

    //@{
    /// for double cylinders: vertical error
    Matrix1D<double> hsvr_FOMs;
    /// for features: error in hot spot
    Matrix1D<double> hsmu_FOMs;
    /// for features: border separability
    Matrix1D<double> hsbr_FOMs;
    /// for features: error in detectability
    Matrix1D<double> hsdt_FOMs;
    /// global: error for double cyilinder elong
    double           hsvr_FOM;
    /// global: error for mean separability
    double           hsmu_FOM;
    /// global: border separability
    double           hsbr_FOM;
    /// global: error in detectability
    double           hsdt_FOM;
    //@}

    //@{
    /// global: Slice histograms
    ImageXmipp      img_histog;
    /// global: Radon Transform FOM
    double           drrt_FOM;
    //@}

    //@{
    /// global: blurring distance
    double           dsbl_FOM;
    /// global: appearance distanc
    double           dsad_FOM;
    //@}
};

/** Effectively compute the FOMs.
    This is the core function in the FOM calculations. The results
    are returned in a structure because they are too many. You may
    cause this function to output some information about the process
    using the \ref Prog_Evaluate_Parameters::tell variable.

    If you set the following flags then these other files are created:
    @code
    SHOW_PROCESS     -->
        w0001_eval_radon.plot:  with the Radon transforms
    SAVE_HISTOGRAMS  -->
        w0001_eval_histog.plot: with the feature histograms
    Without any flag -->
        w0001_eval_shape.plot:  with the shape information
        w0001_eval_slice_histog.img: with the slice histograms as an image
    @endcode */
void compute_FOMs(const Prog_Evaluate_Parameters &prm,
                  EVALUATE_Side_Info &side, EVALUATE_results &results);

/** Show parameters and results.
    This function shows all parameters and results, it also allows you
    to see interactively the voxel values in the phantom and in the
    reconstruction. */
void show_FOMs(const Prog_Evaluate_Parameters &prm,
               EVALUATE_Side_Info &side, const EVALUATE_results &results);

/** Main Evaluation routine.
    This is the main function of the program Evaluate. It takes a
    set of evaluation parameters and returns the results in a structure. */
void ROUT_Evaluate(Prog_Evaluate_Parameters &prm,
                   EVALUATE_results &results);

/** FOM class with all kin of FOMs. */
class FOMs
{
public:
    Matrix1D<double> scL2;
    Matrix1D<double> scL1;
    Matrix1D<double> scL2w;
    Matrix1D<double> scL1w;
    Matrix1D<double> scmu;
    Matrix1D<double> scdev;
    Matrix1D<double> scrange;
    Matrix1D<double> sccorr;
    Matrix1D<double> scinf;
    Matrix1D<double> scresol;
    Matrix1D<double> scL20;
    Matrix1D<double> scL10;
    Matrix1D<double> scmu0;
    Matrix1D<double> scdev0;
    Matrix1D<double> scrange0;
    Matrix1D<double> scL21;
    Matrix1D<double> scL11;
    Matrix1D<double> scmu1;
    Matrix1D<double> scdev1;
    Matrix1D<double> scrange1;
    Matrix1D<double> hsvr;
    Matrix1D<double> hsmu;
    Matrix1D<double> hsbr;
    Matrix1D<double> hsdt;
    Matrix1D<double> drrt;
    Matrix1D<double> dsbl;
    Matrix1D<double> dsad;
public:
    /// Empty constructor. n FOMs are expected
    FOMs(int n);

    /// Set FOMs form measure k from EVALUATE_Results
    void set_FOMs(int k, EVALUATE_results &results);

    /** Compute FOMs stats.
        The mean and stddev of the input FOMs are stored at position i in
        the output FOMs */
    friend void compute_FOMs_stats(const FOMs &foms, int i,
                                   FOMs &fmean, FOMs &fstddev);

    /** Show. */
    friend std::ostream & operator << (std::ostream &out, const FOMs &foms);

    /** Show mean and stddev. Show mean and stddev at position i. */
    friend void show_stats(std::ostream &out, int i, const FOMs &fmean,
                           const FOMs &fstddev);
};

//@}
#endif
