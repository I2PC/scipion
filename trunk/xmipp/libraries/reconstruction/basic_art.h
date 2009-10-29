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
#ifndef _BASIC_ART_HH
#  define _BASIC_ART_HH

#include <fstream>

#include <data/volume.h>
#include <data/selfile.h>

#include "refinement.h"
#include <data/grids.h>
#include "symmetries.h"
#include <data/basis.h>
#include "projection.h"
#include "fourier_filter.h"

struct Recons_info;

/**@defgroup BasicART Basic and common ART
   @ingroup ReconsLibraryPrograms 
    The main difference between ART applied to different cases (single
    particles, crystals, ...) is the single step applied to each case.
    Most of the tasks in the ART are common to all ART processes. All
    these common tasks as well as the common parameters are comprised
    in the Basic_art module. These common tasks are based on the existence
    of an extra_parameter structure containing all the specific information
    for the ART process.

    The user interface program should make a call to the \ref Basic_ROUT_Art
    routine with the corresponding extra_parameter structure.
*/
//@{
/* ART parameters ---------------------------------------------------------- */
/** ART basic parameters.
    This class contains all information needed about the ART process.
    See the user guide for more information about the ART parameters. */
class Basic_ART_Parameters
{
public:
    // Type of the parallel processing
    typedef enum {ART, pCAV, pAVSP, pSART, pBiCAV, pSIRT, pfSIRT, SIRT } t_parallel_mode;

    /* User parameters ...................................................... */
    //@{
    /// Basis function. By default, blobs
    Basis basis;

    /// Number of iterations
    int no_it;

    /// Relaxation parameter
    Matrix1D<double> lambda_list;

    /** Valid methods are ART, pCAV, pAVSP, pSART, pBiCAV, pSIRT and pfSIRT
        for parallel computation. This variable establish the way that particles are
        divided into blocks for parallel processing. If sequential
        processing is wanted, set it to ART or SIRT. This is the default.

        \\Ex: parallel_mode=Basic_ART_Parameters::ART*/
    t_parallel_mode parallel_mode;

    /// Number of projections for each parallel block
    int block_size;

    /** Valid mdoes are ARTK, CAVK and CAV.
        This is the mode of updating a single projection, it has nothing
        to do with the global ART or SIRT mode */
    int eq_mode;

    /// True if random sort of projections
    bool random_sort;

    /// True if no sort must be made
    bool dont_sort;

    /// Sort perpendicular with the last N projections. If -1 with all previous
    int sort_last_N;

    /** Flag for weighted least squares reconstruction (for maximum-likelihood) */
    bool WLS;

    /** Relaxation parameter for WLS residual volume */
    Matrix1D<double> kappa_list;

    /** Vector containing all residual images for wlsART */
    std::vector<Projection> residual_imgs;

    /** Sum ML-weight of all projections*/
    double sum_weight;

    /// Relative size for the grid
    double grid_relative_size;

    /// CC, BCC or FCC (in grids.hh)
    int grid_type;

    /** Projection extension.
        Sometimes for avoiding the box effects produced by the ART and SIRT
        algorithms, projections are extended in size such that the projection
        is supposed to be larger than it really is. The extended part is set
        to 0 telling so that no measure has been done in that area. This
        parameter is measured in pixels, and by default is 0.*/
    int proj_ext;

    /** Interest sphere. If -1 not considered. */
    double R;

    /** Output X size.
        If not given, the input images size is assumed. */
    int Xoutput_volume_size;

    /** Output Y size.
        If not given, the input images size is assumed. */
    int Youtput_volume_size;

    /** Output Z size.
        If not given, the input images size is assumed. */
    int Zoutput_volume_size;

    /// Sampling rate
    double sampling;

    /// File containing symmetries
    FileName fn_sym;

    /// Inpose simetry each sym_each iterations (iteration=projection)
    /// deafult is 0
    int sym_each;

    /// Skip projection with absolute tilt angle greater than max_tilt
    /// deafult is infinite
    double max_tilt;

    /// Refine the translation alignement after n projection presentations
    int ref_trans_after;

    /// Refine the translation alignement after n projection presentations
    double ref_trans_step;

    /// Sparse reconstruction
    double sparseEps;

    /// Tomographic diffussion
    double diffusionWeight;

    /// Force the reconstruction to be symmetric this number of times
    int force_sym;

    /// Do not generate symmetry subgroup
    bool do_not_generate_subgroup;

    /// Do not use symmetrized projections
    bool do_not_use_symproj;

    /// File containing surface mask
    FileName fn_surface_mask;

    /// Selection file with all images to process
    FileName fn_sel;

    /// Goldmask
    double goldmask;

    /// Shifted tomograms
    bool shiftedTomograms;

    /// Root of output filenames
    FileName fn_root;

    /// Grid volume as initial guess
    FileName fn_start;

    /// Stop after this number of images, if 0 then don't use
    int stop_at;

    /// Known volume. If -1, not applied.
    double known_volume;

    /// Selection file with all images to process
    FileName fn_ctf;

    /// Apply unmatched projectors to correct for the CTF
    bool unmatched;

    /** Ray length.
        Basis functions are taken into account only if their distance
        to the projection plane is smaller than this value.
        This value is expressed in basis units. Set it to
        -1 to disable it. Set it to 1 (1 basis away maximum) to
        interpolate a set of planes). */
    double ray_length;

    /// Apply shifts stored in the headers of the 2D-images
    bool apply_shifts;

    /// Apply positivity constraint
    bool positivity;

    /// Print system matrix
    bool print_system_matrix;

    /// Is this a crystal 0 means NO 1 YES
    bool is_crystal;

    /// Variability analysis
    bool variability_analysis;

    /// Noisy reconstruction
    bool noisy_reconstruction;
	
	/// Only for internal purposes, MUST be set when running MPI.
	bool using_MPI;
	
	/// Number of threads to use. Can not be different than 1 when using MPI.
	int threads;

#define TELL_IV                    0x100
#define TELL_ONLY_SYM              0x80
#define TELL_USE_INPUT_BASISVOLUME 0x40
#define TELL_SHOW_ERROR            0x20
#define TELL_MANUAL_ORDER          0x10
#define TELL_SAVE_AT_EACH_STEP     0x8
#define TELL_SAVE_INTERMIDIATE     0x4
#define TELL_SAVE_BASIS            0x2
#define TELL_STATS                 0x1
    /** Debugging level.
        This is a bit valued field, you can set the following bits
        \\TELL_IV: Show intermidiate images if saved,
        \\         Show the reconstructed volume each time the progress
                   bar is updated.
        \\TELL_ONLY_SYM: Skip all the extra projections created using the
           samle symmetry
        \\TELL_USE_INPUT_BASISVOLUME: This flag causes the program to
            not resizing the gridvolume. The same as the input one is
     used as starting point for the algorithm
        \\TELL_SHOW_ERROR: The program will show the error for each projection
        \\TELL_MANUAL_ORDER: The program will ask the number of the
           following projection to process instead of using the
           default perpendicular order
        \\TELL_SAVE_AT_EACH_STEP: At each step (an iteration is
           compound of several steps) the following files are
           written to disk
           @code
              PPPdiff          --> Difference between theoretical and real projections
              PPPtheo          --> Theoretical projection
              PPPread          --> Real projection
              PPPcorr          --> Correction image applied
              PPPbasis.basis   --> Reconstructed volume in the used basis
              PPPvol.vol       --> Reconstructed volume in voxels
           @endcode
        \\TELL_SAVE_INTERMIDIATE: At each iteration a voxel (and possibly
           a basis volume (if the TELL_SAVE_BASIS flag is set) is stored with
           the names fn_root"it"it.vol and .basis (for instance,
           w0001it00.vol and w0001it00.basis).
        \\TELL_SAVE_BASIS: Save basis volume at the end and in the
           intermidiate iterations (if TELL_SAVE_INTERMIDIATE is set).
        \\TELL_STATS: Show image and volume (only of the first grid, of
           the 1, 2 or 4 possible subgrids) statistics.*/
    int tell;

    /// Frequency for saving intermidiate
    int save_intermidiate_every;

    /// Name of file for improved control in parallel jobs
    FileName fn_control;
    //@}

    /* ART Side information ................................................. */
    /** ART Side Information
        The Side information is useful information that needn't be computed
        more than once at the beginning from the ART parameters and that
        is used all over the program. */
    //@{
    /// A list with the symmetry matrices
    SymList         SL;

    /// Projection X dimension
    int             projXdim;

    /// Projection Y dimension
    int             projYdim;

    /// File handler for the history file
    std::ofstream        *fh_hist;

    /// Array with all the sorting information for each projection
    Recons_info     *IMG_Inf;

    /// Order in which projections will be presented to algorithm
    Matrix1D<int>   ordered_list;

    /// Total number of images to process (taking symmetries into account)
    int             numIMG;

    /// Number of different images (without symmetries)
    int             trueIMG;

    /** Volume deformation matrix.
        Samples stored in the basis volume really relate to real space using
        this matrix. This means that a sample which is at index (0,1,0) can
        be placed at a space position given by D*(0,1,0)*grid_size, which could
        be for example (0,1.05,0)*grid_size. This is really useful for crystals.
        If you don't use it, set it to NULL.
        This matrix passes from a to aint, and from b to bint.
        @code
        aint = Dinv*a; a=D*aint;
        bint = Dinv*b; b=D*bint;
        @endcode
        */
    Matrix2D<double> *D;
    /// Just the inverse of D
    Matrix2D<double> *Dinv;

    /** Surface mask.
        The volume is supposed to be 0 where the mask is 1. */
    VolumeXmipp *surface_mask;

    /** POCS frequency.
        POCS restrictions are imposed every (this value) projections.
        By default, 1*/
    int POCS_freq;

    /** CAV equation count.
        This volume contains the number of equations at which each basis
        is involved */
    GridVolumeT<int> *GVNeq;
    //@}
public:
    /** Generate default values for ART parameters.
        Compulsory parameters are not filled and must be given externally.
        See Manual help (ART) to see which ones are compulsory. */
    void default_values();

    /** Read parameters from a command line.
        This function reads the parameters from a command line
        defined by argc and argv. An exception might be thrown by any
        of the internal conversions, this would mean that there is
        an error in the command line and you might show a usage message. */
    void read(int argc, char **argv);

    /** Read parameters from a file.
        An exception is thrown if the file cannot be open */
    void read(const FileName &fn);

    /** Usage message.
        This function shows the way of introdustd::cing the most common parameters. */
    void usage();

    /** Full Usage message.
        This function shows the way of introdustd::cing these parameters. */
    void usage_more();

#define BASIC 0
#define FULL  1
    /** Produce Initial and Side information for ART.
        This function computes from the ART parameters things like
        the pojection size, projection order, symmetry matrices
        list, history handler (don't forget to close it at the end), number
        of images, basis side information and initial basis volume.

        Note: It is supposed that time has been previously configured with
        time_config().

        The info level takes the following values:
         \\ BASIC: Generate basis side information
    \\ FULL: Generate all the rest needed values.

        The rank is a number idetifying the parallel process. If -1 then
        the algorithm is sequential. If 0 then it is the root process.
        */
    void produce_Side_Info(GridVolume &vol_basis0, int level = FULL, int rank = -1);

    /** Compute CAV weights.
        The weights are stored in the GVNeq within this object. If the
        debugging level is greater than 0 then a progress bar is shown
        and the number of equations and unknowns are shown at the end.
        Otherwise nothing is printed (this is the suggested debugging level
        for parallel processing). */
    void compute_CAV_weights(GridVolume &vol_basis0,
                             int numProjs_node, int debug_level = 0);

    /** Lambda for iteration n (first one is iteration 0).
        If the iteration requested is greater than the number of lambdas
        provided then the last lambda in the list is returned. An exception
        is thrown if there are no lambdas in the list. */
    double lambda(int n)
    {
        int imax = XSIZE(lambda_list);
        if (imax == 0)
            REPORT_ERROR(1, "Basic_art: There are no lambdas\n");
        if (n >= imax) return lambda_list(imax -1);
        else         return lambda_list(n);
    }

    /** Kappa (WLS) for iteration n (first one is iteration 0).
        If the iteration requested is greater than the number of lambdas
        provided then the last lambda in the list is returned. An exception
        is thrown if there are no lambdas in the list. */
    double kappa(int n)
    {
        int imax = XSIZE(kappa_list);
        if (imax == 0)
            REPORT_ERROR(1, "Basic_art: There are no kappas\n");
        if (n >= imax) return kappa_list(imax -1);
        else         return kappa_list(n);
    }

    /** Returns X dimension for projections under use. */
    int ProjXdim();

    /** Returns Y dimension for projections under use. */
    int ProjYdim();

};

/** Sort projections orthogonally.
   This function sorts a number of images given by numIMG, whose information
   about their Euler angles are in IMG_inf, into an ordered list which
   gives the indexes. First an image is chosen randomly from the whole
   set. Then all images are compared to the first one, and the most
   perpendicular one is chosen. The remaining set of images are compared
   to this two images and the most perpendicular one to the former two
   is chosen, and so on until no image is left in the set.

   If the result in ordered list is 4, 70, 54, 203, 1, 0, ... it means
   that the first image is the number 4, then goes the 70, then the 54, ...

   If N!=-1 then the product is done only with the last N images. A very
   useful value is N=2*/
void sort_perpendicular(int numIMG, Recons_info *IMG_Inf,
                        Matrix1D<int> &ordered_list, int N = 2);

/** No projection sorting at all.
    This function directly returns the same order as in the selection file */
void no_sort(int numIMG, Matrix1D<int> &ordered_list);

/** Randomize the projections.
   This function sorts randomly a number of images given by numIMG. */
void sort_randomly(int numIMG, Matrix1D<int> &ordered_list);

/** Write first part of ART history.
    This function writes all ART parameters, projection angles, symmetry
    matrices, and grid structure in the history handler provided by
    the side information. At the end the routine makes a call to the
    operator << of the Extra_ART_Parameters to show the specific part
    of the History.

    Basic_ART_Parameters is not constant since things are written in
    \ref Basic_ART_Parameters::fh_hist.*/
template <class Extra_ART_Parameters>
void Basic_ART_Init_history(Basic_ART_Parameters &prm,
                            const Extra_ART_Parameters &eprm, const GridVolume &vol_basis0);

/** Perform all ART iterations.
    This function performs the iterations according to the ART parameters,
    it needs the side information to be fully computed. It throws
    a lot of information to the screen and to the history file (side.fh_hist),
    specially this one must exist.

    The GridVolume must have as input an initial guess for the solution,
    and when this function finishes, it contains the final solution volume
    in basis.

    The rank is the identification number of the process running this function.
    If it is -1, the function is run in seuqential mode. If it is 0, then
    it is the root process.

    See the \ref Basic_ART_Parameters for more information
    about how to generate the iterations.
*/


template <class Extra_ART_Parameters>
void Basic_ART_iterations(Basic_ART_Parameters &prm,
                          Extra_ART_Parameters &eprm, GridVolume &vol_basis, int rank = -1);

/** Main Routine for ART.
    Given any set of Art parameters, this function returns the voxel
    volume which is the solution for this set of projections.
    No initialisation is needed on vol_voxels and vol_basis, they are resized
    to have the same size as the input projections. All output files
    are generated as if the ART program had been called. */
template <class Extra_ART_Parameters>
void Basic_ROUT_Art(Basic_ART_Parameters &prm,
                    Extra_ART_Parameters &eprm, VolumeXmipp &vol_voxels,
                    GridVolume &vol_basis);

/** Reconstruction information.
   This structure contains information for all projections which are
   going to participate in the reconstruction.
   This structure has also got
   information for the symmetry implementation. If there is any symmetry
   then an entry in this table is created using the same projection name
   but different symmetry matrices (only the matrix index is annotated
   in this structure). The Euler angles stored for the symmetrized image
   are the final ones, ie, the original Euler angles symmetrized according
   to the symmetry matrix. Then the symmetry identificator kept in this
   structure is only used to keep some track of what matrix was used to
   symmetrize.
*/
struct Recons_info
{
    /// Projection filename
    FileName fn_proj;
    /// CTF filename
    FileName fn_ctf;
    /// Rotational angle
    float  rot;
    /// Tilting angle
    float  tilt;
    /// Psi angle
    float  psi;
    /** Symmetry number.
        This number express to which symmetry matrix this projection
        is related to (-1: without symmetry, 0: using symmetry matrix 0,
        1: using symmetry matrix 1 ...) */
    int    sym;
    /** Random seed.
        For the reconstruction of pure noise for VSSNR, all images
        coming from the same projection by symmetry relationships should
        have the same random seed. */
    int    seed;
};

/** Build from a Selection File and a Symmetry List.
    The result is stored in the Recons_info array which should point
    to NULL when it is not initialized. */
void build_recons_info(SelFile &selfile, SelFile &selctf,
    const FileName &fn_ctf, const SymList &SL, Recons_info * &IMG_Inf,
    bool do_not_use_symproj);

/* ------------------------------------------------------------------------- */
/** Variability structure */
class VariabilityClass
{
public:
    typedef enum {VAR_none, VAR_measuring, VAR_analyzing} t_VAR_status;
    t_VAR_status VAR_state;
    int Zoutput_volume_size;
    int Youtput_volume_size;
    int Xoutput_volume_size;
    Basic_ART_Parameters *prm;

    /// Vector of training vectors
    std::vector < Matrix3D<double> > VA;

    /// Number of updates so far
    int N;

    /// Constructor
    VariabilityClass(Basic_ART_Parameters *_prm,
                     int _Zoutput_volume_size, int _Youtput_volume_size,
                     int _Xoutput_volume_size);

    /** Start a new ART iteration. */
    void newIteration();

    /** Update data with a new volume.
        The update volume is set to zeros after this function */
    void newUpdateVolume(GridVolume *ptr_vol_out, Projection &read_proj);

    /** Finish analysis. */
    void finishAnalysis();
};

/* ------------------------------------------------------------------------- */
/** POCS structure */
class POCSClass
{
public:
    typedef enum {POCS_measuring, POCS_use, POCS_lowering, POCS_N_measure,
                  POCS_N_use} t_POCS_status;
    t_POCS_status POCS_state;
    double POCS_avg;
    double POCS_stddev;
    double POCS_min;
    double POCS_max;
    double POCS_mean_error;
    double POCS_max_error;
    double POCS_global_mean_error;
    int POCS_freq;
    int POCS_i;
    int POCS_vec_i;
    int POCS_used;
    int POCS_N;
    int Zoutput_volume_size;
    int Youtput_volume_size;
    int Xoutput_volume_size;
    bool apply_POCS;
    Matrix1D<double> POCS_errors;
    Basic_ART_Parameters *prm;

    /// Constructor
    POCSClass(Basic_ART_Parameters *_prm,
              int _Zoutput_volume_size, int _Youtput_volume_size,
              int _Xoutput_volume_size);

    /// Start New ART iteration
    void newIteration();

    /// Start new Projection
    void newProjection();

    /// Apply
    void apply(GridVolume &vol_basis, int it, int images);
};
//@}

#include "basic_art.inc"

#endif
