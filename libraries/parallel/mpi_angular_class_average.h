/***************************************************************************
 * Authors:     AUTHOR_NAME (aerey@cnb.csic.es)
 *
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

#ifndef MPI_ANGULAR_CLASS_AVERAGE_H_
#define MPI_ANGULAR_CLASS_AVERAGE_H_

//mpirun -np 5 xmipp_mpi_angular_class_average --nJobs 70
#include <parallel/xmipp_mpi.h>
#include <data/xmipp_funcs.h>
#include <data/xmipp_program.h>
#include <data/metadata.h>

#include <data/xmipp_fftw.h>
#include <data/args.h>
#include <data/xmipp_image.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/polar.h>
#include <data/basic_pca.h>
#include <data/sampling.h>

//Tags already defined in xmipp
//#define TAG_WORK                     0
//#define TAG_STOP                     1
#define TAG_I_FINISH_WRITTING        12
#define TAG_MAY_I_WRITE              13
#define TAG_YES_YOU_MAY_WRITE        14
#define TAG_DO_NOT_DARE_TO_WRITE     15
#define TAG_I_AM_FREE                16

#define lockWeightIndexesSize 5
#define index_lockIndex 0
#define index_weight 1
#define index_weights1 2
#define index_weights2 3
#define index_ref3d 4

#define ArraySize 8
#define index_DefGroup 0
#define index_2DRef 1
#define index_3DRef 2
#define index_Order 3
#define index_Count 4
#define index_jobId 5
#define index_Rot 6
#define index_Tilt 7

#define AVG_OUPUT_SIZE 10

#define split0 0
#define split1 1
#define split2 2

class MpiProgAngularClassAverage : public XmippMpiProgram
{
public:
    // Seed
    unsigned int master_seed;

    // Metadata with the list of jobs
    MetaData mdJobList;

    /** Input and library docfiles */
    MetaData         DF, DFlib, DFscore;
    /** metadata with classes which have experimental images applied to them */
    MetaData         DFclassesExp;
    /** Output rootnames */
    FileName         fn_out, fn_out1, fn_out2, fn_wien, fn_ref;
    /** Column numbers */
    std::string      col_select;
    /** Upper and lower absolute and relative selection limits */
    double           limit0, limitF, limitRclass, limitRper;
    /** Flags wether to use limit0, limitF and limitR selection */
    bool             do_limit0, do_limitF, do_limitR0class, do_limitRFclass, do_limitR0per, do_limitRFper;
    /** Flag whether to apply mirror operations. By default set to True */
    bool             do_mirrors;
    /** Flag whether also to write out class averages of random halves of the data */
    bool             do_split;
    /** Image dimensions before and after padding (only for Wiener correction) */
    size_t           paddim;
    /** Padding factor */
    double           pad;
    /** One empty image with correct dimensions */
    Image<double>    Iempty;
    /** Do write images assigned to each class */
    bool             do_save_images_assigned_to_classes;
    /** Add output to existing files */
    bool             do_add;
    /** Perform PCA sorting to obtain the average classes */
    bool             do_pcaSorting;
    /** Wiener filter image */
    MultidimArray<double> Mwien;
    /** Selfiles containing all class averages */
    MetaData         SFclasses, SFclasses1, SFclasses2;

    /** Re-alignment of classes */

    /** Input file */
    FileName inFile, refGroup;
    /** Inner and outer radius for rotational alignment */
    int Ri, Ro;
    /** Number of iterations */
    int nr_iter;
    /** Convergence criterion */
    double eps;
    /** Search shift (shifts larger than this will be set to 0)*/
    double max_shift;
    /** Maximum allowed shift in last iteration (shifts larger than this will be set to 0)*/
    double max_shift_change, max_psi_change;
    /** transformers for all rings */
    Polar_fftw_plans global_plans;
    RotationalCorrelationAux rotAux;
    MultidimArray<double> corr;

    /** number of Ctf groups */
    int ctfNum;
    /** Number of 3D references */
    int ref3dNum;
    /** Image dimensions */
    size_t Xdim, Ydim, Zdim, Ndim;

    /** Divide the job in this number block with this number of images */
    size_t mpi_job_size;

    //Lock structure
    MultidimArray<bool> lockArray;
    MultidimArray<double> weightArray;
    MultidimArray<double> weightArrays1;
    MultidimArray<double> weightArrays2;

    MpiProgAngularClassAverage();

    MpiProgAngularClassAverage(int argc, char **argv);

    /** Redefine read */
//    void read(int argc, char** argv);

    void readParams();

    void defineParams();

    void run();

    /** Process a job list (ref3d - ctfGroup - ref2d)
         */
    void mpi_process_loop(double * Def_3Dref_2Dref_JobNo);

    /** Process a single job (ref3d - ctfGroup - ref2d)
         */
    void mpi_process(double * Def_3Dref_2Dref_JobNo);

    /** Initialize
         */
    void mpi_produceSideInfo();

    /** Initialize variables, create infrastructure for mpi job submission, and delete output files.
         */
    void mpi_preprocess();

    /** Read input metadata and filter following the user define constraints
         */
    void filterInputMetadata();

    /** Save discarded images
         */
    void saveDiscardedImages();

    /** Initialize file names.
         */
    void initFileNames();

    /** Get file and stack dimentions, and number of 3d references and number of defocus groups.
         */
    void initDimentions();

    /** Initialize weights matrices.
         */
    void initWeights();

    /** Delete output files if they exist and init stacks.
         */
    void initOutputFiles();

    /**
         */
    void mpi_postprocess();

    /** Compute list of parallel jobs to be performed
             */
    void createJobList();

    /** Write output files
         */
    void mpi_write(
        size_t dirno,
        int ref3dIndex,
        Image<double> avg,
        Image<double> avg1,
        Image<double> avg2,
        MetaData SFclass,
        MetaData SFclass1,
        MetaData SFclass2,
        MetaData SFclassDiscarded,
        double w1,
        double w2,
        double old_w,
        double old_w1,
        double old_w2);

    /** Block output file to avoid concurrent writing
         */
    void mpi_writeController(
            size_t dirno,
            Image<double> avg,
            Image<double> avg1,
            Image<double> avg2,
            MetaData SFclass,
            MetaData SFclass1,
            MetaData SFclass2,
            MetaData SFclassDiscarded,
            MetaData _DF,
            double w1,
            double w2,
            int lockIndex);

    /** Called by mpi_write does the actual writing
         */
    void mpi_writeFile(
    		Image<double> avg,
    		size_t dirno,
    		FileName fileNameStk,
    	    double w_old);

    /**
         */
    void getPolar(
    		MultidimArray<double> &img,
    		Polar<std::complex <double> > &fP,
    		bool conjugated=false,
    		float xoff = 0.,
    		float yoff = 0.);

    /**
         */
    void reAlignClass(
    		Image<double> &avg1,
    		Image<double> &avg2,
    		MetaData &SFclass1,
    		MetaData &SFclass2,
    		std::vector<Image<double> > imgs,
    		std::vector<int> splits,
    		std::vector<int> numbers,
    		size_t dirno,
    		double * my_output);

    /** Apply Wiener filter
         */
    void applyWienerFilter(MultidimArray<double> &img);

};

#endif /* MPI_ANGULAR_CLASS_AVERAGE_H_ */



