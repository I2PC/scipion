/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Jose Roman Bilbao-Castro (jrbcast@ace.ual.es)
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

#ifndef __RECONSTRUCT_FOURIER_H
#define __RECONSTRUCT_FOURIER_H

#include <iostream>
#include <data/fftw.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/projection.h>
#include <data/threads.h>

#include <reconstruction/projection.h>
#include <reconstruction/directions.h>
#include <reconstruction/symmetrize.h>
#include <reconstruction/blobs.h>
#define BLOB_TABLE_SIZE 5000
#define MINIMUMWEIGHT 0.001
#define ACCURACY 0.001

#define EXIT_THREAD 0
#define PROCESS_IMAGE 1
#define PROCESS_WEIGHTS 2
#define PRELOAD_IMAGE 3

/**@defgroup Fourier reconstruction reconstruct_fourier (Fourier reconstruction)
   @ingroup ReconsLibraryPrograms */
//@{
class Prog_RecFourier_prm;
static pthread_mutex_t mutexDocFile= PTHREAD_MUTEX_INITIALIZER;

struct ImageThreadParams
{
	int myThreadID;
	Prog_RecFourier_prm * parent;
        Matrix2D< std::complex<double> > *paddedFourier;
	Matrix2D< std::complex<double> > *localPaddedFourier;
        Matrix2D<double> * symmetry;
        int read;
        Matrix2D<double> * localAInv;
        int imageIndex;
        DocFile * docFile;
        double weight;
        double localweight;

};

/** Fourier reconstruction parameters. */
class Prog_RecFourier_prm
{
public:
    /** Filenames */
    FileName fn_out, fn_sym, fn_sel, fn_doc, fn_control, fn_fsc;

    /** SelFile containing all projections */
    SelFile SF;

    /** DocFile containing all angles */
    DocFile DF;

    /** verbosity flag */
    int verb;

    /** Flag whether to use the weights in the image headers */
    bool do_weights;

    /** Projection padding Factor */
    double padding_factor_proj;

    /** Volume padding Factor */
    double padding_factor_vol;

    /** Sampling rate in Angstroms/pixel */
    double sampling_rate;

    /// Max resolution in Angstroms
    double maxResolution;

    /// Number of iterations for the weight
    int NiterWeight;
    
    /// Number of threads to use in parallel to process a single image
    int numThreads;
    
    /// IDs for the threads
    pthread_t * th_ids;
    
    /// Contains parameters passed to each thread
    ImageThreadParams * th_args;
    
    /// Tells the threads what to do next
    int threadOpCode;
    
    /// Number of rows already processed on an image
    int rowsProcessed;
    
    /// Defines what a thread should do
    static void * processImageThread( void * threadArgs );
    
    /// Controls mutual exclusion on critical zones of code
    pthread_mutex_t workLoadMutex;
    
    /// To create a barrier synchronization for threads
    barrier_t barrier;
    
    /// A status array for each row in an image (processing, processed,etc..)
    int * statusArray;
    
    /// How many image rows are processed at a time by a single thread.
    int thrWidth;
	
public: // Internal members
    // Size of the original images
    int imgSize;
    
    // Padded volume size
    int volPadSizeX;
    int volPadSizeY;
    int volPadSizeZ;

    // Column numbers in the docfile
    int col_rot, col_tilt, col_psi, col_xoff, col_yoff, col_flip, col_weight;

    // Table with blob values
    Matrix1D<double> blob_table, Fourier_blob_table;

    // Inverse of the delta and deltaFourier used in the tables
    double iDelta, iDeltaFourier;

    // Maximum interesting resolution squared
    double maxResolution2;

    // Definition of the blob
    struct blobtype blob; 

    // vector with R symmetry matrices
    std::vector <Matrix2D<double> > R_repository;

    // Fourier transformer for the volume
    XmippFftw transformerVol;

    // Fourier transformer for the images
    XmippFftw transformerImg;

    // An alias to the Fourier transform in transformerVol
    Matrix3D< std::complex<double> > VoutFourier;

    // Volume of Fourier weights
    Matrix3D<double> FourierWeights;

    // Volume of Fourier weights convolved with the kernel
    Matrix3D<double> FourierWeightsConvolved;

    // Kernel used in Fourier
    Matrix3D<double> kernel;

    // Padded image
    Matrix2D<double> paddedImg;

    // Output volume
    VolumeXmipp Vout;
public:
    /// Read arguments from command line
    void read(int argc, char **argv);

    /// Show
    void show();

    /// Usage
    void usage();

    /// Produce side info: fill arrays with relevant transformation matrices
    void produce_Side_info() ;

    /// Get angles (either from reading the header or from a docfile)
    void get_angles_for_image(const FileName &fn, double &rot, double &tilt,
        double &psi, double &xoff, double &yoff, double &flip, double &weight, DocFile * docFile);

    /// Main Loop 
    void run();
    
    void finishComputations( FileName out_name );
    
    /// Process one image
    void processImages( int firstImageIndex, int lastImageIndex, bool saveFSC=false ); //const FileName &fn_img);

    /// Correct weight
    void correctWeight();
};
//@}
#endif
