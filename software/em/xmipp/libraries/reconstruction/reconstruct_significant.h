/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef __RECONSTRUCT_SIGNIFICANT_H
#define __RECONSTRUCT_SIGNIFICANT_H

#include <data/xmipp_program.h>
#include "angular_project_library.h"
#include "volume_initial_simulated_annealing.h"

/**@defgroup ReconstructSignificant Reconstruct multiple volumes analyzing significant correlations
   @ingroup ReconsLibrary */
//@{

/** Significant reconstruction parameters. */
class ProgReconstructSignificant: public XmippProgram
{
public:
    /** Filenames */
    FileName fnIn, fnDir, fnSym, fnInit, fnFirstGallery;

    /** First significance */
    double alpha0;

    /** Last significance */
    double alphaF;

    /** Total number of iterations */
    int Niter;

    /** Keep intermediate volumes */
    bool keepIntermediateVolumes;

    /** Angular sampling */
    double angularSampling;

    /** Maxshift */
    double maxShift;

    /** Minimum tilt */
    double tilt0;

    /** Maximum tilt */
    double tiltF;

    /** Use imed */
    bool useImed;

    /** Strict */
    bool strict;

    /** Neighbourhood in angles */
    double angDistance;

    /** Number of volumes to reconstruct */
    int Nvolumes;

    /** Apply fisher */
    bool applyFisher;

    /** Do reconstruct */
    bool doReconstruct;

    /** Use it for validation */
    bool useForValidation;

    size_t numOrientationsPerParticle;

    bool dontCheckMirrors;


public: // Internal members
    size_t rank, Nprocessors;

    // Metadata with input images and input volumes
    MetaData mdIn;

    // Size of the images
    size_t Xdim;

    // Partial reconstruction metadatas
    std::vector<MetaData> mdReconstructionPartial;

    // Projection matching metadata
    std::vector<MetaData> mdReconstructionProjectionMatching;

    // Set of all correlations
    MultidimArray<double> cc;

    // Set of all weights
    MultidimArray<double> weight;

    // Set of images in the gallery, it should be a metadata but metadatas do not support threads
    std::vector< std::vector<GalleryImage> > mdGallery;

    // Set of input images
    // COSS std::vector<FileName> mdInp;

    // Images
    // COSS Image<double> inputImages;
    std::vector< Image<double> > gallery;
    std::vector< AlignmentTransforms* > galleryTransforms;

	// Current iteration
	int iter;

	// Current alpha
	double currentAlpha;
public:
	/// Empty constructor
	ProgReconstructSignificant();

    /// Read arguments from command line
    virtual void readParams();

    /// Read arguments from command line
    void defineParams();

    /** Show. */
    void show();

    /** Run. */
    void run();

    /// Produce side info: fill arrays with relevant transformation matrices
    void produceSideinfo();

    /// Reconstruct current volume
    void reconstructCurrent();

    /// Generate projections from the current volume
    void generateProjections();

    ///
    void numberOfProjections();

    /// Align images to gallery projections
    void alignImagesToGallery();

    /// Gather alignment
    virtual void gatherAlignment() {}

    /// Synchronize with other processors
    virtual void synchronize() {}
};
//@}
#endif
