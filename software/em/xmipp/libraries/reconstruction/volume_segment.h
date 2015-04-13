/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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
#ifndef _PROG_SEGMENT
#define _PROG_SEGMENT

#include <data/xmipp_funcs.h>
#include <data/multidim_array.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>

///@defgroup VolumeSegment Volume segmentation
///@ingroup ReconsLibrary
//@{
/** Segment parameters. */
class ProgVolumeSegment: public XmippProgram
{
public:
    /// Input volume
    FileName fn_vol;
    /// Segmentation method
    String method;
    /** Desired mass (in voxels), if not given computed from
        the mass and the sampling rate by produce_info. */
    double   voxel_mass;
    /// Desired mass (in Daltons). Not necessary if voxel_mass is provided
    double   dalton_mass;
    /// Desired mass (in aminoacids). Not necessary if voxel_mass is provided
    double   aa_mass;
    /// Sampling rate (in A/pixel). Not necessary if voxel_mass is provided
    double   sampling_rate;
    /// Output mask. If not given it is not written
    FileName fn_mask;
    /// Enable a single threshold measure
    bool     en_threshold;
    /// Threshold
    double   threshold;
    /// Use Otse
    bool otsu;

    /// From here on by Sjors
    // Create probabilistic solvent mask
    bool do_prob;
    /// radius for B.C. Wang-like smoothing procedure
    int wang_radius;

public:
    // Input volume
    Image<double> V;
public:
    /// Read arguments
    void readParams();

    /// Show
    void show() const;

    /// Define parameters
    void defineParams();

    /** Produce side info.
        Read the input volume, and compute the number of voxels
        if not given. An exception is thrown if no way is given to compute
        the voxel mass*/
    void produce_side_info();

    /** Really compute the mask. If a mask name is given then it is
        written to disk.*/
    void segment(Image<double> &mask);

    /** Run */
    void run();
};
//@}
#endif
