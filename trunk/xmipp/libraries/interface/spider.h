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
/*****************************************************************************/
/* INTERACTION WITH SPIDER                                                   */
/*****************************************************************************/

#ifndef _XMIPP_SPIDER_HH
#define _XMIPP_SPIDER_HH

#include <data/metadata.h>

#include <data/funcs.h>
//include <data/selfile.h>
#include <data/docfile.h>
#include <data/volume.h>
#include <data/image.h>

/**@defgroup SpiderInterface Spider
   @ingroup InterfaceLibrary */
//@{
/** Generate a Spider "count" file.
    This function returns a DocFile with (1, 2, 3, ..., imax) */
void generate_Spider_count(int imax, DocFile &DF_out);

/** From a Xmipp selfile to Spider selfile.
    Comments are lost. A Spider Selfile is created in which the -1
    of the Xmipp Selfiles are translated into 0.

    Set new_style to produce the new style of Spider selfiles.*/
void translate_to_Spider_sel(MetaData &SF_in, DocFile &DF_out, bool new_style);

/** Extract angles from a SelFile and store them in a DocFile.
    You can specify the order of the angle extraction by default
    (rot, tilt, psi).
    An exception is thrown if the angles are not one of these.*/
void extract_angles(MetaData &SF_in, DocFile &DF_out,
                    const std::string &ang1 = "rot", const std::string &ang2 = "tilt",
                    const std::string &ang3 = "psi");
void extract_angles(SelFile &SF_in, DocFile &DF_out,
                    const std::string &ang1 = "rot", const std::string &ang2 = "tilt",
                    const std::string &ang3 = "psi");

#ifdef NEVERDEFINED
/** Extract angles from a Docfile and store them in a SelFile.
    You can specify the order of the angle extraction by default
    (rot, tilt, psi).
    An exception is thrown if the angles are not one of these.*/
void write_angles(SelFile &SF_out, DocFile &DF_in,
                  const std::string &ang1 = "rot", const std::string &ang2 = "tilt",
                  const std::string &ang3 = "psi");
#endif

/** Rename ACTIVE images in a selfile.
    The images are renamed with the fn_root given, consecutive numbers
    starting at 1, and the extension given. If no extension is given the
    same one as the input one is used. The correspondence between files
    is shown in stdout. A selFile is returned with the new set of images.
    The renaming is done by calls to cp in the Operating System, so there
    is no restriction on th einput names. */
void rename_for_Spider(MetaData &SF_in, MetaData &SF_out, const FileName &fn_root,
                       const FileName &out_ext);

/** Create empty Spider file.
    Creates a zero filled spider file with the desired dimension. */
void create_empty_Spider_file(const FileName &fn, int Zdim, int Ydim,
                              int Xdim, bool reversed = false, size_t block_size = 102400);

/** 3D Radon transform.
    Creates the 3D radon transform of a volume via Spider.
    An exception is thrown if some intermidiate files (b01.vol, superfeo.vol,
    superfeo2.vol) cannot be created.

    The Delta_rot and tilt are the sampling rates in the rot and tilt space
    measured in degrees.

    The output_size is the size of the output Radon transform. If it is -1
    then it is 1.5*XSIZE(Vol_in).

    The outputis written to file although it is a VolumeXmipp.
**/
void radon_transform(VolumeXmipp &V_in, const FileName &fn_out,
                     double Delta_rot = 2, double Delta_tilt = 2, int output_size = -1);

/** 2D Radon transform.
    Creates the 2D radon transform of an image via Spider.
    An exception is thrown if some intermidiate files (b01.xmp, superfeo.xmp,
    superfeo2.xmp) cannot be created.

    The Delta_ang is the sampling rate in the angular space
    measured in degrees.

    The output_size is the size of the output Radon transform. If it is -1
    then it is 1.5*XSIZE(Vol_in).

    The output is written to file as it cannot be held by
    any of the Xmipp classes.
**/
void radon_transform(ImageXmipp &I_in, const FileName &fn_out,
                     double Delta_ang = 2, int output_size = -1);

/** Fourier Radon transform.
    Makes the Fourier-Radon transform of the volume or image supplied as fn_in.
    The cutoff_freq is normalized to 0.5, while the Fermi
    temperature is a smoothing factor for the filtering.

    An exception is thrown if b01.fft, superfeo.fft or superfeo2.fft cannot be created.
*/
void Fourier_transform_of_Radon_transform(const FileName &fn_in,
        const FileName &fn_out, double cutoff_freq,
        double Fermi_temperature = 0.2);
#ifdef DEPRECATED
/** Angular_refinement via Radon.
    The angular refinement process via the Radon transform is performed.
    All refinements are considered as subsearches within the given range.
    The default range (0,360),(0,90),(0,360) covers the whole
    projection space.

    Input files must be Fourier Radon transforms of the input images.
    The files "peak?????" are created and destroyed again.

    An exception is thrown if the file b01."ext" cannot be created.
*/
void Angular_refinement_Radon(const FileName &fn_vol, const FileName &fn_sel,
                              const FileName &fn_report,
                              double rot0 = 0, double rotF = 360, double rot_step = 2,
                              double tilt0 = 0, double tiltF = 90, double tilt_step = 2,
                              double psi0 = 0, double psiF = 360, double psi_step = 2,
                              double max_shift = 2);

/** Angular refinement via Projection Matching. fn_ext is computed as the
    extension of the volume. The files refangles.(fn_ext) and projlist.(fn_ext)
    are created with the reference angles and the projection list respectively.
    Ideal projections are called ideal****.(fn_ext). The selfile with the
    experimental images (fn_sel) is translated for Spider under the name
    experimentalsel.(fn_ext). The max_shift must be a multiple of the
    shift_step*/
void Angular_refinement_Matching(const FileName &fn_vol, const FileName &fn_sel,
                                 const FileName &fn_report,
                                 double tilt_step = 2,
                                 double max_shift = 2, double shift_step = 1,
                                 double first_ring = 0, double last_ring = -1);

#endif
//@}

#endif
