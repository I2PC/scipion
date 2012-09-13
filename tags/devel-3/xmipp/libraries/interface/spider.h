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

#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
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
void translate_to_Spider_sel(SelFile &SF_in, DocFile &DF_out, bool new_style);

/** Extract angles from a SelFile and store them in a DocFile.
    You can specify the order of the angle extraction by default
    (rot, tilt, psi).
    An exception is thrown if the angles are not one of these.*/
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
void rename_for_Spider(SelFile &SF_in, SelFile &SF_out, const FileName &fn_root,
                       const FileName &out_ext);

//@}

#endif
