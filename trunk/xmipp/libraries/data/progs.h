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

#ifndef PROGS_H
#define PROGS_H

#include "image.h"
#include "volume.h"

/** @defgroup Programs Programs.
 *
 * General program routines.
 *
 * Here you will find some routines which might help you to write programs
 * with a common structure where only a small part changes. These routines
 * make the heavy work for you and you only have to supply the changing
 * part from one to another.
 */

/** Virtual program parameters class.
 * @ingroup Programs
 *
 * Program parameters must be a class inheriting from this class.
 */
class Prog_parameters
{
public:
    /// Output extension
    FileName oext;

    /// Output root
    FileName oroot;

    /// Output file
    FileName fn_out;

    /// Input file
    FileName fn_in;

    /// This flag for application of the transformation as stored in the header
    bool apply_geo;

    /// Use this flag for not writing at every image
    bool each_image_produces_an_output;

    /// Use this flag for not producing a time bar
    bool allow_time_bar;

    /// Empty constructor
    Prog_parameters()
    {
        oroot = oext = fn_out = fn_in = "";
        each_image_produces_an_output = true;
        allow_time_bar = true;
        apply_geo = false;
    }

    /// Read the basic parameters defined for this class
    virtual void read(int argc, char** argv);

    /// Show these parameters
    virtual void show();

    /// Usage
    virtual void usage();

    /// Get input size
    void get_input_size(int& Zdim, int& Ydim, int& Xdim);

    /// Get the number of images to process
    int get_images_to_process();
};

#define IMAGE2IMAGE 1
#define FOURIER2FOURIER 2
#define IMAGE2FOURIER 3
#define FOURIER2IMAGE 4
#define IMAGE2FILE 5
#define FILE2FILE 6

/** Main for SelFiles and individual images without parameters.
 *
 * This routine implements the main of a program which reads either an
 * image/volume or a selfile, processes each image/volume separately
 * (this is the part you have to supply, together with the extra parameters
 * needed by the program), for instance, add some noise, and then save the
 * results. The output images can be the same input ones, or new ones (in case
 * of an image, the user must supply the output name, and for selection files,
 * he must supply an output extension). If an error occurs, the routine shows
 * the error message and exits 1, ie, the program is aborted.
 */
void SF_main(int argc,
             char** argv,
             Prog_parameters* prm,
             void* process_img,
             void* process_vol,
             int operation_mode = IMAGE2IMAGE);
#endif
