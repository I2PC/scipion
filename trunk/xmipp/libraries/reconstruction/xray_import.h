/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2010)
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
#ifndef _PROG_XRAY_IMPORT
#define _PROG_XRAY_IMPORT

#include <data/funcs.h>
#include <data/multidim_array.h>
#include <data/image.h>

///@defgroup XrayImport Xray import
///@ingroup ReconsLibrary
//@{
/** Xray import parameters. */
class Prog_xray_import_prm
{
public:
    /// Input directory
    FileName fnDirData;
    /// Input Flatfield directory
    FileName fnDirFlat;
    /// Output directory
    FileName fnRoot;
    /// Number of pixels to crop from each side. Set to 0 for no cropping
    int cropSize;

public:
    /// Read argument from command line
    void read(int argc, char **argv);

    /// Show
    friend std::ostream & operator << (std::ostream &out,
        const Prog_xray_import_prm &prm);

    /// Usage
    void usage() const;

    // Read an image and correct
    void readAndCorrect(const FileName &fn, Image<double> &I);

    /** Get the darkfield for a directory.
     *  In case there is no darkfield a message is shown and the output image
     *  is empty. */
    void getDarkfield(const FileName &fnDir, Image<double> &IavgDark);

    /** Get the corrected average of a directory.
     *  If there is a darkfield, a file called fnRoot+"_"+fnDir+"_darkfield.xmp"
     *  is saved.
     */
    void getCorrectedAvgOfDirectory(const FileName &fnDir, Image<double> &Iavg);

    /** Really import*/
    void run();
};
//@}
#endif
