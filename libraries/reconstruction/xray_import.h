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

#include <data/xmipp_funcs.h>
#include <data/multidim_array.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>
#include "data/xmipp_hdf5.h"
#include <data/xmipp_threads.h>

///@defgroup XrayImport Xray import
///@ingroup ReconsLibrary
//@{
/** Xray import parameters. */
class ProgXrayImport: public XmippProgram
{
public:
    /// Input directory
    FileName fnInput;
    /// Input Flatfield directory
    FileName fnFlat;
    /// Output directory
    FileName fnRoot;
    /// Output Stack file
    FileName fnOut;
    /// Number of pixels to crop from each side. Set to 0 for no cropping
    int cropSize;
    /// Number of threads
    int thrNum;
    /// bad Pixel filter Mask
    FileName fnBPMask;
    /// Flag to apply flat field correction
    bool flatFix;
    /// Flag to apply dark field correction
    bool darkFix;
    /// Flag to apply log filter
    bool logFilt;
    /// Xray microscopy data origin;
    enum DataSource
    {
        NONE,
        MISTRAL,
        BESSY
    } dSource;
    /// Index number used in Bessy tomograms and flatfield images
    size_t tIni, tEnd, fIni, fEnd;
    /// hdf5 file handler
    XmippH5File H5File;
    double* vCBeam, vExpTime, vslitWidth;
    Matrix1D<double> expTimeArray, cBeamArray, slitWidthArray, anglesArray;

    /// List of input images
    MetaData inMD;
    /// List of output images
    MetaData outMD;

    // Intermediate results
    Image<double> IavgFlat;
    Image<double> IavgDark;
    Image<char>   bpMask;
    std::vector<FileName> filenames;
    std::vector<size_t> objIds;


    // Variables for the threads
    ParallelTaskDistributor *td;
    ThreadManager           *tm;
    Mutex                    mutex;

public:

    /// Constructor
    void init();

    /// Read argument from command line
    void readParams();

    /// Define params
    void defineParams();

    /// Show
    void show() const;

    /** Really import*/
    void run();

    /// Read an image and crop
    void readAndCrop(const FileName &fn, Image<double> &I) const;

    /// Read geometrical info
    void readGeoInfo(const FileName &fn, MDRow &rowGeo) const;

    /// Read related data to normalize the tomogram
    void readCorrectionInfo(const FileName &fn, double &currentBeam,
                            double &expTime, double &slitWidth) const;

    /** Get the darkfield for a directory.
     *  In case there is no darkfield a message is shown and the output image
     *  is empty. */
    void getDarkfield(const FileName &fnDir, Image<double> &IavgDark);

    /** Get the corrected average of a directory.
     *  If there is a darkfield, a file called fnRoot+"_"+fnDir+"_darkfield.xmp"
     *  is saved.
     */
    void getFlatfield(const FileName &fnDir, Image<double> &Iavg);
public:
};
//@}
#endif
