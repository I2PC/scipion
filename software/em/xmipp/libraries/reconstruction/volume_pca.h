/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2015)
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
#ifndef _PROG_VOLUME_PCA
#define _PROG_VOLUME_PCA

#include <data/xmipp_funcs.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>
#include <data/mask.h>
#include <data/basic_pca.h>

///@defgroup VolumePCA Volume PCA
///@ingroup ReconsLibrary
//@{
/** Volume PCA parameters. */
class ProgVolumePCA: public XmippProgram
{
public:
    /// Input set of volumes
    FileName fnVols;
    /// Output set of volumes
    FileName fnVolsOut;
    /// Number of PCA bases
    int NPCA;
    /// Mask
    Mask mask;
    /// Output basis
    FileName fnBasis;
    /// List of percentiles to generate volumes
    StringVector listOfPercentiles;
    /// Average volume
    FileName fnAvgVol;
    /// Output PCA stack
    FileName fnOutStack;
public:
    // Metadata with volumes
    MetaData mdVols;

    // Input volume
    Image<double> V;

    // PCA analyzer
    PCAMahalanobisAnalyzer analyzer;
public:
    /// Read arguments
    void readParams();

    /// Show
    void show() const;

    /// Define parameters
    void defineParams();

    /** Produce side info.*/
    void produce_side_info();

    /** Run */
    void run();
};
//@}
#endif
