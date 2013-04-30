/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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

#ifndef _TRANSFROM_DIMRED
#define _TRANSFROM_DIMRED

#include <data/metadata.h>
#include <data/xmipp_program.h>
#include <data/mask.h>
#include <data/filters.h>

/**@defgroup ProgTransformDimRed Obtain image distances between images in metada
   @ingroup DimRedLibrary */
//@{
/** Analyze cluster parameters. */
class ProgTransformDimRed: public XmippProgram
{
public:
    /** Filename selection file containing the images */
    FileName fnIn;
    FileName fnOut;
    String dimRefMethod;
    int outputDim;
public:
    MetaData SFin;
    Matrix2D<double> X; // Input data
    Mask mask;

    // Auxiliary variables for alignment
    MultidimArray<double> I1, I2;
    Matrix2D<double> M;
    AlignmentAux aux;
    CorrelationAux aux2;
	RotationalCorrelationAux aux3;
public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /// Produce side info
    void produceSideInfo();

    /// Main routine
    void run();

    /// From image to data matrix
    void insertImageInDataMatrix(size_t index, const MultidimArray<double> &mImg);

    /// From image to data matrix
    void extractImageFromDataMatrix(size_t index, MultidimArray<double> &mImg);

    /// Correlation distance between two images
    double progCorrelationDistance(size_t i1, size_t i2);

};
//@}
#endif
