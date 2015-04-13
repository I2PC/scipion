/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano coss@cnb.csic.es
 *              Slavica Jonic        Slavica.Jonic@impmc.jussieu.fr
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

/**@defgroup ART_pseudo reconstruct_art_pseudo (Reconstruct with pseudo atoms)
   @ingroup ReconsLibrary */
//@{

#include "data/xmipp_program.h"
#include "data/xmipp_filename.h"

/** Parameters for reconstructing with pseudoatoms. */
class ProgARTPseudo: public XmippProgram
{
public:
    /// Selfile with the input images
    FileName fnDoc;

    /// Pseudoatom filename
    FileName fnPseudo;

    /// Selfile with the input NMAs
    FileName fnNMA;

    /// Output filename
    FileName fnRoot;

    /// Lambda
    double lambdaART;

    /// Number of iterations
    int Nit;

    /// Sigma of atoms
    double sigma;

    /// Sampling rate
    double sampling;
public:
    /** Define parameters */
    void defineParams();

    /** Read parameters from command line. */
    void readParams();

    /** Show parameters */
    void show() const;

    /** Produce side info */
    void produceSideInfo();

    /** Run */
    void run();

    /** Write Pseudo */
    void writePseudo();

    /** ART single step */
    double ART_single_step(const MultidimArray<double> &Iexp,
        double rot, double tilt, double psi, double shiftX, double shiftY,
        const std::vector<double> &lambda);
public:
    // Input images
    MetaData DF;

    // Atomic positions
    std::vector< Matrix1D<double> > atomPosition;

    // Atomic weights
    std::vector< double > atomWeight;

    // Gaussian projection table
    Matrix1D<double> gaussianProjectionTable;

    // Gaussian projection2 table
    Matrix1D<double> gaussianProjectionTable2;

    // NMA modes
    std::vector < Matrix2D<double> > NMA;
};
//@}
