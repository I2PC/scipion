/***************************************************************************
 *
 * Authors:    Jose Luis Vilas          (jlvilas@cnb.csic.es)
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

#ifndef VALIDATE_TILT_PAIRS_H
#define VALIDATE_TILT_PAIRS_H

#include <data/xmipp_fftw.h>
#include <data/args.h>
#include <data/xmipp_funcs.h>

#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/xmipp_image.h>
#include <data/geometry.h>
#include <data/filters.h>
#include <data/xmipp_program.h>
#include <complex>

/**@defgroup Centilt align_tilt_pairs (Align tilted and untilted images in a random conical tilt experiment)
   @ingroup ReconsLibrary */
//@{
class ProgValidationTiltPairs: public XmippProgram
{
public:
    /**  Filename tiltpar average */
    FileName fntiltimage_In;
    /**  Filename tiltpar average */
    FileName fnuntiltimage_In;
    /**  Filename output document file */
    FileName fnOut;

public:
    /// Define parameters in the command line
    void defineParams();

    /// Transform a set of angles/orientations to a matrix expressed in Pauli basis
    void quaternion2Paulibasis(double rot, double tilt, double psi, std::complex<double>(&L)[4]);

    /// It calculates the product of two matrix expressed in Pauli matrices by their matrix elements A an B
    void Pauliproduct(std::complex<double> A[4], std::complex<double> B[4], std::complex<double> (&P)[4]);

    /// Inverse of a matrix (expressed in Pauli basis)
    void inverse_matrixSU2(std::complex<double> Original[4], std::complex<double> (&Inver)[4]);

    ///A Pauli span is converted to a 2x2 matrix
    void Paulibasis2matrix(std::complex<double> P[4], std::complex<double> (&M)[4]);

    //From a Pauli vecto, its inverse is calculated
    void InversefromPaulibasis(std::complex<double> Original[4],	std::complex<double> (&Inver)[4]);

    /// It takes a 2x2 matrix (where M[0] = m11; M[1]=m12; M[2]=m21; M[3]=m22) and express the matrix M,
    //in the Pauli Basis M=P[0] Sigma0 + P[1] Sigma1 + P[2] Sigma2 + P[3] Sigma3
    void matrix2Paulibasis(std::complex<double> M[4], std::complex<double> (&P)[4]);

    ///Extract angles alpha_x and alpha_y from the transformation matrix
    void extrarotationangles(std::complex<double> R[4], double &alpha_x, double &alpha_y);

    ///It takes two sets of angles (rotu, tiltu, psiu) which define Eu, and (rott, tiltt, psit) which define Et,
    // and returns the tranformation angles, alpha_x an alpha_y for the transformation Et=R(alpha_y)R(alpha_x)Eu.
    //Where Et an Eu, are the tilted and untilted images.
    void angles2tranformation(double untilt_angles[3], double tilt_angles[3], double alpha_x, double alpha_y);

    /// Read parameters from the command line
    void readParams();

    /// Execute de program
    void run();


};
//@}
#endif
