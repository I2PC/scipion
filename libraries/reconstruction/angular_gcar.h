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
#ifndef _PROG_ANGULAR_GCAR
#define _PROG_ANGULAR_GCAR

#include <data/xmipp_program.h>
#include <data/matrix2d.h>
#include <data/sparse_matrix2d.h>
#include <data/matrix1d.h>
#include <data/xmipp_funcs.h>
#include <math.h>
#include <vector>
#include <external/arpack++-2.3/include/blas1c.h>
#include <external/arpack++-2.3/include/lapackc.h>
#include <external/arpack++-2.3/include/arsnsym.h>


class IJpair{
public :
	int i;
	int j;

};

/**@defgroup AngularGCAR Globally Convergent Angular Reconstitution
   @ingroup ReconsLibrary */
//@{
/** Angular GCAR parameters. */
class ProgAngularGCAR: public XmippProgram
{
public:
public:
public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /** Run */
    void run();

    ///  qrand
    void qrand(Matrix2D<double>& q, int K);

	//! Convert a quaternion into a rotation matrix.
    /*!
      \param q is the input quaternion.
      \param rotMatrix is the rotation matrix of the quaternion.
     */
	void qToRot(const Matrix1D<double>& q, Matrix2D<double> &rotMatrix);

	//! Register the estimated orienations to the reference orientations.
	/*!
	  The registration is required to eliminate the arbitrary rotation and
	  reflection in the recovered orientations.
	 */
	/*!
	  \param PHI is the estimated orientation. It is used as output parameter.
	  \param PHIRef is the true orientations.
	 */
	void registerOrientations(Matrix2D<double>& PHI,const Matrix2D<double>& PHIRef);

	//! Compare the computed orientations to reference orientations.
	/*!
	  Each row in PHI and PHIref corresponds to a point on S2, representing the
	  direction vector of some Fourier ray.
	 */
	/*!
	  \param PHI is the estimated orientation.
	  \param PHIRef is the true orientation.
	 */
	void checkOrientations(const Matrix2D<double>& PHI,const Matrix2D<double>& PHIref);

	//! Given a Kx4 array of quaternions, find the common lines between projections k1 and k2.
	/*!
	  \param q is the input matrix.
	  \param k1 is the first projection.
	  \param k2 is the second projection.
	  \param nTheta is the Fourier rays in each projection.
	  \param idx1 is the index of the common line in projection k1.
	  \param idx2 is the index of the common line in projection k2.
	 */
	void commonLineQ(const Matrix2D<double>& q, int k1, int k2, int nTheta, int& idx1, int& idx2);

	void clmatrixCheatQ(const Matrix2D<double>& q, int nTheta, Matrix2D<int>& clmatrix, Matrix2D<int>&  clcorr);

	/*!
	  Use the geometry of 3D Fourier space to improve orientation estimation:
	  	  1) Find the plane of each circle using PCA
	   	  2) Find best equi-spacing between same circle points
	 */
	/*!
	  \param PHI is the input matrix.
	  \param PHI2 is the improved estimation of the orientations PHI.
	  \param K is the number of projections.
	  \param L is the angular resolution of each projection.
	 */
	void cryoCosmetify(const Matrix2D<double>& PHI, Matrix2D<double>& PHI2, int K, int L);

	//! Compute direction vectors of Fourier rays correspoding to given quaternions.
	/*!
	  For each quaternion, which corresponds to the projection
	  orientation of one of the projections, the function computes the
	  direction vectors of the Fourier rays in that projection.
	 */
	/*!
	  \param Q is a 4xK array containing quaternions.
	  \param nTheta is the Fourier rays in each projection.
	  \param PHIRef is the output parameter.
	 */
	void Q2S2(const Matrix2D<double>& Q, Matrix2D<double>& PHIRef, int NTheta);

	//! Recover the orientation of each ray from the matrix W.
	/*!
	  \param W is the spider matrix
	  \param PHI is the coordinates on S2 of each of the KL Fourier rays (KL by 3matrix).
	  \param nEigs is the number of eigenvalues and eigenvectors calculated
	 */
	void cryoOrientations(SparseMatrix2D &W, int nEigs, Matrix2D<double> &PHI);

	/*!
	 * \fn void cryo_S2graph(Matrix2D<int> &clmatrix,int L,int d, SparseMatrix2D &w ,SparseMatrix2D  &D)
	 * \brief  Construct the spiders graph that corresponds to the given common lines matrix
	 * \param  clmatrix   Common lines matrix. If ray l1 in projection k1 and ray l2 in projection k2 are a common line, then entry (k1,k2) is l1 and entry (k2,k1) is l2.
	 * \param  L Number of Fourier rays in each projection.
	 * \param  d Number of samples on each leg of the spider. The total number of samples on each arc of the spider is therefore 2d+1.
	 * \param  W Spider graph. Sparse matrix of size KLxKL. W is normalized to be row-stochastic.
	 */
	void cryo_S2graph(Matrix2D<int> &clmatrix,int L,int d, SparseMatrix2D &w);

private:
	/*!
	 * \fn int double_mod(int a, int b)
	 * \brief  Modulus operation implemented as matlab does(c++ has a different treatment for negative numbers)
	 * \param  a dividend
	 * \param  b divisor
	 */
	int matlab_mod(int a, int b);
	/*!
	 * \fn void sumRows(SparseMatrix2D &w,Matrix1D<double> &resul)
	 * \brief  Sums row by row the colums of a matrix
	 * \param  w The matrix whose colums will be summed
	 * \param  resul colum vector with the results of the operation
	 */
    void sumRows(SparseMatrix2D &w,Matrix1D<double> &resul);
    /*!
     * \fn void sumRows(SparseMatrix2D &w,Matrix1D<double> &resul)
     * \brief Counts the number of times that a pair appears
     * \param elem Vector of sparse elements that stores the number of times that a pair (i,j) appears
     * \param IJ Vector of pairs (i,j)
     * \param indx Size of vector IJ
     */
    void countElems(std::vector<SparseElement> &elem,std::vector<IJpair> &IJ,int indx);
};

class EigElement {
public:
	double eigenvalue;
	int    pos;
};

inline bool operator<(const EigElement& e1, const EigElement& e2)
{
	return e1.eigenvalue>e2.eigenvalue; //Descending order
}



inline bool operator<(const IJpair& _x, const IJpair& _y)
{
	return ( _x.i < _y.i) || (( _x.i== _y.i) && ( _x.j< _y.j));
}

inline bool operator==(const IJpair& _x, const IJpair& _y)
{
	return (_x.j == _y.j) && (_x.i == _y.i);
}


//@}
#endif
