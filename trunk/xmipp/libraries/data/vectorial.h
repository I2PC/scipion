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

#ifndef VECTORIAL_H
#define VECTORIAL_H

#include "multidim_array.h"

/// @defgroup Vectorial Vector volumes
/// @ingroup MultidimensionalArrays

/** @defgroup VectorialSpeedUp Speed up functions
 * @ingroup Vectorial
 */

/** For all elements in a volume fashion.
 * @ingroup VectorialSpeedUp
 *
 * An (k,i,j) position is available inside the loop.
 */
#define FOR_ALL_ELEMENTS_IN_VECTORIAL_MATRIX3D(v) \
    FOR_ALL_ELEMENTS_IN_ARRAY3D((v).__X)

/** For all elements in a multidim array fashion.
 * @ingroup VectorialSpeedUp
 *
 * A single index (i) is available.
 */
#define FOR_ALL_ELEMENTS_IN_MULTIDIM_VECTORIAL_MATRIX3D(v) \
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY((v).__X)

/** Vectorial volume.
 * @ingroup Vectorial
 *
 * A vectorial volume is a "normal" MultidimArray whose elements are vectors
 * instead of single elements are doubles, floats, ... You can access
 * independently to any of the three components as a whole volume, ie, a
 * volume with all the X components, another with all Y components, ... or
 * access to the vector at a given position (logical positions as in ),
 * the X component at that position or the X component at a given
 * multidimensional array position.
 *
 * You can perform arithmetic operations on these vectors, and other common
 * operations such as resize, printShape, ...
 */
class Vectorial_MultidimArray
{
    // The 3 components
    MultidimArray< double > __X;
    MultidimArray< double > __Y;
    MultidimArray< double > __Z;

public:
    /// @defgroup VectorialShape Shape
    /// @ingroup Vectorial

    /** Resize.
     * @ingroup VectorialShape
     */
    void resize(int Zdim, int Ydim, int Xdim)
    {
        __X.resize(Zdim, Ydim, Xdim);
        __Y.resize(Zdim, Ydim, Xdim);
        __Z.resize(Zdim, Ydim, Xdim);
    }

    /** Resize with pattern.
     * @ingroup VectorialShape
     */
    void resize(const MultidimArray< double >& V)
    {
        __X.resize(V);
        __Y.resize(V);
        __Z.resize(V);
    }

    /** Clear.
     * @ingroup VectorialShape
     */
    void clear()
    {
        __X.clear();
        __Y.clear();
        __Z.clear();
    }

    /** Print shape.
     * @ingroup VectorialShape
     */
    void printShape() const
    {
        __X.printShape();
    }

    /** Init zeros.
     * @ingroup VectorialShape
     */
    void initZeros()
    {
        __X.initZeros();
        __Y.initZeros();
        __Z.initZeros();
    }

    /** Set Xmipp origin.
     * @ingroup VectorialShape
     */
    void setXmippOrigin()
    {
        __X.setXmippOrigin();
        __Y.setXmippOrigin();
        __Z.setXmippOrigin();
    }

    /// @defgroup VectorialAccess Component access
    /// @ingroup Vectorial

    /** Vector at a position.
     * @ingroup VectorialAccess
     *
     * The returned vector is a column 3x1 vector. The integer position are
     * logical indexes inside the MultidimArray.
     */
    Matrix1D< double > vector_at(int k, int i, int j) const
    {
        Matrix1D< double > result(3);
        XX(result) = __X(k, i, j);
        YY(result) = __Y(k, i, j);
        ZZ(result) = __Z(k, i, j);
        return result;
    }

    /** Vector at a position, result as argument.
     * @ingroup VectorialAccess
     */
    void vector_at(int k, int i, int j, Matrix1D< double >& result) const
    {
        result.resize(3);
        XX(result) = __X(k, i, j);
        YY(result) = __Y(k, i, j);
        ZZ(result) = __Z(k, i, j);
    }

    /** Constant access to X component.
     * @ingroup VectorialAccess
     *
     * X components are a MultidimArray.
     */
    const MultidimArray< double >& X() const
    {
        return __X;
    }

    /** Access to X components.
     * @ingroup VectorialAccess
     */
    MultidimArray< double >& X()
    {
        return __X;
    }

    /** Get the X components.
     * @ingroup VectorialAccess
     */
    void get_X(MultidimArray< double >& _XXX)
    {
        _XXX = __X;
    }

    /** Set the X components.
     * @ingroup VectorialAccess
     */
    void set_X(const MultidimArray< double >& _XXX)
    {
        __X = _XXX;
    }

    /** Constant access to Y component.
     * @ingroup VectorialAccess
     */
    const MultidimArray< double >& Y() const
    {
        return __Y;
    }

    /** Access to Y components.
     * @ingroup VectorialAccess
     */
    MultidimArray< double >& Y()
    {
        return __Y;
    }

    /** Get the Y components.
     * @ingroup VectorialAccess
     */
    void get_Y(MultidimArray< double >& _Y)
    {
        _Y = __Y;
    }

    /** Set the Y components.
     * @ingroup VectorialAccess
     */
    void set_Y(const MultidimArray< double >& _Y)
    {
        __Y = _Y;
    }

    /** Constant access to Z component.
     * @ingroup VectorialAccess
     */
    const MultidimArray< double >& Z() const
    {
        return __Z;
    }

    /** Access to Z components.
     * @ingroup VectorialAccess
     */
    MultidimArray< double >& Z()
    {
        return __Z;
    }

    /** Get the Z components.
     * @ingroup VectorialAccess
     */
    void get_Z(MultidimArray< double >& _Z)
    {
        _Z = __Z;
    }

    /** Set the Z components.
     * @ingroup VectorialAccess
     */
    void set_Z(const MultidimArray< double >& _Z)
    {
        __Z = _Z;
    }

    /** Constant access to a particular X component.
     * @ingroup VectorialAccess
     */
    double X(int k, int i, int j) const
    {
        return A3D_ELEM(__X, k, i, j);
    }

    /** Access to a particular X component.
     * @ingroup VectorialAccess
     */
    double& X(int k, int i, int j)
    {
        return A3D_ELEM(__X, k, i, j);
    }

    /** Get the X component at (k,i,j).
     * @ingroup VectorialAccess
     */
    double get_X_component(int k, int i, int j)
    {
        return X(k, i, j);
    }

    /** Set the X component at (k,i,j).
     * @ingroup VectorialAccess
     */
    double set_X_component(int k, int i, int j, double val)
    {
        X(k, i, j) = val;
    }

    /** Constant access to a particular Y component.
     * @ingroup VectorialAccess
     */
    double Y(int k, int i, int j) const
    {
        return A3D_ELEM(__Y, k, i, j);
    }

    /** Access to a particular Y component.
     * @ingroup VectorialAccess
     */
    double& Y(int k, int i, int j)
    {
        return A3D_ELEM(__Y, k, i, j);
    }

    /** Get the Y component at (k,i,j).
     * @ingroup VectorialAccess
     */
    double get_Y_component(int k, int i, int j)
    {
        return Y(k, i, j);
    }

    /** Set the Y component at (k,i,j).
     * @ingroup VectorialAccess
     */
    double set_Y_component(int k, int i, int j, double val)
    {
        Y(k, i, j) = val;
    }

    /** Constant access to a particular Z component.
     * @ingroup VectorialAccess
     */
    double Z(int k, int i, int j) const
    {
        return A3D_ELEM(__Z, k, i, j);
    }

    /** Access to a particular Z component.
     * @ingroup VectorialAccess
     */
    double& Z(int k, int i, int j)
    {
        return A3D_ELEM(__Z, k, i, j);
    }

    /** Get the Z component at (k,i,j).
     * @ingroup VectorialAccess
     */
    double get_Z_component(int k, int i, int j)
    {
        return Z(k, i, j);
    }

    /** Set the Z component at (k,i,j).
     * @ingroup VectorialAccess
     */
    double set_Z_component(int k, int i, int j, double val)
    {
        Z(k, i, j) = val;
    }

    /// @defgroup VectorialUtilitites Utilities
    /// @ingroup Vectorial

    /** Substitute each vector by its unit vector.
     * @ingroup VectorialUtilitites
     */
    void normalize_all_vectors()
    {
        FOR_ALL_ELEMENTS_IN_MULTIDIM_VECTORIAL_MATRIX3D(*this)
        {
            double mod = sqrt(
	    	DIRECT_MULTIDIM_ELEM(__X,n) * DIRECT_MULTIDIM_ELEM(__X,n) + 
		DIRECT_MULTIDIM_ELEM(__Y,n) * DIRECT_MULTIDIM_ELEM(__Y,n) + 
		DIRECT_MULTIDIM_ELEM(__Z,n) * DIRECT_MULTIDIM_ELEM(__Z,n));
            DIRECT_MULTIDIM_ELEM(__X,n) /= mod;
            DIRECT_MULTIDIM_ELEM(__Y,n) /= mod;
            DIRECT_MULTIDIM_ELEM(__Z,n) /= mod;
        }
    }

    /** Module of all vectors.
     * @ingroup VectorialUtilitites
     *
     * A volume with all vector modules at each position is returned.
     */
    void module(MultidimArray< double >& result) const
    {
        result.resize(__X);

        FOR_ALL_ELEMENTS_IN_MULTIDIM_VECTORIAL_MATRIX3D(*this)
            DIRECT_MULTIDIM_ELEM(result,  n) = sqrt(
	    	DIRECT_MULTIDIM_ELEM(__X,n) * DIRECT_MULTIDIM_ELEM(__X,n) + 
		DIRECT_MULTIDIM_ELEM(__Y,n) * DIRECT_MULTIDIM_ELEM(__Y,n) + 
		DIRECT_MULTIDIM_ELEM(__Z,n) * DIRECT_MULTIDIM_ELEM(__Z,n));
    }

    /** Write.
     * @ingroup VectorialUtilitites
     *
     * Given a FileName, this function saves 3 volumes named "_X", "_Y" and
     * "_Z".
     */
    void write(const FileName& fn) const
    {
        VolumeXmipp V;

        V = __X;
        V.write(fn.insert_before_extension("_X"));
        V = __Y;
        V.write(fn.insert_before_extension("_Y"));
        V = __Z;
        V.write(fn.insert_before_extension("_Z"));
    }

    /// @defgroup VectorialArithmetic Arithmetic operations
    /// @ingroup Vectorial

#define OPERATION(func, arg1, arg2, result, op) \
    func((arg1).__X, (arg2).__X, (result).__X, op); \
    func((arg1).__Y, (arg2).__Y, (result).__Y, op); \
    func((arg1).__Z, (arg2).__Z, (result).__Z, op);

    /** v3=v1+v2.
     * @ingroup VectorialArithmetic
     */
    Vectorial_MultidimArray operator+(const Vectorial_MultidimArray& op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByArray, *this, op1, temp, '+');
        return temp;
    }

    /** v3=v1-v2.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator-(const Vectorial_MultidimArray& op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByArray, *this, op1, temp, '-');
        return temp;
    }

    /** v3=v1*v2.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator*(const Vectorial_MultidimArray& op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByArray, *this, op1, temp, '*');
        return temp;
    }

    /** v3=v1/v2.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator/(const Vectorial_MultidimArray& op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByArray, *this, op1, temp, '/');
        return temp;
    }

    /** v3=v1^v2.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator^(const Vectorial_MultidimArray& op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByArray, *this, op1, temp, '^');
        return temp;
    }

    /** v3+=v2.
     * @ingroup VectorialAritmetic
     */
    void operator+=(const Vectorial_MultidimArray& op1)
    {
        OPERATION(arrayByArray, *this, op1, *this, '+');
    }

    /** v3-=v2.
     * @ingroup VectorialAritmetic
     */
    void operator-=(const Vectorial_MultidimArray& op1)
    {
        OPERATION(arrayByArray, *this, op1, *this, '-');
    }

    /** v3*=v2.
     * @ingroup VectorialAritmetic
     */
    void operator*=(const Vectorial_MultidimArray& op1)
    {
        OPERATION(arrayByArray, *this, op1, *this, '*');
    }

    /** v3/=v2.
     * @ingroup VectorialAritmetic
     */
    void operator/= (const Vectorial_MultidimArray& op1)
    {
        OPERATION(arrayByArray, *this, op1, *this, '/');
    }

    /** v3^=v2.
     * @ingroup VectorialAritmetic
     */
    void operator^=(const Vectorial_MultidimArray& op1)
    {
        OPERATION(arrayByArray, *this, op1, *this, '^');
    }

#undef OPERATION
#define OPERATION(func, arg1, arg2, result, op) \
    func((arg1).__X, (arg2), (result).__X, op); \
    func((arg1).__Y, (arg2), (result).__Y, op); \
    func((arg1).__Z, (arg2), (result).__Z, op);

    /** v3=v1+k.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator+(double op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByScalar, *this, op1, temp, '+');
        return temp;
    }

    /** v3=v1-k.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator-(double op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByScalar, *this, op1, temp, '-');
        return temp;
    }

    /** v3=v1*k.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator*(double op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByScalar, *this, op1, temp, '*');
        return temp;
    }

    /** v3=v1/k.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator/(double op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByScalar, *this, op1, temp, '/');
        return temp;
    }

    /** v3=v1^k.
     * @ingroup VectorialAritmetic
     */
    Vectorial_MultidimArray operator^(double op1) const
    {
        Vectorial_MultidimArray temp;
        OPERATION(arrayByScalar, *this, op1, temp, '^');
        return temp;
    }

    /** v3+=k.
     * @ingroup VectorialAritmetic
     */
    void operator+=(const double& op1)
    {
        OPERATION(arrayByScalar, *this, op1, *this, '+');
    }

    /** v3-=k.
     * @ingroup VectorialAritmetic
     */
    void operator-=(const double& op1)
    {
        OPERATION(arrayByScalar, *this, op1, *this, '-');
    }

    /** v3*=k.
     * @ingroup VectorialAritmetic
     */
    void operator*=(const double& op1)
    {
        OPERATION(arrayByScalar, *this, op1, *this, '*');
    }

    /** v3/=k.
     * @ingroup VectorialAritmetic
     */
    void operator/=(const double& op1)
    {
        OPERATION(arrayByScalar, *this, op1, *this, '/');
    }

    /** v3^=k.
     * @ingroup VectorialAritmetic
     */
    void operator^=(const double& op1)
    {
        OPERATION(arrayByScalar, *this, op1, *this, '^');
    }

#undef Vectorial_MultidimArray

};


#endif
