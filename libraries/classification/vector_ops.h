/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
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


//-----------------------------------------------------------------------------
// xmippVectorsOps.hh
//-----------------------------------------------------------------------------

#ifndef XMIPP_VECTOR_OPS_H
#define XMIPP_VECTOR_OPS_H

#include <cmath>
#include <stdexcept>
#include <vector>
#include <functional>
#include <algorithm>
#include <sstream>

#include "uniform.h"

/**@defgroup VectorOperations Vector Operations
   @ingroup ClassificationLibrary */
//@{

/**
 *    Adds a scalar "a" to a vector "v"
 *    returns a vector
 *    example: r = v+a
 */
template<class T> std::vector<T> operator+(const std::vector<T>& v, const T& a)
{
    std::vector<T> tmp = v;
    transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(std::plus<T>(), a));
    return tmp;
}

/**
 *    Adds a vector "v" to a scalar "a"
 *    returns a vector
 *    example: r = a+v
 */
template<class T> std::vector<T> operator+(const T& a, const std::vector<T>& v)
{
    std::vector<T> tmp = v;
    transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(std::plus<T>(), a));
    return tmp;
}

/**
 *    Substract a scalar "a" to a vector "v"
 *    returns a vector
 *    example: r = v-a
 */
template<class T> std::vector<T> operator-(const std::vector<T>& v, const T& a)
{
    std::vector<T> tmp = v;
    transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(std::minus<T>(), a));
    return tmp;
}

/**
 *    Substract a vector "v" to a scalar "a"
 *    returns a vector
 *    example: r = a-v
 */
template<class T> std::vector<T> operator-(const T& a, const std::vector<T>& v)
{
    std::vector<T> tmp(v.size(), a);
    transform(tmp.begin(), tmp.end(), v.begin(), tmp.begin(), std::minus<T>());
    return tmp;
}

/**
 *    Multiply a vector "v" by a scalar "a"
 *    returns a vector
 *    example: r = v*a
 */
template<class T> std::vector<T> operator*(const std::vector<T>& v, T a)
{
    std::vector<T> tmp = v;
    transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(std::multiplies<T>(), a));
    return tmp;
}

/**
 *    Multiply a scalar "a" by a vector "v"
 *    returns a vector
 *    example: r = a*v
 */
template<class T> std::vector<T> operator*(T a, const std::vector<T>& v)
{
    std::vector<T> tmp = v;
    transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(std::multiplies<T>(), a));
    return tmp;
}


/**
 *    Divides a vector "v" by a scalar "a"
 *    returns a vector
 *    example: r = v/a
 */
template<class T> std::vector<T> operator/(const std::vector<T>& v, T a)
{
    std::vector<T> tmp = v;
    transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(std::divides<T>(), a));
    return tmp;
}

/**
 *    Divides a scalar "a" by a vector "v"
 *    returns a vector
 *    example: r = a/v
 */
template<class T> std::vector<T> operator/(T a, const std::vector<T>& v)
{
    std::vector<T> tmp(v.size(), a);
    transform(tmp.begin(), tmp.end(), v.begin(), tmp.begin(), std::divides<T>());
    return tmp;
}

/**
 *    *= operator
 *    returns a vector
 *    example: v = v*r
 */
template<class T> std::vector<T>& operator*=(std::vector<T>& v1, const std::vector<T>& v2)
{
    transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::multiplies<T>());
    return v1;
}

/**
 *    /= operator
 *    returns a vector
 *    example: v = v/r
 */
template<class T> std::vector<T>& operator/=(std::vector<T>& v1, const std::vector<T>& v2)
{
    transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::divides<T>());
    return v1;
}


//-----------------------------------------------------------------------------

/**
 * Returns a random initialized vector
 * Parameter: _size   Size of vector
 * Parameter: _lower  Lower value for random elements
 * Parameter: _upper  Upper value for random elements
 */
template <class T>
std::vector<T> randomVector(const unsigned& _size, const T& _lower,
                       const T& _upper)
{
    std::vector<T> v(_size);
    RandomUniformGenerator<T> u(_lower, _upper);
    typename std::vector<T>::iterator i;
    for (i = v.begin() ; i < v.end() ; *i++ = u());

    return v;
}

/**
 * += operator. Ref Stroustrup, p 301.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
std::vector<T>& operator += (std::vector<T>& _v1, const std::vector<T>& _v2)
{
    if (_v1.size() != _v2.size())
        throw std::runtime_error("different size vectors in +=");

    transform(_v1.begin(), _v1.end(), _v2.begin(), _v1.begin(),
              std::plus<T>());

    return _v1;
}

/**
 * -= operator. Ref Stroustrup, p 301.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
std::vector<T>& operator -= (std::vector<T>& _v1, const std::vector<T>& _v2)
{
    if (_v1.size() != _v2.size())
        throw std::runtime_error("different size vectors in -=");

    transform(_v1.begin(), _v1.end(), _v2.begin(), _v1.begin(),
              std::minus<T>());

    return _v1;
}

/**
 * + operator. Ref Stroustrup, p 301.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
std::vector<T> operator + (const std::vector<T>& _v1, const std::vector<T>& _v2)
{
    if (_v1.size() != _v2.size())
        throw std::runtime_error("different size vectors in +");

    // temp. container
    std::vector<T> tmp(_v1.size());

    transform(_v1.begin(), _v1.end(), _v2.begin(), tmp.begin(),
              std::plus<T>());
    return tmp;
}

/**
 * - operator. Ref Stroustrup, p 301.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
std::vector<T> operator - (const std::vector<T>& _v1, const std::vector<T>& _v2)
{
    if (_v1.size() != _v2.size())
        throw std::runtime_error("different size vectors in -");

    // temp. container
    std::vector<T> tmp(_v1.size());

    transform(_v1.begin(), _v1.end(), _v2.begin(), tmp.begin(),
              std::minus<T>());
    return tmp;
}

/**
 * Dot product. Ref Stroustrup, p 301.
 * Parameter: _v1  First argument
 * Parameter: _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
T operator *(const std::vector<T>& _v1, const std::vector<T>& _v2)
{
    if (_v1.size() != _v2.size())
        throw std::runtime_error("different size vectors in *");

    T dotProd = 0;
    typename std::vector<T>::const_iterator i, j;
    for (i = _v1.begin(), j = _v2.begin() ; i != _v1.end() ; i++, j++)
        dotProd += *i * *j;

    return dotProd;
}

/**
 * Multiplies each element times _alpha. Ref Stroustrup, p 301.
 * Parameter: _v      The container
 * Parameter: _alpha  The element to be multiplied by each element in the container
 */
template <class T>
std::vector<T>& operator *= (std::vector<T>& _v, const T _alpha)
{
    typename std::vector<T>::iterator i;
    for (i = _v.begin() ; i != _v.end() ; i++)
        *i *= _alpha;

    return _v;
}

/**
 * Divides each element by _alpha. Ref Stroustrup, p 301.
 * Parameter: _v      The container
 * Parameter: _alpha  The element to be multiplied by each element in the container
 */
template <class T>
std::vector<T>& operator /= (std::vector<T>& _v, const T _alpha)
{
    typename std::vector<T>::iterator i;
    for (i = _v.begin() ; i != _v.end() ; i++)
        *i /= _alpha;

    return _v;
}

/**
 * Standard output for a std::vector<T>
 * Parameter: _os The output stream
 * Parameter: _v  The container to be printed
 */
template <class T>
std::ostream& operator << (std::ostream& _os, const std::vector<T>& _v)
{
    _os << "< ";
    typename std::vector<T>::const_iterator i;
    for (i = _v.begin(); i != _v.end(); i++) _os << *i << " ";
    _os << ">";
    return _os;
}

/**
 * Standard input for a std::vector<T>
 * Parameter: _is  The input stream
 * Parameter: _v   The container to be read
 * @exception  runtime_error  If there are problems reading the vector
 */
template <class T>
std::istream& operator >> (std::istream& _is, std::vector<T>& _v)
{
    _v.clear();

    char c;
    _is >> c;

    if (_is && c != '<')
        _is.setstate(std::ios::failbit);

    bool finish = false;

    while (_is && !finish)
    {
        _is >> c;

        if (!_is)
            return _is;

        if (c == '>')
            finish = true;
        else
        {
            _is.putback(c);
            T item;
            _is >> item;

            if (_is)
                _v.push_back(item);
            else
                return _is;
        }
    }

    if (!_is)
        throw std::runtime_error("Error reading the vector");

    return _is;
}

/** Euclidean distance */
template <class T>
double euclideanDistance(const std::vector<T>& _v1, const std::vector<T>& _v2)
{
    if (_v1.size() != _v2.size())
        throw std::runtime_error("vector of different size in eDist");

    double dist = 0;
    typename std::vector<T>::const_iterator i, j;
    for (i = _v1.begin(), j = _v2.begin() ; i < _v1.end(); i++, j++)
    {
        double tmp = (double)(*i) - (double)(*j);
        dist += tmp * tmp;
    }

    return sqrt(dist);
}

/**
* Norm: norm of a vector (Euclidean distance to origin)
*
*/
template <class T>
T norm(const std::vector<T>& v)
{
    double sum = 0.0;
    typename std::vector<T>::const_iterator i;
    for (i = v.begin(); i != v.end(); i++)
        sum += (double)(*i) * (double)(*i);
    return (T) sqrt(sum);
}
//@}
#endif

