/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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


//-----------------------------------------------------------------------------
// xmippVectorsOps.hh
//-----------------------------------------------------------------------------

#ifndef XMIPP_VECTOR_OPS_H
#define XMIPP_VECTOR_OPS_H


#include <math.h>                // for fabs()
#include <stdexcept>             // runtime_error
#include <vector>                // vector
#include <functional>            // for plus<T>, minus<T>
#include <algorithm>             // for transform()
#include <strstream>             // for ostrstream
#include "xmippUniform.hh"    // for Uniform (Random)

using namespace std;

//-----------------------------------------------------------------------------

/**@name Vector Operations*/
//@{

/**
 *    Adds a scalar "a" to a vector "v"
 *    returns a vector	
 *    example: r = v+a 
 */
template<class T> vector<T> operator+(const vector<T>& v, const T& a)
{
  vector<T> tmp = v;
  transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(plus<T>(), a));
  return tmp;
}

/**
 *    Adds a vector "v" to a scalar "a" 
 *    returns a vector	
 *    example: r = a+v 
 */
template<class T> vector<T> operator+(const T& a, const vector<T>& v)
{
  vector<T> tmp = v;
  transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(plus<T>(), a));
  return tmp;
}

/**
 *    Substract a scalar "a" to a vector "v"
 *    returns a vector	
 *    example: r = v-a 
 */
template<class T> vector<T> operator-(const vector<T>& v, const T& a)
{
  vector<T> tmp = v;
  transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(minus<T>(), a));
  return tmp;
}

/**
 *    Substract a vector "v" to a scalar "a" 
 *    returns a vector	
 *    example: r = a-v 
 */
template<class T> vector<T> operator-(const T& a, const vector<T>& v)
{
  vector<T> tmp(v.size(), a);
  transform(tmp.begin(), tmp.end(), v.begin(), tmp.begin(), minus<T>());
  return tmp;
}

/**
 *    Multiply a vector "v" by a scalar "a" 
 *    returns a vector	
 *    example: r = v*a 
 */
template<class T> vector<T> operator*(const vector<T>& v, T a)
{
  vector<T> tmp = v;
  transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(multiplies<T>(), a));
  return tmp;
}

/**
 *    Multiply a scalar "a" by a vector "v" 
 *    returns a vector	
 *    example: r = a*v 
 */
template<class T> vector<T> operator*( T a, const vector<T>& v)
{
  vector<T> tmp = v;
  transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(multiplies<T>(), a));
  return tmp;
}


/**
 *    Divides a vector "v" by a scalar "a" 
 *    returns a vector	
 *    example: r = v/a 
 */
template<class T> vector<T> operator/(const vector<T>& v, T a)
{
  vector<T> tmp = v;
  transform(tmp.begin(), tmp.end(), tmp.begin(), bind2nd(divides<T>(), a));
  return tmp;
}

/**
 *    Divides a scalar "a" by a vector "v" 
 *    returns a vector	
 *    example: r = a/v
 */
template<class T> vector<T> operator/( T a, const vector<T>& v)
{
  vector<T> tmp(v.size(), a);
  transform(tmp.begin(), tmp.end(), v.begin(), tmp.begin(), divides<T>());
  return tmp;
}

/**
 *    *= operator 
 *    returns a vector	
 *    example: v = v*r 
 */
template<class T> vector<T>& operator*=(vector<T>& v1, const vector<T>& v2)
{
  transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), multiplies<T>());
  return v1;
}

/**
 *    /= operator 
 *    returns a vector	
 *    example: v = v/r 
 */
template<class T> vector<T>& operator/=(vector<T>& v1, const vector<T>& v2)
{
  transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), divides<T>());
  return v1;
}


//-----------------------------------------------------------------------------

/**
 * Returns a random initialized vector
 * @param _size   Size of vector
 * @param _lower  Lower value for random elements
 * @param _upper  Upper value for random elements
 */
template <class T>
vector<T> randomVector(const unsigned& _size, const T& _lower,
     	               const T& _upper)
{
  vector<T> v(_size);
  xmippUniform<T> u(_lower, _upper);
  typename vector<T>::iterator i;
  for (i=v.begin() ; i<v.end() ; *i++ = u());

  return v;
};


/**
 * += operator. Ref Stroustrup, p 301.
 * @param _v1  First argument
 * @param _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
vector<T>& operator += (vector<T>& _v1, const vector<T>& _v2)
{
  if (_v1.size()!=_v2.size())
    throw runtime_error("different size vectors in +=");

  transform( _v1.begin(), _v1.end(), _v2.begin(), _v1.begin(), 
	     plus<T>() );

  return _v1;
};

/**
 * -= operator. Ref Stroustrup, p 301.
 * @param _v1  First argument
 * @param _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
vector<T>& operator -= (vector<T>& _v1, const vector<T>& _v2 )
{
  if (_v1.size()!=_v2.size())
    throw runtime_error("different size vectors in -=");

  transform( _v1.begin(), _v1.end(), _v2.begin(), _v1.begin(),
	     minus<T>() );

  return _v1;
};

/**
 * + operator. Ref Stroustrup, p 301.
 * @param _v1  First argument
 * @param _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T> 
vector<T> operator + ( const vector<T>& _v1, const vector<T>& _v2)
{
  if (_v1.size()!=_v2.size())
    throw runtime_error("different size vectors in +");

  // temp. container
  vector<T> tmp(_v1.size());

  transform( _v1.begin(), _v1.end(), _v2.begin(), tmp.begin(),
	     plus<T>() );
  return tmp;
};

/**
 * - operator. Ref Stroustrup, p 301.
 * @param _v1  First argument
 * @param _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
vector<T> operator - ( const vector<T>& _v1, const vector<T>& _v2 )
{
  if (_v1.size()!=_v2.size())
    throw runtime_error("different size vectors in -");

  // temp. container
  vector<T> tmp(_v1.size());

  transform( _v1.begin(), _v1.end(), _v2.begin(), tmp.begin(),
	     minus<T>() );
  return tmp;
};

/**
 * Dot product. Ref Stroustrup, p 301.
 * @param _v1  First argument
 * @param _v2  Second argument
 * @exception DifferentSize if _v1 and _v2  hasn't the same size
 */
template <class T>
T operator * ( const vector<T>& _v1, const vector<T>& _v2 )
{
  if (_v1.size()!=_v2.size())
    throw runtime_error("different size vectors in *");

  T dotProd = 0;
  typename vector<T>::const_iterator i, j;
  for ( i=_v1.begin(), j=_v2.begin() ; i!=_v1.end() ; i++, j++ )
    dotProd += *i * *j;

  return dotProd;
};

/**
 * Multiplies each element times _alpha. Ref Stroustrup, p 301.
 * @param _v      The container
 * @param _alpha  The element to be multiplied by each element in the container
 */
template <class T>
vector<T>& operator *= ( vector<T>& _v, const T _alpha )
{
  typename vector<T>::iterator i;
  for ( i=_v.begin() ; i!=_v.end() ; i++ )
    *i *= _alpha;

  return _v;
};


/**
 * Divides each element by _alpha. Ref Stroustrup, p 301.
 * @param _v      The container
 * @param _alpha  The element to be multiplied by each element in the container
 */
template <class T>
vector<T>& operator /= ( vector<T>& _v, const T _alpha )
{
  typename vector<T>::iterator i;
  for ( i=_v.begin() ; i!=_v.end() ; i++ )
    *i /= _alpha;

  return _v;
};

/**
 * Standard output for a vector<T>
 * @param _os The output stream
 * @param _v  The container to be printed
 */
template <class T>
ostream& operator << ( ostream& _os, const vector<T>& _v )
{
  _os << "< ";
  typename vector<T>::const_iterator i;
  for (i=_v.begin(); i!=_v.end(); i++) _os << *i << " ";
  //CO: copy( _v.begin(), _v.end(), ostream_iterator<T>( _os, " "));
  _os<< ">";  
  return _os;
};

/**
 * Standard input for a vector<T>
 * @param _is  The input stream
 * @param _v   The container to be read
 * @exception  runtime_error  If there are problems reading the vector
 */
template <class T> 
istream& operator >> ( istream& _is, vector<T>& _v )
{
  _v.clear();

  char c;
  _is >> c;
  
  if (_is && c!='<')
    _is.setstate(ios::failbit);

  bool finish = false;
  
  while (_is && !finish)
  {
    _is >> c;
    
    if (!_is)
      return _is;
    
    if (c=='>')
      finish=true;
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
    throw runtime_error("Error reading the vector");
  
  return _is;
};
//@}

//-----------------------------------------------------------------------------

#endif

