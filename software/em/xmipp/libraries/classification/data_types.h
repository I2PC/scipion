/***************************************************************************
 *
 * Authors:     Alberto Pascual
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
// xmippCDataTypes.h
//-----------------------------------------------------------------------------

#ifndef XMIPPCDATATYPES_H
#define XMIPPCDATATYPES_H

//-----------------------------------------------------------------------------

#include <vector>           // vector
#include <set>              // set
#include <string>           // string

// MAXDOUBLE
#include <limits.h>
#ifndef __CYGWIN__
#ifndef __APPLE__
#include <values.h>
#endif
#endif

#ifndef MAXFLOAT
#define  MAXFLOAT (float)1e38
#endif//MAXFLOAT

#ifndef MAXDOUBLE
#define  MAXDOUBLE  (double)1.79769313486231570e+308
#endif//MAXDOUBLE

#ifndef MAXINT
#define  MAXINT 2147483647
#endif//MAXFLOAT

#ifndef MAXZERO
#define  MAXZERO (float)1e-20   // Useful for defining "zero precision"
#endif//MAXZERO

//-----------------------------------------------------------------------------
/**@defgroup DataTypes Some type definitions
   @ingroup ClassificationLibrary */
//@{
/**
* Feature as float
*/
//typedef double Feature;
//ROB in some compiles Feature classes with feature class
//SO I rename it as floatFeature
typedef float floatFeature;

//-----------------------------------------------------------------------------

/**
* Label as string
*/
typedef std::string Label;

//-----------------------------------------------------------------------------

/**
* FeatureVector as vector of doubles
*/
typedef std::vector<floatFeature> FeatureVector;

//-----------------------------------------------------------------------------

/**
* FeatureSet as a set of xmippfeatures
*/
typedef std::set<floatFeature, std::less<floatFeature> > FeatureSet;
//@}

//-----------------------------------------------------------------------------

#endif//XMIPPCDATATYPES_H
