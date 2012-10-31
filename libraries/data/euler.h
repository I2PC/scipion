/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
/* Follows graphics gems iv and IMATHEULER ideas */
#ifndef _EULER_HH
#define _EULER_HH

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matrix1d.h"
#include "matrix2d.h"
#include "metadata.h"
//@{
/** Routines for computing Euler angles following different
 * approaches
*/

class Euler
{
public:

    typedef enum
    {
        //
        //  All 24 possible orderings
        //

        XYZ = 0x0101, // "usual" orderings
        XZY = 0x0001,
        YZX = 0x1101,
        YXZ = 0x1001,
        ZXY = 0x2101,
        ZYX = 0x2001,

        XZX = 0x0011, // first axis repeated
        XYX = 0x0111,
        YXY = 0x1011,
        YZY = 0x1111,
        ZYZ = 0x2011,
        ZXZ = 0x2111,

        XYZr = 0x2000, // relative orderings -- not common
        XZYr = 0x2100,
        YZXr = 0x1000,
        YXZr = 0x1100,
        ZXYr = 0x0000,
        ZYXr = 0x0100,

        XZXr = 0x2110, // relative first axis repeated
        XYXr = 0x2010,
        YXYr = 0x1110,
        YZYr = 0x1010,
        ZYZr = 0x0110,
        ZXZr = 0x0010,
        //       ||||
        //       VVVV
        //Legend:ABCD
        //  A -> Initial Axis (0==x, 1==y, 2==z)
        //  B -> Parity Even (1==true)
        //  C -> Initial Repeated (1==true)
        //  D -> Frame Static (1==true)
        //

        Legal =   XYZ | XZY | YZX | YXZ | ZXY | ZYX |
                  XZX | XYX | YXY | YZY | ZYZ | ZXZ |
                  XYZr| XZYr| YZXr| YXZr| ZXYr| ZYXr|
                  XZXr| XYXr| YXYr| YZYr| ZYZr| ZXZr,

        //        Min = 0x0000,
        //        Max = 0x2111,
        eulerDefault = ZYZ
    } eulerOrder;

    enum Axis { axisX = 0, axisY = 1, axisZ = 2 };

    enum InputLayout { XYZLayout, IJKLayout };

    //----------------------------------------------------------------
    // Constructors -- all default to ZYZ non-relative
    //   (where there is no argument to specify it)
    //----------------------------------------------------------------
    Euler();
    void init(void);
    Euler(const Euler&);
    Euler(eulerOrder p);
    Euler(const Matrix1D<double> &v,
          eulerOrder o = eulerDefault,
          InputLayout l = IJKLayout);
    Euler(double i, double j,
          double k, eulerOrder o = eulerDefault, InputLayout l = IJKLayout);
    Euler(const Euler &euler, eulerOrder newp);
    Euler(const Matrix2D<double> &, eulerOrder o = eulerDefault);
    //Euler(const Matrix44<T> &, Order o = Default);

    //---------------------------------
    //  Algebraic functions/ Operators
    //---------------------------------

    const Euler& operator=  (const Euler&);
    const Euler& operator=  (const Matrix1D<double>&);

    //--------------------------------------------------------
    // Set the euler value
    //  This does NOT convert the angles, but setXYZVector()
    // does reorder the input vector.
    //--------------------------------------------------------

    static bool  legal(eulerOrder);

    void  setXYZVector(const Matrix1D<double> &);

    eulerOrder  order() const;
    void  setOrder(eulerOrder);

    void  set(Axis initial,
              bool relative,
              bool parityEven,
              bool firstRepeats);

    //---------------------------------------------------------
    // Conversions, toXYZVector() reorders the angles so that
    //  the X rotation comes first, followed by the Y and Z
    //  in cases like XYX ordering, the repeated angle will be
    // in the "z" component
    //---------------------------------------------------------

    void  extract(const Matrix2D<double>&);
    //void  extract(const Matrix44<T>&);
    //void  extract(const Quat<T>&);

    void toMatrix(Matrix2D<double>& M) const;
    //Matrix1D<double>  toMatrix33() const;
    //Matrix44<T>  toMatrix44() const;
    /////Quat<T>  toQuat() const;
    void toXYZVector(Matrix1D<double> v) const;

    //---------------------------------------------------
    // Use this function to unpack angles from ijk form
    //---------------------------------------------------

    void  angleOrder(int &i, int &j, int &k) const;

    //---------------------------------------------------
    // Use this function to determine mapping from xyz to ijk
    // - reshuffles the xyz to match the order
    //---------------------------------------------------

    void  angleMapping(int &i, int &j, int &k) const;

    //----------------------------------------------------------------------
    //
    //  Utility methods for getting continuous rotations. None of these
    //  methods change the orientation given by its inputs (or at least
    //  that is the intent).
    //
    //    angleMod() converts an angle to its equivalent in [-PI, PI]
    //
    //    simpleXYZRotation() adjusts xyzRot so that its components differ
    //                        from targetXyzRot by no more than +-PI
    //
    //    nearestRotation() adjusts xyzRot so that its components differ
    //                      from targetXyzRot by as little as possible.
    //                      Note that xyz here really means ijk, because
    //                      the order must be provided.
    //
    //    makeNear() adjusts "this" Euler so that its components differ
    //               from target by as little as possible. This method
    //               might not make sense for Eulers with different order
    //               and it probably doesn't work for repeated axis and
    //               relative orderings (TODO).
    //
    //-----------------------------------------------------------------------

    static double angleMod (double angle);
    static void  simpleXYZRotation (Matrix1D<double> &xyzRot,
                                    const Matrix1D<double> &targetXyzRot);
    static void  nearestRotation (Matrix1D<double> &xyzRot,
                                  const Matrix1D<double> &targetXyzRot,
                                  eulerOrder order = XYZ);

    void  makeNear (const Euler &target);

    bool  frameStatic() const
    {
        return _frameStatic;
    }
    bool  initialRepeated() const
    {
        return _initialRepeated;
    }
    bool  parityEven() const
    {
        return _parityEven;
    }
    Axis  initialAxis() const
    {
        return _initialAxis;
    }

    Matrix1D <double> vec3;
    double x,y,z;

    void eulerRotate (Matrix2D <double> &M, const Matrix1D <double> &r);

protected:

bool  _frameStatic  :
    1; // relative or static rotations
bool  _initialRepeated :
    1; // init axis repeated as last
bool  _parityEven  :
    1; // "parity of axis permutation"
Axis  _initialAxis  :
    2; // First axis of rotation
};
#define eulerOrderNumber 24
Euler::eulerOrder eulerOrderList[eulerOrderNumber] =
    {
        Euler::XYZ , Euler::XZY , Euler::YZX , Euler::YXZ , Euler::ZXY , Euler::ZYX ,
        Euler::XZX , Euler::XYX , Euler::YXY , Euler::YZY , Euler::ZYZ , Euler::ZXZ ,
        Euler::XYZr, Euler::XZYr, Euler::YZXr, Euler::YXZr, Euler::ZXYr, Euler::ZYXr,
        Euler::XZXr, Euler::XYXr, Euler::YXYr, Euler::YZYr, Euler::ZYZr, Euler::ZXZr
    };

    //@}
#endif
