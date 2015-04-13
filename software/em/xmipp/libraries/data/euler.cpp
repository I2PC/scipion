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
#include "euler.h"

void
Euler::angleOrder(int &i, int &j, int &k) const
{
    i = _initialAxis;
    j = _parityEven ? (i+1)%3 : (i > 0 ? i-1 : 2);
    k = _parityEven ? (i > 0 ? i-1 : 2) : (i+1)%3;
}

void
Euler::angleMapping(int &i, int &j, int &k) const
{
    int m[3];

    m[_initialAxis] = 0;
    m[(_initialAxis+1) % 3] = _parityEven ? 1 : 2;
    m[(_initialAxis+2) % 3] = _parityEven ? 2 : 1;
    i = m[0];
    j = m[1];
    k = m[2];
}

void
Euler::setXYZVector(const Matrix1D<double> &v)
{
    int i,j,k;
    angleMapping(i,j,k);
    vec3(i) = XX(v);
    vec3(j) = YY(v);
    vec3(k) = ZZ(v);
}

//inline Matrix1D<double>
void Euler::toXYZVector(Matrix1D<double> v) const
{
    int i,j,k;
    angleMapping(i,j,k);
    v=vectorR3(vec3(i),vec3(j),vec3(k));
    //    return Vec3<T>((*this)[i],(*this)[j],(*this)[k]);
}

void Euler::init(void)
{
    x=0.;
    y=0.;
    z=0.;
    vec3.initZeros(3);
}

Euler::Euler() :
        _frameStatic(true),
        _initialRepeated(false),
        _parityEven(true),
        _initialAxis(axisX)
{
    init();
}

Euler::Euler(Euler::eulerOrder p) :
        _frameStatic(true),
        _initialRepeated(false),
        _parityEven(true),
        _initialAxis(axisX)
{
    init();
    setOrder(p);
}

Euler::Euler( const Matrix1D <double> &v,
                     Euler::eulerOrder p,
                     Euler::InputLayout l )
{
    init();
    setOrder(p);
    if ( l == XYZLayout )
        setXYZVector(v);
    else
    {
        x = XX(v);
        y = YY(v);
        z = ZZ(v);
    }
}

Euler::Euler(const Euler &euler)
{
    init();
    operator=(euler);
}


Euler::Euler(const Euler &euler,eulerOrder p)
{
    init();
    setOrder(p);
    Matrix2D<double> M(3,3);
    euler.toMatrix(M);
    extract(M);
}


Euler::Euler( double xi, double yi, double zi,
                     Euler::eulerOrder p,
                     Euler::InputLayout l)
{
    setOrder(p);
    if ( l == XYZLayout )
    {
        setXYZVector(vectorR3(xi,yi,zi));
    }

    else
    {
        x = xi;
        y = yi;
        z = zi;
    }
}

Euler::Euler( const Matrix2D<double> &M, Euler::eulerOrder p )
{
    setOrder(p);
    extract(M);
}

void Euler::extract(const Matrix2D<double> &M)
{
    int i,j,k;
    angleOrder(i,j,k);
    if (_initialRepeated)
    {
        // Extract the first angle, x.
        //

        //x = Math::atan2(M[j][i], M[k][i]);
        x = atan2(dMij(M,j,i), dMij(M,k,i));

        //
        // Remove the x rotation from M, so that the remaining
        // rotation, N, is only around two axes, and gimbal lock
        // cannot occur.
        //

        //Vec3<T> r (0, 0, 0);
        Matrix1D <double> r;
        r.initZeros(3);
        VEC_ELEM(r,i) = (_parityEven? -x: x);

        Matrix2D<double> N(4,4);
        N.initIdentity();
        eulerRotate(N,r);
        N = N * M;

        //
        // Extract the other two angles, y and z, from N.
        //

        double sy = sqrt (dMij(N,j,i)*dMij(N,j,i) +
                          dMij(N,k,i)*dMij(N,k,i) );
        y = atan2 (sy, dMij(N,i,i));
        z = atan2 (dMij(N,j,k),dMij(N,j,j));
    }
    else
    {
        //
        // Extract the first angle, x.
        //
        //        x = Math<T>::atan2 (M[j][k], M[k][k]);
        x = atan2 (dMij(M,j,k),dMij(M,k,k));
        //
        // Remove the x rotation from M, so that the remaining
        // rotation, N, is only around two axes, and gimbal lock
        // cannot occur.
        //

        Matrix1D <double> r;
        r.initZeros(3);
        VEC_ELEM(r,i) = (_parityEven? -x: x);
        //std::cerr << "DEBUG_ROB, r:" << r << std::endl;
        Matrix2D<double> N(4,4);
        N.initIdentity();
        eulerRotate(N,r);
        //std::cerr << "DEBUG_ROB, N:" << N << std::endl;
        N = N * M;
        //std::cerr << "DEBUG_ROB, NN:" << N << std::endl;

        //
        // Extract the other two angles, y and z, from N.
        //
        //        T cy = Math<T>::sqrt (N[i][i]*N[i][i] + N[i][j]*N[i][j]);
        double cy = sqrt (dMij(N,i,i)*dMij(N,i,i) +
                          dMij(N,i,j)*dMij(N,i,j) );
        //        y = Math<T>::atan2 (-N[i][k], cy);
        //        z = Math<T>::atan2 (-N[j][i], N[j][j]);
        y = atan2 (-dMij(N,i,k),cy);
        z = atan2 (-dMij(N,j,i),dMij(N,j,j));

    }

    if (!_parityEven)
        //*this *= -1;
    {
        vec3 *= -1;
        x *=-1;
        y *=-1;
        z *=-1;
    }

    if (!_frameStatic)
    {
        double t = x;
        x = z;
        z = t;
    }

}

void Euler::toMatrix(Matrix2D<double>& M) const
{
    int i,j,k;
    angleOrder(i,j,k);

    Matrix1D<double> angles;

    if ( _frameStatic )
        angles = vectorR3(x,y,z);
    else
        angles=vectorR3(z,y,x);

    if ( !_parityEven )
    {
        angles *= -1.0;
    }

    double ci = cos(XX(angles));
    double cj = cos(YY(angles));
    double ch = cos(ZZ(angles));
    double si = sin(XX(angles));
    double sj = sin(YY(angles));
    double sh = sin(ZZ(angles));

    double cc = ci*ch;
    double cs = ci*sh;
    double sc = si*ch;
    double ss = si*sh;

    M.initIdentity(4);

    if ( _initialRepeated )
    {
        dMij(M,i,i) = cj;
        dMij(M,j,i) =  sj*si;
        dMij(M,k,i) =  sj*ci;

        dMij(M,i,j) = sj*sh;
        dMij(M,j,j) = -cj*ss+cc;
        dMij(M,k,j) = -cj*cs-sc;

        dMij(M,i,k) = -sj*ch;
        dMij(M,j,k) =  cj*sc+cs;
        dMij(M,k,k) =  cj*cc-ss;
    }
    else
    {
        dMij(M,i,i) = cj*ch;
        dMij(M,j,i) = sj*sc-cs;
        dMij(M,k,i) = sj*cc+ss;

        dMij(M,i,j) = cj*sh;
        dMij(M,j,j) = sj*ss+cc;
        dMij(M,k,j) = sj*cs-sc;

        dMij(M,i,k) = -sj;
        dMij(M,j,k) = cj*si;
        dMij(M,k,k) = cj*ci;
    }
}


bool
Euler::legal(Euler::eulerOrder order)
{
    return (order & ~Legal) ? false : true;
}

Euler::eulerOrder Euler::order() const
{
    int foo = (_initialAxis == axisZ ?
               0x2000 : (_initialAxis == axisY ? 0x1000 : 0));

    if (_parityEven)
        foo |= 0x0100;
    if (_initialRepeated)
        foo |= 0x0010;
    if (_frameStatic)
        foo++;

    return (eulerOrder)foo;
}

void Euler::setOrder(Euler::eulerOrder p)
{
    set( p & 0x2000 ? axisZ : (p & 0x1000 ? axisY : axisX), // initial axis
         !(p & 0x1),        // static?
         !!(p & 0x100),    // permutation even?
         !!(p & 0x10));    // initial repeats?
}


void Euler::set(Euler::Axis axis,
                bool relative,
                bool parityEven,
                bool firstRepeats)
{
    _initialAxis = axis;
    _frameStatic = !relative;
    _parityEven  = parityEven;
    _initialRepeated = firstRepeats;
}

const Euler& Euler::operator= (const Euler &euler)
{
    x = euler.x;
    y = euler.y;
    z = euler.z;
    _initialAxis = euler._initialAxis;
    _frameStatic = euler._frameStatic;
    _parityEven  = euler._parityEven;
    _initialRepeated = euler._initialRepeated;
    return *this;
}

const Euler& Euler::operator= (const Matrix1D <double> &v)
{
    x = XX(v);
    y = YY(v);
    z = ZZ(v);
    return *this;
}


std::ostream& operator << (std::ostream &o, const Euler &euler)
{
    char a[3] = { 'X', 'Y', 'Z' };

    const char* r = euler.frameStatic() ? "" : "r";
    int i,j,k;
    euler.angleOrder(i,j,k);

    if ( euler.initialRepeated() )
        k = i;

    return o << "("
           << euler.x << " "
           << euler.y << " "
           << euler.z << " "
           << a[i] << a[j] << a[k] << r << ")";
}

double Euler::angleMod (double angle)
{
    angle = fmod( (angle),  (2. * M_PI));

    if (angle < -M_PI)
        angle += 2 * M_PI;
    if (angle > +M_PI)
        angle -= 2 * M_PI;

    return angle;
}

void
Euler::simpleXYZRotation (Matrix1D<double> &xyzRot,
                          const Matrix1D<double> &targetXyzRot)
{
    Matrix1D<double> d  = xyzRot - targetXyzRot;
    XX(xyzRot)  = XX(targetXyzRot) + angleMod(XX(d));
    YY(xyzRot)  = YY(targetXyzRot) + angleMod(YY(d));
    ZZ(xyzRot)  = ZZ(targetXyzRot) + angleMod(ZZ(d));
}

void
Euler::nearestRotation (Matrix1D<double> &xyzRot,
                        const Matrix1D<double> &targetXyzRot,
                        eulerOrder order)
{
    int i,j,k;
    Euler e (0,0,0, order);
    e.angleOrder(i,j,k);

    simpleXYZRotation(xyzRot, targetXyzRot);

    Matrix1D<double> otherXyzRot;
    XX(otherXyzRot) = M_PI+XX(xyzRot);
    YY(otherXyzRot) = M_PI-YY(xyzRot);
    ZZ(otherXyzRot) = M_PI+ZZ(xyzRot);

    simpleXYZRotation(otherXyzRot, targetXyzRot);

    Matrix1D<double> d  = xyzRot - targetXyzRot;
    Matrix1D<double> od = otherXyzRot - targetXyzRot;
    double dMag     = d.dotProduct(d);
    double odMag    = od.dotProduct(od);

    if (odMag < dMag)
    {
        xyzRot = otherXyzRot;
    }
}

void
Euler::makeNear (const Euler &target)
{
    Matrix1D<double> xyzRot ;
    toXYZVector(xyzRot);
    Euler targetSameOrder = Euler(target, order());
    Matrix1D<double> targetXyz;
    targetSameOrder.toXYZVector(targetXyz);

    nearestRotation(xyzRot, targetXyz, order());

    setXYZVector(xyzRot);
}


void Euler::eulerRotate (Matrix2D <double> &M,
                         const Matrix1D <double> &r)
{
    double cos_rz, sin_rz, cos_ry, sin_ry, cos_rx, sin_rx;
    double m00, m01, m02;
    double m10, m11, m12;
    double m20, m21, m22;

    cos_rz = cos (ZZ(r));
    cos_ry = cos (YY(r));
    cos_rx = cos (XX(r));
    sin_rz = sin (ZZ(r));
    sin_ry = sin (YY(r));
    sin_rx = sin (XX(r));

    m00 =  cos_rz *  cos_ry;
    m01 =  sin_rz *  cos_ry;
    m02 = -sin_ry;
    m10 = -sin_rz *  cos_rx + cos_rz * sin_ry * sin_rx;
    m11 =  cos_rz *  cos_rx + sin_rz * sin_ry * sin_rx;
    m12 =  cos_ry *  sin_rx;
    m20 = -sin_rz * -sin_rx + cos_rz * sin_ry * cos_rx;
    m21 =  cos_rz * -sin_rx + sin_rz * sin_ry * cos_rx;
    m22 =  cos_ry *  cos_rx;


    Matrix2D<double> P (M);

    dMij(M,0,0) = dMij(P,0,0) * m00 + dMij(P,1,0) * m01 + dMij(P,2,0) * m02;
    dMij(M,0,1) = dMij(P,0,1) * m00 + dMij(P,1,1) * m01 + dMij(P,2,1) * m02;
    dMij(M,0,2) = dMij(P,0,2) * m00 + dMij(P,1,2) * m01 + dMij(P,2,2) * m02;
    dMij(M,0,3) = dMij(P,0,3) * m00 + dMij(P,1,3) * m01 + dMij(P,2,3) * m02;

    dMij(M,1,0) = dMij(P,0,0) * m10 + dMij(P,1,0) * m11 + dMij(P,2,0) * m12;
    dMij(M,1,1) = dMij(P,0,1) * m10 + dMij(P,1,1) * m11 + dMij(P,2,1) * m12;
    dMij(M,1,2) = dMij(P,0,2) * m10 + dMij(P,1,2) * m11 + dMij(P,2,2) * m12;
    dMij(M,1,3) = dMij(P,0,3) * m10 + dMij(P,1,3) * m11 + dMij(P,2,3) * m12;

    dMij(M,2,0) = dMij(P,0,0) * m20 + dMij(P,1,0) * m21 + dMij(P,2,0) * m22;
    dMij(M,2,1) = dMij(P,0,1) * m20 + dMij(P,1,1) * m21 + dMij(P,2,1) * m22;
    dMij(M,2,2) = dMij(P,0,2) * m20 + dMij(P,1,2) * m21 + dMij(P,2,2) * m22;
    dMij(M,2,3) = dMij(P,0,3) * m20 + dMij(P,1,3) * m21 + dMij(P,2,3) * m22;
}
