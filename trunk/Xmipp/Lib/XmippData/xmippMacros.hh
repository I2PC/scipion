/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#ifndef _XMIPP_MACROS
   #define _XMIPP_MACROS

#ifndef _CYGWIN
   #ifdef __APPLE__
      #include <limits.h>
   #else
      #include <values.h>
   #endif
#endif
#ifndef MINFLOAT
   #define MINFLOAT -1e30
#endif
#ifndef MAXFLOAT
   #define MAXFLOAT  1e30
#endif

//@name Macros */
//@{

// Constants ---------------------------------------------------------------
/**@name Constants*/
//@{
/// PI = 3.1415926535897931
#ifndef PI
#define PI       3.14159265358979323846
#endif

/** In a comparison if two values are closer than this epsilon
    they are said to be the same. Actually set to 1e-6*/
#define XMIPP_EQUAL_ACCURACY 1e-6
//@}

// Numerical Macros --------------------------------------------------------
/**@name Numerical functions */
//@{
/** Absolute value.
    Valid for any kind of number (int, short, float, etc)
    \\ Ex: x=ABS(x); */
#ifndef ABS
#define ABS(x)   (((x)>=0)?(x):(-(x)))
#endif

/** Sign of.
    Valid for any kind of number (int, short, float, etc). It returns +1 or -1
    \\ Ex: if (SGN(x)==-1) cout << "x is negative" << endl; */
#ifndef SGN
#define SGN(x)   (((x)>=0)?1:-1)
#endif

/** Sign of, considering 0 as 0.
    Valid for any kind of number (int, short, float, etc). It returns +1 if
    the number is positive, -1 if the number is negative, and 0 if the
    number is 0.
    \\ Ex: if (SGN0(x)==-1) cout << "x is negative" << endl; */
#ifndef SGN0
#define SGN0(x)   (((x)>=0)? (((x)==0)? 0:1) :-1)
#endif

/** Minimum.
    Valid for any kind of numbers (int, short, float, etc).
    \\ Ex: min_val=MIN(x,y); */
#ifndef MIN
#define MIN(x,y) (((x)>=(y))?(y):(x))
#endif

/** Maximum.
    Valid for any kind of numbers (int, short, float, etc).
    \\ Ex: max_val=MAX(x,y); */
#ifndef MAX
#define MAX(x,y) (((x)>=(y))?(x):(y))
#endif

/** Round to next integer.
    Valid for any kind of numbers (int, short, float, etc). The result
    is of type integer.
    \\ Ex: a=ROUND(-0.8); ---> a=-1
    \\ Ex: a=ROUND(-0.2); ---> a= 0
    \\ Ex: a=ROUND( 0.2); ---> a= 0
    \\ Ex: a=ROUND( 0.8); ---> a= 1*/
#ifndef ROUND
#define ROUND(x) (((x)>0)? (int)((x)+0.5):(int)((x)-0.5))
#endif

/** Round to next larger integer.
    Valid for any kind of numbers (int, short, float, etc). The result
    is of type integer.
    \\ Ex: a=CEIL(-0.8); ---> a= 0
    \\ Ex: a=CEIL(-0.2); ---> a= 0
    \\ Ex: a=CEIL( 0.2); ---> a= 1
    \\ Ex: a=CEIL( 0.8); ---> a= 1*/
#define CEIL(x)  (((x)==(int)(x))? (int)(x):(((x)>0)? (int)((x)+1):(int)(x)))

/** Round to next smaller integer.
    Valid for any kind of numbers (int, short, float, etc). The result
    is of type integer.
    \\ Ex: a=FLOOR(-0.8); ---> a= -1
    \\ Ex: a=FLOOR(-0.2); ---> a= -1
    \\ Ex: a=FLOOR( 0.2); ---> a= 0
    \\ Ex: a=FLOOR( 0.8); ---> a= 0*/
#define FLOOR(x) (((x)==(int)(x))? (int)(x):(((x)>0)? (int)(x):(int)((x)-1)))

/** Return th efractional part of a value.
    The fractional part of 3.7 is 0.7 and of -3.7 is -0.7. */
#define FRACTION(x) ((x)-(int)(x)) 

/** Clip in a saturation fashion.
    CLIP is a macro which acts like a saturation curve, a value x is "clipped"
    to a range defined by x0 and xF, for example the output values for 
    the following x and CLIP(x,-2,2) would be
        x = ... -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 ...
   output = ... -2 -2 -2 -2 -2 -2 -2 -1 0 1 2 2 2 2 2 2 2 ... */
#define CLIP(x,x0,xF) (((x)<(x0))? (x0): (((x)>(xF))? (xF):(x)))

/** Wrapping for integers.
    intWRAP performs a wrapping in the integer set, when the cycle is finsihed
    it begins again. For example, for intWRAP(x,-2,2) would be
         x = ... -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8 ...
    output = ...  2 -2 -1  0  1  2 -2 -1  0  1  2 -2 -1  0  1  2 -2 ... */
#define intWRAP(x,x0,xF) (((x)>=(x0) && (x)<=(xF))? (x): ((x)<(x0))? \
                     ((x)-(int)(((x)-(x0)+1)/((xF)-(x0)+1)-1)*((xF)-(x0)+1)): \
                     ((x)-(int)(((x)-(xF)-1)/((xF)-(x0)+1)+1)*((xF)-(x0)+1)))

/** Wrapping for real numbers.
    realWRAP is used to keep a floating number between a range with a
    wrapping fashion.
    For instance, it is used in trigonometry to say that an angle of
    5*PI is the same as PI, ie, to keep an angle in the range 0...2*PI
    \\ Ex: Corrected_angle=realWRAP(angle,0,2*PI);*/
#define realWRAP(x,x0,xF) (((x)>=(x0) && (x)<=(xF))? (x): ((x)<(x0))? \
                     ((x)-(int)(((x)-(x0))/((xF)-(x0))-1)*((xF)-(x0))): \
                     ((x)-(int)(((x)-(xF))/((xF)-(x0))+1)*((xF)-(x0))))

/** Degrees to radians.
   Ex: angle_in_radians=DEG2RAD(ang_in_degrees);*/
#define DEG2RAD(d)         ((d)*PI/180)

/** Radians to degrees.
   Ex: angle_in_degrees=RAD2DEG(ang_in_radians);*/
#define RAD2DEG(r)         ((r)*180/PI)

/** Cosine in degrees.
   Ex: if (COSD(90)==0) cout << "This is in degrees!\n"; */
#define COSD(x) cos(PI*(x)/180.)

/** ArcCosine in degrees.
   Ex: if (ACOSD(0.5)==60) cout << "This is in degrees!\n"; */
#define ACOSD(x) acos((x))*180./PI

/** Sine in degrees.
   Ex: if (SIND(90)==1) cout << "This is in degrees!\n"; */
#define SIND(x) sin(PI*(x)/180.)

/** ArcSine in degrees.
   Ex: if (ASIND(0.5)==30.) cout << "This is in degrees!\n"; */
#define ASIND(x) asin((x))*180./PI

/** SINC function.
    The sinc function is defined as sin(PI*x)/(PI*x). */
#define SINC(x) (((x) < 0.0001 && (x) > -0.0001)? 1: sin(PI*(x))/(PI*(x)))

/** Returns next positive power of 2.
   It is supposed that the given number is positive although it's not
   needed to be an integer
   \\ Ex: next_power = NEXT_POWER_OF_2(1000);
   \\ ---> next_power=1024; */
#define NEXT_POWER_OF_2(x) pow(2,ceil(log((double)x)/log(2.0)) )

/** Linear interpolation. From low (when a=0) to high (when a=1)
    The following value is returned (equal to (a*h)+((1-a)*l) */
#define LIN_INTERP(a,l,h) ((l)+((h)-(l))*(a))

/** XOR. 
    Logical Xor.*/
#define XOR(a,b) (((a) && !(b)) || (!(a) && (b)))
//@}

// Miscellaneous -----------------------------------------------------------
/**@name Miscellaneous */
//@{
/** Speed up temporary variables.
    The following variables are provided:
    \begin{verbatim}
    float spduptmp0, spduptmp1, spduptmp2;
    int   ispduptmp0, ispduptmp1, ispduptmp2,
          ispduptmp3, ispduptmp4, ispduptmp5;
    \end{verbatim} */
#define SPEED_UP_temps \
    double spduptmp0, spduptmp1, spduptmp2, \
           spduptmp3, spduptmp4, spduptmp5, \
           spduptmp6, spduptmp7, spduptmp8; \
    int   ispduptmp0, ispduptmp1, ispduptmp2, \
          ispduptmp3, ispduptmp4, ispduptmp5;

/** Swap two values.
    It uses a temporal variable which must be of the same type as the two
    parameters */
#define SWAP(a,b,tmp) {\
    tmp=a; \
    a=b; \
    b=tmp; }

/** Starting point for Xmipp volume/image.
    Given a size (in some direction), this function returns the first index
    for a volume/image/array with this size. The formula is 
    \\-(int)((float)(size)/2.0)*/
#define FIRST_XMIPP_INDEX(size) -(int)((float)(size)/2.0)

/** Starting point for Xmipp volume/image.
    Given a size (in some direction), this function returns the first index
    for a volume/image/array with this size. The formula is 
    \\FIRST_XMIPP_INDEX(size)+(size)-1*/
#define LAST_XMIPP_INDEX(size) FIRST_XMIPP_INDEX(size)+(size)-1
//@}
#endif
