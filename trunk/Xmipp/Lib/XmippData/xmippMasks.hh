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
#ifndef _XMIPPMASKS_HH
   #define _XMIPPMASKS_HH
#include "xmippMatrices2D.hh"
#include "xmippMatrices3D.hh"
#include "xmippHistograms.hh"

/**@name Masks
*/
//@{
/*---------------------------------------------------------------------------*/
/* 1D Masks                                                                  */
/*---------------------------------------------------------------------------*/
/**@name 1D masks */
//@{
#define INNER_MASK 1
#define OUTSIDE_MASK 2
#define NO_ACTIVATE 0
#define ACTIVATE 1
/* RaisedCosine ............................................................ */
/** Creates a 1D RaisedCosine mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal RaisedCosines, and outside masks
    are 1-RaisedCosine.
    When entering, the mask is initialiazed to 0 and then the mask is created.
*/
void RaisedCosineMask(matrix1D<double> &mask, 
   double r1, double r2, int mode=INNER_MASK, double x0=0);

/* RaisedCrown ............................................................. */
/** Creates a 1D RaisedCrown mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created within the two
    the radii indicated with an extra region of <pix_width> pixels.
    The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal RaisedCrowns, and outside masks
    are 1-RaisedCrowns.
    When entering, the mask is initialiazed to 0 and then the mask is created.
*/
void RaisedCrownMask(matrix1D<double> &mask, 
   double r1, double r2, double pix_width, int mode=INNER_MASK,
   double x0=0);
//@}

/*---------------------------------------------------------------------------*/
/* 2D Masks                                                                  */
/*---------------------------------------------------------------------------*/
/**@name 2D masks */
//@{
/* Circular ................................................................ */
/** Creates a 2D circular mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void BinaryCircularMask(matrix2D<int> &mask, 
   double radius, int mode=INNER_MASK, double x0=0, double y0=0);

/* DWT Circular ............................................................ */
/** Creates a 2D DWT circular for already sized masks.
    The mask size must be a power of 2. The radius must be expressed in 
    pixel units corresponding to the size of the image. For instance,
    a 64x64 image might have a radius of 32 pixels to concentrate on the
    central part only. The mask is generated only for the desired masks.
    
    If the quadrant="xx" then 01, 10 and 11 are generated together*/
void BinaryDWTCircularMask(matrix2D<int> &mask,
   double radius, int smin, int smax, const string &quadrant);

/* Crown ................................................................... */
/** Creates a 2D crown mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created with the two radii indicated.
    The only two valid modes are INNER_MASK (by default, between the two
    radii) or OUTSIDE_MASK (the negative of the crown). It is supposed that
    R1 is smaller than R2.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void BinaryCrownMask(matrix2D<int> &mask, 
   double R1, double R2, int mode=INNER_MASK, double x0=0, double y0=0);

/* Frame ................................................................... */
/** Creates a 2D frame mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A square placed logically at (x0,y0), by default (0,0),
    is created with the two rectangular dimensions indicated.
    The only two valid modes are INNER_MASK (by default, between the two
    radii) or OUTSIDE_MASK (the negative of the crown). It is supposed that
    R1 is smaller than R2.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void BinaryFrameMask(matrix2D<int> &mask, 
   int Xrect, int Yrect, int mode=INNER_MASK, double x0=0, double y0=0);

/* Gaussian ................................................................ */
/** Creates a 2D gaussian mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal gaussians, and outside masks
    are 1-gaussian.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void GaussianMask(matrix2D<double> &mask, 
   double sigma, int mode=INNER_MASK, double x0=0, double y0=0);

/* RaisedCosine ............................................................ */
/** Creates a 2D RaisedCosine mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal RaisedCosines, and outside masks
    are 1-RaisedCosine.
    When entering, the mask is initialiazed to 0 and then the mask is created.
*/
void RaisedCosineMask(matrix2D<double> &mask, 
   double r1, double r2, int mode=INNER_MASK, double x0=0, double y0=0);

/* RaisedCrown ............................................................. */
/** Creates a 2D RaisedCrown mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created within the two
    the radii indicated with an extra region of <pix_width> pixels.
    The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal RaisedCrowns, and outside masks
    are 1-RaisedCrowns.
    When entering, the mask is initialiazed to 0 and then the mask is created.
*/
void RaisedCrownMask(matrix2D<double> &mask, 
   double r1, double r2, double pix_width, int mode=INNER_MASK,
   double x0=0, double y0=0);

/* Blackman window ......................................................... */
/** 2D blackman window.
    It receives no parameter. */
void BlackmanMask(matrix2D<double> &mask, int mode=INNER_MASK,
   double x0=0, double y0=0);

/* Sinc .................................................................... */
/** Creates a 2D sinc mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal gaussians, and outside masks
    are 1-sinc.
    When entering the mask is initialiazed to 0 and then the mask is created.
    
    Remind that sinc(w*n) is zero at n=1/w;
*/
void SincMask(matrix2D<double> &mask, 
   double omega, int mode=INNER_MASK, double x0=0, double y0=0);

/** Creates a 2D sinc-blackman mask, the mask is resized.
    This function returns a sinc mask windowed by a Blackman window.
    The window is designed to cover a certain power of the sinc
*/
void SincBlackmanMask(matrix2D<double> &mask, 
   double omega, double power_percentage,
   int mode=INNER_MASK, double x0=0, double y0=0);

/* Neighborhood masks ....................................................... */
/** Creates a 3x3 mask with value (1 by default) for those 4-neighbors of the 
    central point (0 otherwise). The parameter center controls whether the center 
    pixel is set to 1 or not 
*/
void mask2D_4neig(matrix2D<int> &mask, int value=1, int center= NO_ACTIVATE);

/** Creates a 3x3 mask with value1 for those 4-neighbors of the central point and
    value2 for the 8 neighbors. The parameter center controls whether
    the center pixel is set to 1 or not */
void mask2D_8neig(matrix2D<int> &mask, int value1=1, int value2=1, 
                   int center= NO_ACTIVATE);
//@}

/*---------------------------------------------------------------------------*/
/* 3D Masks                                                                  */
/*---------------------------------------------------------------------------*/
/**@name 3D masks */
//@{
/* Spherical ............................................................... */
/** Creates a 3D spherical mask for already sized masks.
    The mask is supposed to be already resized and with its logical origin
    defined. A sphere placed logically at (z0,x0,y0),
    by default (0,0,0), is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void BinarySphericalMask(matrix3D<int> &mask,
   double radius, int mode=INNER_MASK, double x0=0, double y0=0, int z0=0);

/* DWT Spherical ............................................................ */
/** Creates a 3D DWT spherical for already sized masks.
    The mask size must be a power of 2. The radius must be expressed in 
    pixel units corresponding to the size of the image. For instance,
    a 64x64x64 image might have a radius of 32 pixels to concentrate on the
    central part only.
    
    If quadrant=xxx then 001,010,011,100,101,110 and 111 are generated
    together*/
void BinaryDWTCircularMask(matrix3D<int> &mask,
   double radius, int smin, int smax, const string &quadrant);

/* Crown ................................................................... */
/** Creates a 3D crown mask for already sized masks.
    The mask is supposed to be already resized and with its logical origin
    defined. A sphere placed logically at (z0,x0,y0),
    by default (0,0,0), is created with
    the radii indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void BinaryCrownMask(matrix3D<int> &mask,
   double R1, double R2, int mode=INNER_MASK, double x0=0, double y0=0, int z0=0);

/* Cylinder ................................................................ */
/** Creates a 3D Cylinder mask for already sized masks.
    The mask is supposed to be already resized and with its logical origin
    defined. A cylinder placed logically at (z0,x0,y0),
    by default (0,0,0), is created with
    the radius and height indicated. The only two valid modes are INNER_MASK
    (by default) or OUTSIDE_MASK.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void BinaryCylinderMask(matrix3D<int> &mask,
   double R, double H, int mode=INNER_MASK, double x0=0, double y0=0, int z0=0);

/* Frame ................................................................... */
/** Creates a 3D frame mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A square placed logically at (x0,y0,z0), by default (0,0,0),
    is created with the two rectangular dimensions indicated.
    The only two valid modes are INNER_MASK (by default, between the two
    radii) or OUTSIDE_MASK (the negative of the crown). It is supposed that
    R1 is smaller than R2.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void BinaryFrameMask(matrix3D<int> &mask, 
   int Xrect, int Yrect, int Zrect, int mode=INNER_MASK,
   double x0=0, double y0=0, double z0=0);

/* Gaussian ................................................................ */
/** Creates a 3D gaussian mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0), by default (0,0),
    is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal gaussians, and outside masks
    are 1-gaussian.
    When entering the mask is initialiazed to 0 and then the mask is created.
*/
void GaussianMask(matrix3D<double> &mask, 
   double sigma, int mode=INNER_MASK, double x0=0, double y0=0, double z0=0);

/* RaisedCosine ............................................................ */
/** Creates a 3D RaisedCosine mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0,z0), by default (0,0,0),
    is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal RaisedCosines, and outside masks
    are 1-RaisedCosine.
    When entering, the mask is initialiazed to 0 and then the mask is created.
*/
void RaisedCosineMask(matrix3D<double> &mask, 
   double r1, double r2, int mode=INNER_MASK, double x0=0, double y0=0,
   double z0=0);

/* RaisedCrown ............................................................. */
/** Creates a 3D RaisedCrown mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0,z0), by default (0,0,0),
    is created within the two
    the radii indicated with an extra region of <pix_width> voxels.
    The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal RaisedCrowns, and outside masks
    are 1-RaisedCrowns.
    When entering, the mask is initialiazed to 0 and then the mask is created.
*/
void RaisedCrownMask(matrix3D<double> &mask, 
   double r1, double r2, double pix_width, int mode=INNER_MASK,
   double x0=0, double y0=0, double z0=0);

/* Blackman window ......................................................... */
/** 3D blackman window.
    It receives no parameter. */
void BlackmanMask(matrix3D<double> &mask, int mode=INNER_MASK,
   double x0=0, double y0=0, double z0=0);

/* Sinc .................................................................... */
/** Creates a 3D sinc mask for already sized masks.
    The mask is supposed to be resized and with its logical origin already
    set. A circle placed logically at (x0,y0,z0), by default (0,0,0),
    is created with
    the radius indicated. The only two valid modes are INNER_MASK (by default)
    or OUTSIDE_MASK. Inner mask are normal gaussians, and outside masks
    are 1-sinc.
    When entering the mask is initialiazed to 0 and then the mask is created.
    
    Remind that sinc(w*t) is zero at t=1/w;
*/
void SincMask(matrix3D<double> &mask, 
   double omega, int mode=INNER_MASK, double x0=0, double y0=0, double z0=0);

/** Creates a 3D sinc-blackman mask, the mask is resized.
    This function returns a sinc mask windowed by a Blackman window.
    The window is designed to cover a certain power of the sinc
*/
void SincBlackmanMask(matrix3D<double> &mask, 
   double omega, double power_percentage,
   int mode=INNER_MASK, double x0=0, double y0=0, double z0=0);

/* Neighborhood masks ....................................................... */
/** Creates a 3x3x3 mask with value (1 by default) for those 6-neighbors of the 
    central point (0 otherwise). The parameter center controls whether the center 
    pixel is set to 1 or not */
    
void mask3D_6neig(matrix3D<int> &mask, int value=1, int center= NO_ACTIVATE);

/** Creates a 3x3x3 mask with value1 (1 by default) for those 6-neighbors and 
    value2 for the  18 neighbors of the central point (0 otherwise). The parameter
    center controls whether the center pixel is set to 1 or not */    
void mask3D_18neig(matrix3D<int> &mask, int value1=1, int value2=1,
 		    int center= NO_ACTIVATE);

/** Creates a 3x3x3 mask with value1 (1 by default) for those 6-neighbors, value2 
    for the  18 neighbors and  value3 for the  26 neighbors of the central point 
    (0 otherwise). The parameter center controls whether the center pixel is set 
    to 1 or not */    
void mask3D_26neig(matrix3D<int> &mask, int value1=1, int value2=1, int value3=1,
 		    int center= NO_ACTIVATE);
//@}

/*---------------------------------------------------------------------------*/
/* Mask Type                                                                 */
/*---------------------------------------------------------------------------*/
/** Parameters for a general Mask.
    This class contains all parameters needed to generate masks. The class
    can read parameters from the command line.
    
    To read a mask from a file within a program do the following
    \begin{verbatim}
    Mask_Params Mask;
    Mask.type=READ_MASK;
    Mask.fn_mask="...";
    Mask.generate_2Dmask();
    
    Mask.apply(input_matrix2D,output_matrix2D);
    \end{verbatim}
    
    To generate a geometric mask within a program do the following:
    \begin{verbatim}
    Mask_Params Mask;
    
    // Define an spherical mask of radius 32 (the active part is
    // within the sphere)
    Mask.type=BINARY_CIRCULAR_MASK;
    Mask.mode=INNER_MASK;
    Mask.R1=32;
    
    // resize the mask after this pattern
    Mask.resize(input_matrix2D);
    
    // Really generate the mask. It is stored internally
    Mask.generate_2Dmask();

    // Apply the mask to some image
    Mask.apply_mask(input_matrix2D,output_matrix2D);
    \end{verbatim}
    */
class Mask_Params {
public:
#define NO_MASK               	  0
#define BINARY_CIRCULAR_MASK  	  1
#define BINARY_CROWN_MASK     	  2
#define BINARY_CYLINDER_MASK  	  3
#define BINARY_FRAME_MASK     	  4
#define GAUSSIAN_MASK         	  5
#define RAISED_COSINE_MASK    	  6
#define BLACKMAN_MASK         	  7
#define SINC_MASK             	  8
#define SINC_BLACKMAN_MASK    	  9
#define READ_MASK             	 10
#define RAISED_CROWN_MASK        11
#define BINARY_DWT_CIRCULAR_MASK 12
/** Mask Type.
    The only valid types are BINARY_CIRCULAR_MASK,
    BINARY_CROWN_MASK, BINARY_CYLINDER_MASK, BINARY_FRAME_MASK,
    GAUSSIAN_MASK, RAISED_COSINE_MASK, BLACKMAN_MASK, SINC_MASK,
    SINC_BLACKMAN_MASK, READ_MASK, RAISED_CROWN_MASK */
    int type;

/** Mode.
    The valid modes are INNER_MASK and OUTSIDE_MASK. */
    int mode;

/** Radius 1.
    Radius for Circular and Cylinder masks and R1 for crowns and raised cosine.
    */
    double R1;

/** Radius 2.
    R2 for crowns and raised cosine. */
    double R2;

/** Pixel width.
    For raised crowns. */
    double pix_width;

/** Height.
    Height for cylinders. */
    double H;

/** Sigma.
    Sigma for gaussians. */
    double sigma;

/** Omega.
    Frequency for sincs. */
    double omega;

/** Rectangular X dimension */
    int Xrect;

/** Rectangular Y dimension */
    int Yrect;

/** Rectangular Z dimension */
    int Zrect;

/** Z origin */
    double z0;
    
/** Y origin */
    double y0;
    
/** X origin */
    double x0;
    
/** Minimum scale for DWT masks. */
    int smin;

/** Maximum scale for DWT masks. */
    int smax;

/** Quadrant.
    If it is empty then all, except 000, are generated. */
    string quadrant;

/** Filename from which the mask is read, if it is the case */
    FileName fn_mask;

/** Geometrix transformation matrix for the mask */
    matrix2D<double> mask_geo;

/** Allowed data types. */
    int allowed_data_types;
/** 1D integer mask. */
    matrix1D<int> imask1D;
/** 2D integer mask. */
    matrix2D<int> imask2D;
/** 3D integer mask. */
    matrix3D<int> imask3D;
/** 1D double mask */
    matrix1D<double> dmask1D;
/** 2D double mask. */
    matrix2D<double> dmask2D;
/** 3D double mask. */
    matrix3D<double> dmask3D;

public:
    #define INT_MASK    1
    #define DOUBLE_MASK 2
    #define ALL_KINDS   INT_MASK | DOUBLE_MASK
/** Constructors. Allowed data types are ALL_KINDS, INT_MASK and DOUBLE_MASK
    used with | .*/
    Mask_Params(int _allowed_data_type=ALL_KINDS);

/** Clear. */
    void clear();

/** Read from command line.
    An exception is thrown if the read mask is not of an allowed type. */
    void read(int argc, char **argv) _THROW;

/** Show. */
    void show() const;

/** Usage. */
    void usage() const;

/** Save 1D mask as a text file */
    void write_1Dmask(const FileName &fn);

/** Save 2D mask as an ImageXmipp. */
    void write_2Dmask(const FileName &fn);

/** Save 3D mask as an VolumeXmipp. */
    void write_3Dmask(const FileName &fn);

/** Return the type of the mask. INT_MASK, DOUBLE_MASK */
    int datatype() {
       if (type==BINARY_CIRCULAR_MASK || type==BINARY_CROWN_MASK ||
           type==BINARY_CYLINDER_MASK || type==BINARY_FRAME_MASK ||
           type==NO_MASK              || type==READ_MASK         ||
	   type==BINARY_DWT_CIRCULAR_MASK)
           return INT_MASK;
       else if (type==GAUSSIAN_MASK   || type==RAISED_COSINE_MASK ||
                type==SINC_MASK       || type==SINC_BLACKMAN_MASK ||
                type==BLACKMAN_MASK   || type==RAISED_CROWN_MASK)
           return DOUBLE_MASK;
    }

/** Resize and set Xmipp origin */
    void resize(int Xdim);
/** Resize and set Xmipp origin. */
    void resize(int Ydim, int Xdim);
/** Resize and set Xmipp origin. */
    void resize(int Zdim, int Ydim, int Xdim);
/** Resize after a pattern */
template <class T>
    void resize(const matrix1D<T> &m);
/** Resize after a pattern */
template <class T>
    void resize(const matrix2D<T> &m);
/** Resize after a pattern */
template <class T>
    void resize(const matrix3D<T> &m);

/** Generate mask for a resized signal.
    It is supposed that the image is already resized and with its logical
    origin set. */
    void generate_1Dmask();

/** Generate mask for an empty signal.*/
    void generate_1Dmask(int Xdim)
       {resize(Xdim); generate_1Dmask();}
    
/** Generate mask for a signal following a pattern. */
template <class T>
   void generate_1Dmask(const matrix1D<T> &m)
       {resize(m); generate_1Dmask();}

/** Generate mask for a resized image.
    It is supposed that the image is already resized and with its logical
    origin set. */
    void generate_2Dmask();

/** Generate mask for an empty image.*/
    void generate_2Dmask(int Ydim, int Xdim)
       {resize(Ydim, Xdim); generate_2Dmask();}
    
/** Generate mask for an image following a pattern. */
template <class T>
   void generate_2Dmask(const matrix2D<T> &m)
       {resize(m); generate_2Dmask();}

/** Generate mask for a resized volume.
    It is supposed that the image is already resized and with its logical
    origin set. */
    void generate_3Dmask();

/** Generate mask for an empty volume. */
    void generate_3Dmask(int Zdim, int Ydim, int Xdim)
       {resize(Zdim, Ydim, Xdim); generate_3Dmask();}
    
/** Generate mask for an image following a pattern. */
template <class T>
   void generate_3Dmask(const matrix3D<T> &m)
       {resize(m); generate_3Dmask();}

/** Apply mask to signal.
    subs_val is the substitute value in case of binary masks*/
template <class T>
   void apply_mask(const matrix1D<T> &I, matrix1D<T> &result, T subs_val=0);

/** Apply mask to image.
    subs_val is the substitute value in case of binary masks*/
template <class T>
   void apply_mask(const matrix2D<T> &I, matrix2D<T> &result, T subs_val=0, const bool &apply_geo=false);

/** Apply mask to volume.
    subs_val is the substitute value in case of binary masks*/
template <class T>
   void apply_mask(const matrix3D<T> &I, matrix3D<T> &result, T subs_val=0);

/** Get binary 1D mask. */
matrix1D<int>    & get_binary_mask1D() {return imask1D;}
/** Get continuous 1D mask. */
matrix1D<double> & get_cont_mask1D()   {return dmask1D;}
/** Get binary 2D mask. */
matrix2D<int>    & get_binary_mask2D() {return imask2D;}
/** Get continuous 2D mask. */
matrix2D<double> & get_cont_mask2D()   {return dmask2D;}
/** Get binary 3D mask. */
matrix3D<int>    & get_binary_mask3D() {return imask3D;}
/** Get continuous 3D mask. */
matrix3D<double> & get_cont_mask3D()   {return dmask3D;}
/** Force to be continuous.
    This function is used when you need a binary mask as a double matrix. */
void force_to_be_continuous() {
   if (datatype()==INT_MASK)
      {type_cast(imask1D,dmask1D);
       type_cast(imask2D,dmask2D);
       type_cast(imask3D,dmask3D);}
}

/** Force to be binary.
    This function is used when you need a double mask as a binary matrix. */
void force_to_be_binary() {
   if (datatype()==DOUBLE_MASK)
      {type_cast(dmask1D,imask1D);
       type_cast(dmask2D,imask2D);
       type_cast(dmask3D,imask3D);}
}
};

/*---------------------------------------------------------------------------*/
/* Mask tools                                                                */
/*---------------------------------------------------------------------------*/
/**@name Tools
   All Mask tools work only in the overlapping area of the given image/volume
   and the mask in logical coordinates. Ie, if you have a mask defined from
   -2 to 2 and you apply it to an image defined from 0 to 63 then only
   those values of the mask between 0 and 2 will be applied. The rest of
   the image will remain untouched. This region where the mask is active
   within the overlapping area will be called in this documentation:
   active area.
 */
//@{
/** Compute statistics in the active area (2D).
    Only the statistics for values in the overlapping between the mask
    and the image for those the mask is not 0 are computed.*/
template <class T>
   void compute_stats_within_binary_mask(const matrix2D<int> &mask,
      const matrix2D<T> &m, T &min_val, T &max_val, double &avg, double &stddev);

/** Apply geometric transformation to a binary mask  */
void apply_geo_binary_2D_mask(matrix2D<int> &mask, const matrix2D<double> &A);

/** Apply geometric transformation to a continuous mask  */
void apply_geo_cont_2D_mask(matrix2D<double> &mask, const matrix2D<double> &A);

/** Apply binary mask to an image (1D).
    The image values for which the input mask is 0 are set to <subs_val>.
    The input and output matrices can be the same ones.
    Only the overlapping values are affected by the mask */
template <class T>
   void apply_binary_mask(const matrix1D<int> &mask, const matrix1D<T> &m_in,
      matrix1D<T> &m_out, T subs_val=(T)0);

/** Apply continuous mask to an image (1D).
    The image is multiplied by the mask.
    The input and output matrices can be the same ones.
    Only the overlapping values are affected by the mask */
template <class T>
   void apply_cont_mask(const matrix1D<double> &mask, const matrix1D<T> &m_in,
      matrix1D<T> &m_out);

/** Apply binary mask to an image (2D).
    The image values for which the input mask is 0 are set to <subs_val>.
    The input and output matrices can be the same ones.
    Only the overlapping values are affected by the mask */
template <class T>
   void apply_binary_mask(const matrix2D<int> &mask, const matrix2D<T> &m_in,
      matrix2D<T> &m_out, T subs_val=(T)0);

/** Apply continuous mask to an image (2D).
    The image is multiplied by the mask.
    The input and output matrices can be the same ones.
    Only the overlapping values are affected by the mask */
template <class T>
   void apply_cont_mask(const matrix2D<double> &mask, const matrix2D<T> &m_in,
      matrix2D<T> &m_out);

/** Compute statistics in the active area (3D).
    Only the statistics for values in the overlapping between the mask
    and the volume for those the mask is not 0 are computed.*/
template <class T>
   void compute_stats_within_binary_mask(const matrix3D<int> &mask,
      const matrix3D<T> &m, T &min_val, T &max_val, double &avg, double &stddev);

/** Apply mask to an volume (3D).
    The volume values for which the input mask is 0 are set to 0. The input
    and output volumes can be the same ones.*/    
template <class T>
   void apply_binary_mask(const matrix3D<int> &mask, const matrix3D<T> &m_in,
      matrix3D<T> &m_out, T subs_val=(T)0);

/** Apply continuous mask to an image (3D).
    The image is multiplied by the mask.
    The input and output matrices can be the same ones.
    Only the overlapping values are affected by the mask */
template <class T>
   void apply_cont_mask(const matrix3D<double> &mask, const matrix3D<T> &m_in,
      matrix3D<T> &m_out);

/** Compute histogram inside mask within its minimum and maximum value (2D).
    Given a matrix as input, this function returns its histogram within
    the minimum and maximum of the matrix inside the mask, in this way all
    the values
    in the matrix are counted. The matrix can be of any numerical type
    (short int, int, double, ...). The number of steps must always be
    given. */
template <class T>
   void compute_hist_within_binary_mask(const matrix2D<int> &mask,
      matrix2D<T> &v, histogram1D &hist, int no_steps);

/** Compute histogram inside mask within two values (2D).
    Given a matrix as input, this function returns the histogram of values
    inside the mask within
    two values, the matrix values outside this range are not counted.
    This can be used to avoid the effect of outliers which causes a
    "compression" in the histogram.
    The matrix can be of any numerical type
    (short int, int, double, ...). The number of steps must always be
    given. */
template <class T>
   void compute_hist_within_binary_mask(const matrix2D<int> &mask,
      const matrix2D<T> &v, histogram1D &hist, T min, T max, int no_steps);

/** Compute histogram inside mask within its minimum and maximum value (3D).
    Given a volume as input, this function returns the histogram of values
    inside the mask within
    the minimum and maximum of the volume, in this way all the values
    in the volume are counted. The volume can be of any numerical type
    (short int, int, double, ...). The number of steps must always be
    given. */
template <class T>
   void compute_hist_within_binary_mask(const matrix3D<int> &mask,
      matrix3D<T> &v, histogram1D &hist, int no_steps);

/** Compute histogram inside mask within two values (3D).
    Given a volume as input, this function returns the histogram of the
    values inside the mask within
    two values, the volume values outside this range are not counted.
    This can be used to avoid the effect of outliers which causes a
    "compression" in the histogram.
    The volume can be of any numerical type
    (short int, int, double, ...). The number of steps must always be
    given. */
template <class T>
   void compute_hist_within_binary_mask(const matrix3D<int> &mask,
      const matrix3D<T> &v, histogram1D &hist, T min, T max, int no_steps);

#define COUNT_ABOVE 1
#define COUNT_BELOW 2
#define COUNT_BETWEEN 3
/** Count pixels/voxels with mask and above a threshold.
    Those pixels within the mask with a value greater or equal than
    a threshold are counted.
    This function makes a call to \Ref{count_with_mask} */
#define count_with_mask_above(mask, m, th) \
    count_with_mask(mask,m,COUNT_ABOVE,th,0);

/** Count pixels/voxels with mask and below a threshold.
    Those pixels within the mask with a value smaller or equal than
    a threshold are counted.
    This function makes a call to \Ref{count_with_mask} */
#define count_with_mask_below(mask, m, th) \
    count_with_mask(mask,m,COUNT_BELOW,th,0);

/** Count pixels/voxels with mask and between two thresholds.
    Those pixels within the mask with a value greater or equal than
    th1 and smaller or equal than th2 are counted.
    This function makes a call to \Ref{count_with_mask} */
#define count_with_mask_between(mask, m, th1, th2) \
    count_with_mask(mask,m,COUNT_BETWEEN,th1,th2);

/** Count pixels with mask and threshold.
    This function returns the number of pixels in the ACTIVE area of an
    image with a value:
    \\COUNT_ABOVE: greater or equal than th1
    \\COUNT_BELOW: smaller or equal than th1
    \\COUNT_BETWEEN: smaller or equal than th1 and greater or equal than th2
    For complex matrices the absolute value is compared. */
template <class T>
   int count_with_mask(const matrix2D<int> &mask,
      const matrix2D<T> &m, int mode, double th1, double th2);
int count_with_mask(const matrix2D<int> &mask,
   const matrix2D<double_complex> &m, int mode, double th1, double th2);

/** Count voxels with mask and threshold.
    This function returns the number of voxels in the ACTIVE area of an
    volume with a value:
    \\COUNT_ABOVE: greater or equal than th1
    \\COUNT_BELOW: smaller or equal than th1
    \\COUNT_BETWEEN: smaller or equal than th1 and greater or equal than th2
    For complex matrices the absolute value is compared. */
template <class T>
   int count_with_mask(const matrix3D<int> &mask,
      const matrix3D<T> &m, int mode, double th1, double th2);
int count_with_mask(const matrix3D<int> &mask,
   const matrix3D<double_complex> &m, int mode, double th1, double th2);

/** Invert binary mask (2D).
    0's are converted in 1's and viceversa*/
template <class T>
   void invert_binary_mask(matrix2D<T> &mask);

/** Invert binary mask (3D).
    0's are converted in 1's and viceversa*/
template <class T>
   void invert_binary_mask(matrix3D<T> &mask);

/** Range adjust within binary mask.
    Make the grey values of m2 fit, in L2 sense, with those in m1. Only the
    voxels within the mask are used to compute the linear transformation.
    If no mask is provided then all voxels are used. */
void range_adjust_within_mask(const matrix2D<double> *mask,
   const matrix2D<double> &m1, matrix2D<double> &m2);

/** Range adjust within binary mask.
    Make the grey values of m2 fit, in L2 sense, with those in m1. Only the
    voxels within the mask are used to compute the linear transformation.
    If no mask is provided then all voxels are used. */
void range_adjust_within_mask(const matrix3D<double> *mask,
   const matrix3D<double> &m1, matrix3D<double> &m2);
//@}
//@}
#endif
