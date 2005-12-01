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
/* ------------------------------------------------------------------------- */
/* VOLUMES                                                                   */
/* ------------------------------------------------------------------------- */
/* When adding functions to this module, don't forget to add them in the
   MultidimInstantiation.cc, too */

#ifndef _XMIPPMATRICES3D_HH
#   define _XMIPPMATRICES3D_HH

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */
#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <complex>

#include "xmippMatrices1D.hh"
#include "xmippMatrices2D.hh"
#include "xmippFuncs.hh"
#include "xmippArgs.hh"

#define maT  matrix3D<T>
#define maT1 matrix3D<T1>
#undef  maTC
#define maTC matrix3D< complex<double> >

/* ************************************************************************* */
/* FORWARD DEFINITIONS                                                       */
/* ************************************************************************* */
#include "Src/MultidimFriends.inc"

template <class T>
   void apply_geom(VT &V2, matrix2D<double> A, const VT &V1, bool inv,
      bool wrap);

template <class T>
   void apply_geom_Bspline(VT &V2, matrix2D<double> A, const VT &V1,
      int Splinedegree, bool inv, bool wrap, T outside=0);

/* ************************************************************************* */
/* CLASS DEFINITION AND PROTOTYPES                                           */
/* ************************************************************************* */
/**@name Xmipp Volumes*/
//@{
/* Speed up ------------------------------------------------------------- */
/**@name Speed up macros
   This macros are defined to allow high speed in critical parts of
   your program. They shouldn't be used systematically as usually there
   is no checking on the correctness of the operation you are performing.
   Speed comes from three facts: first, they are macros and no function
   call is performed (although most of the critical functions are
   inline functions), there is no checking on the correctness of the
   operation (it could be wrong and you are not warned of it), and
   destination vectors are not returned saving time in the copy
   constructor and in the creation/destruction of temporary vectors.*/
//@{
/**@name Size and shape
    Although they are not defined here you can also use STARTINGX and
    FINISHINGX (defined for matrix1D), or STARTINGY and FINISHINGY
    (defined for matrix2D)*/
//@{
/** TRUE if both arrays have the same shape.
    Two arrays have the same shape if they have the same size and the
    same starting point. Be aware that this is a macro which simplifies to
    a boolean. */
#define SAME_SHAPE3D(v1,v2) \
    (XSIZE(v1)==XSIZE(v2) && \
     YSIZE(v1)==YSIZE(v2) && \
     ZSIZE(v1)==ZSIZE(v2) && \
     STARTINGX(v1)==STARTINGX(v2) &&\
     STARTINGY(v1)==STARTINGY(v2) && \
     STARTINGZ(v1)==STARTINGZ(v2))

/** Returns the first valid logical Z index.
    \\Ex: int orgZ=STARTINGZ(V);*/
#define STARTINGZ(m)  ((m).zinit)

/** Returns the last valid logical Z index.
    \\Ex: int finZ=FINISHINGZ(V);*/
#define FINISHINGZ(m) ((m).zinit+(m).zdim-1)

/** Access to Z dimension (size).
    This is a macro equivalent to \Ref{SliNo()}
    \\ Ex:
    \begin{verbatim}
    // Set to 0 1 element out of 8
    for (int k=0; k<ZSIZE(V); k +=2)
       for (int i=0; i<YSIZE(V); i +=2)
           for (int j=0; j<XSIZE(V); j +=2)
               DIRECT_VOL_ELEM(V,k,i,j)=0;
   \end{verbatim}*/
#define ZSIZE(V) ((V).zdim)

/** Access to XY dimension (size).
    Notice that XYSIZE(*this)=XSIZE(*this)*YSIZE(*this). But, this
    value is stored with the volume in order to make voxel access
    faster. */
#define XYSIZE(V) ((V).xydim)

/** For all elements in the array.
    This macro is used to generate loops for the volume in an easy way.
    It defines internal indexes 'k','i' and 'j' which ranges the volume
    using its mathematical definition (ie, logical access).
    \\Ex:
    \begin{verbatim}
    FOR_ALL_ELEMENTS_IN_MATRIX3D(V) {
       cout << m(k,i,j) << " ";
    }
    \end{verbatim} */
#define FOR_ALL_ELEMENTS_IN_MATRIX3D(V) \
    for (int k=STARTINGZ(V); k<=FINISHINGZ(V); k++) \
       for (int i=STARTINGY(V); i<=FINISHINGY(V); i++) \
          for (int j=STARTINGX(V); j<=FINISHINGX(V); j++)

/** For all elements in the array between corners.
    This macro is used to generate loops for a volume in an easy manner.
    Then ZZ(r), YY(r) and XX(r) range from
    (int) ZZ(corner1) to (int)ZZ(corner2),
    (int) YY(corner1) to (int)YY(corner2), (int) XX(corner1) to
    (int) XX(corner2) (included limits) respectively. Notice that corner1
    and corner2 need only be matrix1D. 
    \\Ex:
    \begin{verbatim}
    matrix1D<double> corner1(3), corner2(3), r(3);
    XX(corner1)=-1; XX(corner2)=1;
    YY(corner1)=-2; YY(corner2)=2;
    ZZ(corner1)=-3; ZZ(corner2)=3;
    FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1,corner2) {
       cout << v(r) << " ";
    }
    \end{verbatim}*/
#define FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1,corner2) \
    for (ZZ(r)=ZZ((corner1)); ZZ(r)<=ZZ((corner2)); ZZ(r)++) \
       for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); YY(r)++) \
          for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); XX(r)++)

/** For all elements in common.
    This macro is used to generate loops for all the elements logically
    in common between two volumes in an easy manner.
    Then k, i and j (locally defined) range from
    MAX(STARTINGZ(V1),STARTINGZ(V2)) to MIN(FINISHINGZ(V1),FINISHINGZ(V2)),
    MAX(STARTINGY(V1),STARTINGY(V2)) to MIN(FINISHINGY(V1),FINISHINGY(V2)),
    MAX(STARTINGX(V1),STARTINGX(V2)) to MIN(FINISHINGX(V1),FINISHINGX(V2))
    (included limits) respectively. You need to define SPEED_UP_temps.
    \\Ex:
    \begin{verbatim}
    matrix3D<double> V1(10,10,10), V2(20,20,20);
    V1.set_Xmipp_origin();
    V2.set_Xmipp_origin();
    FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(V1,V2) {
       ...
    }
    \end{verbatim}*/
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(V1,V2) \
    ispduptmp0=MAX(STARTINGZ(V1), STARTINGZ(V2)); \
    ispduptmp1=MIN(FINISHINGZ(V1),FINISHINGZ(V2)); \
    ispduptmp2=MAX(STARTINGY(V1), STARTINGY(V2)); \
    ispduptmp3=MIN(FINISHINGY(V1),FINISHINGY(V2)); \
    ispduptmp4=MAX(STARTINGX(V1), STARTINGX(V2)); \
    ispduptmp5=MIN(FINISHINGX(V1),FINISHINGX(V2)); \
    for (int k=ispduptmp0; k<=ispduptmp1; k++) \
       for (int i=ispduptmp2; i<=ispduptmp3; i++) \
          for (int j=ispduptmp4; j<=ispduptmp5; j++)

/** For all direct elements in the array.
    This macro is used to generate loops for the volume in an easy way.
    It defines internal indexes 'k','i' and 'j' which ranges the volume
    using its physical definition.
    \\Ex:
    \begin{verbatim}
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(V) {
       cout << DIRECT_VOL_ELEM(m,k,i,j) << " ";
    }
    \end{verbatim} */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(V) \
    for (int k=0; k<ZSIZE(V); k++) \
       for (int i=0; i<YSIZE(V); i++) \
          for (int j=0; j<XSIZE(V); j++)
//@}

/**@name Memory access*/
//@{
/** Volume element: Logical access.
    \\ Ex: VOL_ELEM(V,-1,-2,1)=1;
    \\ Ex: val=VOL_ELEM(V,-1,-2,1);*/
#define VOL_ELEM(V,k,i,j) \
   DIRECT_VOL_ELEM(V,(k)-STARTINGZ(V),(i)-STARTINGY(V),(j)-STARTINGX(V))

/** Volume element: Physical access.
    Be careful because this is physical access, usually volumes follow
    the C convention of starting index==0 (X,Y and Z). This function should
    not be used as it goes against the vector library philosophy unless you
    explicitly want to access directly to any value in the volume
    without taking into account its logical position
    \\ Ex: DIRECT_VOL_ELEM(V,0,0,0)=1;
    \\ Ex: val=DIRECT_VOL_ELEM(V,0,0,0);*/
#define DIRECT_VOL_ELEM(V,k,i,j) (V).__m[(k)*XYSIZE(V)+(i)*XSIZE(V)+(j)]

/** A short alias for previous function.
    To avoid writing so much*/
#define dVkij(V,k,i,j) DIRECT_VOL_ELEM(V,k,i,j)

/** Array access.
    This macro gives you access to the array (T *).
    \\ Ex: cout << "This is an int *" << VOL_ARRAY(V) << endl; */
#define VOL_ARRAY(V) MULTIDIM_ARRAY(V)
//@}
//@}

/// Template class for Xmipp volumes
#include "MultidimCommon.hh"
template <class T> class matrix3D {
#include "Src/MultidimBasic.hh"
/**@name Common functions to all multidimensional arrays
   A set of methods are always the same for any multidimensional array.
   Have a look on the more detailed structure. */
public:
//@{
//@Include: Src/MultidimBasic3.hh
//@}
   // This file contains several definitions for statistics and arithmetic
   // operations. To use it we have redirected the internal type maT
   // (multidimensional array<T>) to matrix3D<T>. These definitions are
   // outside because in this way we can reuse the module for other
   // libraries

/* Structure =============================================================== */
// Although the structure is defined as public it should not be used by
// the library user, there are functions enough to handle everything. This is
// done so because C++ does not allow friend classes when the friend class
// is a template.
public:
   int        zdim,ydim,xdim;    // dimensions of array [0...zdim-1]
                                 //                     [0...ydim-1]
                                 //                     [0...xdim-1]
   int        xydim;             // xydim = xdim * ydim
                                 // This is so to make element access faster
   int        zinit,yinit,xinit; // indexes of array  [zinit...zinit+zdim-1]
                                 //                   [yinit...yinit+ydim-1]
                                 //                   [xinit...xinit+xdim-1]

/* Procedures ============================================================== */
public:
   /* Constructors/Destructor ---------------------------------------------- */
   /**@name Constructors*/
   //@{
   /** Empty constructor.
       The empty constructor creates a volume with no memory associated,
       origin=0, size=0, no statistics, ...
       \\Ex: matrix3D<double> V1;*/
   matrix3D() {core_init(); init_shape(); __spcdim=3;}
   
   /** Dimension constructor.
       The dimension constructor creates a volume with memory associated
       (but not assigned to anything, could be full of garbage)
       origin=0, size=the given one, no statistics, ...
       Be careful that first number is the Z dimension (number of
       slices), then the Y dimension (number of rows),
       and at the end the X dimension (number of columns).
       \\Ex: matrix3D<double> V1(3,6,3);*/
   matrix3D(int Zdim, int Ydim, int Xdim)
         {core_init(); init_shape(); resize(Zdim, Ydim, Xdim); __spcdim=3;}

   /** Copy constructor.
       The created volume is a perfect copy of the input volume but
       with a different memory assignment.
       \\Ex: matrix3D<double> V2(V1); */
   matrix3D(const VT &V) {core_init(); init_shape(); *this=V;}

   // Destructor
   ~matrix3D() {core_deallocate();}
   //@}
   
   /* Initialisation ------------------------------------------------------- */
   /**@name Initialisation*/
   //@{
   /** Zero initialisation with a new dimension.
       Be careful to the size order (Zdim, Ydim, Xdim).
       \\Ex: v1.init_zeros(6,3);*/         
   void init_zeros(int Zdim, int Ydim, int Xdim)
         {resize(Zdim,Ydim,Xdim); init_constant((T)0);}
   //@}

   /* Memory related ------------------------------------------------------- */
   /**@name Size and shape
      The shape of a volume is defined by its origin and its size.
      The size is clear, and the origin
      is the logical position of the first real position of the array. For
      instance, if we have a matrix of dimension (2,5,3)=(Zdim,Ydim, Xdim) and
      origin (0,-2,-1), this means
      that the array is representing the logical positions
      \begin{verbatim}
      Slice 0
      [(0,-2,-1) (0,-2,0) (0,-2,1)
       (0,-1,-1) (0,-1,0) (0,-1,1)
       (0, 0,-1) (0, 0,0) (0, 0,1)
       (0, 1,-1) (0, 1,0) (0, 1,1)
       (0, 2,-1) (0, 2,0) (0, 2,1)]
       
      Slice 1
      [(1,-2,-1) (1,-2,0) (1,-2,1)
       (1,-1,-1) (1,-1,0) (1,-1,1)
       (1, 0,-1) (1, 0,0) (1, 0,1)
       (1, 1,-1) (1, 1,0) (1, 1,1)
       (1, 2,-1) (1, 2,0) (1, 2,1)]
      \end{verbatim}
      we could access to any of these positions (Ex: v(0,-2,1)=3;) and actually
      any try to access to a position related to 5 (Ex: v(0,4,1)=3;), although
      it physically exists, is not logically correct and hence it will
      throw an exception. The startingX and finishingX positions for this
      sample vector are -1 and 1 respectively, for Y are -2 and 2 and for Z
      are 0 and 1.
      The "for" iterations through the matrix should include these two
      values if you want to cover the whole matrix.

      \begin{verbatim}
         for (int k=STARTINGZ(V); k<=FINISHINGZ(V); k++)
            for (int i=STARTINGY(V); i<=FINISHINGY(V); i++)
               for (int j=STARTINGX(V); j<=FINISHINGX(V); j++)
                  VOL_ELEM(V,k,i,j) += 1;
      \end{verbatim}*/
      
   //@{
   /** Init shape.
       ydim,xdim=0, startingy,startingx=0.*/
   void init_shape()
      {xinit=yinit=zinit=0; xdim=ydim=zdim=0;}

   /** Copy shape.
       Copy shape variables from a pattern AND THE ARRAY IS RESIZED */
   template <class T1>
      void copy_shape(const maT1 &v)
         {if (XSIZE(*this)!=XSIZE(v) || YSIZE(*this)!=YSIZE(v) ||
              ZSIZE(*this)!=ZSIZE(v) )
              resize(ZSIZE(v), YSIZE(v), XSIZE(v));
          STARTINGX(*this)=STARTINGX(v);
          STARTINGY(*this)=STARTINGY(v);
          STARTINGZ(*this)=STARTINGZ(v);}

   /** Resize to a given size.
       This function resize the actual array to the given size. The origin
       is not modified. If the actual array is larger than the pattern
       then the values outside the new size are lost, if it is smaller
       then 0's are added. An exception is thrown if there is no memory.
       \\Ex: V1.resize(3,3,2);
   */
   void resize(int Zdim, int Ydim, int Xdim) {
      if (Xdim==XSIZE(*this) && Ydim==YSIZE(*this) &&
          Zdim==ZSIZE(*this)) return;
      if (Xdim<=0 || Ydim<=0 || Zdim<=0) {clear(); return;}

      // Ask for memory   
      T* new_m=new T [Zdim*Ydim*Xdim];
      if (new_m==NULL) REPORT_ERROR(1001,"Resize: no memory left");

      // Copy needed elements, fill with 0 if necessary
      int YXdim=Ydim*Xdim;
      for (int k=0; k<Zdim; k++)
         for (int i=0; i<Ydim; i++)
            for (int j=0; j<Xdim; j++) {
               T val;
               if      (k>=ZSIZE(*this)) val=0;
               else if (i>=YSIZE(*this)) val=0;
               else if (j>=XSIZE(*this)) val=0;
               else                      val=DIRECT_VOL_ELEM(*this,k,i,j);
               new_m[k*YXdim+Xdim*i+j]=val;
            }

      // deallocate old vector
      core_deallocate();

      // assign *this vector to the newly created
      MULTIDIM_ARRAY(*this)=new_m;
      XSIZE(*this)=Xdim;
      YSIZE(*this)=Ydim;
      ZSIZE(*this)=Zdim;
      XYSIZE(*this)=YXdim;
      __dim=Zdim*YXdim;
   }

   /** Produce an array suitable for working with Numerical Recipes.
       This function must be used only as
       a preparation for routines which need that the first physical
       index is 1 and not 0 as it usually is in C. New memory
       is needed to hold the new double pointer array.
       Click here to see an \URL[example]{../Extra_Docs/examples.html#LU} */
   T *** adapt_for_numerical_recipes() const
      {T ***m=NULL;
      ask_Tvolume(m, 1, ZSIZE(*this), 1, YSIZE(*this), 1, XSIZE(*this));
      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(*this)
         m[k+1][i+1][j+1]=DIRECT_VOL_ELEM(*this,k,i,j);
      return m;}

   /** Kill an array produced for numerical recipes.
       Nothing needs to be done in fact. */
   void kill_adaptation_for_numerical_recipes(T ***m) const
      {free_Tvolume(m, 1, ZSIZE(*this), 1, YSIZE(*this), 1, XSIZE(*this));}

   /** Intersects.
       TRUE if this array intersects with the box defined by the
       arguments (x0 is the starting X).*/
   bool intersects(double x0, double y0, double z0, double xdim, double ydim,
      double zdim) const;

   /** Outside.
       TRUE if the logical index given is outside the definition region
       of this array. */
   bool outside(int k, int i, int j) const;

   /** isBorder.
       TRUE if the logical index given belong to the border of the matrix */
   bool isBorder(int k,int i,int j);

   /** Set logical origin in Xmipp fashion.
       This function adjust the starting points in the volume such that
       the center of the volume is defined in the Xmipp fashion.
       \\Ex: V1.set_Xmipp_origin();
       @see Conventions */
   void set_Xmipp_origin()
      {zinit=FIRST_XMIPP_INDEX(zdim);
       yinit=FIRST_XMIPP_INDEX(ydim);
       xinit=FIRST_XMIPP_INDEX(xdim);}

   /** Move origin to.
       This function adjust logical indexes such that the Xmipp origin
       of the array moves to the specified position. For instance, an array
       whose x indexes go from -1 to 1, if we move the origin to 4, then
       the x indexes go from 3 to 5. This is very useful for convolution
       operations where you only need to move the logical starting of the
       array. 
       See \Ref{FIRST_XMIPP_INDEX} */
   void move_origin_to(int k, int i, int j) {
      zinit=k+FIRST_XMIPP_INDEX(zdim);
      yinit=i+FIRST_XMIPP_INDEX(ydim);
      xinit=j+FIRST_XMIPP_INDEX(xdim);}

   /** Sets the Z origin.
       The logical position of the first physical Z position is set with
       this function. By default the origin is 0 that is the standard
       convention in C.
       \\Ex: V.startingZ()=0; */
   int& startingZ() {return zinit;}

   /** Another function for setting the Z origin. */
   void set_startingZ(int _zinit) {zinit=_zinit;}

   /** Returns the first valid logical Z index.
       \\Ex: int orgZ=V.startingZ();*/
   int  startingZ() const  {return zinit;}

   /** Returns the last valid logical Z index.
       \\Ex: int finZ=V.finishingZ();*/
   int  finishingZ() const {return zinit+zdim-1;}

   /** Sets the Y origin.
       The logical position of the first physical Y position is set with
       this function. By default the origin is 0 that is the standard
       convention in C.
       \\Ex: V.startingY=-2; */
   int& startingY() {return yinit;}

   /** Another function for setting the Y origin. */
   void set_startingY(int _yinit) {yinit=_yinit;}

   /** Returns the first valid logical Y index.
       \\Ex: int orgY=V.startingY();*/
   int  startingY() const {return yinit;}

   /** Returns the last valid logical Y index.
       \\Ex: int finY=V.finishingY();*/
   int  finishingY() const {return yinit+ydim-1;}

   /** Sets the X origin.
       The logical position of the first physical X position is set with
       this function. By default the origin is 0 that is the standard
       convention in C.
       \\Ex: V.startingX=-1; */
   int& startingX() {return xinit;}

   /** Another function for setting the X origin. */
   void set_startingX(int _xinit) {xinit=_xinit;}

   /** Returns the first valid logical X index.
       \\Ex: int orgX=V.startingX();*/
   int  startingX() const {return xinit;}

   /** Returns the last valid logical X index.
       \\Ex: int finX=V.finishingX();*/
   int  finishingX() const {return xinit+xdim-1;}

   /** Returns the volume dimension.
       Pay attention to the dimension order (Z,Y,X).
       \\ Ex: V.get_dim(Zdim,Ydim,Xdim);*/
   void get_dim(int &Ydim, int &Xdim, int &Zdim) const
        {Xdim=xdim; Ydim=ydim; Zdim=zdim;}

   /** Returns Z dimension.
       \\Ex: int Zdim=V.SliNo();*/
   int  SliNo() const {return zdim;}

   /** Returns Y dimension.
       \\Ex: int Ydim=V.RowNo();*/
   int  RowNo() const {return ydim;}

   /** Returns X dimension.
       \\Ex: int Xdim=V.ColNo();*/
   int  ColNo() const {return xdim;}
   
   /** Same shape.
       Returns true if this object has got the same shape (origin and
       size) than the argument*/
   bool same_shape(const VT &op) const
      {return SAME_SHAPE3D(*this,op);}
   //@}
      
   /* Information Extraction ----------------------------------------------- */
   /**@name Memory access
      This functions allows you to access to the matrix elements.*/
   //@{
   /** Volume element access via index.
       Returns the value of a matrix logical position. In our example we could
       access from v(0,-2,-1) to v(1,2,1). The elements can be used either
       by value or by reference. An exception is thrown if the index is
       outside the logical range. Be careful that the argument order is
       (Z,Y,X).
       \\ Ex: V(0,-2,1)=1;
       \\ Ex: val=V(0,-2,1);*/
   T&   operator () (int k, int i,int j) const {
        if (k<zinit || k>=zinit+zdim)
           REPORT_ERROR(1203,"Matrix3D::operator (): Matrix3D subscript (k) out of range");
        if (i<yinit || i>=yinit+ydim)
           REPORT_ERROR(1203,"Matrix3D::operator (): Matrix3D subscript (i) out of range");
        if (j<xinit || j>=xinit+xdim)
           REPORT_ERROR(1203,"Matrix3D::operator (): Matrix3D subscript (j) out of range");
        return VOL_ELEM(*this, k,i,j);}

   /** Get the voxel at (k,i,j) (logical access)*/
   T get_voxel(int k,int i,int j) const {return (*this)(k,i,j);}

   /** Set the voxel at (k,i,j) (logical access)*/
   void set_voxel(int k,int i,int j, T val) {(*this)(k,i,j)=val;}

   /** Volume element access via double vector.
       Returns the value of a matrix logical position, but this time the
       element position is determined by a R3 vector.
       The elements can be used either by value or by reference.
       An exception is thrown if the index is outside the logical range.
       Pay attention in the following example that we are accessing the
       same element as in the previous function but, now we have to give
       first the X position because we are building
       first a vector of the form (x,y,z).
       \\ Ex: V(vector_R3(1,-2,0))=1;
       \\ Ex: val=V(vector_R3(1,-2,0));*/
   T&   operator () (const matrix1D<double> &v) const
         {return VOL_ELEM((*this),ROUND(ZZ(v)),
            ROUND(YY(v)),ROUND(XX(v)));}

   /** Volume element access via integer vector.*/
   T&   operator () (const matrix1D<int> &v) const
         {return VOL_ELEM((*this),ZZ(v),YY(v),XX(v));}

   /** Interpolates the value of the 3D matrix M at the point (x,y,z).
       (x,y,z) are in logical coordinates. */
   T   interpolated_elem(double x, double y, double z, T outside_value=(T)0) {
       int x0 = FLOOR(x); double fx = x - x0; int x1=x0+1;
       int y0 = FLOOR(y); double fy = y - y0; int y1=y0+1;
       int z0 = FLOOR(z); double fz = z - z0; int z1=z0+1;

       T d000 = (outside(z0,y0,x0)) ? outside_value : VOL_ELEM(*this,z0,y0,x0);
       T d001 = (outside(z0,y0,x1)) ? outside_value : VOL_ELEM(*this,z0,y0,x1);
       T d010 = (outside(z0,y1,x0)) ? outside_value : VOL_ELEM(*this,z0,y1,x0);
       T d011 = (outside(z0,y1,x1)) ? outside_value : VOL_ELEM(*this,z0,y1,x1);
       T d100 = (outside(z1,y0,x0)) ? outside_value : VOL_ELEM(*this,z1,y0,x0);
       T d101 = (outside(z1,y0,x1)) ? outside_value : VOL_ELEM(*this,z1,y0,x1);
       T d110 = (outside(z1,y1,x0)) ? outside_value : VOL_ELEM(*this,z1,y1,x0);
       T d111 = (outside(z1,y1,x1)) ? outside_value : VOL_ELEM(*this,z1,y1,x1);

       double dx00 = LIN_INTERP(fx, (double) d000, (double) d001);
       double dx01 = LIN_INTERP(fx, (double) d100, (double) d101);
       double dx10 = LIN_INTERP(fx, (double) d010, (double) d011);
       double dx11 = LIN_INTERP(fx, (double) d110, (double) d111);
       double dxy0 = LIN_INTERP(fy, (double) dx00, (double) dx10);
       double dxy1 = LIN_INTERP(fy, (double) dx01, (double) dx11);

       return (T) LIN_INTERP(fz, dxy0, dxy1);
   }

   /** Interpolates the value of the 3D matrix M at the point (x,y,z) knowing
       that this image is a set of B-spline coefficients.
       (x,y,z) are in logical coordinates.*/
   T   interpolated_elem_as_Bspline(double x, double y, double z,
       int SplineDegree=3) {
       int SplineDegree_1=SplineDegree-1;

       // Logical to physical
       z-=STARTINGZ(*this);
       y-=STARTINGY(*this);
       x-=STARTINGX(*this);

       int lmax = XSIZE(*this);
       int mmax = YSIZE(*this);
       int nmax = ZSIZE(*this);
       int l1 = CLIP(CEIL(x - SplineDegree_1),0,XSIZE(*this)-1);
       int l2 = CLIP(l1 + SplineDegree,0,XSIZE(*this)-1);
       int m1 = CLIP(CEIL(y - SplineDegree_1),0,YSIZE(*this)-1);
       int m2 = CLIP(m1 + SplineDegree,0,YSIZE(*this)-1);
       int n1 = CLIP(CEIL(z - SplineDegree_1),0,ZSIZE(*this)-1);
       int n2 = CLIP(n1 + SplineDegree,0,ZSIZE(*this)-1);

       double zyxsum = 0.0;
       for (int n=n1; n<=n2; n++)
          if (n<nmax && n>=-1L) {
             int plane_n = YSIZE(*this)*XSIZE(*this)*n;
      	     double yxsum = 0.0;
	     for (int m=m1; m<=m2; m++)
                 if (m<mmax && m>-1L) {
        	    int row_m = plane_n+XSIZE(*this)*m;
      		    double xsum = 0.0;
        	    for (int l=l1; l<=l2; l++) {
        	        double xminusl = x-(double)l;
        	        double Coeff = (double) __m[row_m+l];
		        switch (SplineDegree) {
			   case 2: xsum += Coeff * Bspline02(xminusl); break; 
                	   case 3: xsum += Coeff * Bspline03(xminusl); break;
                	   case 4: xsum += Coeff * Bspline04(xminusl); break;
                	   case 5: xsum += Coeff * Bspline05(xminusl); break;
                	   case 6: xsum += Coeff * Bspline06(xminusl); break;
                	   case 7: xsum += Coeff * Bspline07(xminusl); break;
                	   case 8: xsum += Coeff * Bspline08(xminusl); break;
                	   case 9: xsum += Coeff * Bspline09(xminusl); break;
		        }
        	    }
        	    double yminusm=y-(double)m;
		    switch (SplineDegree) {
        	       case 2: yxsum += xsum * Bspline02(yminusm); break;
        	       case 3: yxsum += xsum * Bspline03(yminusm); break;
        	       case 4: yxsum += xsum * Bspline04(yminusm); break;
        	       case 5: yxsum += xsum * Bspline05(yminusm); break;
        	       case 6: yxsum += xsum * Bspline06(yminusm); break;
        	       case 7: yxsum += xsum * Bspline07(yminusm); break;
        	       case 8: yxsum += xsum * Bspline08(yminusm); break;
        	       case 9: yxsum += xsum * Bspline09(yminusm); break;
		    }
                 }
              double zminusn=z-(double)n;
	      switch (SplineDegree) {
                 case 2: zyxsum += yxsum * Bspline02(zminusn); break;
                 case 3: zyxsum += yxsum * Bspline03(zminusn); break;
                 case 4: zyxsum += yxsum * Bspline04(zminusn); break;
                 case 5: zyxsum += yxsum * Bspline05(zminusn); break;
                 case 6: zyxsum += yxsum * Bspline06(zminusn); break;
                 case 7: zyxsum += yxsum * Bspline07(zminusn); break;
                 case 8: zyxsum += yxsum * Bspline08(zminusn); break;
                 case 9: zyxsum += yxsum * Bspline09(zminusn); break;
              }
         }
       return (T)zyxsum;
   }

   /** Logical to physical index translation.
       This function returns the physical position of a logical one. See
       \URL[Conventions]{../../../Extra_Docs/Conventions.html}
       for more information about these two different accesses.
       \\ Ex: m.logical2physical(k_log,i_log,j_log,k_phys,i_phys,j_phys); */
   void logical2physical(int k_log, int i_log, int j_log,
       int &k_phys, int &i_phys, int &j_phys) const
       {k_phys=k_log-zinit; i_phys=i_log-yinit; j_phys=j_log-xinit;}
   
   /** Physical to logical index translation.
       This function returns the logical position of a physical one. See
       \URL[Conventions]{../../../Extra_Docs/Conventions.html}
       for more information about these two different accesses.
       \\ Ex: m.physical2logical(i_phys,j_phys,i_log,j_log); */
   void physical2logical(int k_phys, int i_phys, int j_phys,
       int &k_log, int &i_log, int &j_log) const
       {k_log=k_phys+zinit; i_log=i_phys+yinit; j_log=j_phys+xinit;}

   /** Get Slice.
       This function returns a slice (a matrix) corresponding to the choosen
       slice inside matrix, the numbering of the slices is also logical not
       physical. By default slices are taken perpendicular to the Z axis,
       but you can specify different axis ('X' and 'Y'). When cutting slices
       the following axes conventions are followed.
       \begin{verbatim}
       Cut along Z axis       Y(2D)=Y(3D)    X(2D)=X(3D)
       Cut along Y axis       Y(2D)=Z(3D)    X(2D)=X(3D)
       Cut along X axis       Y(2D)=Z(3D)    X(2D)=-Y(3D)
       \end{verbatim}
       \\Ex: matrix2D<double> m=V.slice(0);*/
   mT  getSlice(int i, char axis='Z') const
       {mT temp; getSlice(i,temp,axis); return temp;}
   
   /** Slice access for reading.
       This function returns a slice (a matrix) corresponding to the choosen
       slice inside matrix, the numbering of the slices is also logical not
       physical. This function differs from the previous one in that this one
       cuts and assign in a single step instead of in two steps, as in
       the previous example.
       \\Ex: V.slice(0,m);*/
   void getSlice(int k, mT &M, char axis='Z') const {
      if (xdim==0) {M.clear(); return;}   
      switch (axis) {
         case 'Z':
            if (k<zinit || k>=zinit+zdim)
               REPORT_ERROR(1203,"Slice: matrix3D subscript (k) out of range");

            k=k-zinit;
            M.resize(ydim,xdim);
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                  DIRECT_MAT_ELEM(M,i,j)=DIRECT_VOL_ELEM(*this,k,i,j);
            STARTINGX(M)=STARTINGX(*this);
            STARTINGY(M)=STARTINGY(*this);
            break;
         case 'Y':
            if (k<yinit || k>=yinit+ydim)
               REPORT_ERROR(1203,"Slice: matrix3D subscript (i) out of range");

            k=k-yinit;
            M.resize(zdim,xdim);
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                  DIRECT_MAT_ELEM(M,i,j)=DIRECT_VOL_ELEM(*this,i,k,j);
            STARTINGX(M)=STARTINGX(*this);
            STARTINGY(M)=STARTINGZ(*this);
            break;
         case 'X':
            if (k<xinit || k>=xinit+xdim)
               REPORT_ERROR(1203,"Slice: matrix3D subscript (j) out of range");

            k=k-xinit;
            M.resize(zdim,ydim);
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                  DIRECT_MAT_ELEM(M,i,j)=DIRECT_VOL_ELEM(*this,i,j,k);
            STARTINGX(M)=STARTINGY(*this);
            STARTINGY(M)=STARTINGZ(*this);
            break;
         default:
            REPORT_ERROR(1205,(string)"Slice: not supported axis "+axis);
      }
   }
   
   /** Slice access for writing.
       This function sets a matrix corresponding to the choosen
       slice inside volume, the numbering of the slices is also logical not
       physical.
       \\Ex: V.setSlice(1,(V.slice(0)));
       \\--> Copies slice 0 in slice 1 */
    void setSlice(int k, const mT &v) {
      if (xdim==0) return;
      if (k<zinit || k>=zinit+zdim)
         REPORT_ERROR(1203,"setSlice: matrix3D subscript (k) out of range");
      if (v.RowNo()!=ydim || v.ColNo()!=xdim)
         REPORT_ERROR(1202,"setSlice: matrix3D dimensions different from the matrix ones");

      k=k-zinit;
      FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(v)
            DIRECT_VOL_ELEM(*this,k,i,j)=DIRECT_MAT_ELEM(v,i,j);
   }
   //@}

   /* Other utilities ------------------------------------------------------ */
   /**@name Utilities*/
   //@{
   /** This function must take two arrays of the same size, and operate
      element by element according to the operation required. This is the
      function which really implements the operations. Simple calls to it
      perform much faster than calls to the corresponding operators.
      Although it is supposed to be a hidden function not useable by
      normal programmers.
      
      It must be implemented in every Matrix module, this is so because
      of the Matrix2D, for which the multiplication is not a component
      by component multiplication but an algebraic one.*/
   friend void array_by_array(const maT &op1, const maT &op2, maT &result,
      char operation) {
         if (!op1.same_shape(op2))
            REPORT_ERROR(1007,
               (string)"Array_by_array: different shapes ("+operation+")");
         if (operation=='x') operation='*';
         result.resize(op1);
         core_array_by_array(op1, op2, result, operation);
      }

   /** Reverse volume values over X axis.
       Maybe better with an example:

       \begin{verbatim}
         slice 0
         [01 02 03          [07 08 09
          04 05 06           04 05 06
          07 08 09]          01 02 03]
                     ----->
         slice 1
          [11 12 13          [17 18 19
           14 15 16           14 15 16
           17 18 19]          11 12 13]
      \end{verbatim}
       \\Ex: V2=V1.reverseX();*/
   VT reverseX() const {VT temp(*this); temp.self_reverseX(); return temp;}

   /** Reverse matrix values over X axis, keep in this object. */
   void self_reverseX() {
      for (int k=0; k<zdim; k++)
         for (int i=0; i<ydim; i++)
            for (int j=0; j<=(int)(xdim-1)/2; j++) {
               T aux;
               if (k==0 && i==0)
                  cout << "Changing " << j << " " << XSIZE(*this)-1-j << endl;
               SWAP(DIRECT_VOL_ELEM(*this,k,i,j),
                    DIRECT_VOL_ELEM(*this,k,i,XSIZE(*this)-1-j),aux);
            }
      STARTINGX(*this)=-FINISHINGX(*this);
   }
   
   /** Reverse matrix values over Y axis.
       Maybe better with an example:

       \begin{verbatim}
         slice 0
         [01 02 03          [03 02 01
          04 05 06           06 05 04
          07 08 09]          09 08 07]
                     ----->
         slice 1
         [11 12 13          [13 12 11
          14 15 16           16 15 14
          17 18 19]          19 18 17]
       \end{verbatim}
       \\Ex: V2=V1.reverseY();*/
   VT reverseY() const {VT temp(*this); temp.self_reverseY(); return temp;}

   /** Reverse matrix values over Y axis, keep in this object. */
   void self_reverseY(){
      for (int k=0; k<zdim; k++)
         for (int i=0; i<=(int)(ydim-1)/2; i++)
            for (int j=0; j<xdim; j++) {
               T aux;
               SWAP(DIRECT_VOL_ELEM(*this,k,i,j),
                    DIRECT_VOL_ELEM(*this,k,YSIZE(*this)-1-i,j),aux);
            }
      STARTINGY(*this)=-FINISHINGY(*this);
   }
   
   /** Reverse matrix values over Z axis.
       Maybe better with an example:

       \begin{verbatim}
         slice 0
         [01 02 03          [11 12 13
          04 05 06           14 15 16
          07 08 09]          17 18 19]
                     ----->
         slice 1
         [11 12 13          [01 02 03
          14 15 16           04 05 06
          17 18 19]          07 08 09]
       \end{verbatim}
       \\Ex: V2=V1.reverseZ();*/
   VT reverseZ() const {VT temp(*this); temp.self_reverseZ(); return temp;}

   /** Reverse matrix values over Z axis, keep in this object. */
   void self_reverseZ() {
      for (int k=0; k<=(int)(zdim-1)/2; k++)
         for (int i=0; i<ydim; i++)
            for (int j=0; j<xdim; j++) {
               T aux;
               SWAP(DIRECT_VOL_ELEM(*this,k,i,j),
                    DIRECT_VOL_ELEM(*this,ZSIZE(*this)-1-k,i,j),aux);
            }
      STARTINGZ(*this)=-FINISHINGZ(*this);
   }
   
   /** Put a window to volume.
       The volume is windowed within the two positions given to this function.
       Indexes always refer to logical indexes. If a position is outside the
       actual matrix range then the matrix is padded init_value until the
       new position is reached. In the following example suppose that m1
       is the following and that the origin is (-1,-1,-1).

       \begin{verbatim}
         slice 0
         [01 02 03          [
          04 05 06           04 05 06 0
          07 08 09]          07 08 09 0]
                     ----->
         slice 1
         [11 12 13          [
          14 15 16           14 15 16 0
          17 18 19]          17 18 19 0]
       \end{verbatim}
       \\Ex: V1.window(0,0,-1,1,1,2);
       */
   void window(int z0, int y0, int x0, int zF, int yF, int xF, T init_value=0) {
      VT result(zF-z0+1, yF-y0+1, xF-x0+1);
      result.zinit=z0;
      result.yinit=y0;
      result.xinit=x0;

      for (int k=z0; k<=zF; k++)
         for (int i=y0; i<=yF; i++)
            for (int j=x0; j<=xF; j++)
                if ((k>=zinit && k<=zinit+zdim-1) &&
                    (i>=yinit && i<=yinit+ydim-1) &&
                    (j>=xinit && j<=xinit+xdim-1))
                        VOL_ELEM(result,k,i,j)=VOL_ELEM(*this,k,i,j);
                else
                        VOL_ELEM(result,k,i,j)=init_value;
      *this=result;
   }
   //@}

   /* Geometrical transformations ------------------------------------------ */
   /**@name Geometrical Transformations
      In all geometrical transformations a periodic extension of the volume
      is supposed, ie, if a voxel goes out on the left, it is entering on
      the right, ...*/
   //@{
   /** Applies a geometrical transformation.
       Any geometrical transformation defined by the matrix A (double (4x4)!!
       ie, in homogeneous R3 coordinates) is applied to the volume V1.
       The result is stored in V2 (it cannot be the same as the input volume).
       An exception is thrown if the
       transformation matrix is not 4x4.
	   
	   Structure of the transformation matrix: It should have the following 
	   components
	   r11 r12 r13 x
	   r21 r22 r23 y
	   r31 r32 r33 z
	   0   0   0   1
	   where (x,y,z) is the translation desired, and Rij are the components of
	   the rotation matrix R. If you want to apply a scaling factor to the 
	   transformation, then multiply r11, r22 and r33 by it. 
	   
       The result volume is resized to the same dimensions as V1 if V2 is empty
       (0x0) at the beginning, if it is not, ie, if V2 has got some size
       then only those values in the volume are filled, this is very
       useful for resizing the volume, then you manually resize the output
       volume to the desired size and then call this routine.

       The relationship between the output coordinates and the input ones are
       \begin{verbatim}
           out = A * in
       (x,y,z) = A * (x',y',z')
       \end{verbatim}

       This function works independently from the logical indexing of each
       matrix, it sets the logical center and the physical center of the image
       and work with these 2 coordinate spaces. At the end the original logical
       indexing of each matrix is kept.

       The procedure followed goes from coordinates in the output volume
       to the ones in the input one, so the inverse of the A matrix is
       needed. There is a flag telling if the given matrix is already
       the inverse one or the normal one. If it is the normal one internally
       the matrix is inversed. If you are to do many "rotations" then
       some time is spent in inverting the matrix. Normally the
       matrix is the normal one.
       
       There is something else to tell about the geometrical tranformation.
       The value of the voxel in the output volume is computed via
       bilinear interpolation in the input volume. If any of the voxels
       participating in the interpolation falls outside the input volume,
       then automatically the corresponding output voxel is set to 0, unless
       that the wrap flag has been set to 1. In this case if the voxel
       falls out by the right hand then it is "wrapped" and the corresponding
       voxel in the left hand is used. The same is appliable to top-bottom.
       Usually wrap mode is off. Wrap mode is interesting for translations
       but not for rotations, for example.
       
       The inverse mode and wrapping mode should be taken by default by the
       routine, g++ seems to have problems with template functions outside
       a class with default parameters. So, I'm sorry, you will have to
       put them always. The usual combination is
       apply_geom(...,IS_NOT_INV,DONT_WRAP). Although you can also use the
       constants IS_INV, or WRAP.
       \\Ex: matrix2D<double> A(4,4); A.init_identity; apply_geom(V2,A,V1);*/
   friend void apply_geom<>(VT &V2, matrix2D<double> A,
            const VT &V1, bool inv, bool wrap);

   /** Apply geom with B-spline interpolation. */
   friend void apply_geom_Bspline<>(VT &V2, matrix2D<double> A,
       const VT &V1, int Splinedegree, bool inv, bool wrap, 
       T outside);

   /** Self apply geom.
       As apply geometry, but the result is kept in this object */
   void self_apply_geom(matrix2D<double> A, bool inv, bool wrap)
      {VT aux; apply_geom(aux, A,*this,inv,wrap); *this=aux;}

   /** Self apply geom Bspline. */
   void self_apply_geom_Bspline(matrix2D<double> A, int SplineDegree,
      bool inv, bool wrap, T outside=0)
      {VT aux; apply_geom_Bspline(aux, A,*this,SplineDegree,inv,wrap,outside);
          *this=aux;}

   /** Rotate a volume around system axis.
       The rotation angle is in degrees, and the rotational axis is
       either 'X', 'Y' or 'Z'. An exception is thrown if the axis given
       is not one of these.
       \\Ex: V2=V1.rotate(60);*/
   void rotate(double ang, char axis, VT &result, bool wrap=DONT_WRAP) const
      {matrix2D<double> temp=rot3D_matrix(ang, axis);
       apply_geom(result,temp,*this,IS_NOT_INV,wrap);}

   /** Rotate a volume arounf system axis (BSpline). */
   void rotate_Bspline(int Splinedegree, double ang, char axis, VT &result,
      bool wrap=DONT_WRAP, T outside=0) const
      {matrix2D<double> temp=rot3D_matrix(ang, axis);
       apply_geom_Bspline(result,temp,*this,IS_NOT_INV,wrap,outside);}

   /** Rotate a volume around system axis, return result.*/
   VT rotate(double ang, char axis, bool wrap=DONT_WRAP) const
      {VT aux; rotate(ang, axis, aux, wrap); return aux;}

   /** Rotate a volume around system axis, return result (Bspline).*/
   VT rotate_Bspline(int Splinedegree, double ang, char axis,
      bool wrap=DONT_WRAP) const
      {VT aux; rotate_Bspline(Splinedegree, ang, axis, aux, wrap); return aux;}

   /** Rotate a volume around system axis, keep in this object. */
   void self_rotate(double ang, char axis, bool wrap=DONT_WRAP)
      {VT aux; rotate(ang, axis, aux, wrap); *this=aux;}

   /** Rotate a volume around system axis, keep in this object (Bspline). */
   void self_rotate_Bspline(int Splinedegree, double ang, char axis,
      bool wrap=DONT_WRAP)
      {VT aux; rotate_Bspline(Splinedegree, ang, axis, aux, wrap); *this=aux;}

   /** Rotate a volume around any axis.
       The rotation angle is in degrees, and the rotational axis is
       given as a R3 vector. An exception is thrown if the axis is not a
       R3 vector. The axis needs not to be unitary.
       \\Ex: V2=V1.rotate(60,vector_R3(1,1,1));*/
   void rotate(double ang, const matrix1D<double> &axis, VT &result,
      bool wrap=DONT_WRAP) const
      {matrix2D<double> temp=rot3D_matrix(ang,axis);
       apply_geom(result,temp,*this,IS_NOT_INV,wrap);}

   /** Rotate a volume around any axis (Bspline).*/
   void rotate_Bspline(int Splinedegree, double ang,
      const matrix1D<double> &axis, VT &result, 
      bool wrap=DONT_WRAP, T outside=0) const
      {matrix2D<double> temp=rot3D_matrix(ang,axis);
       apply_geom_Bspline(result,temp,*this,Splinedegree,IS_NOT_INV,wrap,outside);}

   /** Rotate a volume around any axis, return result. */
   VT rotate(double ang, const matrix1D<double> v, bool wrap=DONT_WRAP) const
      {VT aux; rotate(ang, v, aux, wrap); return aux;}

   /** Rotate a volume around any axis, return result (Bspline). */
   VT rotate_Bspline(int Splinedegree, double ang, const matrix1D<double> v,
      bool wrap=DONT_WRAP) const
      {VT aux; rotate_Bspline(Splinedegree, ang, v, aux, wrap); return aux;}

   /** Rotate a volume around any axis, keep in this object. */
   void self_rotate(double ang, const matrix1D<double> &v, bool wrap=DONT_WRAP)
      {VT aux; rotate(ang, v, aux, wrap); *this=aux;}

   /** Rotate a volume around any axis, keep in this object (Bspline). */
   void self_rotate_Bspline(int Splinedegree, double ang,
      const matrix1D<double> &v, bool wrap=DONT_WRAP)
      {VT aux; rotate_Bspline(Splinedegree, ang, v, aux, wrap); *this=aux;}

   /** Translate a volume.
       The shift is given as a R3 vector (shift_X, shift_Y, shift_Z);
       An exception is thrown if the displacement is not a R3 vector.
       \\Ex: V2=V1.translate(vector_R3(0,0,2));
       \\--> Displacement of 2 pixels down */
   void translate(const matrix1D<double> &v, VT &result, bool wrap=WRAP) const
      {matrix2D<double> temp=translation3D_matrix(v);
       apply_geom(result,temp,*this,IS_NOT_INV,wrap);}

   /** Translate a volume (Bspline).*/
   void translate_Bspline(int Splinedegree, const matrix1D<double> &v,
      VT &result, bool wrap=WRAP) const
      {matrix2D<double> temp=translation3D_matrix(v);
       apply_geom(result,temp,*this,IS_NOT_INV,wrap);}

   /** Translate a volume, return result.*/
   VT translate(const matrix1D<double> &v, bool wrap=WRAP) const
      {VT aux; translate(v, aux, wrap); return aux;}

   /** Translate a volume, return result (Bspline).*/
   VT translate_Bspline(int Splinedegree, const matrix1D<double> &v,
      bool wrap=WRAP) const
      {VT aux; translate_Bspline(Splinedegree, v, aux, wrap); return aux;}

   /** Translate a volume, keep in this object.*/
   void self_translate(const matrix1D<double> &v, bool wrap=WRAP)
      {VT aux; translate(v, aux, wrap); *this=aux;}

   /** Translate a volume, keep in this object (Bspline).*/
   void self_translate_Bspline(int Splinedegree, const matrix1D<double> &v,
      bool wrap=WRAP)
      {VT aux; translate_Bspline(Splinedegree, v, aux, wrap); *this=aux;}

   /** Translate center of mass to center.
       If the input has very high values, sometimes it is better to
       rescale it to be between 0 and 1. */
   void self_translate_center_of_mass_to_center(bool wrap=WRAP) {
      set_Xmipp_origin();
      matrix1D<double> center;
      center_of_mass(center);
      center*=-1;
      self_translate(center,wrap);
   }

   /** Translate center of mass to center (Bspline).*/
   void self_translate_center_of_mass_to_center_Bspline(
      int Splinedegree, bool wrap=WRAP) {
      set_Xmipp_origin();
      matrix1D<double> center;
      center_of_mass(center);
      center*=-1;
      self_translate_Bspline(Splinedegree,center,wrap);
   }

   /** Scales to a new size.
       The volume is scaled (resampled) to fill a new size. It is not the
       same as "window" in this same class. The size can be larger or 
       smaller than the actual one.
       \\Ex: V2=V1.scale_to_size(128,128,128);*/
   void scale_to_size(int Zdim, int Ydim, int Xdim, VT &result) const
      {matrix2D<double> temp(4,4);
       temp.init_identity();
       DIRECT_MAT_ELEM(temp,0,0)=(double)Xdim/(double)xdim;
       DIRECT_MAT_ELEM(temp,1,1)=(double)Ydim/(double)ydim;
       DIRECT_MAT_ELEM(temp,2,2)=(double)Zdim/(double)zdim;
       result.resize(Zdim,Ydim,Xdim);
       apply_geom(result,temp,*this,IS_NOT_INV,WRAP);}

   /** Scales to a new size (Bspline).*/
   void scale_to_size_Bspline(int Splinedegree, int Zdim, int Ydim, int Xdim,
      VT &result) const
      {matrix2D<double> temp(4,4);
       temp.init_identity();
       DIRECT_MAT_ELEM(temp,0,0)=(double)Xdim/(double)xdim;
       DIRECT_MAT_ELEM(temp,1,1)=(double)Ydim/(double)ydim;
       DIRECT_MAT_ELEM(temp,2,2)=(double)Zdim/(double)zdim;
       result.resize(Zdim,Ydim,Xdim);
       apply_geom_Bspline(result,temp,*this,Splinedegree,IS_NOT_INV,WRAP);}

   /** Scales to a new size, return result*/
   VT scale_to_size(int Zdim, int Ydim, int Xdim) const
      {VT aux; scale_to_size(Zdim,Ydim,Xdim,aux); return aux;}

   /** Scales to a new size, return result (Bspline).*/
   VT scale_to_size_Bspline(int Splinedegree, int Zdim, int Ydim,
      int Xdim) const
      {VT aux; scale_to_size_Bspline(Splinedegree, Zdim,Ydim,Xdim,aux);
       return aux;}

   /** Scales to a new size., keep in this object*/
   void self_scale_to_size(int Zdim, int Ydim, int Xdim)
      {VT aux; scale_to_size(Zdim,Ydim,Xdim,aux); *this=aux;}

   /** Scales to a new size., keep in this object (Bspline)*/
   void self_scale_to_size_Bspline(int Splinedegree,
      int Zdim, int Ydim, int Xdim)
      {VT aux; scale_to_size_Bspline(Splinedegree,Zdim,Ydim,Xdim,aux);
       *this=aux;}

   /** Reduce the image by 2 using a BSpline pyramid. */
   void pyramid_reduce(matrix3D<double> &result, int levels=1) const {
      matrix3D<double> aux, aux2;
      produce_spline_coeffs(aux,3);
      for (int i=0; i<levels; i++) {
         aux.reduce_Bspline(aux2,3);
         aux=aux2;
      }
      aux2.produce_image_from_spline_coeffs(result,3);
   }

   /** Expand the image by 2 using a BSpline pyramid. */
   void pyramid_expand(matrix3D<double> &result, int levels=1) const {
      matrix3D<double> aux, aux2;
      produce_spline_coeffs(aux,3);
      for (int i=0; i<levels; i++) {
         aux.expand_Bspline(aux2,3);
         aux=aux2;
      }
      aux2.produce_image_from_spline_coeffs(result,3);
   }

   /** Produce spline coefficients.*/
   #ifndef DBL_EPSILON
      #define DBL_EPSILON 1e-50
   #endif
   void produce_spline_coeffs(matrix3D<double> &coeffs, int SplineDegree=3)
      const {
      coeffs.init_zeros(ZSIZE(*this),YSIZE(*this),XSIZE(*this));
      STARTINGX(coeffs)=STARTINGX(*this);
      STARTINGY(coeffs)=STARTINGY(*this);
      STARTINGZ(coeffs)=STARTINGZ(*this);
      int Status;
      matrix3D<double> aux;
      type_cast(*this,aux);
      ChangeBasisVolume(MULTIDIM_ARRAY(aux),MULTIDIM_ARRAY(coeffs),
          XSIZE(*this), YSIZE(*this), ZSIZE(*this),
          CardinalSpline, BasicSpline, SplineDegree,
          MirrorOffBounds, DBL_EPSILON, &Status);
      if (Status)
         REPORT_ERROR(1,"matrix3D::produce_spline_coeffs: Error");
   }

   /** Produce image from B-spline coefficients. */
   void produce_image_from_spline_coeffs(
      matrix3D<double> &img, int SplineDegree=3) const {
      img.init_zeros(ZSIZE(*this),YSIZE(*this),XSIZE(*this));
      STARTINGX(img)=STARTINGX(*this);
      STARTINGY(img)=STARTINGY(*this);
      STARTINGZ(img)=STARTINGZ(*this);
      int Status;
      matrix3D<double> aux;
      type_cast(*this,aux);
      ChangeBasisVolume(MULTIDIM_ARRAY(aux),MULTIDIM_ARRAY(img),
          XSIZE(*this), YSIZE(*this), ZSIZE(*this),
          BasicSpline, CardinalSpline, SplineDegree,
          MirrorOnBounds, DBL_EPSILON, &Status);
      if (Status)
         REPORT_ERROR(1,"matrix3D::produce_spline_img: Error");
   }
   #undef DBL_EPSILON

   /** Expand a set of B-spline coefficients.
       Knowing that this matrix is a set of B-spline coefficients,
       produce the expanded set of B-spline coefficients using the
       two-scale relationship. */
   void expand_Bspline(matrix3D<double> &expanded, int SplineDegree=3) const {
      double   g[200];     /* Coefficients of the reduce filter */
      long     ng;         /* Number of coefficients of the reduce filter */
      double   h[200];     /* Coefficients of the expansion filter */
      long     nh;         /* Number of coefficients of the expansion filter */
      short    IsCentered; /* Equal TRUE if the filter is a centered spline, FALSE otherwise */

      // Get the filter
      if (GetPyramidFilter( "Centered Spline", SplineDegree,
            g, &ng, h, &nh, &IsCentered))
	 REPORT_ERROR(1,"Unable to load the filter coefficients");

      matrix3D<double> aux;
      type_cast(*this,aux);
      expanded.resize(2*ZSIZE(aux),2*YSIZE(aux),2*XSIZE(aux));
      Expand_3D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux), ZSIZE(aux),
         MULTIDIM_ARRAY(expanded), h, nh, IsCentered);
   }

   /** Reduce a set of B-spline coefficients.
       Knowing that this matrix is a set of B-spline coefficients,
       produce the reduced set of B-spline coefficients using the
       two-scale relationship. */
   void reduce_Bspline(matrix3D<double> &reduced, int SplineDegree=3) const {
      double   g[200];     /* Coefficients of the reduce filter */
      long     ng;         /* Number of coefficients of the reduce filter */
      double   h[200];     /* Coefficients of the expansion filter */
      long     nh;         /* Number of coefficients of the expansion filter */
      short    IsCentered; /* Equal TRUE if the filter is a centered spline, FALSE otherwise */

      // Get the filter
      if (GetPyramidFilter( "Centered Spline", SplineDegree,
            g, &ng, h, &nh, &IsCentered))
	 REPORT_ERROR(1,"Unable to load the filter coefficients");

      matrix3D<double> aux;
      type_cast(*this,aux);
      if      (XSIZE(aux)%2!=0 && YSIZE(aux)%2!=0 && ZSIZE(aux)%2!=0)
          aux.resize(ZSIZE(aux-1),YSIZE(aux)-1,XSIZE(aux)-1);
      else if (XSIZE(aux)%2!=0 && YSIZE(aux)%2!=0 && ZSIZE(aux)%2==0)
          aux.resize(ZSIZE(aux),YSIZE(aux)-1,XSIZE(aux)-1);
      else if (XSIZE(aux)%2!=0 && YSIZE(aux)%2==0 && ZSIZE(aux)%2!=0)
          aux.resize(ZSIZE(aux)-1,YSIZE(aux),XSIZE(aux)-1);
      else if (XSIZE(aux)%2!=0 && YSIZE(aux)%2==0 && ZSIZE(aux)%2==0)
          aux.resize(ZSIZE(aux),YSIZE(aux),XSIZE(aux)-1);
      else if (XSIZE(aux)%2==0 && YSIZE(aux)%2!=0 && ZSIZE(aux)%2!=0)
          aux.resize(ZSIZE(aux)-1,YSIZE(aux)-1,XSIZE(aux));
      else if (XSIZE(aux)%2==0 && YSIZE(aux)%2!=0 && ZSIZE(aux)%2==0)
          aux.resize(ZSIZE(aux),YSIZE(aux)-1,XSIZE(aux));
      else if (XSIZE(aux)%2==0 && YSIZE(aux)%2==0 && ZSIZE(aux)%2!=0)
          aux.resize(ZSIZE(aux)-1,YSIZE(aux),XSIZE(aux));
      reduced.resize(ZSIZE(aux)/2,YSIZE(aux)/2,XSIZE(aux)/2);
      Reduce_3D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux), ZSIZE(aux),
         MULTIDIM_ARRAY(reduced), g, ng, IsCentered);
   }

   /** Maximum element.
       This function returns the index of the maximum element of an array.
       array(k,i,j). Returns -1 if the array is empty*/
   void max_index(int &kmax, int &imax, int &jmax) const {
      if (XSIZE(*this)==0) {kmax=imax=jmax=-1; return;}
      kmax=imax=jmax=0;
      T   max=VOL_ELEM(*this,kmax,imax,jmax);
      FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
         if (VOL_ELEM(*this,k,i,j)>max)
            {max=VOL_ELEM(*this,k,i,j); kmax=k; imax=i; jmax=j;}
   }

   /** Minimum element.
       This function returns the index of the minimum element of an array.
       array(k,i,j). Returns -1 if the array is empty*/
   void min_index(int &kmin, int &imin, int &jmin) const {
      if (XSIZE(*this)==0) {kmin=imin=jmin=-1; return;}
      kmin=imin=jmin=0;
      T   min=VOL_ELEM(*this,kmin,imin,jmin);
      FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
         if (VOL_ELEM(*this,k,i,j)>min)
            {min=VOL_ELEM(*this,k,i,j); kmin=k; imin=i; jmin=j;}
   }
   //@}

   /* Iterators ------------------------------------------------------------ */
   /**@name Iterators*/
   //@{
   /** Apply the same scalar function to all slices.
       This function must take a matrix and return a single value, a column
       vector with these values is returned.
       \\Ex:T matrix_sum(mT &m) {return m.sum();};
            v1=V.for_all_slices(&matrix_sum);*/
   vT for_all_slices (T (*f)(mT&)) const {
      vT temp;

      if (ZSIZE(*this)==0) {temp.clear(); return temp;};
      temp.resize(ZSIZE(*this));
      STARTINGX(temp)=STARTINGZ(*this);

      for (int k=STARTINGZ(*this); k<=FINISHINGZ(*this); k++) {
         mT aux;
         getSlice(k,aux);
         VEC_ELEM(temp,k)=(*f)(aux);
      }

      return temp;
   }

   /** Apply the same matricial function to all slices.
       This function must take a matrix and return another matrix (of the
       same shape as the input one), a new
       volume with these transformed matrices is returned.
       \\Ex:mT matrix_norm(mT &m) {return m/m.sum();};
            V2=V.for_all_slices(&matrix_norm);*/
   VT for_all_slices (mT (*f)(mT&)) const {
      VT temp;

      if (xdim==0) {temp.clear(); return temp;};
      temp.copy_shape(*this);

      for (int k=STARTINGZ(*this); k<=FINISHINGZ(*this); k++) {
         mT aux;
         getSlice(k,aux);
         aux=(*f)(aux);
         temp.setSlice(k,aux);
      }
      return temp;
   }
   //@}

};

#include "Src/MultidimFriends_implementation.hh"

template <>
void core_array_by_scalar< complex<double> >(const maTC &op1,
   const complex<double> &op2, maTC &result, char operation);

template <>
void core_scalar_by_array< complex<double> >(const complex<double> &op1,
   const maTC &op2, maTC &result, char operation);

template <>
void core_array_by_array< complex<double> >(const maTC &op1, const maTC &op2,
   maTC &result, char operation);

/**@name Related functions
   These functions are not methods of matrix1D */
//@{

/* Other useful functions -------------------------------------------------- */
/**@name Miscellaneous*/
//@{
/** Reduce both volumes to a common size.
    Search the range of logical indexes for which both volumes have got
    valid values, and cut both to that size, the corresponding
    origin is automatically computed.
    \\Ex: matrix3D<double> V1(4,5,3);
    \\V1.startingX()=-2; V1.startingY()=-2; V1.startingZ()=-2;
    \\matrix3D<double> V2(4,2,3);
    \\V2.startingX()=0; V2.startingY()=0; V2.startingZ()=0;
    \\cut_to_common_size(V1,V2);
    \\--> V1 and V2 range from (0,0,0)=(z,y,x) to (1,1,0) */
template <class T>
   void cut_to_common_size(VT &V1, VT &V2) {
      int z0=MAX(STARTINGZ(V1) ,STARTINGZ(V2));
      int zF=MIN(FINISHINGZ(V1),FINISHINGZ(V2));
      int y0=MAX(STARTINGY(V1) ,STARTINGY(V2));
      int yF=MIN(FINISHINGY(V1),FINISHINGY(V2));
      int x0=MAX(STARTINGX(V1) ,STARTINGX(V2));
      int xF=MIN(FINISHINGX(V1),FINISHINGX(V2));
      V1.window(z0,y0,x0,zF,yF,xF);
      V2.window(z0,y0,x0,zF,yF,xF);
   }


/** Does a radial average of a volume, around the voxel where is the origin.
    A vector radial_mean is returned where:
	 - the first element is the mean of the voxels whose
	   distance to the origin is (0-1),
	 - the second element is the mean of the voxels 
       whose distance to the origin is (1-2)
	 - and so on. 
    A second vector radial_count is returned containing the number of 
    voxels over which each radial average was calculated.
    Sjors nov2003: if rounding=true, element=round(distance); 
         - so the first element is the mean of the voxels whose
	   distance to the origin is (0.5-1.5),
	 - the second element is the mean of the voxels 
       whose distance to the origin is (1.5-2.5)
	 - and so on. */
template <class T>
void radial_average(const matrix3D<T> &m, const matrix1D<int> &center_of_rot,
   matrix1D<T> &radial_mean, matrix1D<int> &radial_count, 
   const bool &rounding=false) {
      matrix1D<double> idx(3);

      /* First determine the maximum distance that one should expect,
         to set the dimension of the radial average vector */
      matrix1D<int> distances(8);
      double z=STARTINGZ(m)-ZZ(center_of_rot);
      double y=STARTINGY(m)-YY(center_of_rot);
      double x=STARTINGX(m)-XX(center_of_rot);
      distances(0)=(int)floor(sqrt(x*x+y*y+z*z));
      x=FINISHINGX(m)-XX(center_of_rot);
      distances(1)=(int)floor(sqrt(x*x+y*y+z*z));
      y=FINISHINGY(m)-YY(center_of_rot);
      distances(2)=(int)floor(sqrt(x*x+y*y+z*z));
      x=STARTINGX(m)-XX(center_of_rot);
      distances(3)=(int)floor(sqrt(x*x+y*y+z*z));
      z=FINISHINGZ(m)-ZZ(center_of_rot);
      distances(4)=(int)floor(sqrt(x*x+y*y+z*z));
      x=FINISHINGX(m)-XX(center_of_rot);
      distances(5)=(int)floor(sqrt(x*x+y*y+z*z));
      y=STARTINGY(m)-YY(center_of_rot);
      distances(6)=(int)floor(sqrt(x*x+y*y+z*z));
      x=STARTINGX(m)-XX(center_of_rot);
      distances(7)=(int)floor(sqrt(x*x+y*y+z*z));

      int dim=(int)CEIL(distances.compute_max())+1;
      if (rounding) dim++;

      // Define the vectors
      radial_mean.resize(dim);
      radial_mean.init_zeros();
      radial_count.resize(dim);
      radial_count.init_zeros();

      /* Perform the radial sum and count pixels that contribute to
         every distance */
      FOR_ALL_ELEMENTS_IN_MATRIX3D(m) {
         ZZ(idx)=k-ZZ(center_of_rot);
         YY(idx)=i-YY(center_of_rot);
         XX(idx)=j-XX(center_of_rot);
         // Determine distance to the center
         int distance;
         if (rounding) distance=(int)ROUND(idx.module());
         else distance=(int)floor(idx.module());

         // Sum te value to the pixels with the same distance  
         radial_mean(distance)+=m(k,i,j);
         // Count the pixel
         radial_count(distance)++;
      }

      // Perform the mean
      FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
         radial_mean(i)/=(T)radial_count(i);
   }

//@}
//@}
//@}
#undef maT
#undef maT1

/* Print shape ------------------------------------------------------------- */
template <class T>
   void VT::print_shape(ostream &out) const {
   out << "Size(Z,Y,X): " << ZSIZE(*this) << "x" << YSIZE(*this) << "x"
       << XSIZE(*this)
       << " k=[" << STARTINGZ(*this) << ".." << FINISHINGZ(*this) << "]"
       << " i=[" << STARTINGY(*this) << ".." << FINISHINGY(*this) << "]"
       << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
}

/* Get size--- ------------------------------------------------------------- */
template <class T>
   void VT::get_size(int *size) const
   {size[0]=xdim; size[1]=ydim; size[2]=zdim;}

/* Outside ----------------------------------------------------------------- */
template <class T>
   bool VT::outside(const matrix1D<double> &v) const {
   if (XSIZE(v)<3)
      REPORT_ERROR(1,"Outside: index vector has got not enough components");
   return (XX(v)<STARTINGX(*this) || XX(v)>FINISHINGX(*this) ||
           YY(v)<STARTINGY(*this) || YY(v)>FINISHINGY(*this) ||
           ZZ(v)<STARTINGZ(*this) || ZZ(v)>FINISHINGZ(*this));
}

template <class T>
   bool VT::outside(int k, int i, int j) const {
   return (j<STARTINGX(*this) || j>FINISHINGX(*this) ||
           i<STARTINGY(*this) || i>FINISHINGY(*this) ||
           k<STARTINGZ(*this) || k>FINISHINGZ(*this));
}

/* Intersects -------------------------------------------------------------- */
template <class T>
   bool VT::intersects(const VT &m) const {
      return intersects(STARTINGZ(m),STARTINGY(m),STARTINGX(m),
         XSIZE(m)-1,YSIZE(m)-1,ZSIZE(m)-1);
}

template <class T>
   bool VT::intersects(const matrix1D<double> &corner1,
      const matrix1D<double> &corner2) const {
       if (XSIZE(corner1)!=2 || XSIZE(corner2)!=2)
          REPORT_ERROR(1002,"intersects 1D: corner sizes are not 1");
       return intersects(XX(corner1),YY(corner1),ZZ(corner1),
          XX(corner2)-XX(corner1),YY(corner2)-YY(corner1),
          ZZ(corner2)-ZZ(corner1));
}


template <class T>
   bool VT::intersects(double x0, double y0, double z0, double xdim, double ydim,
      double zdim) const {
   SPEED_UP_temps;
   spduptmp0=MAX(STARTINGZ(*this), z0);
   spduptmp1=MIN(FINISHINGZ(*this),z0+zdim);
   if (spduptmp0>spduptmp1) return false;
   
   spduptmp0=MAX(STARTINGY(*this), y0);
   spduptmp1=MIN(FINISHINGY(*this),y0+ydim);
   if (spduptmp0>spduptmp1) return false;

   spduptmp0=MAX(STARTINGX(*this), x0);
   spduptmp1=MIN(FINISHINGX(*this),x0+xdim);
   if (spduptmp0>spduptmp1) return false;
   return true;
}

/* Corner ------------------------------------------------------------------ */
template <class T>
   bool VT::isCorner(const matrix1D<double> &v) {
   if (XSIZE(v)<3)
      REPORT_ERROR(1,"isCorner: index vector has got not enough components");
   return ((ZZ(v)==STARTINGZ(*this)  && XX(v)==STARTINGX(*this)  && YY(v)==STARTINGY(*this))  ||
           (ZZ(v)==STARTINGZ(*this)  && XX(v)==STARTINGX(*this)  && YY(v)==FINISHINGY(*this)) ||
           (ZZ(v)==STARTINGZ(*this)  && XX(v)==FINISHINGX(*this) && YY(v)==STARTINGY(*this))  ||
           (ZZ(v)==STARTINGZ(*this)  && XX(v)==FINISHINGX(*this) && YY(v)==FINISHINGY(*this)) ||
           (ZZ(v)==FINISHINGZ(*this) && XX(v)==STARTINGX(*this)  && YY(v)==STARTINGY(*this))  ||
           (ZZ(v)==FINISHINGZ(*this) && XX(v)==STARTINGX(*this)  && YY(v)==FINISHINGY(*this)) ||
           (ZZ(v)==FINISHINGZ(*this) && XX(v)==FINISHINGX(*this) && YY(v)==STARTINGY(*this))  ||
           (ZZ(v)==FINISHINGZ(*this) && XX(v)==FINISHINGX(*this) && YY(v)==FINISHINGY(*this)));
}

/* Border ------------------------------------------------------------------ */
template <class T>
   bool VT::isBorder(const matrix1D<int> &v) 
{
   if (XSIZE(v)<3)
      REPORT_ERROR(1,"isBorder: index vector has got not enough components");
   return  isBorder(ZZ(v),YY(v),XX(v));
}

template <class T>
   bool VT::isBorder(int k,int i,int j) 
{
   return (j==STARTINGX(*this)  || j==FINISHINGX(*this)  ||
           k==STARTINGZ(*this)  || k==FINISHINGZ(*this)  ||
           i==STARTINGY(*this)  || i==FINISHINGY(*this));
}

/* Patch ------------------------------------------------------------------- */
template <class T>
   void VT::patch(const VT &patch_array, char operation) {
      SPEED_UP_temps;
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(patch_array,*this)
         switch (operation) {
            case '=': VOL_ELEM(*this,k,i,j) =VOL_ELEM(patch_array,k,i,j); break;
            case '+': VOL_ELEM(*this,k,i,j)+=VOL_ELEM(patch_array,k,i,j); break;
            case '-': VOL_ELEM(*this,k,i,j)-=VOL_ELEM(patch_array,k,i,j); break;
            case '*': VOL_ELEM(*this,k,i,j)*=VOL_ELEM(patch_array,k,i,j); break;
            case '/': VOL_ELEM(*this,k,i,j)/=VOL_ELEM(patch_array,k,i,j); break;
         }
}

/* Output stream ----------------------------------------------------------- */
template <class T>
   ostream& operator << (ostream& ostrm, const VT& v) {   
   if (v.xdim==0)
      ostrm << "NULL matrix3D\n";
   else
      ostrm << endl;
      double max_val=ABS(MULTIDIM_ELEM(v,0));
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(v)
         max_val=MAX(max_val,ABS(MULTIDIM_ELEM(v,i)));
      int prec=best_prec(max_val,10);
      
      for (int k=STARTINGZ(v); k<=FINISHINGZ(v); k++) {
         ostrm << "Slice No. " << k << endl;
         for (int i=STARTINGY(v); i<=FINISHINGY(v); i++) {
            for (int j=STARTINGX(v); j<=FINISHINGX(v); j++) {
               ostrm << FtoA((double)VOL_ELEM(v,k,i,j),10,prec) << ' ';
            }
            ostrm << endl;
         }
      }
      
   return ostrm;
}

/* Apply a geometrical transformation -------------------------------------- */
// It is translated from function transforma in Lib/varios2.c
// We should check which one performs better.
//#define DEBUG
template <class T>
   void apply_geom(VT &V2, matrix2D<double> A, const VT &V1, bool inv,
      bool wrap) {
   int m1, n1, o1, m2, n2, o2;
   double x, y, z, xp, yp, zp;
   double minxp, minyp, maxxp, maxyp, minzp, maxzp;
   int   cen_x, cen_y, cen_z, cen_xp, cen_yp, cen_zp;
   double wx, wy, wz; // Weights in X,Y,Z directions for bilinear interpolation
   
   if ((XSIZE(A)!=4) || (YSIZE(A)!=4))
      REPORT_ERROR(1102,"Apply_geom3D: geometrical transformation is not 4x4");
   if (A.IsIdent()) {V2=V1; return;}
   if (XSIZE(V1)==0)  {V2.clear(); return;}
   
   if (!inv) A=A.inv();
   
   if (XSIZE(V2)==0) V2.resize(V1);
   
   // For scalings the output matrix is resized outside to the final
   // size instead of being resized inside the routine with the
   // same size as the input matrix
   if (XSIZE(V2)==0) V2.resize(V1);
   
   // Find center of matrix3D
   cen_z  =(int)(V2.zdim/2);
   cen_y  =(int)(V2.ydim/2);
   cen_x  =(int)(V2.xdim/2);
   cen_zp =(int)(V1.zdim/2);
   cen_yp =(int)(V1.ydim/2);
   cen_xp =(int)(V1.xdim/2);
   minxp  = -cen_xp;
   minyp  = -cen_yp;
   minzp  = -cen_zp;
   maxxp  = V1.xdim-cen_xp-1;
   maxyp  = V1.ydim-cen_yp-1;
   maxzp  = V1.zdim-cen_zp-1;
   #ifdef DEBUG
      cout << "Geometry 2 center=("
           << cen_z  << "," << cen_y  << "," << cen_x  << ")\n"
           << "Geometry 1 center=("
           << cen_zp << "," << cen_yp << "," << cen_xp << ")\n"
	   << "           min=("
           << minzp  << "," << minyp  << "," << minxp  << ")\n"
	   << "           max=("
           << maxzp  << "," << maxyp  << "," << maxxp  << ")\n"
      ;
   #endif
    
   // Now we go from the output matrix3D to the input matrix3D, ie, for any voxel
   // in the output matrix3D we calculate which are the corresponding ones in
   // the original matrix3D, make an interpolation with them and put this value
   // at the output voxel
   
   // V2 is not initialised to 0 because all its pixels are rewritten
   for (int k=0; k<V2.zdim; k++)
      for (int i=0; i<V2.ydim; i++) {
         // Calculate position of the beginning of the row in the output matrix3D
         x= -cen_x;
         y=i-cen_y;
         z=k-cen_z;

         // Calculate this position in the input image according to the
         // geometrical transformation
         // they are related by
         // coords_output(=x,y) = A * coords_input (=xp,yp)
         xp=x*dMij(A,0,0) + y*dMij(A,0,1) + z*dMij(A,0,2) + dMij(A,0,3);
         yp=x*dMij(A,1,0) + y*dMij(A,1,1) + z*dMij(A,1,2) + dMij(A,1,3);
         zp=x*dMij(A,2,0) + y*dMij(A,2,1) + z*dMij(A,2,2) + dMij(A,2,3);

         for (int j=0; j<V2.xdim; j++) {
            bool interp;
            T tmp;

            #ifdef DEBUG
	       bool show_debug=false;
	       if ((i==0 && j==0 && k==0) ||
	           (i==V2.ydim-1 && j==V2.xdim-1 && k==V2.zdim-1))
		   show_debug=true;
	       if (show_debug)
        	  cout << "(x,y,z)-->(xp,yp,zp)= "
		       << "(" << x  << "," << y  << "," << z  << ") "
		       << "(" << xp << "," << yp << "," << zp << ")\n";
            #endif

            // If the point is outside the volume, apply a periodic extension
            // of the volume, what exits by one side enters by the other
            interp=true;
            if (wrap) {
               if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	           xp>maxxp+XMIPP_EQUAL_ACCURACY) xp=realWRAP(xp, minxp-0.5, maxxp+0.5);
               if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	           yp>maxyp+XMIPP_EQUAL_ACCURACY) yp=realWRAP(yp, minyp-0.5, maxyp+0.5);
               if (zp<minzp-XMIPP_EQUAL_ACCURACY ||
	           zp>maxzp+XMIPP_EQUAL_ACCURACY) zp=realWRAP(zp, minzp-0.5, maxzp+0.5);
            } else {
               if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	           xp>maxxp+XMIPP_EQUAL_ACCURACY) interp=false;
               if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	           yp>maxyp+XMIPP_EQUAL_ACCURACY) interp=false;
               if (zp<minzp-XMIPP_EQUAL_ACCURACY ||
	           zp>maxzp+XMIPP_EQUAL_ACCURACY) interp=false;
            }

            if (interp) {            
               // Calculate the integer position in input volume, be careful
               // that it is not the nearest but the one at the top left corner
               // of the interpolation square. Ie, (0.7,0.7) would give (0,0)
               // Calculate also weights for point m1+1,n1+1
               wx=xp+cen_xp; m1=(int) wx; wx=wx-m1; m2=m1+1;
               wy=yp+cen_yp; n1=(int) wy; wy=wy-n1; n2=n1+1;
               wz=zp+cen_zp; o1=(int) wz; wz=wz-o1; o2=o1+1;
	       
	       #ifdef DEBUG
		  if (show_debug) {
        	     cout << "After wrapping(xp,yp,zp)= "
			  << "(" << xp << "," << yp << "," << zp << ")\n";
	             cout << "(m1,n1,o1)-->(m2,n2,o2)="
			  << "(" << m1 << "," << n1 << "," << o1 << ") "
			  << "(" << m2 << "," << n2 << "," << o2 << ")\n";
		     cout << "(wx,wy,wz)=" 
			  << "(" << wx << "," << wy << "," << wz << ")\n";
		  }
	       #endif

               // Perform interpolation
               // if wx == 0 means that the rightest point is useless for this
               // interpolation, and even it might not be defined if m1=xdim-1
               // The same can be said for wy.
                                              tmp  = (T) ((1-wz)*(1-wy)*(1-wx)*dVkij(V1,o1,n1,m1));
               if (wx!=0 && m2<V1.xdim)       tmp += (T) ((1-wz)*(1-wy)*   wx *dVkij(V1,o1,n1,m2));
               if (wy!=0 && n2<V1.ydim){      tmp += (T) ((1-wz)*   wy *(1-wx)*dVkij(V1,o1,n2,m1));
                  if (wx!=0 && m2<V1.xdim)    tmp += (T) ((1-wz)*   wy *   wx *dVkij(V1,o1,n2,m2));
               }
               if (wz!=0 && o2<V1.zdim){
                                              tmp += (T) (   wz *(1-wy)*(1-wx)*dVkij(V1,o2,n1,m1));
                  if (wx!=0 && m2<V1.xdim)    tmp += (T) (   wz *(1-wy)*   wx *dVkij(V1,o2,n1,m2));
                  if (wy!=0 && n2<V1.ydim){   tmp += (T) (   wz *   wy *(1-wx)*dVkij(V1,o2,n2,m1));
                     if (wx!=0 && m2<V1.xdim) tmp += (T) (   wz *   wy *   wx *dVkij(V1,o2,n2,m2));
                  }
               }
               dVkij(V2,k,i,j) = tmp;
	       #ifdef DEBUG
		  if (show_debug)
	             cout << "tmp1=" << dVkij(V1,o1,n1,m1) << " " << (T) ((1-wz)*(1-wy)*(1-wx)*dVkij(V1,o1,n1,m1)) << endl
        		  << "tmp2=" << dVkij(V1,o1,n1,m2) << " " << (T) ((1-wz)*(1-wy)*   wx *dVkij(V1,o1,n1,m2)) << endl
        		  << "tmp3=" << dVkij(V1,o1,n2,m1) << " " << (T) ((1-wz)*   wy *(1-wx)*dVkij(V1,o1,n2,m1)) << endl
        		  << "tmp4=" << dVkij(V1,o1,n2,m2) << " " << (T) ((1-wz)*   wy *   wx *dVkij(V1,o1,n2,m2)) << endl
        		  << "tmp5=" << dVkij(V1,o2,n1,m1) << " " << (T) (   wz *(1-wy)*(1-wx)*dVkij(V1,o2,n1,m1)) << endl
        		  << "tmp6=" << dVkij(V1,o2,n1,m2) << " " << (T) (   wz *(1-wy)*   wx *dVkij(V1,o2,n1,m2)) << endl
        		  << "tmp7=" << dVkij(V1,o2,n2,m1) << " " << (T) (   wz *   wy *(1-wx)*dVkij(V1,o2,n2,m1)) << endl
        		  << "tmp8=" << dVkij(V1,o2,n2,m2) << " " << (T) (   wz *   wy *   wx *dVkij(V1,o2,n2,m2)) << endl
			  << "tmp= " << tmp << endl;
	       #endif
            }

            // Compute new point inside input image
            xp += dMij(A,0,0);
            yp += dMij(A,1,0);
            zp += dMij(A,2,0);
         }
      }
}
#undef DEBUG

//#define DEBUG
template <class T>
   void apply_geom_Bspline(VT &V2, matrix2D<double> A, const VT &V1,
      int Splinedegree, bool inv, bool wrap, T outside) {
   int m1, n1, o1, m2, n2, o2;
   double x, y, z, xp, yp, zp;
   double minxp, minyp, maxxp, maxyp, minzp, maxzp;
   int   cen_x, cen_y, cen_z, cen_xp, cen_yp, cen_zp;
   
   if ((XSIZE(A)!=4) || (YSIZE(A)!=4))
      REPORT_ERROR(1102,"Apply_geom3D: geometrical transformation is not 4x4");
   if (A.IsIdent()) {V2=V1; return;}
   if (XSIZE(V1)==0)  {V2.clear(); return;}
   
   if (!inv) A=A.inv();
   
   if (XSIZE(V2)==0) V2.resize(V1);
   
   // For scalings the output matrix is resized outside to the final
   // size instead of being resized inside the routine with the
   // same size as the input matrix
   if (XSIZE(V2)==0) V2.resize(V1);
   
   // Find center of matrix3D
   cen_z  =(int)(V2.zdim/2);
   cen_y  =(int)(V2.ydim/2);
   cen_x  =(int)(V2.xdim/2);
   cen_zp =(int)(V1.zdim/2);
   cen_yp =(int)(V1.ydim/2);
   cen_xp =(int)(V1.xdim/2);
   minxp  = -cen_xp;
   minyp  = -cen_yp;
   minzp  = -cen_zp;
   maxxp  = V1.xdim-cen_xp-1;
   maxyp  = V1.ydim-cen_yp-1;
   maxzp  = V1.zdim-cen_zp-1;
   #ifdef DEBUG
      cout << "Geometry 2 center=("
           << cen_z  << "," << cen_y  << "," << cen_x  << ")\n"
           << "Geometry 1 center=("
           << cen_zp << "," << cen_yp << "," << cen_xp << ")\n"
	   << "           min=("
           << minzp  << "," << minyp  << "," << minxp  << ")\n"
	   << "           max=("
           << maxzp  << "," << maxyp  << "," << maxxp  << ")\n"
      ;
   #endif
    
   // Build the B-spline coefficients
   matrix3D<double> Bcoeffs;
   V1.produce_spline_coeffs(Bcoeffs,Splinedegree);
   STARTINGX(Bcoeffs)=(int)minxp;
   STARTINGY(Bcoeffs)=(int)minyp;
   STARTINGZ(Bcoeffs)=(int)minzp;

   // Now we go from the output matrix3D to the input matrix3D, ie, for any voxel
   // in the output matrix3D we calculate which are the corresponding ones in
   // the original matrix3D, make an interpolation with them and put this value
   // at the output voxel
   
   // V2 is not initialised to 0 because all its pixels are rewritten
   for (int k=0; k<V2.zdim; k++)
      for (int i=0; i<V2.ydim; i++) {
         // Calculate position of the beginning of the row in the output matrix3D
         x= -cen_x;
         y=i-cen_y;
         z=k-cen_z;

         // Calculate this position in the input image according to the
         // geometrical transformation
         // they are related by
         // coords_output(=x,y) = A * coords_input (=xp,yp)
         xp=x*dMij(A,0,0) + y*dMij(A,0,1) + z*dMij(A,0,2) + dMij(A,0,3);
         yp=x*dMij(A,1,0) + y*dMij(A,1,1) + z*dMij(A,1,2) + dMij(A,1,3);
         zp=x*dMij(A,2,0) + y*dMij(A,2,1) + z*dMij(A,2,2) + dMij(A,2,3);

         for (int j=0; j<V2.xdim; j++) {
            bool interp;
            T tmp;

            #ifdef DEBUG
	       bool show_debug=false;
	       if ((i==0 && j==0 && k==0) ||
	           (i==V2.ydim-1 && j==V2.xdim-1 && k==V2.zdim-1))
		   show_debug=true;
	       if (show_debug)
        	  cout << "(x,y,z)-->(xp,yp,zp)= "
		       << "(" << x  << "," << y  << "," << z  << ") "
		       << "(" << xp << "," << yp << "," << zp << ")\n";
            #endif

            // If the point is outside the volume, apply a periodic extension
            // of the volume, what exits by one side enters by the other
            interp=true;
            if (wrap) {
               if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	           xp>maxxp+XMIPP_EQUAL_ACCURACY) xp=realWRAP(xp, minxp-0.5, maxxp+0.5);
               if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	           yp>maxyp+XMIPP_EQUAL_ACCURACY) yp=realWRAP(yp, minyp-0.5, maxyp+0.5);
               if (zp<minzp-XMIPP_EQUAL_ACCURACY ||
	           zp>maxzp+XMIPP_EQUAL_ACCURACY) zp=realWRAP(zp, minzp-0.5, maxzp+0.5);
            } else {
               if (xp<minxp-XMIPP_EQUAL_ACCURACY ||
	           xp>maxxp+XMIPP_EQUAL_ACCURACY) interp=false;
               if (yp<minyp-XMIPP_EQUAL_ACCURACY ||
	           yp>maxyp+XMIPP_EQUAL_ACCURACY) interp=false;
               if (zp<minzp-XMIPP_EQUAL_ACCURACY ||
	           zp>maxzp+XMIPP_EQUAL_ACCURACY) interp=false;
            }

            if (interp) {
               dVkij(V2,k,i,j) = (T)Bcoeffs.interpolated_elem_as_Bspline(
	          xp,yp,zp,Splinedegree);
            } else dVkij(V2,k,i,j) = outside;

            // Compute new point inside input image
            xp += dMij(A,0,0);
            yp += dMij(A,1,0);
            zp += dMij(A,2,0);
         }
      }
}

/* Statistics in region ---------------------------------------------------- */
template <class T>
void VT::compute_stats(double &avg, double &stddev, T &min_val, T &max_val,
   const matrix1D<double> &corner1, const matrix1D<double> &corner2) const {
   min_val=max_val=(*this)(corner1);
   matrix1D<double> r(3);
   double N=0, sum=0, sum2=0;
   FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1,corner2) {
      sum+=(*this)(r); sum2+=(*this)(r)*(*this)(r); N++;
      if      ((*this)(r)<min_val) min_val=(*this)(r);
      else if ((*this)(r)>max_val) max_val=(*this)(r);
   }
   if (N!=0) {
      avg=sum/N;
      stddev=sqrt(sum2/N-avg*avg);
   } else {avg=stddev=0;}
}

/* Minimum and maximum in region ------------------------------------------- */
template <class T>
void VT::compute_double_minmax(double &min_val, double &max_val,
   const matrix1D<double> &corner1, const matrix1D<double> &corner2) const {
   min_val=max_val=(*this)(corner1);
   matrix1D<double> r(3);
   FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1,corner2) {
      if      ((*this)(r)<min_val) min_val=(*this)(r);
      else if ((*this)(r)>max_val) max_val=(*this)(r);
   }
}

/* Center of mass ---------------------------------------------------------- */
template <class T>
   void VT::center_of_mass(matrix1D<double> &center, void * mask) {
      center.init_zeros(3);
      double mass=0;
      matrix3D<int> *imask=(matrix3D<int> *) mask;
      FOR_ALL_ELEMENTS_IN_MATRIX3D(*this) {
         if (imask==NULL || VOL_ELEM(*imask,k,i,j)) {
            XX(center)+=j*VOL_ELEM(*this,k,i,j);
            YY(center)+=i*VOL_ELEM(*this,k,i,j);
            ZZ(center)+=k*VOL_ELEM(*this,k,i,j);
	    mass+=VOL_ELEM(*this,k,i,j);
      	 }
      }
      if (mass!=0) center/=mass;
   }

template <>
complex<double> matrix3D< complex<double> >::interpolated_elem(
   double x, double y, double z, complex<double> outside_value);

#endif
