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
#include "../grids.hh"

#include <stdio.h>
#include <string.h> // for memcpy

/*****************************************************************************/
/* Simple Grids                                                              */
/*****************************************************************************/
// Constructor -------------------------------------------------------------
SimpleGrid::SimpleGrid() {
   basis.init_identity(3);
   inv_basis.init_identity(3);
   origin        = vector_R3(0.,0.,0.);
   lowest        = vector_R3(-5.,-5.,-5.);
   highest       = -lowest;
   relative_size = 1;
   R2            = -1;
}

SimpleGrid::SimpleGrid(const SimpleGrid &SG) {
   basis         = SG.basis;
   inv_basis     = SG.inv_basis;
   lowest        = SG.lowest;
   highest       = SG.highest;
   relative_size = SG.relative_size;
   origin        = SG.origin;
   R2            = SG.R2;
}

// Cout --------------------------------------------------------------------
ostream& operator <<(ostream& o, const SimpleGrid &grid) {
   o << "   Simple Grid -----" << endl;
   o << "   Vector 1: " << ((grid.basis).Col(0)).transpose() << endl;
   o << "   Vector 2: " << ((grid.basis).Col(1)).transpose() << endl;
   o << "   Vector 3: " << ((grid.basis).Col(2)).transpose() << endl;
   o << "   Relative size:         " << grid.relative_size       << endl;
   o << "   Interest radius²:      " << grid.R2                  << endl;
   o << "   Origin (univ.coords)   " << grid.origin.transpose()  << endl;
   o << "   Highest (grid. coord)  " << grid.highest.transpose() << endl;
   o << "   Lowest (grid. coord)   " << grid.lowest.transpose()  << endl;
   return o;
}

// Assignment --------------------------------------------------------------
SimpleGrid& SimpleGrid::operator = (const SimpleGrid &SG) {
   if (&SG!=this) {
      basis         = SG.basis;
      inv_basis     = SG.inv_basis;
      lowest        = SG.lowest;
      highest       = SG.highest;
      relative_size = SG.relative_size;
      R2            = SG.R2;
      origin        = SG.origin;
   }
   return *this;
}

// Prepare Grid ------------------------------------------------------------
void SimpleGrid::prepare_grid() _THROW {
   // Compute matrix for inverse basis conversion
   try {
      inv_basis=basis.inv();
   }
   catch (Xmipp_error error) {
      REPORT_ERROR(3001,"The grid vectors are not a true 3D coordinate system");
   }
   lowest.setCol();
   highest.setCol();
   origin.setCol();
}

// Minimum size ------------------------------------------------------------
void Grid::voxel_corners(matrix1D<double> &Gcorner1, matrix1D<double> &Gcorner2,
   const matrix2D<double> *V) const {
   matrix1D<double> SGcorner1(3), SGcorner2(3);     // Subgrid corners
   SPEED_UP_temps;
   
   // Look for the lowest and highest volume coordinate
   Gcorner1.resize(3);  // lowest and highest coord.
   Gcorner2.resize(3);
   for (int n=0; n<GridsNo(); n++) {
      // Find box for this grid
      bool first; first=TRUE;
      for (int k=(int)ZZ(LG[n].lowest); k<=ZZ(LG[n].highest); k++)
	 for (int i=(int)YY(LG[n].lowest); i<=YY(LG[n].highest); i++)
	    for (int j=(int)XX(LG[n].lowest); j<=XX(LG[n].highest); j++) {
	       matrix1D<double> grid_index(3), univ_position(3);
	       VECTOR_R3(grid_index,j,i,k);
	       LG[n].grid2universe(grid_index,univ_position);
	       if (V!=NULL) {M3x3_BY_V3x1(univ_position,*V,univ_position);}
	       if (!LG[n].is_interesting(univ_position)) continue;
	       if (!first) {
        	  XX(SGcorner1)=MIN(XX(SGcorner1),XX(univ_position));
        	  YY(SGcorner1)=MIN(YY(SGcorner1),YY(univ_position));
        	  ZZ(SGcorner1)=MIN(ZZ(SGcorner1),ZZ(univ_position));

        	  XX(SGcorner2)=MAX(XX(SGcorner2),XX(univ_position));
        	  YY(SGcorner2)=MAX(YY(SGcorner2),YY(univ_position));
        	  ZZ(SGcorner2)=MAX(ZZ(SGcorner2),ZZ(univ_position));
	       } else {
	          SGcorner2=SGcorner1=univ_position;
		  first=FALSE;
	       }
      	    }

      // Compare with the rest of the grids
      if (n!=0) {
         XX(Gcorner1)=MIN(XX(Gcorner1),XX(SGcorner1));
         YY(Gcorner1)=MIN(YY(Gcorner1),YY(SGcorner1));
         ZZ(Gcorner1)=MIN(ZZ(Gcorner1),ZZ(SGcorner1));

         XX(Gcorner2)=MAX(XX(Gcorner2),XX(SGcorner2));
         YY(Gcorner2)=MAX(YY(Gcorner2),YY(SGcorner2));
         ZZ(Gcorner2)=MAX(ZZ(Gcorner2),ZZ(SGcorner2));
      } else {
         Gcorner1=SGcorner1;
         Gcorner2=SGcorner2;
      }

      #ifdef DEBUG
      cout << LG[n];
      cout << "SGcorner1 " << SGcorner1.transpose() << endl;
      cout << "SGcorner2 " << SGcorner2.transpose() << endl;
      cout << "Gcorner1  " << Gcorner1.transpose() << endl;
      cout << "Gcorner2  " << Gcorner2.transpose() << endl;
      #endif
   }
}

/*****************************************************************************/
/* Some useful Grids                                                         */
/*****************************************************************************/
/* Create CC Simple grid with a given origin ------------------------------- */
SimpleGrid Create_CC_grid(double relative_size, const matrix1D<double> &corner1,
   const matrix1D<double> &corner2, const matrix1D<double> &origin) {
   SimpleGrid    grid;
   
   // The vectors of the grid are the default ones of (1,0,0), (0,1,0),
   // and (0,0,1), and its inverse matrix is already computed
   grid.relative_size = relative_size;
   grid.origin        = origin;

   // Compute the lowest and highest indexes inside the grid
   grid.universe2grid(corner1,grid.lowest);  grid.lowest.FLOORnD();
   grid.universe2grid(corner2,grid.highest); grid.highest.CEILnD();
   
   grid.R2=-1;
   return grid;
}

/* Create CC grid ---------------------------------------------------------- */
Grid Create_CC_grid(double relative_size, const matrix1D<double> &corner1,
   const matrix1D<double> &corner2) {
   Grid            result;
   SimpleGrid      aux_grid;
   
   matrix1D<double> origin=(corner1+corner2)/2; origin.ROUNDnD();
   aux_grid=Create_CC_grid(relative_size,corner1,corner2,origin);
   result.add_grid(aux_grid);
   return result;
}

Grid Create_CC_grid(double relative_size, int Zdim, int Ydim, int Xdim) {
   Grid            result;
   SimpleGrid      aux_grid;
   
   matrix1D<double> origin=
      vector_R3((double)FLOOR(Xdim/2.0),(double)FLOOR(Ydim/2.0),
         (double)FLOOR(Zdim/2.0));
   aux_grid=Create_CC_grid(relative_size,-origin,
      vector_R3((double)Xdim,(double)Ydim,(double)Zdim)-origin-1,origin);
   result.add_grid(aux_grid);

   return result;
}

/* Create BCC grid --------------------------------------------------------- */
Grid Create_BCC_grid(double relative_size, const matrix1D<double> &corner1,
   const matrix1D<double> &corner2) {
   Grid             result;
   SimpleGrid       aux_grid;
   matrix1D<double> origin=(corner1+corner2)/2; origin.ROUNDnD();
   
   //Even Slice
   //    0 1 2 3 4 5 6 7 8 9 10 11 12 (Col)
   //  0 A   A   A   A   A   A     A
   //  1
   //  2 A   A   A   A   A   A     A
   //  3
   //  4 A   A   A   A   A   A     A
   //  5
   //  6 A   A   A   A   A   A     A
   //  7
   //  8 A   A   A   A   A   A     A
   //  9
   // 10 A   A   A   A   A   A     A
   // 11
   // 12 A   A   A   A   A   A     A
   //(Row)
   //
   //Odd Slice
   //    0 1 2 3 4 5 6 7 8 9 10 11 12 (Col)
   //  0
   //  1   B   B   B   B   B    B
   //  2
   //  3   B   B   B   B   B    B
   //  4
   //  5   B   B   B   B   B    B
   //  6
   //  7   B   B   B   B   B    B
   //  8
   //  9   B   B   B   B   B    B
   // 10
   // 11   B   B   B   B   B    B
   // 12
   //(Row)
   
   // Grid A
   aux_grid=Create_CC_grid(relative_size,corner1,corner2,origin);
   result.add_grid(aux_grid);
   
   // Grid B
   origin=origin+relative_size/2*vector_R3(1.,1.,1.);
   aux_grid=Create_CC_grid(relative_size,corner1,corner2,origin);
   result.add_grid(aux_grid);

   return result;
}

/* Create FCC grid --------------------------------------------------------- */
Grid Create_FCC_grid(double relative_size, const matrix1D<double> &corner1,
   const matrix1D<double> &corner2) {

   Grid             result;
   SimpleGrid       aux_grid;
   matrix1D<double> aux_origin;
   matrix1D<double> cornerb;
   matrix1D<double> origin=(corner1+corner2)/2; origin.ROUNDnD();
   
   //Even Slice
   //    0 1 2 3 4 5 6 7 8 9 10 11 12 (Col)
   //  0 A   A   A   A   A   A     A
   //  1   B   B   B   B   B    B
   //  2 A   A   A   A   A   A     A
   //  3   B   B   B   B   B    B
   //  4 A   A   A   A   A   A     A
   //  5   B   B   B   B   B    B
   //  6 A   A   A   A   A   A     A
   //  7   B   B   B   B   B    B
   //  8 A   A   A   A   A   A     A
   //  9   B   B   B   B   B    B
   // 10 A   A   A   A   A   A     A
   // 11   B   B   B   B   B    B
   // 12 A   A   A   A   A   A     A
   //(Row)
   //
   //Odd Slice
   //    0 1 2 3 4 5 6 7 8 9 10 11 12 (Col)
   //  0   C   C   C   C   C    C
   //  1 D   D   D   D   D   D     D
   //  2   C   C   C   C   C    C
   //  3 D   D   D   D   D   D     D
   //  4   C   C   C   C   C    C
   //  5 D   D   D   D   D   D     D
   //  6   C   C   C   C   C    C
   //  7 D   D   D   D   D   D     D
   //  8   C   C   C   C   C    C
   //  9 D   D   D   D   D   D     D
   // 10   C   C   C   C   C    C
   // 11 D   D   D   D   D   D     D
   // 12   C   C   C   C   C    C
   //(Row)
   
   // Grid A
   aux_grid=Create_CC_grid(relative_size,corner1,corner2,origin);
   result.add_grid(aux_grid);
   // Grid D
   aux_origin=origin+relative_size/2*vector_R3(0.,1.,1.);
   aux_grid=Create_CC_grid(relative_size,corner1,corner2,aux_origin);
   result.add_grid(aux_grid);
   // Grid C
   aux_origin=origin+relative_size/2*vector_R3(1.,0.,1.);
   aux_grid=Create_CC_grid(relative_size,corner1,corner2,aux_origin);
   result.add_grid(aux_grid);
   // Grid B
   cornerb=corner2;
   cornerb(0)=cornerb(0)-1;
   cornerb(1)=cornerb(1)-1;
   aux_origin=origin+relative_size/2*vector_R3(1.,1.,0.);
   aux_grid=Create_CC_grid(relative_size,corner1,cornerb,aux_origin);
   result.add_grid(aux_grid);
   
   return result;
}
#undef MULTIPLY_CC_GRID_BY_TWO

/* CC grid with region of interest ----------------------------------------- */
//#define DEBUG
SimpleGrid Create_grid_within_sphere(double relative_size,
   const matrix1D<double> &origin,
   const matrix1D<double> &X, const matrix1D<double> &Y,
   const matrix1D<double> &Z, double R2) {
   SimpleGrid    grid;
   double R=sqrt(R2);

   grid.set_X(X);
   grid.set_Y(Y);
   grid.set_Z(Z);
   grid.relative_size = relative_size;
   grid.origin        = origin;
   grid.R2=R2;
   grid.lowest.init_zeros(3);
   grid.highest.init_zeros(3);
   grid.prepare_grid();

   // Find grid limits
   int iR=CEIL(R);
   matrix1D<double> univ_position(3), grid_position(3);
   for (int k=-iR; k<=iR; k++)
      for (int i=-iR; i<=iR; i++)
	 for (int j=-iR; j<=iR; j++) {
	    VECTOR_R3(univ_position,j,i,k);
	    if (univ_position.module()>R) continue;
	    grid.universe2grid(univ_position,grid_position);
	    XX(grid.lowest)=MIN(XX(grid.lowest),FLOOR(XX(grid_position)));
	    YY(grid.lowest)=MIN(YY(grid.lowest),FLOOR(YY(grid_position)));
	    ZZ(grid.lowest)=MIN(ZZ(grid.lowest),FLOOR(ZZ(grid_position)));
	    XX(grid.highest)=MAX(XX(grid.highest),CEIL(XX(grid_position)));
	    YY(grid.highest)=MAX(YY(grid.highest),CEIL(YY(grid_position)));
	    ZZ(grid.highest)=MAX(ZZ(grid.highest),CEIL(ZZ(grid_position)));
	 }

   #ifdef DEBUG
      cout << "Sphere radius = " << R << endl
           << "relative size = " << relative_size << endl
	   << "X module      = " << X.module() << endl
	   << "Y module      = " << Y.module() << endl
	   << "Z module      = " << Z.module() << endl
	   << grid
      ;
   #endif
   return grid;
}
#undef DEBUG

/* Create CC grid ---------------------------------------------------------- */
Grid Create_CC_grid(double relative_size, double R) {
   Grid            result;
   SimpleGrid      aux_grid;
   
   matrix1D<double> origin(3); origin.init_zeros();
   matrix1D<double> x(3), y(3), z(3);
   VECTOR_R3(x,1,0,0);
   VECTOR_R3(y,0,1,0);
   VECTOR_R3(z,0,0,1);
   aux_grid=Create_grid_within_sphere(relative_size,origin,x,y,z,R*R);
   result.add_grid(aux_grid);
   return result;
}

/* Create BCC grid --------------------------------------------------------- */
Grid Create_BCC_grid(double relative_size, double R) {
   Grid            result;
   SimpleGrid      aux_grid;
   
   matrix1D<double> origin(3); origin.init_zeros();
   matrix1D<double> x(3), y(3), z(3);
   VECTOR_R3(x,0.5,0.5,-0.5);
   VECTOR_R3(y,0.5,-0.5,0.5);
   VECTOR_R3(z,-0.5,0.5,0.5);
   aux_grid=Create_grid_within_sphere(relative_size,origin,x,y,z,R*R);
   result.add_grid(aux_grid);
   return result;
}

/* Create FCC grid --------------------------------------------------------- */
Grid Create_FCC_grid(double relative_size, double R) {
   Grid            result;
   SimpleGrid      aux_grid;
   
   matrix1D<double> origin(3); origin.init_zeros();
   matrix1D<double> x(3), y(3), z(3);
   VECTOR_R3(x,0.5,0.5,0);
   VECTOR_R3(y,0.5,0,0.5);
   VECTOR_R3(z,0,0.5,0.5);
   aux_grid=Create_grid_within_sphere(relative_size,origin,x,y,z,R*R);
   result.add_grid(aux_grid);
   return result;
}

/*****************************************************************************/
/* Grid Volumes                                                              */
/*****************************************************************************/
/* Assignment -------------------------------------------------------------- */
template <class T>
GridVolumeT<T> & GridVolumeT<T>::operator = (const GridVolumeT<T>& RV) {
   if (this!=&RV) {
      clear();
      G=RV.G;
      for (int i=0; i<RV.VolumesNo(); i++) {
          VolumeT<T>  *V=new VolumeT<T>;
          *V=RV(i);
          LV.push_back(V);
      }
   }
   return *this;
}

/* Clear ------------------------------------------------------------------- */
template <class T>
void GridVolumeT<T>::clear() {
   for (int i=0; i<VolumesNo(); i++) delete LV[i];
   LV.clear();
   G.clear();
}

/* Adapting a Reconstructing Volume to a grid ------------------------------ */
// If it is already adapted to some grid, the old adaptation is forgotten
template <class T>
void GridVolumeT<T>::adapt_to_grid(const Grid &_grid) {
   // Clear old list of volumes
   LV.clear();

   // Keep the incoming Grid at the same time the old grid is forgotten
   G=_grid;

   // Generate a volume for each subgrid
   int                        Zdim,Ydim,Xdim;
   VolumeT<T> *               Vol_aux;
   for (int i=0; i<G.GridsNo(); i++) {
      SimpleGrid & grid=G(i);
      grid.get_size(Zdim,Ydim,Xdim);
      Vol_aux=new VolumeT<T>;
      (*Vol_aux)().resize(Zdim,Ydim,Xdim);  // Using this function
                                            // after empty creation the volume
                                            // is zero-valued.
      STARTINGX((*Vol_aux)())=(int) XX(grid.lowest); // This values are already
      STARTINGY((*Vol_aux)())=(int) YY(grid.lowest); // integer although they
      STARTINGZ((*Vol_aux)())=(int) ZZ(grid.lowest); // are stored as float
      LV.push_back(Vol_aux);
   }
}

// Resize ------------------------------------------------------------------
template <class T>
void GridVolumeT<T>::resize(const matrix1D<double> &corner1,
   const matrix1D<double> &corner2) {
   VolumeT<T> *         Vol_aux;
   vector<VolumeT<T> *> LV_aux;
   
   for (int n=0; n<G.GridsNo(); n++) {
      SimpleGrid &grid=G(n);

      // Resize grid
      grid.universe2grid(corner1,grid.lowest);  grid.lowest.FLOORnD();
      grid.universe2grid(corner2,grid.highest); grid.highest.CEILnD();

      // Resize auxiliary volume
      int Zdim, Ydim, Xdim;
      grid.get_size(Zdim,Ydim,Xdim);
      Vol_aux = new VolumeT<T>;
      (*Vol_aux)().resize(Zdim,Ydim,Xdim);
      STARTINGX((*Vol_aux)())=(int) XX(grid.lowest); // This values are already
      STARTINGY((*Vol_aux)())=(int) YY(grid.lowest); // integer although they
      STARTINGZ((*Vol_aux)())=(int) ZZ(grid.lowest); // are stored as float
      
      // Copy values in common
      VolumeT<T> * origin=LV[n];
      SPEED_UP_temps;
      FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D
         (VOLMATRIX(*Vol_aux),VOLMATRIX(*origin)) {
            VOLVOXEL(*Vol_aux,k,i,j)=VOLVOXEL(*origin,k,i,j);
      }
      
      // Extract old volume and push new one
      delete LV[n];
      LV_aux.push_back(Vol_aux);
   }
   LV=LV_aux;
}

/* Arithmetic operators ---------------------------------------------------- */
#define GRIDVOLUME_BY_SCALAR(op) \
   GridVolumeT<T> result; \
   result.G = G; \
   result.LV.reserve(VolumesNo()); \
   for (int i=0; i<VolumesNo(); i++) \
      array_by_scalar((*this)(i)(),f,result(i)(),op); \
   return result;

template <class T> GridVolumeT<T> GridVolumeT<T>::operator +(T f) const
   {GRIDVOLUME_BY_SCALAR('+');}
template <class T> GridVolumeT<T> GridVolumeT<T>::operator -(T f) const
   {GRIDVOLUME_BY_SCALAR('-');}
template <class T> GridVolumeT<T> GridVolumeT<T>::operator *(T f) const
   {GRIDVOLUME_BY_SCALAR('*');}
template <class T> GridVolumeT<T> GridVolumeT<T>::operator /(T f) const
   {GRIDVOLUME_BY_SCALAR('/');}

template <class T> 
   GridVolumeT<T> operator -(T f, const GridVolumeT<T> &GV) {
   GridVolumeT<T> result;
   VolumeT<T> *Vol_aux;
   result.G = GV.G;
   result.LV.reserve(GV.VolumesNo());
   for (int i=0; i<GV.VolumesNo(); i++) {
      Vol_aux=new VolumeT<T>;
      scalar_by_array(f,GV(i)(),(*Vol_aux)(),'-');
      result.LV.push_back(Vol_aux);
   }
   return result;
}

template <class T> 
GridVolumeT<T> operator /(T f, const GridVolumeT<T> &GV) {
   GridVolumeT<T> result;
   VolumeT<T> *Vol_aux;
   result.G = GV.G;
   result.LV.reserve(GV.VolumesNo());
   for (int i=0; i<GV.VolumesNo(); i++) {
      Vol_aux=new Volume;
      scalar_by_array(f,GV(i)(),(*Vol_aux)(),'/');
      result.LV.push_back(Vol_aux);
   }
   return result;
}

#define GRIDVOL_BY_GRIDVOL(op) \
   GridVolumeT<T> result; \
   VolumeT<T> * Vol_aux; \
   \
   if (VolumesNo()!=GV.VolumesNo()) \
      REPORT_ERROR(3004,(string)"GridVolume::"+op+": Different number of subvolumes");\
   \
   result.G = G;\
   result.LV.reserve(VolumesNo());\
   \
   for (int i=0; i<VolumesNo(); i++) { \
       try { \
          Vol_aux = new VolumeT<T>; \
          array_by_array((*this)(i)(),GV(i)(),(*Vol_aux)(),op); \
          result.LV.push_back(Vol_aux); \
       } catch (Xmipp_error XE) {\
          cout << XE; \
          REPORT_ERROR(3004,(string)"GridVolume::"+op+": Different shape of volume " +\
             ItoA(i)); \
       } \
   } \
   \
   return result;

template <class T> 
GridVolumeT<T> GridVolumeT<T>::operator + (const GridVolumeT<T> &GV) _THROW
   {GRIDVOL_BY_GRIDVOL('+');}

template <class T> 
GridVolumeT<T> GridVolumeT<T>::operator - (const GridVolumeT<T> &GV) _THROW
   {GRIDVOL_BY_GRIDVOL('-');}

template <class T> 
GridVolumeT<T> GridVolumeT<T>::operator * (const GridVolumeT<T> &GV) _THROW
   {GRIDVOL_BY_GRIDVOL('*');}

template <class T> 
GridVolumeT<T> GridVolumeT<T>::operator / (const GridVolumeT<T> &GV) _THROW
   {GRIDVOL_BY_GRIDVOL('/');}

// Write a Grid volume -----------------------------------------------------
template <class T> 
void GridVolumeT<T>::write(const FileName &fn) const {
   VolumeXmippT<T>    V;
   float temp_float;
   size_t floatsize;
   const type_info &typeinfoT = typeid(T); // We need to know what kind
                                           // of variable is T
   const type_info &typeinfoD = typeid(double); 
   const type_info &typeinfoI = typeid(int); 
   
   floatsize= (size_t) sizeof(float);

   if (VolumesNo()==0) return;

   // Create the writing volume ............................................
   int Zdim=0, Ydim=0, Xdim=0;
   for (int v=0; v<VolumesNo(); v++) {
      const VolumeT<T> & this_vol=(*this)(v);
      Zdim += ZSIZE(this_vol());
      Ydim=MAX(Ydim,YSIZE(this_vol()));
      Xdim=MAX(Xdim,XSIZE(this_vol()));
   }
   
   // Check if there is enough space for the control slice
   if (Xdim*Ydim<25) Ydim=(int) CEIL(25.0f/Xdim);

   // A slice is added for control information for each subvolume
   VOLMATRIX(V).init_zeros(Zdim+VolumesNo(),Ydim,Xdim);
   
   // Write Grid volume ....................................................
   #define PACK_DOUBLE(v) \
      {jj=pos%Xdim; ii=pos/Xdim; pos++; VOLVOXEL(V,sli,ii,jj)=(T)(v);}
   #define PACK_INT(v) \
      {jj=pos%Xdim; ii=pos/Xdim; pos++; \
      temp_float = (float) (v); \
      memcpy( &(VOLVOXEL(V,sli,ii,jj)) , &temp_float, floatsize); \
      }
   
   int sli=0;
   for (int v=0; v<VolumesNo(); v++) {
      int pos, ii, jj;           // Position inside the control slice
      int k,i,j;                 // Auxiliar counters

      // Choose grid and volume
      const SimpleGrid & this_grid = grid(v);
      const VolumeT<T> & this_vol = (*this)(v);
      
      // Store Grid data ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      pos=0;

      if(typeinfoT==typeinfoD) {
         for (i=0; i<3; i++)
            for (j=0; j<3; j++) PACK_DOUBLE(MAT_ELEM(this_grid.basis  ,i,j));
         for (i=0; i<3; i++)    PACK_DOUBLE(VEC_ELEM(this_grid.lowest ,i));
         for (i=0; i<3; i++)    PACK_DOUBLE(VEC_ELEM(this_grid.highest,i));
                                PACK_DOUBLE(         this_grid.relative_size);
         for (i=0; i<3; i++)    PACK_DOUBLE(VEC_ELEM(this_grid.origin,i));
                                PACK_DOUBLE(         this_grid.R2);

         // Store volume control ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
         PACK_DOUBLE(ZSIZE(this_vol()));
         PACK_DOUBLE(YSIZE(this_vol()));
         PACK_DOUBLE(XSIZE(this_vol()));
         PACK_DOUBLE(STARTINGZ(this_vol()));
         PACK_DOUBLE(STARTINGY(this_vol()));
         PACK_DOUBLE(STARTINGX(this_vol()));
      }
      else if (typeinfoT==typeinfoI){
          // We use a trick to save the grid information in the volume
          // If the following if is true the trick can not be used   
          if((sizeof(float)!= sizeof(int)))
              REPORT_ERROR(1,
                 "GridVolume is integer and (sizeof(float)!= sizeof(int)");

         for (i=0; i<3; i++)
            for (j=0; j<3; j++) PACK_INT(MAT_ELEM(this_grid.basis  ,i,j));
         for (i=0; i<3; i++)    PACK_INT(VEC_ELEM(this_grid.lowest ,i));
         for (i=0; i<3; i++)    PACK_INT(VEC_ELEM(this_grid.highest,i));
                                PACK_INT(         this_grid.relative_size);
         for (i=0; i<3; i++)    PACK_INT(VEC_ELEM(this_grid.origin,i));
                                PACK_INT(         this_grid.R2);

         // Store volume control ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
         PACK_INT(ZSIZE(this_vol()));
         PACK_INT(YSIZE(this_vol()));
         PACK_INT(XSIZE(this_vol()));
         PACK_INT(STARTINGZ(this_vol()));
         PACK_INT(STARTINGY(this_vol()));
         PACK_INT(STARTINGX(this_vol()));      
      }
      else
         REPORT_ERROR(1,"GridVolume must be double or int\n");

      sli++;

      // Write the whole volume ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      for (k=0; k<ZSIZE(VOLMATRIX(this_vol)); k++) {
         for (i=0; i<YSIZE(VOLMATRIX(this_vol)); i++)
            for (j=0; j<XSIZE(VOLMATRIX(this_vol)); j++)
               DIRECT_VOLVOXEL(V,sli,i,j)=DIRECT_VOLVOXEL(this_vol,k,i,j);
         sli++;
      }
   }
   #undef PACK_DOUBLE
   #undef PACK_INT
   
   // Effectively write the volume .........................................
   V.write(fn);
}

// Read a Grid volume ------------------------------------------------------
//#define DEBUG
template <class T> 
void GridVolumeT<T>::read(const FileName &fn) {
   VolumeXmippT<T>    V;
   VolumeT<T>      * sV;
   SimpleGrid     sG;
   int            sli=0;

   float temp_float;
   size_t floatsize;
   const type_info &typeinfoT = typeid(T); // We need to know what kind
                                          // of variable is T
   const type_info &typeinfoD = typeid(double); 
   const type_info &typeinfoI = typeid(int); 
   
   floatsize= (size_t) sizeof(float);
   // We use a trick to save the grid information in the volume
   // If the following if is true the trick can not be used
   if(  (typeid(T) == typeid(int)) && (sizeof(float)!= sizeof(int) ) )
       {
       cout << "\nError: GridVolume is integer and\n" 
               "(sizeof(float)!= sizeof(int)\n";
       exit(0);
       }

   // Allocate memory ......................................................
   sG.basis.resize(3,3);
   sG.lowest.resize(3);
   sG.highest.resize(3);
   sG.origin.resize(3);

   // Read Reconstructing volume from file .................................
   V.read(fn);
   
   #define UNPACK_DOUBLE(v,cast) \
      {jj=pos%VOLMATRIX(V).xdim; ii=pos/VOLMATRIX(V).xdim; pos++; \
      (v)=(cast)VOLVOXEL(V,sli,ii,jj);}
   #define UNPACK_INT(v,cast) \
      {jj=pos%VOLMATRIX(V).xdim; ii=pos/VOLMATRIX(V).xdim; pos++; \
       memcpy( &temp_float, &(VOLVOXEL(V,sli,ii,jj)),floatsize);\
      (v)=(cast)temp_float;}
   
   while (sli<ZSIZE(V())) {
      int pos, ii, jj;           // Position inside the control slice
      int k,i,j;                 // Auxiliar counters
      int            Zdim, Ydim, Xdim;
      int            Zinit, Yinit, Xinit;

      // Read Grid data ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      pos=0;
      if (typeinfoT==typeinfoD) {
         for (i=0; i<3; i++)
            for (j=0; j<3; j++) UNPACK_DOUBLE(MAT_ELEM(sG.basis  ,i,j),double);
         for (i=0; i<3; i++)    UNPACK_DOUBLE(VEC_ELEM(sG.lowest ,i),int);
         for (i=0; i<3; i++)    UNPACK_DOUBLE(VEC_ELEM(sG.highest,i),int);
                                UNPACK_DOUBLE(         sG.relative_size,double);
         for (i=0; i<3; i++)    UNPACK_DOUBLE(VEC_ELEM(sG.origin,i),double);
                                UNPACK_DOUBLE(         sG.R2,double);
      } else if (typeinfoT==typeinfoI){
         // We use a trick to save the grid information in the volume
         // If the following if is true the trick can not be used   
         if ((sizeof(float)!= sizeof(int)))
            REPORT_ERROR(1,
               "GridVolume is integer and (sizeof(float)!= sizeof(int)");

         for (i=0; i<3; i++)
            for (j=0; j<3; j++) UNPACK_INT(MAT_ELEM(sG.basis  ,i,j),double);
         for (i=0; i<3; i++)    UNPACK_INT(VEC_ELEM(sG.lowest ,i),int);
         for (i=0; i<3; i++)    UNPACK_INT(VEC_ELEM(sG.highest,i),int);
                                UNPACK_INT(         sG.relative_size,double);
         for (i=0; i<3; i++)    UNPACK_INT(VEC_ELEM(sG.origin,i),double);
                                UNPACK_INT(         sG.R2,double);
      }      
      sG.inv_basis=sG.basis.inv();
      
      // Store Grid in the list of the grid volume
      G.add_grid(sG);
   
      // Read Volume Control Information ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      if (typeinfoT==typeinfoD) {
          UNPACK_DOUBLE(Zdim,int);
          UNPACK_DOUBLE(Ydim,int);
          UNPACK_DOUBLE(Xdim,int);
          UNPACK_DOUBLE(Zinit,int);
          UNPACK_DOUBLE(Yinit,int);
          UNPACK_DOUBLE(Xinit,int);
      } else if (typeinfoT==typeinfoI) {
          UNPACK_INT(Zdim,int);
          UNPACK_INT(Ydim,int);
          UNPACK_INT(Xdim,int);
          UNPACK_INT(Zinit,int);
          UNPACK_INT(Yinit,int);
          UNPACK_INT(Xinit,int);
      }

      // Set volume size and origin
      sV=new VolumeT<T>;
      VOLMATRIX(*sV).init_zeros(Zdim,Ydim,Xdim);
      STARTINGZ(VOLMATRIX(*sV))=Zinit;
      STARTINGY(VOLMATRIX(*sV))=Yinit;
      STARTINGX(VOLMATRIX(*sV))=Xinit;
      #ifdef DEBUG
         cout << "The read grid is \n" << sG;
         cout << "Volume dimensions: " << Zdim << " x " << Ydim << " x "
              << Xdim << endl;
         cout << "Volume init: " << Zinit << " x " << Yinit << " x "
              << Xinit << endl;
      #endif
      sli++;

      // Read volume ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      for (k=0; k<ZSIZE(VOLMATRIX(*sV)); k++) {
         for (i=0; i<YSIZE(VOLMATRIX(*sV)); i++)
            for (j=0; j<XSIZE(VOLMATRIX(*sV)); j++) {
               #ifdef DEBUG
                  cout << "Reading from file position (" << sli << "," << i
                       << "," << j << ") to subvolume position ("
                       << k << "," << i << "," << j << ")\n";
               #endif
               DIRECT_VOLVOXEL(*sV,k,i,j)=DIRECT_VOLVOXEL(V,sli,i,j);
            }
         sli++;
      }
      
      // Store volume in the list
      LV.push_back(sV);
   }
   #undef UNPACK_DOUBLE
   #undef UNPACK_INT
}
#undef DEBUG

// Show a grid volume ------------------------------------------------------
template <class T>
ostream& operator << (ostream &o, const GridVolumeT<T> &GV) {
   o << "Grid Volume -----------\n";
   o << GV.G;
   o << "Number of volumes= " << GV.VolumesNo() << endl;
   for (int i=0; i<GV.VolumesNo(); i++) {
      o << "Volume " << i << "------------" << endl;
      o << GV(i)();
   }
   return o;
}

/* Instantiation ----------------------------------------------------------- */
template <class T>
   void instantiate_GridVolume (GridVolumeT<T> Gx) {
   GridVolumeT<double> GVa;// Empty constructor
   GridVolumeT<int> GVb;   // Empty constructor
   GridVolumeT<T> GVc(Gx); // copy constructor
   VolumeT<T>    Vint;
   Grid my_grid;
   T number;
   
   GridVolumeT<T> Gxx(my_grid);// Create using a grid as pattern and basis.
   GVc = Gx;//Assignement;
   Gx.adapt_to_grid(my_grid);
   Gx.init_zeros();
   Gx.clear();
   Gx(0);// Constant access to a volume

   Gx(0)(); 
   Gx.grid(0);// Constant access to a simple grid
   Gx.grid();// Constant access to the whole grid
   Gx+number;Gx-number;Gx*number;Gx/number;
   GVc+Gx; GVc-Gx;GVc*Gx;GVc/Gx;
   Gx.read("kk"); Gx.write("kk");
   cout << Gx;
   Gx.VolumesNo();
   Gx = Gx + (T) 1;      
   matrix1D<double> v; Gx.resize(v,v);
   Gx.resize(GVa);
   Gx.resize(GVb);
}

void instantiateGridVolume() {
   GridVolumeT<double>  GV1; instantiate_GridVolume(GV1);
   GridVolumeT<int>    GV2; instantiate_GridVolume(GV2);
}
