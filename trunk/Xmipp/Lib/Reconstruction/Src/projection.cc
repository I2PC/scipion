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

#include "../projection.hh"

// Some aliases
#define x0   STARTINGX(IMGMATRIX(proj))
#define xF   FINISHINGX(IMGMATRIX(proj))
#define y0   STARTINGY(IMGMATRIX(proj))
#define yF   FINISHINGY(IMGMATRIX(proj))
#define xDim XSIZE(IMGMATRIX(proj))
#define yDim YSIZE(IMGMATRIX(proj))

// Projections from single particles #######################################
// Projections from blob volumes ===========================================

/* Projection of a simple volume ------------------------------------------- */
// The projection is not precleaned (set to 0) before projecting and its
// angles are supposed to be already written (and all Euler matrices
// precalculated
// The projection plane is supposed to pass through the Universal coordinate
// origin

/* Time measures ...........................................................
   This function has been optimized in time, several approaches have been
   tried and here are the improvements in time
   Case: Base (unit time measure)
         The points within the footprint are calculated inside the
         inner loop instead of computing the foot coordinate for the
         first corner and knowing the sampling rate in the oversampled
         image go on the next points. The image access was done directly
         (physical access).

   Notation simplification:
   Case: VOLVOXEL(V,k,i,j) ----> V(k,i,j):   + 38% (Inacceptable)
   Case: DIRECT_IMGPIXEL   ----> IMGPIXEL:   +  5% (Acceptable)
   
   Algorithmic changes:
   Case: project each blob position ---->    +325% (Inacceptable)
         get rid of all beginZ, beginY, prjX, ... and project each
         blob position onto the projection plane
   Case: footprint coordinates outside ----> - 33% (Great!!)
         the footprint coordinates are computed outside the inner
         loop, but you have to use the sampling rate instead to
         move over the oversampled image.
*/

//#define DEBUG_LITTLE
//#define DEBUG

template <class T>
void project_SimpleGrid(VolumeT<T> &vol, const SimpleGrid &grid,
   const struct blobtype &blob,
   const ImageOver &footprint, const ImageOver &footprint2,
   Projection &proj, Projection &norm_proj, int FORW, int eq_mode,
   const VolumeT<int> *VNeq, matrix2D<double> *M,
   VolumeT<T> *vol_var=NULL, double ray_length=-1.0) {
   matrix1D<double> zero(3);                // Origin (0,0,0) 
   matrix1D<double> prjPix(3);              // Position of the pixel within the
                                            // projection (used in VSSNR)
   matrix1D<double> prjX(3);                // Coordinate: Projection of the
   matrix1D<double> prjY(3);                // 3 grid vectors
   matrix1D<double> prjZ(3);
   matrix1D<double> prjOrigin(3);           // Coordinate: Where in the 2D
                                            // projection plane the origin of
                                            // the grid projects

   matrix1D<double> actprj(3);              // Coord: Actual position inside
                                            // the projection plane
   matrix1D<double> beginZ(3);              // Coord: Plane coordinates of the
                                            // projection of the 3D point
                                            // (z0,YY(lowest),XX(lowest))
   matrix1D<double> beginY(3);              // Coord: Plane coordinates of the
                                            // projection of the 3D point
                                            // (z0,y0,XX(lowest))
   double XX_footprint_size;                // The footprint is supposed
   double YY_footprint_size;                // to be defined between
                                            // (-vmax,+vmax) in the Y axis,
                                            // and (-umax,+umax) in the X axis
                                            // This footprint size is the
                                            // R2 vector (umax,vmax).
                                       
   int XX_corner2, XX_corner1;              // Coord: Corners of the
   int YY_corner2, YY_corner1;              // footprint when it is projected
                                            // onto the projection plane
   int           foot_V1, foot_U1;          // Img Coord: coordinate (in
                                            // an image fashion, not in an
                                            // oversampled image fashion)
                                            // inside the blobprint of the
                                            // corner1
   int           foot_V, foot_U;            // Img Coord: coordinate
                                            // corresponding to the blobprint
                                            // point which matches with this
                                            // pixel position
   int           Vsampling, Usampling;      // Sampling rate in Y and X
                                            // directions respectively
                                            // inside the blobprint
   double        vol_corr;                  // Correction for a volum element
   int           N_eq;                      // Number of equations in which
                                            // a blob is involved
   int           i,j,k;                     // volume element indexes

   // Prepare system matrix for printing ...................................
   if (M!=NULL)
      M->init_zeros(YSIZE(proj())*XSIZE(proj()),grid.get_number_of_samples());

   // Project grid axis ....................................................
   // These vectors ((1,0,0),(0,1,0),...) are referred to the grid
   // coord. system and the result is a 2D vector in the image plane
   // The unit size within the image plane is 1, ie, the same as for
   // the universal grid.
   // actprj is reused with a different purpose
   VECTOR_R3(actprj,1,0,0); grid.Gdir_project_to_plane(actprj,proj.euler,prjX);
   VECTOR_R3(actprj,0,1,0); grid.Gdir_project_to_plane(actprj,proj.euler,prjY);
   VECTOR_R3(actprj,0,0,1); grid.Gdir_project_to_plane(actprj,proj.euler,prjZ);
   
   // This time the origin of the grid is in the universal coord system
   Uproject_to_plane(grid.origin, proj.euler, prjOrigin);

   // Footprint size .......................................................
   // The following vectors are integer valued vectors, but they are
   // stored as real ones to make easier operations with other vectors
   XX_footprint_size = footprint.umax();
   YY_footprint_size = footprint.vmax();
   Usampling         = footprint.Ustep();
   Vsampling         = footprint.Vstep();

   // Project the whole grid ...............................................
   // Corner of the plane defined by Z. These coordinates try to be within
   // the valid indexes of the projection (defined between STARTING and
   // FINISHING values, but it could be that they may lie outside.
   // These coordinates need not to be integer, in general, they are
   // real vectors.
   // The vectors returned by the projection routines are R3 but we
   // are only interested in their first 2 components, ie, in the
   // in-plane components
   beginZ=XX(grid.lowest)*prjX + YY(grid.lowest)*prjY + ZZ(grid.lowest)*prjZ
          + prjOrigin;

   // This type conversion gives more speed
   int ZZ_lowest =(int) ZZ(grid.lowest);
   int YY_lowest =(int) YY(grid.lowest);
   int XX_lowest =(int) XX(grid.lowest);
   int ZZ_highest=(int) ZZ(grid.highest);
   int YY_highest=(int) YY(grid.highest);
   int XX_highest=(int) XX(grid.highest);

   #ifdef DEBUG_LITTLE
   int condition;
   condition = FORW && XX(grid.origin)==0;
   if (condition) {
      cout << "Equation mode " << eq_mode << endl;
      cout << "Footprint size " << YY_footprint_size << "x"
           << XX_footprint_size << endl;
      cout << "rot=" << proj.Phi() << " tilt=" << proj.Theta()
           << " psi=" << proj.Psi() << endl;
      cout << "Euler matrix " << proj.euler;
      cout << grid;
      cout << "image limits (" << x0 << "," << y0 << ") (" << xF << ","
           << yF << ")\n";
      cout << "prjX           " << prjX.transpose()      << endl;
      cout << "prjY           " << prjY.transpose()      << endl;
      cout << "prjZ           " << prjZ.transpose()      << endl;
      cout << "prjOrigin      " << prjOrigin.transpose() << endl;
      cout << "beginZ(coord)  " << beginZ.transpose()    << endl;
      cout << "lowest         " << XX_lowest  << " " << YY_lowest
           << " " << XX_lowest  << endl;
      cout << "highest        " << XX_highest << " " << YY_highest
           << " " << XX_highest << endl;
      cout << "Stats of input blob volume "; vol().print_stats(); cout << endl;
      cout.flush();
   }
   #endif

   matrix1D<double> grid_index(3), univ_position(3);
   int number_of_blob=0;
   for (k=ZZ_lowest; k<=ZZ_highest; k++) {
      // Corner of the row defined by Y
      beginY=beginZ;
      for (i=YY_lowest; i<=YY_highest; i++) {
         // First point in the row
	 actprj=beginY;
	 for (j=XX_lowest; j<=XX_highest; j++) {
	    VECTOR_R3(grid_index,j,i,k);
	    grid.grid2universe(grid_index,univ_position);
            
            // Ray length interesting
            bool ray_length_interesting=true;
            if (ray_length!=-1) {
               ray_length_interesting=
                  ABS(point_plane_distance_3D(univ_position, zero,
                      proj.direction))<=ray_length;
            }
            
            
      	    if (grid.is_interesting(univ_position) &&
                ray_length_interesting) {
               // Be careful that you cannot skip any blob, although its
               // value be 0, because it is useful for norm_proj
               #ifdef DEBUG
               condition = FORW && XX(grid.origin)==0 &&
                             ((j==0 && i==0 && k==0));
               if (condition) {
        	  printf("\nProjecting grid coord (%d,%d,%d) ",j,i,k);
        	  cout << "Vol there = " << VOLVOXEL(vol,k,i,j) << endl;
		  printf(" 3D universal position (%f,%f,%f) \n",
		     ZZ(univ_position),YY(univ_position),XX(univ_position));
        	  cout << " Center of the blob proj " << actprj.transpose()  << endl;
		  matrix1D<double> aux;
		  Uproject_to_plane(univ_position, proj.euler, aux);
		  cout << " Center of the blob proj (more accurate) " << aux.transpose() << endl;
               }
               #endif

               // Search for integer corners for this blob
               XX_corner1=CEIL (MAX(x0,XX(actprj)-XX_footprint_size));
               YY_corner1=CEIL (MAX(y0,YY(actprj)-YY_footprint_size));
               XX_corner2=FLOOR(MIN(xF,XX(actprj)+XX_footprint_size));
               YY_corner2=FLOOR(MIN(yF,YY(actprj)+YY_footprint_size));

               #ifdef DEBUG
               if (condition) {
        	  cout << "Clipped and rounded Corner 1 " << XX_corner1
                       << " " << YY_corner1 << " " << endl;
        	  cout << "Clipped and rounded Corner 2 " << XX_corner2
                       << " " << YY_corner2 << " " << endl;
               }
               #endif

               // If the two corners haven't collapsed into a line or a point
               // This happens if the blob falls outside the projection plane
               if (XX_corner1<XX_corner2 && YY_corner1<YY_corner2) {

        	  // Compute the index of the blob for corner1
		  OVER2IMG(footprint,(double)YY_corner1-YY(actprj),
                     (double)XX_corner1-XX(actprj), foot_V1, foot_U1);

        	  if (!FORW) vol_corr=0;

		  // Effectively project this blob
        	  N_eq=(YY_corner2-YY_corner1+1)*(XX_corner2-XX_corner1+1);
        	  foot_V=foot_V1;
		  for (int y=YY_corner1; y<=YY_corner2; y++) {
                     foot_U=foot_U1;
	             for (int x=XX_corner1; x<=XX_corner2; x++) {
                	#ifdef DEBUG
                	if (condition) {
                           cout << "Position in projection (" << x << ","
                        	<< y << ") " << "in footprint ("
                        	<< foot_U << "," << foot_V << ")";
                           double y,x;
                           IMG2OVER(footprint,foot_V,foot_U,y,x);
                           cout << " (d= " << sqrt(y*y+x*x) << ") ";
                           fflush(stdout);
                           cout << "=" << IMGPIXEL(footprint,foot_V,foot_U)
                        	<< " , " << IMGPIXEL(footprint2,foot_V,foot_U);
                	}
                	#endif
                        double a, a2;
                        // Check if volumetric interpolation (i.e., SSNR)
                        if (ray_length==blob.radius) {
                           // This is the VSSNR case
                           // Get the pixel position in the universal coordinate
                           // system
                           SPEED_UP_temps;
                           VECTOR_R3(prjPix,x,y,0);
                           M3x3_BY_V3x1(prjPix,proj.eulert,prjPix);
                           // Get the distance between the blob center
                           // and the pixel
                           V3_MINUS_V3(prjPix,univ_position,prjPix);
                           double r=prjPix.module();
                           a=blob_val(r,blob);
                           a2=a*a;
                        } else {
                           // This is normal reconstruction from projections
                           a =IMGPIXEL(footprint,foot_V,foot_U);
                           a2=IMGPIXEL(footprint2,foot_V,foot_U);
                        }
	        	if (FORW) {
                           switch(eq_mode) {
			      case CAVARTK:
                              case ARTK:
                        	 IMGPIXEL(proj,y,x) += VOLVOXEL(vol,k,i,j)*a;
                        	 IMGPIXEL(norm_proj,y,x) += a2;
				 if (M!=NULL) {
                                    int py, px;
				    proj().logical2physical(y,x,py,px);
				    int number_of_pixel=py*XSIZE(proj())+px;
				    (*M)(number_of_pixel,number_of_blob)=a;
				 }
                        	 break;
                              case CAVK:
                        	 IMGPIXEL(proj,y,x) += VOLVOXEL(vol,k,i,j)*a;
                        	 IMGPIXEL(norm_proj,y,x) += a2*N_eq;
                        	 break;
                              case COUNT_EQ:
                        	 VOLVOXEL(vol,k,i,j)++;
                        	 break;
                              case CAV:
                        	 IMGPIXEL(proj,y,x) += VOLVOXEL(vol,k,i,j)*a;
                        	 IMGPIXEL(norm_proj,y,x) += a2 *
                                    VOLVOXEL(*VNeq,k,i,j);
                        	 break;
                           }

                           #ifdef DEBUG
                              if (condition) {
                        	 cout << " proj= " << IMGPIXEL(proj,y,x)
                                      << " norm_proj=" << IMGPIXEL(norm_proj,y,x)  << endl;
                        	 cout << endl; cout.flush();
                              }
                           #endif
	        	} else {
	                   vol_corr += IMGPIXEL(norm_proj,y,x)*a;
                           #ifdef DEBUG
                           if (condition) {
                              cout << " corr_img= " << IMGPIXEL(norm_proj,y,x)
                                   << " correction=" << vol_corr << endl;
                              cout.flush();
                           }
                           #endif
	        	}

                	// Prepare for next operation
                	foot_U += Usampling;
	             }
                     foot_V += Vsampling;
        	  } // Project this blob

		  if (!FORW) {
                     T correction;
                     switch (eq_mode) {
                	case ARTK: correction=(T) vol_corr;
                           break;
                	case CAVARTK: correction=(T) (vol_corr/N_eq);
                           break;
                	case CAVK:
                	case CAV:  correction=(T) vol_corr;
                           break;
                     }
                     VOLVOXEL(vol,k,i,j)+=correction;
                     if (vol_var!=NULL) {
                        bool interesting=grid.R2==-1 ||
                           (sqrt(XX(univ_position)*XX(univ_position)+
                           YY(univ_position)*YY(univ_position)+
                           ZZ(univ_position)*ZZ(univ_position))+5<
                           sqrt(grid.R2));
                        if (interesting || true)
                           VOLVOXEL(*vol_var,k,i,j)=(T)
                              (0.99*VOLVOXEL(*vol_var,k,i,j)+
                               0.01*ABS(correction));
                               //0.005*ABS(correction));
                     }

                     #ifdef DEBUG
                     if (condition) {
                	printf("\nFinal value at (%d,%d,%d) ",j,i,k);
                	cout << " = " << VOLVOXEL(vol,k,i,j) << endl;
                	cout.flush();
                     }
                     #endif
        	  }
               } // If not collapsed
	       number_of_blob++;
      	    } // If interesting

            // Prepare for next iteration
            V2_PLUS_V2(actprj, actprj, prjX);
      	 }
         V2_PLUS_V2(beginY, beginY, prjY);
      }
      V2_PLUS_V2(beginZ, beginZ, prjZ);
   }
}

#undef DEBUG
#undef DEBUG_LITTLE

/* Project a Grid Volume --------------------------------------------------- */
//#define DEBUG
//#define DEBUG_LITTLE
template <class T>
void project_Volume(
   GridVolumeT<T> &vol,                  // Volume
   const struct blobtype &blob,          // Blob
   const ImageOver  &footprint,          // Footprint of a blob         
   const ImageOver  &footprint2,         // Used for error computation
   Projection       &proj,               // Projection
   Projection       &norm_proj,          // Projection of a unitary volume
   int              Ydim,                // Dimensions of the projection
   int              Xdim,
   double rot, double tilt, double psi,  // Euler angles
   int              FORW,                // 1 if we are projecting a volume
                                         //   norm_proj is calculated
                                         // 0 if we are backprojecting
                                         //   norm_proj must be valid
   int              eq_mode,             // ARTK, CAVARTK, CAVK or CAV
   GridVolumeT<int> *GVNeq,              // Number of equations per blob
   matrix2D<double> *M,                  // System matrix
   GridVolumeT<T>   *vol_var,            // Volume to keep track of the variance
   double            ray_length)         // Ray length of the projection
{
   // If projecting forward initialise projections
   if (FORW) {
      proj.reset(Ydim,Xdim);
      proj.set_angles(rot,tilt,psi);
      norm_proj().resize(proj());
      norm_proj().init_zeros();
   }
   
   #ifdef DEBUG_LITTLE
      if (FORW) {
	 cout << "Number of volumes: " << vol.VolumesNo() << endl
              << "YdimxXdim: " << Ydim << "x" << Xdim << endl;
         for (int i=0; i<vol.VolumesNo(); i++)
	     cout << "Volume " << i << endl << vol.grid(i) << endl;
      }
      
   #endif

   // Project each subvolume
   for (int i=0; i<vol.VolumesNo(); i++) {
      VolumeT<int> *VNeq;
      VolumeT<T>   *Vvar;
      if (GVNeq!=NULL)   VNeq=&((*GVNeq  )(i)); else VNeq=NULL;
      if (vol_var!=NULL) Vvar=&((*vol_var)(i)); else Vvar=NULL;
      project_SimpleGrid(vol(i),vol.grid(i),blob,footprint,footprint2,
         proj, norm_proj, FORW, eq_mode, VNeq, M, Vvar, ray_length);
   #ifdef DEBUG
      ImageXmipp save; save=norm_proj;
      if (FORW) save.write((string)"PPPnorm_FORW"+(char)(48+i));
      else      save.write((string)"PPPnorm_BACK"+(char)(48+i));
   #endif
   }
}
#undef DEBUG_LITTLE
#undef DEBUG

// Projection from a voxel volume ==========================================
/* Project a voxel volume -------------------------------------------------- */
//#define DEBUG
void project_Volume(matrix3D<double> &V, Projection &P, int Ydim, int Xdim,
   double rot, double tilt, double psi) _THROW {
   SPEED_UP_temps;

   // Initialise projection
   P.reset(Ydim, Xdim);
   P.set_angles(rot,tilt,psi);
  
   // Compute the distance for this line crossing one voxel
   int x_0=STARTINGX(V), x_F=FINISHINGX(V);
   int y_0=STARTINGY(V), y_F=FINISHINGY(V);
   int z_0=STARTINGZ(V), z_F=FINISHINGZ(V);

   // Distances in X and Y between the center of the projection pixel begin
   // computed and each computed ray
   double step= 1.0/3.0;
    
   // Avoids divisions by zero and allows orthogonal rays computation
   if(XX(P.direction)==0) XX(P.direction)=XMIPP_EQUAL_ACCURACY;
   if(YY(P.direction)==0) YY(P.direction)=XMIPP_EQUAL_ACCURACY;
   if(ZZ(P.direction)==0) ZZ(P.direction)=XMIPP_EQUAL_ACCURACY;

   // Some precalculated variables
   int x_sign = SGN(XX(P.direction));
   int y_sign = SGN(YY(P.direction));
   int z_sign = SGN(ZZ(P.direction));
   double half_x_sign = 0.5 * x_sign;
   double half_y_sign = 0.5 * y_sign;
   double half_z_sign = 0.5 * z_sign;
   
   FOR_ALL_ELEMENTS_IN_MATRIX2D(P()) {
      matrix1D<double> r_p(3); // r_p are the coordinates of the 
                               // pixel being projected in the
                               // coordinate system attached to the
                               // projection
      matrix1D<double> p1(3);  // coordinates of the pixel in the
      			       // universal space
      double ray_sum=0.0;      // Line integral value
      
      // Computes 4 different rays for each pixel.
      for(int rays_per_pixel = 0; rays_per_pixel<4; rays_per_pixel++)
      { 
	// universal coordinate system                        
	switch(rays_per_pixel){
		case 0:	VECTOR_R3(r_p,j-step,i-step,0);
			break;
		case 1: VECTOR_R3(r_p,j-step,i+step,0);
			break;
		case 2: VECTOR_R3(r_p,j+step,i-step,0);
			break;
		case 3: VECTOR_R3(r_p,j+step,i+step,0);
			break;
	}
			
	// Express r_p in the universal coordinate system
	M3x3_BY_V3x1(p1,P.eulert,r_p);
	
	// Compute the minimum and maximum alpha for the ray
	// intersecting the given volume
	double alpha_xmin=(x_0-0.5-XX(p1))/XX(P.direction);
	double alpha_xmax=(x_F+0.5-XX(p1))/XX(P.direction);
	double alpha_ymin=(y_0-0.5-YY(p1))/YY(P.direction);
	double alpha_ymax=(y_F+0.5-YY(p1))/YY(P.direction);
	double alpha_zmin=(z_0-0.5-ZZ(p1))/ZZ(P.direction);
	double alpha_zmax=(z_F+0.5-ZZ(p1))/ZZ(P.direction);
	
	double alpha_min =MAX(MIN(alpha_xmin,alpha_xmax),
                              MIN(alpha_ymin,alpha_ymax));
               alpha_min =MAX(alpha_min,MIN(alpha_zmin,alpha_zmax));
	double alpha_max =MIN(MAX(alpha_xmin,alpha_xmax),
                              MAX(alpha_ymin,alpha_ymax));
               alpha_max =MIN(alpha_max,MAX(alpha_zmin,alpha_zmax));
	if (alpha_max-alpha_min<XMIPP_EQUAL_ACCURACY)	continue;
	
	#ifdef DEBUG
	   cout << "Pixel:  " << r_p.transpose() << endl
        	<< "Univ:   " << p1.transpose() << endl
        	<< "Dir:    " << P.direction.transpose() << endl
        	<< "Alpha x:" << alpha_xmin << " " << alpha_xmax << endl
        	<< "        " << (p1+alpha_xmin*P.direction).transpose() << endl
        	<< "        " << (p1+alpha_xmax*P.direction).transpose() << endl
        	<< "Alpha y:" << alpha_ymin << " " << alpha_ymax << endl
        	<< "        " << (p1+alpha_ymin*P.direction).transpose() << endl
        	<< "        " << (p1+alpha_ymax*P.direction).transpose() << endl
        	<< "Alpha z:" << alpha_zmin << " " << alpha_zmax << endl
        	<< "        " << (p1+alpha_zmin*P.direction).transpose() << endl
        	<< "        " << (p1+alpha_zmax*P.direction).transpose() << endl
        	<< "alpha  :" << alpha_min  << " " << alpha_max  << endl
        	<< endl;
	#endif

	// Compute the first point in the volume intersecting the ray
	matrix1D<double>  v(3);
	matrix1D<int>    idx(3);
	V3_BY_CT(v,P.direction,alpha_min);
	V3_PLUS_V3(v,p1,v);

	// Compute the index of the first voxel
	XX(idx)=CLIP(ROUND(XX(v)),x_0,x_F);
	YY(idx)=CLIP(ROUND(YY(v)),y_0,y_F);
	ZZ(idx)=CLIP(ROUND(ZZ(v)),z_0,z_F);

	
	#ifdef DEBUG
           cout << "First voxel: " << v.transpose() << endl;
           cout << "   First index: " << idx.transpose() << endl;
           cout << "   Alpha_min: " << alpha_min << endl;
	#endif
   
	// Follow the ray
	double alpha=alpha_min;
	do {  
	   #ifdef DEBUG
	   	cout << " \n\nCurrent Value: " << V(ZZ(idx),YY(idx),XX(idx)) << endl;
	   #endif
	   
	   double alpha_x = (XX(idx)+half_x_sign-XX(p1))/XX(P.direction);
           double alpha_y = (YY(idx)+half_y_sign-YY(p1))/YY(P.direction);
           double alpha_z = (ZZ(idx)+half_z_sign-ZZ(p1))/ZZ(P.direction);

	   // Which dimension will ray move next step into?, it isn't neccesary to be only
	   // one. 
	   double diffx = ABS(alpha-alpha_x);
	   double diffy = ABS(alpha-alpha_y);
	   double diffz = ABS(alpha-alpha_z);
	  
	   double diff_alpha = MIN(MIN(diffx,diffy),diffz);

	   ray_sum+= diff_alpha * V(ZZ(idx),YY(idx),XX(idx)); 
	
	   if(ABS(diff_alpha-diffx) <=XMIPP_EQUAL_ACCURACY){
		  alpha = alpha_x; XX(idx)+=x_sign;
	   }
	   if(ABS(diff_alpha-diffy) <=XMIPP_EQUAL_ACCURACY){
		  alpha = alpha_y; YY(idx)+=y_sign;
	   }
	   if(ABS(diff_alpha-diffz) <=XMIPP_EQUAL_ACCURACY){
		  alpha = alpha_z; ZZ(idx)+=z_sign;
	   }
    
	   #ifdef DEBUG
              cout << "Alpha x,y,z: " << alpha_x << " " << alpha_y
                   << " " << alpha_z << " ---> " << alpha << endl;
           
	      XX(v)+=diff_alpha*XX(P.direction);
              YY(v)+=diff_alpha*YY(P.direction);
              ZZ(v)+=diff_alpha*ZZ(P.direction);
    
	      cout << "    Next entry point: " << v.transpose() << endl
                   << "    Index: " << idx.transpose() << endl
                   << "    diff_alpha: " << diff_alpha << endl
                   << "    ray_sum: " << ray_sum << endl
		   << "    Alfa tot: " << alpha << "alpha_max: " << alpha_max <<
		   endl;
		
           #endif
	} while ((alpha_max-alpha)>XMIPP_EQUAL_ACCURACY);
      } // for

      P(i,j)=ray_sum*0.25;
      #ifdef DEBUG
      	cout << "Assigning P(" << i << "," << j << ")=" << ray_sum << endl;
      #endif
   }
}
#undef DEBUG

// Projection from a voxel volume ==========================================
/* Project a voxel volume -------------------------------------------------- */
//#define DEBUG
void singleWBP(matrix3D<double> &V, Projection &P) _THROW {
   SPEED_UP_temps;

   // Compute the distance for this line crossing one voxel
   int x_0=STARTINGX(V), x_F=FINISHINGX(V);
   int y_0=STARTINGY(V), y_F=FINISHINGY(V);
   int z_0=STARTINGZ(V), z_F=FINISHINGZ(V);

   // Distances in X and Y between the center of the projection pixel begin
   // computed and each computed ray
   double step= 1.0/3.0;
    
   // Avoids divisions by zero and allows orthogonal rays computation
   if(XX(P.direction)==0) XX(P.direction)=XMIPP_EQUAL_ACCURACY;
   if(YY(P.direction)==0) YY(P.direction)=XMIPP_EQUAL_ACCURACY;
   if(ZZ(P.direction)==0) ZZ(P.direction)=XMIPP_EQUAL_ACCURACY;

   // Some precalculated variables
   int x_sign = SGN(XX(P.direction));
   int y_sign = SGN(YY(P.direction));
   int z_sign = SGN(ZZ(P.direction));
   double half_x_sign = 0.5 * x_sign;
   double half_y_sign = 0.5 * y_sign;
   double half_z_sign = 0.5 * z_sign;
   
   FOR_ALL_ELEMENTS_IN_MATRIX2D(P()) {
      matrix1D<double> r_p(3); // r_p are the coordinates of the 
                               // pixel being projected in the
                               // coordinate system attached to the
                               // projection
      matrix1D<double> p1(3);  // coordinates of the pixel in the
      			       // universal space
      double ray_sum=0.0;      // Line integral value
      
      // Computes 4 different rays for each pixel.
/*      for(int rays_per_pixel = 0; rays_per_pixel<4; rays_per_pixel++)
      { 
	// universal coordinate system                        
	switch(rays_per_pixel){
		case 0:	VECTOR_R3(r_p,j-step,i-step,0);
			break;
		case 1: VECTOR_R3(r_p,j-step,i+step,0);
			break;
		case 2: VECTOR_R3(r_p,j+step,i-step,0);
			break;
		case 3: VECTOR_R3(r_p,j+step,i+step,0);
			break;
	}*/
	
	//***********************
	VECTOR_R3(r_p,j,i,0);		
	
	// Express r_p in the universal coordinate system
	M3x3_BY_V3x1(p1,P.eulert,r_p);
	
	// Compute the minimum and maximum alpha for the ray
	// intersecting the given volume
	double alpha_xmin=(x_0-0.5-XX(p1))/XX(P.direction);
	double alpha_xmax=(x_F+0.5-XX(p1))/XX(P.direction);
	double alpha_ymin=(y_0-0.5-YY(p1))/YY(P.direction);
	double alpha_ymax=(y_F+0.5-YY(p1))/YY(P.direction);
	double alpha_zmin=(z_0-0.5-ZZ(p1))/ZZ(P.direction);
	double alpha_zmax=(z_F+0.5-ZZ(p1))/ZZ(P.direction);
	
	double alpha_min =MAX(MIN(alpha_xmin,alpha_xmax),
                              MIN(alpha_ymin,alpha_ymax));
               alpha_min =MAX(alpha_min,MIN(alpha_zmin,alpha_zmax));
	double alpha_max =MIN(MAX(alpha_xmin,alpha_xmax),
                              MAX(alpha_ymin,alpha_ymax));
               alpha_max =MIN(alpha_max,MAX(alpha_zmin,alpha_zmax));
	if (alpha_max-alpha_min<XMIPP_EQUAL_ACCURACY)	continue;
	
	// Compute the first point in the volume intersecting the ray
	matrix1D<double>  v(3);
	matrix1D<int>    idx(3);
	V3_BY_CT(v,P.direction,alpha_min);
	V3_PLUS_V3(v,p1,v);

	// Compute the index of the first voxel
	XX(idx)=CLIP(ROUND(XX(v)),x_0,x_F);
	YY(idx)=CLIP(ROUND(YY(v)),y_0,y_F);
	ZZ(idx)=CLIP(ROUND(ZZ(v)),z_0,z_F);

	
	#ifdef DEBUG
           cout << "First voxel: " << v.transpose() << endl;
           cout << "   First index: " << idx.transpose() << endl;
           cout << "   Alpha_min: " << alpha_min << endl;
	#endif
   
	// Follow the ray
	double alpha=alpha_min;
	do {  
	   #ifdef DEBUG
	   	cout << " \n\nCurrent Value: " << V(ZZ(idx),YY(idx),XX(idx)) << endl;
	   #endif
	   
	   double alpha_x = (XX(idx)+half_x_sign-XX(p1))/XX(P.direction);
           double alpha_y = (YY(idx)+half_y_sign-YY(p1))/YY(P.direction);
           double alpha_z = (ZZ(idx)+half_z_sign-ZZ(p1))/ZZ(P.direction);

	   // Which dimension will ray move next step into?, it isn't neccesary to be only
	   // one. 
	   double diffx = ABS(alpha-alpha_x);
	   double diffy = ABS(alpha-alpha_y);
	   double diffz = ABS(alpha-alpha_z);
	  
	   double diff_alpha = MIN(MIN(diffx,diffy),diffz);

	   V(ZZ(idx), YY(idx), XX(idx)) += diff_alpha * P(i,j);
	
	   if(ABS(diff_alpha-diffx) <=XMIPP_EQUAL_ACCURACY){
		  alpha = alpha_x; XX(idx)+=x_sign;
	   }
	   if(ABS(diff_alpha-diffy) <=XMIPP_EQUAL_ACCURACY){
		  alpha = alpha_y; YY(idx)+=y_sign;
	   }
	   if(ABS(diff_alpha-diffz) <=XMIPP_EQUAL_ACCURACY){
		  alpha = alpha_z; ZZ(idx)+=z_sign;
	   }
    
	} while ((alpha_max-alpha)>XMIPP_EQUAL_ACCURACY);
  //    } // for

       #ifdef DEBUG
      	cout << "Assigning P(" << i << "," << j << ")=" << ray_sum << endl;
      #endif
   }
}
#undef DEBUG

#ifdef NEVER_DEFINED
// Project using spider
void project_Volume(VolumeXmipp &V, Projection &P, int Ydim, int Xdim,
   double rot, double tilt, double psi) _THROW {
   // Initialise projection
   P.reset(Ydim, Xdim);
   P.set_angles(rot,tilt,psi);

   // Project all voxels
   FileName fn_vol=V.name();
   if (fn_vol=="") V.write("superfeo.kkk");
   else system(((string) "ln -s "+fn_vol+" superfeo.kkk").c_str());

   ofstream spider_batch;
   spider_batch.open("b01.kkk");
   if (!spider_batch)
      REPORT_ERROR(1,"Project:: Cannot open file for Spider batch");
   spider_batch
      << "pj 3e\n"
      << "superfeo" << endl
      << Xdim << " " << Ydim << endl
      << "superfeo2\n"
      << rot << endl
      << tilt << endl
      << psi << endl
      << "en\n"
   ;
   spider_batch.close();
   char *spider_prog=getenv("SPIDER");
   if (spider_prog==NULL)
      REPORT_ERROR(1,"Project:: The environment variable SPIDER is not set");
   system(((string)spider_prog+" kkk b01").c_str());
   system("rm LOG*.kkk results*.kkk* b01*.kkk");
   
   P.read("superfeo2.kkk");
   system("rm superfeo*.kkk");
}
#endif
// Projections from crystals particles #####################################
// The projection is not precleaned (set to 0) before projecting and its
// angles are supposed to be already written (and all Euler matrices
// precalculated
// The projection plane is supposed to pass through the Universal coordinate
// origin

/* Algorithm 
Compute Eg, proj(ai), proj(bi) and A
Compute prjX, prjY, prjZ and prjO which are the projections of the origin
   and grid axes
Compute the blob footprint size in the deformed image space
   (=deffootprintsize)
For each blob in the grid
   if it is within the unit cell mask
      // compute actual deformed projection position
      defactprj=A*(k*projZ+i*projY+j*projX+projO)
      // compute corners in the deformed image
      corner1=defactprj-deffootprintsize;
      corner2=defactprj+ deffootprintsize;

      for each point (r) in the deformed projection space between (corner1, corner2)
         compute position (rc) in the undeformed projection space
         compute position within footprint corresponding to rc (=foot)
         if it is inside footprint
            rw=wrap(r);
            update projection at rw with the data at foot  
*/

//#define DEBUG_LITTLE
//#define DEBUG
//#define DEBUG_INTERMIDIATE
#define wrap_as_Crystal(x,y,xw,yw)  \
   xw=(int) intWRAP(x,x0,xF); \
   yw=(int) intWRAP(y,y0,yF);

void project_Crystal_SimpleGrid(Volume &vol, const SimpleGrid &grid,
   const ImageOver &footprint, const ImageOver &footprint2,
   Projection &proj, Projection &norm_proj,
   const matrix1D<double> &shift,
   const matrix1D<double> &aint, const matrix1D<double> &bint,
   const matrix2D<double> &D,  const matrix2D<double> &Dinv,
   const matrix2D<int> &mask, int FORW, int eq_mode) {
   matrix1D<double> prjX(3);                // Coordinate: Projection of the
   matrix1D<double> prjY(3);                // 3 grid vectors
   matrix1D<double> prjZ(3);
   matrix1D<double> prjOrigin(3);           // Coordinate: Where in the 2D
                                            // projection plane the origin of
                                            // the grid projects
   matrix1D<double> prjaint(3);             // Coordinate: Projection of the
   matrix1D<double> prjbint(3);             // 2 deformed lattice vectors
   matrix1D<double> prjDir;                 // Direction of projection
      	             	      	            // in the deformed space

   matrix1D<double> actprj(3);              // Coord: Actual position inside
                                            // the projection plane
   matrix1D<double> defactprj(3);           // Coord: Actual position inside
                                            // the deformed projection plane
   matrix1D<double> beginZ(3);              // Coord: Plane coordinates of the
                                            // projection of the 3D point
                                            // (z0,YY(lowest),XX(lowest))
   matrix1D<double> beginY(3);              // Coord: Plane coordinates of the
                                            // projection of the 3D point
                                            // (z0,y0,XX(lowest))
   matrix1D<double> footprint_size(2);      // The footprint is supposed
      	             	                    // to be defined between
                                            // (-vmax,+vmax) in the Y axis,
                                            // and (-umax,+umax) in the X axis
                                            // This footprint size is the
                                            // R2 vector (umax,vmax).
   matrix1D<double> deffootprint_size(2);   // The footprint size
      	             	      	            // in the deformed space
   int XX_corner2, XX_corner1;              // Coord: Corners of the
   int YY_corner2, YY_corner1;              // footprint when it is projected
                                            // onto the projection plane
   matrix1D<double> rc(2), r(2);            // Position vector which will
      	             	      	            // move from corner1 to corner2.
					    // In rc the wrapping is not
					    // considered, while it is in r
   int           foot_V, foot_U;            // Img Coord: coordinate
                                            // corresponding to the blobprint
                                            // point which matches with this
                                            // pixel position
   double        vol_corr;                  // Correction for a volum element
   int           N_eq;                      // Number of equations in which
                                            // a blob is involved
   int           i,j,k;                     // volume element indexes
   SPEED_UP_temps;   	      	            // Auxiliar variables for
      	             	      	            // fast multiplications

   // Compute the deformed direction of projection .........................
   matrix2D<double> Eulerg;
   Eulerg=proj.euler*D;
   
   // Compute deformation in the projection space ..........................
   // The following two vectors are defined in the deformed volume space
   VECTOR_R3(actprj,XX(aint),YY(aint),0);
      grid.Gdir_project_to_plane(actprj, Eulerg, prjaint);
   VECTOR_R3(actprj,XX(bint),YY(bint),0);
      grid.Gdir_project_to_plane(actprj, Eulerg, prjbint);
   #ifdef DEBUG_LITTLE
      double rot, tilt, psi;
      Euler_matrix2angles(proj.euler,rot,tilt,psi);
      cout << "rot= "  << rot << " tilt= " << tilt
           << " psi= " << psi << endl;
      cout << "D\n" << D << endl;
      cout << "Eulerf\n" << proj.euler << endl;
      cout << "Eulerg\n" << Eulerg << endl;
      cout << "aint    " << aint.transpose() << endl;
      cout << "bint    " << bint.transpose() << endl;
      cout << "prjaint " << prjaint.transpose() << endl;
      cout << "prjbint " << prjbint.transpose() << endl;
      cout.flush();
   #endif
   // Project grid axis ....................................................
   // These vectors ((1,0,0),(0,1,0),...) are referred to the grid
   // coord. system and the result is a 2D vector in the image plane
   // The unit size within the image plane is 1, ie, the same as for
   // the universal grid.
   // Be careful that these grid vectors are defined in the deformed
   // volume space, and the projection are defined in the deformed
   // projections, 
   VECTOR_R3(actprj,1,0,0);
      grid.Gdir_project_to_plane(actprj, Eulerg, prjX);
   VECTOR_R3(actprj,0,1,0);
      grid.Gdir_project_to_plane(actprj, Eulerg, prjY);
   VECTOR_R3(actprj,0,0,1);
      grid.Gdir_project_to_plane(actprj, Eulerg, prjZ);

   // This time the origin of the grid is in the universal coord system
   // but in the deformed space
   Uproject_to_plane(grid.origin,Eulerg,prjOrigin);

   // This is a correction used by the crystallographic symmetries
   prjOrigin += XX(shift)*prjaint + YY(shift)*prjbint;

   #ifdef DEBUG_LITTLE
      cout << "prjX      " << prjX.transpose() << endl;
      cout << "prjY      " << prjY.transpose() << endl;
      cout << "prjZ      " << prjZ.transpose() << endl;
      cout << "prjOrigin " << prjOrigin.transpose() << endl;
      cout.flush();
   #endif

   // Now I will impose that prja becomes (Xdim,0) and prjb, (0,Ydim)
   // A is a matrix such that
   // A*prjaint=(Xdim,0)'
   // A*prjbint=(0,Ydim)'
   matrix2D<double> A(2,2), Ainv(2,2);
   A(0,0)= YY(prjbint)*xDim;   A(0,1)=-XX(prjbint)*xDim;
   A(1,0)=-YY(prjaint)*yDim;   A(1,1)= XX(prjaint)*yDim;
   double nor=1/(XX(prjaint)*YY(prjbint)-XX(prjbint)*YY(prjaint));
   M2x2_BY_CT(A,A,nor);
   M2x2_INV(Ainv,A);

   #ifdef DEBUG_LITTLE
      cout << "A\n" << A << endl;
      cout << "Ainv\n" << Ainv << endl;
      cout << "Check that A*prjaint=(Xdim,0)    "
           << (A*vector_R2(XX(prjaint),YY(prjaint))).transpose() << endl;
      cout << "Check that A*prjbint=(0,Ydim)    "
           << (A*vector_R2(XX(prjbint),YY(prjbint))).transpose() << endl;
      cout << "Check that Ainv*(Xdim,0)=prjaint "
           << (Ainv*vector_R2(xDim,0)).transpose() << endl;
      cout << "Check that Ainv*(0,Ydim)=prjbint "
           << (Ainv*vector_R2(0,yDim)).transpose() << endl;
      cout.flush();
   #endif

   // Footprint size .......................................................
   // The following vectors are integer valued vectors, but they are
   // stored as real ones to make easier operations with other vectors.
   // Look out that this footprint size is in the deformed projection,
   // that is why it is not a square footprint but a parallelogram
   // and we have to look for the smaller rectangle which surrounds it
   XX(footprint_size) = footprint.umax();
   YY(footprint_size) = footprint.vmax();
   N_eq=(2*footprint.umax()+1)*(2*footprint.vmax()+1);
   matrix1D<double> c1(3), c2(3);
   XX(c1)    	      = XX(footprint_size);
   YY(c1)    	      = YY(footprint_size);
   XX(c2)    	      =-XX(c1);
   YY(c2)    	      = YY(c1);
   M2x2_BY_V2x1(c1,A,c1);
   M2x2_BY_V2x1(c2,A,c2);
   XX(deffootprint_size)=MAX(ABS(XX(c1)),ABS(XX(c2)));
   YY(deffootprint_size)=MAX(ABS(YY(c1)),ABS(YY(c2)));

   #ifdef DEBUG_LITTLE
      cout << "Footprint_size " << footprint_size.transpose() << endl;
      cout << "c1: " << c1.transpose() << endl;
      cout << "c2: " << c2.transpose() << endl;
      cout << "Deformed Footprint_size " << deffootprint_size.transpose()
           << endl;
      cout.flush();
   #endif

   // This type conversion gives more speed
   int ZZ_lowest =(int) ZZ(grid.lowest);
   int YY_lowest =MAX((int) YY(grid.lowest),STARTINGY(mask));
   int XX_lowest =MAX((int) XX(grid.lowest),STARTINGX(mask));
   int ZZ_highest=(int) ZZ(grid.highest);
   int YY_highest=MIN((int) YY(grid.highest),FINISHINGY(mask));
   int XX_highest=MIN((int) XX(grid.highest),FINISHINGX(mask));

   // Project the whole grid ...............................................
   // Corner of the plane defined by Z. These coordinates try to be within
   // the valid indexes of the projection (defined between STARTING and
   // FINISHING values, but it could be that they may lie outside.
   // These coordinates need not to be integer, in general, they are
   // real vectors.
   // The vectors returned by the projection routines are R3 but we
   // are only interested in their first 2 components, ie, in the
   // in-plane components
   beginZ=XX_lowest*prjX + YY_lowest*prjY + ZZ_lowest*prjZ + prjOrigin;

   #ifdef DEBUG_LITTLE
      cout << "BeginZ     " << beginZ.transpose()             << endl;
      cout << "Mask       "; mask.print_shape();         cout << endl;
      cout << "Vol        "; vol().print_shape();        cout << endl;
      cout << "Proj       "; proj().print_shape();       cout << endl;
      cout << "Norm Proj  "; norm_proj().print_shape();  cout << endl;
      cout << "Footprint  "; footprint().print_shape();  cout << endl;
      cout << "Footprint2 "; footprint2().print_shape(); cout << endl;
   #endif
   matrix1D<double> grid_index(3);
   for (k=ZZ_lowest; k<=ZZ_highest; k++) {
      // Corner of the row defined by Y
      beginY=beginZ;
      for (i=YY_lowest; i<=YY_highest; i++) {
         // First point in the row
	 actprj=beginY;
	 for (j=XX_lowest; j<=XX_highest; j++) {
	    VECTOR_R3(grid_index,j,i,k);
            #ifdef DEBUG
               cout << "Visiting " << i << " " << j << endl;
            #endif

            // Be careful that you cannot skip any blob, although its
            // value be 0, because it is useful for norm_proj
	    // unless it doesn't belong to the reconstruction mask
	    if (MAT_ELEM(mask,i,j) && grid.is_interesting(grid_index)) {
      	       // Look for the position in the deformed projection
	       M2x2_BY_V2x1(defactprj,A,actprj);

               // Search for integer corners for this blob
               XX_corner1=CEIL (XX(defactprj)-XX(deffootprint_size));
               YY_corner1=CEIL (YY(defactprj)-YY(deffootprint_size));
               XX_corner2=FLOOR(XX(defactprj)+XX(deffootprint_size));
               YY_corner2=FLOOR(YY(defactprj)+YY(deffootprint_size));
	       
               #ifdef DEBUG
                  cout << "  k= " << k << " i= " << i << " j= " << j << endl;
                  cout << "  Actual position: " << actprj.transpose() << endl;
                  cout << "  Deformed position: " << defactprj.transpose() << endl;
                  cout << "  Blob corners: (" << XX_corner1 << "," << YY_corner1
                       << ") (" << XX_corner2 << "," << YY_corner2 << ")\n";
                  cout.flush();
               #endif

      	       // Now there is no way that both corners collapse into a line,
	       // since the projection points are wrapped

               if (!FORW) vol_corr=0;

	       // Effectively project this blob
	       // (xc,yc) is the position of the considered pixel
	       // in the crystal undeformed projection
	       // When wrapping and deformation are considered then x and
	       // y are used
      	       #define xc XX(rc)
	       #define yc YY(rc)
      	       #define x  XX(r)
	       #define y  YY(r)
	       for (y=YY_corner1; y<=YY_corner2; y++)
	          for (x=XX_corner1; x<=XX_corner2; x++) {
      	             // Compute position in undeformed space and
		     // corresponding sample in the blob footprint
		     M2x2_BY_V2x1(rc,Ainv,r);
		     OVER2IMG(footprint,yc-YY(actprj), xc-XX(actprj),
		        foot_V, foot_U);
                     
                     #ifdef DEBUG
                        cout << "    Studying: " << r.transpose()
                             << "--> " << rc.transpose()
                             << "dist (" << xc-XX(actprj)
                             << "," << yc-YY(actprj) << ") --> "
                             << foot_U << " " << foot_V << endl;
                        cout.flush();
                     #endif
                     if (IMGMATRIX(footprint).outside(foot_V,foot_U)) continue;

      	             // Wrap positions
      	             int yw, xw;
		     wrap_as_Crystal(x,y,xw,yw);
                     #ifdef DEBUG
                        cout << "      After wrapping " << xw << " " << yw << endl;
                        cout << "      Value added = " << VOLVOXEL(vol,k,i,j) *
                           IMGPIXEL(footprint,foot_V,foot_U) << " Blob value = "
                           << IMGPIXEL(footprint,foot_V,foot_U) 
                           << " Blob^2 " << IMGPIXEL(footprint2,foot_V,foot_U)
                           << endl;
                        cout.flush();
                     #endif
	             if (FORW) {
                        IMGPIXEL(proj,yw,xw) += VOLVOXEL(vol,k,i,j) *
                           IMGPIXEL(footprint,foot_V,foot_U);
                        switch(eq_mode) {
                           case CAVARTK:
                           case ARTK:
                              IMGPIXEL(norm_proj,yw,xw) +=
                                 IMGPIXEL(footprint2,foot_V,foot_U);
                              break;
                        }
	             } else {
	                vol_corr += IMGPIXEL(norm_proj,yw,xw) *
                           IMGPIXEL(footprint,foot_V,foot_U);
	             }

	          }

               #ifdef DEBUG_INTERMIDIATE
                  proj.write("inter.xmp");
                  norm_proj.write("inter_norm.xmp");
                  cout << "Press any key\n"; char c; cin >> c;
               #endif

	       if (!FORW)
                  switch (eq_mode) {
                     case ARTK:
                        VOLVOXEL(vol,k,i,j) += vol_corr;
                        break;
                     case CAVARTK:
                        VOLVOXEL(vol,k,i,j) += vol_corr/N_eq;
                        break;
                  }
	    }

            // Prepare for next iteration
            V2_PLUS_V2(actprj, actprj, prjX);
      	 }
         V2_PLUS_V2(beginY, beginY, prjY);
      }
      V2_PLUS_V2(beginZ, beginZ, prjZ);
   }
   //#define DEBUG_AT_THE_END
   #ifdef DEBUG_AT_THE_END
      proj.write("inter.xmp");
      norm_proj.write("inter_norm.xmp");
      cout << "Press any key\n"; char c; cin >> c;
   #endif
}
#undef DEBUG
#undef DEBUG_LITTLE
#undef wrap_as_Crystal

/* Project a Grid Volume --------------------------------------------------- */
//#define DEBUG
void project_Crystal_Volume(
   GridVolume &vol,                      // Volume
   const ImageOver  &footprint,          // Footprint of a blob         
   const ImageOver  &footprint2,         // Used for error computation
   Projection       &proj,               // Projection
   Projection       &norm_proj,          // Projection of a unitary volume
   int              Ydim,                // Dimensions of the projection
   int              Xdim,
   double rot, double tilt, double psi,  // Euler angles
   const matrix1D<double> &shift,        // Shift to apply to projections
   const matrix1D<double> &aint,         // First lattice vector (2x1) in voxels
   const matrix1D<double> &bint,         // Second lattice vector (2x1) in voxels
   const matrix2D<double> &D,            // volume deformation matrix
   const matrix2D<double> &Dinv,         // volume deformation matrix
   const matrix2D<int>   &mask,          // volume mask
   int              FORW,                // 1 if we are projecting a volume
                                         //   norm_proj is calculated
                                         // 0 if we are backprojecting
                                         //   norm_proj must be valid
   int eq_mode)                          // ARTK, CAVARTK, CAVK or CAV
{
   // If projecting forward initialise projections
   if (FORW) {
      proj.reset(Ydim,Xdim);
      proj.set_angles(rot,tilt,psi);
      norm_proj().resize(proj());
   }
   
   // Project each subvolume
   for (int i=0; i<vol.VolumesNo(); i++) {
      project_Crystal_SimpleGrid(vol(i),vol.grid(i),footprint,footprint2,
         proj, norm_proj, shift, aint, bint, D, Dinv, mask, FORW, eq_mode);

   #ifdef DEBUG
      ImageXmipp save; save=norm_proj;
      if (FORW) save.write((string)"PPPnorm_FORW"+(char)(48+i));
      else      save.write((string)"PPPnorm_BACK"+(char)(48+i));
   #endif
   }
}
#undef DEBUG

/* Count equations in a grid volume --------------------------------------- */
void count_eqs_in_projection(GridVolumeT<int> &GVNeq,
   const struct blobtype &blob,
   const ImageOver &footprint, const ImageOver &footprint2,
   Projection &read_proj) {
   for (int i=0; i<GVNeq.VolumesNo(); i++)
      project_SimpleGrid(GVNeq(i),GVNeq.grid(i),blob,footprint,footprint2,
         read_proj, read_proj, FORWARD, COUNT_EQ, NULL, NULL);
}

/* Instantiation ----------------------------------------------------------- */
void instatiate_projection() {
   GridVolumeT<int>    GVi;
   GridVolumeT<double> GVd;
   ImageOver           footprint;
   Projection          proj;
   struct blobtype     blob;
   project_Volume(GVi,blob,footprint,footprint,proj,proj,30,30,0,0,0,FORWARD,ARTK);
   project_Volume(GVd,blob,footprint,footprint,proj,proj,30,30,0,0,0,FORWARD,ARTK);
}
