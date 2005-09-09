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
//Fri Nov 19 13:22:15 EST 1999 subsampling strategy modified (R.Marabini)
/* ------------------------------------------------------------------------- */
/* PHANTOMS                                                                  */
/* ------------------------------------------------------------------------- */

#include <stdio.h>
#include "../phantom.hh"
#include <XmippData/xmippGeometry.hh>

/* ######################################################################### */
/* Features                                                                  */
/* ######################################################################### */
/* ------------------------------------------------------------------------- */
/* Prepare                                                                   */
/* ------------------------------------------------------------------------- */
void Sphere::prepare() {
   max_distance=radius;
}

void Blob::prepare() {
   max_distance=radius;
}

void Cylinder::prepare() {
   prepare_Euler();
   max_distance=sqrt(height*height/4+MAX(xradius*xradius,yradius*yradius));
}

void DCylinder::prepare() {
   prepare_Euler();
   max_distance=sqrt((height+separation)*(height+separation)/4+radius*radius);
}

void Cube::prepare() {
   prepare_Euler();
   max_distance=sqrt(xdim*xdim+ydim*ydim+zdim*zdim);
}

void Ellipsoid::prepare() {
   prepare_Euler();
   max_distance=MAX(MAX(xradius,yradius),zradius);
}

void Cone::prepare() {
   prepare_Euler();
   max_distance=sqrt(height*height/4+radius*radius);
}

/* ------------------------------------------------------------------------- */
/* Assignment                                                                */
/* ------------------------------------------------------------------------- */
Feature & Feature::operator = (const Feature &F) {
   if (this==&F) return *this;
   Type         = F.Type;
   Add_Assign   = F.Add_Assign;
   Density      = F.Density;
   Center       = F.Center;
   max_distance = F.max_distance;
   return *this;
}

Oriented_Feature & Oriented_Feature::operator = (const Oriented_Feature &OF) {
   if (this==&OF) return *this;
   rot    = OF.rot;
   tilt   = OF.tilt;
   psi    = OF.psi;
   euler  = OF.euler;
   eulert = OF.eulert;
   Feature::operator = (OF);
   return *this;
}

Sphere & Sphere::operator = (const Sphere &F) {
   if (this==&F) return *this;
   Feature::operator = (F);
   radius = F.radius;
   return *this;
}

Blob & Blob::operator = (const Blob &F) {
   if (this==&F) return *this;
   Feature::operator = (F);
   radius = F.radius;
   alpha  = F.alpha;
   m      = F.m;
   return *this;
}

Cylinder & Cylinder::operator = (const Cylinder &F) {
   if (this==&F) return *this;
   Oriented_Feature::operator = (F);
   xradius = F.xradius;
   yradius = F.yradius;
   height  = F.height;
   return *this;
}

DCylinder & DCylinder::operator = (const DCylinder &F) {
   if (this==&F) return *this;
   Oriented_Feature::operator = (F);
   radius     = F.radius;
   height     = F.height;
   separation = F.separation;
   return *this;
}

Cube & Cube::operator = (const Cube &F) {
   if (this==&F) return *this;
   Oriented_Feature::operator = (F);
   xdim = F.xdim;
   ydim = F.ydim;
   zdim = F.zdim;
   return *this;
}

Ellipsoid & Ellipsoid::operator = (const Ellipsoid &F) {
   if (this==&F) return *this;
   Oriented_Feature::operator = (F);
   xradius = F.xradius;
   yradius = F.yradius;
   zradius = F.zradius;
   return *this;
}

Cone & Cone::operator = (const Cone &F) {
   if (this==&F) return *this;
   Oriented_Feature::operator = (F);
   radius     = F.radius;
   height     = F.height;
   return *this;
}

/* ------------------------------------------------------------------------- */
/* Rotation                                                                  */
/* ------------------------------------------------------------------------- */
void Feature::rotate_center(const matrix2D<double> &E) {
   Center=E*Center;
}

void Feature::rotate(const matrix2D<double> &E) {
   rotate_center(E);
}

void Oriented_Feature::rotate(const matrix2D<double> &E) {
   rotate_center(E);
   prepare();
   euler=euler*E;
   eulert=E.transpose()*eulert;
   Euler_matrix2angles(euler,rot,tilt,psi);
}

/* ------------------------------------------------------------------------- */
/* I/O functions                                                             */
/* ------------------------------------------------------------------------- */
/* Read common part of features -------------------------------------------- */
void Feature::read_common(char *line) _THROW {
   int        stat;
   char       straux[6];
   Center.resize(3);
   stat=sscanf(line,"%s %c %lf %lf %lf %lf",
      straux,
      &(Add_Assign),
      &(Density),
      &(XX(Center)),
      &(YY(Center)),
      &(ZZ(Center)));
   if (stat!=6)
      REPORT_ERROR(3003,
         (string)"Error when reading common part of feature: " + line);
   Type=straux;
}

/* Read a sphere ----------------------------------------------------------- */
void Sphere::read_specific(char *line) _THROW {
   int stat;
   stat=sscanf(line,"%*s %*c %*f %*f %*f %*f %lf",&radius);
   if (stat!=1)
      REPORT_ERROR(3003,(string)"Error when reading a sphere:"+line);
   prepare();
}

/* Read a blob ----------------------------------------------------------- */
void Blob::read_specific(char *line) _THROW {
   int stat;
   stat=sscanf(line,"%*s %*c %*f %*f %*f %*f %lf %lf %d",&radius,&alpha,&m);
   if (stat!=3)
      REPORT_ERROR(3003,(string)"Error when reading a blob:"+line);
   prepare();
}

/* Read a Cylinder --------------------------------------------------------- */
void Cylinder::read_specific(char *line) _THROW {
   int stat;
   stat=sscanf(line,"%*s %*c %*f %*f %*f %*f %lf %lf %lf %lf %lf %lf",&xradius,
      &yradius,&height,&rot, &tilt, &psi);
   if (stat!=6)
      REPORT_ERROR(3003,(string)"Error when reading a cylinder:"+line);
   prepare();
}
 
/* Read a Double Cylinder -------------------------------------------------- */
void DCylinder::read_specific(char *line) _THROW {
   int stat;
   stat=sscanf(line,"%*s %*c %*f %*f %*f %*f %lf %lf %lf %lf %lf %lf",
      &radius,&height, &separation, &rot, &tilt, &psi);
   if (stat!=6)
      REPORT_ERROR(3003,(string)"Error when reading a double cylinder:"+line);
   prepare();
}
 
/* Read a Cube ------------------------------------------------------------- */
void Cube::read_specific(char *line) _THROW {
   int stat;
   stat=sscanf(line,"%*s %*c %*f %*f %*f %*f %lf %lf %lf %lf %lf %lf",
      &xdim,&ydim,&zdim,&rot, &tilt, &psi);
   if (stat!=6)
      REPORT_ERROR(3003,(string)"Error when reading a cube"+line);
   prepare();
}
 
/* Read an Ellipsoid ------------------------------------------------------- */
void Ellipsoid::read_specific(char *line) _THROW {
   int stat;
   stat=sscanf(line,"%*s %*c %*f %*f %*f %*f %lf %lf %lf %lf %lf %lf",
      &xradius,&yradius,&zradius,&rot, &tilt, &psi);
   if (stat!=6)
      REPORT_ERROR(3003,(string)"Error when reading an ellipsoid"+line);
   prepare();
}

/* Read a Cone ------------------------------------------------------------- */
void Cone::read_specific(char *line) _THROW {
   int stat;
   stat=sscanf(line,"%*s %*c %*f %*f %*f %*f %lf %lf %lf %lf %lf",&radius,&height,
      &rot, &tilt, &psi);
   if (stat!=5)
      REPORT_ERROR(3003,(string)"Error when reading a cone"+line);
   prepare();
}
 
/* Show an sphere ---------------------------------------------------------- */
void  Sphere::feat_printf(FILE *fh) const {
   fprintf(fh,"sph    %c     %1.4f    % 7.2f   % 7.2f    % 7.2f    % 7.2f\n",
       Add_Assign, Density, XX(Center), YY(Center), ZZ(Center),
       radius);
}

/* Show a Blob    ---------------------------------------------------------- */
void  Blob::feat_printf(FILE *fh) const {
   fprintf(fh,"blo    %c     %1.4f    % 7.2f   % 7.2f    % 7.2f    % 7.2f"
              "    % 7.2f    %1d\n",
       Add_Assign, Density, XX(Center), YY(Center), ZZ(Center),
       radius, alpha,m);
}

/* Show a cylinder --------------------------------------------------------- */
void  Cylinder::feat_printf(FILE *fh) const {
   fprintf(fh,"cyl    %c     %1.4f    % 7.2f   % 7.2f    % 7.2f    % 7.2f    "
       "% 7.2f    % 7.2f    % 7.2f    % 7.2f    % 7.2f\n",
       Add_Assign, Density, XX(Center), YY(Center), ZZ(Center),
       xradius,yradius,height,
       rot,tilt,psi);
}

/* Show a double cylinder -------------------------------------------------- */
void  DCylinder::feat_printf(FILE *fh) const {
   fprintf(fh,"dcy    %c     %1.4f    % 7.2f   % 7.2f    % 7.2f    % 7.2f    "
       "% 7.2f    % 7.2f    % 7.2f    % 7.2f    % 7.2f\n",
       Add_Assign, Density, XX(Center), YY(Center), ZZ(Center),
       radius,height,separation,
       rot,tilt,psi);
}

/* Show a cube ------------------------------------------------------------- */
void  Cube::feat_printf(FILE *fh) const {
   fprintf(fh,"cub    %c     %1.4f    % 7.2f   % 7.2f    % 7.2f    % 7.2f    "
       "% 7.2f    % 7.2f    % 7.2f    % 7.2f   % 7.2f\n",
       Add_Assign, Density, XX(Center), YY(Center), ZZ(Center),
       xdim,ydim,zdim,
       rot,tilt,psi);
}

/* Show an ellipsoid ------------------------------------------------------- */
void  Ellipsoid::feat_printf(FILE *fh) const {
   fprintf(fh,"ell    %c     %1.4f    % 7.2f   % 7.2f    % 7.2f    % 7.2f    "
       "% 7.2f    % 7.2f    % 7.2f    % 7.2f   % 7.2f\n",
       Add_Assign, Density, XX(Center), YY(Center), ZZ(Center),
       xradius,yradius,zradius,
       rot,tilt,psi);
}

/* Show a cone ------------------------------------------------------------- */
void  Cone::feat_printf(FILE *fh) const {
   fprintf(fh,"con    %c     %1.4f    % 7.2f   % 7.2f    % 7.2f    % 7.2f    "
       "% 7.2f    % 7.2f    % 7.2f    % 7.2f\n",
       Add_Assign, Density, XX(Center), YY(Center), ZZ(Center),
       radius,height,
       rot,tilt,psi);
}

/* Show feat --------------------------------------------------------------- */
ostream& operator << (ostream &o, const Feature *F){
   if (F!=NULL) {
      o << "Feature --------" << endl;
      o << "   Type:        " << F->Type << endl;
      o << "   Add_Assign:  " << F->Add_Assign << endl;
      o << "   Density:     " << F->Density << endl;
      o << "   Center:      " << F->Center.transpose() << endl;
      if      (F->Type=="sph") o << *((Sphere *) F);
      else if (F->Type=="blo") o << *((Blob *) F);
      else if (F->Type=="cyl") o << *((Cylinder *) F);
      else if (F->Type=="dcy") o << *((DCylinder *) F);
      else if (F->Type=="cub") o << *((Cube *) F);
      else if (F->Type=="ell") o << *((Ellipsoid *) F);
      else if (F->Type=="con") o << *((Cone *) F);
   }
   return o;
}

/* Show sphere ------------------------------------------------------------- */
ostream& operator << (ostream &o, const Sphere &f) {
   o << "   Radius: " << f.radius << endl;
   return o;
}

/* Show Blob   ------------------------------------------------------------- */
ostream& operator << (ostream &o, const Blob &f) {
   o << "   Radius: "  << f.radius << endl;
   o << "   Alpha:  "  << f.alpha << endl;
   o << "   m:      "  << f.m << endl;
   return o;
}

/* Show cylinder ----------------------------------------------------------- */
ostream& operator << (ostream &o, const Cylinder &f) {
   o << "   XRadius: " << f.xradius << endl;
   o << "   YRadius: " << f.yradius << endl;
   o << "   Height:  " << f.height << endl;
   o << "   Rot:     " << f.rot << endl;
   o << "   Tilt:    " << f.tilt << endl;
   o << "   Psi:     " << f.psi << endl;
   return o;
}

/* Show double cylinder ---------------------------------------------------- */
ostream& operator << (ostream &o, const DCylinder &f) {
   o << "   Radius: " << f.radius << endl;
   o << "   Height: " << f.height << endl;
   o << "   Separ.: " << f.separation << endl;
   o << "   Rot:    " << f.rot << endl;
   o << "   Tilt:   " << f.tilt << endl;
   o << "   Psi:    " << f.psi << endl;
   return o;
}

/* Show cube --------------------------------------------------------------- */
ostream& operator << (ostream &o, const Cube &f) {
   o << "   Xdim:   " << f.xdim << endl;
   o << "   Ydim:   " << f.ydim << endl;
   o << "   Zdim:   " << f.zdim << endl;
   o << "   Rot:    " << f.rot << endl;
   o << "   Tilt:   " << f.tilt << endl;
   o << "   Psi:    " << f.psi << endl;
   return o;
}

/* Show ellipsoid ---------------------------------------------------------- */
ostream& operator << (ostream &o, const Ellipsoid &f) {
   o << "   Xradius: " << f.xradius << endl;
   o << "   Yradius: " << f.yradius << endl;
   o << "   Zradius: " << f.zradius << endl;
   o << "   Rot:     " << f.rot << endl;
   o << "   Tilt:    " << f.tilt << endl;
   o << "   Psi:     " << f.psi << endl;
   return o;
}

/* Show cone --------------------------------------------------------------- */
ostream& operator << (ostream &o, const Cone &f) {
   o << "   Radius: " << f.radius << endl;
   o << "   Height: " << f.height << endl;
   o << "   Rot:    " << f.rot << endl;
   o << "   Tilt:   " << f.tilt << endl;
   o << "   Psi:    " << f.psi << endl;
   return o;
}

/* ------------------------------------------------------------------------- */
/* Point Inside                                                              */
/* ------------------------------------------------------------------------- */
// For speed reasons an auxiliar vector of length 3 must be supplied to each
// function

#define DEF_Sph_Blob_point_inside {\
  /*Express r in the feature coord. system*/\
  V3_MINUS_V3(aux,r,Center);\
  /*Check if it is inside*/\
  if (XX(aux)*XX(aux) + YY(aux)*YY(aux) +ZZ(aux)*ZZ(aux) <= radius*radius)\
      return 1;\
  return 0;}
  
/* Point inside a sphere --------------------------------------------------- */
int Sphere::point_inside(const matrix1D<double> &r, matrix1D<double> &aux) const {
DEF_Sph_Blob_point_inside
}

/* Point inside a Blob --------------------------------------------------- */
int Blob::point_inside(const matrix1D<double> &r, matrix1D<double> &aux) const {
DEF_Sph_Blob_point_inside
}
#undef DEF_Sph_Blob_point_inside

/* Density inside a Blob --------------------------------------------------- */
double Blob::density_inside(const matrix1D<double> &r, matrix1D<double> &aux) const
{
/*Express r in the feature coord. system*/
V3_MINUS_V3(aux,r,Center);
/*Calculate density*/
return (kaiser_value( aux.module(), radius,  alpha,  m));
}

/* Point inside a cylinder ------------------------------------------------- */
int Cylinder::point_inside(const matrix1D<double> &r, matrix1D<double> &aux) const {
   SPEED_UP_temps;
   double tx,ty;

   // Express r in the feature coord. system
   V3_MINUS_V3(aux,r,Center);
   M3x3_BY_V3x1(aux,euler,aux);

   // Check if it is inside
   tx=XX(aux)/xradius; ty=YY(aux)/yradius;
   if (tx*tx + ty*ty<=1.0
      && ABS(ZZ(aux))<=height/2) return 1;
   return 0;
}

/* Point inside a Double cylinder ------------------------------------------ */
int DCylinder::point_inside(const matrix1D<double> &r, matrix1D<double> &aux) const {
   SPEED_UP_temps;

   // Express r in the feature coord. system
   V3_MINUS_V3(aux,r,Center);
   M3x3_BY_V3x1(aux,euler,aux);

   // Check if inside
   if (XX(aux)*XX(aux)+ YY(aux)*YY(aux)<=radius*radius) {
      double cyl_center=separation/2+height/2;
      if      (ABS(ZZ(aux)-cyl_center)<=height/2) return 1;
      else if (ABS(ZZ(aux)+cyl_center)<=height/2) return 1;
   }
   return 0;
}

/* Point inside a cube ----------------------------------------------------- */
int Cube::point_inside(const matrix1D<double> &r, matrix1D<double> &aux) const {
   SPEED_UP_temps;

   // Express r in the feature coord. system
   V3_MINUS_V3(aux,r,Center);
   M3x3_BY_V3x1(aux,euler,aux);

   // Check if inside
   if (ABS(XX(aux))<=xdim/2 && ABS(YY(aux))<=ydim/2 &&
       ABS(ZZ(aux))<=zdim/2 ) return 1;
   return 0;
}

/* Point inside an ellipsoid ----------------------------------------------- */
int Ellipsoid::point_inside(const matrix1D<double> &r, matrix1D<double> &aux) const {
   SPEED_UP_temps;
   double tx,ty,tz;

   // Express r in the feature coord. system
   V3_MINUS_V3(aux,r,Center);
   M3x3_BY_V3x1(aux,euler,aux);

   // Check if inside
   tx=XX(aux)/xradius; ty=YY(aux)/yradius; tz=ZZ(aux)/zradius;
   if (tx*tx + ty*ty + tz*tz<=1.0) return 1;
   return 0;
}

/* Point inside a cone ----------------------------------------------------- */
int Cone::point_inside(const matrix1D<double> &r, matrix1D<double> &aux) const {
   SPEED_UP_temps;
   double Zradius;

   // Express r in the feature coord. system
   V3_MINUS_V3(aux,r,Center);
   M3x3_BY_V3x1(aux,euler,aux);

   // Check if inside
   if (ABS(ZZ(aux))<=height/2) {
      Zradius=radius*(1-(ZZ(aux)+height/2)/height);
      if (XX(aux)*XX(aux)+ YY(aux)*YY(aux)<=Zradius*Zradius) return 1;
   }
   return 0;
}

/* Voxel inside ------------------------------------------------------------ */
// In all functions the voxelside is supposed to be 1
//#define DEBUG
#ifdef DEBUG
   #define DEBUG_SHOW \
      if (ZZ(r)==0 && YY(r)==0) \
         cout << "Point (z=" << ZZ(aux1) << ",y=" << YY(aux1) << ",x=" \
              << XX(aux1) << ") inside=" << inside << endl;
#else
   #define DEBUG_SHOW
#endif
int Feature::voxel_inside(const matrix1D<double> &r, matrix1D<double> &aux1,
   matrix1D<double> &aux2) const {

   // The subvoxels are visited following a Gray code, so the number
   // of operations is minimized
   XX(aux1)=XX(r)+0.25; YY(aux1)=YY(r)+0.25; ZZ(aux1)=ZZ(r)+0.25; // 000
      int inside = point_inside(aux1,aux2);
   DEBUG_SHOW;
   ZZ(aux1)-=0.5; inside += point_inside(aux1,aux2);              // 001
   DEBUG_SHOW;
   YY(aux1)-=0.5; inside += point_inside(aux1,aux2);              // 011
   DEBUG_SHOW;
   ZZ(aux1)+=0.5; inside += point_inside(aux1,aux2);              // 010
   DEBUG_SHOW;
   XX(aux1)-=0.5; inside += point_inside(aux1,aux2);              // 110
   DEBUG_SHOW;
   ZZ(aux1)-=0.5; inside += point_inside(aux1,aux2);              // 111
   DEBUG_SHOW;
   YY(aux1)+=0.5; inside += point_inside(aux1,aux2);              // 101
   DEBUG_SHOW;
   ZZ(aux1)+=0.5; inside += point_inside(aux1,aux2);              // 100
   DEBUG_SHOW;
   return inside;
}
/* voxel_inside_by_normalized_density ------------------------------------*/
double Feature::voxel_inside_by_normalized_density(
                             const matrix1D<double> &r, 
			     matrix1D<double> &aux1,
   matrix1D<double> &aux2) const {
#ifdef NEVER
   if(Type=="blo")
     {
cout << "den=" <<   density_inside(r, aux2) << endl;  
     return(density_inside(r, aux2));
     }
   else
     return((double)voxel_inside(r,aux1, aux2));
#endif
   // The subvoxels are visited following a Gray code, so the number
   // of operations is minimized
   XX(aux1)=XX(r)+0.25; YY(aux1)=YY(r)+0.25; ZZ(aux1)=ZZ(r)+0.25; // 000
      double inside = (double)point_inside(aux1,aux2)
                             *density_inside(r, aux2);
   DEBUG_SHOW;
   ZZ(aux1)-=0.5; inside += (double)point_inside(aux1,aux2)
                                    *density_inside(r, aux2); // 001
   DEBUG_SHOW;
   YY(aux1)-=0.5; inside += (double)point_inside(aux1,aux2)
                                  *density_inside(r, aux2);   // 011
   DEBUG_SHOW;
   ZZ(aux1)+=0.5; inside += (double)point_inside(aux1,aux2)
                                   *density_inside(r, aux2); // 010
   DEBUG_SHOW;
   XX(aux1)-=0.5; inside += (double)point_inside(aux1,aux2)
                                   *density_inside(r, aux2); // 110
   DEBUG_SHOW;
   ZZ(aux1)-=0.5; inside += (double)point_inside(aux1,aux2)
                                   *density_inside(r, aux2); // 111
   DEBUG_SHOW;
   YY(aux1)+=0.5; inside += (double)point_inside(aux1,aux2)
                                   *density_inside(r, aux2); // 101
   DEBUG_SHOW;
   ZZ(aux1)+=0.5; inside += (double)point_inside(aux1,aux2)
                                   *density_inside(r, aux2); // 100
   DEBUG_SHOW;
   return inside;
     
}
#undef DEBUG

/* Intersects sphere ------------------------------------------------------- */
int Feature::intersects_sphere(const matrix1D<double> &r, double radius,
   matrix1D<double> &aux1, matrix1D<double> &aux2, matrix1D<double> &aux3)
   const {
   double radius2=radius*radius;
   bool intersects=FALSE;
   for (double k=FLOOR(ZZ(r)-radius); k<=CEIL(ZZ(r)+radius) && !intersects; k++)
      for (double i=FLOOR(YY(r)-radius); i<=CEIL(YY(r)+radius) && !intersects; i++)
	 for (double j=FLOOR(XX(r)-radius); j<=CEIL(XX(r)+radius) && !intersects; j++) {
	    if ((k-ZZ(r))*(k-ZZ(r))+(i-YY(r))*(i-YY(r))+(j-XX(r))*(j-XX(r))>
	        radius2) continue;
	    VECTOR_R3(aux3,j,i,k);
	    intersects=voxel_inside(aux3,aux1,aux2);
	 }
   return intersects;
}

/* ------------------------------------------------------------------------- */
/* Draw in                                                                   */
/* ------------------------------------------------------------------------- */
/* Corners ----------------------------------------------------------------- */
void Feature::corners(const Volume *V, matrix1D<double> &corner1,
   matrix1D<double> &corner2) {
   corner1.resize(3);
   corner2.resize(3);
   XX(corner1)=MAX(FLOOR(XX(Center)-max_distance),STARTINGX(VOLMATRIX(*V)));
   YY(corner1)=MAX(FLOOR(YY(Center)-max_distance),STARTINGY(VOLMATRIX(*V)));
   ZZ(corner1)=MAX(FLOOR(ZZ(Center)-max_distance),STARTINGZ(VOLMATRIX(*V)));
   XX(corner2)=MIN(CEIL (XX(Center)+max_distance),FINISHINGX(VOLMATRIX(*V)));
   YY(corner2)=MIN(CEIL (YY(Center)+max_distance),FINISHINGY(VOLMATRIX(*V)));
   ZZ(corner2)=MIN(CEIL (ZZ(Center)+max_distance),FINISHINGZ(VOLMATRIX(*V)));

#ifdef PORSI
   array_by_scalar(Center,max_distance,corner1,'-');
   array_by_scalar(Center,max_distance,corner2,'+');
   corner1=FLOORnD(corner1);
   corner2=CEILnD(corner2);
   XX(corner1)=MAX(XX(corner1),STARTINGX(VOLMATRIX(*V)));
   YY(corner1)=MAX(YY(corner1),STARTINGY(VOLMATRIX(*V)));
   ZZ(corner1)=MAX(ZZ(corner1),STARTINGZ(VOLMATRIX(*V)));
   XX(corner2)=MIN(XX(corner2),FINISHINGX(VOLMATRIX(*V)));
   YY(corner2)=MIN(YY(corner2),FINISHINGY(VOLMATRIX(*V)));
   ZZ(corner2)=MIN(ZZ(corner2),FINISHINGZ(VOLMATRIX(*V)));
#endif
}

/* Draw a feature ---------------------------------------------------------- */
//#define DEBUG
#define Vr VOLVOXEL((*V),(int)ZZ(r),(int)YY(r),(int)XX(r))
void Feature::draw_in(Volume *V, int colour_mode, double colour) {
   matrix1D<double>   aux1(3), aux2(3), corner1(3), corner2(3), r(3);
   int               add;
   double             inside;
   double             final_colour;
   
   if (colour_mode==INTERNAL) {final_colour=Density; add=Add_Assign=='+';}
   else                       {final_colour=colour;  add=0;}

   corners(V,corner1,corner2);
   #ifdef DEBUG
      cout << "Drawing \n";
      cout << this;
      cout << "colour_mode=" << colour_mode << endl;
      cout << "Add_Assign= " << Add_Assign  << endl;
      cout << "add=        " << add         << endl;
      cout << "   Corner 1" << corner1.transpose() << endl;
      cout << "   Corner 2" << corner2.transpose() << endl;
   #endif
   FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1,corner2) {
      inside=voxel_inside_by_normalized_density(r,aux1, aux2);
      #ifdef DEBUG
         //int condition=(ZZ(r)==-12) && (YY(r)==1);
         int condition=1;
         if (condition)
            cout << "   r=" << r.transpose() << " inside= " << inside;
      #endif
      if (inside!=0) {
         double drawing_colour=final_colour*inside/8;
         if (add) Vr +=drawing_colour;
         else     Vr  =MAX(drawing_colour,Vr);
         #ifdef DEBUG
            if (condition)
               cout << "   V(r)=" << VOLVOXEL((*V),(int)ZZ(r),(int)YY(r),(int)XX(r));
         #endif
      }
      #ifdef DEBUG
         if (condition)
            cout << endl;
      #endif
   }
}
#undef DEBUG

/* Sketch a feature -------------------------------------------------------- */
void Feature::sketch_in(Volume *V, double colour) {
   matrix1D<double>   aux1(3),aux2(3), corner1(3), corner2(3), r(3);
   int               inside;
      
   corners(V,corner1,corner2);
   FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1,corner2) {
      inside=voxel_inside(r,aux1,aux2);
      if (inside!=0 && inside!=8)
         VOLVOXEL((*V),(int)ZZ(r),(int)YY(r),(int)XX(r))=colour;
   }
}

/* Shift a feature --------------------------------------------------------- */
void Feature::shift(double shiftX, double shiftY, double shiftZ) {
   XX(Center) += shiftX;
   YY(Center) += shiftY;
   ZZ(Center) += shiftZ;
}

/* Apply a general transformation to a feature ------------------------------ */
void Feature::self_apply_geom(const matrix2D<double> &A) {
   matrix1D<double> r(4);
   XX(r)=XX(Center);
   YY(r)=YY(Center);
   ZZ(r)=ZZ(Center);
   r(3)=1;
   r=A*r;
   XX(Center)=XX(r);
   YY(Center)=YY(r);
   ZZ(Center)=ZZ(r);
}

/* ------------------------------------------------------------------------- */
/* Intersection                                                              */
/* ------------------------------------------------------------------------- */
// A line is supposed to be defined as a direction vector and a passing point
// this way the parametric equation of the line is
// (x,y,z)=(x1,y1,z1)+t(dx,dy,dz)
// where (x,y,z)    is the generic point belonging to this line
//       (x1,y1,z1) is the passing point
//       (dx,dy,dz) is the direction vector
//       t          is a free parameter

double Sphere::intersection(
   const matrix1D<double> &direction,
   const matrix1D<double> &passing_point,
   matrix1D<double> &r,
   matrix1D<double> &u) const {
   // This is done in order to correct the different lengths seen by
   // rays with different "speed". It is related to the jacobian of
   // the transformation from a non-unit direction to a unit one.
   double norm=direction.module();

   // Set the passing point in the ellipsoid coordinate system
   // and normalise to a unit sphere
   V3_MINUS_V3(r,passing_point,Center);
   V3_BY_CT(r, r, 1/radius);
   V3_BY_CT(u, direction, 1/radius);
      return intersection_unit_sphere(u,r)/norm;
}
double Blob::intersection(
   const matrix1D<double> &direction,
   const matrix1D<double> &passing_point,
   matrix1D<double> &r,
   matrix1D<double> &u) const {
#ifdef NEVERDEFINED
   // This is done in order to correct the different lengths seen by
   // rays with different "speed". It is related to the jacobian of
   // the transformation from a non-unit direction to a unit one.
   double norm=direction.module();

   // Set the passing point in the ellipsoid coordinate system
   // and normalise to a unit sphere
   V3_MINUS_V3(r,passing_point,Center);
   V3_BY_CT(r, r, 1/radius);
   V3_BY_CT(u, direction, 1/radius);
   return intersection_unit_sphere(u,r)/norm;
#endif
   return(kaiser_proj(point_line_distance_3D(Center,passing_point,direction), 
                      radius,alpha, m));

}

//#define DEBUG
double Cylinder::intersection(
   const matrix1D<double> &direction,
   const matrix1D<double> &passing_point,
   matrix1D<double> &r,
   matrix1D<double> &u) const {
   double norm=direction.module();
   SPEED_UP_temps;

   // Set the passing point in the cylinder coordinate system
   // and normalise to a unit cylinder
   V3_MINUS_V3(r,passing_point,Center);
   M3x3_BY_V3x1(r,euler,r);
   XX(r) /= xradius;
   YY(r) /= yradius;
   ZZ(r) /= height;

   // Express also the direction in the cyilinder coordinate system
   // and normalise to a unit cylinder
   M3x3_BY_V3x1(u,euler,direction);
   XX(u) /= xradius;
   YY(u) /= yradius;
   ZZ(u) /= height;

   #ifdef DEBUG
   cout << "Intersecting .-.-.-.-.-.-.\n";
   cout << *this;
   cout << "   direction(Univ) = " << direction.transpose() << endl;
   cout << "   passing  (Univ) = " << passing_point.transpose() << endl;
   cout << "   direction(Obj.) = " << u.transpose() << endl;
   cout << "   passing  (Obj.) = " << r.transpose() << endl;
   cout << "   intersection    = " << intersection_unit_cylinder(u,r) << endl;
   #endif

   return intersection_unit_cylinder(u,r)/norm;
}
#undef DEBUG

double DCylinder::intersection(
   const matrix1D<double> &direction,
   const matrix1D<double> &passing_point,
   matrix1D<double> &r,
   matrix1D<double> &u) const {
   double norm=direction.module();
   SPEED_UP_temps;

   // Express also the direction in the cylinder coordinate system
   // and normalise to a unit cylinder
   M3x3_BY_V3x1(u,euler,direction);
   XX(u) /= radius;
   YY(u) /= radius;
   ZZ(u) /= height;

   // Top cylinder
   // Set the passing point in the cylinder coordinate system
   // and normalise to a unit cylinder
   V3_MINUS_V3(r,passing_point,Center);
   M3x3_BY_V3x1(r,euler,r);
   ZZ(r) -= (separation/2+height/2);
   XX(r) /= radius;
   YY(r) /= radius;
   ZZ(r) /= height;
   double i1=intersection_unit_cylinder(u,r);

   // Bottom cylinder
   // Set the passing point in the cylinder coordinate system
   // and normalise to a unit cylinder
   V3_MINUS_V3(r,passing_point,Center);
   M3x3_BY_V3x1(r,euler,r);
   ZZ(r) += (separation/2+height/2);
   XX(r) /= radius;
   YY(r) /= radius;
   ZZ(r) /= height;
   double i2=intersection_unit_cylinder(u,r);
   
   return (i1+i2)/norm;
}

double Cube::intersection(
   const matrix1D<double> &direction,
   const matrix1D<double> &passing_point,
   matrix1D<double> &r,
   matrix1D<double> &u) const {
   double norm=direction.module();
   SPEED_UP_temps;

   // Set the passing point in the cube coordinate system
   // and normalise to a unit cube
   V3_MINUS_V3(r,passing_point,Center);
   M3x3_BY_V3x1(r,euler,r);
   XX(r) /= xdim;
   YY(r) /= ydim;
   ZZ(r) /= zdim;
   
   // Express also the direction in the cube coordinate system
   // and normalise to a unit cube
   M3x3_BY_V3x1(u,euler,direction);
   XX(u) /= xdim;
   YY(u) /= ydim;
   ZZ(u) /= zdim;
   
   return intersection_unit_cube(u,r)/norm;
}

double Ellipsoid::intersection(
   const matrix1D<double> &direction,
   const matrix1D<double> &passing_point,
   matrix1D<double> &r,
   matrix1D<double> &u) const {
   double norm=direction.module();
   SPEED_UP_temps;

   // Set the passing point in the ellipsoid coordinate system
   // and normalise to a unit sphere
   V3_MINUS_V3(r,passing_point,Center);
   M3x3_BY_V3x1(r,euler,r);
   XX(r) /= xradius;
   YY(r) /= yradius;
   ZZ(r) /= zradius;
   
   // Express also the direction in the ellipsoid coordinate system
   // and normalise to a unit sphere
   M3x3_BY_V3x1(u,euler,direction);
   XX(u) /= xradius;
   YY(u) /= yradius;
   ZZ(u) /= zradius;
   
   return intersection_unit_sphere(u,r)/norm;
}

double Cone::intersection(
   const matrix1D<double> &direction,
   const matrix1D<double> &passing_point,
   matrix1D<double> &r,
   matrix1D<double> &u) const {
   return 0;
}

/* ------------------------------------------------------------------------- */
/* Projecting                                                                */
/* ------------------------------------------------------------------------- */
/* Project a feature to a plane -------------------------------------------- */
//#define DEBUG_LITTLE
//#define DEBUG
//#define DEBUG_EVEN_MORE
void Feature::project_to(Projection &P, const matrix2D<double> &VP,
   const matrix2D<double> &PV) const {
   #define SUBSAMPLING 2                  // for every measure 2x2 line
                                          // integrals will be taken to
                                          // avoid numerical errors
   #define SUBSTEP 1/(SUBSAMPLING*2.0)

   matrix1D<double> origin(3);
   matrix1D<double> direction;
      VP.getRow(2,direction);
      direction.self_transpose();
   matrix1D<double> corner1(3), corner2(3);
   matrix1D<double> act(3);
   SPEED_UP_temps;

   // Find center of the feature in the projection plane ...................
   // Step 1). Project the center to the plane, the result is in the
   //          universal coord system
   M3x3_BY_V3x1(origin,VP,Center);

//   matrix1D<double> origin_debug(3);
//   Uproject_to_plane(Center,P.direction,0,origin_debug);
   
//#define DEBUG_LITTLE
   #ifdef DEBUG_LITTLE
      cout << "Actual feature\n"     << this << endl;
      cout << "Center              " << Center.transpose() << endl;
      cout << "VP matrix\n"          << VP << endl;
      cout << "P.direction         " << P.direction.transpose() << endl;
      cout << "direction           " << direction.transpose() << endl;
      cout << "P.euler matrix      " << P.euler << endl;
      cout << "max_distance        " << max_distance << endl;
      cout << "origin              " << origin.transpose() << endl;
//      cout << "origin_debug (Univ.coord) " << origin_debug.transpose() << endl;
   #endif
/*   
   // Step 2). Express this projected center in the projection coord system
   M3x3_BY_V3x1(origin_debug,P.euler,origin_debug);
//   if (A!=NULL) M2x2_BY_V2x1(origin,*A,origin_);
   #ifdef DEBUG_LITTLE
      cout << "origin (Proj.coord) " << origin_debug.transpose() << endl;
   #endif
*/
   
   // Find limits for projection ...........................................
   // Choose corners for the projection of this feature. It is supposed
   // to have at the worst case a projection of size max_distance
   VECTOR_R3(corner1, max_distance, max_distance, max_distance);
   VECTOR_R3(corner2,-max_distance,-max_distance,-max_distance);
   #ifdef DEBUG_LITTLE
      cout << "Corner1 : " << corner1.transpose() << endl
           << "Corner2 : " << corner2.transpose() << endl;
   #endif

   box_enclosing(corner1,corner2,VP,corner1,corner2);
//   if (A!=NULL) {
//      rectangle_enclosing(corner1,corner2,*A,corner1,corner2);
      #ifdef DEBUG_LITTLE
         cout << "Corner1 moves to : " << corner1.transpose() << endl
              << "Corner2 moves to : " << corner2.transpose() << endl;
      #endif
//   }

   V3_PLUS_V3(corner1,origin,corner1);
   V3_PLUS_V3(corner2,origin,corner2);
   #ifdef DEBUG_LITTLE
      cout << "Corner1 finally is : " << corner1.transpose() << endl
           << "Corner2 finally is : " << corner2.transpose() << endl;
   #endif
/*
   matrix1D<double> corner1_debug(2),corner2_debug(2);
   VECTOR_R2(corner1_debug, max_distance, max_distance);
   VECTOR_R2(corner2_debug,-max_distance,-max_distance);
   #ifdef DEBUG_LITTLE
      cout << "Corner1_debug : " << corner1_debug.transpose() << endl
           << "Corner2_debug : " << corner2_debug.transpose() << endl;
   #endif
   V2_PLUS_V2(corner1_debug,origin_debug,corner1_debug);
   V2_PLUS_V2(corner2_debug,origin_debug,corner2_debug);
   #ifdef DEBUG_LITTLE
      cout << "Corner1_debug finally is : " << corner1_debug.transpose() << endl
           << "Corner2_debug finally is : " << corner2_debug.transpose() << endl;
   #endif
*/
   // Discard not necessary components
   corner1.resize(2);
   corner2.resize(2);

   // Clip to image size
   sort_two_vectors(corner1,corner2);
   XX(corner1)=CLIP(ROUND(XX(corner1)),STARTINGX(P()),FINISHINGX(P()));
   YY(corner1)=CLIP(ROUND(YY(corner1)),STARTINGY(P()),FINISHINGY(P()));
   XX(corner2)=CLIP(ROUND(XX(corner2)),STARTINGX(P()),FINISHINGX(P()));
   YY(corner2)=CLIP(ROUND(YY(corner2)),STARTINGY(P()),FINISHINGY(P()));

   #ifdef DEBUG_LITTLE
      cout << "corner1      " << corner1.transpose() << endl;
      cout << "corner2      " << corner2.transpose() << endl;
      cout.flush();
   #endif

   // Check if there is something to project
   if (XX(corner1)==XX(corner2)) return;
   if (YY(corner1)==YY(corner2)) return;
   
   // Study the projection for each point in the projection plane ..........
   // (u,v) are in the deformed projection plane (if any deformation)
   for (int v=(int)YY(corner1); v<=(int)YY(corner2); v++)
      for (int u=(int)XX(corner1); u<=(int)XX(corner2); u++) {
         double length=0;
         #ifdef DEBUG_EVEN_MORE
            cout << "Studying point (" << u << "," << v << ")\n";
            cout.flush();
         #endif

         // Perform subsampling ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
         double u0=u-(int)(SUBSAMPLING/2.0)*SUBSTEP;
         double v0=v-(int)(SUBSAMPLING/2.0)*SUBSTEP;
         double actv=v0;
         for (int subv=0; subv<SUBSAMPLING; subv++) {
            double actu=u0;
            for (int subu=0; subu<SUBSAMPLING; subu++) {
               // Compute the coordinates of point (subu,subv) which is
               // within the plane in the universal coordinate system
               XX(act)=actu;
               YY(act)=actv;
               ZZ(act)=0;
//               if (Ainv!=NULL) M2x2_BY_V2x1(act,*Ainv,act);
//               M3x3_BY_V3x1(act,P.eulert,act);
               M3x3_BY_V3x1(act,PV,act);
               
               // Compute the intersection of a ray which passes through
               // this point and its direction is perpendicular to the
               // projection plane
               double possible_length=intersection(direction,act);
               if (possible_length>0) length += possible_length;
               
               #ifdef DEBUG_EVEN_MORE
               cout << "Averaging at (" << actu << "," << actv << ")\n";
               cout << "   which in univ. coords is " << act.transpose() << endl;
               cout << "   intersection there " << possible_length << endl;
               #endif
               // Prepare for next iteration
               actu +=SUBSTEP*2.0;
               }
            actv +=SUBSTEP*2.0;
         }
         length /=(SUBSAMPLING*SUBSAMPLING);
         #ifdef DEBUG
         cout << "Final value added at position (" << u << "," << v << ")="
              << length << endl;
         #endif
         
         // Add at the correspondant pixel the found intersection ,,,,,,,,,,
         IMGPIXEL(P,v,u) += length*Density;
       }
}
#undef DEBUG_LITTLE
#undef DEBUG
#undef DEBUG_EVEN_MORE

/* ------------------------------------------------------------------------- */
/* Scaling by a factor                                                       */
/* ------------------------------------------------------------------------- */
#define COPY_COMMON_PART \
   f->Type         = Type; \
   f->Add_Assign   = Add_Assign; \
   f->Density      = Density; \
   f->Center       = Center;

#define COPY_ANGLES \
   f->rot          = rot; \
   f->tilt         = tilt; \
   f->psi          = psi;

/* Scale a sphere ---------------------------------------------------------- */
Feature * Sphere::scale(double factor) const {
   Sphere *f;
   f=new Sphere;
   COPY_COMMON_PART;

   f->radius       = factor * radius;
   f->prepare();
   return (Feature *)f;
}

/* Scale a blob ---------------------------------------------------------- */
Feature * Blob::scale(double factor) const {
   Blob *f;
   f=new Blob;
   COPY_COMMON_PART;

   f->radius       = factor * radius;
   f->alpha        = alpha;
   f->m		   = m;
   f->prepare();
   return (Feature *)f;
}
   
/* Scale a cylinder -------------------------------------------------------- */
Feature * Cylinder::scale(double factor) const {
   Cylinder *f;
   f=new Cylinder;
   COPY_COMMON_PART;
   COPY_ANGLES;
   
   f->xradius      = factor * xradius;
   f->yradius      = factor * yradius;
   f->height       = factor * height;
   f->prepare();
   return (Feature *)f;
}
   
/* Scale a double cylinder ------------------------------------------------- */
Feature * DCylinder::scale(double factor) const {
   DCylinder *f;
   f=new DCylinder;
   COPY_COMMON_PART;
   COPY_ANGLES;

   f->radius       = factor * radius;
   f->height       = factor * height;
   f->separation   = separation - 2*(factor-1)*height;
   f->prepare();
   
   return (Feature *)f;
}
   
/* Scale a cube ------------------------------------------------------------ */
Feature * Cube::scale(double factor) const {
   Cube *f;
   f=new Cube;
   COPY_COMMON_PART;
   COPY_ANGLES;

   f->xdim         = factor * xdim;
   f->ydim         = factor * ydim;
   f->zdim         = factor * zdim;
   f->prepare();
   return (Feature *)f;
}
   
/* Scale an ellipsoid ------------------------------------------------------ */
Feature * Ellipsoid::scale(double factor) const {
   Ellipsoid *f;
   f=new Ellipsoid;
   COPY_COMMON_PART;
   COPY_ANGLES;

   f->xradius      = factor * xradius;
   f->yradius      = factor * yradius;
   f->zradius      = factor * zradius;
   f->prepare();
   return (Feature *)f;
}
   
/* Scale a cone ------------------------------------------------------------ */
Feature * Cone::scale(double factor) const {
   Cone *f;
   f=new Cone;
   COPY_COMMON_PART;
   COPY_ANGLES;

   f->radius       = factor * radius;
   f->height       = factor * height;
   f->prepare();
   return (Feature *)f;
}

#undef COPY_COMMON_PART
#undef COPY_ANGLES

/* ------------------------------------------------------------------------- */
/* Backgrounds                                                               */
/* ------------------------------------------------------------------------- */
/* Encircle any feature ---------------------------------------------------- */
Feature *Feature::encircle(double radius) const {
   Sphere *f;
   f=new Sphere;
   
   if (radius==0) radius=1.5*max_distance;
   
   f->Type         = "sph";
   f->Add_Assign   = Add_Assign;
   f->Density      = Density;
   f->Center       = Center;
   f->max_distance = radius;
   f->radius       = radius;
   
   return (Feature *)f;
}

Feature *Feature::background(int back_mode, double back_param) const _THROW {
   switch(back_mode) {
      case ENLARGE_MODE: return scale(back_param); break;
      case SPHERE_MODE:  return encircle(back_param); break;
      default: REPORT_ERROR(3006,"Feature::background: mode not supported");
         break;
   }
}

/* ------------------------------------------------------------------------- */
/* Init random                                                               */
/* ------------------------------------------------------------------------- */
void Sphere::init_rnd(
   double minradius, double maxradius,
   double minden,    double maxden,
   double minx0,     double maxx0,
   double miny0,     double maxy0,
   double minz0,     double maxz0) {

   randomize_random_generator();
   Center.resize(3);
   Type           = "sph";
   Add_Assign     = '+';
   Density        = rnd_unif(minden,maxden);
   XX(Center)     = rnd_unif(minx0,maxx0);
   YY(Center)     = rnd_unif(miny0,maxy0);
   ZZ(Center)     = rnd_unif(minz0,maxz0);

   radius         = rnd_unif(minradius,maxradius);

   max_distance   = radius;
}
void Blob::init_rnd(
   double minradius, double maxradius,
   double minalpha,  double maxalpha,
   double minorder,  double maxorder,
   double minden,    double maxden,
   double minx0,     double maxx0,
   double miny0,     double maxy0,
   double minz0,     double maxz0) {

   randomize_random_generator();
   Center.resize(3);
   Type           = "blo";
   Add_Assign     = '+';
   Density        = rnd_unif(minden,maxden);
   XX(Center)     = rnd_unif(minx0,maxx0);
   YY(Center)     = rnd_unif(miny0,maxy0);
   ZZ(Center)     = rnd_unif(minz0,maxz0);

   radius         = rnd_unif(minradius,maxradius);
   alpha	  = rnd_unif(minalpha,maxalpha);
   m		  = (int)(rnd_unif(minorder,maxorder)+0.5);
   max_distance   = radius;
}

void Cylinder::init_rnd(
   double minxradius,  double maxxradius,
   double minyradius,  double maxyradius,
   double minheight,   double maxheight,
   double minden,      double maxden,
   double minx0,       double maxx0,
   double miny0,       double maxy0,
   double minz0,       double maxz0,
   double minrot,      double maxrot,
   double mintilt,     double maxtilt,
   double minpsi,      double maxpsi) {

   randomize_random_generator();
   Center.resize(3);
   Type           = "cyl";
   Add_Assign     = '+';
   Density        = rnd_unif(minden,maxden);
   XX(Center)     = rnd_unif(minx0,maxx0);
   YY(Center)     = rnd_unif(miny0,maxy0);
   ZZ(Center)     = rnd_unif(minz0,maxz0);

   xradius        = rnd_unif(minxradius,maxxradius);
   yradius        = rnd_unif(minyradius,maxyradius);
   height         = rnd_unif(minheight,maxheight);
   rot            = rnd_unif(minrot,maxrot);
   tilt           = rnd_unif(mintilt,maxtilt);
   psi            = rnd_unif(minpsi,maxpsi);
   Euler_angles2matrix(rot,tilt,psi,euler);
   eulert=euler.transpose();

   max_distance   = sqrt(height*height+MAX(xradius*xradius,yradius*yradius));
}

void DCylinder::init_rnd(
   double minradius,   double maxradius,
   double minheight,   double maxheight,
   double minsep,      double maxsep,
   double minden,      double maxden,
   double minx0,       double maxx0,
   double miny0,       double maxy0,
   double minz0,       double maxz0,
   double minrot,      double maxrot,
   double mintilt,     double maxtilt,
   double minpsi,      double maxpsi) {

   randomize_random_generator();
   Center.resize(3);
   Type           = "dcy";
   Add_Assign     = '+';
   Density        = rnd_unif(minden,maxden);
   XX(Center)     = rnd_unif(minx0,maxx0);
   YY(Center)     = rnd_unif(miny0,maxy0);
   ZZ(Center)     = rnd_unif(minz0,maxz0);

   radius         = rnd_unif(minradius,maxradius);
   height         = rnd_unif(minheight,maxheight);
   separation     = rnd_unif(minsep,maxsep);
   rot            = rnd_unif(minrot,maxrot);
   tilt           = rnd_unif(mintilt,maxtilt);
   psi            = rnd_unif(minpsi,maxpsi);
   Euler_angles2matrix(rot,tilt,psi,euler);
   eulert=euler.transpose();

   max_distance   = sqrt((height+separation)*(height+separation)/4
                         +radius*radius);
}

void Cube::init_rnd(
   double minXdim,     double maxXdim,
   double minYdim,     double maxYdim,
   double minZdim,     double maxZdim,
   double minden,      double maxden,
   double minx0,       double maxx0,
   double miny0,       double maxy0,
   double minz0,       double maxz0,
   double minrot,      double maxrot,
   double mintilt,     double maxtilt,
   double minpsi,      double maxpsi) {

   randomize_random_generator();
   Center.resize(3);
   Type           = "cub";
   Add_Assign     = '+';
   Density        = rnd_unif(minden,maxden);
   XX(Center)     = rnd_unif(minx0,maxx0);
   YY(Center)     = rnd_unif(miny0,maxy0);
   ZZ(Center)     = rnd_unif(minz0,maxz0);

   if (minYdim==0) minYdim=minXdim;
   if (minZdim==0) minZdim=minXdim;
   if (maxYdim==0) maxYdim=maxXdim;
   if (maxZdim==0) maxZdim=maxXdim;
   xdim           = rnd_unif(minXdim,maxXdim);
   ydim           = rnd_unif(minYdim,maxYdim);
   zdim           = rnd_unif(minZdim,maxZdim);
   rot            = rnd_unif(minrot,maxrot);
   tilt           = rnd_unif(mintilt,maxtilt);
   psi            = rnd_unif(minpsi,maxpsi);
   Euler_angles2matrix(rot,tilt,psi,euler);
   eulert=euler.transpose();

   max_distance   = sqrt(xdim*xdim+ydim*ydim+zdim*zdim);
}

void Ellipsoid::init_rnd(
   double minXradius,  double maxXradius,
   double minYradius,  double maxYradius,
   double minZradius,  double maxZradius,
   double minden,      double maxden,
   double minx0,       double maxx0,
   double miny0,       double maxy0,
   double minz0,       double maxz0,
   double minrot,      double maxrot,
   double mintilt,     double maxtilt,
   double minpsi,      double maxpsi) {

   randomize_random_generator();
   Center.resize(3);
   Type           = "ell";
   Add_Assign     = '+';
   Density        = rnd_unif(minden,maxden);
   XX(Center)     = rnd_unif(minx0,maxx0);
   YY(Center)     = rnd_unif(miny0,maxy0);
   ZZ(Center)     = rnd_unif(minz0,maxz0);

   if (minYradius==0) minYradius=minXradius;
   if (minZradius==0) minZradius=minXradius;
   if (maxYradius==0) maxYradius=maxXradius;
   if (maxZradius==0) maxZradius=maxXradius;
   xradius        = rnd_unif(minXradius,maxXradius);
   yradius        = rnd_unif(minYradius,maxYradius);
   zradius        = rnd_unif(minZradius,maxZradius);
   rot            = rnd_unif(minrot,maxrot);
   tilt           = rnd_unif(mintilt,maxtilt);
   psi            = rnd_unif(minpsi,maxpsi);
   Euler_angles2matrix(rot,tilt,psi,euler);
   eulert=euler.transpose();

   max_distance   = MAX(MAX(xradius,yradius),zradius);
}

void Cone::init_rnd(
   double minradius,   double maxradius,
   double minheight,   double maxheight,
   double minden,      double maxden,
   double minx0,       double maxx0,
   double miny0,       double maxy0,
   double minz0,       double maxz0,
   double minrot,      double maxrot,
   double mintilt,     double maxtilt,
   double minpsi,      double maxpsi) {

   randomize_random_generator();
   Center.resize(3);
   Type           = "con";
   Add_Assign     = '+';
   Density        = rnd_unif(minden,maxden);
   XX(Center)     = rnd_unif(minx0,maxx0);
   YY(Center)     = rnd_unif(miny0,maxy0);
   ZZ(Center)     = rnd_unif(minz0,maxz0);

   radius         = rnd_unif(minradius,maxradius);
   height         = rnd_unif(minheight,maxheight);
   rot            = rnd_unif(minrot,maxrot);
   tilt           = rnd_unif(mintilt,maxtilt);
   psi            = rnd_unif(minpsi,maxpsi);
   Euler_angles2matrix(rot,tilt,psi,euler);
   eulert=euler.transpose();

   max_distance   = MAX(radius,height);
}

/* ------------------------------------------------------------------------- */
/* Mean and Variance in a plane                                              */
/* ------------------------------------------------------------------------- */
void Feature::mean_variance_in_plane(Volume *V, double z, double &mean,
   double &var){
   double sum1=0;
   double sum2=0;
   double no_points=0;
   matrix1D<double> r(3), aux1(3), aux2(3);
   
   mean=0;
   var=0;
   if (z<STARTINGZ(VOLMATRIX(*V)) || z>FINISHINGZ(VOLMATRIX(*V))) return;
   
   ZZ(r)=z;
   for (YY(r)=STARTINGY(VOLMATRIX(*V)); YY(r)<=FINISHINGY(VOLMATRIX(*V)); YY(r)++)
      for (XX(r)=STARTINGX(VOLMATRIX(*V)); XX(r)<=FINISHINGX(VOLMATRIX(*V)); XX(r)++) {
         if (voxel_inside(r,aux1,aux2)==8) {
            double voxel=VOLVOXEL(*V,(int)ZZ(r),(int)YY(r),(int)XX(r));
            sum1 += voxel;
            sum2 += voxel*voxel;
            no_points++;
         }
      }
   if (no_points!=0) {
      mean = sum1/no_points;
      var  = sum2/no_points - mean*mean;
   }
}

/* ######################################################################### */
/* Phantoms                                                                  */
/* ######################################################################### */
/* Constructors ------------------------------------------------------------ */
Phantom::Phantom() {
   xdim = ydim = zdim = 0;
   Background_Density = 0;
   fn="";
   current_scale=1;
}

void Phantom::clear() {
   xdim = ydim = zdim = 0;
   Background_Density = 0;
   fn="";
   for (int i=0; i<VF.size(); i++) delete VF[i];
   VF.clear();
}

Phantom & Phantom::operator = (const Phantom &P) {
   if (&P==this) return *this;
   clear();
   fn=P.fn;
   xdim=P.xdim;
   ydim=P.ydim;
   zdim=P.zdim;
   Background_Density=P.Background_Density;
   Sphere     *sph;
   Blob       *blo;
   Cylinder   *cyl;
   DCylinder  *dcy;
   Cube       *cub;
   Ellipsoid  *ell;
   Cone       *con;
   for (int i=0; i<P.VF.size(); i++)
       if      (P.VF[i]->Type=="sph")
            {sph=new Sphere;    *sph=*((Sphere *)    P.VF[i]); add(sph);}
       else if (P.VF[i]->Type=="blo")
            {blo=new Blob;  *blo=*((Blob *)          P.VF[i]); add(blo);}
       else if (P.VF[i]->Type=="cyl")
            {cyl=new Cylinder;  *cyl=*((Cylinder *)  P.VF[i]); add(cyl);}
       else if (P.VF[i]->Type=="dcy")
            {dcy=new DCylinder; *dcy=*((DCylinder *) P.VF[i]); add(dcy);}
       else if (P.VF[i]->Type=="cub")
            {cub=new Cube;      *cub=*((Cube *)      P.VF[i]); add(cub);}
       else if (P.VF[i]->Type=="ell")
            {ell=new Ellipsoid; *ell=*((Ellipsoid *) P.VF[i]); add(ell);}
       else if (P.VF[i]->Type=="con")
            {con=new Cone;      *con=*((Cone      *) P.VF[i]); add(con);}
   return *this;
}

/* Prepare for work -------------------------------------------------------- */
void Phantom::prepare() {
   for (int i=0; i<VF.size(); i++) VF[i]->prepare();
}

/* Maximum distance -------------------------------------------------------- */
double Phantom::max_distance() const {
   double retval=0;
   for (int i=0; i<VF.size(); i++)
       retval=MAX(retval,VF[i]->max_distance+VF[i]->Center.module());
   return retval;
}

/* Volume ------------------------------------------------------------------ */
double Phantom::volume() const {
   double retval=0;
   for (int i=0; i<VF.size(); i++)
       retval+=VF[i]->volume();
   return retval;
}

/* Read Volume Description ------------------------------------------------- */
void Phantom::read(const FileName &fn_phantom, bool apply_scale) _THROW {

FILE *fh_phantom;
char line[201];
int Global_Feature_Read=0; // Indicates if the line with volume dimensions
                           // has been already read
int        stat;
Sphere     *sph;
Blob       *blo;
Cylinder   *cyl;
DCylinder  *dcy;
Cube       *cub;
Ellipsoid  *ell;
Cone       *con;
Feature    *feat, *scaled_feat;
string     feat_type;
double     scale;          // The scale factor is not stored
char       straux[6];

// Clear actual phantom
   clear();

// Open Volume Description File
   if ((fh_phantom = fopen(fn_phantom.c_str(), "r")) == NULL)
      REPORT_ERROR(3003,(string)"Phantom::read: Cannot open the phantom file: "
         +fn_phantom);
    fn = fn_phantom;
 
// Read the file
   while (fgets (line, 200,fh_phantom ) != NULL) {
      if (line[0]==0) continue;
      if (line[0]=='#') continue;
      if (line[0]=='\n') continue;

      // Read volume dimensions and global density .........................
      if (Global_Feature_Read==0) {
         Global_Feature_Read=1;
         stat=sscanf(line, "%d %d %d %lf %lf",&xdim,&ydim,&zdim,
            &Background_Density, &scale);
         if (stat <4)
            REPORT_ERROR(3003,"Phantom::read: God bless us, check the volume"
               " dimensions and global density in volume description file");
         if (stat==4) scale=1;
	 if (apply_scale) {
            xdim = (int) CEIL(scale*xdim);
            ydim = (int) CEIL(scale*ydim);
            zdim = (int) CEIL(scale*zdim);
	    current_scale=1;
	 } else current_scale=scale;
         continue;
      }

      // Read feature description ..........................................
      stat=sscanf(line,"%s",straux);
      feat_type=straux;
      if (stat!=1)
         REPORT_ERROR(3003,
            (string)"Phantom::read: Not correct feature type"+line);

      if        (feat_type=="sph") {
         sph = new Sphere; feat=sph; sph->read_common(line); sph->read_specific(line);
      } else if (feat_type=="blo") {
         blo = new Blob; feat=blo; blo->read_common(line); blo->read_specific(line);
      } else if (feat_type=="cyl") {
         cyl = new Cylinder; feat=cyl; cyl->read_common(line); cyl->read_specific(line);
      } else if (feat_type=="dcy") {
         dcy = new DCylinder; feat=dcy; dcy->read_common(line); dcy->read_specific(line);
      } else if (feat_type=="cub") {
         cub = new Cube; feat=cub; cub->read_common(line); cub->read_specific(line);
      } else if (feat_type=="ell") {
         ell = new Ellipsoid; feat=ell; ell->read_common(line); ell->read_specific(line);
      } else if (feat_type=="con") {
         con = new Cone; feat=con; con->read_common(line); con->read_specific(line);
      } else
         REPORT_ERROR(3003,(string)"Phantom::read: Unknown feature type: "+line);

      // Scale and Store feature
      if (apply_scale) {
	 scaled_feat=feat->scale(scale);
	 scaled_feat->Center = scaled_feat->Center * scale;
	 delete feat;

	 // Store feature
	 VF.push_back(scaled_feat);
      } else
         VF.push_back(feat);
   }
   fclose(fh_phantom);
}

/* Show whole phantom ------------------------------------------------------ */
ostream& operator << (ostream &o, const Phantom &P) {
   cout << "Phantom ---------------------------------------\n";
   cout << "Dimensions: " << P.xdim << " x " << P.ydim << " x " << P.zdim << endl;
   cout << "Background density: " << P.Background_Density << endl;
   for (int i=0; i<P.VF.size(); i++) o << P.VF[i];
   return o;
}

/* Write Volume Description ------------------------------------------------ */
void Phantom::write(const FileName &fn_phantom) _THROW {

FILE *fh_phantom;
char line[201];

// Open Volume Description File
   if ((fh_phantom = fopen(fn_phantom.c_str(), "w")) == NULL)
      REPORT_ERROR(3003,(string)"Phantom::write: Cannot open the phantom file "
         +fn_phantom +" for output");

// Write global comment and size
   fprintf(fh_phantom,"#Phantom Xdim Ydim Zdim Background density\n");
   fprintf(fh_phantom,"       %d    %d   %d   %f",xdim,ydim,zdim,Background_Density);
   if (current_scale!=1) fprintf(fh_phantom,"   %f",current_scale);
   fprintf(fh_phantom,"\n");

// Write description comment and features
   fprintf(fh_phantom,"#Type +/= Density X_Center Y_Center Z_Center\n");
   for (int i=0; i<VF.size(); i++) VF[i]->feat_printf(fh_phantom);
   fclose(fh_phantom);
}

/* Voxel Inside any feature ------------------------------------------------ */
int Phantom::voxel_inside_any_feat(const matrix1D<double> &r,
   matrix1D<double> &aux1, matrix1D<double> &aux2) const {
   int inside;
   for (int i=0; i<VF.size(); i++) {
      inside=VF[i]->voxel_inside(r, aux1, aux2);
      if (inside!=0) return i+1;
   }
   return 0;
}

/* Any feature intersects sphere ------------------------------------------- */
int Phantom::any_feature_intersects_sphere(const matrix1D<double> &r,
   double radius, matrix1D<double> &aux1, matrix1D<double> &aux2,
   matrix1D<double> &aux3) const {
   bool intersects;
   for (int i=0; i<VF.size(); i++) {
      intersects=VF[i]->intersects_sphere(r, radius, aux1, aux2, aux3);
      if (intersects) return i+1;
   }
   return 0;
}

/* Draw a Phantom ---------------------------------------------------------- */
// Always suppose CC grid
void Phantom::draw_in(Volume *V) {
   V->adapt_to_size(zdim,ydim,xdim);
   (*V)().init_constant(Background_Density);
   for (int i=0; i<VF.size(); i++) VF[i]->draw_in(V);
}

/* Label a Phantom --------------------------------------------------------- */
// Always suppose CC grid
void Phantom::label(Volume *V) {
   matrix1D<double> r(3), aux1(3), aux2(3);
   V->adapt_to_size(zdim,ydim,xdim);
   FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*V)) {
      ZZ(r)=k; YY(r)=i; XX(r)=j;
      int sel_feat=voxel_inside_any_feat(r,aux1,aux2);
      // If it is not in the background, check that it is completely
      // inside the feature, if not set it to border.
      if (sel_feat!=0)
         if (VF[sel_feat-1]->voxel_inside(r,aux1,aux2)!=8) sel_feat=-sel_feat;
      VOLVOXEL(*V,k,i,j)=sel_feat;
   }
}

/* Sketch a Phantom -------------------------------------------------------- */
// Always suppose CC grid
void Phantom::sketch_in(Volume *V, int clean, double colour) {
   if (clean) V->adapt_to_size(zdim,ydim,xdim);
   for (int i=0; i<VF.size(); i++) VF[i]->sketch_in(V, colour);
}

/* Shift a phantom --------------------------------------------------------- */
void Phantom::shift(double shiftX, double shiftY, double shiftZ) {
   for (int i=0; i<VF.size(); i++) VF[i]->shift(shiftX,shiftY,shiftZ);
}

/* Rotate a phantom -------------------------------------------------------- */
void Phantom::rotate(const matrix2D<double> &E) {
   for (int i=0; i<VF.size(); i++) VF[i]->rotate(E);
}

/* Apply geometrical transformatio to a phantom ---------------------------- */
void Phantom::self_apply_geom(const matrix2D<double> &A, int inv) _THROW {
   if ((XSIZE(A)!=4) || (YSIZE(A)!=4))
      REPORT_ERROR(1102,"Apply_geom3D: geometrical transformation is not 4x4");
   if (A.IsIdent()) return;
   matrix2D<double> T;
   if (inv==IS_INV) T=A.inv();
   else             T=A;

   for (int i=0; i<VF.size(); i++) VF[i]->self_apply_geom(T);
}

/* Projecting a phantom ---------------------------------------------------- */
//#define DEBUG
void Phantom::project_to(Projection &P, int Ydim, int Xdim,
   double rot, double tilt, double psi, const matrix2D<double> *A) const {
   #ifdef DEBUG
      cout << "Ydim=" << Ydim << " Xdim=" << Xdim << endl
           << "rot=" << rot << " tilt=" << tilt << " psi=" << psi << endl
           << "A\n" << A;
   #endif
   // Initialise projection
   P.adapt_to_size(Ydim, Xdim);
   P.set_angles(rot,tilt,psi);

   // Compute volume to Projection matrix
   matrix2D<double> VP=P.euler;
   if (A!=NULL) VP=(*A)*VP;
   matrix2D<double> PV=VP.inv();
   // Project all features
   for (int i=0; i<VF.size(); i++) VF[i]->project_to(P,VP,PV);
}
#undef DEBUG

void Phantom::project_to(Projection &P, 
   double rot, double tilt, double psi, const matrix2D<double> *A) const {
   P.set_angles(rot,tilt,psi);

   // Compute volume to Projection matrix
   matrix2D<double> VP=P.euler;
   if (A!=NULL) VP=(*A)*VP;
   matrix2D<double> PV=VP.inv();

   // Project all features
   for (int i=0; i<VF.size(); i++) VF[i]->project_to(P,VP,PV);
}

void Phantom::project_to(Projection &P, const matrix2D<double> &VP, double    disappearing_th) const {
   matrix2D<double> PV=VP.inv();

   // Project all features
   for (int i=0; i<VF.size(); i++) 
      {
      if(rnd_unif(0,1)<disappearing_th)
         VF[i]->project_to(P,VP,PV);
      }  
}

/* Surface ----------------------------------------------------------------- */
//#define DEBUG
void Phantom::surface(double z0, double radius, int direction, Image *P)
   const _THROW {
   if (z0!=zdim)
      if (z0<FIRST_XMIPP_INDEX(zdim) || z0>LAST_XMIPP_INDEX(zdim))
         REPORT_ERROR(1,"Phantom::surface: z0 outside phantom");
   #ifdef DEBUG
      cout << "Direction: " << direction << endl;
      cout << "z0:        " << z0        << endl;
      cout << "zdim:      " << zdim      << endl;
   #endif

   matrix1D<double> aux1(3), aux2(3), aux3(3), r(3);
   if (XSIZE((*P)())==0) P->adapt_to_size(ydim,xdim);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(*P)) {
      #ifdef DEBUG
         cout << "Processing (" << i << "," << j << ")" << endl;
      #endif
      // Init ray
      int k;
      if (direction==POS_NEG) k=LAST_XMIPP_INDEX(zdim)+1;
      else                    k=FIRST_XMIPP_INDEX(zdim)-1;
      bool finished;
      finished=FALSE;

      #ifdef DEBUG
         cout << "Initial k=" << k << endl;
      #endif
      // Check that it is not inside and move ray
      // at the end k takes the right value for the image
      while (!finished) {
         VECTOR_R3(r,j,i,k);
         #ifdef DEBUG
            cout << "Checking " << r.transpose() << endl;
         #endif
         // If it is inside move a step backward and finish
         if (any_feature_intersects_sphere(r,radius,aux1,aux2,aux3)) {
            finished=TRUE;
            if (direction==POS_NEG) k++; else k--;
         } else {
         // Else, move until you find z0
            if (z0!=zdim)
               if (direction==POS_NEG) {
                  k--;
                  if (k<z0) {finished=TRUE; k=CEIL(z0);}
               } else {
                  k++;
                  if (k>z0) {finished=TRUE; k=FLOOR(z0);}
               }
            // or you reach the end of the volume
            else
               if (direction==POS_NEG) {
                  k--;
                  if (k<FIRST_XMIPP_INDEX(zdim)) {finished=TRUE;}
               } else {
                  k++;
                  if (k>LAST_XMIPP_INDEX(zdim))  {finished=TRUE;}
               }
         }
      }
      
      IMGPIXEL(*P,i,j)=k;
   }
}
#undef DEBUG
