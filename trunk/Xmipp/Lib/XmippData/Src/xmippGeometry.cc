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

#include <iostream>
#include <math.h>

#include "../xmippGeometry.hh"
#include "../xmippFuncs.hh"

/* ######################################################################### */
/* Geometrical Operations                                                    */
/* ######################################################################### */

/* Project a point to a plane ---------------------------------------------- */
void Uproject_to_plane(const matrix1D<double> &point,
   const matrix1D<double> &direction, double distance,
   matrix1D<double> &result) {
   
   if (XSIZE(result)!=3) result.resize(3);
   double xx=distance-(XX(point)*XX(direction)+YY(point)*YY(direction)+
                      ZZ(point)*ZZ(direction));
   XX(result)=XX(point)+xx*XX(direction);
   YY(result)=YY(point)+xx*YY(direction);
   ZZ(result)=ZZ(point)+xx*ZZ(direction);
}

/* Project a point to a plane ---------------------------------------------- */
void Uproject_to_plane(const matrix1D<double> &r,
   double rot, double tilt, double psi, matrix1D<double> &result) {
   matrix2D<double> euler;
   Euler_angles2matrix(rot, tilt, psi, euler);
   Uproject_to_plane(r,euler,result);
}

/* Project a point to a plane ---------------------------------------------- */
void Uproject_to_plane(const matrix1D<double> &r,
   const matrix2D<double> &euler, matrix1D<double> &result) {
   SPEED_UP_temps;
   if (XSIZE(result)!=3) result.resize(3);
   M3x3_BY_V3x1(result,euler,r);
}

/* Spherical distance ------------------------------------------------------ */
double spherical_distance(const matrix1D<double> &r1, const matrix1D<double> &r2) {
   double r1r2=XX(r1)*XX(r2)+YY(r1)*YY(r2)+ZZ(r1)*ZZ(r2);
   double R1=sqrt(XX(r1)*XX(r1)+YY(r1)*YY(r1)+ZZ(r1)*ZZ(r1));
   double R2=sqrt(XX(r2)*XX(r2)+YY(r2)*YY(r2)+ZZ(r2)*ZZ(r2));
   double ang=acos(r1r2/(R1*R2));
   return ang*R1;
}
/* Point to line distance -------------------------------------------------- */

double point_line_distance_3D(const matrix1D<double> &p, 
                              const matrix1D<double> &a,
			      const matrix1D<double> &v)
			      
{
#ifdef NEVEREVER
double f,g,t;
double x,y,z;
double r;

f= (XX(x1)-XX(x0))*(XX(x2)-XX(x1))+
   (YY(x1)-YY(x0))*(YY(x2)-YY(x1))+
   (ZZ(x1)-ZZ(x0))*(ZZ(x2)-ZZ(x1));

x=XX(x2)-XX(x1);
y=YY(x2)-YY(x1);
z=ZZ(x2)-ZZ(x1);

g= x*x+y*y+z*z;

t=f/g;

x=(XX(x1)-XX(x0)+(XX(x2)-XX(x1))*t);
y=(YY(x1)-YY(x0)+(YY(x2)-YY(x1))*t);
z=(ZZ(x1)-ZZ(x0)+(ZZ(x2)-ZZ(x1))*t);

r=x*x+y*y+z*z; 
if(r<0)
  {
  cout << "Horror: The distance of a line to a point can not be negative"
       << "Congratulation you have found a bug in Xmipp."<< endl;
  exit(0);     
  }
#endif
matrix1D<double> p_a(3);

V3_MINUS_V3(p_a,p,a); 
return (vector_product(p_a,v).module()/v.module());
}			      

/* Point to plane distance ------------------------------------------------- */
double point_plane_distance_3D(const matrix1D<double> &p, 
                               const matrix1D<double> &a,
			       const matrix1D<double> &v)
			      
{
matrix1D<double> p_a(3);
V3_MINUS_V3(p_a,p,a); 
return (dot_product(p_a,v)/v.module());
}			      

/* Rectangle enclosing ----------------------------------------------------- */
void rectangle_enclosing(const matrix1D<double> &v0, const matrix1D<double> &vF,
    const matrix2D<double> &V, matrix1D<double> &corner1,
    matrix1D<double> &corner2) {
   SPEED_UP_temps;
   matrix1D<double> v(2);
   corner1.resize(2);
   corner2.resize(2);
   
   // Store values for reusing input as output vectors
   double XX_v0=XX(v0);
   double YY_v0=YY(v0);
   double XX_vF=XX(vF);
   double YY_vF=YY(vF);

   VECTOR_R2(v,XX_v0,YY_v0);
   M2x2_BY_V2x1(v,V,v);
   XX(corner1)=XX(v); XX(corner2)=XX(v);
   YY(corner1)=YY(v); YY(corner2)=YY(v);

   #define DEFORM_AND_CHOOSE_CORNERS2D \
      M2x2_BY_V2x1(v,V,v); \
      XX(corner1)=MIN(XX(corner1),XX(v)); XX(corner2)=MAX(XX(corner2),XX(v)); \
      YY(corner1)=MIN(YY(corner1),YY(v)); YY(corner2)=MAX(YY(corner2),YY(v));

   VECTOR_R2(v,XX_vF,YY_v0); DEFORM_AND_CHOOSE_CORNERS2D;
   VECTOR_R2(v,XX_v0,YY_vF); DEFORM_AND_CHOOSE_CORNERS2D;
   VECTOR_R2(v,XX_vF,YY_vF); DEFORM_AND_CHOOSE_CORNERS2D;
}

/* Rectangle enclosing ----------------------------------------------------- */
void box_enclosing(const matrix1D<double> &v0, const matrix1D<double> &vF,
    const matrix2D<double> &V, matrix1D<double> &corner1,
    matrix1D<double> &corner2) {
   SPEED_UP_temps;
   matrix1D<double> v(3);
   corner1.resize(3);
   corner2.resize(3);
   
   // Store values for reusing input as output vectors
   double XX_v0=XX(v0);
   double YY_v0=YY(v0);
   double ZZ_v0=ZZ(v0);
   double XX_vF=XX(vF);
   double YY_vF=YY(vF);
   double ZZ_vF=ZZ(vF);

   VECTOR_R3(v,XX_v0,YY_v0,ZZ_v0);
   M3x3_BY_V3x1(v,V,v);
   XX(corner1)=XX(v); XX(corner2)=XX(v);
   YY(corner1)=YY(v); YY(corner2)=YY(v);
   ZZ(corner1)=ZZ(v); ZZ(corner2)=ZZ(v);

   #define DEFORM_AND_CHOOSE_CORNERS3D \
      M3x3_BY_V3x1(v,V,v); \
      XX(corner1)=MIN(XX(corner1),XX(v)); XX(corner2)=MAX(XX(corner2),XX(v)); \
      YY(corner1)=MIN(YY(corner1),YY(v)); YY(corner2)=MAX(YY(corner2),YY(v)); \
      ZZ(corner1)=MIN(ZZ(corner1),ZZ(v)); ZZ(corner2)=MAX(ZZ(corner2),ZZ(v));

   VECTOR_R3(v,XX_vF,YY_v0,ZZ_v0); DEFORM_AND_CHOOSE_CORNERS3D;
   VECTOR_R3(v,XX_v0,YY_vF,ZZ_v0); DEFORM_AND_CHOOSE_CORNERS3D;
   VECTOR_R3(v,XX_vF,YY_vF,ZZ_v0); DEFORM_AND_CHOOSE_CORNERS3D;
   VECTOR_R3(v,XX_v0,YY_v0,ZZ_vF); DEFORM_AND_CHOOSE_CORNERS3D;
   VECTOR_R3(v,XX_vF,YY_v0,ZZ_vF); DEFORM_AND_CHOOSE_CORNERS3D;
   VECTOR_R3(v,XX_v0,YY_vF,ZZ_vF); DEFORM_AND_CHOOSE_CORNERS3D;
   VECTOR_R3(v,XX_vF,YY_vF,ZZ_vF); DEFORM_AND_CHOOSE_CORNERS3D;
}

/* ######################################################################### */
/* Euler Operations                                                          */
/* ######################################################################### */

/* Euler angles --> matrix ------------------------------------------------- */
void Euler_angles2matrix(double alpha, double beta, double gamma,
   matrix2D<double> &A) {
   double ca, sa, cb, sb, cg, sg;
   double cc, cs, sc, ss;
   
   A.resize(3,3);
   alpha = DEG2RAD(alpha);
   beta  = DEG2RAD(beta);
   gamma = DEG2RAD(gamma);
   
   ca = cos(alpha); cb = cos(beta); cg = cos(gamma);
   sa = sin(alpha); sb = sin(beta); sg = sin(gamma);
   cc = cb*ca; cs = cb*sa;
   sc = sb*ca; ss = sb*sa;
      
   A(0,0) =  cg*cc-sg*sa; A(0,1) =  cg*cs+sg*ca; A(0,2) = -cg*sb;
   A(1,0) = -sg*cc-cg*sa; A(1,1) = -sg*cs+cg*ca; A(1,2) = sg*sb;
   A(2,0) =  sc;          A(2,1) =  ss;          A(2,2) = cb;
}

/* Euler direction --------------------------------------------------------- */
void Euler_direction(double alpha, double beta, double gamma,
   matrix1D<double> &v) {
   double ca, sa, cb, sb;
   double cc, cs, sc, ss;
   
   v.resize(3);
   alpha = DEG2RAD(alpha);
   beta  = DEG2RAD(beta);
   
   ca = cos(alpha); cb = cos(beta);
   sa = sin(alpha); sb = sin(beta);
   sc = sb*ca; ss = sb*sa;
      
   VEC_ELEM(v,0)=sc;      VEC_ELEM(v,1)=ss;      VEC_ELEM(v,2)=cb;
}

/* Euler direction2angles ------------------------------- */
//gamma is useless but I keep it for simmetry
//with Euler_direction
void Euler_direction2angles(matrix1D<double> &v,
   double &alpha, double &beta, double &gamma) {
   double abs_ca, abs_sa, sb, cb;
   double aux_alpha;
   double aux_beta;
   double error, newerror;
   matrix1D<double> v_aux;

   v_aux.resize(3);
   cb=VEC_ELEM(v,2);

   if( fabs((cb))>0.999847695)/*one degree */
   {
   cout<< "\nWARMING: Routine Euler_direction2angles is not reliable\n"
          "for small tilt angles. Up to 0.001 deg it should be OK\n"
          "for most applications but you never know";
   }

   if( fabs((cb-1.))<FLT_EPSILON)
      {alpha=0.; beta=0.;}
   else{/*1*/

      aux_beta=acos(cb); /* beta between 0 and PI */


      if(aux_beta>PI/2) aux_beta -= PI; /* beta between -Pi/2 
                                       and PI/2 */
      if(aux_beta< (-PI/2)) aux_beta += PI; /* this should never 
                                happend but who knows*/
      sb = sin(aux_beta);

      abs_ca = fabs(VEC_ELEM(v,0))/sb;
      if(fabs((abs_ca-1.))<FLT_EPSILON)
        aux_alpha = 0.;
      else
        aux_alpha = acos(abs_ca);

      VEC_ELEM(v_aux,0)=sin(aux_beta)*cos(aux_alpha);
      VEC_ELEM(v_aux,1)=sin(aux_beta)*sin(aux_alpha);
      VEC_ELEM(v_aux,2)=cos(aux_beta);

      error = fabs(dot_product(v,v_aux)-1.); 
      alpha=aux_alpha; beta=aux_beta;
      
      VEC_ELEM(v_aux,0)=sin(beta)*cos(-1.*aux_alpha);
      VEC_ELEM(v_aux,1)=sin(beta)*sin(-1.*aux_alpha);
      VEC_ELEM(v_aux,2)=cos(beta);
      newerror = fabs(dot_product(v,v_aux)-1.); 
      if(error>newerror) {alpha = -1.*aux_alpha; 
                          beta  = aux_beta;
			  error = newerror;}

      VEC_ELEM(v_aux,0)=sin(-beta)*cos(-1.*aux_alpha);
      VEC_ELEM(v_aux,1)=sin(-beta)*sin(-1.*aux_alpha);
      VEC_ELEM(v_aux,2)=cos(-beta);
      newerror = fabs(dot_product(v,v_aux)-1.); 
      if(error>newerror) {alpha = -1.*aux_alpha; 
                          beta  = -1.*aux_beta;
			  error = newerror;}
      
      VEC_ELEM(v_aux,0)=sin(-beta)*cos(aux_alpha);
      VEC_ELEM(v_aux,1)=sin(-beta)*sin(aux_alpha);
      VEC_ELEM(v_aux,2)=cos(-beta);
      newerror = fabs(dot_product(v,v_aux)-1.); 
      
      if(error>newerror) {alpha = aux_alpha; 
                          beta  = -1.*aux_beta;
			  error = newerror;}
  }/*else 1 end*/
   gamma = 0.;
   beta  = RAD2DEG(beta);   
   alpha = RAD2DEG(alpha);
}/*Eulerdirection2angles end*/			   

/* Matrix --> Euler angles ------------------------------------------------- */
#define CHECK
//#define DEBUG
void Euler_matrix2angles(matrix2D<double> &A, double &alpha, double &beta,
   double &gamma) _THROW {
   double abs_sb, sign_sb;

   if (XSIZE(A)!=3 || YSIZE(A)!=3)
      REPORT_ERROR(1102,"Euler_matrix2angles: The Euler matrix is not 3x3");

   abs_sb = sqrt(A(0,2)*A(0,2)+A(1,2)*A(1,2));
   if (abs_sb > 16*FLT_EPSILON) {
      gamma = atan2(A(1,2),-A(0,2));
      alpha = atan2(A(2,1),A(2,0));
      if (ABS(sin(gamma))<FLT_EPSILON)
         sign_sb=SGN(-A(0,2)/cos(gamma));
         // if (sin(alpha)<FLT_EPSILON) sign_sb=SGN(-A(0,2)/cos(gamma));
         // else sign_sb=(sin(alpha)>0) ? SGN(A(2,1)):-SGN(A(2,1));
      else
         sign_sb = (sin(gamma)>0) ? SGN(A(1,2)):-SGN(A(1,2));
      beta  = atan2(sign_sb*abs_sb,A(2,2));
   } else { 
      if (SGN(A(2,2))>0) {
         // Let's consider the matrix as a rotation around Z
         alpha = 0;
         beta  = 0;
         gamma = atan2(-A(1,0),A(0,0));
      } else {
         alpha = 0;
         beta  = PI;
         gamma = atan2(A(1,0),-A(0,0));
      }
   }

   gamma = RAD2DEG(gamma);
   beta  = RAD2DEG(beta);
   alpha = RAD2DEG(alpha);

   #ifdef double
      matrix2D<double> Ap;
      Euler_angles2matrix(alpha,beta,gamma,Ap);
      if (A!=Ap) {
         cout << "---\n";
         cout << "Euler_matrix2angles: I have computed angles "
            " which doesn't match with the original matrix\n";
         cout << "Original matrix\n" << A;
         cout << "Computed angles alpha=" << alpha << " beta=" << beta 
              << " gamma=" << gamma << endl;
         cout << "New matrix\n" << Ap;
         cout << "---\n";
      }
   #endif

   #ifdef DEBUG
      cout << "abs_sb " << abs_sb << endl;
      cout << "A(1,2) " << A(1,2) << " A(0,2) " << A(0,2) << " gamma "
           << gamma << endl;
      cout << "A(2,1) " << A(2,1) << " A(2,0) " << A(2,0) << " alpha "
           << alpha << endl;
      cout << "sign sb " << sign_sb << " A(2,2) " << A(2,2)
           << " beta " << beta << endl;
   #endif
}
#undef CHECK
#undef DEBUG

#ifdef NEVERDEFINED
// Michael's method
void Euler_matrix2angles(matrix2D<double> A, double *alpha, double *beta,
   double *gamma) {
   double abs_sb;

   if      (ABS(A(1,1))>FLT_EPSILON) {
      abs_sb=sqrt((-A(2,2)*A(1,2)*A(2,1)-A(0,2)*A(2,0))/A(1,1));
   } else if (ABS(A(0,1))>FLT_EPSILON) {
      abs_sb=sqrt((-A(2,1)*A(2,2)*A(0,2)+A(2,0)*A(1,2))/A(0,1));
   } else if (ABS(A(0,0))>FLT_EPSILON) {
      abs_sb=sqrt((-A(2,0)*A(2,2)*A(0,2)-A(2,1)*A(1,2))/A(0,0));
   } else
      EXIT_ERROR(1,"Don't know how to extract angles");

   if (abs_sb>FLT_EPSILON) {    
      *beta  = atan2(abs_sb,A(2,2));
      *alpha = atan2(A(2,1)/abs_sb, A(2,0)/abs_sb);
      *gamma = atan2(A(1,2)/abs_sb,-A(0,2)/abs_sb);
   } else {
      *alpha = 0;
      *beta  = 0;
      *gamma = atan2(A(1,0),A(0,0));
   }
   
   *gamma = rad2deg(*gamma);
   *beta  = rad2deg(*beta);
   *alpha = rad2deg(*alpha);
}
#endif
void Euler_Angles_after_compresion(const double rot, double tilt, double psi,
   double &new_rot, double &new_tilt, double &new_psi,  matrix2D<double> &D)
{
   int i; 
   matrix1D<double> w(3);
   matrix1D<double> new_w(3);
   matrix2D<double> D_1(3,3);

   double module;
   double newrot, newtilt, newpsi;
   
   //if D has not inverse we are not in business   
   try{ 
     D_1=D.inv();
   } catch(Xmipp_error &XE){
     cout << XE;
     exit(1);
   }

   Euler_direction( rot, tilt, psi, w);
   if( fabs(VEC_ELEM(w,2))>0.999847695)/*cos one degree */
   {
   Euler_direction( rot, 10., psi, w);
   new_w =  (matrix1D<double>)(D_1*w)/((D_1*w).module());
   Euler_direction2angles (new_w, new_rot, new_tilt, new_psi);

   Euler_direction( rot, tilt, psi, w);
   new_w =  (matrix1D<double>)((D_1*w)/((D_1*w).module()));
   new_tilt= SGN(new_tilt)*fabs(ACOSD(VEC_ELEM(new_w,2)));
   new_psi = psi;
   
   // so, for small tilt the value of the rot is not realiable
   // doubleo overcome this problem I first calculate the rot for 
   // any arbitrary large tilt angle and the right rotation
   // and then I calculate the new tilt.
   // Please notice that the new_rotation is not a funcion of
   // the old tilt angle so I can use any arbitrary tilt angle
   }
   else
   {
   new_w =  (matrix1D<double>)(D_1*w)/((D_1*w).module());
   Euler_direction2angles (new_w, new_rot, new_tilt, new_psi);
   new_psi = psi;
   }
}

/* Euler up-down correction ------------------------------------------------ */
void Euler_up_down(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi) {
   newrot  = rot;
   newtilt = tilt+180;
   newpsi  = -(180+psi);
}

/* Same view, differently expressed ---------------------------------------- */
void Euler_another_set(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi) {
   newrot  = rot+180;
   newtilt = -tilt;
   newpsi  = -180+psi;
}

/* Euler mirror Y ---------------------------------------------------------- */
void Euler_mirrorY(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi) {
   newrot  = rot;
   newtilt = tilt+180;
   newpsi  = -psi;
}

/* Euler mirror X ---------------------------------------------------------- */
void Euler_mirrorX(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi) {
   newrot  = rot;
   newtilt = tilt+180;
   newpsi  = 180-psi;
}

/* Euler mirror XY --------------------------------------------------------- */
void Euler_mirrorXY(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi) {
   newrot  = rot;
   newtilt = tilt;
   newpsi  = 180+psi;
}

/* Apply a transformation matrix to Euler angles --------------------------- */
void Euler_apply_transf (const matrix2D<double> &L, const matrix2D<double> &R,
  double rot,     double tilt, double psi,
  double &newrot, double &newtilt, double &newpsi) {

  matrix2D<double> euler(3,3),temp;
  Euler_angles2matrix(rot,tilt,psi,euler);
  temp=L*euler*R;
  Euler_matrix2angles(temp,newrot,newtilt,newpsi);
}

/* Rotate matrix3D with 3 Euler angles ------------------------------------- */
matrix2D<double> Euler_rot3D_matrix(double rot, double tilt, double psi) {
   matrix2D<double> temp;
   Euler_angles2matrix(rot,-tilt,psi,temp);
   temp.resize(4,4);
   DIRECT_MAT_ELEM(temp,3,3)=1;
   return temp;
}

void Euler_rotate(const matrix3D<double> &V, double rot, double tilt, double psi,
   matrix3D<double> &result) {
   apply_geom(result,Euler_rot3D_matrix(rot,tilt,psi),V,IS_NOT_INV,DONT_WRAP);
}

matrix3D<double> Euler_rotate(const matrix3D<double> &V,
   double rot, double tilt, double psi) {
   matrix3D<double> aux; Euler_rotate(V,rot,tilt,psi,aux); return aux;
}


/* ######################################################################### */
/* Intersections                                                             */
/* ######################################################################### */

/* Intersection with a unit sphere ----------------------------------------- */
double intersection_unit_sphere(
   const matrix1D<double> &u,     // direction
   const matrix1D<double> &r)     // passing point
{

   // Some useful constants
   double A=XX(u)*XX(u)+YY(u)*YY(u)+ZZ(u)*ZZ(u);
   double B=XX(r)*XX(u)+YY(r)*YY(u)+ZZ(r)*ZZ(u);
   double C=XX(r)*XX(r)+YY(r)*YY(r)+ZZ(r)*ZZ(r)-1.0;
   double B2_AC=B*B-A*C;
   
   // A degenerate case?
   if (A==0) {
      if (B==0) return -1; // The ellipsoid doesn't intersect
      return 0;            // The ellipsoid is tangent at t=-C/2B
   }
   if (B2_AC<0) return -1;
   
   // A normal intersection
   B2_AC=sqrt(B2_AC);
   double t1=(-B-B2_AC)/A;  // The two parameters within the line for
   double t2=(-B+B2_AC)/A;  // the solution
   return ABS(t2-t1);
}

/* Intersection with a unit cylinder --------------------------------------- */
double intersection_unit_cylinder(
   const matrix1D<double> &u,     // direction
   const matrix1D<double> &r)     // passing point
{
   // Intersect with an infinite cylinder of radius=ry
   double A=XX(u)*XX(u)+YY(u)*YY(u);
   double B=XX(r)*XX(u)+YY(r)*YY(u);
   double C=XX(r)*XX(r)+YY(r)*YY(r)-1;
   
   double B2_AC=B*B-A*C;
   if      (A==0)    {
      if (C>0) return 0;       // Paralell ray outside the cylinder
      else     return 1/ZZ(u); // return height
   } else if (B2_AC<0) return 0;
   B2_AC=sqrt(B2_AC);
   
   // Points at intersection
   double t1=(-B - B2_AC)/A;
   double t2=(-B + B2_AC)/A;
   double z1 = ZZ(r) + t1 * ZZ(u);
   double z2 = ZZ(r) + t2 * ZZ(u);
   
   // Check position of the intersecting points with respect to
   // the finite cylinder, if any is outside correct it to the
   // right place in the top or bottom of the cylinder
   if (ABS(z1)>=0.5) t1=(SGN(z1)*0.5-ZZ(r))/ZZ(u);
   if (ABS(z2)>=0.5) t2=(SGN(z2)*0.5-ZZ(r))/ZZ(u);
   
   return ABS(t1-t2);
}

/* Intersection with a unit cube ------------------------------------------- */
double intersection_unit_cube(
   const matrix1D<double> &u,     // direction
   const matrix1D<double> &r)     // passing point
{
   double t1, t2, t;
   int found_t=0;

   #define ASSIGN_IF_GOOD_ONE \
      if (ABS(XX(r)+t*XX(u))-XMIPP_EQUAL_ACCURACY<=0.5 && \
          ABS(YY(r)+t*YY(u))-XMIPP_EQUAL_ACCURACY<=0.5 && \
          ABS(ZZ(r)+t*ZZ(u))-XMIPP_EQUAL_ACCURACY<=0.5) {\
          if      (found_t==0) {found_t++; t1=t;} \
          else if (found_t==1) {found_t++; t2=t;} \
      }

   // Intersect with x=0.5 and x=-0.5
   if (XX(u)!=0) {
      t=( 0.5-XX(r))/XX(u); ASSIGN_IF_GOOD_ONE;
      t=(-0.5-XX(r))/XX(u); ASSIGN_IF_GOOD_ONE;
   }

   // Intersect with y=0.5 and y=-0.5
   if (YY(u)!=0 && found_t!=2) {
      t=( 0.5-YY(r))/YY(u); ASSIGN_IF_GOOD_ONE;
      t=(-0.5-YY(r))/YY(u); ASSIGN_IF_GOOD_ONE;
   }

   // Intersect with z=0.5 and z=-0.5
   if (ZZ(u)!=0 && found_t!=2) {
      t=( 0.5-ZZ(r))/ZZ(u); ASSIGN_IF_GOOD_ONE;
      t=(-0.5-ZZ(r))/ZZ(u); ASSIGN_IF_GOOD_ONE;
   }
   
   if (found_t==2) return ABS(t1-t2);
   else            return 0;
}
