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

#include "../Prog_project_crystal.hh"
#include "../Prog_art_crystal.hh"
#include <XmippData/xmippArgs.hh>

/* Empty constructor ======================================================= */
Crystal_Projection_Parameters::Crystal_Projection_Parameters() {
   crystal_Xdim=0;
   crystal_Ydim=0;
   orthogonal=false;
   a.clear();
   b.clear();
   Nshift_avg=0;
   Nshift_dev=0;
   disappearing_th=0;
   DF_shift_bool=false;
   DF_shift.clear();
}

/* Read Crystal Projection Parameters ====================================== */
void Crystal_Projection_Parameters::read(FileName fn_crystal,double scale) {
   FILE    *fh_param;
   char    line[201];
   int     lineNo=0;
   char    *auxstr;

   if ((fh_param = fopen(fn_crystal.c_str(), "r")) == NULL)
      REPORT_ERROR(3005,
         (string)"Prog_Project_Parameters::read: There is a problem "
         "opening the file "+fn_crystal);

   while (fgets (line, 200,fh_param ) != NULL) {
      if (line[0]==0)    continue;
      if (line[0]=='#')  continue;
      if (line[0]=='\n') continue;
      switch(lineNo) {
         case 0:
            crystal_Xdim=AtoI(first_token(line),3007,
               "Prog_Project_Crystal::read: Error in Crystal X dimension");
            crystal_Ydim=AtoI(next_token(),3007,
               "Prog_Project_Crystal::read: Error in Crystal Y dimension");
            lineNo++;
	    crystal_Xdim = ROUND(scale*crystal_Xdim);
	    crystal_Ydim = ROUND(scale*crystal_Ydim);
            break;
         case 1:
            a.resize(3);
            XX(a)=scale*AtoF(first_token(line),3007,
               "Prog_Project_Crystal::read: Error in X component of a");
            YY(a)=scale*AtoF(next_token(),3007,
               "Prog_Project_Crystal::read: Error in Y component of a");
            ZZ(a)=0;
            lineNo++;
            break;
         case 2:
            b.resize(3);
            XX(b)=scale*AtoF(first_token(line),3007,
               "Prog_Project_Crystal::read: Error in X component of b");
            YY(b)=scale*AtoF(next_token(),3007,
               "Prog_Project_Crystal::read: Error in Y component of b");
            ZZ(b)=0;
            lineNo++;
            break;
         case 3:
           Nshift_dev=scale*AtoF(first_word(line),3007,
              "Prog_Project_Parameters::read: Error in magnitude shift noise");
           auxstr=next_token();
	   if (auxstr!=NULL)
              Nshift_avg=scale*AtoF(auxstr,3007,
                 "Prog_Project_Parameters::read: Error in magnitude shift bias");
	   else Nshift_avg=0;
           lineNo++;
           break;
         case 4:
            disappearing_th=AtoF(first_token(line),3007,
               "Prog_Project_Crystal::read: Error in disappearing threshold");
            lineNo++;
            break;
         case 5:
            orthogonal=(strcmp(first_token(line),"Yes")==0);
            lineNo++;
            break;
         case 6:
            // shift file
            fn_shift=first_word(line);
            if (strcmp(fn_shift.c_str(),"NULL"))
	        DF_shift_bool=true;
	    	
            lineNo++;
            break;
      } /* switch end */  
   } /* while end */
   if (lineNo!=7 && lineNo!=6)
      REPORT_ERROR(3007,(string)"Prog_Project_Crystal::read: I "
         "couldn't read all parameters from file " + fn_crystal);

   fclose(fh_param);
}

/* Write =================================================================== */
void Crystal_Projection_Parameters::write(FileName fn_crystal) {
   FILE *fh_param;

   if ((fh_param = fopen(fn_crystal.c_str(), "w")) == NULL)
      REPORT_ERROR(3005,
         (string)"Prog_Project_Parameters::write: There is a problem "
         "opening the file "+fn_crystal+" for output");

   fprintf(fh_param,"# Crystal dimensions (X, Y)\n");
   fprintf(fh_param,"%d %d\n",crystal_Xdim,crystal_Ydim);
   fprintf(fh_param,"# Crystal a lattice vector (X, Y)\n");
   fprintf(fh_param,"%f %f\n",XX(a),YY(a));
   fprintf(fh_param,"# Crystal b lattice vector (X, Y)\n");
   fprintf(fh_param,"%f %f\n",XX(b),YY(b));

   fprintf(fh_param,"#     Noise (and bias) applied to the magnitude shift\n");
   fprintf(fh_param,"%f ",Nshift_dev);
   if (Nshift_avg!=0) fprintf(fh_param,"%f \n",Nshift_avg);
   else fprintf(fh_param,"\n");
   
   fprintf(fh_param,"# Disappearing threshold\n");
   fprintf(fh_param,"%f\n",disappearing_th);
      
   fprintf(fh_param,"# Orthogonal Projections\n");
   if (orthogonal) fprintf(fh_param,"Yes\n");
   else            fprintf(fh_param,"No\n");
      
//   fprintf(fh_param,"# Grid relative size\n");

   fprintf(fh_param,"# File with shifts for each unit cell\n");
   fprintf(fh_param,"%s",fn_shift.c_str());
      
   fclose(fh_param);
}

/* Project crystal --------------------------------------------------------- */
//#define DEBUG
//#define DEBUG_MORE
void project_crystal(Phantom &phantom, Projection &P,
   const Projection_Parameters &prm,
   PROJECT_Side_Info &side, const Crystal_Projection_Parameters &prm_crystal,
   float rot, float tilt, float psi) {
   SPEED_UP_temps;
   //Scale crystal output
//   int aux_crystal_Ydim = ROUND(prm_crystal.crystal_Ydim
//   			   *phantom.phantom_scale);
//   int aux_crystal_Xdim = ROUND(prm_crystal.crystal_Xdim
//                                       *phantom.phantom_scale);
   
   // Initialize whole crystal projection
   P.adapt_to_size(prm_crystal.crystal_Ydim,prm_crystal.crystal_Xdim);
   
   // Compute lattice vectors in the projection plane
   P.set_angles(rot,tilt,psi);
   matrix1D<double> proja=P.euler*prm_crystal.a;
   matrix1D<double> projb=P.euler*prm_crystal.b;
   
   // Check if orthogonal projections
   // (projXdim,0)'=A*aproj
   // (0,projYdim)'=A*bproj
   matrix2D<double> Ainv, A,D,Dinv,AuxMat;
   if (prm_crystal.orthogonal) {
      A.resize(2,2);
      A(0,0)= YY(projb)*XSIZE(P()); A(0,1)=-XX(projb)*XSIZE(P());
      A(1,0)=-YY(proja)*YSIZE(P()); A(1,1)= XX(proja)*YSIZE(P());
      double nor=1/(XX(proja)*YY(projb)-XX(projb)*YY(proja));
      M2x2_BY_CT(A,A,nor);
      A.resize(3,3); A(2,2)=1;
      Ainv=A.inv();
      AuxMat.resize(3,3);
      //matrix with ceytal vectors
      AuxMat(0,0)= XX(prm_crystal.a); AuxMat(0,1)= YY(prm_crystal.a); AuxMat(0,2)= 0.;
      AuxMat(1,0)= XX(prm_crystal.b); AuxMat(1,1)= YY(prm_crystal.b); AuxMat(1,2)= 0.;
      AuxMat(2,0)= 0.               ; AuxMat(2,1)= 0.               ; AuxMat(2,2)= 1.;
      D.resize(3,3);
      D(0,0)= XSIZE(P()); D(0,1)= 0.;	      D(0,2)= 0.;
      D(1,0)= 0.;	  D(1,1)= YSIZE(P()); D(1,2)= 0.;
      D(2,0)= 0.;	  D(2,1)= 0.;	      D(2,2)= 1.;
      D.inv(D);
      M3x3_BY_M3x3(D,AuxMat,D);
      Dinv.resize(3,3);
      Dinv=D.inv();
   } else {
      A.init_identity(3);
      Ainv.init_identity(3);
   }
   //#define DEBUG
   #ifdef DEBUG
      cout << "P shape "; P().print_shape(); cout << endl;
      cout << "P.euler " << P.euler; cout << endl;
      cout << "rot= " << rot << " tilt= " << tilt << " psi= " << psi << endl;
      cout << "a " << prm_crystal.a.transpose() << endl;
      cout << "b " << prm_crystal.b.transpose() << endl;
      cout << "proja " << proja.transpose() << endl;
      cout << "projb " << projb.transpose() << endl;
      cout << "proj_Xdim " << prm.proj_Xdim << " proj_Ydim " << prm.proj_Ydim
           << endl;
      cout << "A\n" << A << endl;
      cout << "Ainv\n" << Ainv << endl;
      cout << "D\n" << D << endl;
      cout << "Dinv\n" << Dinv << endl;
   #endif

   // Compute aproj and bproj in the deformed projection space
   matrix1D<double> aprojd=A*proja;
   matrix1D<double> bprojd=A*projb;
   #ifdef DEBUG
      cout << "aprojd " << aprojd.transpose() << endl;
      cout << "bprojd " << bprojd.transpose() << endl;
   #endif

   // Get rid of all unnecessary components
   matrix2D<double> A2D=A; A2D.resize(2,2);
   proja.resize(2);
   projb.resize(2);
   aprojd.resize(2);
   bprojd.resize(2);

   // The unit cell projection size as well
   matrix1D<double> corner1(2), corner2(2);
   // First in the compressed space
   XX(corner1)=FIRST_XMIPP_INDEX(prm.proj_Xdim);
   YY(corner1)=FIRST_XMIPP_INDEX(prm.proj_Ydim);
   XX(corner2)= LAST_XMIPP_INDEX(prm.proj_Xdim);
   YY(corner2)= LAST_XMIPP_INDEX(prm.proj_Ydim);
   // Now deform
   #ifdef DEBUG
      cout << "corner1 before deformation " << corner1.transpose() << endl;
      cout << "corner2 before deformation " << corner2.transpose() << endl;
   #endif
   rectangle_enclosing(corner1,corner2,A2D,corner1,corner2);
   #ifdef DEBUG
      cout << "corner1 after deformation " << corner1.transpose() << endl;
      cout << "corner2 after deformation " << corner2.transpose() << endl;
   #endif

   matrix2D<double> cell_shiftX, cell_shiftY;
   matrix2D<int>    cell_inside;
   matrix2D<double> exp_shifts_matrix_X;
   matrix2D<double> exp_shifts_matrix_Y;

   fill_cell_positions(P, proja, projb, aprojd, bprojd, corner1, corner2,
      prm_crystal, cell_shiftX, cell_shiftY, cell_inside,
      exp_shifts_matrix_X, exp_shifts_matrix_Y);
      
   // Fill a table with all exp shifts 
   init_shift_matrix(prm_crystal,cell_inside, exp_shifts_matrix_X,
                                              exp_shifts_matrix_Y);
   // add the shifts to the already compute values
   FOR_ALL_ELEMENTS_IN_MATRIX2D(exp_shifts_matrix_X) {
	 // Add experimental shifts
	 cell_shiftX(i,j) += exp_shifts_matrix_X(i,j);
	 cell_shiftY(i,j) += exp_shifts_matrix_Y(i,j);
   }
  
   #ifdef DEBUG
      cout << "Cell inside shape "; cell_inside.print_shape(); cout << endl;
      cout << "Cell inside\n" << cell_inside << endl;
      cout << "Cell shiftX\n" << cell_shiftX << endl;
      cout << "Cell shiftY\n" << cell_shiftY << endl;
   #endif

   // Prepare matrices to go from uncompressed space to deformed projection
   matrix2D<double> AE=A*P.euler;   // From uncompressed to deformed
   matrix2D<double> AEinv=AE.inv(); // From deformed to uncompressed

   double density_factor=1.0;
   if (prm_crystal.orthogonal) {
   // Remember to compute de density factor
      matrix1D<double> projection_direction(3);
      (P.euler).getCol(2,projection_direction);
      projection_direction.self_transpose();
      density_factor=(projection_direction*Dinv).module();
      #ifdef DEBUG
      cout << "projection_direction" << projection_direction << endl;
      cout << "projection_direction*A" << projection_direction*A << endl;
      #endif
   }
   #ifdef DEBUG
      cout << "X proyectado=" << (AE*vector_R3(1.0,0.0,0.0)).transpose() << endl;
      cout << "Y proyectado=" << (AE*vector_R3(0.0,1.0,0.0)).transpose() << endl;
      cout << "P.euler_shape="<< endl;
      (P.euler).print_shape();
      cout << "P.euler="      << P.euler << endl;
      cout << "AE="           << AE << endl;
      cout << "AEinv="        << AEinv << endl;
      cout << "Ainv="            << Ainv << endl;
      cout << "density_factor="      << density_factor << endl;
   #endif

   // Project all cells
   FOR_ALL_ELEMENTS_IN_MATRIX2D(cell_inside)
//      if (cell_inside(i,j) && rnd_unif(0,1)<prm_crystal.disappearing_th) {
      if (cell_inside(i,j) ) {
         // Shift the phantom
         // Remind that displacements are defined in the deformed projection
         // that is why they have to be translated to the Universal
         // coordinate system
         matrix1D<double> cell_shift(3);
         VECTOR_R3(cell_shift,cell_shiftX(i,j),cell_shiftY(i,j),0.0f);
	 //SHIFT still pending
	 //cell_shift = cell_shift*phantom.phantom_scale;
         #ifdef DEBUG
            cout << "cell_shift on deformed projection plane "
                 << cell_shift.transpose() << endl;
         #endif
         M3x3_BY_V3x1(cell_shift, AEinv, cell_shift);
         #ifdef DEBUG
            cout << "cell_shift on real space "
                 << cell_shift.transpose() << endl;
         #endif
         
         Phantom aux;
         aux=phantom;
         // Now the cell shift is defined in the uncompressed space
         // Any phantom in the projection line will give the same shift
         // in the projection we are interested in the phantom whose center
         // is in the XYplane
         cell_shift=cell_shift-ZZ(cell_shift)/ZZ(P.direction)*P.direction;
         #ifdef DEBUG
            cout << "cell_shift after moving to ground "
                 << cell_shift.transpose() << endl;
         #endif
         aux.shift(XX(cell_shift),YY(cell_shift),ZZ(cell_shift));

         // Project this phantom
         aux.project_to(P, AE,prm_crystal.disappearing_th);
	 // Multiply by factor 
	 P()=P()*density_factor;
         #ifdef DEBUG_MORE
            cout << "After Projecting ...\n" << aux << endl;
            P.write("inter");
            cout << "Hit any key\n";
            char c; cin >> c;
         #endif
      }
}
#undef DEBUG
#undef DEBUG_MORE

/* Find crystal limits ----------------------------------------------------- */
#define MIN_MODULE 1e-2
void find_crystal_limits(
   const matrix1D<double> &proj_corner1, const matrix1D<double> &proj_corner2, 
   const matrix1D<double> &cell_corner1, const matrix1D<double> &cell_corner2,
   const matrix1D<double> &a, const matrix1D<double> &b,
   int &iamin, int &iamax, int &ibmin, int &ibmax) {
   if (a.module()<MIN_MODULE || b.module()<MIN_MODULE) 
      REPORT_ERROR(1,"find_crystal_limits: one of the lattice vectors is "
         "extremely small");
   
   // Compute area to cover
   double x0=XX(proj_corner1)+XX(cell_corner1);
   double y0=YY(proj_corner1)+YY(cell_corner1);
   double xF=XX(proj_corner2)+XX(cell_corner2);
   double yF=YY(proj_corner2)+YY(cell_corner2);

   // Initialize
   SPEED_UP_temps;
   matrix2D<double> A(2,2), Ainv(2,2);
   A.setCol(0,a); A.setCol(1,b);
   M2x2_INV(Ainv,A);

   // Now express each corner in the a,b coordinate system
   matrix1D<double> r(2);
   VECTOR_R2(r,x0,y0);
   M2x2_BY_V2x1(r,Ainv,r);
   iamin=FLOOR(XX(r)); iamax=CEIL(XX(r));
   ibmin=FLOOR(YY(r)); ibmax=CEIL(YY(r));
   
   #define CHANGE_COORDS_AND_CHOOSE_CORNERS2D \
      M2x2_BY_V2x1(r,Ainv,r); \
      iamin=MIN(FLOOR(XX(r)),iamin); iamax=MAX(CEIL(XX(r)),iamax); \
      ibmin=MIN(FLOOR(YY(r)),ibmin); ibmax=MAX(CEIL(YY(r)),ibmax);

   VECTOR_R2(r,x0,yF); CHANGE_COORDS_AND_CHOOSE_CORNERS2D;
   VECTOR_R2(r,xF,y0); CHANGE_COORDS_AND_CHOOSE_CORNERS2D;
   VECTOR_R2(r,xF,yF); CHANGE_COORDS_AND_CHOOSE_CORNERS2D;
}

/* Move following spiral --------------------------------------------------- */
/* Given a cell, we study the cells in the inmidiate surrounding, this gives
   a structure of the like
   
   xxx
   xXx
   xxx
   
   We will mean by . a non-visited cell, by x a visited one and by o the
   following cell we should visit  given this structure. To identify
   structures we will use the sum of visited cells in rows and columns. The
   following patterns define the spiral movement:
   
   .o.0 ...0 ...0 .xx2 .xx2 xx.2 xo.1 .o.0 ...0
   .x.1 ox.1 .xx2 .xx2 .xo1 xxo2 xx.2 xx.2 ox.1
   ...0 .x.1 .ox1 .o.0 ...0 ...0 ...0 xx.2 xx.2 etc
   010  020  012  022  021  220  210  220  120

   This function supposes that r is inside the matrix
*/
void move_following_spiral(matrix1D<double> &r, const matrix2D<int> &visited) {
int r1=0, r2=0, r3=0, c1=0, c2=0, c3=0;

#define x0 STARTINGX(visited)
#define y0 STARTINGY(visited)
#define xF FINISHINGX(visited)
#define yF FINISHINGY(visited)
#define i  (int)YY(r)
#define j  (int)XX(r)
// Compute row and column sums
if (i>y0 && j>x0) {r1 += visited(i-1,j-1); c1 += visited(i-1,j-1);}
if (i>y0        ) {r1 += visited(i-1,j  ); c2 += visited(i-1,j  );}
if (i>y0 && j<xF) {r1 += visited(i-1,j+1); c3 += visited(i-1,j+1);}
if (        j>x0) {r2 += visited(i  ,j-1); c1 += visited(i  ,j-1);}
                   r2 += visited(i  ,j  ); c2 += visited(i  ,j  );
if (        j<xF) {r2 += visited(i  ,j+1); c3 += visited(i  ,j+1);}
if (i<yF && j>x0) {r3 += visited(i+1,j-1); c1 += visited(i+1,j-1);}
if (i<yF        ) {r3 += visited(i+1,j  ); c2 += visited(i+1,j  );}
if (i<yF && j<xF) {r3 += visited(i+1,j+1); c3 += visited(i+1,j+1);}

#ifdef DEBUG
   cout << r1 << " " << r2 << " " << r3 << " "
        << c1 << " " << c2 << " " << c3 << endl;
#endif

// Decide where to go
if      (r1==0 && r2==1 && r3==0 && c1==0 && c2==1 && c3==0) {         YY(r)--;}
else if (r1==0 && r2==1 && r3==1 && c1==0 && c2==2 && c3==0) {XX(r)--;         }
else if (r1==0 && r2==2 && r3==1 && c1==0 && c2==1 && c3==2) {         YY(r)++;}
else if (r1==2 && r2==2 && r3==0 && c1==0 && c2==2 && c3==2) {         YY(r)++;}
else if (r1==2 && r2==1 && r3==0 && c1==0 && c2==2 && c3==1) {XX(r)++;         }
else if (r1==2 && r2==2 && r3==0 && c1==2 && c2==2 && c3==0) {XX(r)++;         }
else if (r1==1 && r2==2 && r3==0 && c1==2 && c2==1 && c3==0) {         YY(r)--;}
else if (r1==0 && r2==2 && r3==2 && c1==2 && c2==2 && c3==0) {         YY(r)--;}
else if (r1==0 && r2==1 && r3==2 && c1==1 && c2==2 && c3==0) {XX(r)--;         }
else if (r1==1 && r2==2 && r3==2 && c1==3 && c2==2 && c3==0) {         YY(r)--;}
else if (r1==0 && r2==2 && r3==3 && c1==1 && c2==2 && c3==2) {XX(r)--;         }
else if (r1==0 && r2==2 && r3==2 && c1==0 && c2==2 && c3==2) {XX(r)--;         }
else if (r1==2 && r2==2 && r3==1 && c1==0 && c2==2 && c3==3) {         YY(r)++;}
else if (r1==3 && r2==2 && r3==0 && c1==2 && c2==2 && c3==1) {XX(r)++;         }
}
#undef i
#undef j
#undef x0
#undef y0
#undef xF
#undef yF

/* Fill cell rotations and shifts ------------------------------------------ */
//#define DEBUG
void fill_cell_positions(Projection &P,
   matrix1D<double> &aproj,   matrix1D<double> &bproj,
   matrix1D<double> &aprojd,  matrix1D<double> &bprojd,
   matrix1D<double> &corner1, matrix1D<double> &corner2,
   const Crystal_Projection_Parameters &prm_crystal,
   matrix2D<double> &cell_shiftX,
   matrix2D<double> &cell_shiftY,
   matrix2D<int>    &cell_inside,
   matrix2D<double> &exp_shifts_matrix_X,
   matrix2D<double> &exp_shifts_matrix_Y) {

   // Compute crystal limits
   int iamin, iamax, ibmin, ibmax;
   find_crystal_limits(vector_R2(STARTINGX(P()),STARTINGY(P())),
      vector_R2(FINISHINGX(P()),FINISHINGY(P())),
      corner1,corner2,aprojd,bprojd,iamin,iamax,ibmin,ibmax);

   #ifdef DEBUG
      P().print_shape(); cout << endl;
      cout << "aprojd=" << aproj.transpose() << endl;
      cout << "bprojd=" << bproj.transpose() << endl;
      cout << iamin << " " << iamax << " " << ibmin << " " << ibmax << endl;
   #endif
   
   // Compute weight table in the undeformed space
   matrix2D<double> weight(3,3);
   STARTINGX(weight)=-1; STARTINGY(weight)=-1;
   weight(0,0)=0;
   weight(-1, 0)=weight( 1, 0)=1/aproj.module();
   weight( 0,-1)=weight( 0, 1)=1/bproj.module();
   weight(-1, 1)=weight( 1,-1)=1/(aproj-bproj).module();
   weight(-1,-1)=weight( 1, 1)=1/(aproj+bproj).module();

   // Resize all needed matrices
   cell_shiftX.init_zeros(ibmax-ibmin+1,iamax-iamin+1);
   STARTINGX(cell_shiftX)=iamin; STARTINGY(cell_shiftX)=ibmin;
   cell_shiftY.init_zeros(cell_shiftX);
   cell_inside.init_zeros(ibmax-ibmin+1,iamax-iamin+1);
   STARTINGX(cell_inside)=iamin; STARTINGY(cell_inside)=ibmin;

   // Visited has got one cell more than the rest in all directions
   matrix2D<int> visited;
   int visited_size=MAX(iamax-iamin+1,ibmax-ibmin+1)+2;
   visited.init_zeros(visited_size,visited_size);
   STARTINGX(visited)=iamin-(visited_size-(iamax-iamin+1)+1)/2;
   STARTINGY(visited)=ibmin-(visited_size-(ibmax-ibmin+1)+1)/2;
   #ifdef DEBUG
      cout << "weight=" << weight;
      cout << "cell_shiftX shape "; cell_shiftX.print_shape(); cout << endl;
      cout << "cell_inside shape "; cell_inside.print_shape(); cout << endl;
      cout << "visited shape "; visited.print_shape(); cout << endl;
   #endif

   // Perform the crystal deformation starting in the middle of the
   // crystal and going in spirals until the four corners are visited
   // we find one corner
   matrix1D<double> r(2),ri(2),sh(2);
   r.init_zeros();
   #define INDEX(r) (int)YY(r),(int)XX(r)
   while (!visited.isCorner(r)) {
      visited(INDEX(r))=true;
      #ifdef DEBUG
         cout << "   Visiting " << r.transpose() << endl;
      #endif

      // Weight is computed in the undeformed space
      double total_weight=0;
      double total_shiftX=0;
      double total_shiftY=0;
      for (YY(sh)=-1; YY(sh)<=1; YY(sh)++)
          for (XX(sh)=-1; XX(sh)<=1; XX(sh)++) {
              V2_PLUS_V2(ri,r,sh);
              if (!cell_shiftX.outside(ri)) {
                 total_weight += weight(INDEX(sh));
                 total_shiftX += weight(INDEX(sh))*cell_shiftX(INDEX(ri));
                 total_shiftY += weight(INDEX(sh))*cell_shiftY(INDEX(ri));
              }
          }
      if (total_weight==0) total_weight=1;
      if (!cell_shiftX.outside(r)) {
         cell_shiftX(INDEX(r))=rnd_gaus(prm_crystal.Nshift_avg,
            prm_crystal.Nshift_dev)+total_shiftX/total_weight;
         cell_shiftY(INDEX(r))=rnd_gaus(prm_crystal.Nshift_avg,
            prm_crystal.Nshift_dev)+total_shiftY/total_weight;
      }

      // Move to next position
      move_following_spiral(r,visited);
   }
   
   #ifdef DEBUG
      cout << "Cell shift X without absolute displacements" << cell_shiftX;
      cout << "Cell shift Y without absolute displacements" << cell_shiftY;
   #endif

   // The previous shifts are relative to the final position, now
   // express the real final position

   FOR_ALL_ELEMENTS_IN_MATRIX2D(visited) {
      if (!cell_shiftX.outside(i,j)) {
         // Move to final position
         cell_shiftX(i,j) += j*XX(aprojd)+i*XX(bprojd);
         cell_shiftY(i,j) += j*YY(aprojd)+i*YY(bprojd);
         	 
         // Check if there is intersection
         matrix1D<double> auxcorner1(2), auxcorner2(2);
         XX(auxcorner1)=XX(corner1)+cell_shiftX(i,j);
         YY(auxcorner1)=YY(corner1)+cell_shiftY(i,j);
         XX(auxcorner2)=XX(corner2)+cell_shiftX(i,j);
         YY(auxcorner2)=YY(corner2)+cell_shiftY(i,j);
         
         cell_inside(i,j)=P().intersects(auxcorner1,auxcorner2);
//TEMPORAL FIX FOR PHANTOM AS BIG AS THE WHOLE CRYSTAL
         cell_inside(i,j)=1;
//ROBERTO         
         #ifdef DEBUG
            cout << "(i,j)=(" << i << "," << j << ")\n";
            cout << "   Projection shape "; P().print_shape(); cout << endl;
            cout << "   AuxCorner1 " << auxcorner1.transpose() << endl
                 << "   Origin     " << cell_shiftX(i,j) << " "
                                     << cell_shiftY(i,j) << endl
                 << "   AuxCorner2 " << auxcorner2.transpose() << endl;
            cout << "   Inside= " << cell_inside(i,j) << endl;
         #endif

      }
   }
}

/* Fill aux matrix with experimental shifs to add to unit cell 
   projection********************************************************/
   
   void init_shift_matrix(const Crystal_Projection_Parameters &prm_crystal, 
                          matrix2D<int>    &cell_inside,
			  matrix2D<double> &exp_shifts_matrix_X,
 			  matrix2D<double> &exp_shifts_matrix_Y)
  {
   DocFile        aux_DF_shift;//crystal_param is cont
   aux_DF_shift=prm_crystal.DF_shift;
   exp_shifts_matrix_X.resize(cell_inside);  
   exp_shifts_matrix_X.init_zeros(); 
   exp_shifts_matrix_Y.resize(cell_inside);  
   exp_shifts_matrix_Y.init_zeros(); 
   
   //#define DEBUG2
   #ifdef DEBUG2
   cout << aux_DF_shift;
   cout << "exp_shifts_matrix_X shape" << endl;
   exp_shifts_matrix_X.print_shape();
   cout << endl;
   #endif
   #undef DEBUG2 
   //fill matrix with docfile data
   aux_DF_shift.go_first_data_line();
   int max_x, max_y, min_x , min_y;
   
   while (!aux_DF_shift.eof()) {
      //Check that we are not outside the matrix
      if (!exp_shifts_matrix_X.outside(ROUND(aux_DF_shift(1)),
                                      ROUND(aux_DF_shift(0)) 
				     )){
         exp_shifts_matrix_X(ROUND(aux_DF_shift(1)),ROUND(aux_DF_shift(0)))
                            =aux_DF_shift(4);
         exp_shifts_matrix_Y(ROUND(aux_DF_shift(1)),ROUND(aux_DF_shift(0)))
                            =aux_DF_shift(5);
      }		    
      aux_DF_shift.next_data_line();
   }
      //#define DEBUG2
      #ifdef DEBUG2
      cout << exp_shifts_matrix_X;
      #endif
      #undef DEBUG2 
   
   }   

