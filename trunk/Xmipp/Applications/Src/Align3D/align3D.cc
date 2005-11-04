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

#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippFilters.hh>
#include <XmippData/xmippGeometry.hh>
#include <XmippData/xmippMasks.hh>

void Usage(const Mask_Params &m);

int main (int argc, char **argv) {
   FileName fn1, fn2;
   double   rot0, rotF, tilt0, tiltF, psi0, psiF;
   double   step_rot, step_tilt, step_psi;
   double   scale0, scaleF, step_scale;
   double   z0, zF, y0, yF, x0, xF, step_z, step_y, step_x;
   double   grey_scale0, grey_scaleF, step_grey;
   int      tell;
   bool     apply;
   Mask_Params mask(INT_MASK);
   bool     mask_enabled, mean_in_mask;
   
   #define COVARIANCE       1
   #define LEAST_SQUARES    2
   int alignment_method;
   
   #define MINIMUM          1
   #define MAXIMUM          2
   int optimize_criterion;
   
   // Get parameters =======================================================
   try {
      fn1=get_param(argc,argv,"-i1");
      fn2=get_param(argc,argv,"-i2");
      get_3_double_params(argc,argv,"-rot",rot0,rotF,step_rot,0,0,1);
      get_3_double_params(argc,argv,"-tilt",tilt0,tiltF,step_tilt,0,0,1);
      get_3_double_params(argc,argv,"-psi",psi0,psiF,step_psi,0,0,1);
      get_3_double_params(argc,argv,"-scale",scale0,scaleF,step_scale,1,1,1);
      get_3_double_params(argc,argv,"-grey_scale",grey_scale0,grey_scaleF,
         step_grey,1,1,1);
      get_3_double_params(argc,argv,"-z",z0,zF,step_z,0,0,1);
      get_3_double_params(argc,argv,"-y",y0,yF,step_y,0,0,1);
      get_3_double_params(argc,argv,"-x",x0,xF,step_x,0,0,1);
      
      mask_enabled=check_param(argc,argv,"-mask");
      mean_in_mask=check_param(argc,argv,"-mean_in_mask");
      if (mask_enabled)
         mask.read(argc,argv);
      
      if (step_rot  ==0) step_rot  =1;
      if (step_tilt ==0) step_tilt =1;
      if (step_psi  ==0) step_psi  =1;
      if (step_scale==0) step_scale=1;
      if (step_grey ==0) step_grey =1;
      if (step_z    ==0) step_z    =1;
      if (step_y    ==0) step_y    =1;
      if (step_x    ==0) step_x    =1;
      tell=check_param(argc,argv,"-show_fit");
      apply=check_param(argc,argv,"-apply");
      
      if (check_param(argc,argv,"-covariance"))
         {alignment_method=COVARIANCE; optimize_criterion=MAXIMUM;}
      else if (check_param(argc,argv,"-least_squares"))
         {alignment_method=LEAST_SQUARES; optimize_criterion=MINIMUM;}
      else {alignment_method=COVARIANCE; optimize_criterion=MAXIMUM;}
   } catch (Xmipp_error XE) {cout << XE; Usage(mask); exit(1);}
   
   // Main program =========================================================
   //#define DEBUG
   try {
      VolumeXmipp V1(fn1); V1().set_Xmipp_origin();
      VolumeXmipp V2(fn2); V2().set_Xmipp_origin();
      double mean_1=V1().compute_avg();
      double mean_2=V2().compute_avg();
      matrix3D<double> *V, Vaux;
      matrix1D<double> r(3), sc(3);
      matrix2D<double> A;

      // Initialize best_fit
      double best_rot, best_tilt, best_psi;
      double best_sc, best_grey, best_fit, fit;
      matrix1D<double> best_r(3);
      bool first=true;
      
      // Generate mask
      const matrix3D<int> *mask_ptr;
      if (mask_enabled) {
         mask.generate_3Dmask(V1());
	 mask_ptr=&(mask.get_binary_mask3D());
      } else mask_ptr=NULL;

      // Count number of iterations
      int times=1;
      if (!tell) {
         if (grey_scale0!=grey_scaleF)
            times *= FLOOR(1+(grey_scaleF-grey_scale0)/step_grey);
         if (rot0!=rotF) times *= FLOOR(1+(rotF-rot0)/step_rot);
         if (tilt0!=tiltF) times *= FLOOR(1+(tiltF-tilt0)/step_tilt);
         if (psi0!=psiF) times *= FLOOR(1+(psiF-psi0)/step_psi);
         if (scale0!=scaleF) times *= FLOOR(1+(scaleF-scale0)/step_scale);
         if (z0!=zF) times *= FLOOR(1+(zF-z0)/step_z);
         if (y0!=yF) times *= FLOOR(1+(yF-y0)/step_y);
         if (x0!=xF) times *= FLOOR(1+(xF-x0)/step_x);
         init_progress_bar(times);
      } else
         cout << "#Scale Z Y X rot tilt psi grey_factor fitness\n";
      
      // Iterate
      int itime=0;
      int step_time=CEIL((double)times/60.0);
      for (double grey=grey_scale0; grey<=grey_scaleF ; grey+=step_grey)
       for (double rot=rot0; rot<=rotF ; rot+=step_rot)
        for (double tilt=tilt0; tilt<=tiltF ; tilt+=step_tilt)
         for (double psi=psi0; psi<=psiF ; psi+=step_psi)
          for (XX(sc)=scale0; XX(sc)<=scaleF ; XX(sc)+=step_scale)
           for (ZZ(r)=z0; ZZ(r)<=zF ; ZZ(r)+=step_z)
            for (YY(r)=y0; YY(r)<=yF ; YY(r)+=step_y)
             for (XX(r)=x0; XX(r)<=xF ; XX(r)+=step_x) {
                // Rotate?
                if (rot!=0 || tilt!=0 || psi!=0) {
                   Euler_angles2matrix(rot,tilt,psi,A);
                   A.resize(4,4); A(3,3)=1;
                } else A.init_identity(4);

                // Translate?
                if (XX(r)!=0 || YY(r)!=0 || ZZ(r)!=0)
                   A=A*translation3D_matrix(r);

                // Scale?
                if (XX(sc)!=1) {
                   YY(sc)=ZZ(sc)=XX(sc);
                   A=A*scale3D_matrix(sc);
                }

                // Apply geometrical transformation
                if (!A.IsIdent()) {
                   apply_geom_Bspline(Vaux,A,V2(),3,IS_NOT_INV,WRAP);
                   V=&Vaux;
                } else V=&(V2());

                // Scale grey level?
                if (grey!=1) {
                   if (V==&Vaux) Vaux *= grey;
                   else {
                      array_by_scalar(V2(),grey,Vaux,'*');
                      V=&Vaux;
                   }
                   mean_2=Vaux.compute_avg();
                }

      	        // Only within mask?
		if (mean_in_mask && mask_enabled) {
		   double dummy;
		   compute_stats_within_binary_mask(*mask_ptr,
		      *V, dummy,dummy,mean_2,dummy);
		   compute_stats_within_binary_mask(*mask_ptr,
		      V1(),dummy,dummy,mean_1,dummy);
	        }

                // Correlate
		switch (alignment_method) {
		   case (COVARIANCE):
                      fit=correlation_index(V1(),*V,mask_ptr);
		      break;
		   case (LEAST_SQUARES):
		      fit=rms(V1(),*V,mask_ptr);
		      break;
	        }
		#ifdef DEBUG
		   VolumeXmipp save; save()=*V; save.write("PPPV.vol");
		   char c;
		   cout << "Press any key\n"; cin >> c;
		#endif

                // The best?
                if ((fit>best_fit && optimize_criterion==MAXIMUM) ||
		    (fit<best_fit && optimize_criterion==MINIMUM) ||
		     first) {
                   best_fit=fit; best_sc=XX(sc); best_r=r;
                   best_rot=rot; best_tilt=tilt; best_psi=psi;
		   best_grey=grey;
		   first=false;
                }

                // Show fit
                if (tell)
                   cout << XX(sc) << " " << r.transpose() << " "
                        << rot << " " << tilt << " " << psi << " "
			<< grey << " "
                        << fit << endl;
                else
                   if (++itime%step_time==0) progress_bar(itime);
             }
      if (!tell) progress_bar(times);
      if (!first)
	 cout << "The best correlation is for\n"
              << "Scale                  : " << best_sc << endl
              << "Translation (X,Y,Z)    : " << best_r.transpose() << endl
              << "Rotation (rot,tilt,psi): "
        	 << best_rot << " " << best_tilt << " " << best_psi << endl
	      << "Best grey factor       : " << best_grey << endl
              << "Fitness value          : " << best_fit << endl;
      if (apply) {
	 Euler_angles2matrix(best_rot,best_tilt,best_psi,A);
	 A.resize(4,4); A(3,3)=1;
	 A=A*translation3D_matrix(best_r);
	 A=A*scale3D_matrix(vector_R3(best_sc,best_sc,best_sc));
	 V2()*=best_grey;
         apply_geom_Bspline(Vaux,A,V2(),3,IS_NOT_INV,WRAP);
	 V2()=Vaux;
	 V2.write();
      }
   } catch (Xmipp_error XE) {cout << XE;}
}

void Usage(const Mask_Params &m) {
   cerr << "Purpose: Align two volumes varying orientation, position and scale\n";
   cerr << "Usage: align3D [options]\n"
        << "   -i1 <volume1>                           : the first volume to align\n"
        << "   -i2 <volume2>                           : the second one\n"
        << "  [-rot        <rot0>  <rotF>  <step_rot>  : in degrees\n"
        << "  [-tilt       <tilt0> <tiltF> <step_tilt> : in degrees\n"
        << "  [-psi        <psi0>  <psiF>  <step_psi>  : in degrees\n"
        << "  [-scale      <sc0>   <scF>   <step_sc>   : size scale margin\n"
        << "  [-grey_scale <sc0>   <scF>   <step_sc>   : grey scale margin\n"
        << "  [-z          <z0>    <zF>    <step_z>    : Z position in pixels\n"
        << "  [-y          <y0>    <yF>    <step_y>    : Y position in pixels\n"
        << "  [-x          <x0>    <xF>    <step_x>    : X position in pixels\n"
	<< "  [-show_fit]                              : Show fitness values\n"
	<< "  [-apply]                                 : Apply best movement to -i2\n"
	<< "  [-mean_in_mask]                          : Use the means within the mask\n"
	<< "  [-covariance]                            : Covariance fitness criterion\n"
	<< "  [-least_squares]                         : LS fitness criterion\n"
   ;
   m.usage(); cout << endl;
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Align3D {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Align3D/Help/align3D.html";
      help="Align two volumes";
      OPEN MENU menu_align3D;
      COMMAND LINES {
	+ usual: align3D -i1 $FILE1 -i2 $FILE2
                       [-rot        $ROT0  $ROTF  $STEP_ROT]
                       [-tilt       $TILT0 $TILTF $STEP_TILT]
                       [-psi        $PSI0  $PSIF  $STEP_PSI]
                       [-scale      $SC0   $SCF   $STEP_SC]
                       [-grey_scale $GREY0 $GREYF $STEP_GREY]
                       [-z          $Z0    $ZF    $STEP_Z]
                       [-y          $Y0    $YF    $STEP_Y]
                       [-x          $X0    $XF    $STEP_X]
                       [-show_fit] [-apply] [-mean_in_mask]
                       [$FIT_MEASURE]
      }
      PARAMETER DEFINITIONS {
        $FILE1 {
	   label="Input file 1";
	   type=file existing;
	}
        $FILE2 {
	   label="Input file 2";
	   type=file existing;
	}
        OPT(-rot) {label="Rotational angle";}
           $ROT0     {type=float; label="Initial"; by default=0;}
           $ROTF     {type=float; label="Final"; by default=0;}
           $STEP_ROT {type=float; label="Step"; by default=1;}
        OPT(-tilt) {label="Tilting angle";}
           $TILT0     {type=float; label="Initial"; by default=0;}
           $TILTF     {type=float; label="Final"; by default=0;}
           $STEP_TILT {type=float; label="Step"; by default=1;}
        OPT(-psi) {label="In-plane rotational angle";}
           $PSI0     {type=float; label="Initial"; by default=0;}
           $PSIF     {type=float; label="Final"; by default=0;}
           $STEP_PSI {type=float; label="Step"; by default=1;}
        OPT(-scale) {label="Size scale";}
           $SC0     {type=float; label="Initial"; by default=1;}
           $SCF     {type=float; label="Final"; by default=1;}
           $STEP_SC {type=float; label="Step"; by default=1;}
        OPT(-grey_scale) {label="Grey scale";}
           $GREY0     {type=float; label="Initial"; by default=1;}
           $GREYF     {type=float; label="Final"; by default=1;}
           $STEP_GREY {type=float; label="Step"; by default=1;}
        OPT(-z) {label="Z position";}
           $Z0     {type=float; label="Initial"; by default=0;}
           $ZF     {type=float; label="Final"; by default=0;}
           $STEP_Z {type=float; label="Step"; by default=1;}
        OPT(-y) {label="Y position";}
           $Y0     {type=float; label="Initial"; by default=0;}
           $YF     {type=float; label="Final"; by default=0;}
           $STEP_Y {type=float; label="Step"; by default=1;}
        OPT(-x) {label="X position";}
           $X0     {type=float; label="Initial"; by default=0;}
           $XF     {type=float; label="Final"; by default=0;}
           $STEP_X {type=float; label="Step"; by default=1;}
        OPT(-show_fit) {label="Show fitness values";}
        OPT(-apply) {label="Apply best fit to volume 2";}
        OPT(-mean_in_mask) {label="use the mean within a mask for covariance";}
        $FIT_MEASURE {
           label="Fitness method";
           type=Exclusion {
              "Covariance" {-covariance}
              "Least squares" {-least_squares}
           };
        }
      }
   }

   MENU menu_align3D {
      "I/O parameters"
      $FILE1
      $FILE2
      OPT($FIT_MEASURE)
      "Aligning parameters"
      OPT(-rot)
      OPT(-tilt)
      OPT(-psi)
      OPT(-scale)
      OPT(-grey_scale)
      OPT(-z)
      OPT(-y)
      OPT(-x)
      "Miscellaneaous"
      OPT(-show_fit)
      OPT(-apply)
      "Mask"
      OPT(-mean_in_mask)
   }
*/
