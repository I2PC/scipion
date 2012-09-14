/* author: Carlos Oscar Sorzano
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
/* Adjusting --------------------------------------------------------------- */
Volume  vol_phantom_avg0;           // Phantom without average
Volume  vol_recons_avg0;            // Reconstruction with no average
Volume  scaled_vol;                 // Temporary volume
float   scale=-1;                   // Scale for the gray level adjusting
float   scale_m1, scale_p1;         // Variables for descendent gradient
//float   eFOM_m1, eFOM_p1;         // approach
float   fFOM_m1, fFOM_p1;           // approach
float   scale_step;                 // Step to move the scale
float   substracted_mean=-1;        // Mean to substract from the recons.

       } else if (!strcmp (argv[i], "-adjust")) {
	  adjust=1;
          if (i+1<argc) {
             int stat2=sscanf(argv[i+1],"%f",&scale);
             if (stat2==1) i++;
             else scale=-1;
             if (!(i+1<argc)) {
                cout << "The mean to substract is missing\n";
                Usage();
             }
             stat2=sscanf(argv[i+1],"%f",&substracted_mean);
             if (stat2==1) i++;
             else {
                cout << "The mean to substract is missing\n";
                Usage();
             }

/* Adjusting? ============================================================== */
   if (adjust) {
      float sum_cross=0, sum_sq=0;
      int num_voxels=0;
      float diff=0;
      for (k=vol_label().startingZ(); k<=vol_label().finishingZ(); k++)
         for (i=vol_label().startingY(); i<=vol_label().finishingY(); i++)
            for (j=vol_label().startingX(); j<=vol_label().finishingX(); j++) {
               if (vol_label(k,i,j)!=0) {
                  diff += vol_recons(k,i,j) - vol_phantom(k,i,j);
                  num_voxels++;
               }
            }
      float m=diff/num_voxels;
      cout << "m " << m << endl;

      for (k=vol_label().startingZ(); k<=vol_label().finishingZ(); k++)
         for (i=vol_label().startingY(); i<=vol_label().finishingY(); i++)
            for (j=vol_label().startingX(); j<=vol_label().finishingX(); j++) {
               if (vol_label(k,i,j)!=0) {
                  sum_cross    += (vol_phantom(k,i,j)-m)*vol_recons(k,i,j);
                  sum_sq       += (vol_recons(k,i,j)-m)*(vol_recons(k,i,j)-m);
               }
            }
      vol_phantom().compute_stats();
      vol_recons().compute_stats();
      cout << "Before adjusting\n";
      cout << "Phantom Stats: "; vol_phantom().print_stats(); cout << endl;
      cout << "Recons  Stats: "; vol_recons().print_stats(); cout << endl;
      cout << "lambda = " << sum_cross/sum_sq << endl;
      vol_recons()=(vol_recons()-m)*(sum_cross/sum_sq);
      vol_phantom().compute_stats();
      vol_recons().compute_stats();
      cout << "After adjusting\n";
      cout << "Phantom Stats: "; vol_phantom().print_stats(); cout << endl;
      cout << "Recons  Stats: "; vol_recons().print_stats(); cout << endl;
      vol_phantom().compute_stats();
      vol_recons().compute_stats();
      vol_recons.write("inter_global");
      vol_phantom.write("inter_phantom");
   }
#ifdef NUNCA
      float sum_cross=0, sum_sq=0;
      float num_voxels=0;
      float mean_phantom=0, mean_recons=0;
      for (k=vol_label().startingZ(); k<=vol_label().finishingZ(); k++)
         for (i=vol_label().startingY(); i<=vol_label().finishingY(); i++)
            for (j=vol_label().startingX(); j<=vol_label().finishingX(); j++) {
               if (vol_label(k,i,j)!=0) {
                  mean_phantom += vol_phantom(k,i,j);
                  mean_recons  += vol_recons(k,i,j);
                  num_voxels++;
               }
            }
      mean_phantom /= num_voxels;
      mean_recons  /= num_voxels;
      cout << "mean phantom " << mean_phantom << endl;
      cout << "mean recons " << mean_recons << endl;

      float diff = (mean_recons-mean_phantom);
      for (k=vol_label().startingZ(); k<=vol_label().finishingZ(); k++)
         for (i=vol_label().startingY(); i<=vol_label().finishingY(); i++)
            for (j=vol_label().startingX(); j<=vol_label().finishingX(); j++) {
               if (vol_label(k,i,j)!=0) {
                  printf("%d,%d,%d ---> phantom = %f, recons= %f\n",
                     k,i,j,vol_phantom(k,i,j)-diff,vol_recons(k,i,j));
                  sum_cross    += (vol_phantom(k,i,j)-diff)*vol_recons(k,i,j);
                  sum_sq       += (vol_recons(k,i,j)-diff)*
                     (vol_recons(k,i,j)-diff);
               }
            }
      vol_phantom().compute_stats();
      vol_recons().compute_stats();
      cout << "Before adjusting\n";
      cout << "Phantom Stats: "; vol_phantom().print_stats(); cout << endl;
      cout << "Recons  Stats: "; vol_recons().print_stats(); cout << endl;
      cout << "lambda = " << sum_cross/sum_sq << endl;
      vol_recons()=(vol_recons()-diff)*(sum_cross/sum_sq);
      vol_phantom().compute_stats();
      vol_recons().compute_stats();
      cout << "After adjusting\n";
      cout << "Phantom Stats: "; vol_phantom().print_stats(); cout << endl;
      cout << "Recons  Stats: "; vol_recons().print_stats(); cout << endl;
      vol_phantom().compute_stats();
      vol_recons().compute_stats();
      vol_recons.write("inter_global");
      vol_phantom.write("inter_phantom");
   }
   
// Otro metodo --------------------------
   if (adjust) {
      // The reconstruction is tried to be adjusted to the phantom
      // Set the zero mean for both .....................................
      vol_phantom().compute_stats();
      vol_recons().compute_stats();
      cout << "Before adjusting\n";
      cout << "Phantom Stats: "; vol_phantom().print_stats(); cout << endl;
      cout << "Recons  Stats: "; vol_recons().print_stats(); cout << endl;
      vol_recons_avg0()=vol_recons()-vol_recons().avg();
      vol_phantom_avg0()=vol_phantom()-vol_phantom().avg();

      // Now adjust scale factor ........................................
      // The first trial is a scale factor which equalizes deviations
      scale=vol_phantom_avg0().stddev()/vol_recons_avg0().stddev();
      scaled_vol()=vol_recons_avg0()*scale;
      fFOM=compute_feature_FOM(&vol_phantom_avg0,&scaled_vol, &vol_label,
         feature_error);

      i=0;
      scale_step=scale/5;
      while (i<5) {
         cout << "Scale: " << scale << " (step= " << scale_step << 
            ") ---> fFOM=" << fFOM << endl;

         scale_m1=scale-scale_step;
         scaled_vol()=vol_recons_avg0()*scale_m1;
         fFOM_m1=compute_feature_FOM(&vol_phantom_avg0,&scaled_vol, &vol_label,
            feature_error);

         scale_p1=scale+scale_step;
         scaled_vol()=vol_recons_avg0()*scale_p1;
         fFOM_p1=compute_feature_FOM(&vol_phantom_avg0,&scaled_vol, &vol_label,
            feature_error);

         if      (fFOM_p1>fFOM) {scale=scale_p1; fFOM=fFOM_p1;}
         else if (fFOM_m1>fFOM) {scale=scale_m1; fFOM=fFOM_m1;}
         else                   {scale_step /=2; i++;}

      }

      // Scaling the reconstruction volume
      vol_recons()=scale*vol_recons_avg0() + vol_phantom().avg();
   }


// Otra forma -------------------------
   if (adjust) {
      if (scale==-1) {
         // The reconstruction is tried to be adjusted to the phantom
         // Set the zero mean for both .....................................
         vol_phantom().compute_stats();
         vol_recons().compute_stats();
         cout << "Before adjusting\n";
         cout << "Phantom Stats: "; vol_phantom().print_stats(); cout << endl;
         cout << "Recons  Stats: "; vol_recons().print_stats(); cout << endl;
         vol_recons_avg0()=vol_recons()-vol_recons().avg();
         vol_phantom_avg0()=vol_phantom()-vol_phantom().avg();

         // Now adjust scale factor ........................................
         // The first trial is a scale factor which equalizes deviations
         scale=vol_phantom_avg0().stddev()/vol_recons_avg0().stddev();
         scaled_vol()=vol_recons_avg0()*scale;
         eFOM=compute_error_FOM(&vol_phantom_avg0,&scaled_vol, error_radius);

         i=0;
         scale_step=scale/5;
         while (i<5) {
            cout << "Scale: " << scale << " (step= " << scale_step << 
               ") ---> eFOM=" << eFOM << endl;

            scale_m1=scale-scale_step;
            scaled_vol()=vol_recons_avg0()*scale_m1;
            eFOM_m1=compute_error_FOM(&vol_phantom_avg0,&scaled_vol, error_radius);

            scale_p1=scale+scale_step;
            scaled_vol()=vol_recons_avg0()*scale_p1;
            eFOM_p1=compute_error_FOM(&vol_phantom_avg0,&scaled_vol, error_radius);

            if      (eFOM_p1>eFOM) {scale=scale_p1; eFOM=eFOM_p1;}
            else if (eFOM_m1>eFOM) {scale=scale_m1; eFOM=eFOM_m1;}
            else                   {scale_step /=2; i++;}

         }
      // Scaling the reconstruction volume
      vol_recons()=scale*vol_recons_avg0() + vol_phantom().avg();
      } else {
         // The adjusting parameters have been given in the command line
         vol_recons().compute_stats();
         vol_recons().print_stats(); cout << endl;
         vol_recons()=scale*(vol_recons()-substracted_mean)+vol_phantom().avg();
         vol_recons().compute_stats();
         vol_recons().print_stats(); cout << endl;
      }
   }
#endif
   printf("[-adjust [<scale> <mean>]: adjust, and I give or no some parameters\n");

