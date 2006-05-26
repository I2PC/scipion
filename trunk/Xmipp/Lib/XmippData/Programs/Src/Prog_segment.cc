/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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

#include "../Prog_segment.hh"
#include "../../xmippArgs.hh"
#include "../../xmippMorphology.hh"
#include "../../xmippFilters.hh"

// Read arguments ==========================================================
void Prog_segment_prm::read(int argc, char **argv) {
   fn_vol=get_param(argc,argv,"-i");
   voxel_mass=AtoF(get_param(argc,argv,"-voxel_mass","-1"));
   dalton_mass=AtoF(get_param(argc,argv,"-dalton_mass","-1"));
   aa_mass=AtoF(get_param(argc,argv,"-aa_mass","-1"));
   sampling_rate=AtoF(get_param(argc,argv,"-sampling_rate","-1"));
   fn_mask=get_param(argc,argv,"-o","");
   en_threshold=check_param(argc,argv,"-threshold");
   if (en_threshold)
      threshold=AtoF(get_param(argc,argv,"-threshold"));
}

// Show ====================================================================
ostream & operator << (ostream &out, const Prog_segment_prm &prm) {
   out << "Input file   : " << prm.fn_vol        << endl
       << "Voxel mass   : " << prm.voxel_mass    << endl
       << "Dalton mass  : " << prm.dalton_mass   << endl
       << "AA mass      : " << prm.aa_mass       << endl
       << "Sampling rate: " << prm.sampling_rate << endl
       << "Output mask  : " << prm.fn_mask       << endl
       << "Enable thres.: " << prm.en_threshold  << endl
       << "Threshold    : " << prm.threshold     << endl
   ;
   return out;
}

// usage ===================================================================
void Prog_segment_prm::usage() const {
   cerr << "   -i <input volume>       : Volume to segment\n"
        << "  [-voxel_mass  <mass>  |  : Mass in voxels\n"
        << "   [-dalton_mass <mass> |  : Mass in daltons\n"
        << "    -aa_mass     <mass>]   : Mass in aminoacids\n"
        << "   -sampling_rate <Tm>]    : Sampling rate (A/pix)\n"
        << "  [-o <output mask=\"\">]    : Output mask\n"
        << "  [-threshold <th>]        : Apply a single threshold\n"
   ;
}

// Produce side information ================================================
void Prog_segment_prm::produce_side_info() {
   V.read(fn_vol);
   double sampling_rate3=sampling_rate*sampling_rate*sampling_rate;
   if (voxel_mass==-1 && !en_threshold) {
      if ((dalton_mass==-1 && aa_mass==-1) || sampling_rate==-1)
         REPORT_ERROR(1,"Prog_segment_prm: No way to compute voxel mass");
      if (dalton_mass!=-1) voxel_mass=dalton_mass*1.207/sampling_rate3;
      else                 voxel_mass=aa_mass*110*1.207/sampling_rate3;
   }
   cout << endl << "Derived voxel_mass=" << voxel_mass << endl;
}

// Count voxels ============================================================
// Segment with a given threshold and compute the number of voxels in the
// biggest piece
//#define DEBUG
double segment_threshold(const Volume *V_in, Volume *V_out,
   double threshold) {
   Volume aux;

   // Binarize input volume
   (*V_out)()=(*V_in)();
   (*V_out)().threshold("below",threshold,threshold);
   (*V_out)().binarize(threshold);
   
   #ifdef DEBUG
      cout << threshold << endl;
      VolumeXmipp save; save()=(*V_in)(); save.write("PPP0.vol");
      save()=(*V_out)(); save.write("PPP1.vol");
   #endif
   
   // Apply morphological opening to input volume
   aux().resize((*V_out)());
   opening3D((*V_out)(), aux(),18,0,1);
   closing3D(aux(), (*V_out)(),18,0,1);

   #ifdef DEBUG
      save()=(*V_out)(); save.write("PPP2.vol");
   #endif

   // Count the number of different objects
   int no_comp=label_volume((*V_out)(),aux());
   matrix1D<double> count(no_comp+1);
   FOR_ALL_ELEMENTS_IN_MATRIX3D(aux())
      count((int)aux(k,i,j))++;

   #ifdef DEBUG
      cout << count << endl << endl;
      cout << "Press any key\n";
      char c; cin >> c;
   #endif

   // Pick the maximum
   count(0)=0; // We don't want to pick the background
   int imax; count.max_index(imax);

   // Select the mask with only that piece
   FOR_ALL_ELEMENTS_IN_MATRIX3D((*V_out)())
      (*V_out)(k,i,j)=aux(k,i,j)==imax;

   return count(imax);
}

// Really segment ==========================================================
void Prog_segment_prm::segment(VolumeXmipp &mask) {
   double th_min, th_max, val_min, val_max;
   V().compute_double_minmax(val_min,val_max);
   th_min=val_min;
   th_max=val_max;
   double mass_min=MULTIDIM_SIZE(V());
   double mass_max=1;
   
   bool ok=false;
   double th_med;
   if (!en_threshold) {
      // Perform a bracketing search until the mass is 
      // within a 0.1% of the desired mass
      do {
         th_med=(th_min+th_max)*0.5;
         double mass_med=segment_threshold(&V,&mask,th_med);
         cout << "Threshold= " << th_med
              << " mass of the main piece= " << mass_med << endl;
         if (ABS(mass_med-voxel_mass)/voxel_mass<0.001) {ok=true; break;}
         if ((th_max-th_min)/(val_max-val_min)<0.0001) break;
         if (mass_med<voxel_mass) {th_max=th_med; mass_max=mass_med;}
         else                     {th_min=th_med; mass_min=mass_med;}
      } while (true);
   } else {
      // Perform a single thresholding
      double mass_med=segment_threshold(&V,&mask,threshold);
      cout << "Threshold= " << threshold
           << " mass of the main piece= " << mass_med << endl;
      ok=true;
   }
   
   // Save mask if necessary
   if (ok) segment_threshold(&V,&mask,th_med);
   else {
      cout << "An exact threshold has not been found using " << th_min
             << " that produces " << mass_min << " voxels\n";
      segment_threshold(&V,&mask,th_min);
   }
   mask.write(fn_mask);
}
