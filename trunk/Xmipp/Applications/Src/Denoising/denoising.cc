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

#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippWavelets.hh>
#include <XmippData/xmippHistograms.hh>
#include <XmippData/xmippFilters.hh>

#define REMOVE_SCALE      0
#define SOFT_THRESHOLDING 1
#define BAYESIAN          2
#define ADAPTIVE_SOFT     3
#define CENTRAL           4
#define SHAH              5

class Denoising_parameters: public Prog_parameters {
public:
   string DWT_type;
   int    denoising_type;
   int    scale;
   double threshold;
   int    R;
   int    Shah_outer;
   int    Shah_inner;
   int    Shah_refinement;
   matrix1D<double> Shah_weight;
   bool   Shah_edge;
public:
   void read(int argc, char **argv) _THROW {
      Prog_parameters::read(argc,argv);
      DWT_type=get_param(argc,argv,"-type","DAUB12");
      if      (DWT_type=="DAUB4")  set_DWT_type(DAUB4);
      else if (DWT_type=="DAUB12") set_DWT_type(DAUB12);
      else if (DWT_type=="DAUB20") set_DWT_type(DAUB20);
      else REPORT_ERROR(1,"DWT::read: Unknown DWT type");
      string aux=get_param(argc,argv,"-denoising","remove_scales");
      if      (aux=="remove_scales") denoising_type=REMOVE_SCALE;
      else if (aux=="soft_thresholding") denoising_type=SOFT_THRESHOLDING;
      #ifdef NEVER_DEFINED
      else if (aux=="bayesian")          denoising_type=BAYESIAN;
      #endif
      else if (aux=="adaptive_soft")     denoising_type=ADAPTIVE_SOFT;
      else if (aux=="central")           denoising_type=CENTRAL;
      else if (aux=="difussion")         denoising_type=SHAH;
      else                               denoising_type=REMOVE_SCALE;
      scale=AtoI(get_param(argc,argv,"-scale","0"));
      threshold=AtoF(get_param(argc,argv,"-th","50"));
      R=AtoI(get_param(argc,argv,"-R","-1"));
      Shah_outer=AtoI(get_param(argc,argv,"-outer","10"));
      Shah_inner=AtoI(get_param(argc,argv,"-inner","1"));
      Shah_refinement=AtoI(get_param(argc,argv,"-refinement","1"));
      if (check_param(argc,argv,"-Shah_weight"))
         Shah_weight=get_vector_param(argc,argv,"-Shah_weight",4);
      else {
         Shah_weight.init_zeros(4);
	 Shah_weight(1)=Shah_weight(2)=50; Shah_weight(3)=0.02;
      }
      Shah_edge=check_param(argc,argv,"-only_edge");
   }
   
   void show() {
      Prog_parameters::show();
      if (denoising_type!=SHAH) cout << "DWT type: " << DWT_type << endl;
      cout << "Denoising: ";
      switch (denoising_type) {
         case REMOVE_SCALE:
            cout << " Remove scale " << scale << endl;
            break;
         case SOFT_THRESHOLDING:
            cout << " Soft thresholding " << threshold << endl;
            break;
	 #ifdef NEVER_DEFINED
	 case BAYESIAN:
	    cout << " Bayesian\n";
	    break;
	 #endif
	 case ADAPTIVE_SOFT:
	    cout << " Adaptive soft thresholding\n";
	    break;
	 case CENTRAL:
	    cout << " Keeping central part " << R << " pixels\n";
	    break;
	 case SHAH:
	    cout << " Shah difussion\n"
	         << " Outer iterations " << Shah_outer << endl
	         << " Inner iterations " << Shah_inner << endl
		 << " Refinement interations " << Shah_refinement << endl
		 << " Weight " << Shah_weight.transpose() << endl;
	    if (Shah_edge) cout << " Generating edge image\n";
	    break;
      }
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "  [-type <str=DAUB12]        : DWT type. Valid types are:\n"
           << "                               DAUB4, DAUB12, DAUB20\n"
           << "  [-denoising <str=remove_scale] : Denoising action\n"
           << "                               remove_scale\n"
      	   #ifdef NEVER_DEFINED
	   << "                               bayesian\n"
	   #endif
           << "                               soft_thresholding\n"
	   << "                               adaptive_soft\n"
	   << "                               central\n"
	   << "                               diffusion\n"
           << "  [-scale <s=0>]             : scale to remove\n"
           << "  [-threshold <th=50>]       : threshold of values (%) to remove\n"
	   << "  [-R <r=-1>]                : Radius to keep, by default half the size\n"
	   << "  [-outer <it=10>]           : Difussion outer iterations\n"
	   << "  [-innter <it=1>]           : Difussion inner iterations\n"
	   << "  [-refinement <it=1>]       : Difussion refinement iterations\n"
	   << "  [-Shah_weight [w0,w1,w2,w3]]:Diffusion weights\n"
	   << "                               w0=data matching (=0)\n"
	   << "                               w1=1st derivative smooth (=50)\n"
	   << "                               w2=edge strength (=50)\n"
	   << "                               w3=edge smoothness (=0.02)\n"
	   << "  [-only_edge]               : Produce the edge image of the diffusion\n"
      ;
   }
};

void process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Denoising_parameters *eprm=(Denoising_parameters *) prm;
   if (eprm->denoising_type!=SHAH) DWT(img(),img());
   double th;
   matrix2D<double> surface_strength, edge_strength;
   histogram1D hist;
   switch (eprm->denoising_type) {
      case REMOVE_SCALE:
         clean_quadrant(img(),eprm->scale, "01");
         clean_quadrant(img(),eprm->scale, "10");
         clean_quadrant(img(),eprm->scale, "11");
         break;
      case SOFT_THRESHOLDING:
         compute_hist(img(),hist,100);
         soft_thresholding(img(),hist.percentil(eprm->threshold));
         break;
      #ifdef NEVER_DEFINED
      case BAYESIAN:
         bayesian_wiener_filtering(img(),eprm->scale);
	 break;
      #endif
      case ADAPTIVE_SOFT:
      	 adaptive_soft_thresholding(img(),eprm->scale);
	 break;
      case CENTRAL:
         DWT_keep_central_part(img(),eprm->R);
	 break;
      case SHAH:
	 Smoothing_Shah(img(), surface_strength, edge_strength,
	    eprm->Shah_weight, eprm->Shah_outer, eprm->Shah_inner,
	    eprm->Shah_refinement);
	 if (eprm->Shah_edge) img()=edge_strength;
	 else                 img()=surface_strength;
         break;
   }
   if (eprm->denoising_type!=SHAH) IDWT(img(),img());
}

void process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Denoising_parameters *eprm=(Denoising_parameters *) prm;
   DWT(vol(),vol());
   double th;
   histogram1D hist;
   switch (eprm->denoising_type) {
      case REMOVE_SCALE:
         clean_quadrant(vol(),eprm->scale, "001");
         clean_quadrant(vol(),eprm->scale, "010");
         clean_quadrant(vol(),eprm->scale, "011");
         clean_quadrant(vol(),eprm->scale, "100");
         clean_quadrant(vol(),eprm->scale, "101");
         clean_quadrant(vol(),eprm->scale, "110");
         clean_quadrant(vol(),eprm->scale, "111");
         break;
      case SOFT_THRESHOLDING:
         compute_hist(vol(),hist,100);
         soft_thresholding(vol(),hist.percentil(eprm->threshold));
         break;
      case ADAPTIVE_SOFT:
      	 cout << "Adaptive soft-thresholding not implemented for volumes\n";
	 break;
      case CENTRAL:
      	 cout << "Keep central part not implemented for volumes\n";
	 break;
      case SHAH:
         cout << "Shah Difussion not implemented for volumes\n";
	 break;
   }
   IDWT(vol(),vol());
}

int main (int argc, char **argv) {
   Denoising_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
