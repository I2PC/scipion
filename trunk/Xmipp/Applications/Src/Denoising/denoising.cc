/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2000 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippWavelets.hh>
#include <XmippData/xmippHistograms.hh>

#define REMOVE_SCALE      0
#define SOFT_THRESHOLDING 1
#define BAYESIAN          2
#define ADAPTIVE_SOFT     3
#define CENTRAL           4

class Denoising_parameters: public Prog_parameters {
public:
   string DWT_type;
   int    denoising_type;
   int    scale;
   double threshold;
   int    R;
public:
   void read(int argc, char **argv) _THROW {
      Prog_parameters::read(argc,argv);
      DWT_type=get_param(argc,argv,"-type","DAUB12");
      if      (DWT_type=="DAUB4")  set_DWT_type(DAUB4);
      else if (DWT_type=="DAUB12") set_DWT_type(DAUB12);
      else if (DWT_type=="DAUB20") set_DWT_type(DAUB20);
      else REPORT_ERROR(1,"DWT::read: Unknown DWT type");
      string aux=get_param(argc,argv,"-denoising","remove_scales");
      if (aux=="remove_scales") denoising_type=REMOVE_SCALE;
      else if (aux=="soft_thresholding") denoising_type=SOFT_THRESHOLDING;
      #ifdef NEVER_DEFINED
      else if (aux=="bayesian")          denoising_type=BAYESIAN;
      #endif
      else if (aux=="adaptive_soft")     denoising_type=ADAPTIVE_SOFT;
      else if (aux=="central")           denoising_type=CENTRAL;
      else                               denoising_type=REMOVE_SCALE;
      scale=AtoI(get_param(argc,argv,"-scale","0"));
      threshold=AtoF(get_param(argc,argv,"-th","50"));
      R=AtoI(get_param(argc,argv,"-R","-1"));
   }
   
   void show() {
      Prog_parameters::show();
      cout << "DWT type: " << DWT_type << endl;
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
           << "  [-scale <s=0>]             : scale to remove\n"
           << "  [-threshold <th=50>]       : threshold of values (%) to remove\n"
	   << "  [-R <r=-1>]                : Radius to keep, by default half the size\n"
      ;
   }
};

void process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Denoising_parameters *eprm=(Denoising_parameters *) prm;
   DWT(img(),img());
   double th;
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
   }
   IDWT(img(),img());
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
   }
   IDWT(vol(),vol());
}

int main (int argc, char **argv) {
   Denoising_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}
