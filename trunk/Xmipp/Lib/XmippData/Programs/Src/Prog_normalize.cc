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

#include "../Prog_normalize.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippMicrograph.hh>
#include <XmippData/xmippMasks.hh>

/* Read -------------------------------------------------------------------- */
void Normalize_parameters::read(int argc, char **argv) {
   Prog_parameters::read(argc,argv);

   // Get normalizing method
   string aux;
   aux = get_param(argc,argv,"-method","NewXmipp");
   if      (aux=="OldXmipp")      normalizing_method=OLDXMIPP;
   else if (aux=="Near_OldXmipp") normalizing_method=NEAR_OLDXMIPP;
   else if (aux=="NewXmipp")      normalizing_method=NEWXMIPP;
   else if (aux=="NewXmipp2")     normalizing_method=NEWXMIPP2;
   else if (aux=="Michael")       normalizing_method=MICHAEL;
   else if (aux=="Random")        normalizing_method=RANDOM;
   else if (aux=="None")          normalizing_method=NONE;
   else if (aux=="Ramp")          normalizing_method=RAMP;
   else REPORT_ERROR(1,"Normalize: Unknown normalizing method");

   // Normalizing a volume
   normalizing_vol=check_param(argc,argv,"-vol");

   apply_geo=false;
   // Get background mask
   if (!normalizing_vol) {
      if (normalizing_method==NEWXMIPP || normalizing_method==NEWXMIPP2 || 
          normalizing_method==MICHAEL || normalizing_method==NEAR_OLDXMIPP ||
	  normalizing_method==RAMP) {
	enable_mask=check_param(argc,argv,"-mask");
	if (enable_mask) {
	  mask_prm.allowed_data_types=INT_MASK;
	  mask_prm.read(argc,argv);
	} else {
	  enable_mask=false;
	  int i=position_param(argc,argv,"-background");
	  if (i+2>=argc)
	    REPORT_ERROR(1,"Normalize: Not enough parameters after -background");
	  aux=argv[i+1];
	  r=AtoI(argv[i+2]);
	  
	  if (aux=="frame")       background_mode=FRAME;
	  else if (aux=="circle") background_mode=CIRCLE;
	  else
	    REPORT_ERROR(1,"Normalize: Unknown background mode");
	}
	produce_side_info();

	 // Default is to apply inverse transformation from image header to the mask
	 apply_geo=!check_param(argc,argv,"-dont_apply_geo");
      } else
         background_mode=NONE;

      if (normalizing_method==RANDOM) {
         int i=position_param(argc,argv,"-prm");
         if (i+4>=argc)
            REPORT_ERROR(1,"Normalize_parameters::read: Not enough parameters after -prm");
         a0=AtoF(argv[i+1]);
         aF=AtoF(argv[i+2]);
         b0=AtoF(argv[i+3]);
         bF=AtoF(argv[i+4]);
      }
   }
}
   
/* Produce side information ------------------------------------------------ */
void Normalize_parameters::produce_side_info() {
   int Zdim, Ydim, Xdim;
   get_input_size(Zdim, Ydim, Xdim);
   if (!normalizing_vol) {
      if (Zdim!=1)
         REPORT_ERROR(1,"Normalize: This program works only with images");

      if (!enable_mask) {
         bg_mask.resize(Ydim,Xdim); bg_mask.set_Xmipp_origin();
         switch (background_mode) {
            case FRAME:
               BinaryFrameMask(bg_mask,Xdim-2*r,Ydim-2*r,OUTSIDE_MASK);
	       break;
            case CIRCLE:
      	       BinaryCircularMask(bg_mask, r, OUTSIDE_MASK);
	       break;
         }
      } else {
         mask_prm.generate_2Dmask(Ydim, Xdim);
         bg_mask=mask_prm.imask2D;
      }
   }
}

/* Show -------------------------------------------------------------------- */
void Normalize_parameters::show() {
   Prog_parameters::show();
   cout << "Normalizing volumes: " << normalizing_vol << endl;
   if (!normalizing_vol) {
      cout << "Normalizing method: ";
      switch (normalizing_method) {
         case OLDXMIPP:      cout << "OldXmipp\n"; break;
         case NEAR_OLDXMIPP: cout << "Near_OldXmipp\n"; break;
         case NEWXMIPP:      cout << "NewXmipp\n"; break;
         case NEWXMIPP2:     cout << "NewXmipp2\n"; break;
         case MICHAEL:       cout << "Michael\n"; break;
         case NONE:          cout << "None\n"; break;
         case RAMP:          cout << "Ramp\n"; break;
	 case RANDOM:        cout << "Random a=[" << a0 << "," << aF << "], "
	                          << "b=[" << b0 << "," << bF << "]\n"; break;
      }
      if (normalizing_method==NEWXMIPP || normalizing_method==NEWXMIPP2 || 
          normalizing_method==NEAR_OLDXMIPP || normalizing_method==MICHAEL ||
	  normalizing_method==RAMP) {
         cout << "Background mode: ";
         switch (background_mode) {
            case NONE :  cout << "None\n"; break;
            case FRAME:  cout << "Frame, width=" << r << endl; 
                         cout << "Apply transformation to mask: " << apply_geo << endl; 
                         break;
            case CIRCLE: cout << "Circle, radius=" << r << endl; 
                         cout << "Apply transformation to mask: " << apply_geo << endl; 
                         break;
         }
      }
   }
   if (normalizing_method==NEWXMIPP && enable_mask) mask_prm.show();
}

/* Usage ------------------------------------------------------------------- */
void Normalize_parameters::usage() {
   Prog_parameters::usage();
   cerr << "NORMALIZATION OF VOLUMES\n"
        << "  [-vol]                    : Activate this mode\n";
   cerr << "NORMALIZATION OF IMAGES\n"
        << "  [-method <mth=NewXmipp>   : Normalizing method. Valid ones are:\n"
        << "                              OldXmipp, Near_OldXmipp, NewXmipp\n"
        << "                              NewXmipp2, Michael, None, Random, Ramp\n"
        << "                              Methods NewXmipp, Michael, Near_OldXmipp\n"
        << "                              and Ramp need a background mask:\n"
        << "  [-background frame <r>  | : Rectangular background of r pixels\n"
        << "   -background circle <r> | : Circular background outside radius=r\n"
        << "   -mask <options>]           Use an alternative type of background mask\n"
        << "                               (see xmipp_mask for options) \n"
	<< "  [-prm a0 aF b0 bF]        : Only in random mode. y=ax+b\n"
   ;
}

/* Apply geometric transformation to the mask ------------------------------ */
void Normalize_parameters::apply_geo_mask(ImageXmipp &img) {
  matrix2D<double> tmp;
  tmp.resize(bg_mask);
  type_cast(bg_mask,tmp);
  double outside=DIRECT_MAT_ELEM(tmp,0,0);
  // Instead of IS_INV for images use IS_NOT_INV for masks!
  tmp.self_apply_geom_Bspline(img.get_transformation_matrix(),3,
     IS_NOT_INV,DONT_WRAP,outside);
  type_cast(tmp,bg_mask);
}

/* Apply ------------------------------------------------------------------- */
void Normalize_parameters::apply(Image *img) {
   double a, b;
   switch (normalizing_method) {
      case OLDXMIPP:
         normalize_OldXmipp(img);
         break;
      case NEAR_OLDXMIPP:
         normalize_Near_OldXmipp(img,bg_mask);
         break;
      case NEWXMIPP:
         normalize_NewXmipp(img,bg_mask);
         break;
      case NEWXMIPP2:
         normalize_NewXmipp2(img,bg_mask);
         break;
      case RAMP:
       	 normalize_ramp(img,bg_mask);
	 break;
      case MICHAEL:
         normalize_Michael(img,bg_mask);
         break;
      case RANDOM:
         a=rnd_unif(a0,aF);
	 b=rnd_unif(b0,bF);
	 FOR_ALL_ELEMENTS_IN_MATRIX2D((*img)())
	    (*img)(i,j)=a*(*img)(i,j)+b;
         break;
   }
}
