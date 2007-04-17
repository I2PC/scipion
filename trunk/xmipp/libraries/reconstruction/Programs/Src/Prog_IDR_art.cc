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

#include "../Prog_IDR_art.hh"

void Prog_IDR_ART_Parameters::read(int argc, char **argv) {
   fn_exp=get_param(argc,argv,"-exp");
   fn_vol=get_param(argc,argv,"-vol");
   fn_ctf=get_param(argc,argv,"-ctf");
   fn_root=get_param(argc,argv,"-o");
   mu=AtoF(get_param(argc,argv,"-mu","1.8"));
   adjust_gray_levels=check_param(argc,argv,"-adjust_gray_levels");
}

void Prog_IDR_ART_Parameters::produce_side_info() {
   SF_exp.read(fn_exp);
   SF_ctf.read(fn_ctf);
   if (SF_ctf.ImgNo()!=SF_exp.ImgNo())
      REPORT_ERROR(1,"Prog_IDR_ART_Parameters: The number of images in "
         "the ctf and original selfiles do not match");
   V.read(fn_vol);
   V().set_Xmipp_origin();
}

void Prog_IDR_ART_Parameters::show() {
   cout << "Experimental images: " << fn_exp << endl
        << "Input volume: " << fn_vol << endl
        << "CTF Selfile: " << fn_ctf << endl
        << "Relaxation factor: " << mu << endl
        << "Output rootname: " << fn_root << endl
	<< "Adjust gray levels: " << adjust_gray_levels << endl
   ;
	
}

void Prog_IDR_ART_Parameters::Usage() {
   cerr << "Usage: IDR\n"
        << "   -exp <selfile>       : Selfile with experimental images\n"
        << "   -vol <volume>        : Voxel volume with the current reconstruction\n"
        << "   -ctf <selfile>       : Selfile with the CTFs (no phase flipping)\n"
        << "   -o <rootname>        : Output rootname\n"
        << "  [-mu <mu=1.8>]        : Relaxation factor\n"
	<< "  [-adjust_gray_levels] : Adjust gray levels\n"
   ;
}

/* IDR correction ---------------------------------------------------------- */
//#define DEBUG
void Prog_IDR_ART_Parameters::IDR_correction() {
   Projection Ireal, Inorm, Itheo, Itheo_CTF;

   SF_exp.go_first_ACTIVE();
   SF_ctf.go_first_ACTIVE();
   SF_new.clear();
   cerr << "Modifying input data ...\n";
   init_progress_bar(SF_exp.ImgNo());
   int istep=CEIL((double)SF_exp.ImgNo()/60.0);
   int imgs=0;
   while (!SF_exp.eof()) {
      // Read current input image
      FileName fn_img=SF_exp.NextImg();
      Ireal.read(fn_img);

      // Project the volume in the same direction
      project_Volume(V(), Itheo, YSIZE(Ireal()), XSIZE(Ireal()),
	 Ireal.rot(), Ireal.tilt(), Ireal.psi());

      // Copy to theo_CTF and resize
      Itheo_CTF()=Itheo();
      int Ydim=YSIZE(Itheo());
      int Xdim=XSIZE(Itheo());
      Itheo_CTF().set_Xmipp_origin();
      Itheo_CTF().window(FIRST_XMIPP_INDEX(2*Ydim),FIRST_XMIPP_INDEX(2*Xdim),
                         LAST_XMIPP_INDEX (2*Ydim),LAST_XMIPP_INDEX (2*Xdim));

      // Read CTF file
      ctf.FilterBand=CTF;
      ctf.ctf.read(SF_ctf.NextImg());
      ctf.ctf.enable_CTFnoise=false;
      ctf.ctf.Produce_Side_Info();
      ctf.generate_mask(Itheo_CTF());
      ctf.correct_phase();
      
      // Apply CTF
      ctf.apply_mask_Space(Itheo_CTF());
      Itheo_CTF().window(FIRST_XMIPP_INDEX(Ydim),FIRST_XMIPP_INDEX(Xdim),
                         LAST_XMIPP_INDEX (Ydim),LAST_XMIPP_INDEX (Xdim));

      // Center the all images
      Ireal().set_Xmipp_origin();
      Itheo().set_Xmipp_origin();

      // If adjust gray levels
      if (adjust_gray_levels) {
	 double avg, stddev, min_val, max_val;
	 Ireal().compute_stats(avg, stddev, min_val, max_val);
	 Itheo().statistics_adjust(avg, stddev);
      }

      // Apply IDR process
      FOR_ALL_ELEMENTS_IN_MATRIX2D(Ireal())
	 IMGPIXEL(Itheo,i,j)=mu*IMGPIXEL(Ireal,i,j)+
	    (IMGPIXEL(Itheo,i,j)-mu*IMGPIXEL(Itheo_CTF,i,j));

      // Save output image
      FileName fn_out=fn_root+ItoA(fn_img.get_number(),5)+".xmp";
      Itheo.write(fn_out);
      SF_new.insert(fn_out);

      //#define DEBUG
      #ifdef DEBUG
	 Itheo.write(fn_out.add_prefix("PPPtheo_")+".xmp");
	 Itheo_CTF.write(fn_out.add_prefix("PPPtheo_CTF_")+".xmp");
	 Ireal.write(fn_out.add_prefix("PPPreal_")+".xmp");
	 ImageXmipp save;
	 save()=Itheo()-mu*Itheo_CTF();
	 save.write(fn_out.add_prefix("PPPdiff_")+".xmp");
	 cout << "Press any key to continue\n"; char c; cin >> c;
      #endif

      if (imgs++%istep==0) progress_bar(imgs);
   }
   progress_bar(SF_exp.ImgNo());
   SF_new.write(fn_root+".sel");
}
