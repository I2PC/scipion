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

#define OTHER_FILTERS
#ifdef OTHER_FILTERS
   #include "../Prog_FourierFilter.hh"
#endif

#ifdef _HAVE_VTK
void Prog_IDR_ART_Parameters::read(const FileName &fn) {
   if (art_prm!=NULL) delete art_prm;
   art_prm=new Basic_ART_Parameters;
   art_prm->read(fn);

   FILE *fh;
   if ((fh = fopen(fn.c_str(), "r")) == NULL)
      REPORT_ERROR(3005,
         (string)"Prog_IDR_ART_Parameters::read: There is a problem "
         "opening the file "+fn);

   idr_iterations=AtoI(get_param(fh,"idr_iterations",0,"1"));
   it=AtoI(get_param(fh,"first iteration",0,"0"));
   dont_rewrite=check_param(fh,"dont_rewrite");
   fn_blob_volume=get_param(fh,"only_reproject",0,"");
   if (fn_blob_volume!="") idr_iterations=0;
   fn_ctf=get_param(fh,"CTF");
   mu0_list=get_vector_param(fh,"mu",-1);
   if (XSIZE(mu0_list)==0)
      {mu0_list.resize(1); mu0_list.init_constant(1);}
   muF_list=get_vector_param(fh,"muF",-1);
   fclose(fh);
}

void Prog_IDR_ART_Parameters::produce_side_info() _THROW {
   // Generate relaxation parameters
   if (XSIZE(muF_list)==0) mu_list=mu0_list;
   else {
      mu_list.resize(mu0_list);
      for (int i=0; i<XSIZE(mu0_list); i++)
         mu_list(i)=rnd_unif(mu0_list(i),muF_list(i));
   }

   // SelFiles and root name
   SF_current.read(art_prm->fn_sel);
   SF_original.read(art_prm->fn_sel);
   fn_root=art_prm->fn_root;

   // Read CTF
   if (Is_FourierImageXmipp(fn_ctf)) {
      ctf.FilterBand=ctf.FilterShape=FROM_FILE;
      ctf.fn_mask=fn_ctf;
      ctf.generate_mask(NULL);
      multiple_CTFs=FALSE;
   } else {
      SF_ctf.read(fn_ctf);
      if (SF_ctf.ImgNo()!=SF_original.ImgNo())
         REPORT_ERROR(1,"Prog_IDR_ART_Parameters: The number of images in "
            "the ctf and original selfiles do not match");
      multiple_CTFs=TRUE;
   }

   // Read Blob volume
   if (fn_blob_volume!="") vol_blobs.read(fn_blob_volume);
   art_prm->produce_Side_Info(vol_blobs,BASIC);
}

void Prog_IDR_ART_Parameters::show() {
   cout << "IDR parameters =================================================\n"
        << "IDR iterations= " << idr_iterations << endl
	<< "Don't rewrite=  "; print(cout,dont_rewrite);
   cout << "Only reproject= " << fn_blob_volume << endl
        << "CTF=            " << fn_ctf << endl
        << "Mu0=            " << mu0_list.transpose() << endl
        << "MuF=            " << muF_list.transpose() << endl
        << "Mu=             " << mu_list.transpose() << endl
        << "ART parameters =================================================\n"
        << /* art_prm << */ endl
   ;
}

void Prog_IDR_ART_Parameters::Usage() {
   cerr << "[idr_iterations=<no it=1>]            : No. of IDR iterations\n"
	<< "[mu=<mu=1>|\"[\"<mu0>,<mu1>,...\"]\"      : Relaxation parameters\n"
	<< "[muF=<muF>|\"[\"<muF0>,<muF1>,...\"]\"    : Final Relaxation parameters\n"
        << "CTF=<Xmipp Fourier file>|<selfile>    : CTF file\n"
        << "[dont_rewrite=<rw=No>]                : Keep or not the input projs\n"
	<< "[only_reproject=<Blob volume>]        : Do not reconstruct, use this volume\n"
	<< "                                        instead\n";
   cerr << "ART parameters =================================================\n";
   art_prm->usage();
}

/* IDR correction ---------------------------------------------------------- */
void Prog_IDR_ART_Parameters::IDR_correction(GridVolume &vol_blobs, int it) {
   FileName fn_img, fn_out;
   Projection Ireal, Inorm, Itheo, Itheo_CTF;
   vtkImageData *FFTtheo=NULL;

   SF_original.go_first_ACTIVE();
   SF_ctf.go_first_ACTIVE();
   SF_current.clear();
   cerr << "Modifying input data ...\n";
   init_progress_bar(SF_original.ImgNo());
   int istep=CEIL((double)SF_original.ImgNo()/60.0);
   int imgs=0;
   while (!SF_original.eof()) {
      // Read current input image
      fn_img=SF_original.NextImg();
      Ireal.read(fn_img);

      // Read CTF file if necessary
      if (multiple_CTFs) {
         ctf.FilterBand=ctf.FilterShape=FROM_FILE;
         ctf.fn_mask=SF_ctf.NextImg();
         ctf.generate_mask(NULL);
      }

      // Project the volume in the same direction
      project_Volume (vol_blobs, art_prm->blobprint,
	 art_prm->blobprint2, Itheo, Inorm,
	 YSIZE(Ireal()), XSIZE(Ireal()),
	 Ireal.rot(), Ireal.tilt(), Ireal.psi(),
	 FORWARD);

      // Apply CTF
      FFT_VTK(Itheo(),FFTtheo,TRUE);
      ctf.apply_mask(FFTtheo);
      IFFT_VTK(FFTtheo,Itheo_CTF(),TRUE);

      // Center the three images
      Ireal().set_Xmipp_origin();
      Itheo().set_Xmipp_origin();
      Itheo_CTF().set_Xmipp_origin();

      // Produce output filename
      int fn_number=fn_img.get_number();
      fn_out=fn_img.get_root()+"_idr";
      if (dont_rewrite) fn_out += ItoA(it,2);
      fn_out+=ItoA(fn_number,5);

      //#define DEBUG
      #ifdef DEBUG
	 Itheo.write((string)"theo_"+fn_out+".xmp");
	 Itheo_CTF.write((string)"theo_CTF_"+fn_out+".xmp");
	 Ireal.write((string)"real_"+fn_out+".xmp");
	 ImageXmipp save;
	 save()=Itheo()-mu(it)*Itheo_CTF();
	 save.write((string)"diff_"+fn_out+".xmp");
	 // cout << "Press any key to continue\n"; char c; cin >> c;
      #endif

      // Apply IDR process
      FOR_ALL_ELEMENTS_IN_MATRIX2D(Ireal()) {
	 IMGPIXEL(Itheo,i,j)=mu(it)*IMGPIXEL(Ireal,i,j)+
	    (IMGPIXEL(Itheo,i,j)-mu(it)*IMGPIXEL(Itheo_CTF,i,j));
      }

      // Save output image
      fn_out+=".xmp";
      Itheo.write(fn_out);
      SF_current.insert(fn_out);

      if (imgs++%istep==0) progress_bar(imgs);
   }
   progress_bar(SF_original.ImgNo());
   SF_current.write(fn_root+"_idr.sel");
}

/* Core routine ------------------------------------------------------------ */
void Basic_ROUT_IDR_Art(Prog_IDR_ART_Parameters &prm, VolumeXmipp &vol_recons) {
   // Reconstruction variables
   Plain_ART_Parameters plain_art_prm;
   prm.art_prm->tell |= TELL_SAVE_BLOBS;

   do {
      cout << "Running IDR iteration no. " << prm.it << " with mu= "
           << prm.mu(prm.it) << endl;
      // Run ART ...........................................................
      if (prm.fn_blob_volume=="") {
	 prm.art_prm->fn_sel =prm.SF_current.name();
	 prm.art_prm->fn_root=prm.fn_root+"_idr"+ItoA(prm.it,2);
	 if (prm.it>0) prm.art_prm->tell |= TELL_USE_INPUT_BLOBVOLUME;
	 Basic_ROUT_Art(*(prm.art_prm),plain_art_prm,vol_recons,prm.vol_blobs);
      }

      // Apply IDR correction ..............................................
      if (prm.it<prm.idr_iterations || prm.fn_blob_volume!="")
         prm.IDR_correction(prm.vol_blobs,prm.it);

      // Finish this iteration .............................................
      prm.it++;
   } while (prm.it<prm.idr_iterations+1);
}

#endif
