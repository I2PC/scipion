/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
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

#include "../xmippSpider.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippImages.hh>

// Generate Count File -----------------------------------------------------
void generate_Spider_count(int imax, DocFile &DF_out) {
   matrix1D<double>   aux(1);

   DF_out.clear();
   DF_out.append_comment((string)"Count for Spider up to "+ItoA(imax));

   for (aux(0)=1; aux(0)<=imax; aux(0)++)
      DF_out.append_data_line(aux);
}

// Translate to Spider selfile ---------------------------------------------
void translate_to_Spider_sel(SelFile &SF_in, DocFile &DF_out) {
   matrix1D<double>   aux(1);
   int               selline=1;

   DF_out.clear();
   DF_out.append_comment((string)"Translation for Spider of "+SF_in.name());

   SF_in.go_beginning();
   while (!SF_in.eof()) {
      if (!SF_in.Is_COMMENT()) {
         if (SF_in.Is_ACTIVE()) aux(0)=1;
         else                   aux(0)=0;
         DF_out.append_data_line(aux);
      }
      SF_in.next();
   }
}

// Extract angles ----------------------------------------------------------
void extract_angles(SelFile &SF_in, DocFile &DF_out,
   const string &ang1="rot", const string &ang2="tilt",
   const string &ang3="psi") {

   check_angle_descr(ang1);
   check_angle_descr(ang2);
   check_angle_descr(ang3);

   DF_out.clear();
   DF_out.append_comment((string)"Angles for "+SF_in.name()+
      ".   Angle order: "+ang1+" "+ang2+" "+ang3);

   int i=0;
   time_config();
   cerr << "Extracting angles ...\n";
   init_progress_bar(SF_in.ImgNo());
   while (!SF_in.eof()) {
      // Read image
      ImageXmipp P; P.read(SF_in.NextImg());
      if(P.Is_flag_set()==0 || P.Is_flag_set()>2)
      DF_out.append_angles(P.rot(),P.tilt(),P.psi(),
                            ang1,ang2,ang3);
      else if(P.Is_flag_set()==1)
      DF_out.append_angles(P.rot(), P.tilt(), P.psi(),
                           P.rot1(),P.tilt1(),P.psi1(),
                           ang1,ang2,ang3);
      else if(P.Is_flag_set()==2)
      DF_out.append_angles(P.rot(), P.tilt(), P.psi(),
                           P.rot1(),P.tilt1(),P.psi1(),
                           P.rot2(),P.tilt2(),P.psi2(),
                           ang1,ang2,ang3);
      i++;
      if (i%10==0) progress_bar(i);
   }
   progress_bar(SF_in.ImgNo());
}
#ifdef NEVERDEFINED
// write_angles
void write_angles(SelFile &SF_in, DocFile &DF_in,
   const string &ang1="rot", const string &ang2="tilt",
   const string &ang3="psi") {
   double rot, tilt, psi;
   int FirstLine_ColNo;
   
   check_angle_descr(ang1);
   check_angle_descr(ang2);
   check_angle_descr(ang3);
   
//   cout << "FirstLine_ColNo" << DF_in.FirstLine_ColNo();

   int i=0;
   time_config();
   cerr << "Writting new headers ...\n";
   init_progress_bar(SF_in.ImgNo());
   FirstLine_ColNo=DF_in.FirstLine_ColNo();
   while (!SF_in.eof()) {
      // Read image
      ImageXmipp P; P.read(SF_in.NextImg());
      P.clear_fFlag_flag();
      if(FirstLine_ColNo >=3){
         DF_in.get_angles(i+1, rot, tilt, psi, ang1, ang2, ang3);
         P.set_eulerAngles(rot,tilt,psi);
      }   
      if(FirstLine_ColNo >=6){
         DF_in.get_angles1(i+1, rot, tilt, psi, ang1, ang2, ang3);
         P.set_eulerAngles1(rot,tilt,psi);
      }   
      if(FirstLine_ColNo >=9){
         DF_in.get_angles2(i+1, rot, tilt, psi, ang1, ang2, ang3);
         P.set_eulerAngles2(rot,tilt,psi);
      }   
      
      P.write(P.name());
      i++;
      if (i%10==0) progress_bar(i);
   }
   progress_bar(SF_in.ImgNo());

}
#endif
// Rename for Spider -------------------------------------------------------
void rename_for_Spider(SelFile &SF_in, SelFile &SF_out, const FileName &fn_root,
   const FileName &out_ext) {
   FileName fn_in, fn_out;
   int counter=1;

   SF_out.clear();
   while (!SF_in.eof()) {
      fn_in=SF_in.NextImg();
      fn_out=fn_root+ItoA(counter,5);
      if (out_ext=="") fn_out=fn_out.add_extension(fn_in.get_extension());
      else             fn_out=fn_out.add_extension(out_ext);
      SF_out.insert(fn_out);

      cout << "Renaming " << fn_in << " as " << fn_out << endl;
      string command=(string)"cp "+fn_in+" "+fn_out;
      system(command.c_str());
      
      counter++;
   }
}

// Create empty Spider file ------------------------------------------------
void create_empty_Spider_file(const FileName &fn, int Zdim, int Ydim,
   int Xdim, bool reversed, size_t block_size) _THROW {
   unsigned char * buffer = (unsigned char*) calloc(sizeof(unsigned char),
      block_size);
   if (buffer==NULL) 
      REPORT_ERROR(1,"create_empty_Spider_file: No memory left");
   FILE * fd = fopen(fn.c_str(), "w");
   if (fd==NULL)
      REPORT_ERROR(1,(string)"create_empty_Spider_file: Cannot open file" + fn);

   // Write Header
   headerXmipp header;
   if (Zdim==1) header.headerType() = headerXmipp::IMG_XMIPP;
   else         header.headerType() = headerXmipp::VOL_XMIPP;
   header.set_dimension(Ydim, Xdim);
   header.Slices() = 1;
   header.set_header();
   header.set_time();
   header.set_date();
   header.write(fd, reversed);

   // Write file bulk
   int size=Zdim*Ydim*Xdim*sizeof(float);
   for (size_t i=0; i< size/block_size; i++)
       fwrite( buffer, sizeof(unsigned char), block_size, fd);
   fwrite(buffer, sizeof(unsigned char), size%block_size, fd);
   fclose(fd);
}

// 3D Radon ----------------------------------------------------------------
void radon_transform(VolumeXmipp &V_in, const FileName &fn_out,
   double Delta_rot, double Delta_tilt, int output_size) _THROW {
   if (output_size==-1) output_size=CEIL(1.5*XSIZE(V_in()));

   if (V_in.name()=="") V_in.write("superfeo.vol");
   else system(((string)"ln -s "+V_in.name()+" superfeo.vol").c_str());

   // Generate spider batch
   ofstream spider_batch;
   spider_batch.open("b01.vol");
   if (!spider_batch)
      REPORT_ERROR(1,"3D_Radon_transform:: Cannot open file for Spider batch");
   spider_batch
      << "rm 3d\n"
      << "superfeo\n"
      << "n\n"
      << "n\n"
      << output_size << endl
      << "superfeo2\n"
      << Delta_rot << " " << Delta_tilt << endl
      << FLOOR(0.5*XSIZE(V_in())) << endl
      << "n\n"
      << "en\n"
   ;
   spider_batch.close();
   
   char *spider_prog=getenv("SPIDER");
   if (spider_prog==NULL)
      REPORT_ERROR(1,"Project:: The environment variable SPIDER is not set");
   system(((string)spider_prog+" vol b01").c_str());
   system("rm LOG.vol results.vol* b01*.vol superfeo.vol");
   system(((string)"mv superfeo2.vol "+fn_out).c_str());
}

// 2D Radon ----------------------------------------------------------------
void radon_transform(ImageXmipp &I_in, const FileName &fn_out,
   double Delta_ang, int output_size) _THROW {
   if (output_size==-1) output_size=CEIL(1.5*XSIZE(I_in()));

   if (I_in.name()=="") I_in.write("superfeo.xmp");
   else system(((string)"ln -s "+I_in.name()+" superfeo.xmp").c_str());

   // Generate spider batch
   ofstream spider_batch;
   spider_batch.open("b01.xmp");
   if (!spider_batch)
      REPORT_ERROR(1,"2D_Radon_transform:: Cannot open file for Spider batch");
   spider_batch
      << "rm 2d\n"
      << "superfeo\n"
      << "-90 88\n"
      << Delta_ang << endl
      << "superfeo2\n"
      << output_size << endl
      << FLOOR(0.5*XSIZE(I_in())) << endl
      << "0 0\n"
      << "n\n"
      << "en\n"
   ;
   spider_batch.close();

   char *spider_prog=getenv("SPIDER");
   if (spider_prog==NULL)
      REPORT_ERROR(1,"Project:: The environment variable SPIDER is not set");
   system(((string)spider_prog+" xmp b01").c_str());
   system("rm LOG.xmp results.xmp* b01*.xmp superfeo.xmp");
   system(((string)"mv superfeo2.xmp "+fn_out).c_str());
}

// Fourier Radon transform -------------------------------------------------
void Fourier_transform_of_Radon_transform(const FileName &fn_in,
   const FileName &fn_out, double cutoff_freq,
   double Fermi_temperature) _THROW {
   system(((string)"ln -s "+fn_in+" superfeo.fft").c_str());

   // Generate spider batch
   ofstream spider_batch;
   spider_batch.open("b01.fft");
   if (!spider_batch)
      REPORT_ERROR(1,"Fourier_Radon_transform:: Cannot open file for Spider batch");
   spider_batch
      << "rm ft\n"
      << "superfeo\n"
      << "n\n"
      << "*\n"
      << "superfeo2\n"
      << "y\n"
      << "2\n"
      << "m\n"
      << "8\n"
      << "5\n"
      << cutoff_freq << endl
      << Fermi_temperature << endl
      << "a\n"
      << "n\n"
      << "n\n"
      << "en\n"
   ;
   spider_batch.close();

   char *spider_prog=getenv("SPIDER");
   if (spider_prog==NULL)
      REPORT_ERROR(1,"Project:: The environment variable SPIDER is not set");
   system(((string)spider_prog+" fft b01").c_str());
   system("rm LOG.fft results.fft* b01*.fft superfeo.fft");
   system(((string)"mv superfeo2.fft "+fn_out).c_str());
}

// Angular refinement ------------------------------------------------------
void Angular_refinement(const FileName &fn_vol, const FileName &fn_sel,
   const FileName &fn_report,
   double rot0, double rotF, double rot_step,
   double tilt0, double tiltF, double tilt_step,
   double psi0, double psiF, double psi_step,
   double max_shift) _THROW {

   SelFile SF, SF_kk;
   SF.read(fn_sel);
   FileName fn_ext=fn_vol.get_extension();

   // Rename input images
   rename_for_Spider(SF, SF_kk, "kk",fn_ext);
   SF_kk.go_first_ACTIVE();

   FileName fn_first=SF_kk.get_current_file();
   int first_image=fn_first.get_number();
   int last_image=first_image+SF_kk.ImgNo()-1;

   // Generate spider batch
   ofstream spider_batch;
   spider_batch.open(((string)"b01."+fn_ext).c_str());
   if (!spider_batch)
      REPORT_ERROR(1,"Angular refinement:: Cannot open file for Spider batch");
   #ifdef NEVER_DEFINED
      /* This is for spider V03.7 */
   spider_batch
      << "rm orfsf\n"
      << "0 0\n"
      << "*\n"
      << "*\n"
      << fn_vol.without_extension() << endl
      << fn_first.without_extension() << endl
      << first_image << "-" << last_image << endl
      << "0 0\n"
      << "0\n"
      << "peak00001\n"
      << "s\n"
      << max_shift << endl
      << "s\n"
      << "1\n"
      << tilt0 << " " << tiltF << endl
      << tilt_step << endl
      << rot0 << " " << rotF << endl
      << rot_step << endl
      << psi0 << " " << psiF << endl
      << psi_step << endl
      << "n\n"
      << fn_report << endl
      << "en\n"
   ;
   #endif
   spider_batch
      << "rm orfsfs\n"
      << "0 0\n"
      << "E\n"
      << "0\n"
      << "*\n"
      << "*\n"
      << fn_vol.without_extension() << endl
      << fn_first.without_extension() << endl
      << first_image << "-" << last_image << endl
      << "0 0\n"
      << "0\n"
      << "peak00001\n"
      << "s\n"
      << max_shift << endl
      << "s\n"
      << "1\n"
      << tilt0 << " " << tiltF << endl
      << tilt_step << endl
      << rot0 << " " << rotF << endl
      << rot_step << endl
      << psi0 << " " << psiF << endl
      << psi_step << endl
      << "n\n"
      << fn_report << endl
      << "en\n"
   ;
   spider_batch.close();

   char *spider_prog=getenv("SPIDER");
   if (spider_prog==NULL)
      REPORT_ERROR(1,"Project:: The environment variable SPIDER is not set");
   system(((string)"rm "+fn_report+"."+fn_ext).c_str());
   system(((string)spider_prog+" "+fn_ext+" b01").c_str());
   system(((string)"rm LOG."+fn_ext+" results."+fn_ext+
    "* b01*."+fn_ext).c_str());
   system(((string)"for i in peak?????."+fn_ext+" ; do rm $i ; done").
      c_str());

   SF_kk.go_first_ACTIVE();
   while (!SF_kk.eof())
      system(((string)"rm "+SF_kk.NextImg()).c_str());
}
