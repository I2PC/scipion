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
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/Programs/Prog_downsample.hh>
#include <Reconstruction/Programs/Prog_assign_CTF.hh>
#include <stdlib.h>
#include <unistd.h>

void Usage(const Prog_downsample_prm &down_prm);

int main(int argc, char **argv) {
   FileName fn_sel;
   int      input_type;
   #define  TIFF_IMAGES        0
   #define  RAW_IMAGES         1
   #define  DOWNSAMPLED_IMAGES 2
   bool     unzip;
   int      proc; // Number of available processors
   
   FileName output_dir;
   FileName fn_params;
   Prog_downsample_prm down_prm;
   
   // Get input parameters .................................................
   try {
      fn_sel=get_param(argc,argv,"-i");
      if      (check_param(argc,argv,"-tif"))   input_type=TIFF_IMAGES;
      else if (check_param(argc,argv,"-raw"))   input_type=RAW_IMAGES;
      else if (check_param(argc,argv,"-down"))  input_type=DOWNSAMPLED_IMAGES;
      else REPORT_ERROR(1,"No data type for micrographs");
      if (input_type==TIFF_IMAGES || input_type==RAW_IMAGES)
         down_prm.read(argc,argv,TRUE);
      unzip=check_param(argc,argv,"-unzip");
      output_dir=get_param(argc,argv,"-odir");
      fn_params=get_param(argc,argv,"-prm");
      proc=AtoI(get_param(argc,argv,"-proc","1"));
   } catch (Xmipp_error XE) {cout << XE; Usage(down_prm); exit(-1);}
   
   try {
      SelFile SF(fn_sel);
      char current_dir[512]="";
      char *check_current_dir=getcwd(current_dir,512);
      if (check_current_dir==NULL) strcpy(current_dir,".");

      // Read input assign_CTF parameters ..................................
      Prog_assign_CTF_prm assign_CTF_prm;
      assign_CTF_prm.read(fn_params,TRUE);
      
      // Prepare all batches ...............................................
      int N=0; // number of batches
      while (!SF.eof()) {
         FileName fn_img=SF.NextImg();
         FileName fn_root=fn_img.remove_all_extensions();
         
         ofstream fh_batch;
         fh_batch.open((fn_root+".bat").c_str(),ios::out);
         if (!fh_batch)
            REPORT_ERROR(1,(string)"Assign_multiple_CTFs: Cannot open "+
               fn_root+".bat for writing");
         fh_batch << "cd " << output_dir << endl << endl;
         
         FileName remove_file, current_file;
         int current_input_type=input_type;
         // Unzip ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
         if (unzip) {
            current_file=remove_file=fn_img.without_extension();
            fh_batch << "# Unzip input image" << endl
                     << "gunzip -c " << current_dir << "/" << fn_img
                     << " > " << current_file << endl;
            if (current_file.get_extension()!="tif") {
               fh_batch << "cp " << current_dir << "/" << current_file
                        << ".inf .\n";
               remove_file+=(string)" "+current_file+".inf";
            }
            fh_batch << endl;
         } else {
            current_file=(string)current_dir+"/"+fn_img;
         }
         
         // Tiff 2 raw ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
         if (current_input_type==TIFF_IMAGES) {
            fh_batch << "# Generate Raw image" << endl
                     << "xmipp_tiff2raw "<< current_file << " "
                     << fn_root+".raw" << endl;
            if (remove_file!="") {
               fh_batch << "rm -f " << remove_file << endl;
               remove_file="";
            }
            fh_batch << endl;
            current_file=fn_root+".raw";
            remove_file=current_file+" "+current_file+".inf";
            current_input_type=RAW_IMAGES;
         }

         // Downsample raw image ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
         if (current_input_type==RAW_IMAGES) {
            down_prm.fn_micrograph=current_file;
            down_prm.fn_downsampled=fn_root+".down";
            fh_batch << "#Downsample input image\n"
                     << "xmipp_downsample " << down_prm.command_line() << endl;
            if (remove_file!="") {
               fh_batch << "rm -f " << remove_file << endl;
               remove_file="";
            }
            fh_batch << endl;
            current_file=down_prm.fn_downsampled;
            remove_file=current_file+" "+current_file+".inf";
            current_input_type=DOWNSAMPLED_IMAGES;
         }

         // Assign CTF ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
         fh_batch << "# Assign different CTFs\n"
                  << "xmipp_assign_CTF -i " << fn_root+"_assign_CTF.param\n";
         if (remove_file!="") {
            fh_batch << "rm -f " << remove_file << endl;
            remove_file="";
         }
                  
         fh_batch.close();
         system(((string)"chmod 755 "+fn_root+".bat").c_str());
         N++;

         // Prepare assign CTF parameter file ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
         assign_CTF_prm.selfile_fn=(string)current_dir+"/"+fn_root+".sel";
         assign_CTF_prm.picked_fn=(string)current_dir+"/"+fn_root+".pos";
         assign_CTF_prm.ARMAfile=fn_root+"_arma";
         assign_CTF_prm.CTFfile=fn_root+"_ctf";
         assign_CTF_prm.image_fn=current_file;
         assign_CTF_prm.write(output_dir+"/"+fn_root+"_assign_CTF.param");
      }

      // Launch all processes ..............................................
      SF.go_first_ACTIVE();
      int n=0;
      for (int i=0; i<proc; i++) {
         string batches="(";
         for (int j=0; j<FLOOR(N/proc); j++) {
            FileName fn_img=SF.NextImg();
            FileName fn_root=fn_img.remove_all_extensions();
            batches+=fn_root+".bat; ";
            n++;
         }
         if (i<(N%proc)) {
            FileName fn_img=SF.NextImg();
            FileName fn_root=fn_img.remove_all_extensions();
            batches+=fn_root+".bat; ";
            n++;
         }
         batches+=") &";
         cout << "Launching " << batches << endl;
         system(batches.c_str());
      }
   } catch (Xmipp_error XE) {cout << XE;}
}

/* ------------------------------------------------------------------------- */
void Usage(const Prog_downsample_prm &down_prm) {
   cerr << "Usage: compute_CTF [options]\n"
        << "   -i <selfile>           : Selfile with all input micrographs\n"
        << "   -tif|-raw|-down        : Type of input data\n"
        << "  [-unzip]                : gunzip the input image\n"
        << "   -odir <output_dir>     : Output directory\n"
        << "   -prm  <textfile>       : Parameters for Assign_CTF\n"
        << "  [-proc <n=1>]           : Number of processors available\n"
   ;
   down_prm.usage();
}
