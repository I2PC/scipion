/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.uam.es)
 *              Carlos Oscar Sánchez Sorzano
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

#include "../Prog_assign_CTF.hh"

#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippMicrograph.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippHeader.hh>

/* Read parameters ========================================================= */
void Prog_assign_CTF_prm::read(const FileName &fn_prm, bool do_not_read_files) {
   // Read parameters for adjust CTF from input file
   adjust_CTF_prm.read(fn_prm);
   
   // Read specific parameters for this program from input file
   FILE *fh_param;
   if ((fh_param = fopen(fn_prm.c_str(), "r")) == NULL)
      REPORT_ERROR(1,(string)"assign_CTF: There is a problem "
            "opening the file "+fn_prm); 
                                      
   reversed          =check_param(fh_param,"reverse endian");
   N_horizontal      =AtoI(get_param(fh_param,"N_horizontal",0));
   N_vertical        =AtoI(get_param(fh_param,"N_vertical",0,"-1"));
   if (N_vertical==-1) N_vertical=N_horizontal;

   compute_at_particle  =check_param(fh_param,"compute_at_particle");
   micrograph_averaging =check_param(fh_param,"micrograph_averaging");
   piece_averaging      =check_param(fh_param,"piece_averaging");
   Nside_piece          =AtoI(get_param(fh_param,"Nside_piece",0,"5"));
   if (check_param(fh_param,"periodogram"))
        PSD_mode=Periodogram;
   else PSD_mode=ARMA;
   dont_adjust_CTF      =check_param(fh_param,"dont_adjust_CTF");

   if (!do_not_read_files) {
      image_fn          =get_param(fh_param,"image",0,"");
         selfile_mode=image_fn=="";
      if (selfile_mode) micrograph_averaging=true;
      if (selfile_mode) selfile_fn=get_param(fh_param,"selfile",0);
      else              selfile_fn=get_param(fh_param,"selfile",0,"");
      picked_fn         =get_param(fh_param,"picked",0,"");
      FileName fn_root  =image_fn.without_extension();
      if (PSD_mode==ARMA) PSDfn_root=fn_root+"_ARMA";
      else                PSDfn_root=fn_root+"_Periodogram";
   }
   fclose(fh_param);

   // Read ARMA parameters from input file
   if (PSD_mode==ARMA) ARMA_prm.read(fn_prm);
}

/* Write parameters ========================================================= */
void Prog_assign_CTF_prm::write(const FileName &fn_prm,
   string directory) {
   ofstream fh_param;
   fh_param.open(fn_prm.c_str(),ios::out);
   if (!fh_param)
      REPORT_ERROR(1,(string)"assign_CTF: There is a problem "
            "opening the file "+fn_prm+" for write");
   fh_param << "# Assign CTF parameters\n";
   string aux;
   bool remove_directories=directory!="";
   if (!remove_directories) aux=image_fn;
   else                     aux=directory+"/"+image_fn.remove_directories();
   if (!selfile_mode)
      fh_param << "image="                << aux                  << endl;
   fh_param << "N_horizontal="         << N_horizontal         << endl
            << "N_vertical="           << N_vertical           << endl;
   if (!remove_directories) aux=selfile_fn;
   else                     aux=directory+"/"+selfile_fn.remove_directories();
   fh_param << "selfile="              << aux                  << endl;
   if (!remove_directories) aux=picked_fn;
   else                     aux=directory+"/"+picked_fn.remove_directories();
   fh_param << "picked="               << aux                  << endl;
   if (compute_at_particle) fh_param << "compute_at_particle=yes\n";
   if (micrograph_averaging) fh_param << "micrograph_averaging=yes\n";
   if (piece_averaging) {
      fh_param << "piece_averaging=yes\n"
               << "Nside_piece=" << Nside_piece << endl;
   }
   if (PSD_mode==Periodogram) fh_param << "Periodogram=yes\n";
   if (dont_adjust_CTF) fh_param << "dont_adjust_CTF=yes\n";
  
   fh_param << endl;
   fh_param.close();

   adjust_CTF_prm.write(fn_prm,false);
   ARMA_prm.write(fn_prm,false);
}

/* Compute PSD by piece averaging ========================================== */
void Prog_assign_CTF_prm::PSD_piece_by_averaging(matrix2D<double> &piece,
   matrix2D<double> &psd) {
   int small_Ydim=2*YSIZE(piece)/Nside_piece;
   int small_Xdim=2*XSIZE(piece)/Nside_piece;
   matrix2D<double> small_piece(small_Ydim,small_Xdim);
   
   int Xstep=(XSIZE(piece)-small_Xdim)/(Nside_piece-1);
   int Ystep=(YSIZE(piece)-small_Ydim)/(Nside_piece-1);
   psd.init_zeros(small_piece);
   for (int ii=0; ii<Nside_piece; ii++)
      for (int jj=0; jj<Nside_piece; jj++) {
         // Take the corresponding small piece from the piece
         int i0=ii*Xstep;
         int j0=jj*Ystep;
         
         int i,j,ib,jb;
         for (i=0,ib=i0; i<small_Ydim; i++, ib++)
            for (j=0, jb=j0; j<small_Xdim; j++, jb++)
               DIRECT_MAT_ELEM(small_piece,i,j)=
                  DIRECT_MAT_ELEM(piece,ib,jb);

         // Compute the PSD of the small piece   
         matrix2D<double> small_psd;
         small_psd.init_zeros(small_piece);
         if (PSD_mode==ARMA) {
            // Compute the ARMA model
            matrix2D<double> ARParameters,MAParameters;
            double dSigma=CausalARMA(small_piece,ARMA_prm.N_AR,ARMA_prm.M_AR,
                ARMA_prm.N_MA,ARMA_prm.M_MA,ARParameters,MAParameters);
            ARMAFilter(small_piece,small_psd,ARParameters,MAParameters,dSigma);
         } else {
            // Compute the periodogram
            matrix2D< complex<double> > Periodogram;
            FourierTransform(small_piece,Periodogram);
            FFT_magnitude(Periodogram, small_psd);
            small_psd*=small_psd;
            small_psd*=small_Ydim*small_Xdim;
         }
         
         // Add to the average
         psd+=small_psd;
      }

      // Compute the average of all the small pieces and enlarge
      psd/=(Nside_piece*Nside_piece);
      CenterFFT(piece,true);
      psd.self_scale_to_size_Bspline(3,YSIZE(piece),XSIZE(piece));
      CenterFFT(piece,false);
      psd.threshold("below",0,0);
}

/* Main ==================================================================== */
//#define DEBUG
void Prog_assign_CTF_prm::process() {
   // Open input files -----------------------------------------------------
   // Open coordinates
   ifstream PosFile; // File with picked coordinates
   if (picked_fn!="") PosFile.open(picked_fn.c_str());

   // Open selfile with images
   SelFile SF; // Selfile
   if (selfile_fn!="") SF.read(selfile_fn);

   // Open the selfile for the CTFs, if there is a selfile of particles
   FileName fn_root;
   if (!selfile_mode) fn_root=image_fn.remove_all_extensions();
   else               fn_root=selfile_fn.remove_all_extensions();
   ofstream OutputFile_ctf;
   if (selfile_fn!="") OutputFile_ctf.open(
      (selfile_fn.without_extension()+".ctf.sel").c_str());

   // Open the micograph
   Micrograph M_in;
   int bits; // Micrograph depth
   int Ydim, Xdim; // Micrograph dimensions
   if (!selfile_mode) {
      M_in.open_micrograph(image_fn,reversed);
      bits=M_in.depth();
      M_in.size(Xdim, Ydim);
   }

   // Compute the number of divisions --------------------------------------
   int div_Number=0;
   int div_NumberX, div_NumberY;
   if (compute_at_particle) {
      // Check if sel file is empty
      if (SF.LineNo()==0)
         REPORT_ERROR(1,(string)"Prog_assign_CTF_prm: sel file "+SF.name()+
           "is empty ");

      // Count the number of lines in the Position file
      string line;
      PosFile.clear();       
      PosFile.seekg(0, ios::beg);
      while (getline(PosFile,line)) if (line[0]!='#') div_Number++;

      // check that the number of entries in the pos file is the right one
      if (SF.LineNo()!=div_Number) {
         cerr << "Prog_assign_CTF_prm: number of entries in "
              << "pos file: "<< picked_fn.c_str() 
              << "(" << div_Number << ") "
              << " and sel file " 
              << SF.name() << "(" << SF.LineNo() << ") "
              << "is different.\n"
              << "I cannot go any further sorry\n";
         exit(1);
      }
   } else if (micrograph_averaging) {
      if (!selfile_mode) {
         // If averaging, allow overlap among pieces
         div_NumberX=CEIL((double)Xdim/(N_horizontal/2))-1;
         div_NumberY=CEIL((double)Ydim/(N_vertical  /2))-1;
         div_Number=div_NumberX*div_NumberY;
      } else div_Number=SF.ImgNo();
    } else {
      // If not averaging, do not allow overlap
      div_NumberX=CEIL((double)Xdim/N_horizontal);
      div_NumberY=CEIL((double)Ydim/N_vertical  );
      div_Number=div_NumberX*div_NumberY;
    }

   // Process each piece ---------------------------------------------------
   PosFile.clear();       
   PosFile.seekg(0, ios::beg); // Start of file
   SF.go_beginning();
   ImageXmipp psd_avg;
   cerr << "Computing models of each piece ...\n";
   init_progress_bar(div_Number);

   // Prepare these filenames in case they are needed
   FileName fn_avg, fn_avg_model;
   if      (micrograph_averaging && PSD_mode==ARMA)        fn_avg=fn_root+"_ARMAavg.psd";
   else if (micrograph_averaging && PSD_mode==Periodogram) fn_avg=fn_root+"_Periodogramavg.psd";
   int N=1;      // Index of current piece
   int i=0, j=0; // top-left corner of the current piece
   while (N<=div_Number) {
      // Compute the top-left corner of the piece ..........................
      if (compute_at_particle) {
         // Read position of the particle
         string line;
         getline(PosFile,line);
         while (line[0]=='#') getline(PosFile,line);       
         float fi, fj; sscanf(line.c_str(),"%f %f",&fj,&fi);
         i = (int) fi;
         j = (int) fj;
         #ifdef DEBUG
            cout << "line" << line << endl;
            cout << "read from file (j,i)= (" << j << "," << i << ")" <<endl;
            cout << "Particle file name: " << SF.get_current_file() << endl;
         #endif

         // j,i are the window center, we need the top-left corner
         j -= (int) (N_horizontal/2);
         i -= (int) (N_vertical/2);
         if (i<0) i=0;
         if (j<0) j=0;
         if (i>Ydim-N_horizontal) i=Ydim-N_horizontal-1;
         if (j>Xdim-N_vertical)   j=Xdim-N_vertical-1;
      } else {
         if (!selfile_mode) {
            int Xstep=N_horizontal, Ystep=N_vertical;
            if (micrograph_averaging) {Xstep/=2; Ystep/=2;}
            i=((N-1)/div_NumberX)*Ystep;
            j=((N-1)%div_NumberX)*Xstep;
         }
      }
      
      // test if the full piece is inside the micrograph
      if (!selfile_mode) {
         if (i+N_vertical>Ydim)   i=Ydim-N_vertical;
         if (j+N_horizontal>Xdim) j=Xdim-N_horizontal;
      }
      
      // Extract micrograph piece ..........................................
      matrix2D<double> piece(N_vertical,N_horizontal);
      if (!selfile_mode) {
         for (int k=0; k<YSIZE(piece); k++)
             for (int l=0; l<XSIZE(piece); l++)
                 piece(k,l)=M_in(j+l,i+k);
      } else {
         ImageXmipp I; I.read(SF.get_current_file());
         piece=I();
      }
      piece.statistics_adjust(0,1);

      // Estimate the power spectrum .......................................
      ImageXmipp psd; psd().resize(piece);
      if (!piece_averaging)
         if (PSD_mode==ARMA) {
            // Compute the ARMA model
            matrix2D<double> ARParameters,MAParameters;
            double dSigma=CausalARMA(piece,ARMA_prm.N_AR,ARMA_prm.M_AR,
                ARMA_prm.N_MA,ARMA_prm.M_MA,ARParameters,MAParameters);
            ARMAFilter(piece,psd(),ARParameters,MAParameters,dSigma);
         } else {
            // Compute the periodogram
            matrix2D< complex<double> > Periodogram;
            FourierTransform(piece,Periodogram);
            FFT_magnitude(Periodogram, psd());
            psd()*=psd();
            psd()*=N_vertical*N_horizontal;
         }
      else PSD_piece_by_averaging(piece,psd());

      // Perform averaging if applicable ...................................
      if (micrograph_averaging) {
         if (N==1) psd_avg() =psd();
         else      psd_avg()+=psd();
         if (micrograph_averaging && PSD_mode==ARMA) {
            psd_avg.write(fn_avg);
            if (N==1) system(((string)"xmipp_show -psd "+
               fn_avg+" -poll &").c_str());
         }
      }

      // Compute the theoretical model if not averaging ....................
      if (!micrograph_averaging) {
         FileName piece_fn_root;
         if (compute_at_particle) {
            piece_fn_root = SF.get_current_file();
            piece_fn_root = piece_fn_root.get_baseName();
            SF.next();                       
         } else
            piece_fn_root=PSDfn_root+ItoA(N,5);

         psd.write(piece_fn_root+".psd");

      	 if (!dont_adjust_CTF) {
            // Estimate the CTF parameters of this piece 
            adjust_CTF_prm.fn_ctf=piece_fn_root+".psd";
            if (!dont_adjust_CTF) {
               double fitting_error=ROUT_Adjust_CTF(adjust_CTF_prm,false);
               if (compute_at_particle)
        	  OutputFile_ctf << piece_fn_root+".ctfmodel_halfplane" << " 1\n";
            }
	 }
      }

      // Increment the division counter
      progress_bar(++N);
      if (selfile_mode) SF.NextImg();
   }
   M_in.close_micrograph();
   progress_bar(div_Number);

   // If averaging, compute the CTF model ----------------------------------
   if (micrograph_averaging) {
      psd_avg()/=div_Number;
      psd_avg.write(fn_avg);

      if (!dont_adjust_CTF) {
	 // Estimate the CTF parameters
	 cerr << "Adjusting CTF model to the PSD ...\n";
	 adjust_CTF_prm.fn_ctf=fn_avg;
	 double fitting_error=ROUT_Adjust_CTF(adjust_CTF_prm,false);
      }
   }

   // Assign a CTF to each particle ----------------------------------------
   if (!compute_at_particle && selfile_fn!="" && !dont_adjust_CTF) {
      // Process the Selfile
      SF.go_beginning();
      if (!selfile_mode) {
         PosFile.close();
         PosFile.open(picked_fn.c_str());
         if (!PosFile)
            REPORT_ERROR(1,(string)"Prog_assign_CTF_prm::process: Could not open "+
               picked_fn+" for reading");
      }
      while (!SF.eof()) {  
         if (!selfile_mode) {
            string line;
            getline(PosFile,line);
            while (line[0]=='#') getline(PosFile,line);
            if (SF.Is_ACTIVE()) {
               if (!micrograph_averaging) {
                  // Read coordinates of the particle
                  float fX, fY; sscanf(line.c_str(),"%f %f",&fX,&fY);
                  int Y = (int) fY;
                  int X = (int) fX;

                  // Decide which is its piece
                  int idx_X=FLOOR((double)X/N_horizontal);
                  int idx_Y=FLOOR((double)Y/N_vertical);
                  int idx_piece=idx_Y*div_NumberX+idx_X+1;
                  OutputFile_ctf << PSDfn_root+ItoA(idx_piece,5)+".ctfmodel_halfplane"
                                 << " 1\n";
               } else
                  OutputFile_ctf << fn_avg.without_extension()+".ctfmodel_halfplane"
                                 << " 1\n";
            }
         } else {
            OutputFile_ctf << fn_avg.without_extension()+".ctfmodel_halfplane"
                           << " 1\n";
         }
         SF.next();
      }
   }
   if (selfile_fn!="") OutputFile_ctf.close();
   PosFile.close();
}
