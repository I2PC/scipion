/***************************************************************************
 *
 * Authors:    Sjors Scheres
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

/* INCLUDES ---------------------------------------------------------------- */
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippImages.hh>
#include <XmippData/xmippSelFiles.hh>
#include <XmippData/xmippDocFiles.hh>
/* PROTOTYPES -------------------------------------------------------------- */
void Usage();

#define PLUS  1.
#define MINUS -1.

/* MAIN -------------------------------------------------------------------- */
int main (int argc,char *argv[]) {
   float           rot, tilt, psi, xshift, yshift;
   int             i,ncol,col_rot=1,col_tilt=2,col_psi=3,col_xshift=4,col_yshift=5,col_weight=-1,col_mirror=-1;
   float           sign_rot=PLUS,sign_tilt=PLUS,sign_psi=PLUS;
   float           sign_xshift=PLUS,sign_yshift=PLUS;
   int             key_img;
   string          root,ext;
   bool            extract,assign,reset,verb,do_weights=false,do_mirrors=false;
   float           weight,mirror;
   FileName        fn_img,fn_out,fn_tst;
   SelFile         SF;
   DocFile         DF;
   ImageXmipp      img;

// Check command line options ===========================================
   try {

     extract=check_param(argc,argv,"-extract");
     assign=check_param(argc,argv,"-assign");
     reset=check_param(argc,argv,"-reset");
     if (!(extract || assign || reset)) REPORT_ERROR(1,"Please specify: -extract, -reset  or -assign");

     if (extract) {
       SF.read(get_param(argc,argv,"-i"));
       fn_out=get_param(argc,argv,"-o");
     }
     else if (reset) {
       SF.read(get_param(argc,argv,"-i"));
     }
     else if (assign) {
       DF.read(get_param(argc,argv,"-i"));
       fn_out=get_param(argc,argv,"-o","");
       verb=check_param(argc,argv,"-verb");

       // Columns numbers
       if ((i=position_param(argc,argv,"-columns"))!=-1) {
	 if (i+5>=argc) {
	   REPORT_ERROR(1,"Not enough integers after -columns");
	 }
	 col_rot=AtoI(argv[i+1]);
	 col_tilt=AtoI(argv[i+2]);
	 col_psi=AtoI(argv[i+3]);
	 col_xshift=AtoI(argv[i+4]);
	 col_yshift=AtoI(argv[i+5]);

	 // Check colum signs 
	 if (col_rot<0) sign_rot=MINUS;
	 if (col_tilt<0) sign_tilt=MINUS;
	 if (col_psi<0) sign_psi=MINUS;
	 if (col_xshift<0) sign_xshift=MINUS;
	 if (col_yshift<0) sign_yshift=MINUS;

       }

       // Also assign weights or mirror flags?
       ncol=5;
       do_weights=check_param(argc,argv,"-weight");
       if (do_weights) {
	 col_weight=AtoI(get_param(argc,argv,"-weight","6"));
	 ncol++;
       }
       do_mirrors=check_param(argc,argv,"-mirror");
       if (do_mirrors) {
	 col_mirror=AtoI(get_param(argc,argv,"-mirror","7"));
	 ncol++;
       }
     }

   } catch (Xmipp_error XE) {cout << XE; Usage();}

   try {

// Extracting information  ==================================================
     if (extract) {
       DF.reserve(SF.ImgNo());
       matrix1D<double> docline;
       DF.append_comment("Headerinfo columns: rot (1) , tilt (2), psi (3), Xoff (4), Yoff (5), Weight (6), Mirror (7)");

       docline.init_zeros(7);
       int j=0;
       SF.go_beginning();
       while (!SF.eof()) {
	 fn_img=SF.NextImg();
	 img.read(fn_img);
	 docline(0)=img.rot();
	 docline(1)=img.tilt();
	 docline(2)=img.psi();
	 docline(3)=img.Xoff();
	 docline(4)=img.Yoff();
	 docline(5)=img.weight();
	 docline(6)=img.flip();
	 DF.append_comment(fn_img);
	 DF.append_data_line(docline);
	 j++;
       }
       DF.write(fn_out);
     }
// Resetting headers  ==================================================
     else if (reset) {
       SF.go_beginning();
       cerr <<" Resetting all angles, origin offsets, weights and mirror flags to zero ... "<<endl;
       while (!SF.eof()) {
	 img.read(SF.NextImg());
	 img.set_eulerAngles(0.,0.,0.);
	 img.set_originOffsets(0.,0.);
	 img.weight()=0.;
	 img.flip()=0.;
	 img.write(img.name());
       }
     }
// Assigning information  ==================================================
     else if (assign) {

       // Assigning angles and shifts from document file
       DF.go_beginning();
       if (DF.get_current_line().Is_comment()) fn_tst=(DF.get_current_line()).get_text();
       if (strstr(fn_tst.c_str(),"Headerinfo")==NULL) {

	 // Non-NewXmipp type document file
	 cerr << "Warning!! Docfile is of non-NewXmipp type. "<<endl;
	 if (fn_out=="") REPORT_ERROR(1,"Please specify corresponding selfile: -o <selfile>");
	 SF.read(fn_out);
	 if (SF.ImgNo()!=DF.get_last_key()) REPORT_ERROR(1,"docfile and corresponding selfile have unequal (active) entries");
	 else cerr << "Corresponding selfile has expected number of entries"<<endl;
	 SF.go_first_ACTIVE();
	 DF.go_beginning();
	 while (!SF.eof()) {
	   fn_img=SF.get_current_file();
	   img.read(fn_img);
	   img.clear_fFlag_flag();
	   DF.adjust_to_data_line();
	   if (col_rot==0  )  rot=0.;     else rot    = DF(ABS(col_rot)-1)*sign_rot;
	   if (col_tilt==0 ) tilt=0.;     else tilt   = DF(ABS(col_tilt)-1)*sign_tilt;
	   if (col_psi==0  )  psi=0.;     else psi    = DF(ABS(col_psi)-1)*sign_psi;
	   if (col_xshift==0 ) xshift=0.; else xshift = DF(ABS(col_xshift)-1)*sign_xshift;
	   if (col_yshift==0 ) yshift=0.; else yshift = DF(ABS(col_yshift)-1)*sign_yshift;
	   if (do_weights) weight = DF(ABS(col_weight)-1);
	   if (do_mirrors) mirror=DF(ABS(col_mirror)-1);
	 
	   // Assign angles
	   img.set_eulerAngles(rot,tilt,psi);
	   img.set_originOffsets(xshift,yshift);
	   if (do_weights) img.weight()=weight;
	   if (do_mirrors) img.flip()=mirror;
	   if (verb) {
	     cout << fn_img  << " : rot = " << rot << " tilt = " << tilt
		  << " psi = " << psi << " Xoff = " << xshift <<" Yoff = " << yshift;
	     if (do_weights) cout <<" Weight = " <<weight;
	     if (do_mirrors) cout <<" Mirror = " <<DF(ABS(col_mirror)-1);
	     cout << endl;
	   }
	   img.write(fn_img);
	   
	   //Move to next line in SF and document file
	   SF.NextImg();
	   DF.next();
	 }
       } else {

	 // NewXmipp-type document file
	 int n=0;
	 int nmax=DF.dataLineNo();
	 while (n<nmax) {
	   n++;
	   DF.next();
	   if (DF.get_current_line().Is_comment()) fn_img=((DF.get_current_line()).get_text()).erase(0,3);
	   else  REPORT_ERROR(1,"Problem with NewXmipp-type document file");
	   img.read(fn_img);
	   img.clear_fFlag_flag();
	   DF.adjust_to_data_line();
	   if (col_rot==0  )  rot=0.;     else rot    = DF(ABS(col_rot)-1)*sign_rot;
	   if (col_tilt==0 ) tilt=0.;     else tilt   = DF(ABS(col_tilt)-1)*sign_tilt;
	   if (col_psi==0  )  psi=0.;     else psi    = DF(ABS(col_psi)-1)*sign_psi;
	   if (col_xshift==0 ) xshift=0.; else xshift = DF(ABS(col_xshift)-1)*sign_xshift;
	   if (col_yshift==0 ) yshift=0.; else yshift = DF(ABS(col_yshift)-1)*sign_yshift;
	   if (do_weights) weight = DF(ABS(col_weight)-1);
	   if (do_mirrors) mirror=DF(ABS(col_mirror)-1);
	   
	   // Assign angles
	   img.set_eulerAngles(rot,tilt,psi);
	   img.set_originOffsets(xshift,yshift);
	   if (do_weights) img.weight()=weight;
	   if (do_mirrors) img.flip()=mirror;
	   if (verb) {
	     cout << fn_img  << " : rot = " << rot << " tilt = " << tilt
		  << " psi = " << psi << " Xoff = " << xshift <<" Yoff = " << yshift;
	     if (do_weights) cout <<" Weight = " <<weight;
	     if (do_mirrors) cout <<" Mirror = " <<DF(ABS(col_mirror)-1);
	     cout << endl;
	   }
	   img.write(fn_img);
	 }
       }
       if (!verb) cerr <<" done"<<endl;
     }
     
   } catch (Xmipp_error XE) {cout << XE;}
}

/* Usage ------------------------------------------------------------------- */
void Usage() {
    printf("Purpose:\n");
    printf(" Access the geometric transformation (angles & shifts) in the header of 2D-images.\n");
    printf("Usage:\n");
    printf(" There are two modes: getting, resetting or putting information \n");
    printf("   headerinfo -extract \n");
    printf("        -i <selfile>       : input selfile\n");
    printf("        -o <docfile>       : output document file\n");
    printf("   headerinfo -reset \n");
    printf("        -i  <selfile>      : reset all information in the header to zero values \n");
    printf("   headerinfo -assign  \n");
    printf("        -i <docfile>       : input document file\n");
    printf("       [-o <selfile>]      : selfile containing the images to be assessed\n");
    printf("                             (only necessary for non-NewXmipp docfiles) \n");
    printf("       [-verb]             : output assigned information to screen \n");
    printf("       [-columns] <rot=1> <tilt=2> <psi=3> <Xoff=4> <Yoff=5> \n"
           "                           : where the 5 integers are the column numbers for the \n"
           "                           : respective angles and offsets in the docfile\n"
           "                           : Negative column numbers result in a sign change\n"
           "                           : Zeros result in zero values\n");
    printf("       [-weight <col_w=6>] : Set ML-weights (from column number col_w) \n");
    printf("       [-mirror <col_m=7>] : Set mirror-flag (from column col_m) (0=no-flip; 1=flip)\n");
    exit(1);
}
