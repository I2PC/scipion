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

#include <XmippData/xmippDocFiles.hh>
#include <XmippData/xmippArgs.hh>

void Usage();

/* ------------------------------------------------------------------------- */
/* Main                                                                      */
/* ------------------------------------------------------------------------- */
int main (int argc, char *argv[]) {
   FileName fn_DF, fn_SF, fn_SF_out;
   int      col;
   double   limit0, limitF;
   bool     en_limit0, en_limitF;

   // Get input parameters .................................................
   try {
      fn_DF     = get_param(argc, argv, "-doc");
      fn_SF     = get_param(argc, argv, "-sel");
      fn_SF_out = get_param(argc, argv, "-sel_out","");
      col       = AtoI(get_param(argc, argv, "-col","0"));
      en_limit0 = check_param(argc,argv,"-limit0");
      if (en_limit0)
         limit0 = AtoF(get_param(argc,argv,"-limit0"));
      en_limitF = check_param(argc,argv,"-limitF");
      if (en_limitF)
         limitF = AtoF(get_param(argc,argv,"-limitF"));
      if (fn_SF_out=="") fn_SF_out=fn_SF;
   } catch (Xmipp_error XE) {cout << XE; Usage(); exit(0);}
   
   // Process ..............................................................
   try {
      SelFile SF; SF.read(fn_SF);
      DocFile DF; DF.read(fn_DF);
      if (SF.ImgNo(SelLine::ACTIVE)+SF.ImgNo(SelLine::DISCARDED)!=
         DF.dataLineNo())
         REPORT_ERROR(1,"Select images: SelFile and DocFile do not have the "
            "same number of lines");
      if (col>=DF.FirstLine_ColNo())
         REPORT_ERROR(1,"Select images: Column not valid for this DocFile");
      select_images(DF, SF, col, en_limit0, limit0, en_limitF, limitF);
      SF.write(fn_SF_out);
   } catch (Xmipp_error XE) {cout << XE;}   
}

/* Usage ------------------------------------------------------------------- */
void Usage() {
   cerr << "Usage: select_images\n"
        << "   -doc <docfile>                   : Input document file\n"
        << "   -sel <selfile>                   : Corresponding output file\n"
        << "  [-sel_out <selfile_out>=<selfile>]: Output selfile. By default,\n"
        << "                                      the input one\n"
        << "  [-col <col=0>]                    : Column of the docfile\n"
        << "  [-limit0 <limit0>]                : Values below this are discarded\n"
        << "  [-limitF <limitF>]                : Values above this are discarded\n"
   ;
}
