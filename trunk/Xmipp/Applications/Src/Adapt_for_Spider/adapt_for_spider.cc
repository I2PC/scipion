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
#include <XmippInterface/xmippSpider.hh>

void Usage();

int main(int argc, char **argv) {
   string            command;
   FileName          fn_sel;
   FileName          fn_out;
   FileName          out_root;
   FileName          out_ext;
   string            ang1="rot",ang2="tilt",ang3="psi";
   int               maxcount;

   SelFile           SF, SF_out;
   DocFile           DF;
   matrix1D<float>   aux;
   int               selline;
   bool              newsel_style;

// Check command line ------------------------------------------------------
try {
   if (argc<2) REPORT_ERROR(1,"Adapt for Spider: Not enough parameters");
   command=argv[1];

   fn_out=get_param(argc, argv, "-o");
   if        (command=="rename") {
      fn_sel=get_param(argc, argv, "-i");
      out_root=get_param(argc, argv, "-oroot");
      out_ext=get_param(argc, argv, "-oext","");
   } else if (command=="translate_sel") {
      fn_sel=get_param(argc, argv, "-i");
      newsel_style=check_param(argc,argv,"-new_style");
   } else if (command=="generate_count") {
      maxcount=AtoI(get_param(argc, argv, "-max"));
   } else if (command=="extract_angles") {
      int i;
      fn_sel=get_param(argc, argv, "-i");
      if ((i=position_param(argc,argv,"-order"))!=-1) {
         if (i+3>=argc)
            REPORT_ERROR(1,"Adapt for Spider: Not enough parameters behind -ang\n");
         ang1=argv[i+1];
         ang2=argv[i+2];
         ang3=argv[i+3];
      }
   }
} catch (Xmipp_error XE) {cerr << XE; Usage(); exit(1);}

try {

// Perform operation
   if      (command=="rename") {
      SF.read(fn_sel);
      rename_for_Spider(SF,SF_out,out_root,out_ext);
      SF_out.write(fn_out);
   }
   else if (command=="generate_count") {
      generate_Spider_count(maxcount,DF);
      DF.write(fn_out);
   } else if (command=="translate_sel") {
      SF.read(fn_sel);
      translate_to_Spider_sel(SF,DF,newsel_style);
      DF.write(fn_out);
   } else if (command=="extract_angles") {
      SF.read(fn_sel);
      extract_angles(SF,DF,ang1,ang2,ang3);
      DF.write(fn_out);
   }
} catch (Xmipp_error XE) {cerr << XE;}
   exit(0);
}

// Usage -------------------------------------------------------------------
void Usage() {
   cerr << "Usage: adapt_for_Spider rename            : Generate correlative names\n"
        << "            -i <sel_file>                 : Input Xmipp selection file\n"
        << "            -o <output Xmipp SelFile>     : Output Xmipp SelFile\n"
        << "            -oroot <root_name>            : root name for output images\n"
        << "            [-oext <output extension="">  : if nothing is provided the same\n"
        << "                                            as the original images' one is used\n" 
        << "       adapt_for_Spider generate_count    : Generate corresponding Spider SelFile\n"
        << "            -max <max_count>              : Input Xmipp selection file\n"
        << "            -o <DocFile>                  : Output Spider CountFile\n"
        << "       adapt_for_Spider translate_sel     : Generate corresponding Spider SelFile\n"
        << "            -i <Xmipp Selfile>            : Input Xmipp selection file\n"
        << "            -o <Spider Selfile>           : Output Spider SelFile\n"
	<< "           [-new_style]                   : Generate new Spider Selfile style\n"
        << "       adapt_for_Spider extract_angles    : Generate a Docfile with angles\n"
        << "            -i <sel_file>                 : Input Xmipp selection file\n"
        << "            -o <Ang DocFile>              : Output Docfile\n"
        << "           [-order <ang1> <ang2> <ang3>   : order of the angles\n"
        << "                                            by default, psi, tilt, rot\n";
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Adapt_for_spider {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Adapt_for_Spider/Help/adapt_for_spider.html";
      help="Performs several operations for interfacing the SPIDER package";
      OPEN MENU menu_adapt_for_spider;
      COMMAND LINES {
	+ rename: xmipp_adapt_for_spider rename -i $SELFILE_IN -o $SELFILE_OUT
                      -oroot $OUTROOT [-oext $OEXT]
        + translate_sel: xmipp_adapt_for_spider translate_sel -i $SELFILE_IN
                      -o $SELFILE_OUT
        + generate_count: xmipp_adapt_for_spider generate_count -max $MAXCOUNT
                      -o $DOCFILE_OUT
        + extract_angles: xmipp_adapt_for_spider extract_angles -i $SELFILE_IN
                      -o $DOCFILE_OUT
                      [-ang $ANG1 $ANG2 $ANG3] [$ANGL1 $ANGL2 $ANGL3]
      } 
      PARAMETER DEFINITIONS {
        $SELFILE_IN {
	   label="Input SelFile";
	   type=file existing;
	}
        $SELFILE_OUT {
	   label="Output SelFile";
	   type=file;
	}
        $OUTROOT {
            label="Output filename root";
            help="Images are copied with a different filename root";
            type=text;
        }
        $OEXT {
           label="Different output extension";
           type=text;
        }
        $MAXCOUNT {
           label="Maximum count";
           help="Output DocFile goes from 1 to MAXCOUNT";
           type=NATURAL;
        }
        $DOCFILE_OUT {
           label="Document File Out";
           type=file;
        }
        #include "angles.mnu"
      }
   }

   MENU menu_adapt_for_spider {
      "I/O variables"
      $SELFILE_IN
      $SELFILE_OUT
      $OUTROOT
      OPT(-oext)
      $MAXCOUNT
      $DOCFILE_OUT
      OPT($ANGL1)
   }
*/
