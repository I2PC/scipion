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

/* Includes ================================================================ */
#include <XmippData/Programs/Prog_downsample.hh>
#include <XmippData/xmippArgs.hh>
#include "xvsmooth.h"

/* Prototypes ============================================================== */
void Usage(const Prog_downsample_prm &prm);

/* Main program ============================================================ */
int main(int argc, char **argv) {
   Prog_downsample_prm prm;
   bool                smooth;
   bool                reversed;

   // Get input parameters -------------------------------------------------
   try {
      prm.read(argc,argv);
      smooth         = check_param(argc,argv,"-smooth");
      reversed       = check_param(argc,argv,"-reverse_endian");
   } catch (Xmipp_error XE) {cout << XE; Usage(prm); exit(1);}

   try {
      prm.generate_kernel();
      prm.open_input_micrograph();
      prm.create_empty_output_file();
      if (smooth) {
         Micrograph Mp;
         Mp.open_micrograph(prm.fn_downsampled, reversed);
         byte rgb[256]; for (int i=0; i<256; i++) rgb[i]=i;
         byte *result = SmoothResize((byte *) (prm.M.array8()),
            prm.Xdim, prm.Ydim, prm.Xpdim, prm.Ypdim,
            rgb, rgb, rgb, rgb, rgb, rgb, 256); 
         for (int i=0; i<prm.Ypdim; i++)
             for (int j=0; j<prm.Xpdim; j++)
                 Mp.set_val(j,i,result[i*prm.Xpdim+j]);
         Mp.close_micrograph();
      } else prm.Downsample();
      prm.close_input_micrograph();
   } catch (Xmipp_error XE) {cout << XE;}
}

/* Usage =================================================================== */
void Usage(const Prog_downsample_prm &prm) {
   cerr << "Purpose: This file allows you to downsample raw images\n"
        << "Usage: downsample [parameters]\n"
        << "   -i <input_file>        : Raw input file, <input_file>.inf\n"
        << "                            must exist\n"
        << "   -o <output_file>       : Must be different from input one\n"
        << "  [-smooth]               : Use Smoothing for downsampling\n"
   ;
   prm.usage();
}

/* Colimate menu =========================================================== */
/*Colimate:
   PROGRAM DownSample {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/DownSample/Help/downsample.html";
      help="Downsample raw micrographs";
      OPEN MENU Downsample;
      COMMAND LINES {
         + usual: xmipp_downsample -i $FILEIN -o $FILEOUT
                  [-output_bits $OBITS $OBITS_LIST]
                   -Xstep $XSTEP [-Ystep $YSTEP]
                  [-smooth] 
                  [-kernel $KERNEL $KERNEL_LIST
                     [$YDIM $XDIM] [$R] [$SIGMA] [$POWER]]
      }
      PARAMETER DEFINITIONS {
         $FILEIN {
            label="Input micrograph";
            help="Raw file";
            type=file existing;
         }
         $FILEOUT {
            label="Ouput micrograph";
            help="Raw file; it cannot be the same as the input one.";
            type=file;
         }
         $OBITS {shown=no; type=natural;}
         $OBITS_LIST {
            label="Output bits";
            help="By default, 16";
            type=list {
               "16" {$OBITS="16";}
               "8"  {$OBITS="8";}
            };
         }
         $XSTEP {
            label="X step";
            help="Downsampling factor in the X axis; must be integer";
            type=natural;
         }
         $YSTEP {
            label="Y step";
            help="Downsampling factor in the X axis; must be integer";
            type=natural;
         }
         OPT(-smooth) {
            label="Apply a xv-like smoothing";
            help="No kernel must be applied in this case";
         }
         OPT($KERNEL) {label="Kernel definition";}
         $KERNEL {shown=no; type=text;}
         $KERNEL_LIST {
            label="Kernel type";
            help="Notice that the valid parameters are enabled";
            type=list {
               "rectangle" {$KERNEL="rectangle"; OPT($YDIM)=1;
                            OPT($R)=0; OPT($SIGMA)=0; OPT($POWER)=0;}
               "circle"    {$KERNEL="circle"; OPT($YDIM)=0;
                            OPT($R)=1; OPT($SIGMA)=0; OPT($POWER)=0;}
               "gaussian"  {$KERNEL="gaussian"; OPT($YDIM)=0;
                            OPT($R)=1; OPT($SIGMA)=1; OPT($POWER)=0;}
               "pick"      {$KERNEL="pick"; OPT($YDIM)=0;
                            OPT($R)=0; OPT($SIGMA)=0; OPT($POWER)=0;}
               "sinc"      {$KERNEL="sinc"; OPT($YDIM)=0;
                            OPT($R)=0; OPT($SIGMA)=0; OPT($POWER)=1;}
            };
         }
         OPT($YDIM) {label="Rectangle dimensions";}
         $YDIM  {type=natural; label="Kernel Y dimension";}
         $XDIM  {type=natural; label="Kernel X dimension";}
         $R     {type=natural; label="Kernel radius";}
         $SIGMA {type=natural; label="Kernel gaussian sigma";}
         $POWER {type=natural; label="Kernel power"; by default=0.95;}
      }
   }
   MENU Downsample {
      "I/O parameters"
      $FILEIN
      {$FILEOUT OPT($OBITS)}
      "Downsampling factors"
      $XSTEP
      OPT($YSTEP)
      "Kernel shape"
      OPT(-smooth)
      OPT($KERNEL)
   }
*/
