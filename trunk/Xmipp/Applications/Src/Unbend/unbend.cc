/***************************************************************************
 *
 * Authors:     Debora Gil
                Roberto Marabini
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
#include <Reconstruction/Programs/Prog_Umbend.hh>

void Usage(char *argv[]);

int main (int argc, char **argv) {
	
   ImUmbend prm;
   
     
   // Get command line parameters ------------------------------------------
   try {
      prm.inImfile  = get_param(argc, argv, "-i");
      prm.outImfile = get_param(argc, argv, "-o");
      prm.FN_Correlation = get_param(argc,argv,"-cor");
      prm.cc_peak_factor =  AtoF(get_param(argc,argv,"-cc_peak_factor","0.0"));
  
  } catch (Xmipp_error XE) {Usage(argv);}

   // Main program ---------------------------------------------------------
   try {
      ROUT_Umbend(prm);
   } catch (Xmipp_error XE) {cout<<XE;}
   exit(0);
}

void Usage (char *argv[]) {
    cout << "Purpose:\n"
         << "    Computes 2D crystal distorsion AND umbends the crystal\n"
         << "    The input is a cor (MRC) file\n"
    
         << "Usage:" << argv[0] <<" -i InImage -o OutImage -cor filename -cc_peak_factor cc_peak_factor" << endl << endl
         << "\t-i               :  Input Image" << endl
         << "\t-o               :  Output image" << endl
         << "\t-cor             : Correlation file" << endl
         << "\t-cc_peak_factor  : crosscorrelation thershold (0-1)" << endl
    ;    
    exit(1);

}
