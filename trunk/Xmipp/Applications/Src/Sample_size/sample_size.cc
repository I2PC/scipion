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

#include <XmippData/xmippTypes.hh>

void Usage();

// Remove +- from input string, the output is still there
void remove_plus_minus(string &str) {
   for (int i=1; i<=str.length(); i++)
      if (str[i]=='+' || str[i]=='-') str[i]=' ';
}

/* ------------------------------------------------------------------------- */
/* Main                                                                      */
/* ------------------------------------------------------------------------- */
int main (int argc, char **argv) {
   float avg, stddev;
   float alpha;
   float r;
   
// Read command line -------------------------------------------------------
   // Read arguments
   try {
      string str=get_param(argc,argv,"-stats");
      sscanf(str.c_str(),"%f+-%f",&avg,&stddev);
      
      alpha=AtoF(get_param(argc,argv,"-alpha","0.005"));
      r=AtoF(get_param(argc,argv,"-r","0.01"));
   } catch (Xmipp_error XE) {cerr << XE; Usage(); exit(1);}

   // Find t for this alpha using Bolzano's theorem
   float  t1=0,  t2=10;
   float ft1=1, ft2=0;
   while (2*(ft1-ft2)/(ft1+ft2)>0.01) { // The difference is greater than 10%
      float t3=(t1+t2)/2;
      float ft3=gaus_outside_x0(t3);
      if (ft3<alpha) {t2=t3; ft2=ft3;}
      else           {t1=t3; ft1=ft3;}
   }
   float t=(t1+t2)/2;

   // Show number of samples to be estimated
   float n=(t*stddev)/(r*avg); n *=n;
   cout << "You should take at least " << n << " samples\n";
}

/* ------------------------------------------------------------------------- */
/* Usage                                                                     */
/* ------------------------------------------------------------------------- */
void Usage() {
      cerr << "Usage: sample_size -stats <avg1+-dev1> -alpha <alpha=0.005> -r <r=0.01>\n";
}

/* Colimate ---------------------------------------------------------------- */
/*Colimate:
   PROGRAM Sample_size {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Sample_size/Help/sample_size.html";
      help="Compute the sample size so that the true mean is within a certain
            acuracy with a cretain probability of unluckiness";
      OPEN MENU menu_sample_size;
      COMMAND LINES {
         + usual: xmipp_sample_size -stats $STATS1 -alpha $ALPHA
                     -r $R
      }
      PARAMETER DEFINITIONS {
         $STATS1 {
            label="Mean and standard deviation";
            help="Format: avg+-stddev";
            type=text;
         }
         $ALPHA {
            label="Probability of unluckiness";
            help="An unlucky selection is admitted with this probability";
            type=float;
            by default=0.005;
         }
         $R {
            label="Accuracy";
            help="The true mean must be within this percentage of the
                  estimated value";
            type=float;
            by default=0.01;
         }
      }
   }
   
   MENU menu_sample_size {
      $STATS1
      $ALPHA
      $R
   }

*/
