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
   vector<float>  avg, stddev, N;
   int ExpNo;
   
// Read command line -------------------------------------------------------
   // Check number of arguments
   if (argc%2!=1 || argc==1) {Usage();exit(1);}

   // Read arguments
   try {
      ExpNo=(argc-1)/2;
      for (int i=0; i<ExpNo; i++) {
         string str=argv[2*i+2]; remove_plus_minus(str);
         N     .push_back(AtoF(argv[2*i+1]));
         avg   .push_back(AtoF(first_token(str)));
         stddev.push_back(AtoF(next_token()));
      }
   } catch (Xmipp_error XE) {cerr << XE; exit(1);}

// Compute New statistics --------------------------------------------------
   double sum1=0, sum2=0;
   double Ntot=0;
   for (int i=0; i<ExpNo; i++) {
      Ntot += N[i];
      sum1 += N[i]*avg[i];
      sum2 += N[i]*(stddev[i]*stddev[i]+avg[i]*avg[i]);
   }
   
   double newavg=sum1/Ntot;
   double newstddev=sqrt(sum2/Ntot-newavg*newavg);
   cout << "After combining " << Ntot << " " << newavg << "+-" 
        << newstddev << endl;
}

/* ------------------------------------------------------------------------- */
/* Usage                                                                     */
/* ------------------------------------------------------------------------- */
void Usage() {
      cerr << "Usage: combine_stats <N1> <avg1+-dev1> ...\n";
}

/* Menu -------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Combine_stats {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Combine_stats/Help/combine_stats.html";
      help="Combine statistical properties of two populations";
      OPEN MENU Combine_stats;
      COMMAND LINES {
         + usual: combine_stats $N1 $STATS1 $N2 $STATS2
      }
      PARAMETER DEFINITIONS {
         $N1 {label="Experiments on Population 1"; type=NATURAL;}
         $N2 {label="Experiments on Population 2"; type=NATURAL;}
         $STATS1 {label="Statistics of Population 1"; help="avg+-stddev";
            type=TEXT;}
         $STATS2 {label="Statistics of Population 2"; help="avg+-stddev";
            type=TEXT;}
      }
   }
   MENU Combine_stats {
      $N1
      $STATS1
      $N2
      $STATS2
   }
*/
