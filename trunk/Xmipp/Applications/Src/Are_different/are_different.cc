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
// Remade by R. Marabini 12 Junio and 14 Julio 2001
//This program performs:
//  * Small-sample Comparations of two means (based on *independend* 
//    samples from normal distributions with a common variance)
//  * Inference for the mean of a *Paired* difference
/********************************************************************  

Experiments are usually done to compare the average response or the proportion of
responses among two or more groups. For example, you may wish to compare the mean
weight gain of two different feeds given to livestock in order to decide which is
better. In these cases, there is no particular value of interest for each
feed-type; rather the difference in performance is of interest.

Paired t-test: In many applications it is natural to make two measurements of essencially the same kind. For example, measure a FOM over two reconstructions obtained with two different reconstruction algorithms applied to the same data set. This test checks if we can reject the hipothesys that the two algorithms performs equally. NOtice that we do not probe that the two algorithm perfors equally.

Small-sample Comparations...: this method is similar to the above one but does not assume that the data is paired. the input is just the mean and the standard deviation of both populations.

*************************************************************************/
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
   float min_significance, t;
   enum TEST {Two_Ind_Means, Two_Paired_Means};
   TEST test=Two_Ind_Means;
   DocFile  values;    // input file 
   double confidence;
   double mean_data_d=0., sq_sigma_d=0.;
   int number_of_samples;
         
// Read command line -------------------------------------------------------
   if (argc<3) {Usage(); exit(1);}

   try {
      if(check_param(argc, argv, "-i")) test=Two_Paired_Means;
      switch(test)
      {
      case Two_Paired_Means:
         values.read(get_param(argc,argv,"-i"));
         break;
         
      case Two_Ind_Means:  
         for (int i=0; i<2; i++) {
            string str=argv[2*i+2]; remove_plus_minus(str);
            N     .push_back(AtoF(argv[2*i+1]));
            avg   .push_back(AtoF(first_token(str)));
            stddev.push_back(AtoF(next_token()));
         }
         break;
      }//switch end
      min_significance=AtoF(get_param(argc,argv,"-min_sig","0.95"));
   } catch (Xmipp_error XE) {Usage(); exit(1);}

//perform the test
//Inference for tthe Mean of a Paired Difference
   switch(test)
   {
   case Two_Paired_Means:
         number_of_samples=values.dataLineNo();//number of data entries
         values.go_first_data_line();
         for( int i=0; i< number_of_samples; i++)
            {
            mean_data_d += (values(0)-values(1));
            values.next_data_line();            
            }
         mean_data_d /=  number_of_samples;
         values.go_first_data_line();
         for( int i=0; i< number_of_samples; i++)
            {
            sq_sigma_d  += ((values(0)-values(1))-mean_data_d)*
                           ((values(0)-values(1))-mean_data_d);
            values.next_data_line();            
            }
         sq_sigma_d  /=  (number_of_samples-1);
         sq_sigma_d   =  sqrt(sq_sigma_d);
         t = (float)(mean_data_d/ (sq_sigma_d/sqrt((double)number_of_samples)));
         t = ABS(t);
         confidence=student_within_t0(t,number_of_samples-1.);
         break;
//Small-Sample Comparation of Two means (Based on independent Samples 
//from a normal distribution with a common variance
   case Two_Ind_Means:         
// They are different with confidence --------------------------------------
         double var0=stddev[0]*stddev[0];
         double var1=stddev[1]*stddev[1];
         double varp= ((N[0]-1)*var0+(N[1]-1)*var1) / (N[0]+N[1]-2);
         t= (avg[0]-avg[1])/sqrt(varp*(1/N[0]+1/N[1]));

            confidence=student_within_t0(t,(N[0]+N[1]-2));
         break;
    }//switch end
   if (confidence>=min_significance)
      cout << "I would say they are different\n"
           << "Up to a significance level of "
           << confidence << endl;   
   else
      cout << "I would NOT say they are the different\n"
           << "Unless you relax your Significance level to " << confidence
           << endl;   
}

/* ------------------------------------------------------------------------- */
/* Usage                                                                     */
/* ------------------------------------------------------------------------- */
void Usage() {
      cerr << "If t_test is desired:\n" 
           << "\tUsage: are_different <N1 avg1+-dev1> <N2 avg2+-dev2> \n"
           << "\t\t[-min_sig <threshold=0.95>]  : minimum significance\n"
           << "\t\t                               to say they are different\n";
      cerr << "If continuos random variable test is desired:\n" 
           << "\tUsage: -i doc_file_with_data \n"
           << "\t\t[-min_sig <threshold=0.95>]  : minimum significance\n"
           << "\t\t                               to say they are different\n";
}

/* Menu -------------------------------------------------------------------- */
/*Colimate part missing
*/
