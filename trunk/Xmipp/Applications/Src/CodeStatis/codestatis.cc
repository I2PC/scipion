/* author: ALberto Pascual
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
#include <stdio.h>
#include <stdlib.h>
#include <algo.h>
#include <vector>


typedef vector< vector< float > > coordinates;   

double dist( vector<float> c1, vector<float> c2) {
  double dummy, sum2;
  sum2 = 0;
  for (unsigned i=0; i<3; i++) {
    dummy = c1[i]-c2[i];
    sum2 += dummy*dummy;
  }  
  return sqrt(sum2);
} 

int main(int argc, char **argv) {


  FILE *ifp, *ofp, *tmpfp, *ocfp, *varfp; 
  char *iname, *oname, *tmpname, *ocname, *varname;
  char *str;
  int dim, index;
  float rmsf = 0;
  float spacing;
  vector <coordinates> coordSet;
  coordinates tmpCoord, avecoord;
        
  // Read arguments
  
  try { 
    iname = get_param(argc, argv, "-iname");
    oname = get_param(argc, argv, "-soname");
    ocname = get_param(argc, argv, "-aoname");
    varname = get_param(argc, argv, "-varname");
    spacing = AtoF(get_param(argc, argv, "-spacing", "1"));
  } 
  catch (Xmipp_error) {
    cout << "Usage:" << endl;
    cout << "-iname       : Input file name (contain files names)" << endl;
    cout << "-soname      : Output file name (sorted code vectors)" << endl;
    cout << "-aoname      : Output file name (average code vectors)" << endl;
    cout << "-varname     : Output file name (variation file)" << endl;
    cout << "-spacing     : Spacing (in Angstroms/pixel)" << endl;
    exit(1);
   }
    
  cout << "Given parameters are: " << endl;
  cout << "iname = " << iname << endl;
  cout << "soname = " << oname << endl;
  cout << "aoname = " << ocname << endl;
  cout << "varname = " << varname << endl;
  cout << "spacing = " << spacing << endl;


  str = new char [100];
 
  ifp=fopen(iname,"r");      
  ofp=fopen(oname,"w");      
  ocfp=fopen(ocname,"w");      
  varfp=fopen(varname,"w");      

      
  while (feof(ifp) == 0) {   		// Read all files names

     fscanf(ifp, "%s", str);
     cout << "Processing " << str << " file ....." << endl;
          
     tmpfp=fopen(str,"r");                 // Read first file      
      
     fgets(str, 100, tmpfp);               // Read dimensions
     vector <float> v1;  
     read_float_list(str, 1 ,v1);          // Read 1 values from line in v1
     dim = (int) v1[0];
  
     cout << "dimension : " << dim << endl;

     cout << "reading file....." << endl; 
    
     coordinates coord;

     while (fgets(str, 100, tmpfp) != NULL) {
       vector <float> v;
       read_float_list(str,dim,v);            // Read 3 values from line in v
       coord.push_back(v);            
     }
     
     coordSet.push_back(coord);
          
     fclose(tmpfp);
  }    

//   cout << "sorting...." << endl;
//   for (unsigned s = 0; s < coordSet.size(); s++)
 //    sort(coordSet[s].begin(), coordSet[s].end()); 


   cout << "saving new coordinates....." << endl; 
   
   fprintf(ofp, "%d \n", dim);
   fprintf(ocfp, "%d \n", dim);
   
   for(unsigned p = 0; p < coordSet[0].size(); p++) {    
     cout  << "finding point "<< p << endl; 
     double deltaX = 0, deltaY = 0, deltaZ = 0; 

     /* Save First file coordinate (reference) */
     fprintf(ofp, "*****************\n");       
     for (unsigned d=0; d < dim; d++)
       fprintf(ofp, "%f ", coordSet[0][p][d]);
     fprintf(ofp, "\n");       
      
     tmpCoord.resize(0);
     deltaX += coordSet[0][p][0];
     deltaY += coordSet[0][p][1];
     deltaZ += coordSet[0][p][2];
     tmpCoord.push_back(coordSet[0][p]);

     for(unsigned c = 1; c < coordSet.size(); c++) {                 
       
       double minDist = MAXFLOAT;
       int minDistIndex = 1;
//       for(unsigned p2 = 0; p2 < coordSet[0].size(); p2++) {    
       for(unsigned p2 = 0; p2 < coordSet[c].size(); p2++) {    
          double tmp;
          tmp = dist(coordSet[0][p], coordSet[c][p2]); 
          if (tmp < minDist) {
            minDist = tmp; 
	    minDistIndex = p2;
	  }
       }

       deltaX += coordSet[c][minDistIndex][0];
       deltaY += coordSet[c][minDistIndex][1];
       deltaZ += coordSet[c][minDistIndex][2];
       tmpCoord.push_back(coordSet[c][minDistIndex]);
       for (unsigned d=0; d < dim; d++)
         fprintf(ofp, "%f ", coordSet[c][minDistIndex][d]);
       fprintf(ofp, "\n");       
       
       // deletes the found coordinate.
       coordSet[c].erase(coordSet[c].begin() + minDistIndex);       
		 
     }    
   
     deltaX /= coordSet.size();
     deltaY /= coordSet.size();
     deltaZ /= coordSet.size();
     
     fprintf(ocfp, "%f %f %f\n", deltaX, deltaY, deltaZ);
     fprintf(ofp, "Average: %f %f %f\n", deltaX, deltaY, deltaZ);
     
     vector<float> tmp;
     tmp.push_back(deltaX); 
     tmp.push_back(deltaY); 
     tmp.push_back(deltaZ); 
     avecoord.push_back(tmp);
     
     // Calculate standard deviation
     
     float sdX = 0;
     float sdY = 0;
     float sdZ = 0;
     float rmscum = 0, rmsd;
     for (unsigned i = 0; i < tmpCoord.size(); i++) {
        sdX += (tmpCoord[i][0] - deltaX)*(tmpCoord[i][0] - deltaX);
        sdY += (tmpCoord[i][1] - deltaY)*(tmpCoord[i][1] - deltaY);
        sdZ += (tmpCoord[i][2] - deltaZ)*(tmpCoord[i][2] - deltaZ);
     }
     rmscum = (sdX + sdY + sdZ);
     sdX = (sdX/(tmpCoord.size()));
     sdY = (sdY/(tmpCoord.size()));
     sdZ = (sdZ/(tmpCoord.size()));
     rmscum /= (float) (tmpCoord.size());
     rmsd = sqrt(rmscum);

     fprintf(ofp, "Variance   : %f %f %f\n", sdX, sdY, sdZ);
     fprintf(ofp, "SD         : %f %f %f\n", sqrt(sdX), sqrt(sdY), sqrt(sdZ));
     fprintf(ofp, "RMSE       : %f\n", rmsd);
     fprintf(varfp, "%f %f %f %f\n", spacing*sqrt(sdX), spacing*sqrt(sdY), spacing*sqrt(sdZ), spacing*rmsd);
     rmsf += rmsd;
 }  

  fprintf(ofp, "*****************\n");       
  
  cout << "calculating some stats....." << endl; 
  
  rmsf /= (float) coordSet[0].size();
  fprintf(ofp, "RMS Fluctuation in voxels : %f\n", rmsf);
  fprintf(ofp, "RMS Fluctuation in Angstroms : %f\n", spacing*rmsf);
  
  double minDist = MAXFLOAT;
  double maxDist = 0;
  double aveDist = 0;
  long count = 0;
  for (int i = 0; i < avecoord.size() - 1; i++) {
    for (int j = 0; j < avecoord.size(); j++) {
          if (i == j) continue; 
          double tmp;
          tmp = dist(avecoord[i], avecoord[j]); 
          if (tmp < minDist) 
            minDist = tmp;
          if (tmp > maxDist)
            maxDist = tmp;
    }
    aveDist += minDist;   
    count++;
  }
  
  aveDist /= count;  
/*  fprintf(ofp, "Minimum inter-neuron distance in voxels : %f\n", minDist);
  fprintf(ofp, "Minimum inter-neuron distance in Angstroms : %f\n", spacing*minDist);
  fprintf(ofp, "Maximum inter-neuron distance in voxels : %f\n", maxDist);
  fprintf(ofp, "Maximum inter-neuron distance in Angstroms : %f\n", spacing*maxDist);*/
  fprintf(ofp, "Average inter-neuron distance in voxels : %f\n", aveDist);
  fprintf(ofp, "Average inter-neuron distance in Angstroms : %f\n", spacing*aveDist);
    
  cout << "done....." << endl; 
         
  fclose(ifp);       // close input file
  fclose(ofp);       // close output file
  fclose(ocfp);      // close output file
  fclose(varfp);      // close output file
   exit(0);
      
} 


