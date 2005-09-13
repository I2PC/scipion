/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
 *              
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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

#include "../xmippRM.hh"
#include "../xmippDistances.hh"


/**
* Calculate the Random matrix
* @param ts The vectors.
*/
void xmippRM::calculateRM(xmippCTVectors const &ts, int _k){

   int size=ts.size();
   int dim = ts.itemAt(0).size();
   RM.resize(_k, dim);
   randomize_random_generator();

   int verbosity = listener->getVerbosity();

   if (verbosity)
  	listener->OnReportOperation((string) "Creating random matrix....\n");  

   // Define de random matrix      

   if (matdist == "sparse") {
     if (verbosity == 1)    
       listener->OnInitOperation(_k);
     for (int i = 0; i < _k; i++) {
       for (int j = 0; j < dim; j++) {
          double cm = sqrt(3.0);
          float rn;
	  rn = rnd_unif();
	  if (rn < 0.16666) cm *= 1;
	  else if ((rn >= 0.16666) && (rn < 0.8333)) cm = 0;
	  else cm *= -1;
          RM(i,j) = cm;
       }
       if (verbosity == 1)    
         listener->OnProgress(i);
     }
     if (verbosity == 1)    
        listener->OnProgress(_k);
   } else {
      RM.init_random(0, 1, "gaussian");
   }

    // Normalize columns to unit length
   if (verbosity)
  	listener->OnReportOperation((string) "Normalizing random matrix....\n");  
     if (verbosity == 1)    
       listener->OnInitOperation(dim);
    
    // first calculates the lenght of the vector
    for (int z = 0; z < dim; z++) { 	
	xmippFeature lenght = 0;
	for(int it = 0; it < _k; it++) {
		if (finite(RM(it, z))) {
			lenght += RM(it, z)*RM(it, z);
		}
			
	}
	lenght = sqrt(lenght);

	if (lenght != 0) {
	  // Now normalize the vector
	  for(int it = 0; it < _k; it++) {
		if (finite(RM(it, z)))
			RM(it, z) /= lenght;
	  }
	}
	if (verbosity == 1)    
           listener->OnProgress(z);
     } // z 	
     if (verbosity == 1)    
        listener->OnProgress(dim);
}


/* Project ----------------------------------------------------------------- */
void xmippRM::Project(xmippCTVectors &input, xmippCTVectors &output, int _k) {

   // Calculate RM
   calculateRM(input, _k); 

  // Do the projection
   
  int verbosity = listener->getVerbosity();
  if (verbosity)
  	listener->OnReportOperation((string) "Projecting data....\n");  
  if (verbosity == 1)    
       listener->OnInitOperation(_k);
  
  long size = input.size();
  long dim = input.theItems[0].size();
  
  output.theItems.resize(size);
  output.theTargets.resize(size);
  for (int h = 0; h < size; h++) output.theItems[h].resize(_k, 0); 

    for (int i = 0; i < _k; i++) {
       for (int z = 0; z < size; z++) {
         double cum = 0;
         for (int j = 0; j < dim; j++) {
           cum += RM(i,j)*input.theItems[z][j];  
         } // j
         output.theItems[z][i] = cum; 
       }  // z
       if (verbosity == 1)    
           listener->OnProgress(i);
    } // i
    if (verbosity == 1)    
       listener->OnProgress(_k);
}

/* Clear ------------------------------------------------------------------- */
void xmippRM::clear() {
   RM.clear();
   matdist = "gaussian";
}


