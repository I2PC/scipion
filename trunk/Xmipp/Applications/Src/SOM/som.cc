/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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

/* Part of this code were developed by Lorenzo Zampighi and Nelson Tang
   of the department of Physiology of the David Geffen School of Medicine,
   University of California, Los Angeles   
*/


// To avoid problems with long template names
#pragma warning(disable:4786)

#include <fstream>
#include <Classification/xmippSOM.hh>

/* Prototypes -============================================================= */

void Usage (char **argv);

/* Main function -============================================================= */

main(int argc, char** argv) {


/* Input Parameters ======================================================== */

FileName       fn_in;    	// input file
FileName       fn_out;   	// output file
FileName       cb_in = "";   	// Code vectors input file
FileName       tmpN;		// Temporary variable
unsigned       iter = 10000;	// Iteration number
unsigned       verb = 0;	// Verbosity level
bool           norm = 1;	// Normalize?
unsigned       xdim = 10;	// X-dimension (-->)
unsigned       ydim = 5;	// Y-dimension 
float 	       alpha_0 = 0.1;	// Initial alpha value
float 	       radius_0;	// Initial radius value
string         layout = "HEXA";	// Type of layout (Topology)
bool 	       gaussian = true; // Gaussian kernel
bool 	       bubble = false;  // Bubble kernel
bool 	       saveClusters = false;    // Save clusters in separate files
bool 	       use_rand_cvs = false; // NT: flag to truly randomize codevectors or not


/* Parameters ============================================================== */
   try {

       fn_in = get_param(argc, argv, "-din");

       if (check_param(argc, argv, "-cout"))
          fn_out = get_param(argc, argv, "-cout");
       else {
         Usage(argv); 
	 exit(EXIT_FAILURE);
       } 

       if (check_param(argc, argv, "-cvin"))
       	  cb_in = get_param(argc, argv, "-cvin");

       ydim = AtoI(get_param(argc, argv, "-ydim", "5"));
       xdim = AtoI(get_param(argc, argv, "-xdim", "10"));

       if (check_param(argc, argv, "-hexa")) {
       	  if (check_param(argc, argv, "-rect")) {
	        cout << "Error: you can not define two topologies" << endl;
	 	exit(EXIT_FAILURE);	  	
	  }
	  layout = "HEXA";
       } else if (check_param(argc, argv, "-rect"))
          layout = "RECT";
         

       if (check_param(argc, argv, "-gaussian")) {
       	  if (check_param(argc, argv, "-bubble")) {
	        cout << "Error: you can not define two kernels" << endl;
	 	exit(EXIT_FAILURE);	  	
	  }
	  gaussian = true;
          bubble = false;
       } else if (check_param(argc, argv, "-bubble")) {
          bubble = true;
	  gaussian = false;
       }


       alpha_0 = AtoF(get_param(argc, argv, "-alpha", "0.1"));

       if (check_param(argc, argv, "-radius"))
       	 radius_0 = AtoF(get_param(argc, argv, "-radius"));
       else  
         if (xdim > ydim)
	   radius_0 = xdim;
	 else 
	   radius_0 = ydim;
       
       iter = AtoI(get_param(argc, argv, "-iter", "10000"));
       verb = AtoI(get_param(argc, argv, "-verb", "0"));

       if (check_param(argc, argv, "-norm"))
          norm = true;
       else norm = false;
       
       if (check_param(argc, argv, "-saveclusters"))
          saveClusters = true;
       else saveClusters = false;

       if (check_param(argc, argv, "-randomcodevectors"))
	 use_rand_cvs = true;
       else use_rand_cvs = false;

       if (argc == 1) {Usage(argv);}
       
   }
   catch (Xmipp_error XE) {cout << XE; Usage(argv);}


/* Some validations ===================================================== */
  

   if (iter < 1) {
     cerr << argv[0] << ": invalid value for iter (must be > 1): " << iter << endl;
     exit(EXIT_FAILURE);
   }
   
   if (verb < 0 || verb > 2) {
     cerr << argv[0] << ": invalid value for verbosity (must be between 0 and 2): " << verb << endl;
     exit(EXIT_FAILURE);
   }

   if (alpha_0 < 0) {
     cerr << argv[0] << ": invalid value for alpha (must be > 0): " << alpha_0 << endl;
     exit(EXIT_FAILURE);
   }

   if (radius_0 < 1) {
     cerr << argv[0] << ": invalid value for radius (must be > 1): " << radius_0 << endl;
     exit(EXIT_FAILURE);
   }
   
   if (xdim < 1) {
     cerr << argv[0] << ": invalid value for xdim (must be > 1): " << xdim << endl;
     exit(EXIT_FAILURE);
   }

   if (ydim < 1) {
     cerr << argv[0] << ": invalid value for ydim (must be > 1): " << ydim << endl;
     exit(EXIT_FAILURE);
   }


/* Shows parameters ===================================================== */

   cout << endl << "Parameters used: " << endl;
   cout << "Input data file : " << fn_in << endl;
   cout << "Output file name : " << fn_out << endl;
   if (cb_in != "")
   	cout << "Input code vectors file name : " << cb_in << endl;   
   else if (use_rand_cvs)
     cout << "Using randomized code vectors" << endl;
   if (saveClusters)
     cout << "Save clusters in separate files: " << fn_out << ".(cluster number)"<< endl;
   cout << "Horizontal dimension (Xdim) = " << xdim << endl;
   cout << "Vertical dimension (Ydim) = " << ydim << endl;
   if (layout == "HEXA")
     cout << "Hexagonal topology " << endl;
   else
     cout << "Rectangular topology " << endl;
 
   if (gaussian)
     cout << "Gaussian neighborhood function " << endl;
   else
     cout << "Bubble neighborhood function" << endl;
   cout << "Initial learning rate (alpha) = " << alpha_0 << endl;
   cout << "Initial neighborhood radius (radius) = " << radius_0 << endl;
 
   cout << "Total number of iterations = " << iter << endl;
   cout << "verbosity level = " << verb << endl;
   if (norm)
     cout << "Normalize input data" << endl;
   else
     cout << "Do not normalize input data " << endl;


/* Open training vector ================================================= */

  
  ifstream inStream(fn_in.c_str());
  if (!inStream) {
      cerr << argv[0] << ": can't open file " << fn_in << endl;
      exit(EXIT_FAILURE);
  }

  xmippCTVectors ts(0, true);

  cout << endl << "Reading data file " << fn_in << "....." << endl;     
  inStream >> ts;


/* Real stuff ============================================================== */


 try {
   if (norm) {
   	cout << "Normalizing....." << endl;
   	ts.normalize();  		    // Normalize input data        
   }	
      
   xmippMap* myMap;
   
   if (cb_in != "") {
        cout << "Reading codevectors file " << cb_in << "....." << endl;   
        ifstream codeStream(cb_in.c_str());
        if (!codeStream) {
          cerr << argv[0] << ": can't open file " << cb_in << endl;
          exit(EXIT_FAILURE);
        }
        myMap = new xmippMap(codeStream);
   } else    
        myMap = new xmippMap(layout, xdim, ydim, ts, use_rand_cvs);


   xmippSOM *thisSOM;
   Descent alpha(alpha_0, 0);         // alpha decreases linearly to 0
   Descent radius(radius_0, 1);	      // radius decreases linearly to 1			
   
   xmippSOM::neighType neigh;
   if (gaussian) 
   	neigh = xmippSOM::GAUSSIAN;
   else 
   	neigh = xmippSOM::BUBBLE;   
   thisSOM = new xmippSOM(alpha, radius, neigh, iter);    // Creates SOM Algorithm

   xmippTextualListener myListener;	    // Define the listener class
   myListener.setVerbosity() = verb;	    // Set verbosity level
   thisSOM->setListener(&myListener);       // Set Listener

   thisSOM->train(*myMap, ts);              // Train algorithm

   // Test algorithm
   cout << endl;
   double dist = thisSOM->test(*myMap, ts);
   cout << endl << "Quantization error : " <<  dist << endl;

   // Classifying
   cout << "Classifying....." << endl;
   myMap->classify(&ts);

   // Calibrating
   cout << "Calibrating....." << endl;
   myMap->calibrate(ts);

   if (norm) {
   	cout << "Denormalizing code vectors....." << endl;
   	myMap->unNormalize(ts.getNormalizationInfo()); // de-normalize codevectors        
   }	

  /******************************************************* 
      Saving all kind of Information 
  *******************************************************/

   cout << "Saving code vectors as " << fn_out << ".cod ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".cod"; 
   ofstream codS(tmpN.c_str());
   codS << *myMap;
   codS.flush();    

   cout << "Saving algorithm information as " << fn_out << ".inf ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".inf"; 
   ofstream infS(tmpN.c_str());
   infS << "Kohonen SOM algorithm" << endl << endl;
   infS << "Input data file : " << fn_in << endl;
   if (cb_in != "")
     infS << "Input code vectors file : " << cb_in << endl;
   else if (use_rand_cvs)
     infS << "Using randomized code vectors" << endl;
   infS << "Code vectors output file : " << fn_out <<  ".cod" << endl;
   infS << "Algorithm information output file : " << fn_out <<  ".inf" << endl;
   infS << "Number of feature vectors: " << ts.size() << endl;
   infS << "Number of variables: " << ts.itemAt(0).size() << endl;
   infS << "Horizontal dimension (Xdim) = " << xdim << endl;
   infS << "Vertical dimension (Ydim) = " << ydim << endl;
   if (layout == "HEXA")
     infS << "Hexagonal topology " << endl;
   else
     infS << "Rectangular topology " << endl;
   if (gaussian)
     infS << "Gaussian neighborhood function " << endl;
   else
     infS << "Bubble neighborhood function" << endl;
   infS << "Initial learning rate (alpha) = " << alpha_0 << endl;
   infS << "Initial neighborhood radius (radius) = " << radius_0 << endl;
   infS << "Total number of iterations = " << iter << endl;
   if (norm)
     infS << "Input data normalized" << endl;
   else
     infS << "Input data not normalized" << endl;
   infS << "Quantization error : " <<  dist << endl;

   infS.flush();    

   // assign data to neurons
   if (saveClusters) {
   	cout << "Saving neurons assigments ....." << endl;  
   	for (unsigned i= 0; i < myMap->size(); i++) {
		tmpN = fn_out.c_str() + (string) "."  + ItoA(i); 
   		ofstream cStream(tmpN.c_str());
		for (int j = 0; j < myMap->classifAt(i).size(); j++)
   			cStream << myMap->classifAt(i)[j] << endl;
   		cStream.flush();    
   	}
   }

   // save .vs file to be compatible with SOM_PAK
   cout << "Saving visual file as " << fn_out << ".vs ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".vs"; 
   ofstream vsStream(tmpN.c_str());
   vsStream << ts.theItems[0].size() << " " << myMap->layout() << " " << myMap->width() << " " << myMap->height();
   if (gaussian)
      vsStream  << " gaussian" << endl;
   else    
      vsStream  << " bubble" << endl;
   for (int i= 0; i < ts.size(); i++) {
   	int j = myMap->winner(ts, i);
   	vsStream << myMap->indexToPos(j).first << " " << myMap->indexToPos(j).second << " " << eDist(myMap->theItems[j], ts.theItems[i]) << " " << ts.theTargets[i] << endl;
   }   
   vsStream.flush();    

   // save .his file (Histogram)
   cout << "Saving code vectors histogram file as " << fn_out << ".his ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".his"; 
   ofstream hisStream(tmpN.c_str());
   myMap->printHistogram(hisStream);
   hisStream.flush();    

   // save .err file (Average Quantization Error)
   cout << "Saving code vectors average quantization error file as " << fn_out << ".err ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".err"; 
   ofstream errStream(tmpN.c_str());
   myMap->printQuantError(errStream);
   errStream.flush();    


   cout << endl;
   
   delete myMap;
   
 } catch ( const exception& e ) {
    cout << e.what() << endl;
 }
  return 0;
}


/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage (char **argv) {
  printf (
     "\nUsage: %s [Purpose and Parameters]"
     "\nPurpose: Kohonen Self-Organizing Feature Maps"
     "\n"           
     "\nParameter Values: (note space before value)"
     "\n"
     "\n    -din    file_in           Input data file (plain data)"
     "\n    -cout   file_out          Base name for output data files:"
     "\n    -cvin   file_in           Codevectors input file"
     "\n    -saveclusters    	      save clusters in separate files (Default = No)"     
     "\n    -xdim   H-dimension	      Horizontal size of the map (default = 10)"
     "\n    -ydim   V-dimension       Vertical size of the map (default = 5)"
     "\n    -hexa    	   	      Hexagonal topology (default)"     
     "\n    -rect    	   	      Rectangular topology"     
     "\n    -gaussian   	      Gaussian neighborhood learning kernel (Default)"     
     "\n    -bubble    		      Bubble neighborhood learning kernel"     
     "\n    -alpha     learning rate  Initial learning rate value (default = 0.1)"
     "\n    -radius    radius 	      Initial neighborhood radius"     
     "\n    			      (default = max(xdim, ydim))"     
     "\n    -iter      iterations     Number of iterations (default = 10000)"     
     "\n    -norm            	      Normalize training data (default)"     
     "\n    -verb      verbosity      Information level while running: "     
     "\n    			      0: No information (default)"     
     "\n    			      1: Progress bar"     
     "\n    			      2: Code vectors change between iterations"     
     "\n    -randomcodevectors        Use truly randomized codevectors (default FALSE)"
     "\n			   \n"
     ,argv[0]);
     exit(0);
}
