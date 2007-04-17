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
#include <Classification/xmippFuzzySOM.hh>

/* Prototypes -============================================================= */

void Usage (char **argv);

/* Main function -============================================================= */

main(int argc, char** argv) {


/* Input Parameters ======================================================== */

FileName       fn_in;    	// input file
FileName       fn_out;   	// output file
FileName       cb_in = "";   	// Code vectors input file
FileName       fn_algo_in = ""; // input algorithm file
FileName       tmpN;		// Temporary variable
double         eps = 1e-7;	// Stopping criteria
unsigned       iter = 1000;	// Iteration number
unsigned       verb = 0;	// Verbosity level
bool           norm = 1;	// Normalize?
unsigned       xdim;		// X-dimension (-->)
unsigned       ydim;		// Y-dimension 
double	       m0 = 2.0;	// Initial m
double         m1 = 1.01;	// Final m
double         reg;		// Regularization (smoothness) parameter
string         layout = "RECT"; // topology (layout)
unsigned       annSteps = 1000; // Deterministic Annealing steps
bool 	       saveClusters = false;    // Save clusters in separate files
bool use_rand_cvs = false; // NT: flag to truly randomize codevectors or not

/* Parameters ============================================================== */
   try {

       if (check_param(argc, argv, "-din"))
         fn_in = get_param(argc, argv, "-din");
       else {
         Usage(argv); 
	 exit(EXIT_FAILURE);
       } 

       if (check_param(argc, argv, "-cout"))
          fn_out = get_param(argc, argv, "-cout");
       else {
         Usage(argv); 
	 exit(EXIT_FAILURE);
       } 

       if (check_param(argc, argv, "-cvin"))
       	  cb_in = get_param(argc, argv, "-cvin");


       if (check_param(argc, argv, "-xdim"))
	 xdim= AtoI(get_param(argc, argv, "-xdim"));
       else {
         Usage(argv); 
	 exit(EXIT_FAILURE);
       } 

       if (check_param(argc, argv, "-ydim"))
	 ydim= AtoI(get_param(argc, argv, "-ydim"));
       else {
         Usage(argv); 
	 exit(EXIT_FAILURE);
       } 

       if (check_param(argc, argv, "-hexa")) {
       	  if (check_param(argc, argv, "-rect")) {
	        cout << "Error: you can not define two topologies" << endl;
	 	exit(EXIT_FAILURE);	  	
	  }
	  layout = "HEXA";
       } else if (check_param(argc, argv, "-rect"))
          layout = "RECT";
         
       m0 =  AtoF(get_param(argc, argv, "-m0", "2.0"));
       m1 =  AtoF(get_param(argc, argv, "-m1", "1.01"));
       reg =  AtoF(get_param(argc, argv, "-reg", "0.5"));

       eps = AtoF(get_param(argc, argv, "-eps", "1e-7"));
       iter = AtoI(get_param(argc, argv, "-iter", "1000"));
       verb = AtoI(get_param(argc, argv, "-verb", "0"));

       if (check_param(argc, argv, "-norm"))
          norm = true;
       else norm = false;

       annSteps = AtoI(get_param(argc, argv, "-steps", "1000"));
       
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

   if (m0 <= 1) {
     cerr << argv[0] << ": invalid value for m0 (must be > 1): " << m0 << endl;
     exit(EXIT_FAILURE);
   }

   if (m1 <= 1) {
     cerr << argv[0] << ": invalid value for m1 (must be > 1): " << m1 << endl;
     exit(EXIT_FAILURE);
   }

   if ((annSteps != 0) && (m0 <= m1) ) {
     cerr << argv[0] << ": invalid value for m0 and m1 (m0 must be > m1) " << endl;
     exit(EXIT_FAILURE);
   }

   if ((annSteps < 0) || (annSteps == 1)) {
     cerr << argv[0] << ": invalid value for annSteps (must be > 1): " << annSteps << endl;
     exit(EXIT_FAILURE);
   }

   if (reg < 0) {
     cerr << argv[0] << ": invalid value for smoothness parameter (must be > 0): " << reg << endl;
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

  if (reg == 0) reg = 1.23456789e-7;


/* Shows parameters ===================================================== */

   cout << endl << "Parameters used: " << endl;
   cout << "Input data file : " << fn_in << endl;
   cout << "Output code vector file : " << fn_out << ".cod" << endl;
   if (cb_in != "")
   	cout << "Input code vectors file name : " << cb_in << endl;
   else if (use_rand_cvs)
     cout << "Using randomized code vectors" << endl;
   cout << "Horizontal dimension (Xdim) = " << xdim << endl;
   cout << "Vertical dimension (Ydim) = " << ydim << endl;
   if (layout == "HEXA")
     cout << "Hexagonal topology " << endl;
   else
     cout << "Rectangular topology " << endl;
   cout << "Initial Fuzzy Constant (m0) = " << m0 << endl;
   cout << "Final Fuzzy Constant (m1) = " << m1 << endl;
   cout << "Regularization Constant (reg) = " << reg << endl;
   cout << "Deterministic annealing steps = " << annSteps << endl;
   cout << "Total number of iterations = " << iter << endl;
   cout << "Stopping criteria (eps) = " << eps << endl;
   cout << "verbosity level = " << verb << endl;
   if (norm)
     cout << "Normalize input data" << endl;
   else
     cout << "Do not normalize input data " << endl;


/* Open training vector ================================================= */

  cout << endl << "Reading file " << fn_in << "....." << endl;
  
  ifstream inStream(fn_in.c_str());
  if (!inStream) {
      cerr << argv[0] << ": can't open file " << fn_in << endl;
      exit(EXIT_FAILURE);
  }

  xmippCTVectors ts(0, true);
  try 
    {
      inStream >> ts;
    }
  catch (exception& e)
    {
      cerr << argv[0] << ": can't read file " << fn_in  << " because " << e.what() << endl;
      exit(EXIT_FAILURE);
    }


/* Real stuff ============================================================== */


 try {

   if (norm) {
   	cout << "Normalizing....." << endl;
   	ts.normalize();  		    // Normalize input data        
   }	

   xmippFuzzyMap *myMap;
   
   if (cb_in != "") {
        cout << "Reading fuzzy codevectors file " << cb_in << "....." << endl;   
        ifstream codeStream(cb_in.c_str());
        if (!codeStream) {
          cerr << argv[0] << ": can't open file " << cb_in << endl;
          exit(EXIT_FAILURE);
        }
	myMap = new xmippFuzzyMap(codeStream, ts.size(), true);
   } else    
        myMap = new xmippFuzzyMap(layout, xdim, ydim, ts, use_rand_cvs);

   
   xmippFuzzySOM *thisSOM;
   if (fn_algo_in == "") {
   	thisSOM = new xmippFuzzySOM(m0, m1, annSteps, reg, eps, iter);    // Creates FSOM Algorithm
   } else {
        cout << "Reading algorithm file " << fn_algo_in << "....." << endl << endl;   
        ifstream algoStream(fn_algo_in.c_str());
        if (!algoStream) {
          cerr << argv[0] << ": can't open file " << fn_algo_in << endl;
          exit(EXIT_FAILURE);
        }
   }

   xmippTextualListener myListener;	    // Define the listener class
   myListener.setVerbosity() = verb;	    // Set verbosity level
   thisSOM->setListener(&myListener);       // Set Listener

   if (cb_in != "") {
   	if (ts.isNormalized()) {
   		cout << "Normalizing code vectors....." << endl;
   		myMap->Normalize(ts.getNormalizationInfo()); 	     // normalize code vectors
   	}	
   	thisSOM->train(*myMap, ts);            		     // Train algorithm	
   } else
   	thisSOM->train(*myMap, ts);        		     	     // Train algorithm

   // Test algorithm
   double dist = thisSOM->test(*myMap, ts);
   cout << endl << "Quantization error : " <<  dist << endl;

   // Calculates functional value
   double functional, fidelity, penalty;
   functional = thisSOM->functional(ts, *myMap, m1, reg, fidelity, penalty);
   cout << "Functional : " <<  functional << " (fidelity = " << fidelity << " penalty = " << penalty << " )" << endl << endl;


   // Classifying
   cout << "Classifying....." << endl;
   myMap->classify(&ts);

   // Calibrating
   cout << "Calibrating....." << endl;
   myMap->calibrate(ts);

  /******************************************************* 
      Saving all kind of Information 
  *******************************************************/

   cout << "Saving algorithm information as " << fn_out << ".inf ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".inf"; 
   ofstream infS(tmpN.c_str());
   infS << "Fuzzy SOM algorithm" << endl << endl;
   infS << "Input data file : " << fn_in << endl;
   if (cb_in != "")
     infS << "Input code vectors file : " << cb_in << endl;
   else if (use_rand_cvs)
     infS << "Using randomized code vectors" << endl;
   infS << "Code vectors output file : " << fn_out <<  ".cod" << endl;
   infS << "Algorithm information output file : " << fn_out <<  ".inf" << endl;
   infS << "Number of feature vectors: " << ts.size() << endl;
   infS << "Number of variables: " << ts.theItems[0].size() << endl;
   infS << "Horizontal dimension (Xdim) = " << xdim << endl;
   infS << "Vertical dimension (Ydim) = " << ydim << endl;
   if (layout == "HEXA")
     infS << "Hexagonal topology " << endl;
   else
     infS << "Rectangular topology " << endl;
   if (norm)
     infS << "Input data normalized" << endl;
   else
     infS << "Input data not normalized" << endl;
   infS << "Initial Fuzzy constant (m0) = " << m0 << endl;
   infS << "Final Fuzzy constant (m1) = " << m1 << endl;
   infS << "Smoothness factor (reg) = " << reg << endl;
   infS << "Deterministic annealing steps = " << annSteps << endl;
   infS << "Total number of iterations = " << iter << endl;
   infS << "Stopping criteria (eps) = " << eps << endl;
   infS << "Quantization error : " <<  dist << endl;
   infS << "Functional : " <<  functional << " (fidelity = " << fidelity << " penalty = " << penalty << " )" << endl << endl;
   infS.flush();    

   // assign data to clusters according to fuzzy threshold
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
   vsStream << ts.theItems[0].size() << " " << myMap->layout() << " " << myMap->width() << " " << myMap->height() << " gaussian" << endl;
   for (int i= 0; i < ts.size(); i++) {
   	int j = myMap->fuzzyWinner(i);
   	vsStream << myMap->indexToPos(j).first << " " << myMap->indexToPos(j).second << " " << myMap->memb[i][j] << " " << ts.theTargets[i] << endl;
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

   if (norm) {
   	cout << "Denormalizing code vectors....." << endl;
   	myMap->unNormalize(ts.getNormalizationInfo()); // de-normalize codevectors        
   }	

   cout << "Saving code vectors as " << fn_out << ".cod ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".cod"; 
   ofstream codS(tmpN.c_str());
   codS << *myMap;
   codS.flush();    


   cout << endl;
   
   delete myMap;
   delete thisSOM;

   
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
     "\nPurpose: Fuzzy Self-Organizing Feature Map"
     "\n"           
     "\nParameter Values: (note space before value)"
     "\n"
     "\n    -din    file_in           Input data file"
     "\n    -cout   file_out          Base name for output data files "
     "\n    -cvin   file_in           Codevectors input file"
     "\n    -saveclusters    	      save clusters in separate files (Default = No)"     
     "\n    -xdim   H-dimension	      Horizontal size of the map"
     "\n    -ydim   V-dimension       Vertical size of the map"
     "\n    -hexa    	   	      Hexagonal topology"     
     "\n    -rect    	   	      Rectangular topology (default)"     
     "\n    -steps  steps    	      Deterministic annealing steps (default = 1000)"     
     "\n    -m0     Initial m         Initial Fuzzy constant (default = 2)"
     "\n    -m1     Final m 	      Final Fuzzy constant (default = 1.02)"     
     "\n    -reg    smoothness        Smoothness factor (default = 0.5)"     
     "\n    -eps    Epsilon 	      Stopping criteria (default = 1e-7)"
     "\n    -iter   iterations        Number of iterations (default = 1000)"     
     "\n    -norm            	      Normalize training data (default)"     
     "\n    -verb   verbosity         Information level while running: "     
     "\n    			      0: No information (default)"     
     "\n    			      1: Progress bar"     
     "\n    			      2: Code vectors change between iterations"     
     "\n    -randomcodevectors        Use truly randomized codevectors (default: FALSE)"
     "\n			   \n"
     ,argv[0]);
}
