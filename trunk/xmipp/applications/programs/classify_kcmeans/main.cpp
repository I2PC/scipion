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

// To avoid problems with long template names
#pragma warning(disable:4786)

#include <fstream>

#include <classification/tstudent_kerdensom.h>
#include <classification/gaussian_kerdensom.h>

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
unsigned       iter = 200;	// Iteration number
unsigned       verb = 0;	// Verbosity level
bool           norm = 1;	// Normalize?
unsigned       c;		// Number of clusters
int	       df = 3;		// Degrees of freedom	
bool           useCBook = false;   	// Use codebook
bool 	       saveClusters = false;    // Save clusters in separate files
bool 	       saveCodebook = false;    // Save codebook in a separate file
bool 	       gaussian = true;         // Gaussian Kernel
bool 	       tStudent = false;        // tStudent Kernel

/* Parameters ============================================================== */

   try {

       fn_in = get_param(argc, argv, "-din");

       if (check_param(argc, argv, "-cout"))
          fn_out = get_param(argc, argv, "-cout");
       else {
         Usage(argv);
	 exit(EXIT_FAILURE);
       }

       if (check_param(argc, argv, "-cvin")) {
       	  if (check_param(argc, argv, "-cbin")) {
	        cout << "Error: you can not use two code vectors files" << endl;
	 	exit(EXIT_FAILURE);	  	
	  }
       	  cb_in = get_param(argc, argv, "-cvin");
	  useCBook = false;
       }


       if (check_param(argc, argv, "-cbin")) {
       	  if (check_param(argc, argv, "-cvin")) {
	        cout << "Error: you can not use two code vectors files" << endl;
	 	exit(EXIT_FAILURE);	  	
	  }
       	  cb_in = get_param(argc, argv, "-cbin");
	  useCBook = true;
       }


       c = AtoI(get_param(argc, argv, "-c"));

       if (check_param(argc, argv, "-gaussian")) {
       	  if (check_param(argc, argv, "-tStudent")) {
	        cout << "Error: you can not define two kernels functions" << endl;
	 	exit(EXIT_FAILURE);	  	
	  }
	  gaussian = true;
          tStudent = false;
       } else if (check_param(argc, argv, "-tStudent")) {
          gaussian = false;
	  tStudent = true;
       }

       df = (int) AtoI(get_param(argc, argv, "-df", "3"));

       eps = AtoF(get_param(argc, argv, "-eps", "1e-7"));
       iter = AtoI(get_param(argc, argv, "-iter", "200"));
       verb = AtoI(get_param(argc, argv, "-verb", "0"));

       if (check_param(argc, argv, "-norm"))
          norm = true;
       else norm = false;

       if (check_param(argc, argv, "-saveclusters"))
          saveClusters = true;
       else saveClusters = false;

       if (check_param(argc, argv, "-savecodebook"))
          saveCodebook = true;
       else saveCodebook = false;

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

   if (c < 2) {
     cerr << argv[0] << ": invalid value for c (must be >= 2): " << c << endl;
     exit(EXIT_FAILURE);
   }

   if (df < 2) {
     cerr << argv[0] << ": invalid value for df (must be > 1): " << df << endl;
     exit(EXIT_FAILURE);
   }



/* Shows parameters ===================================================== */

   cout << endl << "Parameters used: " << endl;
   cout << "Input data file : " << fn_in << endl;
   cout << "Output file name : " << fn_out << endl;
   if (cb_in != "")
   	cout << "Input code vectors file name : " << cb_in << endl;
   if (saveClusters)
     cout << "Save clusters in separate files: " << fn_out << ".(cluster number)"<< endl;
   cout << "Number of clusters (c) = " << c << endl;
   if (gaussian)
     cout << "Gaussian Kernel function " << endl;
   else {
     cout << "t-Student Kernel function" << endl;
     cout << "Degrees of freedom (df) = " << df << endl;
   }
   cout << "Total number of iterations = " << iter << endl;
   cout << "Stopping criteria (eps) = " << eps << endl;
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
  cout << endl << "Reading input data file " << fn_in << "....." << endl;
  inStream >> ts;



/* Real stuff ============================================================== */


 try {

   if (norm) {
   	cout << "Normalizing....." << endl;
   	ts.normalize();  		    // Normalize input data
   }	

   xmippFuzzyMap *myMap;

   if (cb_in != "") {
     if (useCBook) {
        cout << "Reading fuzzy codebook file " << cb_in << "....." << endl;
        ifstream codeStream(cb_in.c_str());
        if (!codeStream) {
          cerr << argv[0] << ": can't open file " << cb_in << endl;
          exit(EXIT_FAILURE);
        }
        myMap = new xmippFuzzyMap(codeStream, ts.size(), false);
     } else {
        cout << "Reading fuzzy codevectors file " << cb_in << "....." << endl;
        ifstream codeStream(cb_in.c_str());
        if (!codeStream) {
          cerr << argv[0] << ": can't open file " << cb_in << endl;
          exit(EXIT_FAILURE);
        }
        myMap = new xmippFuzzyMap(codeStream, ts.size(), true);
     }
   } else
        myMap = new xmippFuzzyMap("RECT", c, 1, ts);


   xmippKerDenSOM *thisSOM;
   if (fn_algo_in == "") {
        if (gaussian)
           thisSOM = new xmippGaussianKerDenSOM(0, 0, 0, eps, iter);        // Creates KerDenSOM Algorithm
        else 		
   	   thisSOM = new xmippTStudentKerDenSOM(0, 0, 0, eps, iter, df);    // Creates KerDenSOM Algorithm
   } else {
        cout << "Reading algorithm file " << fn_algo_in << "....." << endl << endl;
        ifstream algoStream(fn_algo_in.c_str());
        if (!algoStream) {
          cerr << argv[0] << ": can't open file " << fn_algo_in << endl;
          exit(EXIT_FAILURE);
        }
   }

   xmippTextualListener myListener;	      // Define the listener class
   myListener.setVerbosity() = verb;	      // Set verbosity level
   thisSOM->setListener(&myListener);         // Set Listener

   if (cb_in != "") {
   	if (ts.isNormalized()) {
   		cout << "Normalizing code vectors....." << endl;
   		myMap->Normalize(ts.getNormalizationInfo()); 	     // normalize code vectors
   	}	
   	thisSOM->train(*myMap, ts, fn_out, true);            // Train algorithm
   } else
   	thisSOM->train(*myMap, ts, fn_out);        		     // Train algorithm
	

   // Test algorithm
   xmippFeature dist = thisSOM->test(*myMap, ts);
   cout << endl << "Quantization error : " <<  dist << endl;

   // Classifying
   cout << "Classifying....." << endl;
   myMap->classify(&ts);

   // Calibrating
   cout << "Calibrating....." << endl;
   myMap->calibrate(ts);

  /*******************************************************
      Saving all kind of Information
  *******************************************************/

   if (saveCodebook) {
   	cout << "Saving whole codebook as " << fn_out << ".cbk ....." << endl;
   	tmpN = fn_out.c_str() + (string) ".cbk";
   	ofstream cbkS(tmpN.c_str());
   	myMap->saveObject(cbkS);
   	cbkS.flush();
   }

   cout << "Saving algorithm information as " << fn_out << ".inf ....." << endl;
   tmpN = fn_out.c_str() + (string) ".inf";
   ofstream infS(tmpN.c_str());
   infS << "Kernel Probability Density Estimator clustering algorithm" << endl;
   infS << "                 Kernel c-Means" << endl << endl;
   infS << "Input data file : " << fn_in << endl;
   if (cb_in != "")
     infS << "Input code vectors file : " << cb_in << endl;
   infS << "Code vectors output file : " << fn_out <<  ".cod" << endl;
   infS << "Whole codebook output file : " << fn_out <<  ".cbk" << endl;
   infS << "Algorithm information output file : " << fn_out <<  ".inf" << endl;
   infS << "Number of feature vectors: " << ts.size() << endl;
   infS << "Number of variables: " << ts.theItems[0].size() << endl;
   infS << "Number of clusters = " << c << endl;
   if (norm)
     infS << "Input data normalized" << endl;
   else
     infS << "Input data not normalized" << endl;

   if (gaussian)
     infS << "Gaussian Kernel function " << endl;
   else {
     infS << "t-Student Kernel function" << endl;
     infS << "Degrees of freedom (df) = " << df << endl;
   }

   infS << "Total number of iterations = " << iter << endl;
   infS << "Stopping criteria (eps) = " << eps << endl;
   infS << "Final Sigma = " << thisSOM->getSigma() << endl;
   infS << "Quantization error : " <<  dist << endl;
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
     "\nPurpose: Kernel Density Estimator Clustering algorithm"
     "\nParameter Values: (note space before value)"
     "\n    -din    file_in           Input data file (plain data)"
     "\n    -cout   file_out          Base name for output data files"
     "\n    -cvin   file_in           Codevectors input file"
     "\n    -saveclusters    	      save clusters in separate files (Default = No)"
     "\n    -c      Clusters	      Number of clusters"
     "\n    -gaussian    	      Gaussian Kernel Function (default)"
     "\n    -tStudent    	      t-Student Kernel Function "
     "\n    -df     df 	      	      t-Student degrees of freedom (default = 3)"
     "\n    -eps    Epsilon 	      Stopping criteria (default = 1e-7)"
     "\n    -iter   iterations        Number of iterations (default = 200)"
     "\n    -norm            	      Normalize training data (default = No)"
     "\n    -verb   verbosity         Information level while running: "
     "\n    			      0: No information (default)"
     "\n    			      1: Progress bar"
     "\n    			      2: Changes between iterations"
     "\n"
     ,argv[0]);
     exit(0);
}
