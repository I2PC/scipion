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
#include <Classification/xmippRM.hh>

/* Prototypes -============================================================= */

void Usage (char **argv);

/* Main function -============================================================= */

main(int argc, char** argv) {


/* Input Parameters ======================================================== */

FileName       fn_in;    	// input file
FileName       fn_out;   	// output file
FileName       cb_in = "";   	// Code vectors input file
FileName       tmpN;		// Temporary variable
unsigned       verb = 0;	// Verbosity level
unsigned       k;               // Dimension of projected subspace  
string         matdist = "sparse"; // Type of random matrix distribution

/* Parameters ============================================================== */
   try {

       fn_in = get_param(argc, argv, "-din");

       if (check_param(argc, argv, "-cout"))
          fn_out = get_param(argc, argv, "-cout");
       else {
         Usage(argv); 
	 exit(EXIT_FAILURE);
       } 

       k = AtoI(get_param(argc, argv, "-k"));

       if (check_param(argc, argv, "-sparse")) {
       	  if (check_param(argc, argv, "-gaussian")) {
	        cout << "Error: you can not define two types of distributions" << endl;
	 	exit(EXIT_FAILURE);	  	
	  }
	  matdist = "sparse";
       } else if (check_param(argc, argv, "-gaussian"))
          matdist = "gaussian";

       verb = AtoI(get_param(argc, argv, "-verb", "1"));

       if (argc == 1) {Usage(argv);}
       
   }
   catch (Xmipp_error XE) {cout << XE; Usage(argv);}


/* Some validations ===================================================== */
  
   if (k < 1) {
     cerr << argv[0] << ": invalid dimension (-k) (must be > 0): " << k << endl;
     exit(EXIT_FAILURE);
   }

   if (verb < 0 || verb > 1) {
     cerr << argv[0] << ": invalid value for verbosity (must be between 0 and 1): " << verb << endl;
     exit(EXIT_FAILURE);
   }


/* Shows parameters ===================================================== */

   cout << endl << "Parameters used: " << endl;
   cout << "Input data file : " << fn_in << endl;
   cout << "Output file name : " << fn_out << endl;
   cout << "Dimension of projected subspace : " << k << endl;
   if (matdist == "sparse")
     cout << "Sparsed distributed random matrix" << endl;
   else
     cout << "Gaussian distributed random matrix" << endl;
   cout << "verbosity level = " << verb << endl;


/* Open training vector ================================================= */

  
  ifstream inStream(fn_in.c_str());
  if (!inStream) {
      cerr << argv[0] << ": can't open file " << fn_in << endl;
      exit(EXIT_FAILURE);
  }

  xmippCTVectors ts(0, false);

  cout << endl << "Reading data file..." << fn_in << "....." << endl << endl;     
  inStream >> ts;

  xmippCTVectors projectedTs(0, true);


/* Real stuff ============================================================== */
 try {
 
   xmippRM RM;

   xmippTextualListener myListener;	    // Define the listener class
   myListener.setVerbosity() = verb;	    // Set verbosity level
   RM.setListener(&myListener);       	    // Set Listener

   RM.Project(ts, projectedTs, k);

  /******************************************************* 
      Saving all kind of Information 
  *******************************************************/

   cout << "Saving projected vectors as " << fn_out << ".dat ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".dat"; 
   ofstream projS(tmpN.c_str());
   projS << k << " " << projectedTs.size() << endl;
   for (int j = 0; j < projectedTs.size() ; j++) {
     for (int i = 0; i < k; i++)
       projS << " " << projectedTs.theItems[j][i];
     projS << " " << ts.theTargets[j] << endl;  
   }
   projS.flush();    

   cout << "Saving algorithm information as " << fn_out << ".inf ....." << endl;  
   tmpN = fn_out.c_str() + (string) ".inf"; 
   ofstream infS(tmpN.c_str());
   infS << "Random Mapping" << endl << endl;
   infS << "Input data file : " << fn_in << endl;
   infS << "Dimension of projected subspace : " << k << endl;
   if (matdist == "sparse")
     infS << "Sparsed distributed random matrix" << endl;
   else
     infS << "Gaussian distributed random matrix" << endl;
   infS << "Projected vectors output file : " << fn_out <<  ".dat" << endl;
   infS << "Algorithm information output file : " << fn_out <<  ".inf" << endl;
   infS << "Number of feature vectors: " << ts.size() << endl;
   infS << "Number of variables: " << ts.itemAt(0).size() << endl;
   infS.flush();    

   cout << endl;
   
 } catch ( Xmipp_error &e ) {
    cout << e << endl;
 }
  return 0;
}


/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage (char **argv) {
  printf (
     "\nUsage: %s [Purpose and Parameters]"
     "\nPurpose: Random Mapping (RM)"
     "\n"           
     "\nParameter Values: (note space before value)"
     "\n"
     "\n    -din    file_in           Input data file (plain data)"
     "\n    -cout   file_out          Base name for output data files:"
     "\n    -k      k                 Dimension of projected subspace" 
     "\n    -sparse                   Sparse distribution of random matrix (default)" 
     "\n    -gaussian                 Gaussian distribution of random matrix" 
     "\n    -verb      verbosity      Information level while running: "     
     "\n    			      0: No information"     
     "\n    			      1: Progress bar (default)"     
     "\n			   \n"
     ,argv[0]);
     exit(0);
}
