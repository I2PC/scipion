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

#include <classification/fcmeans.h>

/* Prototypes -============================================================= */

void Usage(char **argv);

/* Main function -============================================================= */

main(int argc, char** argv)
{


    /* Input Parameters ======================================================== */

    FileName       fn_in;     // input file
    FileName       fn_out;    // output file
    FileName       cb_in = "";    // Code vectors input file
    FileName       tmpN;  // Temporary variable
    unsigned       c;  // Number of clusters
    xmippFeature   m;  // Fuzzy membership
    xmippFeature   eps = 1e-7; // Stopping criteria
    unsigned       iter = 1000; // Iteration number
    unsigned       verb = 0; // Verbosity level
    bool           norm = 1; // Normalize?
    bool         saveClusters = false;    // Save clusters in separate files

    /* Parameters ============================================================== */
    try
    {

        if (checkParameter(argc, argv, "-i"))
            fn_in = getParameter(argc, argv, "-i");
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-o"))
            fn_out = getParameter(argc, argv, "-o");
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-cvin"))
            cb_in = getParameter(argc, argv, "-cvin");

        if (checkParameter(argc, argv, "-c"))
            c = AtoI(getParameter(argc, argv, "-c"));
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-saveclusters"))
            saveClusters = true;
        else saveClusters = false;

        m = AtoF(getParameter(argc, argv, "-m", "2.0"));
        eps = AtoF(getParameter(argc, argv, "-eps", "1e-7"));
        iter = AtoI(getParameter(argc, argv, "-iter", "1000"));
        verb = AtoI(getParameter(argc, argv, "-verb", "0"));

        if (checkParameter(argc, argv, "-norm"))
            norm = true;
        else
            norm = false;

        if (argc == 1)
        {
            Usage(argv);
        }

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage(argv);
    }

    /* Some validations ===================================================== */

    if (c <= 1)
    {
        cerr << "Number of clusters c = " << c << " must be > 1" << endl;
        exit(EXIT_FAILURE);
    }

    if (m <= 1)
    {
        cerr << "Fuzzy constant m = " << m << " must be > 1" << endl;
        exit(EXIT_FAILURE);
    }

    if (iter < 1)
    {
        cerr << argv[0] << ": invalid value for iter (must be > 1): " << iter << endl;
        exit(EXIT_FAILURE);
    }

    if (verb < 0 || verb > 2)
    {
        cerr << argv[0] << ": invalid value for verbosity (must be between 0 and 2): " << verb << endl;
        exit(EXIT_FAILURE);
    }


    /* Shows parameters ===================================================== */

    cout << endl << "Parameters used: " << endl;
    cout << "Input data file : " << fn_in << endl;
    cout << "Output file name : " << fn_out << endl;
    if (cb_in != "")
        cout << "Input cluster centers file name : " << cb_in << endl;
    cout << "Number of clusters = " << c << endl;
    cout << "Fuzzy constant = " << m << endl;
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
    if (!inStream)
    {
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

    if (c >= ts.size())
    {
        cerr << "Number of clusters c = " << c << " must be < " << ts.size() << endl;
        exit(EXIT_FAILURE);
    }


    /* Real stuff ============================================================== */


    try
    {
        if (norm)
        {
            cout << "Normalizing....." << endl;
            ts.normalize();        // Normalize input data
        }

        xmippFCMeans thisFCmeans(m, eps, iter);  // Creates Fuzzy c-means Algorithm

        xmippFCB* thisFCB;
        if (cb_in != "")
        {
            cout << "Reading fuzzy cluster centers file " << cb_in << "....." << endl;
            ifstream codeStream(cb_in.c_str());
            if (!codeStream)
            {
                cerr << argv[0] << ": can't open file " << cb_in << endl;
                exit(EXIT_FAILURE);
            }
            thisFCB = new xmippFCB(codeStream, ts.size()); // Reads FuzzyCodeBook from file
        }
        else
            thisFCB = new xmippFCB(c, ts); // initialize Fuzzy codebook randomly

        xmippTextualListener myListener;     // Define the listener class
        myListener.setVerbosity() = verb;     // Set verbosity level
        thisFCmeans.setListener(&myListener);    // Set Listener

        thisFCmeans.train(*thisFCB, ts);          // Train algorithm

        // Test algorithm
        xmippFeature qerror = thisFCmeans.test(*thisFCB, ts);
        cout << "Quantization error : " <<  qerror << endl;

        // Classifying
        cout << "Classifying....." << endl;
        thisFCB->classify(&ts);

        // Calibrating
        cout << "Calibrating....." << endl;
        thisFCB->fuzzyCalibrate(ts);

        // Shows Validity functionals (If applicable)

        xmippFeature F = thisFCmeans.F(*thisFCB);
        xmippFeature H = thisFCmeans.H(*thisFCB);
        xmippFeature NFI = thisFCmeans.NFI(*thisFCB);
        xmippFeature S = thisFCmeans.S(*thisFCB, ts);
        cout << endl << "Validity Functionals : " << endl;
        cout << "Partition coefficient (max) (F) : " << F << endl;
        cout << "Partition entropy (min) (H) : " << H << endl;
        cout << "Non-Fuzzy Index (max) (NFI) : " << NFI << endl;
        cout << "Compactness and Separation index (min) (S) : " << S << endl;
        cout << endl;

        // assign data to clusters according to fuzzy threshold
        if (saveClusters)
        {
            cout << "Saving clusters assigments ....." << endl;
            for (unsigned i = 0; i < thisFCB->size(); i++)
            {
                tmpN = fn_out.c_str() + (string) "."  + ItoA(i);
                ofstream cStream(tmpN.c_str());
                for (int j = 0; j < thisFCB->classifAt(i).size(); j++)
                    cStream << thisFCB->classifAt(i)[j] << endl;
                cStream.flush();
            }
        }

        // save .his file (Histogram)
        cout << "Saving clusters histogram file as " << fn_out << ".his ....." << endl;
        tmpN = fn_out.c_str() + (string) ".his";
        ofstream hisStream(tmpN.c_str());
        thisFCB->printHistogram(hisStream);
        hisStream.flush();

        // save .err file (Average Quantization Error)
        cout << "Saving cluster centers average quantization error file as " << fn_out << ".err ....." << endl;
        tmpN = fn_out.c_str() + (string) ".err";
        ofstream errStream(tmpN.c_str());
        thisFCB->printQuantError(errStream);
        errStream.flush();

        // save .vs file to be compatible with SOM_PAK
        cout << "Saving visual file as " << fn_out << ".vs ....." << endl;
        tmpN = fn_out.c_str() + (string) ".vs";
        ofstream vsStream(tmpN.c_str());
        vsStream << ts.theItems[0].size() << " " << "FCMeans" << " " << thisFCB->membClusters() << " 1" << " gaussian" << endl;
        for (int i = 0; i < ts.size(); i++)
        {
            int j = thisFCB->fuzzyWinner(i);
            vsStream << j << " 0 " << thisFCB->memb[i][j] << " " << ts.theTargets[i] << endl;
        }
        vsStream.flush();

        cout << "Saving algorithm information as " << fn_out << ".inf ....." << endl;
        tmpN = fn_out.c_str() + (string) ".inf";
        ofstream infS(tmpN.c_str());
        infS << "Fuzzy c-Means Clustering Algorithm (FCMeans)" << endl << endl;
        infS << "Input data file : " << fn_in << endl;
        if (cb_in != "")
            infS << "Input cluster centers file : " << cb_in << endl;
        infS << "Cluster centers output file : " << fn_out <<  ".cod" << endl;
        if (cb_in != "")
            infS << "Input cluster centers file name : " << cb_in << endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << endl;
        infS << "Number of feature vectors: " << ts.size() << endl;
        infS << "Number of variables: " << ts.theItems[0].size() << endl;
        infS << "Number of clusters: " << thisFCB->membClusters() << endl;
        if (norm)
            infS << "Input data normalized" << endl;
        else
            infS << "Input data not normalized" << endl;
        infS << "Fuzzy constant (m) = " << m << endl;
        infS << "Total number of iterations = " << iter << endl;
        infS << "Stopping criteria (eps) = " << eps << endl;
        infS << "Quantization error : " <<  qerror << endl;
        infS << endl << "Validity Functionals : " << endl;
        infS << "Partition coefficient (max) (F) : " << F << endl;
        infS << "Partition entropy (min) (H) : " << H << endl;
        infS << "Non-Fuzzy Index (max) (NFI) : " << NFI << endl;
        infS << "Compactness and Separation index (min) (S) : " << S << endl;
        infS.flush();

        if (norm)
        {
            cout << "Denormalizing cluster centers....." << endl;
            thisFCB->unNormalize(ts.getNormalizationInfo()); // de-normalize cluster centers
        }

        // Save codevectors
        cout << "Saving cluster centers in " << fn_out << ".cod....." << endl;
        tmpN = fn_out.c_str() + (string) ".cod";
        ofstream outStream(tmpN.c_str());
        outStream << ts.theItems[0].size() << " " << thisFCB->size() << endl;
        outStream << *thisFCB;
        outStream.flush();

        delete thisFCB;

    }
    catch (const exception& e)
    {
        cout << e.what() << endl;
    }
    return 0;
}


/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage(char **argv)
{
    printf(
        "\nUsage: %s [Purpose and Parameters]"
        "\nPurpose: Fuzzy partition (clustering) using Fuzzy-cmeans algorithm"
        "\n"
        "\nParameter Values: (note space before value)"
        "\n"
        "\n    -i    file_in        Input data file"
        "\n    -o    file_out       Base name for output data files "
        "\n    -cvin   file_in      Cluster centers input file"
        "\n    -saveclusters      Save clusters in separate files (Default = No)"
        "\n    -c    Clusters   Number of clusters"
        "\n    -m    Fuzzy constant Fuzzy constant (default = 2.0)"
        "\n    -eps  Epsilon   Stopping criteria (default = 1e-7)"
        "\n    -iter iterations    Number of iterations (default = 1000)"
        "\n    -norm       Normalize training data (default: No)"
        "\n    -verb verbosity    Information level while running: "
        "\n        0: No information (default)"
        "\n        1: Progress bar"
        "\n        2: Changes between iterations"
        "\n    \n"
        , argv[0]);
}
