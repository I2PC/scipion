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
            c = textToInteger(getParameter(argc, argv, "-c"));
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-saveclusters"))
            saveClusters = true;
        else saveClusters = false;

        m = textToFloat(getParameter(argc, argv, "-m", "2.0"));
        eps = textToFloat(getParameter(argc, argv, "-eps", "1e-7"));
        iter = textToInteger(getParameter(argc, argv, "-iter", "1000"));
        verb = textToInteger(getParameter(argc, argv, "-verb", "0"));

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
        std::cout << XE;
        Usage(argv);
    }

    /* Some validations ===================================================== */

    if (c <= 1)
    {
        std::cerr << "Number of clusters c = " << c << " must be > 1" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (m <= 1)
    {
        std::cerr << "Fuzzy constant m = " << m << " must be > 1" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (iter < 1)
    {
        std::cerr << argv[0] << ": invalid value for iter (must be > 1): " << iter << std::endl;
        exit(EXIT_FAILURE);
    }

    if (verb < 0 || verb > 2)
    {
        std::cerr << argv[0] << ": invalid value for verbosity (must be between 0 and 2): " << verb << std::endl;
        exit(EXIT_FAILURE);
    }


    /* Shows parameters ===================================================== */

    std::cout << std::endl << "Parameters used: " << std::endl;
    std::cout << "Input data file : " << fn_in << std::endl;
    std::cout << "Output file name : " << fn_out << std::endl;
    if (cb_in != "")
        std::cout << "Input cluster centers file name : " << cb_in << std::endl;
    std::cout << "Number of clusters = " << c << std::endl;
    std::cout << "Fuzzy constant = " << m << std::endl;
    std::cout << "Total number of iterations = " << iter << std::endl;
    std::cout << "Stopping criteria (eps) = " << eps << std::endl;
    std::cout << "verbosity level = " << verb << std::endl;
    if (norm)
        std::cout << "Normalize input data" << std::endl;
    else
        std::cout << "Do not normalize input data " << std::endl;


    /* Open training vector ================================================= */

    std::cout << std::endl << "Reading file " << fn_in << "....." << std::endl;

    std::ifstream inStream(fn_in.c_str());
    if (!inStream)
    {
        std::cerr << argv[0] << ": can't open file " << fn_in << std::endl;
        exit(EXIT_FAILURE);
    }

    xmippCTVectors ts(0, true);
    try
    {
        inStream >> ts;
    }
    catch (std::exception& e)
    {
        std::cerr << argv[0] << ": can't read file " << fn_in  << " because " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    if (c >= ts.size())
    {
        std::cerr << "Number of clusters c = " << c << " must be < " << ts.size() << std::endl;
        exit(EXIT_FAILURE);
    }


    /* Real stuff ============================================================== */


    try
    {
        if (norm)
        {
            std::cout << "Normalizing....." << std::endl;
            ts.normalize();        // Normalize input data
        }

        xmippFCMeans thisFCmeans(m, eps, iter);  // Creates Fuzzy c-means Algorithm

        xmippFCB* thisFCB;
        if (cb_in != "")
        {
            std::cout << "Reading fuzzy cluster centers file " << cb_in << "....." << std::endl;
            std::ifstream codeStream(cb_in.c_str());
            if (!codeStream)
            {
                std::cerr << argv[0] << ": can't open file " << cb_in << std::endl;
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
        std::cout << "Quantization error : " <<  qerror << std::endl;

        // Classifying
        std::cout << "Classifying....." << std::endl;
        thisFCB->classify(&ts);

        // Calibrating
        std::cout << "Calibrating....." << std::endl;
        thisFCB->fuzzyCalibrate(ts);

        // Shows Validity functionals (If applicable)

        xmippFeature F = thisFCmeans.F(*thisFCB);
        xmippFeature H = thisFCmeans.H(*thisFCB);
        xmippFeature NFI = thisFCmeans.NFI(*thisFCB);
        xmippFeature S = thisFCmeans.S(*thisFCB, ts);
        std::cout << std::endl << "Validity Functionals : " << std::endl;
        std::cout << "Partition coefficient (max) (F) : " << F << std::endl;
        std::cout << "Partition entropy (min) (H) : " << H << std::endl;
        std::cout << "Non-Fuzzy Index (max) (NFI) : " << NFI << std::endl;
        std::cout << "Compactness and Separation index (min) (S) : " << S << std::endl;
        std::cout << std::endl;

        // assign data to clusters according to fuzzy threshold
        if (saveClusters)
        {
            std::cout << "Saving clusters assigments ....." << std::endl;
            for (unsigned i = 0; i < thisFCB->size(); i++)
            {
                tmpN = fn_out.c_str() + (std::string) "."  + integerToString(i);
                std::ofstream cStream(tmpN.c_str());
                for (int j = 0; j < thisFCB->classifAt(i).size(); j++)
                    cStream << thisFCB->classifAt(i)[j] << std::endl;
                cStream.flush();
            }
        }

        // save .his file (Histogram)
        std::cout << "Saving clusters histogram file as " << fn_out << ".his ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".his";
        std::ofstream hisStream(tmpN.c_str());
        thisFCB->printHistogram(hisStream);
        hisStream.flush();

        // save .err file (Average Quantization Error)
        std::cout << "Saving cluster centers average quantization error file as " << fn_out << ".err ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".err";
        std::ofstream errStream(tmpN.c_str());
        thisFCB->printQuantError(errStream);
        errStream.flush();

        // save .vs file to be compatible with SOM_PAK
        std::cout << "Saving visual file as " << fn_out << ".vs ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".vs";
        std::ofstream vsStream(tmpN.c_str());
        vsStream << ts.theItems[0].size() << " " << "FCMeans" << " " << thisFCB->membClusters() << " 1" << " gaussian" << std::endl;
        for (int i = 0; i < ts.size(); i++)
        {
            int j = thisFCB->fuzzyWinner(i);
            vsStream << j << " 0 " << thisFCB->memb[i][j] << " " << ts.theTargets[i] << std::endl;
        }
        vsStream.flush();

        std::cout << "Saving algorithm information as " << fn_out << ".inf ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".inf";
        std::ofstream infS(tmpN.c_str());
        infS << "Fuzzy c-Means Clustering Algorithm (FCMeans)" << std::endl << std::endl;
        infS << "Input data file : " << fn_in << std::endl;
        if (cb_in != "")
            infS << "Input cluster centers file : " << cb_in << std::endl;
        infS << "Cluster centers output file : " << fn_out <<  ".cod" << std::endl;
        if (cb_in != "")
            infS << "Input cluster centers file name : " << cb_in << std::endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
        infS << "Number of feature vectors: " << ts.size() << std::endl;
        infS << "Number of variables: " << ts.theItems[0].size() << std::endl;
        infS << "Number of clusters: " << thisFCB->membClusters() << std::endl;
        if (norm)
            infS << "Input data normalized" << std::endl;
        else
            infS << "Input data not normalized" << std::endl;
        infS << "Fuzzy constant (m) = " << m << std::endl;
        infS << "Total number of iterations = " << iter << std::endl;
        infS << "Stopping criteria (eps) = " << eps << std::endl;
        infS << "Quantization error : " <<  qerror << std::endl;
        infS << std::endl << "Validity Functionals : " << std::endl;
        infS << "Partition coefficient (max) (F) : " << F << std::endl;
        infS << "Partition entropy (min) (H) : " << H << std::endl;
        infS << "Non-Fuzzy Index (max) (NFI) : " << NFI << std::endl;
        infS << "Compactness and Separation index (min) (S) : " << S << std::endl;
        infS.flush();

        if (norm)
        {
            std::cout << "Denormalizing cluster centers....." << std::endl;
            thisFCB->unNormalize(ts.getNormalizationInfo()); // de-normalize cluster centers
        }

        // Save codevectors
        std::cout << "Saving cluster centers in " << fn_out << ".cod....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".cod";
        std::ofstream outStream(tmpN.c_str());
        outStream << ts.theItems[0].size() << " " << thisFCB->size() << std::endl;
        outStream << *thisFCB;
        outStream.flush();

        delete thisFCB;

    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
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
