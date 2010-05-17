/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

// To avoid problems with long template names
#pragma warning(disable:4786)

#include <fstream>

#include <classification/tstudent_kerdensom.h>
#include <classification/gaussian_kerdensom.h>

/* Prototypes -============================================================= */

void Usage(char **argv);

/* Main function -============================================================= */

main(int argc, char** argv)
{


    /* Input Parameters ======================================================== */

    FileName       fn_in;     // input file
    FileName       fn_out;    // output file
    FileName       cb_in = "";    // Code vectors input file
    FileName       fn_algo_in = ""; // input algorithm file
    FileName       tmpN;  // Temporary variable
    double         eps = 1e-7; // Stopping criteria
    unsigned       iter = 200; // Iteration number
    unsigned       verb = 0; // Verbosity level
    bool           norm = 1; // Normalize?
    unsigned       c;  // Number of clusters
    int        df = 3;  // Degrees of freedom
    bool           useCBook = false;    // Use codebook
    bool         saveClusters = false;    // Save clusters in separate files
    bool         saveCodebook = false;    // Save codebook in a separate file
    bool         gaussian = true;         // Gaussian Kernel
    bool         tStudent = false;        // tStudent Kernel

    /* Parameters ============================================================== */

    try
    {

        fn_in = getParameter(argc, argv, "-i");

        if (checkParameter(argc, argv, "-o"))
            fn_out = getParameter(argc, argv, "-o");
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-cvin"))
        {
            if (checkParameter(argc, argv, "-cbin"))
            {
                std::cout << "Error: you can not use two code vectors files" << std::endl;
                exit(EXIT_FAILURE);
            }
            cb_in = getParameter(argc, argv, "-cvin");
            useCBook = false;
        }


        if (checkParameter(argc, argv, "-cbin"))
        {
            if (checkParameter(argc, argv, "-cvin"))
            {
                std::cout << "Error: you can not use two code vectors files" << std::endl;
                exit(EXIT_FAILURE);
            }
            cb_in = getParameter(argc, argv, "-cbin");
            useCBook = true;
        }


        c = textToInteger(getParameter(argc, argv, "-c"));

        if (checkParameter(argc, argv, "-gaussian"))
        {
            if (checkParameter(argc, argv, "-tStudent"))
            {
                std::cout << "Error: you can not define two kernels functions" << std::endl;
                exit(EXIT_FAILURE);
            }
            gaussian = true;
            tStudent = false;
        }
        else if (checkParameter(argc, argv, "-tStudent"))
        {
            gaussian = false;
            tStudent = true;
        }

        df = (int) textToInteger(getParameter(argc, argv, "-df", "3"));

        eps = textToFloat(getParameter(argc, argv, "-eps", "1e-7"));
        iter = textToInteger(getParameter(argc, argv, "-iter", "200"));
        verb = textToInteger(getParameter(argc, argv, "-verb", "0"));

        if (checkParameter(argc, argv, "-norm"))
            norm = true;
        else norm = false;

        if (checkParameter(argc, argv, "-saveclusters"))
            saveClusters = true;
        else saveClusters = false;

        if (checkParameter(argc, argv, "-savecodebook"))
            saveCodebook = true;
        else saveCodebook = false;

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

    if (c < 2)
    {
        std::cerr << argv[0] << ": invalid value for c (must be >= 2): " << c << std::endl;
        exit(EXIT_FAILURE);
    }

    if (df < 2)
    {
        std::cerr << argv[0] << ": invalid value for df (must be > 1): " << df << std::endl;
        exit(EXIT_FAILURE);
    }



    /* Shows parameters ===================================================== */

    std::cout << std::endl << "Parameters used: " << std::endl;
    std::cout << "Input data file : " << fn_in << std::endl;
    std::cout << "Output file name : " << fn_out << std::endl;
    if (cb_in != "")
        std::cout << "Input code vectors file name : " << cb_in << std::endl;
    if (saveClusters)
        std::cout << "Save clusters in separate files: " << fn_out << ".(cluster number)" << std::endl;
    std::cout << "Number of clusters (c) = " << c << std::endl;
    if (gaussian)
        std::cout << "Gaussian Kernel function " << std::endl;
    else
    {
        std::cout << "t-Student Kernel function" << std::endl;
        std::cout << "Degrees of freedom (df) = " << df << std::endl;
    }
    std::cout << "Total number of iterations = " << iter << std::endl;
    std::cout << "Stopping criteria (eps) = " << eps << std::endl;
    std::cout << "verbosity level = " << verb << std::endl;
    if (norm)
        std::cout << "Normalize input data" << std::endl;
    else
        std::cout << "Do not normalize input data " << std::endl;


    /* Open training vector ================================================= */


    std::ifstream inStream(fn_in.c_str());
    if (!inStream)
    {
        std::cerr << argv[0] << ": can't open file " << fn_in << std::endl;
        exit(EXIT_FAILURE);
    }

    xmippCTVectors ts(0, true);
    std::cout << std::endl << "Reading input data file " << fn_in << "....." << std::endl;
    inStream >> ts;



    /* Real stuff ============================================================== */


    try
    {

        if (norm)
        {
            std::cout << "Normalizing....." << std::endl;
            ts.normalize();        // Normalize input data
        }

        xmippFuzzyMap *myMap;

        if (cb_in != "")
        {
            if (useCBook)
            {
                std::cout << "Reading fuzzy codebook file " << cb_in << "....." << std::endl;
                std::ifstream codeStream(cb_in.c_str());
                if (!codeStream)
                {
                    std::cerr << argv[0] << ": can't open file " << cb_in << std::endl;
                    exit(EXIT_FAILURE);
                }
                myMap = new xmippFuzzyMap(codeStream, ts.size(), false);
            }
            else
            {
                std::cout << "Reading fuzzy codevectors file " << cb_in << "....." << std::endl;
                std::ifstream codeStream(cb_in.c_str());
                if (!codeStream)
                {
                    std::cerr << argv[0] << ": can't open file " << cb_in << std::endl;
                    exit(EXIT_FAILURE);
                }
                myMap = new xmippFuzzyMap(codeStream, ts.size(), true);
            }
        }
        else
            myMap = new xmippFuzzyMap("RECT", c, 1, ts);


        xmippKerDenSOM *thisSOM;
        if (fn_algo_in == "")
        {
            if (gaussian)
                thisSOM = new xmippGaussianKerDenSOM(0, 0, 0, eps, iter);        // Creates KerDenSOM Algorithm
            else
                thisSOM = new xmippTStudentKerDenSOM(0, 0, 0, eps, iter, df);    // Creates KerDenSOM Algorithm
        }
        else
        {
            std::cout << "Reading algorithm file " << fn_algo_in << "....." << std::endl << std::endl;
            std::ifstream algoStream(fn_algo_in.c_str());
            if (!algoStream)
            {
                std::cerr << argv[0] << ": can't open file " << fn_algo_in << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        xmippTextualListener myListener;       // Define the listener class
        myListener.setVerbosity() = verb;       // Set verbosity level
        thisSOM->setListener(&myListener);         // Set Listener

        if (cb_in != "")
        {
            if (ts.isNormalized())
            {
                std::cout << "Normalizing code vectors....." << std::endl;
                myMap->Normalize(ts.getNormalizationInfo());       // normalize code vectors
            }
            thisSOM->train(*myMap, ts, fn_out, true);            // Train algorithm
        }
        else
            thisSOM->train(*myMap, ts, fn_out);               // Train algorithm


        // Test algorithm
        xmippFeature dist = thisSOM->test(*myMap, ts);
        std::cout << std::endl << "Quantization error : " <<  dist << std::endl;

        // Classifying
        std::cout << "Classifying....." << std::endl;
        myMap->classify(&ts);

        // Calibrating
        std::cout << "Calibrating....." << std::endl;
        myMap->calibrate(ts);

        /*******************************************************
            Saving all kind of Information
        *******************************************************/

        if (saveCodebook)
        {
            std::cout << "Saving whole codebook as " << fn_out << ".cbk ....." << std::endl;
            tmpN = fn_out.c_str() + (std::string) ".cbk";
            std::ofstream cbkS(tmpN.c_str());
            myMap->saveObject(cbkS);
            cbkS.flush();
        }

        std::cout << "Saving algorithm information as " << fn_out << ".inf ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".inf";
        std::ofstream infS(tmpN.c_str());
        infS << "Kernel Probability Density Estimator clustering algorithm" << std::endl;
        infS << "                 Kernel c-Means" << std::endl << std::endl;
        infS << "Input data file : " << fn_in << std::endl;
        if (cb_in != "")
            infS << "Input code vectors file : " << cb_in << std::endl;
        infS << "Code vectors output file : " << fn_out <<  ".cod" << std::endl;
        infS << "Whole codebook output file : " << fn_out <<  ".cbk" << std::endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
        infS << "Number of feature vectors: " << ts.size() << std::endl;
        infS << "Number of variables: " << ts.theItems[0].size() << std::endl;
        infS << "Number of clusters = " << c << std::endl;
        if (norm)
            infS << "Input data normalized" << std::endl;
        else
            infS << "Input data not normalized" << std::endl;

        if (gaussian)
            infS << "Gaussian Kernel function " << std::endl;
        else
        {
            infS << "t-Student Kernel function" << std::endl;
            infS << "Degrees of freedom (df) = " << df << std::endl;
        }

        infS << "Total number of iterations = " << iter << std::endl;
        infS << "Stopping criteria (eps) = " << eps << std::endl;
        infS << "Final Sigma = " << thisSOM->getSigma() << std::endl;
        infS << "Quantization error : " <<  dist << std::endl;
        infS.flush();

        // assign data to clusters according to fuzzy threshold
        if (saveClusters)
        {
            std::cout << "Saving neurons assigments ....." << std::endl;
            for (unsigned i = 0; i < myMap->size(); i++)
            {
                tmpN = fn_out.c_str() + (std::string) "."  + integerToString(i);
                std::ofstream cStream(tmpN.c_str());
                for (int j = 0; j < myMap->classifAt(i).size(); j++)
                    cStream << myMap->classifAt(i)[j] << std::endl;
                cStream.flush();
            }
        }

        // save .vs file to be compatible with SOM_PAK
        std::cout << "Saving visual file as " << fn_out << ".vs ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".vs";
        std::ofstream vsStream(tmpN.c_str());
        vsStream << ts.theItems[0].size() << " " << myMap->layout() << " " << myMap->width() << " " << myMap->height() << " gaussian" << std::endl;
        for (int i = 0; i < ts.size(); i++)
        {
            int j = myMap->fuzzyWinner(i);
            vsStream << myMap->indexToPos(j).first << " " << myMap->indexToPos(j).second << " " << myMap->memb[i][j] << " " << ts.theTargets[i] << std::endl;
        }
        vsStream.flush();

        // save .his file (Histogram)
        std::cout << "Saving code vectors histogram file as " << fn_out << ".his ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".his";
        std::ofstream hisStream(tmpN.c_str());
        myMap->printHistogram(hisStream);
        hisStream.flush();

        // save .err file (Average Quantization Error)
        std::cout << "Saving code vectors average quantization error file as " << fn_out << ".err ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".err";
        std::ofstream errStream(tmpN.c_str());
        myMap->printQuantError(errStream);
        errStream.flush();


        if (norm)
        {
            std::cout << "Denormalizing code vectors....." << std::endl;
            myMap->unNormalize(ts.getNormalizationInfo()); // de-normalize codevectors
        }

        std::cout << "Saving code vectors as " << fn_out << ".cod ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".cod";
        std::ofstream codS(tmpN.c_str());
        codS << *myMap;
        codS.flush();


        std::cout << std::endl;

        delete myMap;
        delete thisSOM;

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
        "\nPurpose: Kernel Density Estimator Clustering algorithm"
        "\nParameter Values: (note space before value)"
        "\n    -i      file_in           Input data file (plain data)"
        "\n    -o      file_out          Base name for output data files"
        "\n    -cvin   file_in           Codevectors input file"
        "\n    -saveclusters           save clusters in separate files (Default = No)"
        "\n    -c      Clusters       Number of clusters"
        "\n    -gaussian           Gaussian Kernel Function (default)"
        "\n    -tStudent           t-Student Kernel Function "
        "\n    -df     df               t-Student degrees of freedom (default = 3)"
        "\n    -eps    Epsilon        Stopping criteria (default = 1e-7)"
        "\n    -iter   iterations        Number of iterations (default = 200)"
        "\n    -norm                   Normalize training data (default = No)"
        "\n    -verb   verbosity         Information level while running: "
        "\n             0: No information (default)"
        "\n             1: Progress bar"
        "\n             2: Changes between iterations"
        "\n"
        , argv[0]);
    exit(0);
}
