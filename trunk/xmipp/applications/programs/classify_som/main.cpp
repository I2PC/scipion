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

/* Part of this code were developed by Lorenzo Zampighi and Nelson Tang
   of the department of Physiology of the David Geffen School of Medicine,
   University of California, Los Angeles
*/


// To avoid problems with long template names
#pragma warning(disable:4786)

#include <fstream>

#include <classification/som.h>

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
    unsigned       iter = 10000; // Iteration number
    unsigned       verb = 0; // Verbosity level
    bool           norm = 1; // Normalize?
    unsigned       xdim = 10; // X-dimension (-->)
    unsigned       ydim = 5; // Y-dimension
    float          alpha_0 = 0.1; // Initial alpha value
    float          radius_0; // Initial radius value
    std::string    layout = "HEXA"; // Type of layout (Topology)
    bool           gaussian = true; // Gaussian kernel
    bool           bubble = false;  // Bubble kernel
    bool           saveClusters = false;    // Save clusters in separate files
    bool           use_rand_cvs = false; // NT: flag to truly randomize codevectors or not


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
            cb_in = getParameter(argc, argv, "-cvin");

        ydim = textToInteger(getParameter(argc, argv, "-ydim", "5"));
        xdim = textToInteger(getParameter(argc, argv, "-xdim", "10"));

        if (checkParameter(argc, argv, "-hexa"))
        {
            if (checkParameter(argc, argv, "-rect"))
            {
                std::cout << "Error: you can not define two topologies" << std::endl;
                exit(EXIT_FAILURE);
            }
            layout = "HEXA";
        }
        else if (checkParameter(argc, argv, "-rect"))
            layout = "RECT";


        if (checkParameter(argc, argv, "-gaussian"))
        {
            if (checkParameter(argc, argv, "-bubble"))
            {
                std::cout << "Error: you can not define two kernels" << std::endl;
                exit(EXIT_FAILURE);
            }
            gaussian = true;
            bubble = false;
        }
        else if (checkParameter(argc, argv, "-bubble"))
        {
            bubble = true;
            gaussian = false;
        }


        alpha_0 = textToFloat(getParameter(argc, argv, "-alpha", "0.1"));

        if (checkParameter(argc, argv, "-radius"))
            radius_0 = textToFloat(getParameter(argc, argv, "-radius"));
        else if (xdim > ydim)
            radius_0 = xdim;
        else
            radius_0 = ydim;

        iter = textToInteger(getParameter(argc, argv, "-iter", "10000"));
        verb = textToInteger(getParameter(argc, argv, "-verb", "0"));

        if (checkParameter(argc, argv, "-norm"))
            norm = true;
        else norm = false;

        if (checkParameter(argc, argv, "-saveclusters"))
            saveClusters = true;
        else saveClusters = false;

        if (checkParameter(argc, argv, "-randomcodevectors"))
            use_rand_cvs = true;
        else use_rand_cvs = false;

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

    if (alpha_0 < 0)
    {
        std::cerr << argv[0] << ": invalid value for alpha (must be > 0): " << alpha_0 << std::endl;
        exit(EXIT_FAILURE);
    }

    if (radius_0 < 1)
    {
        std::cerr << argv[0] << ": invalid value for radius (must be > 1): " << radius_0 << std::endl;
        exit(EXIT_FAILURE);
    }

    if (xdim < 1)
    {
        std::cerr << argv[0] << ": invalid value for xdim (must be > 1): " << xdim << std::endl;
        exit(EXIT_FAILURE);
    }

    if (ydim < 1)
    {
        std::cerr << argv[0] << ": invalid value for ydim (must be > 1): " << ydim << std::endl;
        exit(EXIT_FAILURE);
    }


    /* Shows parameters ===================================================== */

    std::cout << std::endl << "Parameters used: " << std::endl;
    std::cout << "Input data file : " << fn_in << std::endl;
    std::cout << "Output file name : " << fn_out << std::endl;
    if (cb_in != "")
        std::cout << "Input code vectors file name : " << cb_in << std::endl;
    else if (use_rand_cvs)
        std::cout << "Using randomized code vectors" << std::endl;
    if (saveClusters)
        std::cout << "Save clusters in separate files: " << fn_out << ".(cluster number)" << std::endl;
    std::cout << "Horizontal dimension (Xdim) = " << xdim << std::endl;
    std::cout << "Vertical dimension (Ydim) = " << ydim << std::endl;
    if (layout == "HEXA")
        std::cout << "Hexagonal topology " << std::endl;
    else
        std::cout << "Rectangular topology " << std::endl;

    if (gaussian)
        std::cout << "Gaussian neighborhood function " << std::endl;
    else
        std::cout << "Bubble neighborhood function" << std::endl;
    std::cout << "Initial learning rate (alpha) = " << alpha_0 << std::endl;
    std::cout << "Initial neighborhood radius (radius) = " << radius_0 << std::endl;

    std::cout << "Total number of iterations = " << iter << std::endl;
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

    std::cout << std::endl << "Reading data file " << fn_in << "....." << std::endl;
    inStream >> ts;


    /* Real stuff ============================================================== */


    try
    {
        if (norm)
        {
            std::cout << "Normalizing....." << std::endl;
            ts.normalize();        // Normalize input data
        }

        xmippMap* myMap;

        if (cb_in != "")
        {
            std::cout << "Reading codevectors file " << cb_in << "....." << std::endl;
            std::ifstream codeStream(cb_in.c_str());
            if (!codeStream)
            {
                std::cerr << argv[0] << ": can't open file " << cb_in << std::endl;
                exit(EXIT_FAILURE);
            }
            myMap = new xmippMap(codeStream);
        }
        else
            myMap = new xmippMap(layout, xdim, ydim, ts, use_rand_cvs);


        xmippSOM *thisSOM;
        Descent alpha(alpha_0, 0);         // alpha decreases linearly to 0
        Descent radius(radius_0, 1);       // radius decreases linearly to 1

        xmippSOM::neighType neigh;
        if (gaussian)
            neigh = xmippSOM::GAUSSIAN;
        else
            neigh = xmippSOM::BUBBLE;
        thisSOM = new xmippSOM(alpha, radius, neigh, iter);    // Creates SOM Algorithm

        xmippTextualListener myListener;     // Define the listener class
        myListener.setVerbosity() = verb;     // Set verbosity level
        thisSOM->setListener(&myListener);       // Set Listener

        thisSOM->train(*myMap, ts);              // Train algorithm

        // Test algorithm
        std::cout << std::endl;
        double dist = thisSOM->test(*myMap, ts);
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

        std::cout << "Saving algorithm information as " << fn_out << ".inf ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".inf";
        std::ofstream infS(tmpN.c_str());
        infS << "Kohonen SOM algorithm" << std::endl << std::endl;
        infS << "Input data file : " << fn_in << std::endl;
        if (cb_in != "")
            infS << "Input code vectors file : " << cb_in << std::endl;
        else if (use_rand_cvs)
            infS << "Using randomized code vectors" << std::endl;
        infS << "Code vectors output file : " << fn_out <<  ".cod" << std::endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
        infS << "Number of feature vectors: " << ts.size() << std::endl;
        infS << "Number of variables: " << ts.itemAt(0).size() << std::endl;
        infS << "Horizontal dimension (Xdim) = " << xdim << std::endl;
        infS << "Vertical dimension (Ydim) = " << ydim << std::endl;
        if (layout == "HEXA")
            infS << "Hexagonal topology " << std::endl;
        else
            infS << "Rectangular topology " << std::endl;
        if (gaussian)
            infS << "Gaussian neighborhood function " << std::endl;
        else
            infS << "Bubble neighborhood function" << std::endl;
        infS << "Initial learning rate (alpha) = " << alpha_0 << std::endl;
        infS << "Initial neighborhood radius (radius) = " << radius_0 << std::endl;
        infS << "Total number of iterations = " << iter << std::endl;
        if (norm)
            infS << "Input data normalized" << std::endl;
        else
            infS << "Input data not normalized" << std::endl;
        infS << "Quantization error : " <<  dist << std::endl;

        infS.flush();

        // assign data to neurons
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
        vsStream << ts.theItems[0].size() << " " << myMap->layout() << " " << myMap->width() << " " << myMap->height();
        if (gaussian)
            vsStream  << " gaussian" << std::endl;
        else
            vsStream  << " bubble" << std::endl;
        for (int i = 0; i < ts.size(); i++)
        {
            int j = myMap->winner(ts, i);
            vsStream << myMap->indexToPos(j).first << " " << myMap->indexToPos(j).second << " " << eDist(myMap->theItems[j], ts.theItems[i]) << " " << ts.theTargets[i] << std::endl;
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
        "\nPurpose: Kohonen Self-Organizing Feature Maps"
        "\n"
        "\nParameter Values: (note space before value)"
        "\n"
        "\n    -i      file_in           Input data file (plain data)"
        "\n    -o      file_out          Base name for output data files:"
        "\n    -cvin   file_in           Codevectors input file"
        "\n    -saveclusters           save clusters in separate files (Default = No)"
        "\n    -xdim   H-dimension       Horizontal size of the map (default = 10)"
        "\n    -ydim   V-dimension       Vertical size of the map (default = 5)"
        "\n    -hexa               Hexagonal topology (default)"
        "\n    -rect               Rectangular topology"
        "\n    -gaussian          Gaussian neighborhood learning kernel (Default)"
        "\n    -bubble            Bubble neighborhood learning kernel"
        "\n    -alpha     learning rate  Initial learning rate value (default = 0.1)"
        "\n    -radius    radius        Initial neighborhood radius"
        "\n             (default = max(xdim, ydim))"
        "\n    -iter      iterations     Number of iterations (default = 10000)"
        "\n    -norm                   Normalize training data (default)"
        "\n    -verb      verbosity      Information level while running: "
        "\n             0: No information (default)"
        "\n             1: Progress bar"
        "\n             2: Code vectors change between iterations"
        "\n    -randomcodevectors        Use truly randomized codevectors (default FALSE)"
        "\n      \n"
        , argv[0]);
    exit(0);
}
