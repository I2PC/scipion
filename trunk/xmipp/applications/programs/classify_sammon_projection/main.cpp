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

// to avoid warnings with long class names
#pragma warning(disable:4786)

#include <fstream>

#include <data/args.h>
#include <classification/training_vector.h>
#include <classification/sammon.h>

/**
 * This test only reads a training set from a file and prints it
 */

//-----------------------------------------------------------------------------

typedef xmippSammon Sammon;

/* Prototypes -============================================================= */

void Usage(char **argv);

//-----------------------------------------------------------------------------

int main(int argc, char** argv)
{

    /* Input Parameters ======================================================== */
    FileName       fn_in;     // input file
    FileName       fn_out;    // output file
    FileName       tmpN;     // temporary file name
    unsigned       dim = 2;  // Dimension of output space
    unsigned       iter = 100000; // Iteration number
    unsigned       verb = 0; // Verbosity level
    bool           norm = 1; // Normalize?

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

        dim = textToInteger(getParameter(argc, argv, "-dim", "2"));
        iter = textToInteger(getParameter(argc, argv, "-iter", "100000"));
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


    std::cout << std::endl << "Reading file " << fn_in << "....." << std::endl;

    std::ifstream inStream(fn_in.c_str());
    if (!inStream)
    {
        std::cerr << argv[0] << ": can't open file " << fn_in << std::endl;
        exit(EXIT_FAILURE);
    }


    /* Some validations ===================================================== */

    Sammon::In in(0, true);
    try
    {
        inStream >> in;
    }
    catch (std::exception& e)
    {
        std::cerr << argv[0] << ": can't read file " << fn_in  << " because " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    if (dim < 0 || dim >= in.theItems[0].size())
    {
        std::cerr << argv[0] << ": invalid value for dim " << dim << std::endl;
        exit(EXIT_FAILURE);
    }

    if (iter < 1)
    {
        std::cerr << argv[0] << ": invalid value for iter " << iter << std::endl;
        exit(EXIT_FAILURE);
    }

    if (verb < 0 || verb > 2)
    {
        std::cerr << argv[0] << ": invalid value for verbosity " << verb << std::endl;
        exit(EXIT_FAILURE);
    }

    /* Real stuff ============================================================== */


    std::cout << "Parameters used: " << std::endl;
    std::cout << "Input data file : " << fn_in << std::endl;
    std::cout  << "Output file : " << fn_out <<  ".sam" << std::endl;
    std::cout << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
    std::cout << "Output space dimension = " << dim << std::endl;
    if (norm)
        std::cout << "Input data normalized" << std::endl;
    else
        std::cout << "Input data not normalized" << std::endl;
    std::cout << "Total number of iterations = " << iter << std::endl << std::endl;

    if (norm)
    {
        std::cout << "Normalizing....." << std::endl;
        in.normalize();        // Normalize input data
    }
    Sammon::Out out(in.theItems[0].size()); // create output data set
    Sammon sammon(dim, iter, 0.3);        // Crate algorithm
    xmippTextualListener myListener; // Define the listener class
    myListener.setVerbosity() = verb; // Set verbosity level
    sammon.setListener(&myListener); // Set Listener
    sammon(in, out);   // do mapping

    std::cout << "Saving mapped data in " << fn_out << ".sam....." << std::endl;

    tmpN = fn_out.c_str() + (std::string) ".sam";
    std::ofstream outStream(tmpN.c_str());
    outStream << out;
    outStream.flush();
    xmippFeature stress = sammon.getStress();
    std::cout << "Sammon stress (Distance error) : " << stress << std::endl;

    std::cout << "Saving algorithm information as " << fn_out << ".inf ....." << std::endl;
    tmpN = fn_out.c_str() + (std::string) ".inf";
    std::ofstream infS(tmpN.c_str());
    infS << "Sammon algorithm" << std::endl << std::endl;
    infS << "Input data file : " << fn_in << std::endl;
    infS << "Output file : " << fn_out <<  ".cod" << std::endl;
    infS << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
    infS << "Number of feature vectors: " << in.size() << std::endl;
    infS << "Number of variables: " << in.theItems[0].size() << std::endl;
    infS << "Output space dimension = " << dim << std::endl;
    if (norm)
        infS << "Input data normalized" << std::endl;
    else
        infS << "Input data not normalized" << std::endl;
    infS << "Total number of iterations = " << iter << std::endl;
    infS << "Sammon stress (mapping error) : " <<  stress << std::endl;
    infS.flush();


    return 0;

}


/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage(char **argv)
{
    printf(
        "\nUsage: %s [Purpose and Parameters]"
        "\nPurpose: Non-linear mapping of high dimensional input space"
        "\n  into a low-dimensional one"
        "\n"
        "\nParameter Values: (note space before value)"
        "\n"
        "\n    -i    file_in        input data file"
        "\n    -o    file_out       output data file (mapped data)"
        "\n    -dim  dimension      dimension of the output space (default = 2)"
        "\n    -iter iterations    number of iterations (default = 100000)"
        "\n    -norm       Normalize training data (default: No)"
        "\n    -verb verbosity    information level while running: "
        "\n        0: No information (default)"
        "\n        1: Progress bar"
        "\n        2: sammon stress (distance error)"
        "\n    \n"
        , argv[0]);
}
