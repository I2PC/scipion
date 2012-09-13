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

#include <classification/pca.h>

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
    unsigned       verb = 0; // Verbosity level

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

        verb = textToInteger(getParameter(argc, argv, "-verb", "1"));


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


    if (verb < 0 || verb > 2)
    {
        std::cerr << argv[0] << ": invalid value for verbosity (must be between 0 and 2): " << verb << std::endl;
        exit(EXIT_FAILURE);
    }


    /* Shows parameters ===================================================== */

    std::cout << std::endl << "Parameters used: " << std::endl;
    std::cout << "Input data file : " << fn_in << std::endl;
    std::cout << "Output file name : " << fn_out << std::endl;
    std::cout << "verbosity level = " << verb << std::endl;


    /* Open training vector ================================================= */


    std::ifstream inStream(fn_in.c_str());
    if (!inStream)
    {
        std::cerr << argv[0] << ": can't open file " << fn_in << std::endl;
        exit(EXIT_FAILURE);
    }

    xmippCTVectors ts(0, false);

    std::cout << std::endl << "Reading data file " << fn_in << "....." << std::endl;
    inStream >> ts;


    /* Real stuff ============================================================== */


    try
    {

        // Define PCA class and do the projection
        xmippPC myPCA;

        xmippTextualListener myListener;     // Define the listener class
        myListener.setVerbosity() = verb;     // Set verbosity level
        myPCA.setListener(&myListener);          // Set Listener

        myPCA.reset(ts);                         // find eigenvectors and eigenvalues

        /*******************************************************
            Saving all kind of Information
        *******************************************************/

        std::cout << std::endl << "Saving eigenvalues as " << fn_out << ".eval ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".eval";
        std::ofstream evalS(tmpN.c_str());
        evalS << "3  " << myPCA.eigenvec.size() << std::endl;
        double cum = 0;
        for (int i = 0; i < myPCA.eigenval.size(); i++) cum += myPCA.eigenval[i];
        if (cum == 0) cum = 1;
        for (int i = 0; i < myPCA.eigenval.size(); i++)
            evalS << i << " " << myPCA.eigenval[i] << " " << myPCA.eigenval[i] / cum << std::endl;
        evalS.flush();

        std::cout << "Saving eigenvectors as " << fn_out << ".evec ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".evec";
        std::ofstream evecS(tmpN.c_str());
        evecS << myPCA.eigenvec[0].size() << " " << myPCA.eigenvec.size() << std::endl;
        for (int j = 0; j < myPCA.eigenvec.size(); j++)
        {
            for (int i = 0; i < myPCA.eigenvec[j].size(); i++)
                evecS << " " << myPCA.eigenvec[j][i];
            evecS << std::endl;
        }
        evecS.flush();

        std::cout << "Saving algorithm information as " << fn_out << ".inf ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".inf";
        std::ofstream infS(tmpN.c_str());
        infS << "PCA" << std::endl << std::endl;
        infS << "Input data file : " << fn_in << std::endl;
        infS << "Eigenvalues output file : " << fn_out <<  ".eval" << std::endl;
        infS << "Eigenvectors output file : " << fn_out <<  ".evec" << std::endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
        infS << "Number of feature vectors: " << ts.size() << std::endl;
        infS << "Number of variables: " << ts.itemAt(0).size() << std::endl;

        infS.flush();


        std::cout << std::endl;

    }
    catch (Xmipp_error &e)
    {
        std::cout << e << std::endl;
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
        "\nPurpose: Principal Component Analysis (PCA)"
        "\n"
        "\nParameter Values: (note space before value)"
        "\n"
        "\n    -i      file_in           Input data file (plain data)"
        "\n    -o      file_out          Base name for output data files:"
        "\n    -verb   level       Show progress bar: "
        "\n             0: Do not show the progress bar"
        "\n             1: Show the progress bar (default)"
        "\n      \n"
        , argv[0]);
    exit(0);
}
