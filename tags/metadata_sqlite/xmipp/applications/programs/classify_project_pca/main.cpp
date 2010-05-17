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

#include <classification/pca.h>

/* Prototypes -============================================================= */

void Usage(char **argv);

/* Main function -============================================================= */

main(int argc, char** argv)
{


    /* Input Parameters ======================================================== */

    FileName       fn_in;     // input file
    FileName       fn_out;    // output file
    FileName       fn_ein;     // eigen vectors input file
    FileName       fn_evin;     // eigen values input file
    FileName       cb_in = "";    // Code vectors input file
    FileName       tmpN;  // Temporary variable
    unsigned       k;               // Dimension of projected subspace
    float          p;               // Percent of explained variance
    bool           met;             // project by k?
    bool           recon;           // Reconstruct original data using k components?

    /* Parameters ============================================================== */
    try
    {

        fn_in = getParameter(argc, argv, "-i");
        fn_ein = getParameter(argc, argv, "-ein");
        fn_evin = getParameter(argc, argv, "-evin");

        if (checkParameter(argc, argv, "-o"))
            fn_out = getParameter(argc, argv, "-o");
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-recon"))
            recon = true;
        else
            recon = false;

        if (checkParameter(argc, argv, "-k"))
        {
            if (checkParameter(argc, argv, "-p"))
            {
                std::cerr << argv[0] << ": Invalid option. You can not select number of dimensions and percent at the same time" << std::endl;
                exit(EXIT_FAILURE);
            }
            k = textToInteger(getParameter(argc, argv, "-k"));
            met = true;
        }

        if (checkParameter(argc, argv, "-p"))
        {
            if (checkParameter(argc, argv, "-k"))
            {
                std::cerr << argv[0] << ": Invalid option. You can not select number of dimensions and percent at the same time" << std::endl;
                exit(EXIT_FAILURE);
            }
            p = textToFloat(getParameter(argc, argv, "-p"));
            met = false;
        }



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


    if (!met && (p <= 0 || p > 100))
    {
        std::cerr << argv[0] << ": invalid value for percent of explained variance (must be > 0 and <= 100): " << p << std::endl;
        exit(EXIT_FAILURE);
    }

    /* Shows parameters ===================================================== */

    std::cout << std::endl << "Parameters used: " << std::endl;
    std::cout << "Input data file : " << fn_in << std::endl;
    std::cout << "Input eigen vectors file : " << fn_ein << std::endl;
    std::cout << "Input eigen values file : " << fn_evin << std::endl;
    std::cout << "Output file name : " << fn_out << std::endl;
    if (met)
        std::cout << "Dimension of projected subspace : " << k << std::endl;
    else
        std::cout << "Percent of explained variance : " << p << std::endl;
    if (recon)
        std::cout << "Reconstruct original data" << std::endl;


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

    if (met)
    {
        if (k > ts.theItems[0].size() || k <= 0)
        {
            std::cerr << argv[0] << ": invalid value for dimension of projected subspace (must be > 0 and <= dimension of the data): " << k << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    /* Open eigen vectors================================================= */


    std::ifstream einStream(fn_ein.c_str());
    if (!einStream)
    {
        std::cerr << argv[0] << ": can't open file " << fn_ein << std::endl;
        exit(EXIT_FAILURE);
    }

    xmippCTVectors ev(0, false);

    std::cout << "Reading eigen vectors file " << fn_ein << "....." << std::endl;
    einStream >> ev;

    /* Open eigen values================================================= */

    std::ifstream evinStream(fn_evin.c_str());
    if (!evinStream)
    {
        std::cerr << argv[0] << ": can't open file " << fn_evin << std::endl;
        exit(EXIT_FAILURE);
    }

    xmippCTVectors eval(0, false);

    std::cout << "Reading eigen values file " << fn_evin << "....." << std::endl;
    evinStream >> eval;

    /* Real stuff ============================================================== */

    if (!met)
    {
// Calculate the number of components that satisfy p% of the variance
        double cum = 0;
        int bestIndex;
        for (int i = 0; i < eval.size(); i++)
        {
            cum += eval.theItems[i][2] * 100.0;
            if (cum >= p)
            {
                bestIndex = i;
                p = cum;
                break;
            }
        }
        k = bestIndex + 1;
        std::cout << "Explained " << p << "% of the variance with " << k << " components" << std::endl;
    }
    else
    {
// Calculate the percent of variance explained by k components
        double cum = 0;
        for (int i = 0; i < k; i++)
            cum += eval.theItems[i][2] * 100.0;
        p = cum;
        std::cout << "The " << k << " components explain " << p << "% of the variance." << std::endl;
    }

    try
    {

        // Do the projection

// Calculate mean of the vectors
        xmippCTVectors statVec(0, true);
        statVec = ts.getStatVector();

        std::cout << "projecting into " << k << " dimensions..." << std::endl;
        int size = ts.size();
        int dim = ts.itemAt(0).size();

        xmippCTVectors projectedTs(0, true);
        projectedTs.theItems.resize(size);
        projectedTs.theTargets.resize(size);
        for (int h = 0; h < size; h++) projectedTs.theItems[h].resize(k, 0);

        for (int i = 0; i < k; i++)
        {
            for (int z = 0; z < size; z++)
            {
                double cum = 0;
                for (int j = 0; j < dim; j++)
                {
                    cum += ev.theItems[i][j] * (ts.theItems[z][j] - statVec.theItems[0][j]);
                } // j
                projectedTs.theItems[z][i] = cum;
            }  // z
        } // i

        xmippCTVectors reconsTs(0, true);
        if (recon)
        {
            std::cout << "Estimating original data using " << k << " components..." << std::endl;
            reconsTs.theItems.resize(size);
            reconsTs.theTargets.resize(size);
            for (int h = 0; h < size; h++) reconsTs.theItems[h].resize(dim, 0);

            for (int i = 0; i < dim; i++)
            {
                for (int z = 0; z < size; z++)
                {
                    double cum = 0;
                    for (int j = 0; j < k; j++)
                    {
                        cum += ev.theItems[j][i] * projectedTs.theItems[z][j];
                    } // j
                    cum += statVec.theItems[0][i];
                    reconsTs.theItems[z][i] = cum;
                }  // z
            } // i
        } // if recon

        /*******************************************************
            Saving all kind of Information
        *******************************************************/

        std::cout << "Saving projected vectors as " << fn_out << ".dat ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".dat";
        std::ofstream projS(tmpN.c_str());
        projS << k << " " << size << std::endl;
        for (int j = 0; j < size ; j++)
        {
            for (int i = 0; i < k; i++)
                projS << " " << projectedTs.theItems[j][i];
            projS << " " << ts.theTargets[j] << std::endl;
        }
        projS.flush();

        if (recon)
        {
            std::cout << "Saving reconstructed vectors as " << fn_out << ".recon ....." << std::endl;
            tmpN = fn_out.c_str() + (std::string) ".recon";
            std::ofstream rS(tmpN.c_str());
            rS << dim << " " << size << std::endl;
            for (int j = 0; j < size ; j++)
            {
                for (int i = 0; i < dim; i++)
                    rS << " " << reconsTs.theItems[j][i];
                rS << " " << ts.theTargets[j] << std::endl;
            }
            rS.flush();
        }

        std::cout << "Saving algorithm information as " << fn_out << ".inf ....." << std::endl;
        tmpN = fn_out.c_str() + (std::string) ".inf";
        std::ofstream infS(tmpN.c_str());
        infS << "PCA Projection" << std::endl << std::endl;
        infS << "Input data file : " << fn_in << std::endl;
        infS << "Projected vectors output file : " << fn_out <<  ".dat" << std::endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
        infS << "Number of feature vectors: " << ts.size() << std::endl;
        infS << "Number of variables: " << ts.itemAt(0).size() << std::endl;
        if (met)
        {
            infS << "Dimension of projected subspace : " << k << std::endl;
            infS << "Percent of explained variance : " << p << std::endl;
        }
        else
        {
            infS << "Percent of explained variance : " << p << std::endl;
            infS << "Number of final components : " << k << std::endl;
        }
        if (recon)
            infS << "Reconstruct original data" << std::endl;
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
        "\nPurpose: Project data into the first k Principal Components (PCA projection)"
        "\n"
        "\nParameter Values: (note space before value)"
        "\n"
        "\n    -i      file_in           Input data file (plain data)"
        "\n    -o      file_out          Base name for output data files"
        "\n    -ein    file_in           Eigen vectors file"
        "\n    -evin   file_in           Eigen values file"
        "\n    -k      k                 Number of components"
        "\n    -p      p                 Percent of the explained variance"
        "\n    -recon  true/false        If true, reconstruct original data"
        "\n             using the k principal components"
        "\n             (default: false)"
        "\n      \n"
        , argv[0]);
    exit(0);
}

