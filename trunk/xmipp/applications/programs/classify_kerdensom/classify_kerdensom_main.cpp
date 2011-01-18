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

#include <data/program.h>
#include <classification/tstudent_kerdensom.h>
#include <classification/gaussian_kerdensom.h>

/* Parameters class ======================================================= */
class ProgKenderSOM: public XmippProgram
{
public:
    /* Input Parameters ======================================================== */
    FileName       fn_in;        // input file
    FileName       fn_out;       // output file
    FileName       cb_in;        // Code vectors input file
    FileName       fn_algo_in;   // input algorithm file
    FileName       tmpN;         // Temporary variable
    double         eps;          // Stopping criteria
    unsigned       iter;         // Iteration number
    bool           norm;         // Normalize?
    unsigned       xdim;         // X-dimension (-->)
    unsigned       ydim;         // Y-dimension
    double         reg0;         // Initial reg
    double         reg1;         // Final reg
    int            df;           // Degrees of freedom
    std::string    layout;       // layout (Topology)
    unsigned       annSteps;     // Deterministic Annealing steps
    bool           useCBook;     // Use codebook
    bool           saveClusters; // Save clusters in separate files
    bool           saveCodebook; // Save codebook in a separate file
    bool           gaussian;     // Gaussian Kernel
    bool           tStudent;     // tStudent Kernel
public:
    // Define parameters
    void defineParams()
    {
        addUsageLine("Purpose: Kernel Density Estimator Self-Organizing Map");
        addParamsLine("  -i <file_in>                : Input data file (plain data)");
        addParamsLine(" [-o <file_out>]              : Base name for output data files");
        addParamsLine(" [--cvin <file_in>]           : Codevectors input file");
        addParamsLine(" [--cbin <file_in>]           : Codebook input file");
        addParamsLine(" [--saveclusters]             : Save clusters in separate files");
        addParamsLine(" [--savecodebook]             : Save code book");
        addParamsLine(" [--xdim <Hdimension=10>]     : Horizontal size of the map");
        addParamsLine(" [--ydim <Vdimension=5>]      : Vertical size of the map");
        addParamsLine(" [--topology <topology=RECT>] : Lattice topology");
        addParamsLine("    where <topology> RECT HEXA");
        addParamsLine(" [--steps <steps=10>]         : Deterministic annealing steps");
        addParamsLine(" [--reg0  <Initial_reg=1000>] : Initial smoothness factor");
        addParamsLine(" [--reg1  <Final_reg=100>]    : Final  smoothness factor");
        addParamsLine(" [--kernel <kernel=gaussian>] : Kernel function");
        addParamsLine("    where <kernel> gaussian tstudent");
        addParamsLine(" [--df <df=3>]                : t-Student degrees of freedom");
        addParamsLine(" [--ain <algorithmFile>]      : algorithm input file");
        addParamsLine(" [--eps <epsilon=1e-7>]       : Stopping criteria");
        addParamsLine(" [--iter <N=200>]             : Number of iterations");
        addParamsLine(" [--norm]                     : Normalize training data");
    }

    // Read parameters
    void readParams()
    {
        fn_in = getParam("-i");
        if (checkParam("-o"))
            fn_out=getParam("-o");
        if (checkParam("--cvin") && checkParam("--cbin"))
            REPORT_ERROR(ERR_ARG_INCORRECT,"Cannot provide --cvin and --cbin");
        if (checkParam("--cvin"))
        {
            cb_in=getParam("--cvin");
            useCBook = false;
        }
        else if (checkParam("--cbin"))
        {
            cb_in=getParam("--cbin");
            useCBook = true;
        }
        ydim = getIntParam("--ydim");
        xdim = getIntParam("--xdim");
        layout = getParam("--topology");
        std::string kernel=getParam("--kernel");
        if (kernel=="gaussian")
        {
            gaussian = true;
            tStudent = false;
        }
        else if (kernel=="tstudent")
        {
            gaussian = false;
            tStudent = true;
        }
        reg0 = getDoubleParam("--reg0");
        reg1 = getDoubleParam("--reg1");
        df = getIntParam("--df");
        if (checkParam("--ain"))
            fn_algo_in=getParam("--ain");
        eps = getDoubleParam("--eps");
        iter = getIntParam("--iter");
        norm = checkParam("--norm");
        saveClusters = checkParam("--saveclusters");
        saveCodebook = checkParam("--savecodebook");
        annSteps = getIntParam("--steps");

        // Some checks
        if (iter < 1)
            REPORT_ERROR(ERR_ARG_INCORRECT,"iter must be > 1");

        if ((reg0 <= reg1) && (reg0 != 0) && (annSteps > 1))
            REPORT_ERROR(ERR_ARG_INCORRECT,"reg0 must be > reg1");
        if (reg0 == 0)
            annSteps = 0;
        if (reg0 < 0)
            REPORT_ERROR(ERR_ARG_INCORRECT,"reg0 must be > 0");
        if (reg1 < 0)
            REPORT_ERROR(ERR_ARG_INCORRECT,"reg1 must be > 0");
        if (xdim < 1)
            REPORT_ERROR(ERR_ARG_INCORRECT,"xdim must be >= 1");
        if (ydim < 1)
            REPORT_ERROR(ERR_ARG_INCORRECT,"ydim must be >= 1");
        if (df < 2)
            REPORT_ERROR(ERR_ARG_INCORRECT,"df must be > 1");
    }

    void show()
    {
        std::cout << "Input data file : " << fn_in << std::endl;
        std::cout << "Output file name : " << fn_out << std::endl;
        if (cb_in != "")
            std::cout << "Input code vectors file name : " << cb_in << std::endl;
        if (saveClusters)
            std::cout << "Save clusters in separate files: " << fn_out << ".(cluster number)" << std::endl;
        std::cout << "Horizontal dimension (Xdim) = " << xdim << std::endl;
        std::cout << "Vertical dimension (Ydim) = " << ydim << std::endl;
        if (layout == "HEXA")
            std::cout << "Hexagonal topology " << std::endl;
        else
            std::cout << "Rectangular topology " << std::endl;
        std::cout << "Initial smoothness factor (reg0) = " << reg0 << std::endl;
        std::cout << "Final smoothness factor (reg1) = " << reg1 << std::endl;
        if (gaussian)
            std::cout << "Gaussian Kernel function " << std::endl;
        else
        {
            std::cout << "t-Student Kernel function" << std::endl;
            std::cout << "Degrees of freedom (df) = " << df << std::endl;
        }
        std::cout << "Deterministic annealing steps = " << annSteps << std::endl;
        std::cout << "Total number of iterations = " << iter << std::endl;
        std::cout << "Stopping criteria (eps) = " << eps << std::endl;
        if (norm)
            std::cout << "Normalize input data" << std::endl;
        else
            std::cout << "Do not normalize input data " << std::endl;
    }

    // Run
    void run()
    {
        /* Open training vector ================================================= */
        std::ifstream inStream(fn_in.c_str());
        if (!inStream)
        {
            std::cerr << argv[0] << ": can't open file " << fn_in << std::endl;
            exit(EXIT_FAILURE);
        }

        ClassicTrainingVectors ts(0, true);
        std::cout << std::endl << "Reading input data file " << fn_in << "....." << std::endl;
        inStream >> ts;

        /* Real stuff ============================================================== */
        if (norm)
        {
            std::cout << "Normalizing....." << std::endl;
            ts.normalize();        // Normalize input data
        }

        FuzzyMap *myMap;

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
                myMap = new FuzzyMap(codeStream, ts.size(), false);
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
                myMap = new FuzzyMap(codeStream, ts.size(), true);
            }
        }
        else
            myMap = new FuzzyMap(layout, xdim, ydim, ts);


        KerDenSOM *thisSOM;
        if (fn_algo_in == "")
        {
            if (gaussian)
                thisSOM = new GaussianKerDenSOM(reg0, reg1, annSteps, eps, iter);        // Creates KerDenSOM Algorithm
            else
                thisSOM = new TStudentKerDenSOM(reg0, reg1, annSteps, eps, iter, df);    // Creates KerDenSOM Algorithm
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

        TextualListener myListener;       // Define the listener class
        myListener.setVerbosity() = verbose;       // Set verbosity level
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
        if ((annSteps == 0) || (annSteps == 1))
        {
            infS << "Kernel Probability Density Estimator clustering algorithm" << std::endl;
            infS << "                 Kernel c-Means" << std::endl << std::endl;
        }
        else
        {
            infS << "Kernel Probability Density Estimator SOM algorithm" << std::endl;
            infS << "                     KerDenSOM" << std::endl << std::endl;
        }
        infS << "Input data file : " << fn_in << std::endl;
        if (cb_in != "")
            infS << "Input code vectors file : " << cb_in << std::endl;
        infS << "Code vectors output file : " << fn_out <<  ".cod" << std::endl;
        infS << "Whole codebook output file : " << fn_out <<  ".cbk" << std::endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << std::endl;
        infS << "Number of feature vectors: " << ts.size() << std::endl;
        infS << "Number of variables: " << ts.theItems[0].size() << std::endl;
        infS << "Horizontal dimension (Xdim) = " << xdim << std::endl;
        infS << "Vertical dimension (Ydim) = " << ydim << std::endl;
        if (norm)
            infS << "Input data normalized" << std::endl;
        else
            infS << "Input data not normalized" << std::endl;

        if (annSteps > 1)
        {
            if (layout == "HEXA")
                infS << "Hexagonal topology " << std::endl;
            else
                infS << "Rectangular topology " << std::endl;
        }

        if (gaussian)
            infS << "Gaussian Kernel function " << std::endl;
        else
        {
            infS << "t-Student Kernel function" << std::endl;
            infS << "Degrees of freedom (df) = " << df << std::endl;
        }

        if (annSteps > 1)
        {
            infS << "Initial smoothness factor (reg0) = " << reg0 << std::endl;
            infS << "Final smoothness factor (reg1) = " << reg1 << std::endl;
            infS << "Deterministic annealing steps = " << annSteps << std::endl;
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
};

/* Main function -============================================================= */
int main(int argc, char** argv)
{
    try
    {
        ProgKenderSOM prm;
        prm.read(argc,argv);
        prm.run();
    }
    catch (XmippError XE)
    {
        std::cout << XE << std::endl;
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}
