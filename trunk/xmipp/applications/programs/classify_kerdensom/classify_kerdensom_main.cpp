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
    bool           gaussian;     // Gaussian Kernel
    bool           tStudent;     // tStudent Kernel
public:
    // Define parameters
    void defineParams()
    {
        addUsageLine("Purpose: Kernel Density Estimator Self-Organizing Map");
        addParamsLine("  -i <file_in>                : Input data file (plain data)");
        addParamsLine(" [-o <file_out>]              : Base name for output data files");
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
        eps = getDoubleParam("--eps");
        iter = getIntParam("--iter");
        norm = checkParam("--norm");
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
        ClassicTrainingVectors ts(0, true);
        std::cout << std::endl << "Reading input data file " << fn_in << "....." << std::endl;
        ts.read(fn_in);

        /* Real stuff ============================================================== */
        if (norm)
        {
            std::cout << "Normalizing....." << std::endl;
            ts.normalize();        // Normalize input data
        }

        FuzzyMap *myMap = new FuzzyMap(layout, xdim, ydim, ts);

        KerDenSOM *thisSOM;
        if (gaussian)
            thisSOM = new GaussianKerDenSOM(reg0, reg1, annSteps, eps, iter);        // Creates KerDenSOM Algorithm
        else
            thisSOM = new TStudentKerDenSOM(reg0, reg1, annSteps, eps, iter, df);    // Creates KerDenSOM Algorithm

        TextualListener myListener;       // Define the listener class
        myListener.setVerbosity() = verbose;       // Set verbosity level
        thisSOM->setListener(&myListener);         // Set Listener
        thisSOM->train(*myMap, ts, fn_out);        // Train algorithm

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
        // assign data to clusters according to fuzzy threshold
        std::cout << "Saving neurons assigments ....." << std::endl;
        for (unsigned i = 0; i < myMap->size(); i++)
        {
            tmpN = fn_out.c_str() + (std::string) "."  + integerToString(i);
            std::ofstream cStream(tmpN.c_str());
            for (int j = 0; j < myMap->classifAt(i).size(); j++)
                cStream << myMap->classifAt(i)[j] << std::endl;
            cStream.flush();
        }

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
    ProgKenderSOM prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
