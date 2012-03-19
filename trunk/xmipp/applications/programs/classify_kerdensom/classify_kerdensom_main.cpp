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

#include <data/xmipp_program.h>
#include <classification/gaussian_kerdensom.h>

/* Parameters class ======================================================= */
class ProgKenderSOM: public XmippProgram
{
public:
    /* Input Parameters ======================================================== */
    FileName       fn_in;        // input file
    FileName       fn_root;      // output file
    FileName       tmpN;         // Temporary variable
    double         eps;          // Stopping criteria
    unsigned       iter;         // Iteration number
    bool           norm;         // Normalize?
    int            xdim;         // X-dimension (-->)
    int            ydim;         // Y-dimension
    double         reg0;         // Initial reg
    double         reg1;         // Final reg
    std::string    layout;       // layout (Topology)
    unsigned       annSteps;     // Deterministic Annealing steps
public:
    // Define parameters
    void defineParams()
    {
        addUsageLine("Purpose: Kernel Density Estimator Self-Organizing Map");
        addUsageLine("+KerDenSOM stands for Kernel Probability Density Estimator Self-Organizing Map.");
        addUsageLine("+It maps a set of high dimensional input vectors into a two-dimensional grid. ");
        addUsageLine("+For more information, please see the following [[http://www.ncbi.nlm.nih.gov/pubmed/11472094][reference]].");
        addUsageLine("+");
        addUsageLine("+The topology of the network can be hexagonal or rectangular (see below). ");
        addUsageLine("+It is advised to design maps with one of its sides larger than the other (e.g. 10x5).");
        addUsageLine("+   Xdim is ------>",true);
        addUsageLine("+   HEXAGONAL:",true);
        addUsageLine("+ O O O O O O O O O",true);
        addUsageLine("+O O O & & & O O O",true);
        addUsageLine("+ O O & @ @ & O O O",true);
        addUsageLine("+O O & @ + @ & O O",true);
        addUsageLine("+ O O & @ @ & O O O",true);
        addUsageLine("+O O O & & & O O O",true);
        addUsageLine("+ O O O O O O O O O",true);
        addUsageLine("+   RECTANGULAR:",true);
        addUsageLine("+O O O O O O O O O",true);
        addUsageLine("+0 O O O & O O O O",true);
        addUsageLine("+O O O & @ & O O O",true);
        addUsageLine("+O O & @ + @ & O O",true);
        addUsageLine("+O O O & @ & O O O",true);
        addUsageLine("+O O O O & O O O O",true);
        addUsageLine("+O O O O O O O O O",true);
        addSeeAlsoLine("image_vectorize");
        addParamsLine("  -i <file_in>                : Input data file");
        addParamsLine("                              : This file is generated by [[image_vectorize_v3][image_vectorize]]");
        addParamsLine(" --oroot <rootname>           : rootname_classes.xmd, rootname_images.xmd and rootname_vectors.xmd will be created");
        addParamsLine("                              : This file mst be read by [[image_vectorize_v3][image_vectorize]]");
        addParamsLine(" [--xdim <Hdimension=10>]     : Horizontal size of the map");
        addParamsLine(" [--ydim <Vdimension=5>]      : Vertical size of the map");
        addParamsLine(" [--topology <topology=RECT>] : Lattice topology");
        addParamsLine("                              :+The following picture will help in understanding the differences ");
        addParamsLine("                              :+between both topologies and the map axis convention:");
        addParamsLine("    where <topology> RECT HEXA");
        addParamsLine(" [--deterministic_annealing <steps=10> <Initial_reg=1000> <Final_reg=100>] : Deterministic annealing");
        addParamsLine(" 							 : controls the smoothness regularization");
        addParamsLine("                              : Set it to 0 0 0 for Kernel C-means");
        addParamsLine(" [--eps <epsilon=1e-7>]       : Stopping criteria");
        addParamsLine(" [--iter <N=200>]             : Number of iterations");
        addParamsLine(" [--norm]                     : Normalize input data");
        addExampleLine("xmipp_image_vectorize -i images.stk -o vectors.xmd");
        addExampleLine("xmipp_classify_kerdensom -i vectors.xmd -o kerdensom.xmd");
    }

    // Read parameters
    void readParams()
    {
        fn_in = getParam("-i");
        if (checkParam("--oroot"))
            fn_root=getParam("--oroot");
        ydim = getIntParam("--ydim");
        xdim = getIntParam("--xdim");
        layout = getParam("--topology");
        annSteps = getIntParam("--deterministic_annealing",0);
        reg0 = getDoubleParam("--deterministic_annealing",1);
        reg1 = getDoubleParam("--deterministic_annealing",2);
        eps = getDoubleParam("--eps");
        iter = getIntParam("--iter");
        norm = checkParam("--norm");

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
    }

    void show()
    {
        std::cout << "Input data file : " << fn_in << std::endl;
        std::cout << "Output rootname : " << fn_root << std::endl;
        std::cout << "Horizontal dimension (Xdim) = " << xdim << std::endl;
        std::cout << "Vertical dimension (Ydim) = " << ydim << std::endl;
        if (layout == "HEXA")
            std::cout << "Hexagonal topology " << std::endl;
        else
            std::cout << "Rectangular topology " << std::endl;
        std::cout << "Initial smoothness factor (reg0) = " << reg0 << std::endl;
        std::cout << "Final smoothness factor (reg1) = " << reg1 << std::endl;
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

        KerDenSOM *thisSOM= new GaussianKerDenSOM(reg0, reg1, annSteps, eps, iter);        // Creates KerDenSOM Algorithm

        TextualListener myListener;       // Define the listener class
        myListener.setVerbosity() = verbose;       // Set verbosity level
        thisSOM->setListener(&myListener);         // Set Listener
        FileName fnClasses=fn_root+"_classes.xmd";
        thisSOM->train(*myMap, ts, fnClasses); // Train algorithm

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
        // Save map size
        MetaData MDkerdensom;
        size_t id=MDkerdensom.addObject();
        MDkerdensom.setValue(MDL_XSIZE,xdim,id);
        MDkerdensom.setValue(MDL_YSIZE,ydim,id);
        MDkerdensom.setValue(MDL_MAPTOPOLOGY,layout,id);
        MDkerdensom.write(formatString("KerDenSOM_Layout@%s",fnClasses.c_str()),MD_APPEND);

        // save intracluster distance and number of vectors in each cluster
        MetaData MDsummary;
        for (unsigned i = 0; i < myMap->size(); i++)
        {
        	id=MDsummary.addObject();
        	MDsummary.setValue(MDL_REF,(int)i+1,id);
        	MDsummary.setValue(MDL_CLASSIFICATION_INTRACLASS_DISTANCE,myMap->aveDistances[i],id);
        	MDsummary.setValue(MDL_COUNT,(size_t)myMap->classifSizeAt(i),id);
        }
        MDsummary.write(formatString("classes@%s",fnClasses.c_str()),MD_APPEND);

        // assign data to clusters according to fuzzy threshold
        std::cout << "Saving neurons assigments ....." << std::endl;
        MetaData vectorContentIn, MDimages;
        vectorContentIn.read(formatString("vectorContent@%s",fn_in.c_str()));
        FileName fn;
        std::vector<size_t> objIds;
        vectorContentIn.findObjects(objIds);
        for (unsigned i = 0; i < myMap->size(); i++)
        {
        	MetaData MD;
            for (int j = 0; j < myMap->classifAt(i).size(); j++)
            {
            	size_t order=myMap->classifAt(i)[j];
            	vectorContentIn.getValue(MDL_IMAGE,fn,objIds[order]);
            	MD.setValue(MDL_IMAGE,fn,MD.addObject());
            	size_t id=MDimages.addObject();
            	MDimages.setValue(MDL_IMAGE,fn,id);
            	MDimages.setValue(MDL_REF,(int) i+1,id);
            }
            if (MD.size()>0)
            	MD.write(formatString("class_%06d@%s",i+1,fnClasses.c_str()),MD_APPEND);
        }
        MDimages.write(fn_root+"_images.xmd");

        // Save code vectors
        if (norm)
        {
            std::cout << "Denormalizing code vectors....." << std::endl;
            myMap->unNormalize(ts.getNormalizationInfo()); // de-normalize codevectors
        }
        MetaData vectorHeaderIn, vectorHeaderOut, vectorContentOut;
        vectorHeaderIn.read(formatString("vectorHeader@%s",fn_in.c_str()));
        vectorHeaderOut.setColumnFormat(false);
        int size, vectorSize;
        size_t idIn=vectorHeaderIn.firstObject();
        size_t idOut=vectorHeaderOut.addObject();
        vectorHeaderIn.getValue(MDL_XSIZE,size,idIn);
        vectorHeaderOut.setValue(MDL_XSIZE,size,idOut);
        vectorHeaderIn.getValue(MDL_YSIZE,size,idIn);
        vectorHeaderOut.setValue(MDL_YSIZE,size,idOut);
        vectorHeaderIn.getValue(MDL_ZSIZE,size,idIn);
        vectorHeaderOut.setValue(MDL_ZSIZE,size,idOut);
        vectorHeaderOut.setValue(MDL_COUNT,(size_t)myMap->size(),idOut);
        vectorHeaderIn.getValue(MDL_CLASSIFICATION_DATA_SIZE,vectorSize,idIn);
        vectorHeaderOut.setValue(MDL_CLASSIFICATION_DATA_SIZE,vectorSize,idOut);
        vectorHeaderOut.write(formatString("vectorHeader@%s_vectors.xmd",fn_root.c_str()),MD_APPEND);
        FileName fnVectorsRaw=fn_root+"_vectors.vec";
        std::ofstream fhVectorsRaw(fnVectorsRaw.c_str(),std::ios::binary);
        if (!fhVectorsRaw)
            REPORT_ERROR(ERR_IO_NOWRITE,fnVectorsRaw);
        for (size_t i = 0; i < myMap->size(); i++)
        {
        	id=vectorContentOut.addObject();
        	vectorContentOut.setValue(MDL_REF,(int)i+1,id);
        	vectorContentOut.setValue(MDL_ORDER,i,id);
        	fhVectorsRaw.write((char*)&(myMap->theItems[i][0]),vectorSize*sizeof(float));
        }
        fhVectorsRaw.close();
        vectorContentOut.write(formatString("vectorContent@%s_vectors.xmd",fn_root.c_str()),MD_APPEND);
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
