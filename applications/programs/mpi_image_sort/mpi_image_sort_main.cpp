 /***************************************************************************
 *
 * Authors:    Carlos Oscar           coss@cnb.csic.es (2010)
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

//#include <mpi.h>
#include <parallel/xmipp_mpi.h>
#include <data/filters.h>
#include <data/xmipp_image.h>
#include <data/mask.h>
#include <data/metadata.h>
#include <data/xmipp_program.h>


class ProgSortImages: public XmippProgram
{
public:
    /** Filename selection file containing the images */
    FileName fnSel;

    /**  Filename output rootname */
    FileName fnRoot;
public:
    MpiNode *node;

    /** Filename of the output stack */
    FileName fnStack;

    // Output selfile
    MetaData SFout;

    // SelFile images
    std::vector< MDRow > toClassify;

    // Control vector to specify which ones have already been done
    Matrix1D<int> stillToDo;

    // Image holding current reference
    Image<double> lastImage;

    // Mask of the background
    MultidimArray<int> mask;
public:
    /// MPI constructor
    ProgSortImages(int argc, char **argv)
    {
        node=new MpiNode(argc,argv);
        if (!node->isMaster())
            verbose=0;
    }

    /// MPI destructor
    ~ProgSortImages()
    {
        delete node;
    }

    /// Read argument from command line
    void readParams()
    {
        fnSel  = getParam("-i");
        fnRoot = getParam("--oroot");
    }

    /// Show
    void show()
    {
        if (!verbose)
            return;
        std::cerr << "Input selfile:    " << fnSel           << std::endl
        << "Output rootname:  " << fnRoot          << std::endl
        ;
    }

    /// Usage
    void defineParams()
    {
        addUsageLine("Sort a set of images by local similarity");
        addUsageLine("+The program takes the first image of the input selfile as reference.");
        addUsageLine("+Then, it looks for the image among the remaining ones that is most similar to it once they have been aligned.");
        addUsageLine("+Now, the most similar image is aligned to the reference image (the first one) and is added to the list of sorted images.");
        addUsageLine("+The second image acts now as the reference, and the most similar image among the remaining images is added to the list of sorted images.");
        addUsageLine("+This process is repeated until no image is left.");
        addUsageLine("+");
        addUsageLine("+Note that only the MPI version of this program exists");
        addParamsLine("   -i <selfile>        : selfile of images");
        addParamsLine("   --oroot <rootname>  : output rootname");
        addParamsLine("                       : rootname.stk contains the list of aligned images.");
        addParamsLine("                       : rootname.xmd contains the correspondence between aligned images ");
        addParamsLine("                       : and the original images as well as the correlation coefficient ");
        addParamsLine("                       : between the aligned image and its predecessor in the list of ");
        addParamsLine("                       : aligned images.");
        addExampleLine("mpirun -np 8 `which xmipp_mpi_image_sort` -i images.xmd --oroot sorted");
    }

    /// Produce side info
    void produceSideInfo(int rank)
    {
        fnStack = fnRoot + ".stk";

        if (rank == 0)
            fnStack.deleteFile();

        // Read input selfile and reference
        MetaData SF;
        SF.read(fnSel);
        FileName fnImg;
        toClassify.reserve(SF.size());
        FileName fnImageStack;
        CorrelationAux aux;
        RotationalCorrelationAux aux2;
        MDRow row;
        bool thereIsClassCount=SF.containsLabel(MDL_CLASS_COUNT);
        bool firstSelected=false;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            SF.getRow(row,__iter.objId);
            bool proceed=true;
            if (thereIsClassCount)
            {
                size_t classCount;
                row.getValue(MDL_CLASS_COUNT,classCount);
                if (classCount==0)
                    proceed=false;
            }
            if (proceed)
            {
                toClassify.push_back(row);
                row.getValue(MDL_IMAGE,fnImg);
                if (!firstSelected)
                {
                    row.getValue(MDL_IMAGE,fnImg);
                    lastImage.read(fnImg);
                    centerImage(lastImage(),aux,aux2);
                    if (rank==0)
                        lastImage.write(fnStack,0,true,WRITE_APPEND);
                    fnImageStack.compose(1,fnStack);
                    if (rank==0)
                    {
                        row.setValue(MDL_IMAGE_ORIGINAL,fnImg);
                        row.setValue(MDL_MAXCC,1.0);
                        SFout.addRow(row);
                    }
                    firstSelected=true;
                }
            }
        }
        if (toClassify.size()==0)
            return;
        stillToDo.resizeNoCopy(toClassify.size());
        stillToDo.initConstant(1);
        VEC_ELEM(stillToDo,0)=0;

        // Prepare mask
        mask.resize(lastImage());
        mask.setXmippOrigin();
        BinaryCircularMask(mask,XSIZE(lastImage())/2, INNER_MASK);
    }

    /// Choose next image
    void chooseNextImage(int rank, int Nproc)
    {
        Image<double> bestImage, I;
        Matrix2D<double> M;
        double bestCorr=-1;
        int bestIdx=-1;
        int count=0;
        FileName fnImageStack, fnImg;
        AlignmentAux aux;
        CorrelationAux aux2;
        RotationalCorrelationAux aux3;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(stillToDo)
        {
            if (!VEC_ELEM(stillToDo,i))
                continue;
            if ((count+1)%Nproc!=rank)
            {
                ++count;
                continue;
            }
            toClassify[i].getValue(MDL_IMAGE,fnImg);
            I.read(fnImg);
            I().setXmippOrigin();
            double corr=alignImagesConsideringMirrors(lastImage(),I(),M,aux,aux2,aux3,&mask);
            centerImage(I(),aux2,aux3);
            if (corr>bestCorr)
            {
                bestCorr=corr;
                bestImage()=I();
                bestIdx=i;
            }
            ++count;
        }

        // Rank 0 receives from the other nodes their best image
        double buffer[2];
        if (rank==0)
        {
            MPI_Status status;
            for (int n=1; n<Nproc; ++n)
            {
                MPI_Recv(buffer, 2, MPI_DOUBLE, n, 0, MPI_COMM_WORLD, &status);
                if (buffer[1]>bestCorr)
                {
                    bestIdx=(int)buffer[0];
                    bestCorr=buffer[1];
                }
            }
            std::cout << "Images to go=" << stillToDo.sum()-1 << " current correlation= " << bestCorr << std::endl;
        }
        else
        {
            buffer[0]=bestIdx;
            buffer[1]=bestCorr;
            MPI_Send(buffer, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        // Now rank 0 redistributes the best image
        if (rank==0)
        {
            buffer[0]=bestIdx;
            MPI_Bcast(buffer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Bcast(buffer, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            bestIdx=(int)buffer[0];
        }

        // All compute the best image
        toClassify[bestIdx].getValue(MDL_IMAGE,fnImg);
        I.read(fnImg);
        I().setXmippOrigin();
        bestCorr=alignImagesConsideringMirrors(lastImage(),I(),M,aux,aux2,aux3,&mask);
        centerImage(I(),aux2,aux3);
        bestImage()=I();
        lastImage()=bestImage();

        if (rank==0)
        {
            int idxStack=SFout.size();
            fnImageStack.compose(idxStack+1,fnStack);
            bestImage.write(fnStack,idxStack,true,WRITE_APPEND);
            MDRow &row=toClassify[bestIdx];
            FileName fnImgOrig;
            row.getValue(MDL_IMAGE,fnImgOrig);
            row.setValue(MDL_IMAGE_ORIGINAL,fnImgOrig);
            row.setValue(MDL_IMAGE,fnImageStack);
            row.setValue(MDL_MAXCC,bestCorr);
            SFout.addRow(row);
        }
        VEC_ELEM(stillToDo,bestIdx)=0;
    }

    /// Main routine
    void run()
    {
        show();
        produceSideInfo(node->rank);
        while (stillToDo.sum()>0)
            chooseNextImage(node->rank,node->size);
        if (node->rank==0)
            SFout.write(fnRoot+".xmd");
    }
};

int main(int argc, char **argv)
{
    // Read input parameters
    ProgSortImages program(argc, argv);
    program.read(argc, argv);
    return program.tryRun();
}
