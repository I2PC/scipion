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

#include <mpi.h>
#include <data/filters.h>
#include <data/funcs.h>
#include <data/image.h>
#include <data/image_collection.h>
#include <data/mask.h>
#include <data/metadata.h>
#include <data/program.h>

class ProgSortImages: public XmippProgram
{
public:
    /** Filename selection file containing the images */
    FileName fnSel;

    /**  Filename output rootname */
    FileName fnRoot;
public:
    /** Filename of the output stack */
    FileName fnStack;

    // Output selfile
    MetaData SFout;

    // SelFile images
    std::vector< FileName > toClassify;

    // Control vector to specify which ones have already been done
    Matrix1D<int> stillToDo;

    // Image holding current reference
    Image<double> lastImage;

    // Mask of the background
    MultidimArray<int> mask;
public:
    /// Read argument from command line
    void readParams()
    {
        fnSel  = getParam("-i");
        fnRoot = getParam("--oroot");
    }

    /// Show
    void show()
    {
        std::cerr << "Input selfile:    " << fnSel           << std::endl
        << "Output rootname:  " << fnRoot          << std::endl
        ;
    }

    /// Usage
    void defineParams()
    {
        addUsageLine("Sort a set of images by local similarity");
        addParamsLine("   -i <selfile>        : selfile of images");
        addParamsLine("   --oroot <rootname>  : output rootname");
    }

    /// Produce side info
    void produceSideInfo(int rank)
    {
        fnStack=fnRoot+".stk";
        if (exists(fnStack) && rank==0)
            unlink(fnStack.c_str());

        // Read input selfile and reference
        ImageCollection SF;
        SF.read(fnSel);
        int idx=0;
        FileName fnImg;
        toClassify.reserve(SF.size());
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            if (idx==0)
            {
                SF.getValue(MDL_IMAGE,fnImg);
                toClassify.push_back(fnImg);
                lastImage.read(fnImg);
                centerImage(lastImage());
                if (rank==0)
                    lastImage.write(fnStack,idx,true,WRITE_APPEND);
                FileName fnImageStack;
                fnImageStack.compose(idx,fnStack);
                if (rank==0)
                {
                    SFout.addObject();
                    SFout.setValue(MDL_IMAGE,fnImageStack);
                    SFout.setValue(MDL_IMAGE_ORIGINAL,fnImg);
                    SFout.setValue(MDL_MAXCC,1.0);
                }
            }
            else
            {
                SF.getValue(MDL_IMAGE,fnImg);
                toClassify.push_back(fnImg);
            }
            idx++;
        }
        stillToDo.resizeNoCopy(SF.size());
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
        FOR_ALL_ELEMENTS_IN_MATRIX1D(stillToDo)
        {
            if (!VEC_ELEM(stillToDo,i))
                continue;
            if ((count+1)%Nproc!=rank)
            {
            	++count;
            	continue;
            }
            I.read(toClassify[i]);
            I().setXmippOrigin();
            double corr=alignImagesConsideringMirrors(lastImage(),I(),M,&mask);
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
        	bestIdx=buffer[0];
        }

        // All compute the best image
        I.read(toClassify[bestIdx]);
        I().setXmippOrigin();
        bestCorr=alignImagesConsideringMirrors(lastImage(),I(),M,&mask);
        bestImage()=I();

        lastImage()=bestImage();
        if (rank==0)
        {
            int idxStack=SFout.size();
            FileName fnImageStack;
            fnImageStack.compose(idxStack,fnStack);
            bestImage.write(fnStack,idxStack,true,WRITE_APPEND);
            SFout.addObject();
            SFout.setValue(MDL_IMAGE,fnImageStack);
            SFout.setValue(MDL_IMAGE_ORIGINAL,toClassify[bestIdx]);
            SFout.setValue(MDL_MAXCC,bestCorr);
        }
        VEC_ELEM(stillToDo,bestIdx)=0;
    }

    /// Main routine
    void run(int rank, int Nproc)
    {
        while (stillToDo.sum()>0)
            chooseNextImage(rank,Nproc);
        if (rank==0)
            SFout.write(fnRoot+".sel");
    }
};

int main(int argc, char **argv)
{
    // Read input parameters
    ProgSortImages program;
    program.read(argc, argv);

    // Start MPI communications
    MPI_Init(&argc, &argv);
    int rank, Nprocessors;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Nprocessors);

    try
    {
        if (rank==0)
            program.show();
        program.produceSideInfo(rank);
        program.run(rank,Nprocessors);
    }
    catch (XmippError XE)
    {
        std::cerr << XE << std::endl;
    }
    return 0;
}
