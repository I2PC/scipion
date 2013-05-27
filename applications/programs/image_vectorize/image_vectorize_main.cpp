/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <data/xmipp_program.h>
#include <data/mask.h>
#include <data/metadata_extension.h>

class ProgImageVectorize: public XmippProgram
{
public:
    FileName fnIn;
    FileName fnOut;
    Mask     mask;
    bool     apply_mask;
    void defineParams()
    {
        addUsageLine("Convert an image/volume to data and viceversa.");
        addUsageLine("Data is the file format used by KerDenSOM and other classification programs");
        addParamsLine("   -i <file> : If Image->Vector: Input Single image/volume/stack or selfile");
        addParamsLine("             : If Vector->Image: Input Vector metadata");
        addParamsLine("   -o <file> : If Image->Vector: Output vector metadata");
        addParamsLine("             : If Vector->Image: Output stack");
        mask.defineParams(this,INT_MASK,NULL,"Extract pixel values from a mask area.");
        addExampleLine("Produce vectors from images or volumes",false);
        addExampleLine("xmipp_image_vectorize -i projections.sel -o vectors.xmd");
        addExampleLine("Produce images from vectors",false);
        addExampleLine("xmipp_image_vectorize -i vectors.xmd -o images.stk");
        addExampleLine("Produce vectors from images or volumes with a mask",false);
        addExampleLine("xmipp_image_vectorize -i projections.sel -o vectors.xmd --mask circular -16");
        addExampleLine("Produce images from vectors that were defined using a mask",false);
        addExampleLine("xmipp_image_vectorize -i vectors.xmd -o images.stk --mask circular -16");
    }

    void readParams()
    {
        fnIn = getParam("-i");
        fnOut = getParam("-o");
        mask.allowed_data_types = INT_MASK;
        if ((apply_mask = checkParam("--mask")))
            mask.readParams(this);
    }

    void show()
    {
        if (verbose==0)
            return;
        std::cout
        << "Input:   " << fnIn    << std::endl
        << "Output:  " << fnOut << std::endl
        ;
    }


    void run()
    {
        StringVector blockList;
        getBlocksInMetaDataFile(fnIn, blockList);
        bool headerFound=false;
        for (size_t i=0; i<blockList.size(); i++)
            if (blockList[i]=="vectorHeader")
            {
                headerFound=true;
                break;
            }
        if (!headerFound)
        {
            // Image -> Vector
            MetaData SF(fnIn);
            MetaData vectorContent, vectorHeader;
            vectorHeader.setColumnFormat(false);
            Image<double> img;
            bool first=true;
            size_t order=0;
            FileName fnImg;
            float *buffer=NULL;
            FileName fnOutRaw=formatString("%s.vec",fnOut.withoutExtension().c_str());
            std::ofstream fhOutRaw(fnOutRaw.c_str(),std::ios::binary);
            if (!fhOutRaw)
                REPORT_ERROR(ERR_IO_NOWRITE,fnOutRaw);
            size_t vectorSize;
            MDRow row;
            FOR_ALL_OBJECTS_IN_METADATA(SF)
            {
                img.readApplyGeo(SF, __iter.objId);

                // Create header
                if (first)
                {
                    if (apply_mask)
                    {
                        mask.generate_mask(ZSIZE(img()),YSIZE(img()),XSIZE(img()));
                        vectorSize=(int)mask.get_binary_mask().sum();
                    }
                    else
                        vectorSize=MULTIDIM_SIZE(img());
                    size_t outId=vectorHeader.addObject();
                    vectorHeader.setValue(MDL_XSIZE,XSIZE(img()),outId);
                    vectorHeader.setValue(MDL_YSIZE,YSIZE(img()),outId);
                    vectorHeader.setValue(MDL_ZSIZE,ZSIZE(img()),outId);
                    vectorHeader.setValue(MDL_COUNT,SF.size(),outId);
                    vectorHeader.setValue(MDL_CLASSIFICATION_DATA_SIZE,vectorSize,outId);
                    vectorHeader.write(formatString("vectorHeader@%s",fnOut.c_str()));
                    buffer=new float[vectorSize];
                    first=false;
                }

                // Save this image in the output metadata
                SF.getRow(row,__iter.objId);
                row.setValue(MDL_ORDER,order++);
                vectorContent.addRow(row);

                // Save raw values
                const MultidimArray<double> &mimg=img();
                const MultidimArray<int> &mmask=mask.get_binary_mask();
                int ii=0;
                if (apply_mask)
                {
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mimg)
                    if (DIRECT_MULTIDIM_ELEM(mmask,n))
                        buffer[ii++]=(float)DIRECT_MULTIDIM_ELEM(mimg,n);
                }
                else
                {
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mimg)
                    buffer[ii++]=(float)DIRECT_MULTIDIM_ELEM(mimg,n);
                }
                fhOutRaw.write((char*)buffer,vectorSize*sizeof(float));
            }
            fhOutRaw.close();
            vectorContent.write(formatString("vectorContent@%s",fnOut.c_str()),MD_APPEND);
            delete []buffer;
        }
        else
        {
            // Vector -> Image
            fnOut.deleteFile();

            MetaData vectorHeader(formatString("vectorHeader@%s",fnIn.c_str()));
            MetaData vectorContent(formatString("vectorContent@%s",fnIn.c_str()));

            // Read header
            size_t headerId=vectorHeader.firstObject();
            size_t Xdim, Ydim, Zdim, vectorSize, imgNo;
            vectorHeader.getValue(MDL_XSIZE,Xdim,headerId);
            vectorHeader.getValue(MDL_YSIZE,Ydim,headerId);
            vectorHeader.getValue(MDL_ZSIZE,Zdim,headerId);
            vectorHeader.getValue(MDL_COUNT,imgNo,headerId);
            vectorHeader.getValue(MDL_CLASSIFICATION_DATA_SIZE,vectorSize,headerId);
            if (apply_mask)
                mask.generate_mask(Zdim,Ydim,Xdim);
            const MultidimArray<int> &mmask=mask.get_binary_mask();

            // Process all images
            Image<double> img;
            img().resizeNoCopy(Zdim,Ydim,Xdim);
            const MultidimArray<double> &mimg=img();
            FileName fnImg, fnIdx;
            float *buffer=new float[vectorSize];
            FileName fnInRaw=formatString("%s.vec",fnIn.withoutExtension().c_str());
            std::ifstream fhInRaw(fnInRaw.c_str(),std::ios::binary);
            if (!fhInRaw)
                REPORT_ERROR(ERR_IO_NOTEXIST,fnInRaw);
            size_t order;
            size_t idx=1;
            FOR_ALL_OBJECTS_IN_METADATA(vectorContent)
            {
                vectorContent.getValue(MDL_IMAGE,fnImg,__iter.objId);
                vectorContent.getValue(MDL_ORDER,order,__iter.objId);

                // Read raw values
                fhInRaw.seekg(order*vectorSize*sizeof(float));
                fhInRaw.read((char*)buffer,vectorSize*sizeof(float));
                if (!fhInRaw)
                    REPORT_ERROR(ERR_IO_NOREAD,
                                 formatString("Could not read image %lu from %s",
                                              order,fnInRaw.c_str()));
                int ii=0;
                if (apply_mask)
                {
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mimg)
                    if (DIRECT_MULTIDIM_ELEM(mmask,n))
                        DIRECT_MULTIDIM_ELEM(mimg,n)=buffer[ii++];
                }
                else
                {
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mimg)
                    DIRECT_MULTIDIM_ELEM(mimg,n)=buffer[ii++];
                }

                // Write image
                fnIdx.compose(idx++,fnOut);
                img.write(fnIdx);
            }
            fhInRaw.close();
            delete []buffer;
        }
    }
}
;

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgImageVectorize program;
    program.read(argc, argv);
    return program.tryRun();
}

