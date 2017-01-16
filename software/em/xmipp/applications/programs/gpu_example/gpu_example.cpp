#include <data/metadata_extension.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>

class ProgGPUExample: public XmippProgram
{
public:
    /** Filename selection file containing the images */
    FileName fnSel;
    /** Filename reference image */
    FileName fnRef;
    /**  Filename output root */
    FileName fnRoot;
public:
    /// Read argument
    void readParams()
    {
        fnSel = getParam("-i");
        fnRef = getParam("--ref");
        fnRoot = getParam("--oroot");
    }

    /// Show
    void show()
    {
        if (verbose==0)
            return;
        std::cerr
        << "Input selfile:        " << fnSel       << std::endl
        << "Input reference:      " << fnRef       << std::endl
        << "Output rootname:       " << fnRoot      << std::endl
        ;
    }

    /// Define parameters
    void defineParams()
    {
        addUsageLine("Aligns a set of images");
        addParamsLine("  -i <selfile>             : Selfile containing images to be aligned");
        addParamsLine("  --oroot <rootname>       : Output rootname");
        addParamsLine(" [--ref <image=\"\">]      : reference image; if none: pyramidal combination of subset of images");
    }

    /// Main routine
    void run()
    {
        Image<double> Iref, I;
        Iref.read(fnRef);
        
        MetaData md;
        md.read(fnSel);
        size_t Xdim, Ydim, Zdim, Ndim;
        getImageSize(md,Xdim,Ydim,Zdim,Ndim);
            
        FileName fnImg;
        FileName fnOutputStack=fnRoot+".stk";
        createEmptyFile(fnOutputStack,Xdim,Ydim,1,md.size(),true,WRITE_REPLACE);
        size_t imageNo=1;
        FOR_ALL_OBJECTS_IN_METADATA(md)
        {
            md.getValue(MDL_IMAGE,fnImg,__iter.objId);
            I.read(fnImg);
            
            I.write(fnOutputStack,imageNo,true,WRITE_REPLACE);
            imageNo+=1;
        }
    }
};

RUN_XMIPP_PROGRAM(ProgGPUExample)
