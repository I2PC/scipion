/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

class ProgVolumeToWeb: public XmippProgram
{
public:
    /// Volume to convert
    FileName fnVol;

    /// Output central slices
    FileName fnCentralSlices;

    /// Output projections
    FileName fnProjections;

    /// Number of slices
    int Nslices;

    /// Maximum Width
    int maxWidth;

    /// Separation
    int sep;
public:
    /// Read parameters from command line
    void readParams()
    {
    	fnVol=getParam("-i");
    	if (checkParam("--central_slices"))
    	{
    		fnCentralSlices=getParam("--central_slices");
    		Nslices=getIntParam("--central_slices",1);
    	}
    	if (checkParam("--projections"))
    		fnProjections=getParam("--projections");
    	maxWidth=getIntParam("--maxWidth");
    	sep=getIntParam("--separation");
    }

    /// define parameters
    void defineParams()
    {
        addUsageLine("Creates a representation for the web of the input volume");
        addParamsLine("   -i <volume>                        : Input volume");
        addParamsLine("  [--central_slices <imgFile> <n=-1>] : Output central slices (normally jpg); n=-1 for all slices");
        addParamsLine("  [--projections <imgFile>]           : Output projections (normally jpg)");
        addParamsLine("  [--maxWidth <w=800>]                : Maximum image width");
        addParamsLine("  [--separation <s=2>]                : Separation in pixels between slices");
        addExampleLine("xmipp_volume_to_web -i volume.vol --central_slices central_slices.jpg --projections projections.jpg");
    }

    /// Run
    void run()
    {
    	Image<double> V;
    	V.read(fnVol);
    	MultidimArray<double> &mV=V();

    	Image<double> I;
    	if (!fnCentralSlices.empty())
    	{
    		// Choose the starting and finishing slice
    		if (Nslices==-1)
    			Nslices=ZSIZE(mV);
    		int kmiddle=ZSIZE(mV)/2;
    		int k0=kmiddle-Nslices/2;
    		int kF=k0+Nslices-1;

    		// Resize the output image
    		int NslicesPerRow=maxWidth/XSIZE(mV);
    		if (NslicesPerRow==0)
    			NslicesPerRow=1;
    		else if (NslicesPerRow>Nslices)
    			NslicesPerRow=Nslices;
    		int Nrows=(int)ceil(((float)Nslices)/NslicesPerRow);
    		if (Nrows==0)
    			Nrows=1;
    		I().resizeNoCopy(Nrows*YSIZE(mV)+(Nrows-1)*sep,NslicesPerRow*XSIZE(mV)+(NslicesPerRow-1)*sep);
    		I().initConstant(mV.computeMax());

    		// Copy the slices into the output image
    		int outI=0, outJ=0;
    		for (int k=k0; k<=kF; ++k)
    		{
    			int outi=outI*(YSIZE(mV)+sep);
    			int outj=outJ*(XSIZE(mV)+sep);
    			for (size_t i=0; i<YSIZE(mV); ++i, ++outi)
    				memcpy(&IMGPIXEL(I,outi,outj),&A3D_ELEM(mV,k,i,0),XSIZE(mV)*sizeof(double));
    			outJ++;
    			if (outJ==NslicesPerRow)
    			{
    				outJ=0;
    				outI++;
    			}
    		}
    		I().rangeAdjust(0,255);
    		I.write(fnCentralSlices);
    	}
    	if (fnProjections!="")
    	{
    		// Take projections
    		MultidimArray<double> pXY, pYZ, pXZ;
    		pXY.initZeros(YSIZE(mV),XSIZE(mV));
    		pXZ.initZeros(ZSIZE(mV),XSIZE(mV));
    		pYZ.initZeros(ZSIZE(mV),YSIZE(mV));
    		double maxVal=-1e38;
    		FOR_ALL_ELEMENTS_IN_ARRAY3D(mV)
    		{
    			double v=A3D_ELEM(mV,k,i,j);
    			A2D_ELEM(pXY,i,j)+=v;
    			A2D_ELEM(pXZ,k,j)+=v;
    			A2D_ELEM(pYZ,k,i)+=v;
    			maxVal=std::max(A2D_ELEM(pXY,i,j),maxVal);
    			maxVal=std::max(A2D_ELEM(pXZ,k,j),maxVal);
    			maxVal=std::max(A2D_ELEM(pYZ,k,i),maxVal);
    		}

    		// Now produce output image
    		I().resizeNoCopy(XMIPP_MAX(YSIZE(mV),ZSIZE(mV)),2*(XSIZE(mV)+sep)+YSIZE(mV));
    		I().initConstant(maxVal);
    		int outi=0, outj=0;
			for (size_t i=0; i<YSIZE(pXY); ++i, ++outi)
				memcpy(&IMGPIXEL(I,outi,outj),&A2D_ELEM(pXY,i,0),XSIZE(pXY)*sizeof(double));

			outi=0;
			outj=XSIZE(mV)+sep;
			for (size_t i=0; i<YSIZE(pXZ); ++i, ++outi)
				memcpy(&IMGPIXEL(I,outi,outj),&A2D_ELEM(pXZ,i,0),XSIZE(pXZ)*sizeof(double));

			outi=0;
			outj=2*(XSIZE(mV)+sep);
			for (size_t i=0; i<YSIZE(pYZ); ++i, ++outi)
				memcpy(&IMGPIXEL(I,outi,outj),&A2D_ELEM(pYZ,i,0),XSIZE(pYZ)*sizeof(double));

    		I().rangeAdjust(0,255);
    		I.write(fnProjections);
    	}
    }
};

int main(int argc, char **argv)
{
    ProgVolumeToWeb prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
