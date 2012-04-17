/***************************************************************************
 *
 * Authors:    Javier Vargas            jvargas@cnb.csic.es
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

#include <data/xmipp_polynomials.h>
#include <data/xmipp_program.h>
#include <data/multidim_array.h>
#include <data/fringe_processing.h>
#include <data/xmipp_image.h>
#include <data/xmipp_fftw.h>

class ProgCTFEstimateFromPSDZernike: public XmippProgram
{
public:
    /// CTF filename
    FileName fn_psd;
    double lambda;
    int size;
    int x;
    int y;
    double R;
    double S;
    double thrs;

public:
    /// CTF amplitude to model
    Image<double> ctftomodel;
public:
    /// Read parameters
    void readParams()
    {
    	fn_psd=getParam("--psd");
    	lambda = getDoubleParam("--lambda");
    	R = getDoubleParam("--freq");
    	S = getDoubleParam("--var");
    	size = getIntParam("--size");
    	x = getIntParam("--x");
    	y = getIntParam("--y");
    	thrs = getDoubleParam("--thrs");

    }

    /// Show parameters
    void show()
    {
    	std::cout << "PSD file:    " << fn_psd << std::endl;
    }

    /// Define Parameters
    void defineParams()
    {
        addUsageLine("Adjust a parametric model to a PSD file.");
        addUsageLine("The PSD is enhanced ([[http://www.ncbi.nlm.nih.gov/pubmed/16987671][See article]]). ");
        addUsageLine("And finally, the CTF is fitted to the PSD, being guided by the enhanced PSD ");
        addUsageLine("([[http://www.ncbi.nlm.nih.gov/pubmed/17911028][See article]]).");
        addParamsLine("--psd <PSDfile> : PSD file");
        addParamsLine("[--lambda <v=2>] : regularization parameter");
        addParamsLine("--freq <R>: Rough estimation of the fringe frequencies");
        addParamsLine("--var <S>: variance of the fringe frequency along the pattern");
        addParamsLine("[--size <size=5>] : regularization window");
        addParamsLine("--x <x> : x coordinate of the direction starting point");
        addParamsLine("--y <y> : x coordinate of the direction starting point");
        addParamsLine("[--thrs <thrs=5>] : intensity thresholding parameter to cutoff the psd intensity above a athis value");
    }

    /// Produce side information
    void produce_side_info()
    {
    	ctftomodel.read(fn_psd);
    }

    /** Run */
    void run()
    {

    	produce_side_info();

    	FringeProcessing fp;
    	MultidimArray<double> mod, in, phase;
    	MultidimArray<double> & im =  ctftomodel();

    	mod.resizeNoCopy(im);
    	in.resizeNoCopy(im);
    	phase.resizeNoCopy(im);

    	CenterFFT(im,true);
    	im.threshold("abs_above", thrs, 0);

        fp.demodulate(im,R,S,lambda,size,x,y,phase,mod,verbose);
    }
};

int main (int argc,char *argv[])
{
	ProgCTFEstimateFromPSDZernike prog_prm;
	prog_prm.read(argc,argv);
	return prog_prm.tryRun();
}
