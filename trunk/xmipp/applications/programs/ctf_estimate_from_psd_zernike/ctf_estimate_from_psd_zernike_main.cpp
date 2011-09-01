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

class ProgCTFEstimateFromPSDZernike: public XmippProgram
{
public:
    /// CTF filename
    FileName             fn_psd;
public:
    /// CTF amplitude to model
    Image<double>        ctftomodel;
public:
    /// Read parameters
    void readParams()
    {
    	fn_psd=getParam("--psd");
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
        addParamsLine("   --psd <PSDfile> : PSD file");
        addSeeAlsoLine("ctf_estimate_from_micrograph, ctf_estimate_from_psd");
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
    	PolyZernikes polynom;
    	Matrix1D<int> coefs(6);
    	VEC_ELEM(coefs,0)=1;
    	VEC_ELEM(coefs,1)=1;
    	VEC_ELEM(coefs,2)=1;
    	VEC_ELEM(coefs,3)=1;
    	VEC_ELEM(coefs,4)=1;
    	VEC_ELEM(coefs,5)=1;
    	VEC_ELEM(coefs,6)=1;
    	//VEC_ELEM(coefs,7)=1;
    	//VEC_ELEM(coefs,8)=1;
    	//VEC_ELEM(coefs,9)=1;

    	//polynomials.create(coefs,ctftomodel());
    	polynom.create(coefs);
    	//polynom.fit(coefs,ctftomodel());
    }
};

int main (int argc,char *argv[])
{
	ProgCTFEstimateFromPSDZernike prog_prm;
	prog_prm.read(argc,argv);
	return prog_prm.tryRun();
}
