/***************************************************************************
 *
 * Authors:  	David Maluenda (dmaluenda@cnb.csic.es)
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

#include "coordinates_carbon_filter.h"

void ProgCoordinatesCarbonFilter::defineParams()
{
    // Parameters
    addParamsLine(" -c <coordinates> : Input coordinates");
    addParamsLine(" -m <micrograph> : Reference volume");
    addParamsLine(" [--patchSize <n=20>] : Patch size for the variance filter");
    addParamsLine(" -o <coordinates> : Ouput coordinates");
}

void ProgCoordinatesCarbonFilter::readParams()
{
	fnInCoord = getParam("-c");
	fnInMic = getParam("-m");
	fnOutCoord = getParam("-o");
	patchSize = getIntParam("--patchSize");
}

void ProgCoordinatesCarbonFilter::show()
{
    if (verbose)
		std::cout
		<< "Input Coordinates:  " << fnInCoord  << std::endl
		<< "Input Micrograph:   " << fnInMic  << std::endl
		<< "Output coordinates: " << fnOutCoord << std::endl
		<< "Patch size:         " << patchSize << std::endl
		;
}

//#define DEBUG
void ProgCoordinatesCarbonFilter::run()
{
    show();

    Micrograph mic;
    mic.open_micrograph(fnInMic);
    mic.read_coordinates(0, fnInCoord);
	mic.add_label("");

	int nPart = mic.ParticleNo();

	std::cout << "Number of particles: " << nPart << std::endl; 

	Image<double> im;
	im.read(fnInMic);

    MultidimArray<double> &matrixMic = im();
    MultidimArray<double> matrixAvg(YSIZE(matrixMic),XSIZE(matrixMic));
    matrixAvg.setXmippOrigin();
    MultidimArray<double> matrixVar = matrixAvg;
    
    std::cout << "Matrix size: " << XSIZE(matrixMic) << "x" << YSIZE(matrixMic) << std::endl;

    varianceFilter(matrixMic, matrixAvg, matrixVar, patchSize);
    matrixAvg.selfABS();

    std::cout << "Resulting matrix size: " << XSIZE(matrixMic) << "x" << YSIZE(matrixMic) << std::endl;

    // OtsuSegmentation(matrixMic);


    
    Image<double> imAvg(matrixAvg), imVar(matrixVar);

    imAvg.write("outputAvg.mrc");
    imVar.write("outputVar.mrc");
    











/*







    MultidimArray<double> patch;
    int patchSize_2=patchSize/2;
	patch.resize(patchSize,patchSize,patchSize);
    patch.setXmippOrigin();

	V.read(fnIn);
    MultidimArray<double> &mV=V();
    R.read(fnRef);
    MultidimArray<double> &mR=R();

    std::cout << "Volume size = " <<  mV.xdim << "x" << mV.ydim
    		  << "x" << mV.zdim   << std::endl  ;



	std::vector< MultidimArray<double> > patchList;
	std::vector< MultidimArray<double> > patchListR;
	int patchListLength=0;
    for (int k=patchSize_2; k<(int)ZSIZE(mV)-patchSize_2; k+=patchSize)
        for (int i=patchSize_2; i<(int)YSIZE(mV)-patchSize_2; i+=patchSize)
            for (int j=patchSize_2; j<(int)XSIZE(mV)-patchSize_2; j+=patchSize)
		        {
		        	mV.window(patch,
		            		k-patchSize_2,i-patchSize_2,j-patchSize_2,
		            	 	k+patchSize_2-1,i+patchSize_2-1,j+patchSize_2-1);
		            patchList.push_back(patch);

		            mR.window(patch,
		            		k-patchSize_2,i-patchSize_2,j-patchSize_2,
		            	 	k+patchSize_2-1,i+patchSize_2-1,j+patchSize_2-1);
		            patchListR.push_back(patch);

		            patchListLength++;
				}


	std::cout << "Vector length = " << patchListLength  << std::endl
			  << "Patch  size  = " << patchList[0].xdim
			  << "x" << patchList[0].ydim << "x" << patchList[0].zdim 
			  << std::endl ;

	CorrelationAux aux;

	MultidimArray< double> textureCorr, noiseCorr, reffCorr;
	textureCorr.resize(patchSize,patchSize,patchSize);
    textureCorr.setXmippOrigin();
	noiseCorr.resize(patchSize,patchSize,patchSize);
    noiseCorr.setXmippOrigin();
    reffCorr.resize(patchSize,patchSize,patchSize);
    reffCorr.setXmippOrigin();

    std::cout << "noiseCorr size = " <<  noiseCorr.xdim << "x" 
     		  << noiseCorr.ydim << "x" << noiseCorr.zdim << std::endl
			  << "textureCorr size = " <<  textureCorr.xdim << "x" 
     		  << textureCorr.ydim << "x" << textureCorr.zdim << std::endl  ;

    int jj;
	double cumTextureCorr=0, cumNoiseCorr=0, autoCorrNoise=0, autoCorrTexture=0,
		   crossCorrNoise=0, crossCorrTexture=0;
	for (int i=0; i<patchListLength ; i++)
		for (int j=0; j<patchListLength ; j++)
		{
			// jj is the j counter in a ciclic version (jj = jj + N)
			jj = i+j;
			if(jj>=patchListLength) jj-=patchListLength;

			correlation_matrix(patchList[i],patchList[jj],noiseCorr,aux,false);
			correlation_matrix(patchList[i],patchListR[jj],textureCorr,aux,false);
			correlation_matrix(patchListR[i],patchListR[jj],reffCorr,aux,false);
			cumTextureCorr += textureCorr.sum2();
			cumNoiseCorr += noiseCorr.sum2();
			if(i==jj)
			{
				autoCorrNoise += noiseCorr.sum2();
				autoCorrTexture += textureCorr.sum2();
			}else{
				crossCorrNoise += noiseCorr.sum2()/patchListLength;
				crossCorrTexture += textureCorr.sum2()/patchListLength;
			}
		}
	

	std::cout 
		<< "GlobalCrossCorrelation    -> " << cumTextureCorr << std::endl
		<< "GlobalAutoCorrelation     -> " << cumNoiseCorr << std::endl 
		<< "Texture/Noise GlobalCorr. -> " << cumTextureCorr/cumNoiseCorr 
		<< std::endl << "	Noise:" << std::endl
		<< "InnerCroosCorrelation     -> " << crossCorrNoise << std::endl
		<< "InnerAutoCorrelation      -> " << autoCorrNoise << std::endl 
		<< "Cross/Auto InnerCorr.     -> " << crossCorrNoise/autoCorrNoise 
		<< std::endl << "	Texture:"<< std::endl
		<< "InnerCroosCorrelation     -> " << crossCorrTexture << std::endl
		<< "InnerAutoCorrelation      -> " << autoCorrTexture << std::endl 
		<< "Cross/Auto InnerCorr.     -> " << crossCorrTexture/autoCorrTexture 
		<< std::endl << "  -->  Texture/Noise InnerCorr: "
		<< (crossCorrTexture/autoCorrTexture)/(crossCorrNoise/autoCorrNoise)
		<< std::endl ;
*/


}
#undef DEBUG
