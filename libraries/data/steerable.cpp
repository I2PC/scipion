/***************************************************************************
 *
 * Authors:     Manuel Sanchez Pau 
 *              Carlos Oscar Sanchez Sorzano
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

#include "steerable.h"
#include "xmipp_fftw.h"
#include "histogram.h"
#include "filters.h"
#include "morphology.h"

// Remove wedge ------------------------------------------------------------
void MissingWedge::removeWedge(MultidimArray<double> &V) const
{
    Matrix2D<double> Epos, Eneg;
    Euler_angles2matrix(rotPos,tiltPos,0,Epos);
    Euler_angles2matrix(rotNeg,tiltNeg,0,Eneg);
    
    Matrix1D<double> freq(3), freqPos, freqNeg;
    Matrix1D<int> idx(3);

    FourierTransformer transformer;
    MultidimArray< std::complex<double> > Vfft;
    transformer.FourierTransform(V,Vfft,false);

    FOR_ALL_ELEMENTS_IN_ARRAY3D(Vfft)
    {
        // Frequency in the coordinate system of the volume
        VECTOR_R3(idx,j,i,k);
        FFT_idx2digfreq(V,idx,freq);

        // Frequency in the coordinate system of the plane
        freqPos=Epos*freq;
        freqNeg=Eneg*freq;
        if (ZZ(freqPos)<0 || ZZ(freqNeg)>0)
            Vfft(k,i,j)=0;
    }
    transformer.inverseFourierTransform();
}

// Constructor -------------------------------------------------------------
Steerable::Steerable(double sigma, MultidimArray<double> &Vtomograph,
    double deltaAng, const std::string &filterType, const MissingWedge *_MW) 
{
    MW=_MW;
    buildBasis(Vtomograph,sigma);
 
    // Choose a,b,c parameters as a function of the filterType
    double a; 
    double b; 
    double c;
    if (filterType == "wall") {
        // for wall structures
        a = -(1.0/4.0);
        b = 5.0/4.0;
        c = 5.0/2.0;
    }
    else{      
        // for filament structures
        a = 1.0;
        b = -(5.0/3.0);
        c = 10.0/3.0;
    }

    // Filter along tilt=0 and 180 => u=(1,0,0)
    double u0=1;
    double u1=0;
    double u2=0;                
    FOR_ALL_ELEMENTS_IN_ARRAY3D(Vtomograph){
        Vtomograph(k,i,j) = basis[0](k,i,j) * (a+b*u0*u0) +
                            basis[1](k,i,j) * (a+b*u1*u1) +
                            basis[2](k,i,j) * (a+b*u2*u2) +
                         c*(basis[3](k,i,j) * u0*u1 + 
                            basis[4](k,i,j) * u0*u2 + 
                            basis[5](k,i,j) * u1*u2);
    }

    // Filter the rest of directions and keep the maximum
    double Ntilt = round(180.0/deltaAng);  
    for (int i=1;i<Ntilt;i++){
        double tilt = deltaAng*i;
        double deltaRoti = deltaAng/SIND(tilt);
        double NrotP = round(360.0/deltaRoti);        
        for (int j=0;j<NrotP;j++){
            double rot = j*deltaRoti;
            double u0 = SIND(rot)*COSD(tilt);
            double u1 = SIND(rot)*SIND(tilt);
            double u2 = COSD(rot);
            FOR_ALL_ELEMENTS_IN_ARRAY3D(Vtomograph)
            {
                double filterval =
                    basis[0](k,i,j) * (a+b*u0*u0) +
                    basis[1](k,i,j) * (a+b*u1*u1) +
                    basis[2](k,i,j) * (a+b*u2*u2) +
                 c*(basis[3](k,i,j) * u0*u1 + 
                    basis[4](k,i,j) * u0*u2 + 
                    basis[5](k,i,j) * u1*u2);

                if(filterval>Vtomograph(k,i,j))
                    Vtomograph(k,i,j) = filterval;
            }
        }
    }
}

/* Build basis ------------------------------------------------------------- */
void Steerable::buildBasis(const MultidimArray<double> &Vtomograph, double sigma)
{
    std::vector< MultidimArray<double> > hx1, hy1, hz1;
    generate1DFilters(sigma, Vtomograph, hx1, hy1, hz1);
    for (int n=0; n<6; n++)
    {
        MultidimArray<double> aux;
        singleFilter(Vtomograph,hx1[n],hy1[n],hz1[n],aux);
        basis.push_back(aux);
    }    
}        

void Steerable::singleFilter(const MultidimArray<double>& Vin,
    MultidimArray<double> &hx1, MultidimArray<double> &hy1, MultidimArray<double> &hz1,
    MultidimArray<double> &Vout){

    MultidimArray< std::complex<double> > H, Aux;
    Vout.initZeros(Vin);

    // Filter in X
    #define MINUS_ONE_POWER(n) (((n)%2==0)? 1:-1)
    FourierTransformer transformer;
    transformer.FourierTransform(hx1,H);
    
    FOR_ALL_ELEMENTS_IN_ARRAY1D(H)
          H(i)*= MINUS_ONE_POWER(i);

    FourierTransformer transformer2;
    
    MultidimArray<double> aux(XSIZE(Vin));
        	   
    transformer2.setReal(aux);		   
		   
    for (size_t k=0; k<ZSIZE(Vin); k++)
        for (size_t i=0; i<YSIZE(Vin); i++)
        {
            for (size_t j=0; j<XSIZE(Vin); j++)
                DIRECT_A1D_ELEM(aux,j)=DIRECT_A3D_ELEM(Vin,k,i,j);
			    
	    transformer2.FourierTransform( );	    
	    transformer2.getFourierAlias( Aux );
	    Aux*=H;
	    transformer2.inverseFourierTransform( );
            	    
	    for (size_t j=0; j<XSIZE(Vin); j++)
                DIRECT_A3D_ELEM(Vout,k,i,j)=XSIZE(aux)*DIRECT_A1D_ELEM(aux,j);
        }

    // Filter in Y
    transformer.FourierTransform(hy1,H);
    
    FOR_ALL_ELEMENTS_IN_ARRAY1D(H)
          H(i)*= MINUS_ONE_POWER(i);

    aux.initZeros(YSIZE(Vin));
    transformer2.setReal(aux);		   
    
    for (size_t k=0; k<ZSIZE(Vin); k++)
        for (size_t j=0; j<XSIZE(Vin); j++)
        {
            for (size_t i=0; i<YSIZE(Vin); i++)
                DIRECT_A1D_ELEM(aux,i)=DIRECT_A3D_ELEM(Vout,k,i,j);

	    transformer2.FourierTransform( );	    
	    transformer2.getFourierAlias( Aux );
	    Aux*=H;
	    transformer2.inverseFourierTransform( );
            
	    for (size_t i=0; i<YSIZE(Vin); i++)
                DIRECT_A3D_ELEM(Vout,k,i,j)=XSIZE(aux)*DIRECT_A1D_ELEM(aux,i);
        }

    // Filter in Z

    transformer.FourierTransform(hz1,H);

    FOR_ALL_ELEMENTS_IN_ARRAY1D(H)
          H(i)*= MINUS_ONE_POWER(i);

    aux.initZeros(ZSIZE(Vin));    
    transformer2.setReal(aux);		   

    for (size_t i=0; i<YSIZE(Vin); i++)
        for (size_t j=0; j<XSIZE(Vin); j++)
        {
            for (size_t k=0; k<ZSIZE(Vin); k++)
                DIRECT_A1D_ELEM(aux,k)=DIRECT_A3D_ELEM(Vout,k,i,j);

	    transformer2.FourierTransform( );	    
	    transformer2.getFourierAlias( Aux );
	    Aux*=H;
	    transformer2.inverseFourierTransform( );

            for (size_t k=0; k<ZSIZE(Vin); k++)
                DIRECT_A3D_ELEM(Vout,k,i,j)=XSIZE(aux)*DIRECT_A1D_ELEM(aux,k);
        }
    
    // If Missing wedge
    if (MW!=NULL)
        MW->removeWedge(Vout);
}

/* Filter generation ------------------------------------------------------- */
void Steerable::generate1DFilters(double sigma,
    const MultidimArray<double> &Vtomograph,
    std::vector< MultidimArray<double> > &hx1,
    std::vector< MultidimArray<double> > &hy1,
    std::vector< MultidimArray<double> > &hz1){

    // Initialization 
    MultidimArray<double> aux;
    aux.initZeros(XSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hx1.push_back(aux);
    
    aux.initZeros(YSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hy1.push_back(aux);

    aux.initZeros(ZSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hz1.push_back(aux);

    double sigma2=sigma*sigma;       
    double k1 =  1.0/pow((2.0*PI*sigma),(3.0/2.0));
    double k2 = -1.0/(sigma2);
    
    FOR_ALL_ELEMENTS_IN_ARRAY1D(hx1[0])
    {        
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
	hx1[0](i) = k1*k2*g*(1.0-(i2/sigma2));
	hx1[1](i) = k1*k2*g;
	hx1[2](i) = k1*k2*g;
	hx1[3](i) = k1*k2*k2*g*i;
	hx1[4](i) = k1*k2*k2*g*i;
	hx1[5](i) = k1*k2*k2*g;
    }    
    FOR_ALL_ELEMENTS_IN_ARRAY1D(hy1[0])
    {
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
        hy1[0](i) = g;
        hy1[1](i) = g*(1.0-(i2/sigma2));
        hy1[2](i) = g;
        hy1[3](i) = g*i;
        hy1[4](i) = g;
        hy1[5](i) = g*i;
    }
    FOR_ALL_ELEMENTS_IN_ARRAY1D(hz1[0])
    {
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
	hz1[0](i) = g;
	hz1[1](i) = g;
	hz1[2](i) = g*(1.0-(i2/sigma2));
	hz1[3](i) = g;
	hz1[4](i) = g*i;
	hz1[5](i) = g*i;
    }
}

void Steerable::generate3DFilter(MultidimArray<double>& h3D,
    std::vector< MultidimArray<double> > &hx1,
    std::vector< MultidimArray<double> > &hy1,
    std::vector< MultidimArray<double> > &hz1)
{
    h3D.initZeros(XSIZE(hz1[0]),XSIZE(hy1[0]),XSIZE(hx1[0]));
    h3D.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY3D(h3D)
        for (int n=0; n<6; n++)
            h3D(k,i,j)+=(hz1[n](k)*hy1[n](i)*hx1[n](j));    
}

