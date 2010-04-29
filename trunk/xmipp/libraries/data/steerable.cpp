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
#include "fftw.h"
#include "volume.h"
#include "histogram.h"
#include "filters.h"
#include "morphology.h"

// Remove wedge ------------------------------------------------------------
void MissingWedge::removeWedge(Matrix3D<double> &V) const
{
    Matrix2D<double> Epos, Eneg;
    Euler_angles2matrix(rotPos,tiltPos,0,Epos);
    Euler_angles2matrix(rotNeg,tiltNeg,0,Eneg);

    Matrix1D<double> freq(3), freqPos, freqNeg;
    Matrix1D<int> idx(3);

    XmippFftw transformer;
    Matrix3D< std::complex<double> > Vfft;
    transformer.FourierTransform(V,Vfft,false);

    FOR_ALL_ELEMENTS_IN_MATRIX3D(Vfft)
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
Steerable::Steerable(double sigma, Matrix3D<double> &Vtomograph,
                     double deltaAng, filterType type, const MissingWedge *_MW, unsigned int numThreads )
{
    MW =_MW;

    time_t start, finish;

    start = time(NULL);
    buildBasis(Vtomograph,sigma, numThreads);
    finish = time(NULL);

    printf( "buildBasis takes %ld secs. \n", finish-start);
    fflush(stdout);

    pthread_t * th_ids = (pthread_t *)malloc( numThreads * sizeof( pthread_t));
    FilterThreadArgs * th_args_filter = (FilterThreadArgs *) malloc ( numThreads * sizeof( FilterThreadArgs ) );

    start = time(NULL);
    for( int nt = 0; nt < numThreads; nt++ )
    {
        // Passing parameters to each thread
        th_args_filter[nt].myThreadID = nt;
        th_args_filter[nt].numThreads = numThreads;
        th_args_filter[nt].Vtomograph = &Vtomograph;
        th_args_filter[nt].basis = &basis;
        th_args_filter[nt].filter = type;
        th_args_filter[nt].deltaAng = deltaAng;
        pthread_create( (th_ids+nt) , NULL, filterThread, (void *)(th_args_filter+nt) );
    }

    for( int nt = 0; nt < numThreads; nt++ )
    {
        pthread_join( *(th_ids+nt), NULL );
    }

    finish = time(NULL);
    printf( "filtering takes %ld secs. \n", finish-start);
    fflush(stdout);

    free( th_ids );
    free( th_args_filter );
}

void * Steerable::filterThread( void * parameters )
{
    FilterThreadArgs * args = (FilterThreadArgs *) parameters;
    unsigned int numThreads = args->numThreads;
    unsigned int myID = args->myThreadID;
    filterType ftType = args->filter;
    Matrix3D<double> * Vtomograph = args->Vtomograph;
    std::vector< Matrix3D<double> > * basis = args->basis;
    double deltaAng = args->deltaAng;
    int index;

    // Choose a,b,c parameters as a function of the filterType
    double a;
    double b;
    double c;

    if ( ftType == FT_WALLS )
    {
        // for wall structures
        a = -(1.0/4.0);
        b = 5.0/4.0;
        c = 5.0/2.0;
    }
    else
    {
        // for filament structures
        a = 1.0;
        b = -(5.0/3.0);
        c = 10.0/3.0;
    }

    // Filter along tilt=0 and 180 => u=(1,0,0)
    double u0=1;
    double u1=0;
    double u2=0;

    int initk, endk;

    if( numThreads == 1 )
    {
        initk = 0;
        endk = ZSIZE(*Vtomograph) - 1;
    }
    else
    {
        int items_per_thread = ZSIZE(*Vtomograph) / numThreads;
        initk = myID * items_per_thread;
        endk = initk + items_per_thread - 1;
        if( myID == numThreads - 1 )
            endk = ZSIZE(*Vtomograph) - 1;
    }

    const Matrix3D<double> &basis0 = (*basis)[0];
    const Matrix3D<double> &basis1 = (*basis)[1];
    const Matrix3D<double> &basis2 = (*basis)[2];
    const Matrix3D<double> &basis3 = (*basis)[3];
    const Matrix3D<double> &basis4 = (*basis)[4];
    const Matrix3D<double> &basis5 = (*basis)[5];

    for (int k=initk; k<=endk; k++)
        for (int i=0; i< YSIZE(*Vtomograph); i++)
            for (int j=0; j< XSIZE(*Vtomograph); j++)
            {
                DIRECT_VOL_ELEM(*Vtomograph,k,i,j)= DIRECT_VOL_ELEM( basis0, k, i, j ) * (a+b) +
                                                    DIRECT_VOL_ELEM( basis1, k, i, j ) * a +
                                                    DIRECT_VOL_ELEM( basis2, k, i, j ) * a;
            }

    // Filter the rest of directions and keep the maximum

    double Ntilt = round(180.0/deltaAng);

    for (int v=1; v<Ntilt; v++)
    {
        double tilt = deltaAng*v;
        double deltaRoti = deltaAng/SIND(tilt);
        double NrotP = round(360.0/deltaRoti);

        for (int w=0; w<NrotP; w++)
        {
            double rot = w*deltaRoti;
            double u0 = SIND(rot)*COSD(tilt);
            double u1 = SIND(rot)*SIND(tilt);
            double u2 = COSD(rot);
            double a_bu0u0 = a+b*u0*u0;
            double a_bu1u1 = a+b*u1*u1;
            double a_bu2u2 = a+b*u2*u2;
            double cu0u1 = c*u0*u1;
            double cu0u2 = c*u0*u2;
            double cu1u2 = c*u1*u2;

            for (int k=initk; k<=endk; k++)
            {
                for (int i=0; i<YSIZE(*Vtomograph); i++)
                {
                    for (int j=0; j<XSIZE(*Vtomograph); j++)
                    {
                        double filterval =
                            DIRECT_VOL_ELEM( basis0, k, i, j ) * a_bu0u0 +
                            DIRECT_VOL_ELEM( basis1, k, i, j ) * a_bu1u1 +
                            DIRECT_VOL_ELEM( basis2, k, i, j ) * a_bu2u2 +
                            DIRECT_VOL_ELEM( basis3, k, i, j ) * cu0u1 +
                            DIRECT_VOL_ELEM( basis4, k, i, j ) * cu0u2 +
                            DIRECT_VOL_ELEM( basis5, k, i, j ) * cu1u2;

                        /*(*basis)[0](k,i,j) * a_bu0u0 +
                        (*basis)[1](k,i,j) * a_bu1u1 +
                        (*basis)[2](k,i,j) * a_bu2u2 +
                        (*basis)[3](k,i,j) * cu0u1 +
                        (*basis)[4](k,i,j) * cu0u2 +
                        (*basis)[5](k,i,j) * cu1u2;*/

                        if(filterval> DIRECT_VOL_ELEM(*Vtomograph,k,i,j)) // (*Vtomograph)(k,i,j))
                            DIRECT_VOL_ELEM(*Vtomograph,k,i,j)= filterval;
                        //(*Vtomograph)(k,i,j) = filterval;
                    }
                }
            }
        }
    }

    return NULL;
}

/* Build basis ------------------------------------------------------------- */
void Steerable::buildBasis(const Matrix3D<double> &Vtomograph, double sigma, unsigned int numThreads)
{
    std::vector< Matrix1D<double> > hx, hy, hz;
    generate1DFilters(sigma, Vtomograph, hx, hy, hz);

    // For this first threads operation, the number of threads is limited
    // to a maximum of 6
    int auxThreads;

    if( numThreads > 6 )
        auxThreads = 6;
    else
        auxThreads = numThreads;

    pthread_t * th_ids = (pthread_t *)malloc( auxThreads * sizeof( pthread_t));
    SteerableThreadArgs * th_args = (SteerableThreadArgs *) malloc ( auxThreads * sizeof( SteerableThreadArgs ) );

    for (int n=0; n<6; n++)
    {
        status_array.push_back( n );
    }

    basis.resize(6);

    for( int nt = 0; nt < auxThreads; nt++ )
    {
        // Passing parameters to each thread
        th_args[nt].parent = this;
        th_args[nt].myThreadID = nt;
        th_args[nt].numThreads = auxThreads;
        th_args[nt].Vin = &Vtomograph;
        th_args[nt].basis = &basis;
        th_args[nt].hx = &hx;
        th_args[nt].hy = &hy;
        th_args[nt].hz = &hz;

        pthread_create( (th_ids+nt) , NULL, processSteerableThread, (void *)(th_args+nt) );
    }

    for( int nt = 0; nt < auxThreads; nt++ )
    {
        pthread_join( *(th_ids+nt), NULL );
    }

    free( th_ids );
    free( th_args );
}

void * Steerable::processSteerableThread( void * parameters )
{
    SteerableThreadArgs * args = (SteerableThreadArgs *) parameters;
    unsigned int numThreads = args->numThreads;
    unsigned int myID = args->myThreadID;
    Steerable * parent = args->parent;
    const Matrix3D<double> * Vin = args->Vin;
    std::vector< Matrix1D<double> > * hx = args->hx;
    std::vector< Matrix1D<double> > * hy = args->hy;
    std::vector< Matrix1D<double> > * hz = args->hz;
    std::vector< Matrix3D<double> > * basis = args->basis;

    int index;

    Matrix3D<double> Vout;

    do
    {
        pthread_mutex_lock( &mutex );

        if( status_array.empty( ) )
        {
            pthread_mutex_unlock( &mutex );
            return NULL;
        }
        else
        {
            index = status_array.back( );
            status_array.pop_back( );
        }

        pthread_mutex_unlock( &mutex );

        Matrix1D< std::complex<double> > H, Aux;
        Vout.initZeros(*Vin);

        // Filter in X
#define MINUS_ONE_POWER(n) (((n)%2==0)? 1:-1)
        XmippFftw transformer;
        transformer.FourierTransform((*hx)[index],H);

        FOR_ALL_ELEMENTS_IN_MATRIX1D(H)
        H(i)*= MINUS_ONE_POWER(i);

        XmippFftw transformer2;

        Matrix1D<double> aux(XSIZE(*Vin));

        transformer2.setReal(aux);

        for (int k=0; k<ZSIZE(*Vin); k++)
            for (int i=0; i<YSIZE(*Vin); i++)
            {
                for (int j=0; j<XSIZE(*Vin); j++)
                    DIRECT_VEC_ELEM(aux,j)=DIRECT_VOL_ELEM(*Vin,k,i,j);

                transformer2.FourierTransform( );
                transformer2.getFourierAlias( Aux );
                Aux*=H;
                transformer2.inverseFourierTransform( );

                for (int j=0; j<XSIZE(*Vin); j++)
                    DIRECT_VOL_ELEM(Vout,k,i,j)=XSIZE(aux)*DIRECT_VEC_ELEM(aux,j);
            }

        // Filter in Y
        transformer.FourierTransform((*hy)[index],H);

        FOR_ALL_ELEMENTS_IN_MATRIX1D(H)
        H(i)*= MINUS_ONE_POWER(i);

        aux.initZeros(YSIZE(*Vin));
        transformer2.setReal(aux);

        for (int k=0; k<ZSIZE(*Vin); k++)
            for (int j=0; j<XSIZE(*Vin); j++)
            {
                for (int i=0; i<YSIZE(*Vin); i++)
                    DIRECT_VEC_ELEM(aux,i)=DIRECT_VOL_ELEM(Vout,k,i,j);

                transformer2.FourierTransform( );
                transformer2.getFourierAlias( Aux );
                Aux*=H;
                transformer2.inverseFourierTransform( );

                for (int i=0; i<YSIZE(*Vin); i++)
                    DIRECT_VOL_ELEM(Vout,k,i,j)=XSIZE(aux)*DIRECT_VEC_ELEM(aux,i);
            }

        // Filter in Z

        transformer.FourierTransform((*hz)[index],H);

        FOR_ALL_ELEMENTS_IN_MATRIX1D(H)
        H(i)*= MINUS_ONE_POWER(i);

        aux.initZeros(ZSIZE(*Vin));
        transformer2.setReal(aux);

        for (int i=0; i<YSIZE(*Vin); i++)
            for (int j=0; j<XSIZE(*Vin); j++)
            {
                for (int k=0; k<ZSIZE(*Vin); k++)
                    DIRECT_VEC_ELEM(aux,k)=DIRECT_VOL_ELEM(Vout,k,i,j);

                transformer2.FourierTransform( );
                transformer2.getFourierAlias( Aux );
                Aux*=H;
                transformer2.inverseFourierTransform( );

                for (int k=0; k<ZSIZE(*Vin); k++)
                    DIRECT_VOL_ELEM(Vout,k,i,j)=XSIZE(aux)*DIRECT_VEC_ELEM(aux,k);
            }

        // If Missing wedge
        if (parent->MW!=NULL)
            (parent->MW)->removeWedge(Vout);

        (*basis)[index] = Vout;

    }
    while(1);
}

/* Filter generation ------------------------------------------------------- */
void Steerable::generate1DFilters(double sigma,
                                  const Matrix3D<double> &Vtomograph,
                                  std::vector< Matrix1D<double> > &hx,
                                  std::vector< Matrix1D<double> > &hy,
                                  std::vector< Matrix1D<double> > &hz)
{

    // Initialization
    Matrix1D<double> aux;
    aux.initZeros(XSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hx.push_back(aux);

    aux.initZeros(YSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hy.push_back(aux);

    aux.initZeros(ZSIZE(Vtomograph));
    aux.setXmippOrigin();
    for (int i=0; i<6; i++) hz.push_back(aux);

    double sigma2=sigma*sigma;
    double k1 =  1.0/pow((2.0*PI*sigma),(3.0/2.0));
    double k2 = -1.0/(sigma2);

    FOR_ALL_ELEMENTS_IN_MATRIX1D(hx[0])
    {
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
        hx[0](i) = k1*k2*g*(1.0-(i2/sigma2));
        hx[1](i) = k1*k2*g;
        hx[2](i) = k1*k2*g;
        hx[3](i) = k1*k2*k2*g*i;
        hx[4](i) = k1*k2*k2*g*i;
        hx[5](i) = k1*k2*k2*g;
    }
    FOR_ALL_ELEMENTS_IN_MATRIX1D(hy[0])
    {
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
        hy[0](i) = g;
        hy[1](i) = g*(1.0-(i2/sigma2));
        hy[2](i) = g;
        hy[3](i) = g*i;
        hy[4](i) = g;
        hy[5](i) = g*i;
    }
    FOR_ALL_ELEMENTS_IN_MATRIX1D(hz[0])
    {
        double i2=i*i;
        double g = -exp(-i2/(2.0*sigma2));
        hz[0](i) = g;
        hz[1](i) = g;
        hz[2](i) = g*(1.0-(i2/sigma2));
        hz[3](i) = g;
        hz[4](i) = g*i;
        hz[5](i) = g*i;
    }
}

void Steerable::generate3DFilter(Matrix3D<double>& h3D,
                                 std::vector< Matrix1D<double> > &hx,
                                 std::vector< Matrix1D<double> > &hy,
                                 std::vector< Matrix1D<double> > &hz)
{
    h3D.initZeros(XSIZE(hz[0]),XSIZE(hy[0]),XSIZE(hx[0]));
    h3D.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX3D(h3D)
    for (int n=0; n<6; n++)
        h3D(k,i,j)+=(hz[n](k)*hy[n](i)*hx[n](j));
}

