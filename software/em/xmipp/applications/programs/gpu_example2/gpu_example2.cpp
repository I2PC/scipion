#include <data/multidim_array.h>
#include <data/matrix2d.h>
#include <data/xmipp_fftw.h>

int main(void)
{
    MultidimArray<double> m1(2,2);
    MultidimArray<double> m2(2,2);
    MultidimArray<double> m3(2,2);


    A2D_ELEM(m1,0,0) = 1.;
    A2D_ELEM(m1,1,0) = 2.;
    A2D_ELEM(m1,0,1) = 3.;
    A2D_ELEM(m1,1,1) = 4.;

    A2D_ELEM(m2,0,0) = 11.;
    A2D_ELEM(m2,1,0) = 22.;
    A2D_ELEM(m2,0,1) = 33.;
    A2D_ELEM(m2,1,1) = 44.;

    m3 = m1 + m2;

    std::cerr << "m1" << m1 << std::endl;
    std::cerr << "m2" << m2 << std::endl;
    std::cerr << "m3" << m3 << std::endl;

    MultidimArray< std::complex< double > > FFT1, coeffs;
    FourierTransformer transformer1;

    //void FourierTransform(T& v, T1& V, bool getCopy=true)
    transformer1.FourierTransform(m2, FFT1, true);
    //void getFourierCopy(T& V)

    std::cerr << "FFT1" << FFT1 << std::endl;

    //void inverseFourierTransform(T& V, T1& v)
    MultidimArray<double> recover(2,2);
    transformer1.inverseFourierTransform(FFT1, recover);

    std::cerr << "FFT1" << FFT1 << std::endl;
    //std::cerr << "What????" << std::endl;
    std::cerr << "recover" << recover << std::endl;

    transformer1.cleanup();

/*
    for(int i=0; i< 2;i++)
       for(int j=0; j< 2;j++)
           A2D_ELEM(m3,i,j) = A2D_ELEM(m1,i,j)+A2D_ELEM(m2,i,j);
*/


    return 0;
}

