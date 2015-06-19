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

#include "fourier_filter.h"

#include <data/args.h>
#include <data/xmipp_image.h>
#include <data/mask.h>
#include <data/xmipp_fft.h>

/* Clear ------------------------------------------------------------------- */
void FourierFilter::init()
{
    FilterShape = RAISED_COSINE;
    FilterBand = LOWPASS;
    w2 = w1 = 0;
    raised_w = 0;
    ctf.clear();
    ctf.enable_CTFnoise = false;
    do_correct_phase = false;
    do_generate_3dmask = false;
    sampling_rate = -1.;
}

/* Empty constructor ---------------------------------------------------------*/
FourierFilter::FourierFilter()
{
    init();
}

/* Define params ------------------------------------------------------------------- */
void FourierFilter::defineParams(XmippProgram *program)
{
    program->addParamsLine("== Fourier ==");
    program->addParamsLine("  [ --fourier <filter_type>]    : Filter in Fourier space");
    program->addParamsLine("         where <filter_type>");
    program->addParamsLine("            low_pass  <w1> <raisedw=0.02>      : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            high_pass <w1> <raisedw=0.02>      : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            band_pass <w1> <w2> <raisedw=0.02> : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            stop_band <w1> <w2> <raisedw=0.02> : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            stop_lowbandx <w1> <raisedw=0.02>  : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            stop_lowbandy <w1> <raisedw=0.02>  : Cutoff freq (<1/2 or A)");
    program->addParamsLine("            wedge <th0> <thF> <rot=0> <tilt=0> <psi=0>  : Missing wedge (along y) for data between th0-thF ");
    program->addParamsLine("                                             : y is rotated by euler angles degrees");
    program->addParamsLine("            cone <th0>                       : Missing cone for tilt angles up to th0 ");
    program->addParamsLine("                                             : do not use mask type for wedge or cone filters");
    program->addParamsLine("            real_gaussian <w1>               : Gaussian in real space with sigma = w1");
    program->addParamsLine("            gaussian <w1>                    : Gaussian in Fourier space with sigma = w1");
    program->addParamsLine("            sparsify <p=0.975>               : Delete p percent of the smallest Fourier coefficients");
    program->addParamsLine("            ctf <ctfile>                     : Provide a .ctfparam file");
    program->addParamsLine("            ctfpos <ctfile>                  : Provide a .ctfparam file");
    program->addParamsLine("                                             : The CTF phase will be corrected before applying");
    program->addParamsLine("            bfactor <B>                      : Exponential filter (positive values for decay) ");
    program->addParamsLine("               requires --sampling;                                                         ");
    program->addParamsLine("            fsc <metadata>                   : Filter with the FSC profile contained in the metadata");
    program->addParamsLine("               requires --sampling;                                                         ");
    program->addParamsLine("            binary_file <file>               : Binary file with the filter");
    program->addParamsLine("                                             :+The filter must be defined in the whole Fourier space (not only the nonredundant part).");
    program->addParamsLine("                                             :+This filter is produced, for instance, by xmipp_ctf_group");
    program->addParamsLine("            astigmatism <sigma=120>          : Filter images according to the astigmatism of their CTF");
    program->addParamsLine("                                             :+sigma is the standard deviation of a Gaussian in degrees for the CTF phase");
    program->addParamsLine("               requires --sampling;                                                         ");
    program->addParamsLine("         alias -f;");
    program->addParamsLine("  [--sampling <sampling_rate>]               : If provided pass frequencies are taken in Ang and for the CTF case");
    program->addParamsLine("                                             : and for the CTF case this is the sampling used to compute the CTF");
    program->addParamsLine("         alias -s;");
    program->addParamsLine("         requires --fourier;");
    program->addParamsLine("  [--save <filename=\"\"> ]                  : Do not apply just save the mask");
    program->addParamsLine("         requires --fourier;");
}

/* Read parameters from command line. -------------------------------------- */
void FourierFilter::readParams(XmippProgram *program)
{
    init();
    if (program->checkParam("--save"))
        maskFn = program->getParam("--save");
    if (program->checkParam("--sampling"))
        sampling_rate = program->getDoubleParam("--sampling");

    // Filter shape .........................................................
    String filter_type;
    filter_type = program->getParam("--fourier");

    // Read frequencies cuttoff options
    if (filter_type == "low_pass")
    {
        w1 = program->getDoubleParam("--fourier", "low_pass");
        raised_w = program->getDoubleParam("--fourier", "low_pass",1);
        FilterBand = LOWPASS;
        FilterShape = RAISED_COSINE;
    }
    else if (filter_type == "high_pass")
    {
        w1 = program->getDoubleParam("--fourier", "high_pass");
        raised_w = program->getDoubleParam("--fourier", "high_pass",1);
        FilterBand = HIGHPASS;
        FilterShape = RAISED_COSINE;
    }
    else if (filter_type == "band_pass")
    {
        w1 = program->getDoubleParam("--fourier", "band_pass");
        w2 = program->getDoubleParam("--fourier", "band_pass", 1);
        raised_w = program->getDoubleParam("--fourier", "band_pass",2);
        FilterBand = BANDPASS;
        FilterShape = RAISED_COSINE;
    }
    else if (filter_type == "stop_band")
    {
        w1 = program->getDoubleParam("--fourier", "stop_band");
        w2 = program->getDoubleParam("--fourier", "stop_band", 1);
        raised_w = program->getDoubleParam("--fourier", "stop_band",2);
        FilterBand = STOPBAND;
        FilterShape = RAISED_COSINE;
    }
    else if (filter_type == "stop_lowbandx")
    {
        w1 = program->getDoubleParam("--fourier", "stop_lowbandx");
        raised_w = program->getDoubleParam("--fourier", "stop_lowbandx",1);
        FilterBand = STOPLOWBANDX;
        FilterShape = RAISED_COSINE;
    }
    else if (filter_type == "stop_lowbandy")
    {
        w1 = program->getDoubleParam("--fourier", "stop_lowbandy");
        raised_w = program->getDoubleParam("--fourier", "stop_lowbandy",1);
        FilterBand = STOPLOWBANDY;
        FilterShape = RAISED_COSINE;
    }
    else if (filter_type == "wedge")
    {
        t1 = program->getDoubleParam("--fourier", "wedge", 0);
        t2 = program->getDoubleParam("--fourier", "wedge", 1);
        rot  = program->getDoubleParam("--fourier", "wedge", 2);
        tilt = program->getDoubleParam("--fourier", "wedge", 3);
        psi  = program->getDoubleParam("--fourier", "wedge", 4);
        FilterShape = FilterBand = WEDGE;
    }
    else if (filter_type == "cone")
    {
        t1 = program->getDoubleParam("--fourier", "cone", 0);
        FilterShape = FilterBand = CONE;
    }
    else if (filter_type == "gaussian")
    {
        w1 = program->getDoubleParam("--fourier", "gaussian");
        FilterShape = GAUSSIAN;
        FilterBand = LOWPASS;
    }
    else if (filter_type == "sparsify")
    {
        percentage = program->getDoubleParam("--fourier", "sparsify");
        FilterShape = SPARSIFY;
        FilterBand = SPARSIFY;
    }
    else if (filter_type == "real_gaussian")
    {
        w1 = program->getDoubleParam("--fourier", "real_gaussian");
        FilterShape = REALGAUSSIAN;
        FilterBand = LOWPASS;
    }
    else if (filter_type == "ctf")
    {
        FilterShape = FilterBand = CTF;
        ctf.enable_CTFnoise = false;
        ctf.read( program->getParam("--fourier", "ctf") );
        if (sampling_rate > 0.)
        {
            std::cerr << "CTF was optaned at sampling rate: " << ctf.Tm << std::endl;
            std::cerr << "Image sampling rate is: " << sampling_rate << std::endl;
            ctf.Tm=sampling_rate;
        }
        ctf.produceSideInfo();
    }
    else if (filter_type == "ctfpos")
    {
        FilterShape = FilterBand = CTFPOS;
        ctf.enable_CTFnoise = false;
        ctf.read( program->getParam("--fourier", "ctfpos") );
        ctf.produceSideInfo();
        do_correct_phase = true;
    }
    else if (filter_type == "bfactor")
    {
        FilterShape = FilterBand = BFACTOR;
        w1 = program->getDoubleParam("--fourier", "bfactor");
    }
    else if (filter_type == "fsc")
    {
        FilterShape = FilterBand = FSCPROFILE;
        fnFSC = program->getParam("--fourier", "fsc");
    }
    else if (filter_type == "binary_file")
    {
        FilterShape = FilterBand = BINARYFILE;
        fnFilter = program->getParam("--fourier", "binary_file");
    }
    else if (filter_type == "astigmatism")
    {
        FilterShape = FilterBand = ASTIGMATISMPROFILE;
        w1 = program->getDoubleParam("--fourier", 1);
        w1 = DEG2RAD(w1);
    }
    else
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE, "This couldn't happen, check argument parser or params definition");

    if (!(FilterBand == BFACTOR || FilterBand==ASTIGMATISMPROFILE) && sampling_rate > 0.)
    {
        if (w1 != 0)
            w1 = sampling_rate / w1;
        if (w2 != 0)
            w2 = sampling_rate / w2;
    }
}

/* Show -------------------------------------------------------------------- */
void FourierFilter::show()
{
    {
        std::cout << "Filter Band: ";
        switch (FilterBand)
        {
        case LOWPASS:
            std::cout << "Lowpass before " << w1 << std::endl;
            break;
        case HIGHPASS:
            std::cout << "Highpass after " << w1 << std::endl;
            break;
        case BANDPASS:
            std::cout << "Bandpass between " << w1 << " and " << w2 << std::endl;
            break;
        case STOPBAND:
            std::cout << "Stopband between " << w1 << " and " << w2 << std::endl;
            break;
        case STOPLOWBANDX:
            std::cout << "Stop low band X up to " << w1 << std::endl;
            break;
        case STOPLOWBANDY:
            std::cout << "Stop low band Y up to " << w1 << std::endl;
            break;
        case CTF:
            std::cout << "CTF\n";
            break;
        case CTFPOS:
            std::cout << "CTFPOS\n";
            break;
        case BFACTOR:
            std::cout << "Bfactor "<< w1 << std::endl
                      << "Sampling rate " << sampling_rate << std::endl;
            break;
        case WEDGE:
            std::cout << "Missing wedge keep data between tilting angles of " << t1 << " and " << t2 << " deg\n";
            break;
        case CONE:
            std::cout << "Missing cone for RCT data with tilting angles up to " << t1 << " deg\n";
            break;
        case FSCPROFILE:
            std::cout << "FSC file " << fnFSC << std::endl
                      << "Sampling rate " << sampling_rate << std::endl;
            break;
        case BINARYFILE:
            std::cout << "Filter file " << fnFilter << std::endl;
            break;
        case ASTIGMATISMPROFILE:
            std::cout << "Astigmatism sigma " << RAD2DEG(w1) << std::endl;
            break;
        }
        if(FilterShape!=CONE && FilterShape!=WEDGE)
            std::cout << "Filter Shape: ";
        switch (FilterShape)
        {
        case RAISED_COSINE:
            std::cout << "Raised cosine with " << raised_w
            << " raised frequencies\n";
            break;
        case GAUSSIAN:
            std::cout << "Gaussian\n";
            break;
        case SPARSIFY:
            std::cout << "Sparsify\n";
            break;
        case REALGAUSSIAN:
            std::cout << "Real Gaussian\n";
            break;
        case CTF:
            std::cout << "CTF\n" << ctf;
            break;
        case CTFPOS:
            std::cout << "CTFPOS\n" << ctf;
            break;
        }
        if (maskFn != "")
            std::cout << "Save mask in file: " << maskFn
            << " disable actual masking"<< std::endl;
    }
}

void FourierFilter::apply(MultidimArray<double> &img)
{
    static bool firstTime = true;
    do_generate_3dmask = (img.zdim > 1);
    if (firstTime)
    {
        generateMask(img);
        firstTime = false;
    }
    if (maskFn != "")
        if (do_generate_3dmask==0 && MULTIDIM_SIZE(img)>1024*1024)
            REPORT_ERROR(ERR_IO_SIZE,"Cannot save 2D mask with xdim*ydim  > 1M");
        else
        {
            Image<int> I;
            Image<double> D;
            if ( XSIZE(maskFourier) !=0 )
            {
                I()=maskFourier;
                I.write(maskFn);
            }
            else if (XSIZE(maskFourierd)!=0)
            {
                D()=maskFourierd;
                D.write(maskFn);
            }
            else
                REPORT_ERROR(ERR_MATRIX_EMPTY,"Cannot save  mask file");

        }
    else
        applyMaskSpace(img);
}

/* Get mask value ---------------------------------------------------------- */
double FourierFilter::maskValue(const Matrix1D<double> &w)
{
    double absw = w.module();
    double wx=fabs(XX(w));
    double wy=fabs(YY(w));

    // Generate mask
    switch (FilterBand)
    {
    case LOWPASS:
        switch (FilterShape)
        {
        case RAISED_COSINE:
            if (absw<w1)
                return 1;
            else if (absw<w1+raised_w)
                return (1+cos(PI/raised_w*(absw-w1)))/2;
            else
                return 0;
            break;
        case GAUSSIAN:
            return 1/sqrt(2*PI*w1)*exp(-0.5*absw*absw/(w1*w1));
            break;
        case REALGAUSSIAN:
            return exp(-PI*PI*absw*absw*w1*w1);
            break;
        }
        break;
    case HIGHPASS:
        switch (FilterShape)
        {
        case RAISED_COSINE:
            if (absw>w1)
                return 1;
            else if (absw>w1-raised_w)
                return (1+cos(PI/raised_w*(w1-absw)))/2;
            else
                return 0;
            break;
        }
        break;
    case BANDPASS:
        switch (FilterShape)
        {
        case RAISED_COSINE:
            if (absw>=w1 && absw<=w2)
                return 1;
            else if (absw>w1-raised_w && absw<w1)
                return (1+cos(PI/raised_w*(w1-absw)))/2;
            else if (absw<w2+raised_w && absw>w2)
                return (1+cos(PI/raised_w*(w2-absw)))/2;
            else
                return 0;
            break;
        }
        break;
    case STOPBAND:
        switch (FilterShape)
        {
        case RAISED_COSINE:
            if (absw>=w1 && absw<=w2)
                return 0;
            else if (absw>w1-raised_w && absw<w1)
                return 1-(1+cos(PI/raised_w*(w1-absw)))/2;
            else if (absw<w2+raised_w && absw>w2)
                return 1-(1+cos(PI/raised_w*(w2-absw)))/2;
            else
                return 1;
            break;
        }
        break;
	case STOPLOWBANDX:
		switch (FilterShape)
		{
		case RAISED_COSINE:
			if (wx<w1)
				return 0;
			else if (wx<w1+raised_w)
				return (1+cos(PI/raised_w*(w1-wx)))/2;
			else
				return 1.0;
			break;
		}
		break;
	case STOPLOWBANDY:
		switch (FilterShape)
		{
		case RAISED_COSINE:
			if (wy<w1)
				return 0;
			else if (wy<w1+raised_w)
				return (1+cos(PI/raised_w*(w1-wy)))/2;
			else
				return 1.0;
			break;
		}
		break;
    case CTF:
        ctf.precomputeValues(XX(w)/ctf.Tm,YY(w)/ctf.Tm);
        return ctf.getValueAt();
        break;
    case CTFPOS:
        ctf.precomputeValues(XX(w)/ctf.Tm,YY(w)/ctf.Tm);
        return ABS(ctf.getValueAt());
        break;
    case BFACTOR:
        {
            double R = absw / sampling_rate;
            return exp( - (w1 / 4.)  * R * R);
        }
        break;
    case FSCPROFILE:
		{
			double R = absw / sampling_rate;
			int idx=-1;
			size_t imax=freqContFSC.size();
			double *ptrFreq=&freqContFSC[0];
			for (size_t i=0; i<imax; ++i)
				if (ptrFreq[i]>R)
				{
					idx=(int)i;
					break;
				}
			if (idx>=1)
			{
				double x0=freqContFSC[idx-1];
				double x1=freqContFSC[idx];
				double y0=FSC[idx-1];
				double y1=FSC[idx];
				return y0+(y1-y0)*(R-x0)/(x1-x0);
			}
			return 0;
		}
    case ASTIGMATISMPROFILE:
    	{
			double R = absw / sampling_rate;
	        ctf.precomputeValues(R,0.0);
	        double phaseU=ctf.getPhaseAt();
	        ctf.precomputeValues(0.0,R);
	        double phaseV=ctf.getPhaseAt();
	        double diff=(phaseV-phaseU)*0.5;
	        return exp(-0.5*diff*diff/(w1*w1));
    	}
        break;
    default:
        REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown mask type");
    }
    return 0;
}

/* Generate mask ----------------------------------------------------------- */
void FourierFilter::generateMask(MultidimArray<double> &v)
{
    if (FilterShape==SPARSIFY)
        return;

    if (FilterShape==FSCPROFILE)
    {
    	MetaData mdFSC(fnFSC);
        mdFSC.getColumnValues(MDL_RESOLUTION_FREQ,freqContFSC);
        mdFSC.getColumnValues(MDL_RESOLUTION_FRC,FSC);
    }
    // std::cout << "Generating mask " << do_generate_3dmask << std::endl;

    if (do_generate_3dmask)
    {
        transformer.setReal(v);
        MultidimArray< std::complex<double> > Fourier;
        transformer.getFourierAlias(Fourier);

        if (FilterShape==WEDGE || FilterShape==CONE)
        {
            maskFourier.initZeros(Fourier);
            maskFourier.setXmippOrigin();
            switch (FilterShape)
            {
            case WEDGE:
                {
                    Matrix2D<double> A;
                    Euler_angles2matrix(rot,tilt,psi,A,false);
                    BinaryWedgeMask(maskFourier, t1, t2, A,true);
                    break;
                }
            case CONE:
                BinaryConeMask(maskFourier, 90. - fabs(t1),INNER_MASK,true);
                break;
            }
        }
        else if (FilterShape==BINARYFILE)
        {
        	Image<double> filter;
        	filter.read(fnFilter);
            scaleToSize(BSPLINE3, maskFourierd, filter(), XSIZE(v), YSIZE(v), ZSIZE(v));
            maskFourierd.resize(Fourier);
            return;
        }

        else
        {
            maskFourierd.initZeros(Fourier);
            maskFourierd.setXmippOrigin();

            w.resizeNoCopy(3);
            for (size_t k=0; k<ZSIZE(Fourier); k++)
            {
                FFT_IDX2DIGFREQ(k,ZSIZE(v),ZZ(w));
                for (size_t i=0; i<YSIZE(Fourier); i++)
                {
                    FFT_IDX2DIGFREQ(i,YSIZE(v),YY(w));
                    for (size_t j=0; j<XSIZE(Fourier); j++)
                    {
                        FFT_IDX2DIGFREQ(j,XSIZE(v),XX(w));
                        DIRECT_A3D_ELEM(maskFourierd,k,i,j)=maskValue(w);
                    }
                }
            }
        }
    }
    else if (MULTIDIM_SIZE(v)<=1024*1024)
    {
        transformer.setReal(v);
        MultidimArray< std::complex<double> > Fourier;
        transformer.getFourierAlias(Fourier);
        if (FilterShape==BINARYFILE)
        {
        	Image<double> filter;
        	filter.read(fnFilter);
            scaleToSize(BSPLINE3, maskFourierd, filter(), XSIZE(v), YSIZE(v), ZSIZE(v));
            maskFourierd.resize(Fourier);
            return;
        }
        maskFourierd.initZeros(Fourier);
        w.resizeNoCopy(3);
        for (size_t k=0; k<ZSIZE(Fourier); k++)
        {
            FFT_IDX2DIGFREQ(k,ZSIZE(v),ZZ(w));
            for (size_t i=0; i<YSIZE(Fourier); i++)
            {
                FFT_IDX2DIGFREQ(i,YSIZE(v),YY(w));
                for (size_t j=0; j<XSIZE(Fourier); j++)
                {
                    FFT_IDX2DIGFREQ(j,XSIZE(v),XX(w));
                    DIRECT_A3D_ELEM(maskFourierd,k,i,j)=maskValue(w);
                }
            }
        }
    }
}

void FourierFilter::applyMaskSpace(MultidimArray<double> &v)
{
    MultidimArray< std::complex<double> > aux3D;
    transformer.FourierTransform(v, aux3D, false);
    applyMaskFourierSpace(v, aux3D);
    transformer.inverseFourierTransform();
}

void FourierFilter::applyMaskFourierSpace(const MultidimArray<double> &v, MultidimArray<std::complex<double> > &V)
{
    if (XSIZE(maskFourier)!=0)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V)
        DIRECT_MULTIDIM_ELEM(V,n)*=DIRECT_MULTIDIM_ELEM(maskFourier,n);
    }
    else if (XSIZE(maskFourierd)!=0)
    {
        double *ptrV=(double*)&DIRECT_MULTIDIM_ELEM(V,0);
        double *ptrMask=(double*)&DIRECT_MULTIDIM_ELEM(maskFourierd,0);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V)
        {
            *ptrV++ *= *ptrMask;
            *ptrV++ *= *ptrMask++;
        }
    }
    else if (FilterShape==SPARSIFY)
    {
        FFT_magnitude(V,vMag);
        vMag.resize(1,1,1,MULTIDIM_SIZE(vMag));
        vMag.sort(vMagSorted);
        double minMagnitude=A1D_ELEM(vMagSorted,(int)(percentage*XSIZE(vMag)));
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(V)
        if (DIRECT_MULTIDIM_ELEM(vMag,n)<minMagnitude)
        {
            double *ptr=(double*)&DIRECT_MULTIDIM_ELEM(V,n);
            *ptr=0;
            *(ptr+1)=0;
        }
    }
    else
    {
        w.resizeNoCopy(3);
        for (size_t k=0; k<ZSIZE(V); k++)
        {
            FFT_IDX2DIGFREQ(k,ZSIZE(v),ZZ(w));
            for (size_t i=0; i<YSIZE(V); i++)
            {
                FFT_IDX2DIGFREQ(i,YSIZE(v),YY(w));
                for (size_t j=0; j<XSIZE(V); j++)
                {
                    FFT_IDX2DIGFREQ(j,XSIZE(v),XX(w));
                    DIRECT_A3D_ELEM(V,k,i,j)*=maskValue(w);
                }
            }
        }
    }
}

/* Mask power -------------------------------------------------------------- */
double FourierFilter::maskPower()
{
    if (XSIZE(maskFourier) != 0)
        return maskFourier.sum2()/MULTIDIM_SIZE(maskFourier);
    else if (XSIZE(maskFourierd) != 0)
        return maskFourierd.sum2()/MULTIDIM_SIZE(maskFourierd);
    else
        return 0;
}

// Correct phase -----------------------------------------------------------
void FourierFilter::correctPhase()
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(maskFourierd)
    if (DIRECT_MULTIDIM_ELEM(maskFourierd,n)< 0)
        DIRECT_MULTIDIM_ELEM(maskFourierd,n)*= -1;
}

// Bandpass -----------------------------------------------------------------
void bandpassFilter(MultidimArray<double> &img, double w1, double w2, double raised_w)
{
    FourierFilter Filter;
    if (w1==0)
    {
        Filter.FilterBand=LOWPASS;
        Filter.w1=w2;
    }
    else if (w2==0.5)
    {
        Filter.FilterBand=HIGHPASS;
        Filter.w1=w1;
    }
    else
    {
        Filter.FilterBand=BANDPASS;
        Filter.w1=w1;
        Filter.w2=w2;
    }
    Filter.FilterShape = RAISED_COSINE;
    Filter.raised_w=raised_w;
    img.setXmippOrigin();
    Filter.generateMask(img);
    Filter.applyMaskSpace(img);
}

void gaussianFilter(MultidimArray<double> &img, double w1)
{
    FourierFilter Filter;
    Filter.FilterShape = GAUSSIAN;
    Filter.FilterBand = LOWPASS;
    Filter.w1=w1;
    img.setXmippOrigin();
    Filter.generateMask(img);
    Filter.applyMaskSpace(img);
}

void realGaussianFilter(MultidimArray<double> &img, double sigma)
{
    FourierFilter Filter;
    Filter.FilterShape = REALGAUSSIAN;
    Filter.FilterBand = LOWPASS;
    Filter.w1=sigma;
    img.setXmippOrigin();
    Filter.generateMask(img);
    Filter.applyMaskSpace(img);
}
