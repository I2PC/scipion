/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "common_lines.h"

#include <data/args.h>
#include <data/mask.h>
#include <data/filters.h>
#include <reconstruction/fourier_filter.h>
#include <reconstruction/radon.h>
#include <fstream>

/* Common line ------------------------------------------------------------- */
CommonLine::CommonLine()
{
    angi=angj=0;
    distanceij=-1;
}

/* Read parameters --------------------------------------------------------- */
void CommonLine_Parameters::read(int argc, char **argv)
{
    fn_sel = getParameter(argc, argv, "-i");
    fn_out = getParameter(argc, argv, "-o", "commonlines.txt");
    lpf    = textToFloat(getParameter(argc, argv, "-lpf", "0.01"));
    hpf    = textToFloat(getParameter(argc, argv, "-hpf", "0.35"));
    stepAng= textToFloat(getParameter(argc, argv, "-stepAng", "3"));
    mem    = textToFloat(getParameter(argc, argv, "-mem", "1"));
    Nthr   = textToInteger(getParameter(argc, argv, "-thr", "1"));
    Nmpi   = 1;
    distance=CORRENTROPY;
    if (checkParameter(argc,argv,"-correlation"))
        distance=CORRELATION;
}

/* Usage ------------------------------------------------------------------- */
void CommonLine_Parameters::usage()
{
    std::cout
        << "Usage: common_lines [Purpose and Parameters]\n"
        << "Purpose: For every pair of images in the selfile, find the common line\n"
        << "Parameter Values: (notice space before value)\n"
        << "    -i <file_in>        : input selfile\n"
        << "   [-o <file_out>]      : if no name is given, commonlines.txt\n"
        << "   [-correlation]       : use correlation instead of correntropy\n"
        << "   [-lpf <f=0.01>]      : low pass frequency (<0.5)\n"
        << "   [-hpf <f=0.35>]      : high pass frequency (<0.5)\n"
        << "   [-stepAng <s=3>]     : angular step\n"
        << "   [-mem <m=1>]         : float number with the memory available in Gb\n"
        << "   [-thr <thr=1>]       : number of threads for each process\n"
    ;
}

/* Side info --------------------------------------------------------------- */
void CommonLine_Parameters::produceSideInfo()
{
    SF.read(fn_sel);
    Nimg=SF.ImgNo();

    // Compute the number of images in each block
    int Xdim, Ydim;
    SF.ImgSize(Ydim,Xdim);
    Nblock=FLOOR(sqrt(mem*pow(2.0,30.0)/(Ydim*Xdim*sizeof(double))));
    Nblock=XMIPP_MIN(Nblock,CEIL(((float)Nimg)/Nmpi));
    
    // Ask for memory for the common line matrix
    CommonLine dummy;
    for (int i=0; i<Nimg*Nimg; i++)
        CLmatrix.push_back(dummy);
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &out, const CommonLine_Parameters &prm)
{
    out << "File in:       " << prm.fn_sel  << std::endl
        << "File out:      " << prm.fn_out  << std::endl
        << "Lowpass:       " << prm.lpf     << std::endl
        << "Highpass:      " << prm.hpf     << std::endl
        << "StepAng:       " << prm.stepAng << std::endl
        << "Memory(Gb):    " << prm.mem     << std::endl
        << "Block size:    " << prm.Nblock  << std::endl
        << "N.Threads:     " << prm.Nthr    << std::endl;
    if (prm.distance==CORRENTROPY)
        out << "Distance:      Correntropy\n";
    else if (prm.distance==CORRELATION)
        out << "Distance:      Correlation\n";
    return out;
}

/* Get and prepare block --------------------------------------------------- */
struct ThreadPrepareImages
{
    int myThreadID;
    CommonLine_Parameters * parent;
    SelFile *SFi;
    std::vector< Matrix2D<double> > *blockRTs;
    
    double sigma;
};

void * threadPrepareImages( void * args )
{
    ThreadPrepareImages * master = (ThreadPrepareImages *) args;
    CommonLine_Parameters * parent = master->parent;
    SelFile SFi=*(master->SFi);
    SFi.go_beginning();
    int Ydim, Xdim;
    SFi.ImgSize(Ydim,Xdim);

    Matrix2D<int> mask;
    mask.resize(Ydim,Xdim);
    mask.setXmippOrigin();
    BinaryCircularMask(mask,Xdim/2, OUTSIDE_MASK);

    FourierMask Filter;
    Filter.FilterBand=BANDPASS;
    Filter.w1=parent->lpf;
    Filter.w2=parent->hpf;
    Filter.raised_w=Filter.w1/3;

    int i=0;
    bool first=true;
    master->sigma=0;
    int Nsigma=0;
    while (!SFi.eof())
    {
        if ((i+1)%parent->Nthr==master->myThreadID)
        {
            ImageXmipp I;
            I.read(SFi.NextImg());
            I().setXmippOrigin();

            if (first)
            {
                Filter.generate_mask(I());
                first=false;
            }
            Filter.apply_mask_Space(I());
            Matrix2D<double> RT;
            Radon_Transform(I(),parent->stepAng,RT);
            (*(master->blockRTs))[i]=RT;

            double min_val, max_val, avg, stddev;
            computeStats_within_binary_mask(mask,I(),
                min_val, max_val, avg, stddev);
            master->sigma+=stddev;
            Nsigma++;
        }
        else
            SFi.next();
        i++;
    }
    master->sigma/=Nsigma;
}

void CommonLine_Parameters::getAndPrepareBlock(int i,
    std::vector< Matrix2D<double> > &blockRTs)
{
    // Get the selfile
    SelFile SFi;
    SF.chooseSubset(i*Nblock,XMIPP_MIN((i+1)*Nblock,Nimg)-1,SFi);
    
    // Ask for space for all the block images
    int jmax=SFi.ImgNo();
    Matrix2D<double> dummy;
    for (int j=0; j<jmax; j++)
        blockRTs.push_back(dummy);
    
    // Read and preprocess the images
    pthread_t * th_ids = new pthread_t[Nthr];
    ThreadPrepareImages * th_args = new ThreadPrepareImages[Nthr];
    for( int nt = 0 ; nt < Nthr ; nt ++ )
    {
        // Passing parameters to each thread
        th_args[nt].parent = this;
        th_args[nt].myThreadID = nt;
        th_args[nt].SFi = &SFi;
        th_args[nt].blockRTs = &blockRTs;
        pthread_create( (th_ids+nt) , NULL, threadPrepareImages,
            (void *)(th_args+nt) );
    }
		
    // Waiting for threads to finish
    sigma=0;
    for( int nt = 0 ; nt < Nthr ; nt ++ )
    {
        pthread_join(*(th_ids+nt), NULL);
        sigma+=th_args[nt].sigma;
    }
    sigma/=Nthr;
    int Ydim, Xdim;
    SF.ImgSize(Ydim,Xdim);
    sigma*=sqrt(2*Xdim);
    
    // Threads structures are not needed any more
    delete( th_ids );
    delete( th_args );
}

/* Process block ----------------------------------------------------------- */
void commonLineTwoImages(std::vector< Matrix2D<double> > &RTsi, int idxi,
    std::vector< Matrix2D<double> >&RTsj, int idxj,
    CommonLine_Parameters *parent, CommonLine &result)
{
    Matrix2D<double> &RTi=RTsi[idxi];
    Matrix2D<double> &RTj=RTsj[idxj];
    
    result.distanceij=2;
    result.angi=-1;
    result.angj=-1;
    for (int ii=0; ii<YSIZE(RTi); ii++)
    {
        Matrix1D<double> linei;
        linei.initZeros(XSIZE(RTi));
        linei.setXmippOrigin();
        FOR_ALL_ELEMENTS_IN_MATRIX1D(linei) linei(i) = RTi(ii,i);

        for (int jj=0; jj<YSIZE(RTj); jj++)
        {
            Matrix1D<double> linej;
            linej.initZeros(XSIZE(RTi));
            linej.setXmippOrigin();
            FOR_ALL_ELEMENTS_IN_MATRIX1D(linej) linej(i) = RTj(jj,i);
            
            // Compute distance between the two lines
            double distance;
            if (parent->distance==CORRELATION)
                distance=1-(correlation_index(linei,linej)+1)/2;
            else if (parent->distance==CORRENTROPY)
                distance=1-correntropy(linei,linej,parent->sigma);

            // Check if this is the best match
            if (distance<result.distanceij)
            {
                result.distanceij=distance;
                result.angi=ii;
                result.angj=jj;
            }
        }
    }
    result.angi*=parent->stepAng;
    result.angj*=parent->stepAng;
}

struct ThreadCompareImages
{
    int myThreadID;
    CommonLine_Parameters * parent;
    int i;
    int j;
    std::vector< Matrix2D<double> > *RTsi;
    std::vector< Matrix2D<double> > *RTsj;
};

void * threadCompareImages( void * args )
{
    ThreadCompareImages * master = (ThreadCompareImages *) args;
    CommonLine_Parameters * parent = master->parent;

    int blockIsize=master->RTsi->size();
    int blockJsize=master->RTsj->size();
    for (int i=0; i<blockIsize; i++)
    {
        long int ii=parent->Nblock*master->i+i;
        for (int j=0; j<blockJsize; j++)
        {
            // Check if this two images have to be compared
            long int jj=parent->Nblock*master->j+j;
            if (ii>=jj) continue;
            if ((ii*blockJsize+jj+1)%parent->Nthr!=master->myThreadID)
                continue;
            
            // Effectively compare the two images
            commonLineTwoImages(*(master->RTsi),i,*(master->RTsj),j,
                parent,parent->CLmatrix[ii*parent->Nimg+jj]);
        }
    }
}

void CommonLine_Parameters::processBlock(int i, int j)
{
    if (i>j) return;

    // Preprocess each one of the selfiles
    std::vector< Matrix2D<double> > RTsi, RTsj;
    getAndPrepareBlock(i,RTsi);
    if (i!=j) getAndPrepareBlock(j,RTsj);

    // Compare all versus all
    // Read and preprocess the images
    pthread_t * th_ids = new pthread_t[Nthr];
    ThreadCompareImages * th_args = new ThreadCompareImages[Nthr];
    for( int nt = 0 ; nt < Nthr ; nt ++ )
    {
        // Passing parameters to each thread
        th_args[nt].parent = this;
        th_args[nt].myThreadID = nt;
        th_args[nt].i = i;
        th_args[nt].j = j;
        th_args[nt].RTsi = &RTsi;
        if (i!=j)
            th_args[nt].RTsj = &RTsj;
        else
            th_args[nt].RTsj = &RTsi;
        pthread_create( (th_ids+nt) , NULL, threadCompareImages,
            (void *)(th_args+nt) );
    }
		
    // Waiting for threads to finish
    for( int nt = 0 ; nt < Nthr ; nt ++ )
        pthread_join(*(th_ids+nt), NULL);
    
    // Threads structures are not needed any more
    delete( th_ids );
    delete( th_args );
}

/* Write results ----------------------------------------------------------- */
void CommonLine_Parameters::writeResults() {
    // Look for the minimum and maximum of the common line matrix
    double minVal=2, maxVal=-2;
    for (int i=0; i<Nimg; i++)
        for (int j=0; j<Nimg; j++)
        {
            double val=CLmatrix[i*Nimg+j].distanceij;
            if (val>0)
            {
                minVal=XMIPP_MIN(minVal,val);
                maxVal=XMIPP_MAX(maxVal,val);
            }
        }

    // Write the common line matrix
    std::ofstream fh_out;
    fh_out.open(fn_out.c_str());
    if (!fh_out)
        REPORT_ERROR(1,(std::string)"Cannot open "+fn_out+" for writing");
    for (int i=0; i<Nimg; i++)
        for (int j=0; j<Nimg; j++)
        {
            int ii=i*Nimg+j;
            if (CLmatrix[ii].distanceij>0)
                fh_out << j << " " << i << " "
                       << ROUND(65535*(CLmatrix[ii].distanceij-minVal)/
                             (maxVal-minVal)) << " "
                       << CLmatrix[ii].angi << " "
                       << CLmatrix[ii].angj << std::endl;
        }
    fh_out.close();
}

/* Main program ------------------------------------------------------------ */
void CommonLine_Parameters::run(int rank)
{
    // Process all blocks
    int maxBlock=CEIL(((float)Nimg)/Nblock);
    if (rank==0)
        init_progress_bar(maxBlock*maxBlock);
    for (int i=0; i<maxBlock; i++)
        for (int j=0; j<maxBlock; j++)
        {
            int numBlock=i*maxBlock+j;
            if ((numBlock+1)%Nmpi==rank)
                processBlock(i,j);
            if (rank==0) progress_bar(i*maxBlock+j);
        }
    if (rank==0) progress_bar(maxBlock*maxBlock);
}
