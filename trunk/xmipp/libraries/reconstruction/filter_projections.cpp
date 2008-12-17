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

#include "filter_projections.h"
#include "projection.h"
#include <data/args.h>
#include <data/filters.h>
#include <data/histogram.h>
#include <pthread.h>

/* Read parameters --------------------------------------------------------- */
void Prog_Filter_Projections_Parameters::read(int argc, char **argv)
{
    if (argc==1) REPORT_ERROR(1,"No filtering criteria provided");

    fn_in=getParameter(argc,argv,"-i");
    fn_out=getParameter(argc,argv,"-o");

    int i=paremeterPosition(argc,argv,"-filter_score");
    if (i!=-1)
    {
        if (i+2<argc)
        {
            fn_score=argv[i+1];
            percentil_score=textToFloat(argv[i+2]);
        }
        else
            REPORT_ERROR(1,"Not enough arguments after -filter_score");
    }
    else
    {
        fn_score="";
        percentil_score=-1;
    }

    i=paremeterPosition(argc,argv,"-filter_cost");
    if (i!=-1)
    {
        if (i+2<argc)
        {
            fn_cost=argv[i+1];
            percentil_cost=textToFloat(argv[i+2]);
        }
        else
            REPORT_ERROR(1,"Not enough arguments after -filter_cost");
    }
    else
    {
        fn_cost="";
        percentil_cost=-1;
    }

    i=paremeterPosition(argc,argv,"-filter_movement");
    if (i!=-1)
    {
        if (i+3<argc)
        {
            fn_movement0=argv[i+1];
            angleLimit=textToFloat(argv[i+2]);
            shiftLimit=textToFloat(argv[i+3]);
        }
        else
            REPORT_ERROR(1,"Not enough arguments after -filter_movement");
    }
    else
    {
        fn_movement0="";
        angleLimit=shiftLimit=-1;
    }

    i=paremeterPosition(argc,argv,"-filter_normalization");
    if (i!=-1)
    {
        if (i+4<argc)
        {
            fn_vol=argv[i+1];
            r1=textToFloat(argv[i+2]);
            r2=textToFloat(argv[i+3]);
            percentil_normalization=textToFloat(argv[i+4]);
        }
        else
            REPORT_ERROR(1,"Not enough arguments after -filter_normalization");
    }
    else
    {
        fn_vol="";
        r1=r2=percentil_normalization=-1;
    }
    
    Nthreads=textToInteger(getParameter(argc,argv,"-thr","1"));
}

/* Usage ------------------------------------------------------------------- */
void Prog_Filter_Projections_Parameters::usage()
{
    std::cerr
        << "Usage:\n"
        << "  -o <rootname>                         : Name of the output docfile and selfile\n"
        << "  [-filter_score <docfile> <percentil>] : The lower percentil is removed\n"
        << "                                        : This file is produced by xmipp_angular_discrete_assign\n"
        << "  [-filter_cost <docfile> <percentil>]  : The higher percentil is removed\n"
        << "                                        : This file is produced by xmipp_angular_continuous_assign\n"
        << "  [-filter_movement <docfile0> <angleLimit> <shiftLimit>]:\n"
        << "                                        : Particles moving more than a certain limit are not considered\n"
        << "  [-filter_normalization <volume> <r0> <rF> <percentil> : The worse fitting particles will be removed\n"
        << "  [-thr <N=1>]                          : Number of threads available\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Filter_Projections_Parameters::show()
{
    std::cout
        << "Filter score docfile: " << fn_score
        << " percentil=" << percentil_score << std::endl
        << "Filter cost docfile: " << fn_cost
        << " percentil=" << percentil_cost << std::endl
        << "Filter movement docfile0: " << fn_movement0
        << " angleLimit=" << angleLimit
        << " shiftLimit=" << shiftLimit << std::endl
        << "Filter normalization volume: " << fn_vol
        << " r1=" << r1
        << " r2=" << r2
        << " percentil=" << percentil_normalization << std::endl
        << "Number of threads=" << Nthreads << std::endl
    ;
}

/* Produce side information ------------------------------------------------ */
void Prog_Filter_Projections_Parameters::produce_side_info()
{
    DF_in.read(fn_in);
    int Nimg=DF_in.dataLineNo();

    if (percentil_score>0)
    {
        DF_score.read(fn_score);
        if (Nimg!=DF_score.dataLineNo())
            REPORT_ERROR(1,"The number of images in score docfile is not the same as in the input docfile");
    }

    if (percentil_cost>0)
    {
        DF_cost.read(fn_cost);
        if (Nimg!=DF_cost.dataLineNo())
            REPORT_ERROR(1,"The number of images in cost docfile is not the same as in the input docfile");
    }

    if (angleLimit>0)
    {
        DF_movement0.read(fn_movement0);
        angleLimit=DEG2RAD(angleLimit);
        if (Nimg!=DF_movement0.dataLineNo())
            REPORT_ERROR(1,"The number of images in movement docfile is not the same as in the input docfile");
    }

    if (percentil_normalization>0)
    {
        V.read(fn_vol);
        V().setXmippOrigin();
        
        ImageXmipp I(DF_in.get_imagename(1));
        I().setXmippOrigin();
        
        Mask_Params aux;
        aux.type = BINARY_CIRCULAR_MASK;
        aux.mode = INNER_MASK;
        aux.R1 = r1;
        aux.resize(I());
        aux.generate_2Dmask();
        typeCast(aux.imask2D,dhardMask);
        ihardMask=aux.imask2D;
        dhardMask.setXmippOrigin();
        ihardMask.setXmippOrigin();
        
        Mask_Params aux2;
        aux2.type = RAISED_COSINE_MASK;
        aux2.mode = INNER_MASK;
        aux2.R1 = r1;
        aux2.R2 = r2;
        aux2.resize(I());
        aux2.generate_2Dmask();
        softMask=aux2.dmask2D;
        softMask.setXmippOrigin();
    }
    
    valid.resize(Nimg);
    correlations.resize(Nimg);
    for (int i=0; i<Nimg; i++)
    {
        valid[i]=true;
        correlations[i]=-2;
    }
}

/* Filter by normalization ------------------------------------------------- */
struct FilterByNormalizationParams {
    Prog_Filter_Projections_Parameters *prm;
    int idThread;
};

static pthread_mutex_t correlationMutex = PTHREAD_MUTEX_INITIALIZER;

void * filterByNormalizationThread(void *args)
{
    // Pick the input parameters
    FilterByNormalizationParams * argsprm=
        (FilterByNormalizationParams *)args;
    Prog_Filter_Projections_Parameters *prm=argsprm->prm;
    int idThread=argsprm->idThread;
    
    // Make a local copy of DF_in (this is because of the threads)
    DocFile DF_in_local=prm->DF_in;
    
    // Create a local copy of the correlations observed
    std::vector<double> local_correlations;
    local_correlations.resize(prm->valid.size());
    
    // Evaluate all projections corresponding to this thread
    int Nimg=prm->valid.size();
    if (idThread==0)
    {
        std::cout << "Filtering by normalization ...\n";
        init_progress_bar(Nimg);
    }
    for (int i=0; i<Nimg; i++)
    {
        local_correlations[i]=-2;
        if (i%prm->Nthreads==idThread && prm->valid[i])
        {
            // Get the experimental image
            ImageXmipp I;
            DF_in_local.get_image(i+1,I);
            I().selfTranslateBSpline(3,vectorR2(I.Xoff(),I.Yoff()));
        
            // Get the corresponding theoretical projection
            Projection P;
            project_Volume(prm->V(), P, YSIZE(I()), XSIZE(I()),
                I.rot(), I.tilt(), I.psi());

            // Adjust the two projections and compute correlation within mask
            rangeAdjust_within_mask(&(prm->dhardMask),P(),I());
            local_correlations[i]=correlation_index(P(),I(),&(prm->ihardMask));
            
            // Mask the image with a raised cosine
            I()*=prm->softMask;

            // Write the adjusted image
            I.write();
            if (idThread==0 && i%60==0) progress_bar(i);
        }
    }
    if (idThread==0) progress_bar(Nimg);
    
    // Dump the results onto the parent correlation vector
    pthread_mutex_lock( &correlationMutex );
    for (int i=0; i<Nimg; i++)
        if (local_correlations[i]>-2)
            prm->correlations[i]=local_correlations[i];
    pthread_mutex_unlock( &correlationMutex );        
}

/* Run --------------------------------------------------------------------- */
void Prog_Filter_Projections_Parameters::run()
{
    int Nimg=valid.size();
    
    // Filter by score .....................................................
    if (percentil_score>0)
    {
        // Compute the histogram of the scores
        int col_score=DF_score.getColNumberFromHeader("Score")  - 1;
        if (col_score<0)
            REPORT_ERROR(1,"Column Score not found in the score docfile");
        Matrix1D<double> score=DF_score.col(col_score);
        histogram1D Hscore;
        compute_hist(score,Hscore,200);
        double threshold=Hscore.percentil(percentil_score);
        std::cout << "Score threshold=" << threshold << std::endl;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(score)
            if (score(i)<threshold)
                valid[i]=false;
    }
    
    // Filter by cost .....................................................
    if (percentil_cost>0)
    {
        // Compute the histogram of the costs
        int col_cost=DF_cost.getColNumberFromHeader("Cost")  - 1;
        if (col_cost<0)
            REPORT_ERROR(1,"Column Cost not found in the cost docfile");
        Matrix1D<double> cost=DF_cost.col(col_cost);
        histogram1D Hcost;
        compute_hist(cost,Hcost,200);
        double threshold=Hcost.percentil(100-percentil_cost);
        std::cout << "Cost threshold=" << threshold << std::endl;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(cost)
            if (cost(i)>threshold)
                valid[i]=false;
    }
    
    // Filter by movement ..................................................
    if (angleLimit>0)
    {
        DF_movement0.go_first_data_line();
        DF_in.go_first_data_line();
        int i=0;
        while (!DF_movement0.eof())
        {
            double rot0    = DF_movement0(0);
            double tilt0   = DF_movement0(1);
            double psi0    = DF_movement0(2);
            double shiftX0 = DF_movement0(3);
            double shiftY0 = DF_movement0(4);

            double rotF    = DF_in(0);
            double tiltF   = DF_in(1);
            double psiF    = DF_in(2);
            double shiftXF = DF_in(3);
            double shiftYF = DF_in(4);
            
            double diffX=shiftXF-shiftX0;
            double diffY=shiftYF-shiftY0;
            
            Matrix2D<double> E0, EF;
            Euler_angles2matrix(rot0,tilt0,psi0,E0);
            Euler_angles2matrix(rotF,tiltF,psiF,EF);
            if (Euler_distanceBetweenMatrices(E0,EF)>angleLimit ||
                sqrt(diffX*diffX+diffY*diffY)>shiftLimit)
                valid[i]=false;

            i++;
            DF_movement0.next_data_line();
            DF_in.next_data_line();
        }
    }

    // Filter by normalization .............................................
    if (fn_vol!="")
    {
        pthread_t * th_ids = new pthread_t[Nthreads];
        FilterByNormalizationParams * th_args=
            new FilterByNormalizationParams[Nthreads];
        for( int nt = 0 ; nt < Nthreads ; nt ++ )
        {
            th_args[nt].prm = this;
            th_args[nt].idThread = nt;
            pthread_create( &th_ids[nt], NULL,
                filterByNormalizationThread, &th_args[nt]);
        }

        // Waiting for threads to finish
        for( int nt=0; nt<Nthreads; nt++)
            pthread_join(th_ids[nt], NULL);

        // Thread structures are not needed any more
        delete [] th_ids;
        delete [] th_args;

        // Check that all correlations are significantly different from 0
        std::vector<double> validCorrelations;
        double eps=-1.96/sqrt(dhardMask.sum()-3);
        for (int i=0; i<Nimg; i++)
            if (valid[i])
            {
                if (correlations[i]<0) valid[i]=false;
                else if (tanh(correlations[i]+eps)<=0) valid[i]=false;
                else validCorrelations.push_back(correlations[i]);
            }
        
        // Filter by low correlations
        Matrix1D<double> vcorrelations;
        vcorrelations.resize(validCorrelations.size());
        FOR_ALL_ELEMENTS_IN_MATRIX1D(vcorrelations)
            vcorrelations(i)=validCorrelations[i];
        histogram1D Hcorr;
        compute_hist(vcorrelations,Hcorr,200);
        double threshold=Hcorr.percentil(percentil_normalization);
        std::cout << "Normalization threshold=" << threshold << std::endl;
        for (int i=0; i<Nimg; i++)
            if (valid[i] && correlations[i]<threshold)
                valid[i]=false;
    }

    // Produce the output docfile ..........................................
    DocFile DF_out;
    SelFile SF_out;
    DF_out.append_comment("Headerinfo columns: rot (1) , tilt (2), psi (3), Xoff (4), Yoff (5)");
    for (int i=0; i<Nimg; i++)
        if (valid[i])
        {
            FileName imagename=DF_in.get_imagename(i+1);
            DF_out.append_comment(imagename);
            SF_out.insert(imagename);
            Matrix1D<double> v(5);
            v(0)=DF_in(i+1,0);
            v(1)=DF_in(i+1,1);
            v(2)=DF_in(i+1,2);
            if (fn_vol!="")
            {
                v(3)=0;
                v(4)=0;
            }
            else
            {
                v(3)=DF_in(i+1,3);
                v(4)=DF_in(i+1,4);
            }
            DF_out.append_data_line(v);
        }
    DF_out.write(fn_out+".doc");
    SF_out.write(fn_out+".sel");
    std::cout << DF_out.dataLineNo() << " images have been kept out of "
              << Nimg << std::endl;
}
