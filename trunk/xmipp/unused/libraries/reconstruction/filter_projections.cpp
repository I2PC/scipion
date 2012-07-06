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

#include "filter_projections.h"
#include <data/projection.h>
#include <data/args.h>
#include <data/filters.h>
#include <data/histogram.h>
#include <pthread.h>

/* Read parameters --------------------------------------------------------- */
void Prog_Filter_Projections_Parameters::read(int argc, char **argv)
{
    if (argc==1)
        REPORT_ERROR(ERR_ARG_MISSING,"No filtering criteria provided");

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
            REPORT_ERROR(ERR_ARG_MISSING,"Not enough arguments after -filter_score");
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
            REPORT_ERROR(ERR_ARG_MISSING,"Not enough arguments after -filter_cost");
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
            REPORT_ERROR(ERR_ARG_MISSING,"Not enough arguments after -filter_movement");
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
            REPORT_ERROR(ERR_ARG_MISSING,"Not enough arguments after -filter_normalization");
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
    << "  -i <docfile>                          : Name of the input docfile\n"
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
    int Nimg =  DF_in.size();

    if (percentil_score>0)
    {
        DF_score.read(fn_score);
        if (Nimg!=DF_score.size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                         "The number of images in score docfile is not the same as in the input docfile");
    }

    if (percentil_cost>0)
    {
        DF_cost.read(fn_cost);
        if (Nimg!=DF_cost.size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                         "The number of images in cost docfile is not the same as in the input docfile");
    }

    if (angleLimit>0)
    {
        DF_movement0.read(fn_movement0);
        angleLimit=DEG2RAD(angleLimit);
        if (Nimg!=DF_movement0.size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                         "The number of images in movement docfile is not the same as in the input docfile");
    }

    if (percentil_normalization>0)
    {
        V.read(fn_vol);
        V().setXmippOrigin();

        Image<double> I;
        FileName fnAux;
        DF_in.getValue(MDL_IMAGE, fnAux, DF_in.firstObject());
        I.read(fnAux);
        I().setXmippOrigin();

        Mask aux;
        aux.type = BINARY_CIRCULAR_MASK;
        aux.mode = INNER_MASK;
        aux.R1 = r1;
        aux.resize(I());
        aux.generate_mask();
        typeCast(aux.get_binary_mask(),dhardMask);
        ihardMask=aux.get_binary_mask();
        dhardMask.setXmippOrigin();
        ihardMask.setXmippOrigin();

        Mask aux2;
        aux2.type = RAISED_COSINE_MASK;
        aux2.mode = INNER_MASK;
        aux2.R1 = r1;
        aux2.R2 = r2;
        aux2.resize(I());
        aux2.generate_mask();
        softMask=aux2.get_cont_mask();
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
struct FilterByNormalizationParams
{
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
    MetaData DF_in_local=prm->DF_in;

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
    Image<double> I;
    FileName fnAux;
    for (int i=0; i<Nimg; i++)
    {
        local_correlations[i]=-2;
        if (i%prm->Nthreads==idThread && prm->valid[i])
        {
            // Get the experimental image
            DF_in_local.getValue(MDL_IMAGE,fnAux,i+1);
            I.readApplyGeo(DF_in_local, i+1);
            //selfTranslate(BSPLINE3,I(),vectorR2(I.Xoff(),I.Yoff()));
            //I.setShifts(0.0,0.0);

            // Get the corresponding theoretical projection
            Projection P;
            projectVolume(prm->V(), P, YSIZE(I()), XSIZE(I()),
                           I.rot(), I.tilt(), I.psi());

            // Compute correlation within mask
            local_correlations[i]=correlationIndex(P(),I(),
                                                    &(prm->ihardMask));

            // Make the noise within the mask to be zero mean and with
            // standard deviation 1
            P()-=I();
            double min_val, max_val, avg, stddev;
            computeStats_within_binary_mask(prm->ihardMask, P(),
                                            min_val, max_val, avg, stddev);
            I()+=avg;
            I()*=1/stddev;

            // Mask the image with a raised cosine
            I()*=prm->softMask;

            // Write the adjusted image
            I.write();
            if (idThread==0 && i%60==0)
                progress_bar(i);
        }
    }
    if (idThread==0)
        progress_bar(Nimg);

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
        std::vector<double> score;
        DF_score.getColumnValues(MDL_MAXCC,score);
        Histogram1D Hscore;
        compute_hist(score,Hscore,200);
        double threshold=Hscore.percentil(percentil_score);
        std::cout << "Score threshold=" << threshold << std::endl;
        int imax=score.size();
        for (int i=0; i<imax; i++)
            if (score[i]<threshold)
                valid[i]=false;
    }

    // Filter by cost .....................................................
    if (percentil_cost>0)
    {
        // Compute the histogram of the costs
        std::vector<double> cost;
        DF_cost.getColumnValues(MDL_COST,cost);
        Histogram1D Hcost;
        compute_hist(cost,Hcost,200);
        double threshold=Hcost.percentil(100-percentil_cost);
        std::cout << "Cost threshold=" << threshold << std::endl;
        int imax=cost.size();
        for (int i=0; i<imax; i++)
            if (cost[i]>threshold)
                valid[i]=false;
    }

    // Filter by movement ..................................................
    if (angleLimit>0)
    {
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA2(DF_movement0, DF_in)
        {
            double rot0, tilt0, psi0, shiftX0, shiftY0;
            DF_movement0.getValue(MDL_ANGLE_ROT,rot0,__iter.objId);
            DF_movement0.getValue(MDL_ANGLE_TILT,tilt0,__iter.objId);
            DF_movement0.getValue(MDL_ANGLE_PSI,psi0,__iter.objId);
            DF_movement0.getValue(MDL_SHITF_X,shiftX0,__iter.objId);
            DF_movement0.getValue(MDL_SHITF_Y,shiftY0,__iter.objId);

            double rotF, tiltF, psiF, shiftXF, shiftYF;
            DF_in.getValue(MDL_ANGLE_ROT,rotF,__iter2.objId);
            DF_in.getValue(MDL_ANGLE_TILT,tiltF,__iter2.objId);
            DF_in.getValue(MDL_ANGLE_PSI,psiF,__iter2.objId);
            DF_in.getValue(MDL_SHITF_X,shiftXF,__iter2.objId);
            DF_in.getValue(MDL_SHITF_Y,shiftYF,__iter2.objId);

            double diffX=shiftXF-shiftX0;
            double diffY=shiftYF-shiftY0;

            Matrix2D<double> E0, EF;
            Euler_angles2matrix(rot0,tilt0,psi0,E0);
            Euler_angles2matrix(rotF,tiltF,psiF,EF);
            double angularDistance=acos(Euler_distanceBetweenMatrices(E0,EF));
            if (angularDistance>angleLimit ||
                sqrt(diffX*diffX+diffY*diffY)>shiftLimit)
                valid[i]=false;

            i++;
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
                if (correlations[i]<0)
                    valid[i]=false;
                else if (tanh(correlations[i]+eps)<=0)
                    valid[i]=false;
                else
                    validCorrelations.push_back(correlations[i]);
            }

        // Filter by low correlations
        MultidimArray<double> vcorrelations;
        vcorrelations.resize(validCorrelations.size());
        FOR_ALL_ELEMENTS_IN_ARRAY1D(vcorrelations)
        vcorrelations(i)=validCorrelations[i];
        Histogram1D Hcorr;
        compute_hist(vcorrelations,Hcorr,200);
        double threshold=Hcorr.percentil(percentil_normalization);
        std::cout << "Normalization threshold=" << threshold << std::endl;
        for (int i=0; i<Nimg; i++)
            if (valid[i] && correlations[i]<threshold)
                valid[i]=false;
    }

    // Produce the output docfile ..........................................
    MetaData DF_out;
    //TODO: CHECK?????
    //std::vector< std::string > imagenames;
    //DF_in.getColumnValues(MDL_IMAGE,imagenames);
    FileName imgFn;
    size_t id;
    int i = 0;
    //for (int i=0; i<Nimg; i++)
    FOR_ALL_OBJECTS_IN_METADATA(DF_in)
    {
        if (valid[i])
        {

            double rotF, tiltF, psiF, shiftXF, shiftYF;
            DF_in.getValue(MDL_IMAGE,imgFn,__iter.objId);
            DF_in.getValue(MDL_ANGLE_ROT,rotF,__iter.objId);
            DF_in.getValue(MDL_ANGLE_TILT,tiltF,__iter.objId);
            DF_in.getValue(MDL_ANGLE_PSI,psiF,__iter.objId);
            if (fn_vol!="")
            {
                shiftXF=0;
                shiftYF=0;
            }
            else
            {
                DF_in.getValue(MDL_SHITF_X,shiftXF,__iter.objId);
                DF_in.getValue(MDL_SHITF_Y,shiftYF,__iter.objId);
            }
            id = DF_out.addObject();
            DF_out.setValue(MDL_IMAGE,imgFn,id);
            DF_out.setValue(MDL_ANGLE_ROT,rotF,id);
            DF_out.setValue(MDL_ANGLE_TILT,tiltF,id);
            DF_out.setValue(MDL_ANGLE_PSI,psiF,id);
            DF_out.setValue(MDL_SHITF_X,shiftXF,id);
            DF_out.setValue(MDL_SHITF_Y,shiftYF,id);
        }
        i++;
    }
    DF_out.write(fn_out+".doc");
    std::cout << DF_out.size() << " images have been kept out of "
    << Nimg << std::endl;
}
