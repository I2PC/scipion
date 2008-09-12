/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.uam.es)
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

/* INCLUDES ---------------------------------------------------------------- */
#include <reconstruction/mlf_newtomo.h>

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char **argv)
{

    int c, nn, imgno, opt_refno;
    double aux, LL=0., sumw_allrefs=0., avePmax=0.;
    SelFile SFa;
    DocFile DFo;
    FileName fn_tmp;
    double * dataRefs;
    double * oldDataRefs;
    double * dataSigma;
    double * dataWsumRefs;
    double * dataWsumWedsPerRef;
    double * dataWsumWedsPerGroup;
    double * dataWsumDist;
    double * dataSumWRefs;
    double * dataSumWGroups;
    double * imgsPmax;
    int    * imgsOptRefNos;

    Prog_mlf_tomo_prm prm;

    // Get input parameters
    try
    {
        prm.read(argc, argv);
        prm.produceSideInfo();
        prm.show();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }

    // Reserve memory for all the required data structures
    try
    {
        dataSumWRefs         = new double[prm.nr_ref];
        dataSumWGroups       = new double[prm.nr_group];
        dataRefs             = new double[prm.nr_ref * prm.size];
        oldDataRefs          = new double[prm.nr_ref * prm.size];
        dataWsumRefs         = new double[prm.nr_ref * prm.size];
        dataWsumWedsPerRef   = new double[prm.nr_ref * prm.hsize];
        dataSigma            = new double[prm.nr_group * prm.hsize];
        dataWsumWedsPerGroup = new double[prm.nr_group * prm.hsize];
        dataWsumDist         = new double[prm.nr_group * prm.hsize];
        imgsPmax             = new double[prm.SFi.ImgNo()];
        imgsOptRefNos        = new int   [prm.SFi.ImgNo()];
    }
    catch (std::bad_alloc&)
    {
        REPORT_ERROR(1,"Error allocating memory in main");
    }

    try
    {
        // For later parallellization: split here prm.SFi
        prm.produceSideInfo2();
        prm.calculateAllFFTWs();
        prm.generateInitialReferences(dataRefs, dataWsumWedsPerRef);
        prm.estimateInitialNoiseSpectra(dataSigma);
        prm.regularise(dataRefs, dataSigma);
        prm.writeOutputFiles(0,-1,dataRefs,dataWsumWedsPerRef,DFo,sumw_allrefs,LL,avePmax);

        // Loop over all regularization steps
        for (int step = 1; step <= prm.reg_steps; step++)
        {

            if (prm.verb > 0) init_progress_bar(prm.Niter);

            // Initiaze oldDataRefs for convergence check to all-zeros 
            for (int ii=0; ii < prm.nr_ref * prm.size; ii++)
                oldDataRefs[ii] = 0.;

            // Albertop-like decrease
            prm.reg = exp(log(prm.reg0) - (log(prm.reg0) - log(XMIPP_MAX(prm.regF,0.001)))*(step-1)/(prm.reg_steps-1));
            if (prm.verb > 0) 
                std::cerr<<std::endl;
                std::cerr << "  Sub-tomogram classification:  step "<<step<<" of "<<prm.reg_steps<<" with regularisation= "<<prm.reg<<std::endl;

            // Iterate until convergence (or Niter) is reached
            for (int iter = prm.istart; iter <= prm.Niter; iter++)
            {

                //if (prm.verb > 0) std::cerr << "  Sub-tomogram classification:  step "<<step<<" iteration " << iter << " of " << prm.Niter << std::endl;

                prm.expectation(dataRefs,
                                dataSigma,
                                dataWsumRefs, 
                                dataWsumWedsPerRef,
                                dataWsumWedsPerGroup,
                                dataWsumDist,
                                dataSumWRefs,
                                LL, 
                                avePmax,
                                DFo);

                prm.maximization(dataRefs,
                                 dataSigma,
                                 dataWsumRefs,
                                 dataWsumWedsPerRef,
                                 dataWsumWedsPerGroup,
                                 dataWsumDist,
                                 dataSumWRefs,
                                 sumw_allrefs, 
                                 avePmax);


                if (prm.checkConvergence(dataRefs,
                                         oldDataRefs) )
                    break;

                if (prm.verb > 0) progress_bar(iter);

            } // end loop iterations
            if (prm.verb > 0) progress_bar(prm.Niter);

            prm.writeOutputFiles(step, 
                                 -1,
                                 dataRefs,
                                 dataWsumWedsPerRef,
                                 DFo,
                                 sumw_allrefs, 
                                 LL, 
                                 avePmax);

            // linear decrease of regularisation parameter
            //prm.reg -= prm.delta_reg;
            //prm.reg = XMIPP_MAX(prm.reg,0.);

        } // end loop regularization steps

        // Write out converged structures
        prm.writeOutputFiles(-1,
                             -1, 
                             dataRefs,
                             dataWsumWedsPerRef,
                             DFo,
                             sumw_allrefs, 
                             LL, 
                             avePmax);

        // Write out docfile with optimal transformation & references
        fn_tmp=prm.fn_root+".doc";
        DFo.write(fn_tmp);

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}



