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
    double aux, LL, sumw_allrefs, avePmax;
    std::vector<double> conv;
    SelFile SFa;
    double * dataRefs;
    double * dataSigma;
    double * dataWsumRefs;
    double * dataWsumWedsPerRef;
    double * dataWsumWedsPerGroup;
    double * dataWsumDist;
    double * dataSumWRefs;
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
        dataRefs             = new double[prm.nr_ref * prm.size];
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
        // Read references in memory and store their FTs
        prm.readAndFftwAllReferences(dataRefs);
        // For later parallellization: split here prm.SFi and prm.SFw!
        //SFa = prm.SFi + prm.SFw;
        prm.calculateAllFFTWs(prm.SFi);
        prm.estimateInitialNoiseSpectra(dataSigma);

        // Loop over all iterations
        for (int iter = prm.istart; iter <= prm.Niter; iter++)
        {

            if (prm.verb > 0) std::cerr << "  Sub-tomogram classification:  iteration " << iter << " of " << prm.Niter << std::endl;

            prm.expectation(prm.SFi, 
                            prm.SFw, 
                            dataRefs,
                            dataSigma,
                            dataWsumRefs, 
                            dataWsumWedsPerRef,
                            dataWsumWedsPerGroup,
                            dataWsumDist,
                            dataSumWRefs,
                            imgsPmax,
                            imgsOptRefNos,
                            LL, 
                            avePmax);

            prm.maximization(dataRefs,
                             dataSigma,
                             dataWsumRefs,
                             dataWsumWedsPerRef,
                             dataWsumWedsPerGroup,
                             dataWsumDist,
                             dataSumWRefs,
                             imgsPmax,
                             imgsOptRefNos,
                             sumw_allrefs, 
                             avePmax );

            prm.writeOutputFiles(iter, 
                                 dataRefs,
                                 sumw_allrefs, 
                                 LL, 
                                 avePmax,
                                 conv);

        } // end loop iterations

        // Write out converged structures
        prm.writeOutputFiles(-1, 
                             dataRefs,
                             sumw_allrefs, 
                             LL, 
                             avePmax,
                             conv);

        // Free the memory (is this necessary?)
        delete [] dataSumWRefs;
        delete [] dataWsumRefs;
        delete [] dataWsumWedsPerRef;
        delete [] dataWsumWedsPerGroup;
        delete [] dataWsumDist;
        delete [] imgsPmax;
        delete [] imgsOptRefNos;
        delete [] dataRefs;
        delete [] dataSigma;

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(0);
    }
}



