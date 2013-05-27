/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

#include "mpi_ml_align2d.h"

/* Some constast to message passing tags */
#define TAG_SEED 1
#define TAG_DOCFILESIZE 2
#define TAG_DOCFILE 3

MpiML2DBase::MpiML2DBase(XmippProgram * prm)
{
    node = NULL;
    program = prm;
}

MpiML2DBase::MpiML2DBase(XmippProgram * prm, MpiNode * mpinode)
{
    node = mpinode;
    created_node = false;
    program = prm;
}

MpiML2DBase::~MpiML2DBase()
{
    if (created_node)
        delete node;
}

void MpiML2DBase::readMpi(int argc, char** argv)
{
    if (node == NULL)
    {
        node = new MpiNode(argc, argv);
        created_node = true;
    }
    //The following makes the asumption that 'this' also
    //inherits from an XmippProgram
    if (!node->isMaster())
        program->verbose = 0;
    // Read subsequently to avoid problems in restart procedure
    for (size_t proc = 0; proc < node->size; ++proc)
    {
        if (proc == node->rank)
            program->read(argc, (const char **)argv);
        node->barrierWait();
    }
    //Send "master" seed to slaves for same randomization
    MPI_Bcast(&program->seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void MpiML2DBase::sendDocfile(const MultidimArray<double> &docfiledata)
{

    // Write intermediate files
    if (!node->isMaster())
    {
        // All slaves send docfile data to the master
        int s_size = MULTIDIM_SIZE(docfiledata);
        MPI_Send(&s_size, 1, MPI_INT, 0, TAG_DOCFILESIZE,
                 MPI_COMM_WORLD);
        MPI_Send(MULTIDIM_ARRAY(docfiledata), s_size, MPI_DOUBLE,
                 0, TAG_DOCFILE, MPI_COMM_WORLD);
    }
    else
    {
        // Master fills docfile
        // Master's own contribution
        ML2DBaseProgram * ml2d = (ML2DBaseProgram*)program;
        ml2d->addPartialDocfileData(docfiledata, ml2d->myFirstImg, ml2d->myLastImg);
        int s_size;
        size_t first_img, last_img;
        MPI_Status status;

        for (size_t docCounter = 1; docCounter < node->size; ++docCounter)
        {
            // receive in order
            MPI_Recv(&s_size, 1, MPI_INT, docCounter, TAG_DOCFILESIZE,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(MULTIDIM_ARRAY(docfiledata), s_size,
                     MPI_DOUBLE, docCounter, TAG_DOCFILE,
                     MPI_COMM_WORLD, &status);
            divide_equally(ml2d->nr_images_global, node->size, docCounter, first_img, last_img);
            ml2d->addPartialDocfileData(docfiledata, first_img, last_img);
        }
    }
}

MpiProgML2D::MpiProgML2D():MpiML2DBase(this)
{}

MpiProgML2D::MpiProgML2D(MpiNode * node):MpiML2DBase(this, node)
{}



void MpiProgML2D::setNumberOfLocalImages()
{
    nr_images_local = divide_equally(nr_images_global, node->size, node->rank, myFirstImg,
                                     myLastImg);
}

void MpiProgML2D::produceSideInfo2()
{
    node->barrierWait();
    ProgML2D::produceSideInfo2();
    //Also sync after finishing produceSideInfo2
    node->barrierWait();
}

void MpiProgML2D::expectation()
{
    MultidimArray<double> Maux;
    double aux;
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();

    ProgML2D::expectation();
    //After expectation, collect data from all nodes
    // Here MPI_allreduce of all wsums,LL and sumfracweight !!!
    MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    LL = aux;
    MPI_Allreduce(&sumfracweight, &aux, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    sumfracweight = aux;
    MPI_Allreduce(&wsum_sigma_noise, &aux, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    wsum_sigma_noise = aux;
    MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
    wsum_sigma_offset = aux;
    for (int refno = 0; refno < model.n_ref * factor_nref; refno++)
    {
        MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[refno]),
                      MULTIDIM_ARRAY(Maux),
                      MULTIDIM_SIZE(wsum_Mref[refno]), MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        wsum_Mref[refno] = Maux;
        MPI_Allreduce(&sumw[refno], &aux, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        sumw[refno] = aux;
        MPI_Allreduce(&sumwsc2[refno], &aux, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        sumwsc2[refno] = aux;
        MPI_Allreduce(&sumw_mirror[refno], &aux, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        sumw_mirror[refno] = aux;
        MPI_Allreduce(&sumw2[refno], &aux, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        sumw2[refno] = aux;
        MPI_Allreduce(&sumwsc[refno], &aux, 1, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        sumwsc[refno] = aux;
    }
}//end of expectation

void MpiProgML2D::endIteration()
{
    // Write output files
    sendDocfile(docfiledata);
    writeOutputFiles(model, OUT_ITER);
}

void MpiProgML2D::writeOutputFiles(const ModelML2D &model, OutputType outputType)
{
    //All nodes should arrive to writeOutput files at same time
    node->barrierWait();
    //Only master write files
    if (node->isMaster())
        ProgML2D::writeOutputFiles(model, outputType);
    else if (outputType == OUT_REFS)
        outRefsMd = FN_CLASSES_MD(getIterExtraPath(fn_root, iter));
    //All nodes wait until files are written
    node->barrierWait();
}

//Just for debuging
void MpiProgML2D::printModel(const String &msg, const ModelML2D & model)
{
    if (node->isMaster())
        ProgML2D::printModel(msg, model);
}

void MpiProgML2D::usage(int verb) const
{
    if (node->isMaster())
        ProgML2D::usage();
}

#define SET_RANK_AND_SIZE() rank = node->rank; size = node->size;

MpiProgMLRefine3D::MpiProgMLRefine3D(int argc, char ** argv, bool fourier):MpiML2DBase(this)
{
    //create mpi node, which will be passed to ml2d
    node = new MpiNode(argc, argv);
    created_node = true;
    fourier_mode = fourier;

    if (fourier)
        ml2d = new MpiProgMLF2D(node);
    else
        ml2d = new MpiProgML2D(node);

    SET_RANK_AND_SIZE();
}

void MpiProgMLRefine3D::copyVolumes()
{
    //only master copy volumes before processing
    if (node->isMaster())
        ProgMLRefine3D::copyVolumes();
    //all nodes waiting until volumes are copied
    node->barrierWait();
}

void MpiProgMLRefine3D::reconstructVolumes()
{
    //code is already parallelized through the rank variable
    LOG("           MpiProgMLRefine3D::reconstructVolumes");
    ProgMLRefine3D::reconstructVolumes();
    //all nodes needs to wait until reconstruction is done
    node->barrierWait();
}

void MpiProgMLRefine3D::postProcessVolumes()
{
    //only master post process
    if (node->isMaster())
        ProgMLRefine3D::postProcessVolumes();
    //all nodes waiting until volumes are copied
    node->barrierWait();
}

void MpiProgMLRefine3D::makeNoiseImages()
{
    //only master post process
    if (node->isMaster())
        ProgMLRefine3D::makeNoiseImages();
    //all nodes waiting until volumes are copied
    node->barrierWait();
}

/// Calculate 3D SSNR, only master and broadcast result
void MpiProgMLRefine3D::calculate3DSSNR(MultidimArray<double> &spectral_signal)
{
    //only master calculate
    if (node->isMaster())
        ProgMLRefine3D::calculate3DSSNR(spectral_signal);
    MPI_Bcast(MULTIDIM_ARRAY(spectral_signal), MULTIDIM_SIZE(spectral_signal), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

/// Convergency check, only master and broadcast result
bool MpiProgMLRefine3D::checkConvergence()
{
    int result = -1;
    LOG("           MpiProgMLRefine3D::checkConvergence: inside");
    //only master check
    //if (node->isMaster())
    if (node->rank == 1)
    {
      LOG("           MpiProgMLRefine3D::checkConvergence: master checkConvergence");
        result = ProgMLRefine3D::checkConvergence() ? 1 : 0;
    }
    node->barrierWait();
    LOG("           MpiProgMLRefine3D::checkConvergence: before Bcast");
    MPI_Bcast(&result, 1, MPI_INT, 1, MPI_COMM_WORLD);

    LOG("           MpiProgMLRefine3D::checkConvergence: leaving...");
    return result == 1;


}

void MpiProgMLRefine3D::createEmptyFiles(int type)
{
    //only master create empty files
    if (node->isMaster())
        ProgMLRefine3D::createEmptyFiles(type);
    //all nodes waiting until volumes are projected
    node->barrierWait();
}

void MpiProgMLRefine3D::projectVolumes(MetaData &mdProj)
{
    ProgMLRefine3D::projectVolumes(mdProj);
    //all nodes waiting until volumes are projected
    node->barrierWait();
}

MpiProgMLRefine3D::~MpiProgMLRefine3D()
{}

MpiProgMLF2D::MpiProgMLF2D():MpiML2DBase(this)
{}

MpiProgMLF2D::MpiProgMLF2D(MpiNode * node):MpiML2DBase(this, node)
{}

void MpiProgMLF2D::expectation()
{
    MultidimArray<double> Maux, Vaux;
    double aux;
    Maux.resize(dim, dim);
    Maux.setXmippOrigin();

    ProgMLF2D::expectation();
    MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    LL = aux;
    MPI_Allreduce(&sumcorr, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sumcorr = aux;
    MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    wsum_sigma_offset = aux;
    if (do_kstest)
    {
        Vaux.resize(sumhist);
        MPI_Allreduce(MULTIDIM_ARRAY(sumhist), MULTIDIM_ARRAY(Vaux),
                      MULTIDIM_SIZE(sumhist), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vaux)
        DIRECT_MULTIDIM_ELEM(sumhist, n) = DIRECT_MULTIDIM_ELEM(Vaux, n);
        for (size_t ires = 0; ires < hdim; ires++)
        {
            Vaux.resize(sumhist);
            MPI_Allreduce(MULTIDIM_ARRAY(resolhist[ires]), MULTIDIM_ARRAY(Vaux),
                          MULTIDIM_SIZE(resolhist[ires]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vaux)
            DIRECT_MULTIDIM_ELEM(resolhist[ires], n) = DIRECT_MULTIDIM_ELEM(Vaux, n);
        }
    }
    for (int refno = 0;refno < model.n_ref; refno++)
    {
        MPI_Allreduce(MULTIDIM_ARRAY(wsum_Mref[refno]), MULTIDIM_ARRAY(Maux),
                      MULTIDIM_SIZE(wsum_Mref[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        wsum_Mref[refno] = Maux;
        if (do_ctf_correction)
        {
            MPI_Allreduce(MULTIDIM_ARRAY(wsum_ctfMref[refno]), MULTIDIM_ARRAY(Maux),
                          MULTIDIM_SIZE(wsum_ctfMref[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsum_ctfMref[refno] = Maux;
        }
        MPI_Allreduce(&sumw[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sumw[refno] = aux;
        MPI_Allreduce(&sumw2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sumw2[refno] = aux;
        MPI_Allreduce(&sumwsc2[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sumwsc2[refno] = aux;
        MPI_Allreduce(&sumwsc[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sumwsc[refno] = aux;
        MPI_Allreduce(&sumw_mirror[refno], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sumw_mirror[refno] = aux;
    }
    for (size_t ifocus = 0;ifocus < nr_focus;ifocus++)
    {
        for (size_t ii = 0; ii <  Mwsum_sigma2[ifocus].size(); ii++)
        {
            MPI_Allreduce(&Mwsum_sigma2[ifocus][ii], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            Mwsum_sigma2[ifocus][ii] = aux;
        }
        if (do_student)
        {
            MPI_Allreduce(&sumw_defocus[ifocus], &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sumw_defocus[ifocus] = aux;
        }
    }
}//end of expectation

void MpiProgMLF2D::produceSideInfo2()
{
    node->barrierWait();
    ProgMLF2D::produceSideInfo2();
    //Also sync after finishing produceSideInfo2
    node->barrierWait();
}

void MpiProgMLF2D::produceSideInfo()
{
    SET_RANK_AND_SIZE();
    ProgMLF2D::produceSideInfo();
}


void MpiProgMLF2D::endIteration()
{
    // Write output files
    LOG("MpiProgMLF2D::endIteration : before sendDocfile");
    sendDocfile(docfiledata);
    LOG("MpiProgMLF2D::endIteration : before writeOutputFiles");
    writeOutputFiles(model, OUT_ITER);
    LOG("MpiProgMLF2D::endIteration : before updateWienerFilters");
    updateWienerFilters(spectral_signal, sumw_defocus, iter);
    LOG("MpiProgMLF2D::endIteration : before barrierWait");
    node->barrierWait();
}

void MpiProgMLF2D::writeOutputFiles(const ModelML2D &model, OutputType outputType)
{
    //All nodes should arrive to writeOutput files at same time
    LOG("MpiProgMLF2D::writeOutputFiles : waiting before writing");
    node->barrierWait();
    //Only master write files
    if (node->isMaster())
    {
        ProgMLF2D::writeOutputFiles(model, outputType);
        LOG("MpiProgMLF2D::writeOutputFiles : master writing ");
    }
    else if (outputType == OUT_REFS)
    {
      outRefsMd = FN_CLASSES_MD(getIterExtraPath(fn_root, iter));
      LOG(formatString("MpiProgMLF2D::writeOutputFiles : slave setting outRefsMd = %s", outRefsMd.c_str()).c_str());
    }
    //All nodes wait until files are written
    LOG("MpiProgMLF2D::writeOutputFiles : waiting after writing");
    node->barrierWait();
}


/** Constructor */
MpiProgMLTomo::MpiProgMLTomo()
{
    node = NULL;
}
/** Destructor */
MpiProgMLTomo::~MpiProgMLTomo()
{
    delete node;
}

/** Redefine the basic Program read to do it sequentially */
void MpiProgMLTomo::read(int argc, char ** argv, bool reportErrors)
{
    if (node == NULL)
        node = new MpiNode(argc, argv);
    //The following makes the asumption that 'this' also
    //inherits from an XmippProgram
    if (!node->isMaster())
        verbose = 0;
    // Read subsequently to avoid problems in restart procedure
    for (size_t proc = 0; proc < node->size; ++proc)
    {
        if (proc == node->rank)
            ProgMLTomo::read(argc, (const char **)argv);
        node->barrierWait();
    }
}
/// Only master will generate initial references
void MpiProgMLTomo::setNumberOfLocalImages()
{
    nr_images_local = divide_equally(nr_images_global, node->size, node->rank, myFirstImg, myLastImg);
}

/// Only master will generate initial references
void MpiProgMLTomo::generateInitialReferences()
{
    if (node->isMaster())
        ProgMLTomo::generateInitialReferences();
    else
        fn_ref = FN_ITER_MD(0);
    node->barrierWait();
}

/// Integrate over all experimental images
void MpiProgMLTomo::expectation(MetaData &MDimg, std::vector< Image<double> > &Iref, int iter,
                                double &LL, double &sumfracweight,
                                std::vector<MultidimArray<double> > &wsumimgs,
                                std::vector<MultidimArray<double> > &wsumweds,
                                double &wsum_sigma_noise, double &wsum_sigma_offset,
                                MultidimArray<double> &sumw)
{
    ProgMLTomo::expectation(MDimg,Iref,iter,LL,sumfracweight,wsumimgs,wsumweds,wsum_sigma_noise,wsum_sigma_offset,sumw);

    double aux;
    MultidimArray<double> Vaux; //1D
    MultidimArray<double> Maux, Maux2; //3D

    Maux.resize(dim, dim, dim);
    Maux.setXmippOrigin();
    Maux2.resize(dim, dim, hdim+1);

    // Here MPI_allreduce of all wsums,LL and sumcorr !!!
    MPI_Allreduce(&LL, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    LL = aux;
    MPI_Allreduce(&sumfracweight, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sumfracweight = aux;
    MPI_Allreduce(&wsum_sigma_noise, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    wsum_sigma_noise = aux;
    MPI_Allreduce(&wsum_sigma_offset, &aux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    wsum_sigma_offset = aux;
    Vaux.resize(nr_ref);
    MPI_Allreduce(MULTIDIM_ARRAY(sumw), MULTIDIM_ARRAY(Vaux),
                  MULTIDIM_SIZE(sumw), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sumw = Vaux;
    for (int refno = 0; refno < 2*nr_ref; ++refno)
    {
        MPI_Allreduce(MULTIDIM_ARRAY(wsumimgs[refno]), MULTIDIM_ARRAY(Maux),
                      MULTIDIM_SIZE(wsumimgs[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        wsumimgs[refno] = Maux;
        if (do_missing)
        {
            MPI_Allreduce(MULTIDIM_ARRAY(wsumweds[refno]), MULTIDIM_ARRAY(Maux2),
                          MULTIDIM_SIZE(wsumweds[refno]), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            wsumweds[refno] = Maux2;
        }
    }
}

///Add info of some processed images to later write to files
void MpiProgMLTomo::addPartialDocfileData(const MultidimArray<double> &data, size_t first, size_t last)
{
    // Write intermediate files
    if (!node->isMaster())
    {
        // All slaves send docfile data to the master
        int s_size = MULTIDIM_SIZE(docfiledata);
        MPI_Send(&s_size, 1, MPI_INT, 0, TAG_DOCFILESIZE, MPI_COMM_WORLD);
        MPI_Send(MULTIDIM_ARRAY(docfiledata), s_size, MPI_DOUBLE, 0, TAG_DOCFILE, MPI_COMM_WORLD);
    }
    else
    {
        // Master fills metadata and add it's contribution
        ProgMLTomo::addPartialDocfileData(docfiledata, myFirstImg, myLastImg);
        int s_size;
        size_t first_img, last_img;
        MPI_Status status;

        for (size_t docCounter = 1; docCounter < node->size; ++docCounter)
        {
            // receive in order
            MPI_Recv(&s_size, 1, MPI_INT, docCounter, TAG_DOCFILESIZE, MPI_COMM_WORLD, &status);
            MPI_Recv(MULTIDIM_ARRAY(docfiledata), s_size, MPI_DOUBLE, docCounter, TAG_DOCFILE, MPI_COMM_WORLD, &status);
            divide_equally(nr_images_global, node->size, docCounter, first_img, last_img);
            ProgMLTomo::addPartialDocfileData(docfiledata, first_img, last_img);
        }
    }
}


/// Only master write output files
void MpiProgMLTomo::writeOutputFiles(const int iter,
                                     std::vector<MultidimArray<double> > &wsumweds,
                                     double &sumw_allrefs, double &LL, double &avefracweight,
                                     std::vector<double> &conv, std::vector<MultidimArray<double> > &fsc)
{

    if (node->isMaster())
        ProgMLTomo::writeOutputFiles(iter, wsumweds, sumw_allrefs, LL, avefracweight, conv, fsc);
    node->barrierWait();
}



