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

/** Some constast to message passing tags */
#define TAG_SEED 1
#define TAG_DOCFILESIZE 2
#define TAG_DOCFILE 3

MpiProgML2D::~MpiProgML2D()
{
    delete node;
}

void MpiProgML2D::read(int argc, char** argv)
{
    node = new MpiNode(argc, argv);
    // Read subsequently to avoid problems in restart procedure
    for (int proc = 0; proc < node->size; ++proc)
    {
        if (proc == node->rank)
            ProgML2D::read(argc, argv);
        node->barrierWait();
    }

    //Send "master" seed to slaves for same randomization
    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!node->isMaster())
      verbose = 0;
}

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

void MpiProgML2D::addPartialDocfileData(const MultidimArray<double> &data, int first, int last)
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
        ProgML2D::addPartialDocfileData(docfiledata, myFirstImg, myLastImg);
        int s_size, first_img, last_img;
        MPI_Status status;

        for (int docCounter = 1; docCounter < node->size; ++docCounter)
        {
            // receive in order
            MPI_Recv(&s_size, 1, MPI_INT, docCounter, TAG_DOCFILESIZE,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(MULTIDIM_ARRAY(docfiledata), s_size,
                     MPI_DOUBLE, docCounter, TAG_DOCFILE,
                     MPI_COMM_WORLD, &status);
            divide_equally(nr_images_global, node->size, docCounter, first_img, last_img);
            ProgML2D::addPartialDocfileData(docfiledata, first_img, last_img);
        }
    }
}

void MpiProgML2D::writeOutputFiles(const ModelML2D &model, OutputType outputType)
{
    //All nodes should arrive to writeOutput files at same time
    node->barrierWait();
    //Only master write files
    if (node->isMaster())
        ProgML2D::writeOutputFiles(model, outputType);
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

