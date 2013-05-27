/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
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

/*-----------------------------------------------------------------------------
 BatchSOM.cc
 Implements Kohonen Self-Organizing Feature Maps by using Batch SOM
-----------------------------------------------------------------------------*/

#include "batch_som.h"

/**
 * Construct a BatchSOM from the code vectors in a stream
 * Parameter: _is  The stream
 */
BatchSOM::BatchSOM(std::istream& _is): SOM(_is)
{
    readSelf(_is);
}

/**
 * Trains the SOM
 * Parameter: _som  The som to train
 * Parameter: _ts   The training set
 */
void BatchSOM::train(ClassificationMap& _som, const ClassicTrainingVectors& _ts) const
{


    unsigned long t = 0;

    int verbosity = listener->getVerbosity();
    if (verbosity)
        listener->OnReportOperation((std::string) "Batch Training Kohonen SOM....\n");
    if (verbosity == 1 || verbosity == 3)
        listener->OnInitOperation(somNSteps);

    SomIn aveVector(_som.theItems[0].size());
    std::vector<unsigned> tmpVector;

    while (t < somNSteps)
    {
        _som.classify(&_ts);
        // Check for each SOM unit
        for (unsigned it = 0; it < _som.size(); it++)
        {
            for (unsigned a = 0; a < aveVector.size(); a++)
                aveVector[a] = 0.;
            long total = 0;
            // Collects a list of the input vectors assigned to the neighborhood
            std::vector<unsigned> neig = _som.neighborhood(_som.indexToPos(it), ceil(somRadius(t, somNSteps)));
            for (std::vector<unsigned>::iterator itt = neig.begin();itt < neig.end();itt++)
            {
                tmpVector =  _som.classifAt(*itt);
                for (unsigned j = 0 ; j < tmpVector.size() ; j++)
                {
                    SomIn v = _ts.theItems[tmpVector[j]];
                    for (unsigned a = 0; a < v.size(); a++)
                        aveVector[a] += v[a];
                    total++;
                }

            }
            if (total != 0)
            {
                for (unsigned a = 0; a < aveVector.size(); a++)
                    aveVector[a] /= (floatFeature) total;
                _som.theItems[it] = aveVector;
            }
        }

        if (verbosity == 1 || verbosity == 3)
            listener->OnProgress(t);
        if (verbosity >= 2)
        {
            char s[100];
            sprintf(s, "Iteration %d of %d.\n", (int)(t + 1), (int)somNSteps);
            listener->OnReportOperation((std::string) s);
        }
        t++;

    } // while t < somSteps

    if (verbosity == 1 || verbosity == 3)
        listener->OnProgress(somNSteps);
}


