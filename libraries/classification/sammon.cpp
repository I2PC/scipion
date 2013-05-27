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
//-----------------------------------------------------------------------------
// Sammon.cc
// Sammon Projection Algorithm
//-----------------------------------------------------------------------------

#include "sammon.h"
#include <data/xmipp_funcs.h>

//-----------------------------------------------------------------------------
// Sammon: Sammon Maps
//-----------------------------------------------------------------------------

Sammon::Sammon(const unsigned _mapped,
                         const unsigned _num_iterations ,
                         const double _learning_rate):
        mapped(_mapped),
        num_iterations(_num_iterations),
        learning_rate(_learning_rate),
        stress(-1),
        verbosity(0)
{
    if (mapped < 1)
        throw std::invalid_argument("xmippSammon: mapped space must be > 1");
    if (num_iterations < 1)
        throw std::invalid_argument("xmippSammon: number of iterations must be > 0");
    if (learning_rate <= 0.0)
        throw std::invalid_argument("xmippSammon: learning rate must be > 0");
}

//-----------------------------------------------------------------------------

void Sammon::operator()(const In& in, Out& out)
{
    // clean the mapped space
    out.clear();
    out.calibrated(in.calibrated());
    // initialization of mapped space
    RandomUniformGenerator<double> uniform(-0.5, 0.5);
    FeatureVector v(mapped);
    unsigned i;
    for (i = 0; i < in.size(); i++)
    {
        generate(v.begin(), v.end(), uniform);
        transform(v.begin(), v.end(), v.begin(),
                  bind2nd(std::divides<double>(), norm(v)));
        out.add(v, in.theTargets[i]);
    }

    // calculate distances in original space
    std::vector<floatFeature> distances(in.size() *(in.size() - 1) / 2);
    std::vector<floatFeature>::iterator distance = distances.begin();
    for (i = 1; i < in.size(); i++)
        for (unsigned j = 0; j < i; j++)
            *distance++ = max(0.001, euclideanDistance(in.theItems[i], in.theItems[j]));

    // centroids of mapped samples
    std::vector<floatFeature> centroid(mapped);

    // first derivative and second derivative of mapping error
    std::vector<floatFeature> dE(in.size());
    std::vector<floatFeature> d2E2(in.size());

    // copy of the samples for each pattern loop p
    std::vector<std::vector<floatFeature> > out2(in.size(), std::vector<floatFeature>(mapped));
    int p;
    unsigned q;

    verbosity = listener->getVerbosity();
    if (verbosity)
        listener->OnReportOperation((std::string) "mapping....\n");
    if (verbosity == 1 || verbosity == 3)
        listener->OnInitOperation(num_iterations);


    for (unsigned iteration = 0; iteration < num_iterations; iteration++)
    {
        // reset centroids
        fill(centroid.begin(), centroid.end(), 0.0);

        // for each pattern ...
        for (p = 0; p < (int)in.size(); p++)
        {
            // reset error
            fill(dE.begin(), dE.end(), 0.0);
            fill(d2E2.begin(), d2E2.end(), 0.0);

            for (int j = 0; j < (int)in.size(); j++)
                if (j != p)
                {
                    unsigned m;
                    if (j < p)
                        m = p * (p - 1) / 2 + j;
                    else
                        m = j * (j - 1) / 2 + p;

                    // distance between two mapped samples
                    double old_dist = distances[m];
                    double new_dist = euclideanDistance(out.theItems[p], out.theItems[j]);
                    double dist_diff = old_dist - new_dist;
                    double dist_prod = old_dist * new_dist;

                    // calculate delta for each vector in mapped space
                    for (q = 0; q < mapped; q++)
                    {
                        double out_diff = out.theItems[p][q] - out.theItems[j][q];

                        if (dist_prod != 0.0)
                        {
                            // first order derivative
                            dE[q] += (dist_diff / dist_prod) * out_diff;
                            // second order derivative
                            d2E2[q] += (1.0 / dist_prod) * (dist_diff - (out_diff * out_diff / new_dist) * (1.0 + dist_diff / new_dist));
                        }
                    }
                }

            // adapt each patterm in the mapped space
            for (q = 0; q < mapped; q++)
                // adjust the q-th dimension of the p-th mapped vector
                centroid[q] += out2[p][q] =
                                   out.theItems[p][q] + learning_rate * dE[q] / fabs(d2E2[q]);
        }

        // shift the mapped vectors to their centroid
        for (p = 0; p < (int)out.size(); p++)
            for (q = 0; q < mapped; q++)
                out.theItems[p][q] = out2[p][q] - (centroid[q] / out.size());

        // shift the mapped vectors to their centroid
        for (p = 0; p < (int)out.size(); p++)
            for (q = 0; q < mapped; q++)
                out.theItems[p][q] = out2[p][q] - (centroid[q] / out.size());


        if (verbosity == 1 || verbosity == 3)
            listener->OnProgress(iteration);

        if ((verbosity > 1) || (iteration == num_iterations - 1))
        {
            // Calculates sammon stress (error in distances)

            stress = 0.;
            double tot = 0.;
            int mutual = 0;
            for (int pp = 0; pp < (int)in.size(); pp++)
                for (int jj = 0; jj < pp; jj++)
                {
                    double d = distances[mutual];
                    tot += d;
                    double ee = d - euclideanDistance(out.theItems[pp], out.theItems[jj]);
                    stress += (ee * ee / d);
                    mutual++;
                }
            stress /= tot;

            // Show it
            if (verbosity >= 2)
            {
                char s[100];
                sprintf(s, "Iteration %d of %d. Sammon Stress: %f\n", iteration + 1, num_iterations, stress);
                listener->OnReportOperation((std::string) s);
            }

        } // if verbosity

    }
    if (verbosity == 1 || verbosity == 3)
        listener->OnProgress(num_iterations);

}

