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
// SOM.cc
// Implements Kohonen Self-Organizing Feature Maps
//-----------------------------------------------------------------------------

#include "som.h"

/**
 * Construct a SOM from the code vectors in a stream
 * Parameter: _is  The stream
 */
SOM::SOM(std::istream& _is): ClassificationAlgorithm<ClassificationMap>()
{
    readSelf(_is);
}


/**
 * Sets the alpha function
 * Parameter: _alpha  alpha(t)
 */
void SOM::alpha(Descent _alpha)
{
    somAlpha = _alpha;
}

/**
 * Sets the radius function
 * Parameter: _radius  radius(t)
 */
void SOM::radius(Descent _radius)
{
    somRadius = _radius;
}

/**
 * Sets the number of training steps
 * Parameter: _nSteps  Number of training steps
 */
void SOM::nSteps(const unsigned long& _nSteps)
{
    somNSteps = _nSteps;
}


/**
 * Trains the SOM
 * Parameter: _som  The som to train
 * Parameter: _ts   The training set
 */
void SOM::train(ClassificationMap& _som, ClassicTrainingVectors& _ts) const
{
    unsigned long t = 0;

    int verbosity = listener->getVerbosity();
    if (verbosity)
        listener->OnReportOperation((std::string) "Training Kohonen SOM....\n");
    if (verbosity == 1 || verbosity == 3)
        listener->OnInitOperation(somNSteps);


    while (t < somNSteps)
    {
        for (std::vector<SomIn>::iterator i = _ts.theItems.begin();
             t < somNSteps && i < _ts.theItems.end() ; i++, t++)
        {
            // get the best matching.
            SomIn& theBest = _som.test(*i)
                             ;
            if (somNeigh == BUBBLE)
            { // Bubble
                // update the neighborhood around the best one
                std::vector<unsigned> neig = _som.neighborhood(_som.codVecPos(theBest),
                                             ceil(somRadius(t, somNSteps)));
                for (std::vector<unsigned>::iterator it = neig.begin();it < neig.end();it++)
                {
                    SomIn& v = _som.theItems[*it]
                               ;
                    for (unsigned j = 0; j < v.size(); j++)
                        v[j] += ((*i)[j] - v[j]) * somAlpha(t, somNSteps);
                }
            }
            else
            { // Gaussian
                // update all neighborhood convoluted by a gaussian
                double radius = somRadius(t, somNSteps);
                double alpha = somAlpha(t, somNSteps);
                for (unsigned it = 0 ; it < _som.size(); it++)
                {
                    double dist = _som.neighDist(_som.codVecPos(theBest), _som.indexToPos(it));
                    double alp = alpha * (double) exp((double)(-dist * dist / (2.0 * radius * radius)));
                    SomIn& v = _som.theItems[it];
                    for (unsigned j = 0; j < v.size(); j++)
                        v[j] += ((*i)[j] - v[j]) * alp;
                }
            } // else

        } // for examples

        if (verbosity == 1 || verbosity == 3)
            listener->OnProgress(t);
        if (verbosity >= 2)
        {
            char s[100]
            ;
            sprintf(s, "Iteration %d of %d.\n", (int)t, (int)somNSteps);
            listener->OnReportOperation((std::string) s);
        }
    } // while t < somSteps


    if (verbosity == 1 || verbosity == 3)
        listener->OnProgress(somNSteps);

}


/**
 * Tests the SOM
 * Parameter: _som        The som to test
 * Parameter: _examples   The training set of examples
 */
double SOM::test(const ClassificationMap& _som, const TS& _examples) const
{

    // Defines verbosity level
    int verbosity = listener->getVerbosity();
    if (verbosity)
    {
        listener->OnReportOperation((std::string) "Estimating quantization error....\n");
        listener->OnInitOperation(_examples.size());
    }


    /* Scan all data entries */
    double qerror = 0.0;
    for (size_t i = 0; i < _examples.size(); i++)
    {
        SomIn& theBest = _som.test(_examples.theItems[i]); // get the best
        qerror += euclideanDistance(theBest, _examples.theItems[i]);
        if (verbosity)
        {
            int tmp = (int)((_examples.size() * 5) / 100);
            if ((tmp == 0) && (i != 0))
                tmp = i;
            else
                tmp = 1;
            if ((i % tmp) == 0)
                listener->OnProgress(i);
        }
    }
    if (verbosity)
        listener->OnProgress(_examples.size());
    return (qerror / (double) _examples.size());
}

/**
* Clears the Algorithm
*/

void SOM::clear()
{
    somNeigh = GAUSSIAN;
    somNSteps = 0;
    listener = NULL; // it can not be deleted here
}

/**
* Standard output for a SOM algorithm
* Parameter: _os The output stream
*/
void SOM::printSelf(std::ostream& _os) const
{
    _os << (int) somNeigh << std::endl;
    _os << somNSteps << std::endl;
    somAlpha.printSelf(_os);
    somRadius.printSelf(_os);
}

/**
* Standard input for a som algorithm
* Parameter: _is The input stream
*/
void SOM::readSelf(std::istream& _is)
{
    clear();
    try
    {
        if (_is)
            _is >> (int&) somNeigh;
        if (_is)
            _is >> somNSteps;
        somAlpha.readSelf(_is);
        somRadius.readSelf(_is);

    }
    catch (std::exception& e)
    {
        std::ostrstream msg;
        msg << e.what() << std::endl << "Error reading the SOM algorithm";
        throw std::runtime_error(msg.str());
    }

}


/**
 * Saves the SOM class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void SOM::saveObject(std::ostream& _os) const
{
    printSelf(_os);
}


/**
 * Loads the SOM class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void SOM::loadObject(std::istream& _is)
{
    readSelf(_is);
}




//---------------------------------------------------------------------------
// Class Descent
//---------------------------------------------------------------------------
/**
 * Returns the function value associated a step if the transition from
 * the initial value to the final value es made in _nSteps steps
 * Parameter: _step    The actual step
 * Parameter: _nsteps  The number of steps to reach the final val from the
 *       initial one
 */
double Descent::operator()(const unsigned _step, const unsigned _nSteps) const
{
    if (_nSteps == 0 || initial == final || _step >= _nSteps)
        return final;
    else
        return final + ((initial - final) *
                        ((double)(_nSteps - _step) / (double)_nSteps));
}


/**
* Standard output for a Descent class
* Parameter: _os The output stream
*/
void Descent::printSelf(std::ostream& _os) const
{
    _os << initial << std::endl;
    _os << final << std::endl;
}

/**
* Standard input for a Descent class
* Parameter: _is The input stream
*/
void Descent::readSelf(std::istream& _is)
{
    try
    {
        if (_is)
            _is >> initial;
        if (_is)
            _is >> final;
    }
    catch (std::exception& e)
    {
        std::ostrstream msg;
        msg << e.what() << std::endl << "Error reading Descent class";
        throw std::runtime_error(msg.str());
    }

}


/**
 * Saves the Descent class into a stream.
 * this method can be used to save the status of the class.
 * Parameter: _os The output stream
 */
void Descent::saveObject(std::ostream& _os) const
{
    printSelf(_os);
}


/**
 * Loads the Descent class from a stream.
 * this method can be used to load the status of the class.
 * Parameter: _is The output stream
 */
void Descent::loadObject(std::istream& _is)
{
    readSelf(_is);
}

