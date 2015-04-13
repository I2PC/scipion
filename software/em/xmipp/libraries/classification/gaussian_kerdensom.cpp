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
// GaussianKerDenSOM.cc
// Implements Smoothly Distributed Kernel Probability Density Estimator Self-Organizing Map
// Uses a Gaussian Kernel Function.
//-----------------------------------------------------------------------------

#include <fstream>
#include <ctime>

#include "gaussian_kerdensom.h"
#include <data/metadata.h>

#define  MAXZ -11282

//-----------------------------------------------------------------------------
/**
 * Trains the GaussianDenSOM
 * Parameter: _som  The KerDenSom to train
 * Parameter: _ts   The training set
 * Parameter: _update True if uses _som as starting point for training.
 * Parameter: _sigma If update = true, uses this sigma for the training.
 */
void GaussianKerDenSOM::train(FuzzyMap& _som, TS& _examples, FileName& _fn,
                              bool _update, double _sigma, bool _saveIntermediate)
{
	MetaData MDconvergence;
	FileName tmpN;
    numNeurons = _som.size();
    numVectors = _examples.size();
    dim = _examples.theItems[0].size();
    tmpV.resize(dim, 0.);
    tmpDens.resize(numNeurons, 0.);
    tmpMap.resize(numNeurons);
    for (size_t i = 0; i < numNeurons; i++)
        tmpMap[i].resize(dim, 0.);
    tmpD.resize(numNeurons);
    tmpD1.resize(numNeurons);
    double stopError;

    int verbosity = listener->getVerbosity();
    if (verbosity)
        listener->OnReportOperation((std::string) "\nTraining KerDenSOM....\n");

    // Initialization
    if (verbosity)
        listener->OnReportOperation((std::string) "\nInitializing....\n");
    if (!_update)
    {
        initU(&_som);
        updateV1(&_som, &_examples);
        sigma = updateSigmaI(&_som, &_examples);
    }
    else
    {
        double alpha;
        sigma = _sigma;
        if (_sigma == 0.0)
        {
            updateU1(&_som, &_examples);
            sigma = updateSigmaI(&_som, &_examples);
        }
        updateU(&_som, &_examples, sigma, alpha);
    }


    double regtmp, tmpregMax, tmpregMin, pen, lkhood;
    tmpregMax = log(reg0);
    tmpregMin = log(reg1);

    if (annSteps == 0)
        annSteps = 1;
    for (size_t iter = 0; iter < annSteps; iter++)
    {

        if (verbosity)
        {
            if (annSteps > 1)
            {
                char s[100];
                sprintf(s, "\nTraining Deterministic Annealing step %d of %d....\n", (int)(iter + 1), (int)annSteps);
                listener->OnReportOperation((std::string) s);
            }
            else
                listener->OnReportOperation((std::string) "Training ....\n");
        }

        if (annSteps == 1)
            regtmp = reg1;
        else
            regtmp = exp(tmpregMax - iter * (tmpregMax - tmpregMin) / (annSteps - 1));
        stopError = mainIterations(&_som, &_examples, sigma, regtmp);
        if (verbosity)
            listener->OnReportOperation((std::string) "Calculating cost function....\n");
        double funct = functional(&_examples, &_som, sigma, reg1, lkhood, pen);

        size_t id=MDconvergence.addObject();
        MDconvergence.setValue(MDL_KERDENSOM_REGULARIZATION,regtmp,id);
        MDconvergence.setValue(MDL_KERDENSOM_FUNCTIONAL,funct,id);
        MDconvergence.setValue(MDL_KERDENSOM_SIGMA,sigma,id);

        if (verbosity)
        {
            char s[100];
            sprintf(s, "Code vectors variation: %g\n", stopError);
            listener->OnReportOperation((std::string) s);
        }


        if (annSteps > 1)
        {
            // Classifying
            if (verbosity)
                listener->OnReportOperation((std::string) "Classifying....\n");
            _som.classify(&_examples);

            // Calibrating
            if (verbosity)
                listener->OnReportOperation((std::string) "Calibrating....\n");
            _som.calibrate(_examples);

            if (_examples.isNormalized())
            {
                if (verbosity)
                    listener->OnReportOperation((std::string) "Denormalizing code vectors....\n");
                _som.unNormalize(_examples.getNormalizationInfo()); // de-normalize codevectors
            }

            if (_saveIntermediate)
            {
                // Saves each codebook (for all iterations)
                tmpN = _fn.c_str() + (std::string) "_" + integerToString(iter) + (std::string) ".cod";
                if (verbosity)
                    listener->OnReportOperation((std::string) "Saving code vectors....\n");
                std::ofstream codS(tmpN.c_str());
                codS << _som;
                codS.flush();

                // save .his file (Histogram)
                if (verbosity)
                    listener->OnReportOperation((std::string) "Saving histogram file....\n");
                tmpN = _fn.c_str() + (std::string) "_" + integerToString(iter) + (std::string) ".his";
                std::ofstream hisStream(tmpN.c_str());
                _som.printHistogram(hisStream);
                hisStream.flush();

                // save .err file (Average Quantization Error)
                if (verbosity)
                    listener->OnReportOperation((std::string) "Saving Average Quantization Error file....\n");
                tmpN = _fn.c_str() + (std::string) "_" + integerToString(iter) + (std::string) ".err";
                std::ofstream errStream(tmpN.c_str());
                _som.printQuantError(errStream);
                errStream.flush();

                // save .vs file to be compatible with SOM_PAK
                if (verbosity)
                    listener->OnReportOperation((std::string) "Saving visual file....\n");
                tmpN = _fn.c_str() + (std::string) "_" + integerToString(iter) + (std::string) ".vs";
                std::ofstream vsStream(tmpN.c_str());
                vsStream << _examples.theItems[0].size() << " " << _som.layout() << " " << _som.width() << " " << _som.height() << " gaussian" << std::endl;
                for (size_t i = 0; i < _examples.size(); i++)
                {
                    int j = _som.fuzzyWinner(i);
                    vsStream << _som.indexToPos(j).first << " " << _som.indexToPos(j).second << " " << _som.memb[i][j] << " " << _examples.theTargets[i] << std::endl;
                }
                vsStream.flush();

                // save .inf file
                if (verbosity)
                    listener->OnReportOperation((std::string) "Saving inf file....\n");
                tmpN = _fn.c_str() + (std::string) "_" + integerToString(iter) + (std::string) ".inf";
                std::ofstream infS(tmpN.c_str());
                infS << "Kernel Probability Density Estimator SOM algorithm" << std::endl << std::endl;
                infS << "Deterministic annealing step " << iter + 1 << " out of " << annSteps << std::endl;
                infS << "Number of feature vectors: " << _examples.size() << std::endl;
                infS << "Number of variables: " << _examples.theItems[0].size() << std::endl;
                infS << "Horizontal dimension (Xdim) = " << _som.width() << std::endl;
                infS << "Vertical dimension (Ydim) = " << _som.height() << std::endl;
                if (_examples.isNormalized())
                    infS << "Input data normalized" << std::endl;
                else
                    infS << "Input data not normalized" << std::endl;
                if (_som.layout() == "rect")
                    infS << "Rectangular topology " << std::endl;
                else
                    infS << "Hexagonal topology " << std::endl;
                infS << "Gaussian Kernel function " << std::endl;
                infS << "Total number of iterations = " << somNSteps << std::endl;
                infS << "Stopping criteria (eps) = " << epsilon << std::endl << std::endl ;

                infS << "Smoothness factor (regularization) = " << regtmp << std::endl;
                infS << "Functional value = " << funct << std::endl;
                infS << "Sigma (Kernel width) = " << sigma << std::endl;
                infS.flush();
            }

            if (_examples.isNormalized())
            {
                if (verbosity)
                    listener->OnReportOperation((std::string) "Normalizing code vectors....\n");
                _som.Normalize(_examples.getNormalizationInfo());       // normalize code vectors
            }
        } // if annSteps

    } // for
    MDconvergence.write(formatString("KerDenSOM_Convergence@%s",_fn.c_str()),MD_APPEND);

    tmpV.clear();
    tmpDens.clear();
    tmpMap.clear();
    tmpD.clear();
    tmpD1.clear();

}

//-----------------------------------------------------------------------------
/**
 * Update the U (Membership)
 */
double GaussianKerDenSOM::updateU(FuzzyMap* _som, const TS* _examples,
		                          const double& _sigma, double& _alpha)
{
    // Create auxiliar stuff
    double auxDist;
    double rr2, max1, d1, tmp, r1;

    double irr1 =1.0/( 2.0 * _sigma);
    _alpha = 0;

    // Update Membership matrix
    double *ptrTmpD=&tmpD[0];
    double *ptrTmpD1=&tmpD1[0];
    double idim=1.0/dim;
    for (size_t k = 0; k < numVectors; k++)
    {
        max1 = -MAXFLOAT;
        const floatFeature *ptrExample=&(_examples->theItems[k][0]);
        for (size_t i = 0; i < numNeurons; i ++)
        {
            auxDist = 0;
            const floatFeature *ptrCodeVector=&(_som->theItems[i][0]);
            for (size_t j = 0; j < dim; j++)
            {
                double tmp=((double)ptrExample[j] - (double)ptrCodeVector[j]);
                auxDist += tmp * tmp;
            }
            auxDist*=idim;
            ptrTmpD[i] = auxDist;
            rr2 = -auxDist * irr1;
            ptrTmpD1[i] = rr2;
            if (max1 < rr2)
                max1 = rr2;
        }
        r1 = 0;
        for (size_t j = 0; j < numNeurons; j ++)
        {
            rr2 = ptrTmpD1[j] - max1;
            if (rr2 < MAXZ)
                d1 = 0;
            else
                d1 = (double)exp(rr2);
            r1 += d1;
            ptrTmpD1[j] = d1;
        }
        double ir1=1.0/r1;

        floatFeature *ptrSomMembK=&(_som->memb[k][0]);
        for (size_t j = 0; j < numNeurons; j ++)
        {
            tmp = ptrTmpD1[j] * ir1;
            ptrSomMembK[j] = (floatFeature) tmp;
            _alpha += tmp * ptrTmpD[j];
        }
    } // for k
    return 0.0;
}

//-----------------------------------------------------------------------------

// Estimate Sigma part II
double GaussianKerDenSOM::updateSigmaII(FuzzyMap* _som, const TS* _examples, const double& _reg, const double& _alpha)
{
    size_t cc, j;

    if (_reg == 0)
        return(_alpha / (double)(numVectors*dim));

    // Computing Sigma (Part II)
    double q = 0.;
    for (cc = 0; cc < numNeurons; cc++)
    {
        _som->localAve(_som->indexToPos(cc), tmpV);
        double t1 = (double) _som->getLayout().numNeig(_som, (SomPos) _som->indexToPos(cc));
        for (j = 0; j < dim; j++)
            q += (t1 * ((double)_som->theItems[cc][j]) - (double)(tmpV[j])) * ((double)_som->theItems[cc][j]);
    }
    return (double)((_alpha + _reg*q) / (double)(numVectors*dim));
}


//-----------------------------------------------------------------------------

/**
 * Estimate the PD (Method 1: Using the code vectors)
 */
double GaussianKerDenSOM::codeDens(const FuzzyMap* _som, const FeatureVector* _example, double _sigma) const
{
    double s = 0;
    double K=-1.0/(2*_sigma);
    const floatFeature *ptrExample=&((*_example)[0]);
    for (size_t cc = 0; cc < numNeurons; cc++)
    {
        double t = 0;
        const floatFeature *ptrCodeVector=&(_som->theItems[cc][0]);
        for (size_t j = 0; j < dim; j++)
        {
            double diff=(double)(ptrExample[j]) - (double)(ptrCodeVector[j]);
            t += diff*diff;
        }
        t *= K;
        if (t < MAXZ)
            t = 0;
        else
            t = exp(t);
        s += t;
    }
    return std::pow(2*PI*_sigma, -0.5*dim)*s / numNeurons;
}

//-----------------------------------------------------------------------------

/**
 * Estimate the PD (Method 2: Using the data)
 */
double GaussianKerDenSOM::dataDens(const TS* _examples, const FeatureVector* _example, double _sigma) const
{
    unsigned j, vv;
    double t = 0, s = 0;
    for (vv = 0; vv < numVectors; vv++)
    {
        t = 0;
        for (j = 0; j < dim; j++)
            t += ((double)((*_example)[j]) - (double)(_examples->theItems[vv][j])) * ((double)((*_example)[j]) - (double)(_examples->theItems[vv][j]));
        t = -t / (2.0 * _sigma);
        if (t < MAXZ)
            t = 0;
        else
            t = exp(t);
        s += t;
    }
    double tmp = (double) dim / 2.0;
    return (double)(pow((double)(2*PI*_sigma), -tmp)*s / (double)numVectors);
}


//-----------------------------------------------------------------------------

/**
 * Determines the functional value
 * Returns the likelihood and penalty parts of the functional
 */
double GaussianKerDenSOM::functional(const TS* _examples, const FuzzyMap* _som,
		                             double _sigma, double _reg, double& _likelihood,
		                             double& _penalty)
{
    unsigned j, vv, cc;
    double t;
    _likelihood = 0;
    for (vv = 0; vv < numVectors; vv++)
    {
        t = codeDens(_som, &(_examples->theItems[vv]), _sigma);
        if (t == 0)
        {
            t = 1e-300;
        }
        _likelihood += log(t);
    }
    _likelihood = -_likelihood;
    _penalty = 0;

    if (_reg != 0)
    {
        for (cc = 0; cc < numNeurons; cc++)
        {
            _som->localAve(_som->indexToPos(cc), tmpV);
            t = _som->getLayout().numNeig(_som, (SomPos) _som->indexToPos(cc));
            for (j = 0; j < dim; j++)
            {
                _penalty += (t * (double)(_som->theItems[cc][j]) -
                		         (double)(tmpV[j])) * (double)(_som->theItems[cc][j]);
            }
        }
        _penalty = _reg * _penalty / (2.0 * _sigma);
    }
    return (_likelihood + _penalty);
}

//-----------------------------------------------------------------------------
