/***************************************************************************
 *
 * Authors:     Enrique Recarte Llorens   (erecallo@hotmail.com)
 *              Carlos Oscar S. Sorzano   (coss@cnb.csic.es)
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

#include "naive_bayes.h"
#include <string>
#include <queue>

bool debugging = true;
//#define DEBUG_SPLITTING_USING_ENTROPY
//#define DEBUG_WEIGHTS
//#define DEBUG_CLASSIFICATION
//#define DEBUG_FINE_CLASSIFICATION

/* ------------------------------------------------------------------------- */
/* Leaf nodes                                                                */
/* ------------------------------------------------------------------------- */
/* Split histograms -------------------------------------------------------- */
// Split several histograms within the indexes l0 and lF so that
// the entropy after division is maximized
int splitHistogramsUsingEntropy(const std::vector<histogram1D> &hist,
    int l0, int lF)
{
    // Number of classes
    int K = hist.size();
    
    // Set everything outside l0 and lF to zero, and make it a PDF
    std::vector<histogram1D> histNorm;
    for (int k = 0; k < K; k++)
    {
        histogram1D histaux = hist[k];
        for (int l = 0; l < XSIZE(histaux); l++)
            if (l < l0 || l > lF) 
                histaux(l) = 0;
        histaux /= histaux.sum();
        histNorm.push_back(histaux);
    }

    // Compute for each class the probability of being l<=l0 and l>l0
    Matrix2D<double> p(K, 2);
    for (int k = 0; k < K; k++)
    {
        p(k, 0) = histNorm[k](l0);
        p(k, 1) = 0;
        for (int l = l0 + 1; l <= lF; l++)
            p(k, 1) += histNorm[k](l);
    }

    // Compute the splitting l giving maximum entropy
    double maxEntropy = 0;
    int lmaxEntropy = -1;
    int l = l0;
    while (l < lF)
    {
        // Compute the entropy of the clases if we split by l
        double entropy = 0;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(p)
            if (p(i, j) != 0) entropy -= p(i, j) * log10(p(i, j));

        #ifdef DEBUG_SPLITTING_USING_ENTROPY
            std::cout << "Splitting at "  << l << " entropy=" << entropy
                      << std::endl;
        #endif

        // Check if this is the maximum
        if (entropy > maxEntropy)
        {
            maxEntropy = entropy;
            lmaxEntropy = l;
        }
         
        // Move to next split point
        l++;
        
        // Update probabilities of being l<=l0 and l>l0
        for (int k = 0; k < K; k++)
        {
            p(k, 0) += histNorm[k](l);
            p(k, 1) -= histNorm[k](l);
        }
    }
    
    #ifdef DEBUG_SPLITTING_USING_ENTROPY
        std::cout << "Finally in l=[" << l0 << "," << lF
                  << " Max Entropy:" << maxEntropy
                  << " lmax=" << lmaxEntropy << std::endl;
    #endif

    // If the point giving the maximum entropy is too much on the extreme,
    // substitute it by the middle point
    if (lmaxEntropy<=2 || lmaxEntropy>=lF-2)
        lmaxEntropy = (int) CEIL((lF + l0)/2.0);

    return lmaxEntropy;
}

/* Constructor ------------------------------------------------------------- */
LeafNode::LeafNode(const std::vector < Matrix1D<double> > &leafFeatures, 
    int discrete_levels)
{
    __discreteLevels = discrete_levels;
    K = leafFeatures.size();

    // Compute the minimum and maximum of each class
    double minval, maxval;
    for(int k=0; k<K; k++)
    {
        double minvalk, maxvalk;
        leafFeatures[k].computeDoubleMinMax(minvalk, maxvalk);
        if (k==0)
        {
            minval=minvalk;
            maxval=maxvalk;
        }
        else
        {
            minval=XMIPP_MIN(minval,minvalk);
            maxval=XMIPP_MAX(maxval,maxvalk);
        }
    }

    // Compute the PDF of each class
    std::vector<histogram1D> hist(K);
    for (int k=0; k<K; k++)
    {
        compute_hist(leafFeatures[k], hist[k], minval, maxval, 100);
        hist[k] += 1; // Apply Laplace correction
    }

    // Split the histograms into discrete_level (power of 2) bins
    std::queue< Matrix1D<int> > intervals, splittedIntervals;
    Matrix1D<int> limits(2);
    VECTOR_R2(limits,0,99);
    intervals.push(limits);
    for (int i=0; i<ROUND(log2(__discreteLevels)); i++)
    {
        // Split all the intervals in the queue
        while (!intervals.empty())
        {
            Matrix1D<int> currentInterval = intervals.front();
            intervals.pop();
            int lsplit = splitHistogramsUsingEntropy(hist,
                currentInterval(0), currentInterval(1));
            VECTOR_R2(limits,currentInterval(0),lsplit);
            splittedIntervals.push(limits);
            VECTOR_R2(limits,lsplit+1, currentInterval(1));
            splittedIntervals.push(limits);
        }
        
        // Copy the splitted intervals to the interval list
        while (!splittedIntervals.empty())
        {
            intervals.push(splittedIntervals.front());
            splittedIntervals.pop();
        }
    }

    // Compute the bins of the split
    Matrix1D<int> newBins(__discreteLevels);
    int imax=intervals.size();
    for (int i=0; i<imax; i++)
    {
        newBins(i) = intervals.front()(1);
        intervals.pop(); 
    }

    // Compute now the irregular histograms
    for (int k=0; k<K; k++)
    {
        IrregularHistogram1D irregHist;
        irregHist.init(hist[k], newBins);
        irregHist.selfNormalize();
        __leafPDF.push_back(irregHist);    
    }
}

/* Get the number of levels ------------------------------------------------ */
int LeafNode::getNumberOfLevels()
{
    return __discreteLevels;
}

/* Assign probability ------------------------------------------------------ */
double LeafNode::assignProbability(double value, int k)
{
    int index = __leafPDF[k].val2Index(value);
    return __leafPDF[k](index);
}

/* Compute weight ---------------------------------------------------------- */
double LeafNode::computeWeight() const
{
    double retval=0;
    for (int k1=0; k1<K; k1++)
        for (int k2=0; k2<K; k2++)
        {
            if (k1==k2) continue;
            retval+=KLDistance(
                __leafPDF[k1].getHistogram(),
                __leafPDF[k2].getHistogram());
        }
    return (retval/(K*K-K));
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out, const LeafNode &leaf)
{
    for (int k=0; k<leaf.K; k++)
        _out << "Histogram of class " << k << "=\n"
             << leaf.__leafPDF[k] << std::endl;
    _out << "Classification power=" << leaf.computeWeight() << std::endl;
    return _out;
}

/* ------------------------------------------------------------------------- */
/* Naive Bayes classifier                                                    */
/* ------------------------------------------------------------------------- */
NaiveBayes::NaiveBayes(
    const std::vector< Matrix2D<double> > &features,
    const Matrix1D<double> &priorProbs,
    int discreteLevels)
{ 
    K = features.size();
    Nfeatures=XSIZE(features[0]);
    __priorProbsLog10 = priorProbs;
    __priorProbsLog10.selfLog10();
    __priorProbsLog10.resize(K);

    // Build a leafnode for each feature and assign a weight
    __weights.initZeros(Nfeatures);
    std::vector < Matrix1D<double> > aux(K);	
    for (int f=0; f<Nfeatures; f++)
    {
        for (int k=0; k<K; k++)
            features[k].getCol(f, aux[k]);
        __leafs.push_back(new LeafNode(aux,discreteLevels));
        __weights(f)=__leafs[f]->computeWeight();
        #ifdef DEBUG_WEIGHTS
            if(debugging == true)
            {
                std::cout << "Node " << f << std::endl
                          << *(__leafs[f]) << std::endl;
            }
        #endif            
    }
    __weights /= __weights.computeMax();
    
    // Set default cost matrix
    __cost.resize(K,K);
    __cost.initConstant(1);
    for (int i=0; i<K; i++) __cost(i,i)=0;
}

/* Destructor -------------------------------------------------------------- */
NaiveBayes::~NaiveBayes()
{
   for (int i = 0; i < __leafs.size(); i++)
      delete __leafs[i];
}

/* Set cost matrix --------------------------------------------------------- */
void NaiveBayes::setCostMatrix(const Matrix2D<double> &cost)
{
    if (XSIZE(cost)!=K || YSIZE(cost)!=K)
        REPORT_ERROR(1,"Cost matrix does not have the apropriate size");
    __cost=cost;
}

/* Do inference ------------------------------------------------------------ */
int NaiveBayes::doInference(const Matrix1D<double> &newFeatures,
    double &cost)
{
    Matrix1D<double> classesProbs;
    classesProbs = __priorProbsLog10;
    for(int f=0; f<Nfeatures; f++)
        for (int k=0; k<K; k++)
        {
      	    double p = __leafs[f]->assignProbability(newFeatures(f), k);
            
            if (ABS(p) < 1e-2) classesProbs(k) += -2*__weights(f);
            else               classesProbs(k) += __weights(f)*std::log10(p);

            #ifdef DEBUG_FINE_CLASSIFICATION
                if(debugging == true)
                {
                    std::cout << "Feature " << f
                              << " Probability for class " << k << " = "
                              << classesProbs(k) << " increase= " << p 
                              << std::endl;
                    char c;
// COSS                    std::cin >> c;
//                    if (c=='q') debugging = false;
                }
            #endif        	        
        }
    classesProbs-=classesProbs.computeMax();
//    std::cout << "classesProbs " << classesProbs.transpose() << std::endl;
    
    for (int k=0; k<K; k++)
        classesProbs(k)=pow(10.0,classesProbs(k));
    classesProbs/=classesProbs.sum();
//    std::cout << "classesProbs norm " << classesProbs.transpose() << std::endl;
    
    Matrix1D<double> allCosts;
    allCosts=__cost*classesProbs;
    for (int k=0; k<K; k++)
        allCosts(k)=log10(allCosts(k));
//    std::cout << "allCosts " << allCosts.transpose() << std::endl;

    int bestk=0;
    cost=allCosts(0);
    for (int k=1; k<K; k++)
        if (allCosts(k)<cost)
        {
            cost=allCosts(k);
            bestk=k;
        }

    #ifdef DEBUG_CLASSIFICATION
        if(debugging == true)
        {
            for (int k=0; k<K; k++)
                classesProbs(k)=log10(classesProbs(k));
            std::cout << "Class probababilities=" << classesProbs.transpose()
                      << "\n  costs=" << allCosts.transpose()
                      << "  best class=" << bestk << std::endl;
            char c;
            // COSS std::cin >> c;
            // if (c=='q') debugging = false;
        }
    #endif
    return bestk;
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &_out, const NaiveBayes &naive)
{
    for (int f=0; f<naive.Nfeatures; f++)
        _out << "Node " << f << std::endl
             << *(naive.__leafs[f]) << std::endl;
    return _out;
}

/* Ensemble constructor ---------------------------------------------------- */
EnsembleNaiveBayes::EnsembleNaiveBayes(
    const std::vector < Matrix2D<double> >  &features,
    const Matrix1D<double> &priorProbs,
    int discreteLevels, int numberOfClassifiers, 
    double samplingFeatures, double samplingIndividuals,
    const std::string &newJudgeCombination)
{
    int NFeatures=XSIZE(features[0]);
    int NsubFeatures=CEIL(NFeatures*samplingFeatures);
    K=features.size();
    judgeCombination=newJudgeCombination;
    for (int n=0; n<numberOfClassifiers; n++)
    {
        // Produce the set of features for this subclassifier
        Matrix1D<int> subFeatures(NsubFeatures);
        FOR_ALL_ELEMENTS_IN_MATRIX1D(subFeatures)
            subFeatures(i)=ROUND(rnd_unif(0,NFeatures-1));

        // Container for the new training sample
        std::vector< Matrix2D<double> >  newFeatures;

        // Produce the data set for each class
        for (int k=0; k<K; k++)
        {
            int NIndividuals=YSIZE(features[k]);
            int NsubIndividuals=CEIL(NIndividuals*samplingIndividuals);
            Matrix1D<int> subIndividuals(NsubIndividuals);
            FOR_ALL_ELEMENTS_IN_MATRIX1D(subIndividuals)
                subIndividuals(i)=ROUND(rnd_unif(0,NsubIndividuals-1));

            Matrix2D<double> newFeaturesK;
            newFeaturesK.initZeros(NsubIndividuals,NsubFeatures);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(newFeaturesK)
                newFeaturesK(i,j)=features[k](subIndividuals(i),subFeatures(j));
            
            newFeatures.push_back(newFeaturesK);
        }
        
        // Create a Naive Bayes classifier with this data
        NaiveBayes *nb=new NaiveBayes(newFeatures, priorProbs, discreteLevels);
        ensemble.push_back(nb);
        ensembleFeatures.push_back(subFeatures);
    }
}

/* Destructor -------------------------------------------------------------- */
EnsembleNaiveBayes::~EnsembleNaiveBayes()
{
    int nmax=ensemble.size();
    for (int n=0; n<nmax; n++)
        delete ensemble[n];
}

/* Set cost matrix --------------------------------------------------------- */
void EnsembleNaiveBayes::setCostMatrix(const Matrix2D<double> &cost)
{
    int nmax=ensemble.size();
    for (int n=0; n<nmax; n++)
        ensemble[n]->setCostMatrix(cost);
}

/* Do inference ------------------------------------------------------------ */
int EnsembleNaiveBayes::doInference(const Matrix1D<double> &newFeatures,
    double &cost, Matrix1D<int> &votes)
{

    int nmax=ensemble.size();
    Matrix1D<double> minCost, maxCost;
    votes.initZeros(K);
    minCost.initZeros(K); minCost.initConstant(1);
    maxCost.initZeros(K); maxCost.initConstant(1);
    double bestMinCost=0;
    int bestClass;
    Matrix1D<double> newFeaturesn;
    for (int n=0; n<nmax; n++)
    {
        double costn;
        newFeaturesn.initZeros(XSIZE(ensembleFeatures[n]));
        FOR_ALL_ELEMENTS_IN_MATRIX1D(newFeaturesn)
            newFeaturesn(i)=newFeatures(ensembleFeatures[n](i));
        int k=ensemble[n]->doInference(newFeaturesn, costn);
        votes(k)++;
        if (minCost(k)>0 || minCost(k)>costn) minCost(k)=costn;
        if (maxCost(k)>0 || maxCost(k)<costn) maxCost(k)=costn;
        if (minCost(k)<bestMinCost)
        {
            bestMinCost=minCost(k);
            bestClass=k;
        }
    }
    if      (judgeCombination[bestClass]=='m') cost=minCost(bestClass);
    else if (judgeCombination[bestClass]=='M') cost=maxCost(bestClass);
    else    cost=minCost(bestClass);
    return bestClass;
}
