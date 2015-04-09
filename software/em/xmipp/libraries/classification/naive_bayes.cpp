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
 *  e-mail address 'xmipp@cnb.csic.es'                                  
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
int splitHistogramsUsingEntropy(const std::vector<Histogram1D> &hist,
                                size_t l0, size_t lF)
{
    // Number of classes
    int K = hist.size();

    // Set everything outside l0 and lF to zero, and make it a PDF
    std::vector<Histogram1D> histNorm;
    for (int k = 0; k < K; k++)
    {
        Histogram1D histaux = hist[k];
        for (size_t l = 0; l < XSIZE(histaux); l++)
            if (l < l0 || l > lF)
                DIRECT_A1D_ELEM(histaux,l) = 0;
        histaux *= 1.0/histaux.sum();
        histNorm.push_back(histaux);
    }

    // Compute for each class the probability of being l<=l0 and l>l0
    MultidimArray<double> p(K, 2);
    for (int k = 0; k < K; k++)
    {
        const Histogram1D& histogram=histNorm[k];
        DIRECT_A2D_ELEM(p,k, 0) = DIRECT_A1D_ELEM(histogram,l0);
        DIRECT_A2D_ELEM(p,k, 1) = 0;
        for (size_t l = l0 + 1; l <= lF; l++)
            DIRECT_A2D_ELEM(p,k, 1) += DIRECT_A1D_ELEM(histogram,l);
    }

    // Compute the splitting l giving maximum entropy
    double maxEntropy = 0;
    int lmaxEntropy = -1;
    size_t l = l0;
    while (l < lF)
    {
        // Compute the entropy of the classes if we split by l
        double entropy = 0;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(p)
        {
            double aux=DIRECT_MULTIDIM_ELEM(p,n);
            if (aux != 0)
                entropy -= aux * log10(aux);
        }

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
        ++l;

        // Update probabilities of being l<=l0 and l>l0
        for (int k = 0; k < K; k++)
        {
            const Histogram1D& histogram=histNorm[k];
            double aux=DIRECT_A1D_ELEM(histogram,l);
            DIRECT_A2D_ELEM(p,k, 0) += aux;
            DIRECT_A2D_ELEM(p,k, 1) -= aux;
        }
    }

#ifdef DEBUG_SPLITTING_USING_ENTROPY
    std::cout << "Finally in l=[" << l0 << "," << lF
    << " Max Entropy:" << maxEntropy
    << " lmax=" << lmaxEntropy << std::endl;
#endif

    // If the point giving the maximum entropy is too much on the extreme,
    // substitute it by the middle point
    if (lmaxEntropy<=2 || lmaxEntropy>=(int)lF-2)
        lmaxEntropy = (int)ceil((lF + l0)/2.0);

    return lmaxEntropy;
}

/* Constructor ------------------------------------------------------------- */
LeafNode::LeafNode(const std::vector < MultidimArray<double> > &leafFeatures,
                   int discrete_levels)
{
    __discreteLevels = discrete_levels;
    K = leafFeatures.size();
    if (__discreteLevels==0)
    {
        // This is a dummy node for features that cannot classify
        MultidimArray<int> newBins(1);
        A1D_ELEM(newBins,0)=0;
        Histogram1D hist;
        hist.resize(1);
        A1D_ELEM(hist,0)=1;
        IrregularHistogram1D irregHist;
        for (int k=0; k<K; k++)
        {
            irregHist.init(hist, newBins);
            irregHist.selfNormalize();
            __leafPDF.push_back(irregHist);
        }
    }
    else
    {
        // Compute the minimum and maximum of each class
        double minval=0., maxval=0.;
        for(int k=0; k<K; k++)
        {
            double minvalk=0., maxvalk=0.;
            leafFeatures[k].computeDoubleMinMax(minvalk, maxvalk);
            if (k==0)
            {
                minval=minvalk;
                maxval=maxvalk;
            }
            else
            {
                minval=std::min(minval,minvalk);
                maxval=std::max(maxval,maxvalk);
            }
        }
        if (minval==maxval)
        {
            __discreteLevels=0;
            return;
        }

        // Compute the PDF of each class
        std::vector<Histogram1D> hist(K);
        for (int k=0; k<K; k++)
        {
            // There is variation of this feature for this class
            compute_hist(leafFeatures[k], hist[k], minval, maxval, 100);
            hist[k] += 1; // Apply Laplace correction
        }

        // Split the histograms into discrete_level (power of 2) bins
        std::queue< Matrix1D<int> > intervals, splittedIntervals;
        Matrix1D<int> limits(2);
        VECTOR_R2(limits,0,99);
        intervals.push(limits);
        int imax=ROUND(log2(__discreteLevels));
        for (int i=0; i<imax; i++)
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
        MultidimArray<int> newBins(__discreteLevels);
        imax=intervals.size();
        for (int i=0; i<imax; i++)
        {
            A1D_ELEM(newBins,i) = intervals.front()(1);
            intervals.pop();
        }

        // Compute now the irregular histograms
        IrregularHistogram1D irregHist;
        for (int k=0; k<K; k++)
        {
            irregHist.init(hist[k], newBins);
            irregHist.selfNormalize();
            __leafPDF.push_back(irregHist);
        }
    }
}

/* Get the number of levels ------------------------------------------------ */
int LeafNode::getNumberOfLevels()
{
    return __discreteLevels;
}

/* Assign probability ------------------------------------------------------ */
double LeafNode::assignProbability(double value, int k) const
{
    const IrregularHistogram1D& hist=__leafPDF[k];
    int index = hist.val2Index(value);
    return DIRECT_A1D_ELEM(hist.__hist,index);
}

/* Compute weight ---------------------------------------------------------- */
double LeafNode::computeWeight() const
{
    double retval=0;
    for (int k1=0; k1<K; k1++)
        for (int k2=0; k2<K; k2++)
        {
            if (k1==k2)
                continue;
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
    const std::vector< MultidimArray<double> > &features,
    const Matrix1D<double> &priorProbs,
    int discreteLevels)
{
    K = features.size();
    Nfeatures=XSIZE(features[0]);
    __priorProbsLog10.initZeros(K);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(__priorProbsLog10)
    VEC_ELEM(__priorProbsLog10,i)=log10(VEC_ELEM(priorProbs,i));

    // Create a dummy leaf for features that cannot classify
    std::vector < MultidimArray<double> > aux(K);
    dummyLeaf=new LeafNode(aux,0);

    // Build a leafnode for each feature and assign a weight
    __weights.initZeros(Nfeatures);
    for (int f=0; f<Nfeatures; f++)
    {
        for (int k=0; k<K; k++)
            features[k].getCol(f, aux[k]);
        LeafNode *leaf=new LeafNode(aux,discreteLevels);
        if (leaf->__discreteLevels>0)
        {
            __leafs.push_back(leaf);
            DIRECT_A1D_ELEM(__weights,f)=__leafs[f]->computeWeight();
        }
        else
        {
            __leafs.push_back(dummyLeaf);
            DIRECT_A1D_ELEM(__weights,f)=0;
            delete leaf;
        }
#ifdef DEBUG_WEIGHTS

        if(debugging == true)
        {
            std::cout << "Node " << f << std::endl;
            std::cout << *(__leafs[f]) << std::endl;
            //char c;
            //std::cin >> c;
        }
#endif

    }
    double norm=__weights.computeMax();
    if (norm>0)
    	__weights *= 1.0/norm;

    // Set default cost matrix
    __cost.resizeNoCopy(K,K);
    __cost.initConstant(1);
    for (int i=0; i<K; i++)
        MAT_ELEM(__cost,i,i)=0;
}

/* Destructor -------------------------------------------------------------- */
NaiveBayes::~NaiveBayes()
{
    int imax=__leafs.size();
    for (int i = 0; i < imax; i++)
        if (__leafs[i]!=dummyLeaf)
            delete __leafs[i];
    delete dummyLeaf;
}

/* Set cost matrix --------------------------------------------------------- */
void NaiveBayes::setCostMatrix(const Matrix2D<double> &cost)
{
    size_t iK=(size_t) K;
    if (MAT_XSIZE(cost)!=iK || MAT_YSIZE(cost)!=iK)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,"Cost matrix does not have the apropriate size");
    __cost=cost;
}

/* Do inference ------------------------------------------------------------ */
int NaiveBayes::doInference(const MultidimArray<double> &newFeatures, double &cost,
                            Matrix1D<double> &classesProbs, Matrix1D<double> &allCosts)
{
    classesProbs=__priorProbsLog10;
    for(int f=0; f<Nfeatures; f++)
    {
        const LeafNode &leaf_f=*(__leafs[f]);
        double newFeatures_f=DIRECT_A1D_ELEM(newFeatures,f);
        for (int k=0; k<K; k++)
        {
            double p = leaf_f.assignProbability(newFeatures_f, k);

            if (fabs(p) < 1e-2)
                VEC_ELEM(classesProbs,k) += -2*DIRECT_A1D_ELEM(__weights,f);
            else
                VEC_ELEM(classesProbs,k) += DIRECT_A1D_ELEM(__weights,f)*std::log10(p);

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
    }

    classesProbs-=classesProbs.computeMax();
    //    std::cout << "classesProbs " << classesProbs.transpose() << std::endl;

    for (int k=0; k<K; k++)
        VEC_ELEM(classesProbs,k)=pow(10.0,VEC_ELEM(classesProbs,k));
    classesProbs*=1.0/classesProbs.sum();
    //    std::cout << "classesProbs norm " << classesProbs.transpose() << std::endl;

    allCosts=__cost*classesProbs;
    //    std::cout << "allCosts " << allCosts.transpose() << std::endl;

    int bestk=0;
    cost=VEC_ELEM(allCosts,0)=std::log10(VEC_ELEM(allCosts,0));
    for (int k=1; k<K; k++)
    {
        VEC_ELEM(allCosts,k)=std::log10(VEC_ELEM(allCosts,k));
        if (VEC_ELEM(allCosts,k)<cost)
        {
            cost=VEC_ELEM(allCosts,k);
            bestk=k;
        }
    }

#ifdef DEBUG_CLASSIFICATION
    if(debugging == true)
    {
        for (int k=0; k<K; k++)
            classesProbs(k)=log10(classesProbs(k));
        std::cout << "Class probababilities=" << classesProbs.transpose()
        << "\n  costs=" << allCosts.transpose()
        << "  best class=" << bestk << " cost=" << cost << std::endl;
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
    {
        _out << "Node " << f << std::endl;
        _out << *(naive.__leafs[f]) << std::endl;
    }
    return _out;
}

/* Ensemble constructor ---------------------------------------------------- */
#define WEIGHTED_SAMPLING
EnsembleNaiveBayes::EnsembleNaiveBayes(
    const std::vector < MultidimArray<double> >  &features,
    const Matrix1D<double> &priorProbs,
    int discreteLevels, int numberOfClassifiers,
    double samplingFeatures, double samplingIndividuals,
    const std::string &newJudgeCombination)
{
    int NFeatures=XSIZE(features[0]);
    int NsubFeatures=CEIL(NFeatures*samplingFeatures);
    K=features.size();
    judgeCombination=newJudgeCombination;

#ifdef WEIGHTED_SAMPLING
    // Measure the classification power of each variable
    NaiveBayes *nb_weights=new NaiveBayes(features, priorProbs, discreteLevels);
    MultidimArray<double> weights=nb_weights->__weights;
    delete nb_weights;
    double sumWeights=weights.sum();
#endif

    for (int n=0; n<numberOfClassifiers; n++)
    {
        // Produce the set of features for this subclassifier
        MultidimArray<int> subFeatures(NsubFeatures);
        FOR_ALL_ELEMENTS_IN_ARRAY1D(subFeatures)
        {
#ifdef WEIGHTED_SAMPLING
            double random_sum_weight=rnd_unif(0,sumWeights);
            int j=0;
            do
            {
                double wj=DIRECT_A1D_ELEM(weights,j);
                if (wj<random_sum_weight)
                {
                    random_sum_weight-=wj;
                    j++;
                    if (j==NFeatures)
                    {
                        j=NFeatures-1;
                        break;
                    }
                }
                else
                    break;
            }
            while (true);
            DIRECT_A1D_ELEM(subFeatures,i)=j;
#else

            DIRECT_A1D_ELEM(subFeatures,i)=round(rnd_unif(0,NFeatures-1));
#endif

        }

        // Container for the new training sample
        std::vector< MultidimArray<double> >  newFeatures;

        // Produce the data set for each class
        for (int k=0; k<K; k++)
        {
            int NIndividuals=YSIZE(features[k]);
            int NsubIndividuals=CEIL(NIndividuals*samplingIndividuals);
            MultidimArray<int> subIndividuals(NsubIndividuals);
            FOR_ALL_ELEMENTS_IN_ARRAY1D(subIndividuals)
            subIndividuals(i)=ROUND(rnd_unif(0,NsubIndividuals-1));

            MultidimArray<double> newFeaturesK;
            newFeaturesK.initZeros(NsubIndividuals,NsubFeatures);
            const MultidimArray<double>& features_k=features[k];
            FOR_ALL_ELEMENTS_IN_ARRAY2D(newFeaturesK)
            DIRECT_A2D_ELEM(newFeaturesK,i,j)=DIRECT_A2D_ELEM(features_k,
                                              DIRECT_A1D_ELEM(subIndividuals,i),
                                              DIRECT_A1D_ELEM(subFeatures,j));

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
                                    double &cost, MultidimArray<int> &votes,
                                    Matrix1D<double> &classesProbs, Matrix1D<double> &allCosts)
{
    int nmax=ensemble.size();
    MultidimArray<double> minCost, maxCost;
    votes.initZeros(K);
    minCost.initZeros(K);
    minCost.initConstant(1);
    maxCost.initZeros(K);
    maxCost.initConstant(1);
    double bestMinCost=0;
    int bestClass=0;
    MultidimArray<double> newFeaturesn;
    for (int n=0; n<nmax; n++)
    {
        double costn;
        newFeaturesn.initZeros(XSIZE(ensembleFeatures[n]));
        FOR_ALL_ELEMENTS_IN_ARRAY1D(newFeaturesn)
        newFeaturesn(i)=newFeatures(ensembleFeatures[n](i));
        int k=ensemble[n]->doInference(newFeaturesn, costn, classesProbs, allCosts);
        votes(k)++;
        if (minCost(k)>0 || minCost(k)>costn)
            minCost(k)=costn;
        if (maxCost(k)>0 || maxCost(k)<costn)
            maxCost(k)=costn;
        if (minCost(k)<bestMinCost)
        {
            bestMinCost=minCost(k);
            bestClass=k;
        }
    }
    if      (judgeCombination[bestClass]=='m')
        cost=minCost(bestClass);
    else if (judgeCombination[bestClass]=='M')
        cost=maxCost(bestClass);
    else
        cost=minCost(bestClass);
    return bestClass;
}

/* Do inference for class ------------------------------------------------- */
int EnsembleNaiveBayes::doInferenceForClass(int classNumber, const Matrix1D<double> &newFeatures, double &cost,
        Matrix1D<double> &classesProbs, Matrix1D<double> &allCosts)
{
    int nmax=ensemble.size();
    double minCost=1, maxCost=1;
    int votes=0;
    MultidimArray<double> newFeaturesn;
    for (int n=0; n<nmax; n++)
    {
        double costn;
        const MultidimArray<int> &ensembleFeatures_n=ensembleFeatures[n];
        newFeaturesn.resizeNoCopy(XSIZE(ensembleFeatures_n));
        FOR_ALL_ELEMENTS_IN_ARRAY1D(newFeaturesn)
        {
            int idx=A1D_ELEM(ensembleFeatures_n,i);
            A1D_ELEM(newFeaturesn,i)=VEC_ELEM(newFeatures,idx);
        }

        int k=ensemble[n]->doInference(newFeaturesn, costn, classesProbs, allCosts);
        if (k==classNumber)
        {
            votes++;
            if (minCost>0 || minCost>costn)
                minCost=costn;
            if (maxCost>0 || maxCost<costn)
                maxCost=costn;
        }
    }
    if      (judgeCombination[classNumber]=='m')
        cost=minCost;
    else if (judgeCombination[classNumber]=='M')
        cost=maxCost;
    else
        cost=minCost;
    return votes;
}
