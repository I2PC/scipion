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

#ifndef XMIPP__SVM_CLASSIFIER_HH__
#define XMIPP__SVM_CLASSIFIER_HH__

/* Includes ---------------------------------------------------------------- */
#include <data/xmipp_program.h>
#include "svm.h"

/**@defgroup SVMClassifier SVM Classifier
   @ingroup ClassificationLibrary */
//@{

/** SVM classifier class.
 *  This class use the SVMLIB Library in order to classify
 *  the data using the svm method
 */
class SVMClassifier
{
public:
    svm_parameter param;
    svm_problem prob;
    svm_model *model;
public:

    SVMClassifier();
    ~SVMClassifier();
    void SVMTrain(MultidimArray<double> &trainSet,MultidimArray<double> &lable);
    double  predict(MultidimArray<double> &featVec,double &score);
    void SaveModel(const FileName &fnModel);
    void LoadModel(const FileName &fnModel);
};
//@}
#endif
