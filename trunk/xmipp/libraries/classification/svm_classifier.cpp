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
#include "svm_classifier.h"

SVMClassifier::SVMClassifier()
{
    param.svm_type = C_SVC;
    param.kernel_type = LINEAR1;
    param.degree = 3;
    param.gamma = 0;
    param.coef0 = 0;
    param.nu = 0.5;
    param.cache_size = 100;
    param.C = 1;
    param.eps = 0.001;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
}
SVMClassifier::~SVMClassifier()
{
    svm_free_and_destroy_model(&model);
    svm_destroy_param(&param);
    free(prob.y);
    free(prob.x);
}
void SVMClassifier::SVMTrain(MultidimArray<double> &trainSet,MultidimArray<double> &lable)
{
    prob.l = YSIZE(trainSet);
    prob.y = new double[prob.l];
    prob.x = new svm_node *[prob.l+1];
    const char *error_msg;
    for (int i=0;i<YSIZE(trainSet);i++)
    {
        prob.x[i]=new svm_node[XSIZE(trainSet)+1];
        int cnt = 0;
        for (int j=0;j<XSIZE(trainSet);j++)
        {
            if (trainSet(i,j)==0)
                continue;
            else
            {
                prob.x[i][cnt].value=DIRECT_A2D_ELEM(trainSet,i,j);
                prob.x[i][cnt].index=j+1;
                cnt++;
            }
        }
        prob.x[i][cnt].index=-1;
        prob.x[i][cnt].value=2;
        prob.y[i] = DIRECT_A1D_ELEM(lable,i);
    }
    std::ofstream fh_auto;
    fh_auto.open("trainvectors.txt");
    for (int i=0;i<YSIZE(trainSet);i++)
    {
        for (int j=0;j<=XSIZE(trainSet);j++)
        {
            fh_auto <<"("<<prob.x[i][j].value<<","<<prob.x[i][j].index<<")"<<" ";
        }

        fh_auto<<" "<<prob.y[i]<<std::endl;
    }
    error_msg = svm_check_parameter(&prob,&param);
    if(error_msg)
    {
        fprintf(stderr,"ERROR: %s\n",error_msg);
        exit(1);
    }
    model=svm_train(&prob,&param);
}
int SVMClassifier::predict(MultidimArray<double> &featVec)
{
    svm_node *x_space;
    int cnt=0;
    int nr_class=svm_get_nr_class(model);
    double *prob_estimates=(double *) malloc(nr_class*sizeof(double));
    x_space=new svm_node[XSIZE(featVec)+1];

    for (int i=0;i<XSIZE(featVec);i++)
    {
        if (DIRECT_A1D_ELEM(featVec,i)==0)
            continue;
        else
        {
            x_space[cnt].value=DIRECT_A1D_ELEM(featVec,i);
            x_space[cnt].index=i+1;
            cnt++;
        }
    }
    x_space[cnt].index=-1;
    return svm_predict_probability(model,x_space,prob_estimates);
}
void SVMClassifier::SaveModel(const FileName &fnModel)
{
    std::cerr<<svm_get_nr_class(model);
    svm_save_model(fnModel.c_str(),model);
}
void SVMClassifier::LoadModel(const FileName &fnModel)
{
    model=svm_load_model(fnModel.c_str());
}
