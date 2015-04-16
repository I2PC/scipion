/***************************************************************************
 *
 * Authors:     Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
 *
 * Universidad San Pablo C.E.U.
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

#include "angular_commonline.h"
#include <data/filters.h>
#include <data/geometry.h>
#include <data/xmipp_image.h>
#include <reconstruction/radon.h>

/* Constructor ------------------------------------------------------------- */
EulerSolver::EulerSolver(int dim, int pop,
                         const Matrix1D<int> &newAlreadyOptimized,
                         const Matrix1D<double> &newCurrentSolution,
                         const MultidimArray<int> &newImgIdx,
                         const Prog_Angular_CommonLine *newParent): DESolver(dim, pop)
{
    parent = newParent;
    Ndim = dim;
    NToSolve=dim/3;
    Nimg=VEC_XSIZE(newAlreadyOptimized);
    imgAvgCorrelation.initZeros(Nimg);
    imgMinCorrelation.initZeros(Nimg);
    correlationMatrix.initZeros(Nimg,Nimg);
    commonline.initZeros(3);
    alreadyOptimized=&newAlreadyOptimized;
    currentSolution=&newCurrentSolution;
    imgIdx=&newImgIdx;
    show=false;
}

void EulerSolver::setShow(bool newShow)
{
    show=newShow;
}

/* Energy function for the solver ------------------------------------------ */
double EulerSolver::EnergyFunction(double trial[], bool &bAtSolution)
{
    // Check limits
    int i,j;

    for (i=0, j=0; i<Ndim; i++,j=(j+1)%3)
    {
        if (j==1)
            trial[i]=realWRAP(trial[i],0,180);
        else
            trial[i]=realWRAP(trial[i],0,360);
    }

    // Get number of symmetries
    int Nsym=parent->L.size();

    // Evaluate goal function for this trial
    double retval=0;
    double Ncomparisons=0;
    double worseSimilarity=2;
    int idx;
    static Matrix1D<int> imgCorrelationN;
    static double minval=0;
    imgAvgCorrelation.initZeros(Nimg);
    imgMinCorrelation.initZeros(Nimg);
    imgMinCorrelation.initConstant(2);
    correlationMatrix.resize(Nimg,Nimg);
    imgCorrelationN.initZeros(Nimg);
    for (int imgi=0; imgi<Nimg; imgi++)
    {
        // Get the right angles for this image
        if      ((*alreadyOptimized)(imgi)==0)
            continue;
        else if ((*alreadyOptimized)(imgi)==1)
        {
            for (idx=0; idx<NToSolve; idx++)
                if ((*imgIdx)(idx)==imgi)
                    break;
            idx*=3;
            /*
            trial[idx]   = parent->initialSolution(3*imgi);
            trial[idx+1] = parent->initialSolution(3*imgi+1);
            trial[idx+2] = parent->initialSolution(3*imgi+2);
            */
            roti  = trial[idx++];
            tilti = trial[idx++];
            psii  = trial[idx];
        }
        else
        {
            idx=3*imgi;
            roti  = (*currentSolution)(idx++);
            tilti = (*currentSolution)(idx++);
            psii  = (*currentSolution)(idx);
        }

        // Loop for the second image
        for (int imgj=imgi+1; imgj<Nimg; imgj++)
        {
            // Get the right angles for this image
            if      ((*alreadyOptimized)(imgj)==0)
                continue;
            else if ((*alreadyOptimized)(imgj)==1)
            {
                for (idx=0; idx<NToSolve; idx++)
                    if ((*imgIdx)(idx)==imgj)
                        break;
                idx*=3;
                /*
                trial[idx]   = parent->initialSolution(3*imgj);
                trial[idx+1] = parent->initialSolution(3*imgj+1);
                trial[idx+2] = parent->initialSolution(3*imgj+2);
                */
                rotj  = trial[idx++];
                tiltj = trial[idx++];
                psij  = trial[idx];
            }
            else
            {
                idx=3*imgj;
                rotj  = (*currentSolution)(idx++);
                tiltj = (*currentSolution)(idx++);
                psij  = (*currentSolution)(idx);
            }

            // Check that at least one of the two images is new
            if ((*alreadyOptimized)(imgi)==2 && (*alreadyOptimized)(imgj)==2)
                continue;

            double similarity=similarityBetweenTwoLines(imgi,imgj);
            if (similarity>0)
            {
                retval -= similarity;
                Ncomparisons++;
                worseSimilarity=XMIPP_MIN(worseSimilarity,similarity);
                imgAvgCorrelation(imgi)+=similarity;
                imgAvgCorrelation(imgj)+=similarity;
                imgMinCorrelation(imgi)=XMIPP_MIN(similarity,
                                                  imgMinCorrelation(imgi));
                imgMinCorrelation(imgj)=XMIPP_MIN(similarity,
                                                  imgMinCorrelation(imgj));
                correlationMatrix(imgi,imgj)=correlationMatrix(imgj,imgi)=
                                                 similarity;
                imgCorrelationN(imgi)++;
                imgCorrelationN(imgj)++;
            }

            double original_rotj=rotj;
            double original_tiltj=tiltj;
            double original_psij=psij;

            for (int sym=0; sym<Nsym; sym++)
            {
                Euler_apply_transf(parent->L[sym],parent->R[sym],
                                   original_rotj,original_tiltj,original_psij,
                                   rotj, tiltj, psij);
                similarity=similarityBetweenTwoLines(imgi,imgj);
                if (similarity>0)
                {
                    retval -= similarity;
                    Ncomparisons++;
                    worseSimilarity=XMIPP_MIN(worseSimilarity,similarity);
                    imgAvgCorrelation(imgi)+=similarity;
                    imgAvgCorrelation(imgj)+=similarity;
                    imgMinCorrelation(imgi)=XMIPP_MIN(similarity,
                                                      imgMinCorrelation(imgi));
                    imgMinCorrelation(imgj)=XMIPP_MIN(similarity,
                                                      imgMinCorrelation(imgj));
                    correlationMatrix(imgi,imgj)=correlationMatrix(imgj,imgi)=
                                                     similarity;
                    imgCorrelationN(imgi)++;
                    imgCorrelationN(imgj)++;
                }
            }
        }
    }
    if (Ncomparisons>0)
        retval/=Ncomparisons;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(imgAvgCorrelation)
    if (imgCorrelationN(i)!=0)
        imgAvgCorrelation(i)/=imgCorrelationN(i);
    if (show)
        std::cout
        << "Average Distance  = " << retval << std::endl
        << "Worse   Distance  = " << -worseSimilarity << std::endl
        << "imgAvgCorrelation = " << imgAvgCorrelation.transpose() << std::endl
        << "imgMinCorrelation = " << imgMinCorrelation.transpose() << std::endl
        << std::endl
        << std::endl
        ;
    if (retval<minval)
    {
        minval=retval;
        std::cout << "MinEnergy=" << minval << std::endl;
    }
    return retval;
    //    return 0.5*(retval-worseSimilarity);
}

/* Distance between two common lines --------------------------------------- */
//#define DEBUG
double EulerSolver::similarityBetweenTwoLines(int imgi, int imgj)
{
    // Compute the direction of the common line in the
    // universal coordinate system
    Euler_direction(roti, tilti, psii, normali);
    Euler_direction(rotj, tiltj, psij, normalj);
    vectorProduct(normali,normalj,commonline);
    if (commonline.module()<XMIPP_EQUAL_ACCURACY)
    {
        // They do not share a common line but a common plane
        if (show)
        {
            std::cout
            << "imgi=" << imgi << " (rot,tilt,psi)=("
            << roti << "," << tilti << "," << psii << ") normali="
            << normali.transpose() << std::endl
            << "imgj=" << imgj << " (rot,tilt,psi)=("
            << rotj << "," << tiltj << "," << psij << ") normalj="
            << normalj.transpose() << std::endl
            ;
            std::cout << "These two images are taken from the same direction\n";
            std::cout << "Press any key\n";
            char c;
            std::cin >> c;
        }
        return -1;
    }
    commonline.selfNormalize();

    // Compute the direction of the common line in each
    // of the projection coordinate systems
    Uproject_to_plane(commonline, roti, tilti, psii, commonlinei);
    Uproject_to_plane(commonline, rotj, tiltj, psij, commonlinej);

    // Compute the angle of the common line in i and j images
    double angi=RAD2DEG(atan2(YY(commonlinei),XX(commonlinei)));
    double angj=RAD2DEG(atan2(YY(commonlinej),XX(commonlinej)));

    int idxAngi = (int)intWRAP(-((int)angi),0,359);
    int idxAngj = (int)intWRAP(-((int)angj),0,359);

    double retval1=
        correlationIndex(parent->radon[imgi][idxAngi],
                          parent->radon[imgj][idxAngj]);
    int idxBesti=parent->bestLine(imgi,imgj);
    int idxBestj=parent->bestLine(imgj,imgi);
    int retvali=ABS(idxBesti-idxAngi);
    if (retvali>180)
        retvali=ABS(retvali-360);
    int retvalj=ABS(idxBestj-idxAngj);
    if (retvalj>180)
        retvalj=ABS(retvalj-360);
    retvali=XMIPP_MIN(retvali,10);
    retvalj=XMIPP_MIN(retvalj,10);

    const double i180=0.5*1.0/180.0;
    double retval2=1-(retvali+retvalj)*i180;

    if (show)
    {
        std::cout
        << "imgi=" << imgi << " (rot,tilt,psi)=("
        << roti << "," << tilti << "," << psii << ") normali="
        << normali.transpose() << std::endl
        << "imgj=" << imgj << " (rot,tilt,psi)=("
        << rotj << "," << tiltj << "," << psij << ") normalj="
        << normalj.transpose() << std::endl
        << "commonline= " << commonline.transpose() << std::endl
        << "in imgi=" << commonlinei.transpose() << " anglei=" << angi
        << " (" << idxAngi << ")\n"
        << "in imgj=" << commonlinej.transpose() << " anglej=" << angj
        << " (" << idxAngj << ")\n"
        << "Distance between lines = " << retval1 << " " << retval2
        << std::endl
        ;
        parent->radon[imgi][idxAngi].write("PPPradoni1.txt");
        parent->radon[imgj][idxAngj].write("PPPradonj1.txt");
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
    }

#ifdef DEBUG
    for (idxAngi=0; idxAngi<360; idxAngi++)
    {
        double retval3=
            correlationIndex(parent->radon[imgi][idxAngi],
                              parent->radon[imgj][idxAngj]);
        ;

        if (show)
        {
            std::cout
            << " (" << idxAngi << ")\n"
            << "Distance between lines = " << retval3 << std::endl
            << std::endl
            ;
            parent->radon[imgi][idxAngi].write("PPPradoni2.txt");
            parent->radon[imgj][idxAngj].write("PPPradonj2.txt");
            parent->radonDerivative[imgi][idxAngi].write("PPPradonDerivativei2.txt");
            parent->radonDerivative[imgj][idxAngj].write("PPPradonDerivativej2.txt");
            std::cout << "Press any key\n";
            char c;
            std::cin >> c;
            if (c=='q')
                break;
        }
    }
#endif
    return retval1*retval2;
}
#undef DEBUG

/* Wrapper for Powell ------------------------------------------------------ */
EulerSolver *global_Eulersolver=NULL;

double wrapperSolverEnergy(double trial[], void *prm)
{
    bool bAtSolution;
    return global_Eulersolver->EnergyFunction(trial+1,bAtSolution);
}

/* Parameters -------------------------------------------------------------- */
void Prog_Angular_CommonLine::read(int argc, const char **argv)
{
    fnSel = getParameter(argc,argv,"-i");
    fnOut = getParameter(argc,argv,"-oang");
    fnSym = getParameter(argc,argv,"-sym","");
    NGen   = textToInteger(getParameter(argc,argv,"-NGen","50000"));
    NGroup = textToInteger(getParameter(argc,argv,"-NGroup","10"));
    tryInitial = checkParameter(argc,argv,"-tryInitial");
}

void Prog_Angular_CommonLine::show() const
{
    std::cout << "Selfile:     " << fnSel      << std::endl
    << "Output:      " << fnOut      << std::endl
    << "Symmetry:    " << fnSym      << std::endl
    << "Generations: " << NGen       << std::endl
    << "Groups:      " << NGroup     << std::endl
    << "Try Initial: " << tryInitial << std::endl
    ;
}

void Prog_Angular_CommonLine::usage() const
{
    std::cerr << "angular_commonline\n"
    << "   -i <selfile>     : SelFile with input images\n"
    << "   -oang <docfile>  : Docfile with the angular assignment\n"
    << "  [-NGen <g=50000>] : Number of generations\n"
    << "  [-NGroup <N=10>]  : Number of group comparisons\n"
    << "  [-tryInitial]     : Do not optimize\n"
    << "  [-sym <symfile>]  : Symmetry\n"
    ;
}

/* Produce side info ------------------------------------------------------- */
void Prog_Angular_CommonLine::produceSideInfo()
{
    // Read selfile images and the initial angles
    SF.read(fnSel);
    int Nimg=SF.size();
    initialSolution.resize(3*Nimg);
    int idx=0;
    Image<double> I;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        I.readApplyGeo(SF,__iter.objId);
        img.push_back(I());
        initialSolution(3*idx)=I.rot();
        initialSolution(3*idx+1)=I.tilt();
        initialSolution(3*idx+2)=I.psi();
        idx++;
    }

    // Set current solution
    // The first image is already optimized and its angles are 0,0,0
    currentSolution.initZeros(3*Nimg);
    alreadyOptimized.initZeros(Nimg);
    alreadyOptimized(0)=2;

    // Symmetry List
    if (fnSym!="")
    {
        SL.readSymmetryFile(fnSym);
        for (int sym=0; sym<SL.symsNo(); sym++)
        {
            Matrix2D<double> auxL, auxR;
            SL.getMatrices(sym,auxL,auxR);
            auxL.resize(3,3);
            auxR.resize(3,3);
            L.push_back(auxL);
            R.push_back(auxR);
        }
    }

    // Compute the Radon transform of each image
    const double deltaRot=1.0;
    MultidimArray<double> RT;
    MultidimArray<double> projection;
    std::cout << "Preprocessing images ...\n";
    init_progress_bar(Nimg);
    for (int n=0; n<Nimg; n++)
    {
        // Initialize list of projection and derivatives for this image
        std::vector< MultidimArray<double> > dummyList;
        radon.push_back(dummyList);

        // Compute Radon transform of each image
        Radon_Transform(img[n],deltaRot,RT);

        // Separate each projection line and compute derivative
        for (size_t i=0; i<YSIZE(RT); i++)
        {
            RT.getRow(i,projection);
            radon[n].push_back(projection);
        }
        progress_bar(n);
    }
    progress_bar(Nimg);

    // Compute the best line in image i that correlates with image j
    bestLine.initZeros(Nimg,Nimg);
    bestLineCorrelation.initZeros(Nimg,Nimg);
    std::cout << "Computing common lines ...\n";
    init_progress_bar(Nimg);
    for (int n1=0; n1<Nimg; n1++)
    {
        bestLineCorrelation(n1,n1)=1;
        bestLine(n1,n1)=-1;
        for (int n2=n1+1; n2<Nimg; n2++)
        {
            int lmax=radon[n1].size();
            for (int l1=0; l1<lmax; l1++)
                for (int l2=0; l2<lmax; l2++)
                {
                    double corrl1l2=correlationIndex(radon[n1][l1],
                                                      radon[n2][l2]);
                    if (corrl1l2>bestLineCorrelation(n1,n2))
                    {
                        bestLineCorrelation(n1,n2)=
                            bestLineCorrelation(n2,n1)=corrl1l2;
                        bestLine(n1,n2)=l1;
                        bestLine(n2,n1)=l2;
                    }
                }
        }
        progress_bar(n1);
    }
    progress_bar(Nimg);
}

/* OptimizeGroup ----------------------------------------------------------- */
double Prog_Angular_CommonLine::optimizeGroup(const Matrix1D<int> &imgIdx,
        Matrix1D<double> &solution, bool show)
{
    int    count_repetitions = 0;
    double current_energy = 0;
    double previous_energy = 2;
    double first_energy = 2;
    int    NGenStep=NGen/100;
    int    NToSolve=VEC_XSIZE(imgIdx);

    // Optimize with Differential Evolution
    // Setup solver
    Matrix1D<double> minAllowed(3*NToSolve), maxAllowed(3*NToSolve);
    int idx=0;
    for (int i=0; i<NToSolve; i++)
    {
        maxAllowed(idx++)=360;
        maxAllowed(idx++)=180;
        maxAllowed(idx++)=360;
    }
    solver=new EulerSolver(3*NToSolve,30*NToSolve,
                           alreadyOptimized, currentSolution, imgIdx, this);
    solver->Setup(MATRIX1D_ARRAY(minAllowed), MATRIX1D_ARRAY(maxAllowed),
                  stBest2Bin, 0.5, 0.8);
    global_Eulersolver=solver;

    // Really optimize
    bool done = false;
    int count=0;
    while (!done)
    {
        // Go NGenStep generations ahead
        solver->Solve(NGenStep);
        current_energy = solver->Energy();
        count++;

        if (current_energy<0 && first_energy>0)
            first_energy=current_energy;

        // Check if rather slow convergence
        if (first_energy<0 &&
            (ABS((previous_energy-current_energy)/
                 (previous_energy-first_energy)))<0.01)
            count_repetitions++;
        else
            count_repetitions=0;

        // Show
        if (show)
            std::cout << "Iteration: " << count
            << "   Energy: " << current_energy
            << "  ( " << count_repetitions << ")\n";


        // Stopping criterion
        if (count_repetitions>=10 && (count>=50 || count==100))
            done=true;

        previous_energy = current_energy;
    }

    // Optimize with Powell
    Matrix1D<double> steps(3*NToSolve);
    steps.initConstant(1);

    solution.initZeros(steps);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(solution)
    VEC_ELEM(solution,i) = solver->Solution()[i];
    int iter;
    double retval;
    powellOptimizer(solution,1,3*NToSolve,wrapperSolverEnergy,NULL,
                    0.001,retval,iter,steps,show);

    delete solver;
    return retval;
}

/* Optimize ---------------------------------------------------------------- */
void Prog_Angular_CommonLine::optimize(Matrix1D<double> &solution)
{
    size_t Nimg = SF.size();
    solution.initZeros(3*Nimg);
    Matrix1D<int> assigned, tabuPenalization;
    assigned.initZeros(Nimg);
    assigned(0)=1;
    tabuPenalization.initZeros(Nimg);

    // Assign all images
    int removalCounter=0;
    while (assigned.sum()<Nimg)
    {
        // Initialize the list of Euler vectors
        // There is a std::vector< Matrix1D<double> > for each image
        // Each Matrix1D<double> is an orientation found for this image
        // in an independent experiment
        std::vector < std::vector< Matrix1D<double> > > eulerAngles;
        std::vector < std::vector< double > > correlations;
        eulerAngles.clear();
        correlations.clear();
        for (size_t i=0; i<Nimg; i++)
        {
            std::vector< Matrix1D<double> > dummy1;
            std::vector< double > dummy2;
            eulerAngles.push_back(dummy1);
            correlations.push_back(dummy2);
        }

        // For each image try different pairs and get different alignments
        std::cout << "Aligning image pairs, "
        << Nimg-assigned.sum() << " images remaining ...\n";
        init_progress_bar((int)(Nimg-assigned.sum()));
        for (size_t i=0; i<Nimg; i++)
        {
            if (assigned(i) || tabuPenalization(i)>0)
                continue;
            Matrix1D<int> backupAlreadyOptimized=alreadyOptimized;

            // Perform NGroup experiments
            for (int n=0; n<NGroup; n++)
            {
                // Prepare the vector of images to optimize
                alreadyOptimized=backupAlreadyOptimized;
                alreadyOptimized(i)=1;

                // Align these two images
                Matrix1D<int> imgIdx(1);
                imgIdx(0)=i;
                Matrix1D<double> auxSolution, anglesi(3);
                double energy=optimizeGroup(imgIdx,auxSolution,false);
                anglesi(0)=auxSolution(0);
                anglesi(1)=auxSolution(1);
                anglesi(2)=auxSolution(2);

                // Keep results
                eulerAngles[i].push_back(anglesi);
                correlations[i].push_back(-energy);
                alreadyOptimized=backupAlreadyOptimized;
            }
            progress_bar(i);
        }
        progress_bar((int)(Nimg-assigned.sum()));

        // Compute for each image the variance in the top assignment
        size_t Nsym=SL.symsNo();
        double bestDistance=-2;
        int besti=-1;
        size_t topN=NGroup;
        for (size_t i=0; i<Nimg; i++)
            if (eulerAngles[i].size()<2*topN && !assigned(i) &&
                tabuPenalization(i)==0)
                topN=(size_t)ceil(eulerAngles[i].size()/2);
        for (size_t i=0; i<Nimg; i++)
        {
            // If the particle has  already been assigned skip it
            if (assigned(i) || tabuPenalization(i)>0)
                continue;

            // Sort the images by ascending correlation
            MultidimArray<double> aux;
            aux.initZeros(eulerAngles[i].size());
            for (size_t n=0; n<eulerAngles[i].size(); n++)
                aux(n)=correlations[i][n];
            MultidimArray<int> idx;
            aux.indexSort(idx);

            // Among the top, compute the distance between
            // the different Euler angles
            double distance=0;
            /*
                        if (topN!=1)
                        {
                            for (int j1=XSIZE(idx)-topN; j1<XSIZE(idx); j1++)
                            {
                                Matrix2D<double> E1;
                                int m=idx(j1)-1;
                                Euler_angles2matrix(eulerAngles[i][m](0),
                                    eulerAngles[i][m](1),eulerAngles[i][m](2),E1);
                                for (int j2=j1+1; j2<XSIZE(idx); j2++)
                                {
                                    Matrix2D<double> E2;
                                    m=idx(j2)-1;
                                    double rot=eulerAngles[i][m](0);
                                    double tilt=eulerAngles[i][m](1);
                                    double psi=eulerAngles[i][m](2);
                                    double otherrot, othertilt, otherpsi;
                                    Euler_angles2matrix(rot,tilt,psi,E2);
                                    double maxDistance=Euler_distanceBetweenMatrices(E1,E2);
                                    for (int sym=0; sym<Nsym; sym++)
                                    {
                                        Euler_apply_transf(L[sym],R[sym],
                                            rot, tilt, psi, otherrot, othertilt, otherpsi);
                                        Euler_angles2matrix(otherrot,othertilt,otherpsi,E2);
                                        double symDistance=Euler_distanceBetweenMatrices(E1,E2);;
                                        if (symDistance>maxDistance)
                                            maxDistance=symDistance;
                                    }
                                    distance+=maxDistance;
                                }
                            }
                            distance/=topN*(topN-1)/2.0;
                        }
                        else
                            distance=1.0;
            */
            distance=aux.computeMax();
            std::cout << "Image " << i << " distance=" << distance << " "
            << " Ncomparisons= " << eulerAngles[i].size()
            << " topN=" << topN << std::endl;
            std::cout.flush();

            if (distance>bestDistance)
            {
                bestDistance=distance;
                besti=i;
            }
        }

        if (bestDistance>0)
        {
            // Sort the images by ascending correlation in the best cluster
        	MultidimArray<double> aux;
            aux.initZeros(eulerAngles[besti].size());
            for (size_t n=0; n<eulerAngles[besti].size(); n++)
                aux(n)=correlations[besti][n];
            MultidimArray<int> idx;
            aux.indexSort(idx);

            // Keep only the topN
            std::vector< Matrix1D<double> > bestEulerAngles;
            std::vector< double > bestCorrelations;
            for (size_t n=XSIZE(idx)-topN; n<XSIZE(idx); n++)
            {
                bestEulerAngles.push_back(eulerAngles[besti][idx(n)-1]);
                bestCorrelations.push_back(correlations[besti][idx(n)-1]);
            }

            // Show best resolved image
            std::cout << "Candidates for image " << besti << " distance="
            << bestDistance << std::endl;
            for (size_t n=0; n<bestEulerAngles.size(); n++)
            {
                Matrix1D<double> direction;
                Euler_direction(bestEulerAngles[n](0),
                                bestEulerAngles[n](1),bestEulerAngles[n](2),
                                direction);
                std::cout << bestEulerAngles[n].transpose() << " corr= "
                << bestCorrelations[n] << " dir= "
                << direction.transpose() << std::endl;
            }

            // Cluster the different solutions
            std::vector< int > clusterBestAssignment;
            std::vector< double > clusterBestCorrelation;
            Matrix1D<int> alreadyClustered;
            alreadyClustered.initZeros(topN);
            for (size_t j1=0; j1<topN; j1++)
            {
                if (!alreadyClustered(j1))
                {
                    alreadyClustered(j1)=1;

                    Matrix2D<double> E1;
                    Euler_angles2matrix(bestEulerAngles[j1](0),
                                        bestEulerAngles[j1](1),bestEulerAngles[j1](2),
                                        E1);
                    double bestCorrelation=bestCorrelations[j1];
                    int bestAssignment=j1;
                    for (size_t j2=j1+1; j2<topN; j2++)
                    {
                        if (!alreadyClustered(j2))
                        {
                            Matrix2D<double> E2;
                            Euler_angles2matrix(bestEulerAngles[j2](0),
                                                bestEulerAngles[j2](1),
                                                bestEulerAngles[j2](2),E2);
                            double bestCorrE1E2=
                                Euler_distanceBetweenMatrices(E1,E2);
                            for (size_t sym=0; sym<Nsym; sym++)
                            {
                                double otherrot, othertilt, otherpsi;
                                Euler_apply_transf(L[sym],R[sym],
                                                   bestEulerAngles[j2](0), bestEulerAngles[j2](1),
                                                   bestEulerAngles[j2](2), otherrot, othertilt,
                                                   otherpsi);
                                Euler_angles2matrix(otherrot,othertilt,otherpsi,E2);
                                double aux=Euler_distanceBetweenMatrices(E1,E2);
                                if (aux>bestCorrE1E2)
                                    bestCorrE1E2=aux;
                            }
                            if (bestCorrE1E2>0.97)
                            {
                                alreadyClustered(j2)=1;
                                if (bestCorrelation<bestCorrelations[j2])
                                {
                                    bestCorrelation=bestCorrelations[j2];
                                    bestAssignment=j2;
                                }
                            }
                        }
                    }

                    clusterBestAssignment.push_back(bestAssignment);
                    clusterBestCorrelation.push_back(bestCorrelation);

                    std::cout << "Cluster headed by "
                    << bestEulerAngles[bestAssignment].transpose()
                    << " corr= " << bestCorrelation << std::endl;
                }
            }

            if (clusterBestAssignment.size()>=5 && assigned.sum()>5)
            {
                alreadyOptimized(besti) = 0;
                assigned(besti)         = 0;
                tabuPenalization(besti)+=10;
            }
            else
            {
                // Set the status of this image to optimized
                alreadyOptimized(besti) = 2;
                assigned(besti)         = 1;

                // Try all the solutions in the cluster just to make sure
                int bestCluster=-1;
                double bestEnergy=0;
                Matrix1D<double> bestCurrentSolution,
                bestCurrentImgAvgCorrelation,
                bestCurrentImgMinCorrelation;
                MultidimArray<double> bestCurrentCorrelationMatrix;
                if (clusterBestAssignment.size()>1)
                {
                    Matrix1D<double> backupCurrentSolution;
                    backupCurrentSolution=currentSolution;
                    for (size_t n=0; n<clusterBestAssignment.size(); n++)
                    {
                        std::cout << "Trying solution of cluster " << n << std::endl;
                        currentSolution(3*besti)   = bestEulerAngles[
                                                         clusterBestAssignment[n]](0);
                        currentSolution(3*besti+1) = bestEulerAngles[
                                                         clusterBestAssignment[n]](1);
                        currentSolution(3*besti+2) = bestEulerAngles[
                                                         clusterBestAssignment[n]](2);
                        double energy=realignCurrentSolution();
                        if (energy<bestEnergy)
                        {
                            bestEnergy=energy;
                            bestCluster=n;
                            bestCurrentSolution=currentSolution;
                            bestCurrentImgAvgCorrelation=currentImgAvgCorrelation;
                            bestCurrentImgMinCorrelation=currentImgMinCorrelation;
                            bestCurrentCorrelationMatrix=currentCorrelationMatrix;
                        }
                        currentSolution=backupCurrentSolution;
                    }
                }
                else
                {
                    currentSolution(3*besti)   = bestEulerAngles[
                                                     clusterBestAssignment[0]](0);
                    currentSolution(3*besti+1) = bestEulerAngles[
                                                     clusterBestAssignment[0]](1);
                    currentSolution(3*besti+2) = bestEulerAngles[
                                                     clusterBestAssignment[0]](2);
                    double energy=realignCurrentSolution();
                    bestEnergy=energy;
                    bestCluster=0;
                    bestCurrentSolution=currentSolution;
                    bestCurrentImgAvgCorrelation=currentImgAvgCorrelation;
                    bestCurrentImgMinCorrelation=currentImgMinCorrelation;
                    bestCurrentCorrelationMatrix=currentCorrelationMatrix;
                }

                // Realign the current solution
                std::cout << "Cluster chosen " << bestCluster << " angles="
                <<  bestEulerAngles[clusterBestAssignment[
                                        bestCluster]].transpose() << std::endl;
                currentSolution=bestCurrentSolution;
                currentImgAvgCorrelation=bestCurrentImgAvgCorrelation;
                currentImgMinCorrelation=bestCurrentImgMinCorrelation;
                currentCorrelationMatrix=bestCurrentCorrelationMatrix;
                std::cout << "Image Avg Correlation: " << currentImgAvgCorrelation
                << "\nImage Min Correlation: " << currentImgMinCorrelation
                << std::endl;
            }

            // Cleaning of the "garbage"
            int totalAssigned=(int)assigned.sum();
            std::cout << "removal=" << removalCounter
            << " totalAssigned=" << totalAssigned
            << std::endl;
            removeViaClusters(currentCorrelationMatrix);
            if (removalCounter!=0 && totalAssigned>3)
            {
                std::vector<int> imgsToRemove;
                int imin=-1;
                double worseCorrelation=2;
                double bestCorrelation=-2;
                if (removalCounter==0)
                {
                    double meanCorr=0;
                    double stdCorr=0;
                    double Ncorr=0;

                    FOR_ALL_ELEMENTS_IN_MATRIX1D(currentImgAvgCorrelation)
                    {
                        if (currentImgAvgCorrelation(i)==1)
                        {
                            meanCorr+=currentImgAvgCorrelation(i);
                            stdCorr+=currentImgAvgCorrelation(i)*
                                     currentImgAvgCorrelation(i);
                            Ncorr++;
                            switch (removalCounter)
                            {
                            case 1:
                                if (currentImgAvgCorrelation(i)<worseCorrelation)
                                {
                                    worseCorrelation=currentImgAvgCorrelation(i);
                                    imin=i;
                                }
                                break;
                            case 2:
                                if (currentImgAvgCorrelation(i)>bestCorrelation)
                                {
                                    bestCorrelation=currentImgAvgCorrelation(i);
                                    imin=i;
                                }
                                break;
                            }
                        }
                    }
                    imgsToRemove.push_back(imin);
                    if (Ncorr>0)
                    {
                        meanCorr/=Ncorr;
                        stdCorr=sqrt(stdCorr/Ncorr-meanCorr*meanCorr);
                        std::cout << "Mean=" << meanCorr << " std="
                        << stdCorr << std::endl;
                        if (ABS(worseCorrelation-meanCorr)>3*stdCorr)
                            imgsToRemove.push_back(imin);
                    }
                }
                else
                    imgsToRemove=removeViaClusters(currentCorrelationMatrix);

                if (imgsToRemove.size()!=0)
                {
                    for (size_t n=0; n<imgsToRemove.size(); n++)
                    {
                        int imin=imgsToRemove[n];
                        std::cout << "Image " << imin << " removed from the "
                        << "current assignment corr="
                        << currentImgAvgCorrelation(imin)
                        << std::endl;

                        alreadyOptimized(imin)=assigned(imin)=0;
                        currentSolution(3*imin)=currentSolution(3*imin+1)=
                                                    currentSolution(3*imin+2)=0;
                        currentImgAvgCorrelation(imin)=0;
                        currentImgMinCorrelation(imin)=2;
                        for (size_t i=0; i<XSIZE(currentCorrelationMatrix); i++)
                            currentCorrelationMatrix(i,imin)=
                                currentCorrelationMatrix(imin,i)=0;
                        tabuPenalization(imin)+=10;
                    }
                }

            }
            removalCounter=(removalCounter+1)%3;
        }

        // Remove one from the penalization of every image
        FOR_ALL_ELEMENTS_IN_MATRIX1D(tabuPenalization)
        if (tabuPenalization(i)>0)
            tabuPenalization(i)--;
    }

    solution=currentSolution;
}

#ifdef NEVERDEFINED
void Prog_Angular_CommonLine::optimize(Matrix1D<double> &solution)
{
    int Nimg = SF.ImgNo();
    alreadyOptimized.initZeros(Nimg);
    alreadyOptimized.initConstant(1);
    alreadyOptimized(0)=2;
    Matrix1D<int> imgIdx(Nimg-1);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(imgIdx)
    imgIdx(i)=i+1;
    double energy=optimizeGroup(imgIdx,solution,true);
}
#endif

// Compute clusters
double Prog_Angular_CommonLine::computeClusters(
    const MultidimArray<double> &correlationMatrix,
    std::vector< std::vector<int> > &clusters,
    MultidimArray<double> &worseCorrelationMatrix,
    MultidimArray<double> &bestCorrelationMatrix, bool show) const
{
    // Extract the list of images currently optimized and a smaller
    // correlation matrix
    std::vector<int> idxImgs;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
    if (alreadyOptimized(i)!=0)
        idxImgs.push_back(i);
    int Noptimized=idxImgs.size();

    // Initially every element is a cluster
    worseCorrelationMatrix.initZeros(Noptimized,Noptimized);
    clusters.clear();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(worseCorrelationMatrix)
    worseCorrelationMatrix(i,j)=correlationMatrix(idxImgs[i],idxImgs[j]);
    bestCorrelationMatrix=worseCorrelationMatrix;
    for (int i=0; i<Noptimized; i++)
    {
        std::vector<int> singleElement;
        singleElement.clear();
        singleElement.push_back(i);
        clusters.push_back(singleElement);
    }
    std::cout << "Correlation matrix\n"
    << worseCorrelationMatrix << std::endl;

    // Compute the average distance between all elements
    double avgDistance=0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(worseCorrelationMatrix)
    if (i!=j)
        avgDistance+=worseCorrelationMatrix(i,j);
    else
        avgDistance+=1;
    avgDistance/=Noptimized*Noptimized;

    // Start the real clustering
    while (clusters.size()>2)
    {
        // Look for the two closest clusters
    	size_t besti=0, bestj=0;
        double bestCorr=-1;
        for (size_t i=0; i<YSIZE(worseCorrelationMatrix); i++)
            for (size_t j=i+1; j<XSIZE(worseCorrelationMatrix); j++)
                if (worseCorrelationMatrix(i,j)>bestCorr)
                {
                    bestCorr=worseCorrelationMatrix(i,j);
                    besti=i;
                    bestj=j;
                }

        // std::cout << "Joining " << besti << " and " << bestj << std::endl;
        // std::cout << "Worse Correlation Matrix:\n" << worseCorrelationMatrix << std::endl;
        // std::cout << "Best Correlation Matrix:\n" << bestCorrelationMatrix << std::endl;

        // Join them in the list
        for (size_t n=0; n<clusters[bestj].size(); n++)
            clusters[besti].push_back(clusters[bestj][n]);
        clusters.erase(clusters.begin()+bestj);

        // Readjust the distance between this new cluster and the rest
        // of the existing
        for (size_t i=0; i<YSIZE(worseCorrelationMatrix); i++)
        {
            if (i!=besti)
            {
                worseCorrelationMatrix(besti,i)=
                    worseCorrelationMatrix(i,besti)=XMIPP_MIN(
                                                        worseCorrelationMatrix(i,besti),
                                                        worseCorrelationMatrix(i,bestj));
                bestCorrelationMatrix(besti,i)=
                    bestCorrelationMatrix(i,besti)=XMIPP_MAX(
                                                       bestCorrelationMatrix(i,besti),
                                                       bestCorrelationMatrix(i,bestj));
            }
            else
            {
                worseCorrelationMatrix(besti,besti)=
                    worseCorrelationMatrix(bestj,bestj)=
                        worseCorrelationMatrix(besti,bestj);
                bestCorrelationMatrix(besti,besti)=
                    bestCorrelationMatrix(bestj,bestj)=
                        bestCorrelationMatrix(besti,bestj);
            }
        }

        // Move everything from bestj to the left
        for (size_t i=0; i<YSIZE(worseCorrelationMatrix); i++)
            for (size_t j=bestj; j<XSIZE(worseCorrelationMatrix)-1; j++)
            {
                worseCorrelationMatrix(i,j)=worseCorrelationMatrix(i,j+1);
                bestCorrelationMatrix(i,j)=bestCorrelationMatrix(i,j+1);
            }

        // Move everything from bestj to the top
        for (size_t i=bestj; i<YSIZE(worseCorrelationMatrix)-1; i++)
            for (size_t j=0; j<XSIZE(worseCorrelationMatrix); j++)
            {
                worseCorrelationMatrix(i,j)=worseCorrelationMatrix(i+1,j);
                bestCorrelationMatrix(i,j)=bestCorrelationMatrix(i+1,j);
            }

        // Remove the last row and column of worseCorrelation
        worseCorrelationMatrix.resize(YSIZE(worseCorrelationMatrix)-1,
                                      XSIZE(worseCorrelationMatrix)-1);
        bestCorrelationMatrix.resize(YSIZE(bestCorrelationMatrix)-1,
                                     XSIZE(bestCorrelationMatrix)-1);
    }

    // Substitute the cluster indexes by image indexes
    for (size_t n=0; n<2; n++)
        for (size_t i=0;i<clusters[n].size(); i++)
            clusters[n][i]=idxImgs[clusters[n][i]];

    // Compute the separability
    Matrix1D<double> avgDistancek;
    avgDistancek.initZeros(2);
    for (size_t n=0; n<2; n++)
    {
        for (size_t i=0;i<clusters[n].size(); i++)
            for (size_t j=0;j<clusters[n].size(); j++)
                if (i!=j)
                    avgDistancek(n)+=correlationMatrix(
                                         clusters[n][i],clusters[n][j]);
                else
                    avgDistancek(n)+=1;
        avgDistancek(n)/=clusters[n].size()*clusters[n].size();
    }
    double mergeDistance1=avgDistance/(0.5*(avgDistancek(0)+avgDistancek(1)));
    double mergeDistance2;
    if (worseCorrelationMatrix(0,0)>0 && worseCorrelationMatrix(1,1)>0)
        mergeDistance2=worseCorrelationMatrix(0,1)/
                       (0.5*(avgDistancek(0)+avgDistancek(1)));
    else
        mergeDistance2=0;

    if (show)
    {
        for (int n=0; n<2; n++)
        {
            std::cout << "Cluster " << n << ": ";
            for (size_t i=0;i<clusters[n].size(); i++)
                std::cout << clusters[n][i] << " ";
            std::cout << std::endl;
        }
        std::cout << "Merge distance=" << mergeDistance1 << " ="
        << avgDistance << "/(0.5*(" << avgDistancek(0) << "+"
        << avgDistancek(1) << ")\n";
        std::cout << "Merge distance=" << mergeDistance2 << " ="
        << worseCorrelationMatrix(0,1) << "/(0.5*("
        << avgDistancek(0) << "+"
        << avgDistancek(1) << ")\n";
        std::cout << "Worse Correlation matrix\n" << worseCorrelationMatrix << std::endl;
        std::cout << "Best Correlation matrix\n" << bestCorrelationMatrix << std::endl;
    }
    return mergeDistance2;
}

// Performs a HCA to see which image to remove
std::vector<int> Prog_Angular_CommonLine::removeViaClusters(
    const MultidimArray<double> &correlationMatrix)
{
    std::vector<int> retval;
    std::vector< std::vector<int> > clusters;
    MultidimArray<double> worseCorrelationMatrix;
    MultidimArray<double> bestCorrelationMatrix;
    computeClusters(correlationMatrix,clusters,
                    worseCorrelationMatrix, bestCorrelationMatrix, true);
    if (clusters[0].size()+clusters[1].size()<=4)
        return retval;

    // Find the less populated cluster
    int nmin=0;
    int nmax=1;
    if (clusters[0].size()>clusters[1].size())
    {
        nmin=1;
        nmax=0;
    }

    // Choose the element or elements to remove
    double diameter0=(bestCorrelationMatrix(0,0)-worseCorrelationMatrix(0,0));
    double diameter1=(bestCorrelationMatrix(1,1)-worseCorrelationMatrix(1,1));
    double diameter01=(bestCorrelationMatrix(0,1)-worseCorrelationMatrix(0,1));
    bool separated=(diameter01>(diameter0+diameter1));

    if (separated && clusters[nmin].size()<=3 && clusters[nmax].size()>=6)
    {
        for (size_t i=0;i<clusters[nmin].size(); i++)
        {
            int imin=clusters[nmin][i];
            retval.push_back(imin);
        }
    }
    else if (clusters[nmax].size()>=4 && clusters[nmin].size()>=3)
    {
        // Look for the worse image within this cluster
        int imin=-1;
        double worseCorrelation=2;
        for (size_t i=0;i<clusters[nmin].size(); i++)
        {
            int imgIndex=clusters[nmin][i];
            if (currentImgAvgCorrelation(imgIndex)<worseCorrelation)
            {
                worseCorrelation=currentImgAvgCorrelation(imgIndex);
                imin=imgIndex;
            }
        }
        retval.push_back(imin);
    }

    // Return which image to remove
    std::cout << "Remove images: ";
    for (size_t i=0; i< retval.size(); i++)
        std::cout << retval[i] << " ";
    std::cout << std::endl;
    return retval;
}

double Prog_Angular_CommonLine::realignCurrentSolution()
{
    // Realign all images that have already been optimized
    FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
    if (alreadyOptimized(i)==2)
        alreadyOptimized(i)=1;
    int NToSolve=(int)alreadyOptimized.sum();
    Matrix1D<int> imgIdx(NToSolve);
    int idx=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
    if (alreadyOptimized(i)==1)
        imgIdx(idx++)=i;

    solver=new EulerSolver(3*NToSolve,30*NToSolve,
                           alreadyOptimized, currentSolution, imgIdx, this);
    global_Eulersolver=solver;

    Matrix1D<double> steps(3*NToSolve), solution;
    steps.initConstant(1);
    solution.initZeros(steps);
    idx=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(imgIdx)
    {
        solution(idx++)=currentSolution(3*imgIdx(i));
        solution(idx++)=currentSolution(3*imgIdx(i)+1);
        solution(idx++)=currentSolution(3*imgIdx(i)+2);
    }
    int iter;
    double energy;
    powellOptimizer(solution,1,3*NToSolve,wrapperSolverEnergy,NULL,
                    0.001,energy,iter,steps,true);

    idx=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(imgIdx)
    {
        currentSolution(3*imgIdx(i))   = solution(idx++);
        currentSolution(3*imgIdx(i)+1) = solution(idx++);
        currentSolution(3*imgIdx(i)+2) = solution(idx++);
    }

    FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
    if (alreadyOptimized(i)==1)
        alreadyOptimized(i)=2;

    currentImgAvgCorrelation=solver->imgAvgCorrelation.transpose();
    currentImgMinCorrelation=solver->imgMinCorrelation.transpose();
    currentCorrelationMatrix=solver->correlationMatrix;
    delete solver;
    return energy;
}

/* Try solution ------------------------------------------------------------ */
double Prog_Angular_CommonLine::trySolution(const Matrix1D<double> &solution)
{
    Matrix1D<double> minAllowed(VEC_XSIZE(solution)), maxAllowed(VEC_XSIZE(solution));
    int idx=0;
    int Nimg=VEC_XSIZE(solution)/3;
    for (int i=0; i<Nimg; i++)
    {
        maxAllowed(idx++)=360;
        maxAllowed(idx++)=180;
        maxAllowed(idx++)=360;
    }

    alreadyOptimized.initZeros(Nimg);
    alreadyOptimized.initConstant(1);
    MultidimArray<int> imgIdx(Nimg);
    imgIdx.initLinear(0,Nimg-1,1,"incr");
    bool bAtSolution;
    solver=new EulerSolver(VEC_XSIZE(solution),1,
                           alreadyOptimized,currentSolution,imgIdx,this);
    solver->Setup(MATRIX1D_ARRAY(minAllowed), MATRIX1D_ARRAY(maxAllowed),
                  stBest2Bin, 0.5, 0.8);
    solver->setShow(true);
    double energy=solver->EnergyFunction(MATRIX1D_ARRAY(solution),bAtSolution);
    std::cout << "Correlation matrix\n" << solver->correlationMatrix
    << std::endl;
    return energy;
}

/* Run -------------------------------------------------------------------- */
void Prog_Angular_CommonLine::run()
{
    if (tryInitial)
    {
        trySolution(initialSolution);
    }
    else
    {
        // Look for the solution
        Matrix1D<double> solution;
        optimize(solution);

        MetaData DF;

        int idx=0;
        Image<double> I;
        size_t id;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            I.readApplyGeo(SF,__iter.objId);
            I.setEulerAngles(currentSolution(idx),currentSolution(idx+1),
            		currentSolution(idx+2));
            idx+=3;
            I.write();

            id=DF.addObject();
            DF.setValue(MDL_ANGLE_ROT,I.rot(),id);
            DF.setValue(MDL_ANGLE_TILT,I.tilt(),id);
            DF.setValue(MDL_ANGLE_PSI,I.psi(),id);
        }
        DF.write(fnOut);
    }
}
