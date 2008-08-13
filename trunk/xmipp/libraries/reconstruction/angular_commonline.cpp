/***************************************************************************
 *
 * Authors:     Ignacio Fernández Moreno (hellmoon666@gmail.com)
 *              Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

#include "angular_commonline.h"
#include <data/filters.h>
#include <data/geometry.h>
#include <data/image.h>
#include <data/docfile.h>
#include <reconstruction/radon.h>

/* Constructor ------------------------------------------------------------- */
EulerSolver::EulerSolver(int dim, int pop,
    const Matrix1D<int> &newAlreadyOptimized,
    const Matrix1D<double> &newCurrentSolution,
    const Matrix1D<int> &newImgIdx,
    const Prog_Angular_CommonLine *newParent): DESolver(dim, pop)
{
    parent = newParent;
    Ndim = dim;
    NToSolve=dim/3;
    Nimg=XSIZE(newAlreadyOptimized);
    imgCorrelation.initZeros(Nimg);
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

    for (i=0, j=0; i<Ndim; i++,j=(j+1)%3) {
        if (j==1) trial[i]=realWRAP(trial[i],0,180);
        else      trial[i]=realWRAP(trial[i],0,360);
    }

    // Get number of symmetries
    int Nsym=parent->L.size();

    // Evaluate goal function for this trial
    double retval=0;
    double Ncomparisons=0;
    double worseSimilarity=2;
    int idx;
    static Matrix1D<int> imgCorrelationN;
    imgCorrelation.initZeros(Nimg);
    imgCorrelationN.initZeros(Nimg);
    for (int imgi=0; imgi<Nimg; imgi++)
    {
        // Get the right angles for this image
        if      ((*alreadyOptimized)(imgi)==0) continue;
        else if ((*alreadyOptimized)(imgi)==1)
        {
            for (idx=0; idx<NToSolve; idx++)
                if ((*imgIdx)(idx)==imgi) break;
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
            if      ((*alreadyOptimized)(imgj)==0) continue;
            else if ((*alreadyOptimized)(imgj)==1)
            {
                for (idx=0; idx<NToSolve; idx++)
                    if ((*imgIdx)(idx)==imgj) break;
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
            if ((*alreadyOptimized)(imgj)==2 && (*alreadyOptimized)(imgj)==2)
                continue;

            double similarity=similarityBetweenTwoLines(imgi,imgj);
            if (similarity>0)
            {
                retval -= similarity;
                Ncomparisons++;
                worseSimilarity=XMIPP_MIN(worseSimilarity,similarity);
                imgCorrelation(imgi)+=similarity;
                imgCorrelation(imgj)+=similarity;
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
                    imgCorrelation(imgi)+=similarity;
                    imgCorrelation(imgj)+=similarity;
                    imgCorrelationN(imgi)++;
                    imgCorrelationN(imgj)++;
                }
            }
	}	
    }
    if (Ncomparisons>0) retval/=Ncomparisons;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(imgCorrelation)
        if (imgCorrelationN(i)!=0) imgCorrelation(i)/=imgCorrelationN(i);
    if (show)
        std::cout
            << "Average Distance = " << retval << std::endl
            << "Worse   Distance = " << -worseSimilarity << std::endl
            << "Final   Distance = " << 0.5*(retval-worseSimilarity)
            << "imgCorrelation = " << imgCorrelation.transpose() << std::endl
            << std::endl
            << std::endl
        ;
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
            char c; std::cin >> c;
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
    
    double retval1=0.5*(
        correlation_index(parent->radon[imgi][idxAngi],
                          parent->radon[imgj][idxAngj])+
        correlation_index(parent->radonDerivative[imgi][idxAngi],
                          parent->radonDerivative[imgj][idxAngj]));

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
            << "Distance between lines = " << retval1 << std::endl
        ;
        parent->radon[imgi][idxAngi].write("PPPradoni1.txt");
        parent->radon[imgj][idxAngj].write("PPPradonj1.txt");
        parent->radonDerivative[imgi][idxAngi].write("PPPradonDerivativei1.txt");
        parent->radonDerivative[imgj][idxAngj].write("PPPradonDerivativej1.txt");
    }

    // Try now with the opposite direction
    idxAngi = (int)intWRAP(-((int)angi)+180,0,359);
    double retval2=0.5*(
        correlation_index(parent->radon[imgi][idxAngi],
                          parent->radon[imgj][idxAngj])+
        correlation_index(parent->radonDerivative[imgi][idxAngi],
                          parent->radonDerivative[imgj][idxAngj]));

    if (show)
    {
        Matrix2D<double> Euleri, Eulerj;
        Euler_angles2matrix(roti,tilti,psii,Euleri);
        Euler_angles2matrix(rotj,tiltj,psij,Eulerj);
        std::cout
            << "in imgi=" << commonlinei.transpose() << " anglei=" << angi
            << " (" << idxAngi << ")\n"
            << "Distance between lines = " << retval2 << std::endl
            << "Euler i" << Euleri << std::endl
            << "Euler j" << Eulerj << std::endl
            << std::endl
        ;
        parent->radon[imgi][idxAngi].write("PPPradoni2.txt");
        parent->radon[imgj][idxAngj].write("PPPradonj2.txt");
        parent->radonDerivative[imgi][idxAngi].write("PPPradonDerivativei2.txt");
        parent->radonDerivative[imgj][idxAngj].write("PPPradonDerivativej2.txt");
        std::cout << "Press any key\n";
        char c; std::cin >> c;
    }

    #ifdef DEBUG
    for (idxAngi=0; idxAngi<360; idxAngi++)
    {
        double retval3=0.5*(
            correlation_index(parent->radon[imgi][idxAngi],
                              parent->radon[imgj][idxAngj])+
            correlation_index(parent->radonDerivative[imgi][idxAngi],
                              parent->radonDerivative[imgj][idxAngj]));

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
            char c; std::cin >> c;
            if (c=='q') break;
        }
    }
    #endif
    return XMIPP_MAX(retval1,retval2);
}
#undef DEBUG

/* Wrapper for Powell ------------------------------------------------------ */
EulerSolver *global_Eulersolver=NULL;

double wrapperSolverEnergy(double trial[]) {
   bool bAtSolution;
   return global_Eulersolver->EnergyFunction(trial+1,bAtSolution);
}

/* Parameters -------------------------------------------------------------- */
void Prog_Angular_CommonLine::read(int argc, char **argv)
{
    fnSel = getParameter(argc,argv,"-i");
    fnOut = getParameter(argc,argv,"-o");
    fnSym = getParameter(argc,argv,"-sym","");
    NGen   = textToInteger(getParameter(argc,argv,"-NGen","50000"));
    NGroup = textToInteger(getParameter(argc,argv,"-NGroup","7"));
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
              << "   -o <docfile>     : Docfile with the angular assignment\n"
              << "  [-NGen <g=50000>] : Number of generations\n"
              << "  [-NGroup <g=7>]   : Number of generations\n"
              << "  [-tryInitial]     : Do not optimize\n"
              << "  [-sym <symfile>]  : Symmetry\n"
    ;
}

/* Produce side info ------------------------------------------------------- */
void Prog_Angular_CommonLine::produceSideInfo()
{
    // Read selfile images and the initial angles
    SF.read(fnSel);
    int Nimg=SF.ImgNo();
    initialSolution.resize(3*Nimg);
    int idx=0;
    while (!SF.eof())
    {
        ImageXmipp I;
        I.read(SF.NextImg());
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
        SL.read_sym_file(fnSym);
        for (int sym=0; sym<SL.SymsNo(); sym++)
        {
            Matrix2D<double> auxL, auxR;
            SL.get_matrices(sym,auxL,auxR);
            auxL.resize(3,3);
            auxR.resize(3,3);
            L.push_back(auxL);
            R.push_back(auxR);
        }
    }

    // Compute the Radon transform and its derivative of each image
    const double deltaRot=1.0;
    Matrix2D<double> RT;
    Matrix1D<double> projection, projectionDerivative;
    std::cout << "Preprocessing images ...\n";
    init_progress_bar(Nimg);
    for (int n=0; n<Nimg; n++)
    {
        // Initialize list of projection and derivatives for this image
        std::vector< Matrix1D<double> > dummyList;
        radon.push_back(dummyList);
        radonDerivative.push_back(dummyList);

        // Compute Radon transform of each image
        Radon_Transform(img[n],deltaRot,RT);

        // Separate each projection line and compute derivative
        for (int i=0; i<YSIZE(RT); i++)
        {
            RT.getRow(i,projection);
            projection.numericalDerivative(projectionDerivative);
            radon[n].push_back(projection);
            radonDerivative[n].push_back(projectionDerivative);
        }
        progress_bar(n);
    }
    progress_bar(Nimg);
}

/* OptimizeGroup ----------------------------------------------------------- */
double Prog_Angular_CommonLine::optimizeGroup(const Matrix1D<int> &imgIdx,
    Matrix1D<double> &solution, bool show) {
    int    count_repetitions = 0;
    double current_energy = 0;
    double previous_energy = 2;
    double first_energy = 2;
    int    NGenStep=NGen/100;
    int    Nimg=SF.ImgNo();
    int    NToSolve=XSIZE(imgIdx);

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
    solver->Setup(MULTIDIM_ARRAY(minAllowed), MULTIDIM_ARRAY(maxAllowed),
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

        if (current_energy<0 && first_energy>0) first_energy=current_energy;

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
	if (count_repetitions>=10 && count>=50 || count==100)
            done=true;

	previous_energy = current_energy;
    }

    // Optimize with Powell
    Matrix1D<double> steps(3*NToSolve);
    steps.initConstant(1);

    solution.initZeros(steps);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(solution)
   	DIRECT_VEC_ELEM(solution,i) = solver->Solution()[i];
    int iter;
    double retval;
    powellOptimizer(solution,1,3*NToSolve,wrapperSolverEnergy,
        0.001,retval,iter,steps,show);

    delete solver;
    return retval;
}

/* Optimize ---------------------------------------------------------------- */
void Prog_Angular_CommonLine::optimize(Matrix1D<double> &solution)
{
    int Nimg = SF.ImgNo();
    solution.initZeros(3*Nimg);
    Matrix1D<int> assigned;
    assigned.initZeros(Nimg);
    assigned(0)=1;

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
        for (int i=0; i<Nimg; i++)
        {
            std::vector< Matrix1D<double> > dummy1;
            std::vector< double > dummy2;
            eulerAngles.push_back(dummy1);
            correlations.push_back(dummy2);
        }

        // For each image try different pairs and get different alignments
        std::cout << "Aligning image pairs, "
                  << Nimg-assigned.sum() << " images remaining ...\n";
        init_progress_bar(Nimg);
        for (int i=1; i<Nimg; i++)
        {
            if (assigned(i)) continue;
            Matrix1D<int> backupAlreadyOptimized=alreadyOptimized;

            Matrix1D<int> comparedTo;
            comparedTo=assigned;
            comparedTo(i)=1;
            bool allowComparisonWithFixedImages=false;
            // Perform NGroup experiments
            for (int n=0; n<NGroup; n++)
            {
                if (comparedTo.sum()==XSIZE(comparedTo))
                    allowComparisonWithFixedImages=true;

                // Look for an image to which we have not compared
                int j;
                bool tryAgain=true;
                do {
                    j=ROUND(rnd_unif(0,Nimg-1));
                    if (!allowComparisonWithFixedImages && comparedTo(j)==0)
                        tryAgain=false;
                    else if (allowComparisonWithFixedImages && comparedTo(j)==1)
                        tryAgain=false;
                } while (tryAgain);
                comparedTo(j)++;

                // Prepare the vector of images to optimize
                alreadyOptimized=backupAlreadyOptimized;
                alreadyOptimized(i)=alreadyOptimized(j)=1;

                // Align these two images
                Matrix1D<int> imgIdx(2);
                VECTOR_R2(imgIdx,i,j);
                Matrix1D<double> auxSolution, anglesi(3), anglesj(3);
                double energy=optimizeGroup(imgIdx,auxSolution,false);
                anglesi(0)=auxSolution(0);
                anglesi(1)=auxSolution(1);
                anglesi(2)=auxSolution(2);
                anglesj(0)=auxSolution(3);
                anglesj(1)=auxSolution(4);
                anglesj(2)=auxSolution(5);

                // Keep results
                eulerAngles[i].push_back(anglesi);
                correlations[i].push_back(-energy);
                if (!allowComparisonWithFixedImages)
                {
                    eulerAngles[j].push_back(anglesj);
                    correlations[j].push_back(-energy);
                }
                alreadyOptimized=backupAlreadyOptimized;
            }
            progress_bar(i);
        }
        progress_bar(Nimg);

        // Compute for each image the variance in the top assignment
        int Nsym=SL.SymsNo();
        double bestDistance=0;
        int besti=-1;
        int topN=NGroup;
        for (int i=1; i<Nimg; i++)
            if (eulerAngles[i].size()<2*topN && !assigned(i))
                topN=CEIL(eulerAngles[i].size()/2);
        for (int i=1; i<Nimg; i++)
        {
            // If the particle has  already been assigned skip it
            if (assigned(i)) continue;
            
            // Sort the images by ascending correlation
            Matrix1D<double> aux;
            aux.initZeros(eulerAngles[i].size());
            for (int n=0; n<eulerAngles[i].size(); n++)
                aux(n)=correlations[i][n];
            Matrix1D<int> idx=aux.indexSort();

            // Among the top, compute the distance between
            // the different Euler angles
            double distance=0;
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
                        double maxDistance=ABS(E1(2,0)*E2(2,0)+
                            E1(2,1)*E2(2,1)+E1(2,2)*E2(2,2));
                        for (int sym=0; sym<Nsym; sym++)
                        {
                            Euler_apply_transf(L[sym],R[sym],
                                rot, tilt, psi, otherrot, othertilt, otherpsi);
                            Euler_angles2matrix(otherrot,othertilt,otherpsi,E2);
                            double symDistance=ABS(E1(2,0)*E2(2,0)+
                                E1(2,1)*E2(2,1)+E1(2,2)*E2(2,2));
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

        // Sort the images by ascending correlation in the best cluster
        Matrix1D<double> aux;
        aux.initZeros(eulerAngles[besti].size());
        for (int n=0; n<eulerAngles[besti].size(); n++)
            aux(n)=correlations[besti][n];
        Matrix1D<int> idx=aux.indexSort();

        // Keep only the topN
        std::vector< Matrix1D<double> > bestEulerAngles;
        std::vector< double > bestCorrelations;
        for (int n=XSIZE(idx)-topN; n<XSIZE(idx); n++)
        {
            bestEulerAngles.push_back(eulerAngles[besti][idx(n)-1]);
            bestCorrelations.push_back(correlations[besti][idx(n)-1]);
        }

        // Show best resolved image
        std::cout << "Candidates for image " << besti << " distance="
                  << bestDistance << std::endl;
        for (int n=0; n<bestEulerAngles.size(); n++)
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
        for (int j1=0; j1<topN; j1++)
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
                for (int j2=j1+1; j2<topN; j2++)
                {
                    if (!alreadyClustered(j2))
                    {
                        Matrix2D<double> E2;
                        Euler_angles2matrix(bestEulerAngles[j2](0),
                            bestEulerAngles[j2](1),
                            bestEulerAngles[j2](2),E2);
                        double bestCorrE1E2=
                            Euler_distanceBetweenMatrices(E1,E2);
                        for (int sym=0; sym<Nsym; sym++)
                        {
                            double otherrot, othertilt, otherpsi;
                            Euler_apply_transf(L[sym],R[sym],
                                bestEulerAngles[j2](0), bestEulerAngles[j2](1),
                                bestEulerAngles[j2](2), otherrot, othertilt,
                                otherpsi);
                            Euler_angles2matrix(otherrot,othertilt,otherpsi,E2);
                            double aux=Euler_distanceBetweenMatrices(E1,E2);
                            if (aux>bestCorrE1E2) bestCorrE1E2=aux;
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

        // Set the status of these images to optimized 
        alreadyOptimized(besti) = 2;
        assigned(besti)         = 1;
        
        // Try all the solutions in the cluster just to make sure
        int bestCluster=-1;
        double bestEnergy=0;
        Matrix1D<double> bestCurrentSolution, bestCurrentImageCorrelation;
        if (clusterBestAssignment.size()>1)
        {
            Matrix1D<double> backupCurrentSolution;
            backupCurrentSolution=currentSolution;
            for (int n=0; n<clusterBestAssignment.size(); n++)
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
                    bestCurrentImageCorrelation=currentImageCorrelation;
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
            bestCurrentImageCorrelation=currentImageCorrelation;
        }
        
        // Realign the current solution
        std::cout << "Cluster chosen " << bestCluster << " angles="
                  <<  bestEulerAngles[clusterBestAssignment[
                        bestCluster]].transpose() << std::endl;
        currentSolution=bestCurrentSolution;
        currentImageCorrelation=bestCurrentImageCorrelation;
        std::cout << "Image Correlation: " << currentImageCorrelation
                  << std::endl;
        
        // Every two images remove the worse one in the list
        int totalAssigned=assigned.sum();
        removalCounter=(removalCounter+1)%2;
        if (totalAssigned>4 && removalCounter==0)
        {
            int imin=-1;
            double worseCorrelation=2;
            FOR_ALL_ELEMENTS_IN_MATRIX1D(currentImageCorrelation)
                if (currentImageCorrelation(i)>0 &&
                    currentImageCorrelation(i)<worseCorrelation &&
                    i>0)
                {
                    worseCorrelation=currentImageCorrelation(i);
                    imin=i;
                }
            if (imin!=-1)
            {
                alreadyOptimized(imin)=assigned(imin)=0;
                currentSolution(3*imin)=currentSolution(3*imin+1)=
                    currentSolution(3*imin+2)=0;
                std::cout << "Image " << imin << " removed from the "
                          << "current assignment because its correlation was "
                          << worseCorrelation << std::endl;
            }
        }
    }
    
    solution=currentSolution;
}

double Prog_Angular_CommonLine::realignCurrentSolution()
{
    // Realign all images that have already been optimized
    FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
        if (alreadyOptimized(i)==2) alreadyOptimized(i)=1;
    alreadyOptimized(0)=2;
    int NToSolve=alreadyOptimized.sum()-2;
    Matrix1D<int> imgIdx(NToSolve);
    int idx=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
        if (alreadyOptimized(i)==1) imgIdx(idx++)=i;

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
    powellOptimizer(solution,1,3*NToSolve,wrapperSolverEnergy,
        0.001,energy,iter,steps,true);

    idx=0;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(imgIdx)
    {
   	 currentSolution(3*imgIdx(i))   = solution(idx++);
   	 currentSolution(3*imgIdx(i)+1) = solution(idx++);
   	 currentSolution(3*imgIdx(i)+2) = solution(idx++);
    }

    FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
        if (alreadyOptimized(i)==1) alreadyOptimized(i)=2;

    currentImageCorrelation=solver->imgCorrelation.transpose();
    delete solver;
    return energy;
}

#ifdef NEVER_DEFINED
void Prog_Angular_CommonLine::optimize(Matrix1D<double> &solution)
{
    int Nimg = SF.ImgNo();
    solution.initZeros(3*Nimg);

    // Compute the number of groups to try
    std::vector< Matrix1D<int> > groupsToTry=groups;

    // Align all groups
    while (groupsToTry.size()!=0)
    {
        // Align each group independently and annotate the best score
        Matrix1D<double> scores;
        scores.initZeros(groupsToTry.size());
        Matrix1D<int> backupAlreadyOptimized=alreadyOptimized;
        std::vector< Matrix1D<double> > individualSolutions;
        individualSolutions.clear();
        std::vector< Matrix1D<int> > individualImgIdx;
        individualImgIdx.clear();
        FOR_ALL_ELEMENTS_IN_MATRIX1D(scores)
        {
            // Prepare the set of images to optimize
            Matrix1D<int> imgIdx;
            imgIdx=groupsToTry[i];
            int NToSolve=XSIZE(imgIdx);
            std::cout << "Processing group " << imgIdx.transpose() << std::endl;

            // Prepare the vector of images to optimize
            alreadyOptimized=backupAlreadyOptimized;
            FOR_ALL_ELEMENTS_IN_MATRIX1D(imgIdx)
                alreadyOptimized(imgIdx(i))=1;

            // Really optimize
            Matrix1D<double> auxSolution;
            
            std::cout << "Already Optimized " << alreadyOptimized.transpose()
                      << std::endl;
            std::cout << "Current solution" << currentSolution.transpose()
                      << std::endl;
            
            scores(i)=optimizeGroup(imgIdx,auxSolution,true);
            individualSolutions.push_back(auxSolution);
            individualImgIdx.push_back(imgIdx);
            alreadyOptimized=backupAlreadyOptimized;
        }
        
        // Take the best aligned group
        int imin;
        std::cout << "Scores=" << scores.transpose() << std::endl;
        scores.minIndex(imin);
        std::cout << "The best aligned group is "
                  << individualImgIdx[imin].transpose()
                  << " with a score of "
                  << scores(imin) << std::endl
                  << "Solution " << individualSolutions[imin].transpose()
                  << std::endl;

        // Set the status of these images to optimized and copy
        // the solution
        FOR_ALL_ELEMENTS_IN_MATRIX1D(individualImgIdx[imin])
        {
            alreadyOptimized(individualImgIdx[imin](i))    = 2;
            currentSolution(3*individualImgIdx[imin](i))   =
                individualSolutions[imin](3*i);
            currentSolution(3*individualImgIdx[imin](i)+1) =
                individualSolutions[imin](3*i+1);
            currentSolution(3*individualImgIdx[imin](i)+2) =
                individualSolutions[imin](3*i+2);
        }
        
        // Remove the best group from the list of groups to try
        std::vector < Matrix1D<int> >::iterator groupIterator=
            groupsToTry.begin();
        for (int i=0; i<imin; i++)
            groupIterator++;
        groupsToTry.erase(groupIterator);

        // Realign all images that have already been optimized
        FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
            if (alreadyOptimized(i)==2) alreadyOptimized(i)=1;
        alreadyOptimized(0)=2;
        int NToSolve=alreadyOptimized.sum()-2;
        Matrix1D<int> imgIdx(NToSolve);
        int idx=0;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
            if (alreadyOptimized(i)==1) imgIdx(idx++)=i;
        
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
        powellOptimizer(solution,1,3*NToSolve,wrapperSolverEnergy,
            0.001,energy,iter,steps,true);

        idx=0;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(imgIdx)
        {
   	     currentSolution(3*imgIdx(i))   = solution(idx++);
   	     currentSolution(3*imgIdx(i)+1) = solution(idx++);
   	     currentSolution(3*imgIdx(i)+2) = solution(idx++);
        }

        FOR_ALL_ELEMENTS_IN_MATRIX1D(alreadyOptimized)
            if (alreadyOptimized(i)==1) alreadyOptimized(i)=2;

        delete solver;
    }
    solution=currentSolution;
}
#endif

/* Try solution ------------------------------------------------------------ */
double Prog_Angular_CommonLine::trySolution(const Matrix1D<double> &solution)
{
    Matrix1D<double> minAllowed(XSIZE(solution)), maxAllowed(XSIZE(solution));
    int idx=0;
    int Nimg=XSIZE(solution)/3;
    for (int i=0; i<Nimg; i++)
    {
        maxAllowed(idx++)=360;
        maxAllowed(idx++)=180;
        maxAllowed(idx++)=360;
    }

    Matrix1D<int> imgIdx(Nimg);
    imgIdx.initLinear(1,Nimg-1,1,"incr");
    bool bAtSolution;
    solver=new EulerSolver(XSIZE(solution),1,
        alreadyOptimized,currentSolution,imgIdx,this);
    solver->Setup(MULTIDIM_ARRAY(minAllowed), MULTIDIM_ARRAY(maxAllowed),
        stBest2Bin, 0.5, 0.8);
    solver->setShow(true);
    double energy=solver->EnergyFunction(MULTIDIM_ARRAY(solution),bAtSolution);
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
        
        DocFile DF;
        DF.append_comment("Rot Tilt Psi");
        
        int idx=0;
        SF.go_first_ACTIVE();
        while (!SF.eof())
        {
            ImageXmipp I;
            I.read(SF.NextImg());
            I.set_rot ((float)(currentSolution(idx++)));
            I.set_tilt((float)(currentSolution(idx++)));
            I.set_psi ((float)(currentSolution(idx++)));
            I.write();
            
            Matrix1D<double> assignment(5);
            assignment(0)=I.rot();
            assignment(1)=I.tilt();
            assignment(2)=I.psi();
            DF.append_data_line(assignment);
        }
        DF.write(fnOut);
    }
}
