/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <stdexcept>
#include <iostream>

#include <votca/xtp/eigen.h>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <votca/xtp/bse_operator.h>
#include <votca/xtp/davidsonsolver.h>



using boost::format;
using std::flush;

namespace votca { 
    namespace xtp {

using namespace std;

DavidsonSolver::DavidsonSolver(ctp::Logger &log) : _log(log) { }


Eigen::ArrayXd DavidsonSolver::_sort_index(Eigen::VectorXd &V) const
{
    Eigen::ArrayXd idx = Eigen::ArrayXd::LinSpaced(V.rows(),0,V.rows()-1);
    std::sort(idx.data(),idx.data()+idx.size(),
              [&](int i1, int i2){return V[i1]<V[i2];});
    return idx; 
}

Eigen::MatrixXd DavidsonSolver::_get_initial_eigenvectors(Eigen::VectorXd &d, int size_initial_guess) const
{
    Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(d.size(),size_initial_guess);
    Eigen::ArrayXd idx = DavidsonSolver::_sort_index(d);

    for (int j=0; j<size_initial_guess;j++) {
        guess(idx(j),j) = 1.0; 
    }

    return guess;
}

Eigen::MatrixXd DavidsonSolver::_solve_linear_system(Eigen::MatrixXd &A, Eigen::VectorXd &r) const
{
    Eigen::MatrixXd w;

    if (this->jacobi_linsolve == "CG") {
        Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower|Eigen::Upper> cg;
        cg.compute(A);
        w = cg.solve(r);
    }   

    else if (this->jacobi_linsolve == "GMRES") {
        Eigen::GMRES<Eigen::MatrixXd, Eigen::IdentityPreconditioner> gmres;
        gmres.compute(A);
        w = gmres.solve(r);
    }

    else if (this->jacobi_linsolve == "LLT") {
        w = A.llt().solve(r);
    }

    // default 
    else {
        throw std::runtime_error( " Linear Solver " + this->jacobi_linsolve + " not recognized");
    }
    
    return w;
}

template <class MatrixReplacement>
Eigen::MatrixXd DavidsonSolver::_jacobi_orthogonal_correction(MatrixReplacement &A, Eigen::VectorXd &r, Eigen::VectorXd &u, double lambda) const
{
    //! Solve the linear system of the jacobi correction
    /*!
        \param A Matrix of the system
        \param r residue vector
        \param u eigenvetor
        \param lambda eigenvalue
    */

    // form the projector P = 1 - u * u.T
    Eigen::MatrixXd P = -u*u.transpose();
    P.diagonal().array() += 1.0;

    // project the matrix P * (A - lambda*I) * P^T
    Eigen::MatrixXd projA = A*P.transpose();
    projA -= lambda*P.transpose();
    projA = P * projA;

    return DavidsonSolver::_solve_linear_system(projA,r);
}

template <class MatrixReplacement>
void DavidsonSolver::solve(MatrixReplacement &A, int neigen, int size_initial_guess)
{
        
    if(this->jacobi_correction)
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Jacobi-Davidson - linsolve : " << this->jacobi_linsolve << flush;
    else
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Davidson (DPR)" << flush;

    double norm;
    int size = A.rows();

    //. search space exeeding the system size
    if (max_search_space > size) {
        CTP_LOG(ctp::logDEBUG, _log)
                << ctp::TimeStamp() << " Warning Max search space (" << max_search_space << ") larger than system size (" << size << ")" << flush;   
        max_search_space = size;
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Max search space set to " << size << flush;                            
    }
    
    // initial guess size
    if (size_initial_guess == 0) {
        size_initial_guess = 2*neigen;
    }
    int search_space = size_initial_guess;

    // initialize the guess eigenvector
    Eigen::VectorXd Adiag = A.diagonal();    
    Eigen::MatrixXd V = DavidsonSolver::_get_initial_eigenvectors(Adiag, size_initial_guess);
    std::sort(Adiag.data(),Adiag.data()+Adiag.size());

    Eigen::MatrixXd thinQ; // thin matrix for QR 
    Eigen::VectorXd lambda; // eigenvalues hodlers
    
    // temp varialbes 
    Eigen::MatrixXd T, U, q;
    Eigen::VectorXd w, tmp;

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " iter\tSearch Space\tNorm" << flush;    
    for (int iiter = 0; iiter < iter_max; iiter++ ) {
        
        // orthogonalise the vectors
        if (iiter>0) {
            thinQ = Eigen::MatrixXd::Identity(V.rows(),V.cols());
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(V);
            V = qr.householderQ() * thinQ;
        }

        // project the matrix on the trial subspace and diagonalize
        T = A * V;
        T = V.transpose()*T;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
        lambda = es.eigenvalues();
        U = es.eigenvectors();
        
        // Ritz eigenvectors
        q = V.block(0,0,V.rows(),search_space)*U;
        norm = 0.0;

        // correction vectors
        for (int j=0; j<size_initial_guess; j++) {   

            // residue vector
            w = A*q.col(j) - lambda(j)*q.col(j);
            norm += w.norm();

            // jacobi-davidson correction
            if (this->jacobi_correction) {
                tmp = q.col(j);
                w = DavidsonSolver::_jacobi_orthogonal_correction<MatrixReplacement>(A,w,tmp,lambda(j));
            } else  {
                w = w / ( lambda(j) - Adiag(j) );  
            }

            // append the correction vector to the search space
            V.conservativeResize(Eigen::NoChange,V.cols()+1);
            V.col(V.cols()-1) = w;
        }

        // normalized the norm of the residue vector
        norm /= search_space;
        CTP_LOG(ctp::logDEBUG, _log)  << ctp::TimeStamp() 
           << format(" %1$4d \t %2$12d \t %3$4.2e/%4$.0e") % iiter % search_space % norm % tol << flush;         
        
        // break if converged, update otherwise
        if (norm < tol) {
            break;
        }
        else {
            search_space += size_initial_guess;
        }

        // restart
        if (search_space > max_search_space) {
            V = q.block(0,0,V.rows(),size_initial_guess);
            search_space = size_initial_guess;
        }
    }

    // store the eigenvalues/eigenvectors
    this->_eigenvalues = lambda.head(neigen);
    this->_eigenvectors = q.block(0,0,q.rows(),neigen);

    // normalize the eigenvectors
    for (int i=0; i<neigen; i++){
        this->_eigenvectors.col(i).normalize();
    }
}

template void DavidsonSolver::solve<Eigen::MatrixXd>(Eigen::MatrixXd &A, int neigen, int size_initial_guess=0);
template void DavidsonSolver::solve<BSE_OPERATOR>(BSE_OPERATOR &A, int neigen, int size_initial_guess=0);

}}