/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/aomatrix.h>
#include <votca/xtp/bsecoupling.h>
#include <votca/tools/constants.h>
#include <boost/format.hpp>



namespace votca { namespace xtp {

using boost::format;
using namespace tools;

void BSECoupling::Initialize(Property& options){
    
    #if (GWBSE_DOUBLE)
        CTP_LOG(ctp::logDEBUG, *_pLog) <<  " Compiled with full double support" << flush;   
    #else
        CTP_LOG(ctp::logDEBUG, *_pLog) <<  " Compiled with float/double mixture (standard)" << flush;   
    #endif
    
    std::string key = Identify(); 
    _doSinglets=false;
    _doTriplets=false;
   _output_perturbation=false;
    
    
    _openmp_threads = 0;
    
    if ( options.exists( key + ".openmp") ) {
                 _openmp_threads = options.get(key + ".openmp").as<int> ();
            }
    
    string spintype   = options.get(key + ".spin").as<string> ();
        if(spintype=="all"){
            _doSinglets=true;
            _doTriplets=true;
        }
        else if(spintype=="triplet"){
            _doTriplets=true;
        }
        else if(spintype=="singlet"){
            _doSinglets=true;
        }
        else{
            throw std::runtime_error((boost::format("Choice % for type not known. Available singlet,triplet,all") % spintype).str());
        }
     
     
   if ( options.exists( key + ".algorithm") ) {
                string algorithm = options.get(key + ".algorithm").as<string> ();
                 if(algorithm=="perturbation"){
                    _output_perturbation=true;
                 }
   }

        _levA  = options.get(key + ".moleculeA.states").as<int> ();
        _levB  = options.get(key + ".moleculeB.states").as<int> ();
        _occA  = options.get(key + ".moleculeA.occLevels").as<int> ();
        _occB  = options.get(key + ".moleculeB.occLevels").as<int> ();
        _unoccA  = options.get(key + ".moleculeA.unoccLevels").as<int> ();
        _unoccB  = options.get(key + ".moleculeB.unoccLevels").as<int> ();
  
}

void BSECoupling::WriteToProperty(const Orbitals& orbitalsA, const Orbitals& orbitalsB, 
                        Property& triplet_summary, int stateA, int stateB, double JAB){
  Property &coupling_summary = triplet_summary.add("coupling", (format("%1$1.6e") % JAB).str()); 
  double energyA = orbitalsA.BSETripletEnergies()(stateA)*conv::hrt2ev;
  double energyB = orbitalsB.BSETripletEnergies()(stateB)*conv::hrt2ev;
  coupling_summary.setAttribute("excitonA", stateA);
  coupling_summary.setAttribute("excitonB", stateB);
  coupling_summary.setAttribute("energyA", (format("%1$1.6e") % energyA).str());
  coupling_summary.setAttribute("energyB", (format("%1$1.6e") % energyB).str());
  coupling_summary.setAttribute("pert", (format("%1$1.6e") % getTripletCouplingElement( stateA , stateB, 0)).str());
  coupling_summary.setAttribute("diag", (format("%1$1.6e") % getTripletCouplingElement( stateA , stateB, 1)).str());
}

void BSECoupling::Addoutput(Property &type_summary,Orbitals& orbitalsA, 
                               Orbitals& orbitalsB){
   
    string algorithm="full_diag";
    int methodindex=1;
    if (_output_perturbation){
        algorithm="perturbation";
        methodindex=0;
    }
    type_summary.setAttribute("algorithm",algorithm);
    if (_doSinglets){
        Property &singlet_summary = type_summary.add("singlets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB <_levB ; ++stateB ) {
               double JAB = getSingletCouplingElement( stateA , stateB, methodindex);
               WriteToProperty(orbitalsA, orbitalsB, singlet_summary, stateA, stateB, JAB);    
           } 
        }
    }
    
    
    if ( _doTriplets){
        Property &triplet_summary = type_summary.add("triplets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB < _levA ; ++stateB ) {
               double JAB = getTripletCouplingElement( stateA , stateB, methodindex);
               WriteToProperty(orbitalsA, orbitalsB, triplet_summary, stateA, stateB, JAB);           
           } 
        }
    }       
}


double BSECoupling::getSingletCouplingElement( int levelA, int levelB, int methodindex) {
    return JAB_singlet[methodindex]( levelA  , levelB +  _levA ) * votca::tools::conv::hrt2ev;
}



double BSECoupling::getTripletCouplingElement( int levelA, int levelB, int methodindex) {
    return JAB_triplet[methodindex]( levelA  , levelB + _levA ) * votca::tools::conv::hrt2ev;
}


/**
 * \brief evaluates electronic couplings  
 *   
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 * @return false if failed
 */
bool BSECoupling::CalculateCouplings(const Orbitals& orbitalsA, const Orbitals& orbitalsB, const Orbitals& orbitalsAB) {
       CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Calculating exciton couplings" << flush;
     // set the parallelization 
    #ifdef _OPENMP
    
    if ( _openmp_threads > 0 ) omp_set_num_threads(_openmp_threads);      
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << " Using "<< omp_get_max_threads()<<" threads" << flush;
    #endif
    
    
    
    orbitalsAB.setCoupledExcitonsA(_levA);
    orbitalsAB.setCoupledExcitonsB(_levB);
    //check to see if ordering of atoms agrees
    const std::vector<QMAtom*> atomsA=orbitalsA.QMAtoms();
    const std::vector<QMAtom*> atomsB=orbitalsB.QMAtoms();
    const std::vector<QMAtom*> atomsAB=orbitalsAB.QMAtoms();
    
    for (int i=0;i<atomsAB.size();i++){
        QMAtom* dimer=atomsAB[i];
        QMAtom* monomer=NULL;
        if (i<atomsA.size()){
            monomer=atomsA[i];
        }
        else if (i<atomsB.size()+atomsA.size() ){
            monomer=atomsB[i-atomsA.size()];
        }
        else{
            throw runtime_error((boost::format("Number of Atoms in dimer %3i and the two monomers A:%3i B:%3i does not agree") %atomsAB.size() %atomsA.size() %atomsB.size()).str());
        }
        
        if(monomer->getType() != dimer->getType()){
            throw runtime_error("\nERROR: Atom types do not agree in dimer and monomers\n");
        }
        if(tools::abs(monomer->getPos()-dimer->getPos())>0.001){
            CTP_LOG(ctp::logERROR,*_pLog) << "======WARNING=======\n Coordinates of monomers and dimer atoms do not agree, do you know what you are doing?\n " << flush;
            break;
        }
        
    }
    
    
    // constructing the direct product orbA x orbB
    int _basisA = orbitalsA.getBasisSetSize();
    int _basisB = orbitalsB.getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        CTP_LOG(ctp::logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }

    // number of levels stored in monomers
    int _levelsA = orbitalsA.getNumberOfLevels();
    int _levelsB = orbitalsB.getNumberOfLevels();
    
        
    // get exciton information of molecule A
    int _bseA_cmax        = orbitalsA.getBSEcmax();
    int _bseA_cmin        = orbitalsA.getBSEcmin();
    int _bseA_vmax        = orbitalsA.getBSEvmax();
    int _bseA_vmin        = orbitalsA.getBSEvmin();
    int _bseA_vtotal      = _bseA_vmax - _bseA_vmin +1 ;
    int _bseA_ctotal      = _bseA_cmax - _bseA_cmin +1 ;
    int _bseA_size        = _bseA_vtotal * _bseA_ctotal;
    int _bseA_singlet_exc = orbitalsA.BSESingletCoefficients().cols();
    int _bseA_triplet_exc = orbitalsA.BSETripletCoefficients().cols();

    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   molecule A has " << _bseA_singlet_exc << " singlet excitons with dimension " << _bseA_size << flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   molecule A has " << _bseA_triplet_exc << " triplet excitons with dimension " << _bseA_size << flush;
    
    // now, two storage assignment matrices for two-particle functions
    Eigen::MatrixXi _combA;
    _combA.resize(_bseA_size,2);
    int _cnt = 0;
    for ( int _v = 0; _v < _bseA_vtotal; _v++){
        for ( int _c = 0; _c < _bseA_ctotal; _c++){
            _combA(_cnt,0) = _v;
            _combA(_cnt,1) = _bseA_vtotal + _c;
            _cnt++;
        }
    }
    
    // get exciton information of molecule B
    int _bseB_cmax        = orbitalsB.getBSEcmax();
    int _bseB_cmin        = orbitalsB.getBSEcmin();
    int _bseB_vmax        = orbitalsB.getBSEvmax();
    int _bseB_vmin        = orbitalsB.getBSEvmin();
    int _bseB_vtotal      = _bseB_vmax - _bseB_vmin +1 ;
    int _bseB_ctotal      = _bseB_cmax - _bseB_cmin +1 ;
    int _bseB_size        = _bseB_vtotal * _bseB_ctotal;
    int _bseB_singlet_exc = orbitalsB.BSESingletCoefficients().cols();
    int _bseB_triplet_exc = orbitalsB.BSETripletCoefficients().cols();
    

    
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   molecule B has " << _bseB_singlet_exc << " singlet excitons with dimension " << _bseB_size << flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   molecule B has " << _bseB_triplet_exc << " triplet excitons with dimension " << _bseB_size << flush;
    
    // now, two storage assignment matrices for two-particle functions
    Eigen::MatrixXi _combB;
    _combB.resize(_bseB_size,2);
    _cnt = 0;
    for ( int _v = 0; _v < _bseB_vtotal; _v++){
        for ( int _c = 0; _c < _bseB_ctotal; _c++){
            _combB(_cnt,0) = _bseA_vtotal + _bseA_ctotal + _v;
            _combB(_cnt,1) = _bseA_vtotal + _bseA_ctotal + _bseB_vtotal + _c;
            _cnt++;
        }
    }
    
    if(_levA>_bseA_singlet_exc){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of excitons you want is greater than stored for molecule A. Setting to max number available" << flush; 
        _levA=_bseA_singlet_exc;
    }
    if(_levB>_bseB_singlet_exc){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of excitons you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levB=_bseB_singlet_exc;
    }
    
    
    if(_levA>_bseA_singlet_exc){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of Frenkel states you want is greater than stored for molecule A. Setting to max number available" << flush; 
        _levA=_bseA_singlet_exc;
    }
    if(_levB>_bseB_singlet_exc){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of Frenkel states you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levB=_bseB_singlet_exc;
    }
    
    if(_unoccA>_bseA_ctotal){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of occupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccA=_bseA_ctotal;
    }
    else if (_unoccA<0){
        _unoccA=_bseA_ctotal;
    }
    if(_unoccB>_bseB_ctotal){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of occupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccB=_bseB_ctotal;
    }
    else if (_unoccB<0){
        _unoccB=_bseB_ctotal;
    }
    
    if(_occA>_bseA_vtotal){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of unoccupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occA=_bseA_vtotal;
    }
    else if (_occA<0){
        _occA=_bseA_vtotal;
    }
    if(_occB>_bseB_vtotal){
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Number of unoccupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occB=_bseB_vtotal;
    }else if (_occB<0){
        _occB=_bseB_vtotal;
    }
    
    
   
    
    // get exciton information of pair AB
    int _bseAB_cmax = orbitalsAB.getBSEcmax();
    int _bseAB_cmin = orbitalsAB.getBSEcmin();
    int _bseAB_vmax = orbitalsAB.getBSEvmax();
    int _bseAB_vmin = orbitalsAB.getBSEvmin();
    int _bseAB_vtotal = _bseAB_vmax - _bseAB_vmin +1 ;
    int _bseAB_ctotal = _bseAB_cmax - _bseAB_cmin +1 ;
    int _bseAB_size   = _bseAB_vtotal * _bseAB_ctotal;
    // check if electron-hole interaction matrices are stored
    if ( ! orbitalsAB.hasEHinteraction_triplet() && _doTriplets){
        CTP_LOG(ctp::logERROR,*_pLog) << "BSE EH for triplets not stored " << flush;
        return false;
    }
    if ( ! orbitalsAB.hasEHinteraction_singlet() && _doSinglets){
        CTP_LOG(ctp::logERROR,*_pLog) << "BSE EH for singlets not stored " << flush;
        return false;
    }
    const MatrixXfd&    _eh_t = _orbitalsAB.eh_t(); 
    const MatrixXfd&    _eh_s = _orbitalsAB.eh_s(); 
    if(_doTriplets){
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   dimer AB has BSE EH interaction triplet with dimension " << _eh_t.rows() << " x " <<  _eh_t.cols() << flush;
    }
    if(_doSinglets){
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   dimer AB has BSE EH interaction singlet with dimension " << _eh_s.rows() << " x " <<  _eh_s.cols() << flush;
    }
    // now, two storage assignment matrices for two-particle functions
    Eigen::MatrixXi _combAB;
    _combAB.resize(_bseAB_size,2);
    _cnt = 0;
    for ( int _v = 0; _v < _bseAB_vtotal; _v++){
        for ( int _c = 0; _c < _bseAB_ctotal; _c++){
            //_combAB(_cnt,0) = _v;
            //_combAB(_cnt,1) = _bseAB_vtotal + _c;
            
            _combAB(_cnt,0) = _bseAB_vmin + _v;
            _combAB(_cnt,1) = _bseAB_vmin + _bseAB_vtotal + _c;
            
            _cnt++;
        }
    }
    
    

    
    // DFT levels of monomers can be reduced to those used in BSE
    _levelsA = _bseA_vtotal + _bseA_ctotal;
    _levelsB = _bseB_vtotal + _bseB_ctotal;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   levels used in BSE of molA: " << _bseA_vmin << " to " << _bseA_cmax << " total: " << _bseA_vtotal + _bseA_ctotal <<  flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   levels used in BSE of molB: " << _bseB_vmin << " to " << _bseB_cmax << " total: " << _bseB_vtotal + _bseB_ctotal <<  flush;
    
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        CTP_LOG(ctp::logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
    
    //       | Orbitals_A          0 |      | Overlap_A |     
    //       | 0          Orbitals_B |.T  X   | Overlap_B |  X  ( Orbitals_AB )
    
    Eigen::MatrixXd _psi_AxB =Eigen::MatrixXd::Zero( _levelsA + _levelsB, _basisA + _basisB  );
    
  
    
    // constructing merged orbitals
    _psi_AxB.block(0,0,_levelsA , _basisA) = orbitalsA.MOCoefficients().block(_bseA_vmin,0, _bseA_cmax+1-_bseA_vmin, _basisA );
    _psi_AxB.block(_levelsA, _basisA,_levelsB,_basisB) =orbitalsB.MOCoefficients().block(_bseB_vmin,0,_bseB_cmax+1-_bseA_vmin,_basisB); 
    
    // psi_AxB * S_AB * psi_AB
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   projecting monomer onto dimer orbitals" << flush; 
    
     Eigen::MatrixXd _overlapAB;
    if ( orbitalsAB.hasAOOverlap() ) {
            CTP_LOG(ctp::logDEBUG,*_pLog) << "Reading overlap matrix from orbitals" << flush; 
           _overlapAB= orbitalsAB.AOOverlap();
    }else{
        CTP_LOG(ctp::logDEBUG,*_pLog) << "Calculating overlap matrix for basisset: "<< orbitalsAB.getDFTbasis()<< flush; 
        BasisSet _dftbasisset;
        AOBasis _dftbasis;
        _dftbasisset.LoadBasisSet(orbitalsAB.getDFTbasis());

        _dftbasis.AOBasisFill(_dftbasisset, orbitalsAB.QMAtoms());
        AOOverlap _dftAOoverlap;
        _dftAOoverlap.Fill(_dftbasis);
        _overlapAB=_dftAOoverlap.Matrix();
    }
    
  
    Eigen::MatrixXd _psi_AxB_dimer_basis = _psi_AxB.transpose()*_overlapAB*orbitalsAB.MOCoefficients();  
    _overlapAB.resize(0,0);
    
    
    //cout<< "_psi_AxB_dimer"<<endl;
    int LevelsA = _levelsA;
    for (int i=0;i<_psi_AxB_dimer_basis.rows();i++){
        double mag=_psi_AxB_dimer_basis.row(i).squaredNorm();
        if (mag<0.95){
            int monomer = 0;
            int level = 0;
            if ( i < LevelsA ) {
                monomer = 1;
                level   = _bseA_vmin + i;
            } else {
                monomer = 2;
                level   = _bseB_vmin + i -_levelsA;
                
            }
            CTP_LOG(ctp::logERROR,*_pLog) << "\nERROR: " << i << " Projection of orbital " << level << " of monomer " << monomer << " on dimer is insufficient,mag="<<mag<<" maybe the orbital order is screwed up, otherwise increase dimer basis.\n"<<flush;
        }
    }
   
    
    //notation AB is CT states with A+B-, BA is the counterpart
    //Setting up CT-states:
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Setting up CT-states" << flush; 
    //Number of A+B- states
    int noAB=_occA*_unoccB;
    //Number of A-B+ states
    int noBA=_unoccA*_occB;
    
    
    Eigen::MatrixXi comb_CTAB=Eigen::MatrixXi::Zero(noAB,2);
    
    _cnt = 0;   
    // iterate A over occupied, B over unoccupied
    int v_start=_bseA_vtotal-_occA;
    for ( int _v = v_start; _v < _bseA_vtotal; _v++){
        for ( int _c = 0; _c <_unoccB; _c++){            
            comb_CTAB(_cnt,0) =_v;
            comb_CTAB(_cnt,1) = _bseA_vtotal+_bseA_ctotal+_bseB_vtotal + _c;
           
            _cnt++;
        }
    }
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  <<"  "<<noAB <<" CT states A+B- created" << flush;
 
    Eigen::MatrixXi  comb_CTBA=Eigen::MatrixXi::Zero(noBA,2);
    _cnt = 0;
    // iterate A over unoccupied, B over occupied
    v_start=_bseB_vtotal-_occB;
    for ( int _v = v_start; _v < _bseB_vtotal; _v++){
        for ( int _c = 0; _c <_unoccA; _c++){            
            comb_CTBA(_cnt,0) =_bseA_vtotal+_bseA_ctotal+_v;
            comb_CTBA(_cnt,1) = _bseA_vtotal+ _c;
            
            _cnt++;
        }
    }
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  <<"  "<<noBA <<" CT states B+A- created" << flush;
    
    
    
    // these 4 matrixes, matrix(i,j) contains the j-th dimer MO component of the i-th excitation
   
    ctAB.resize(noAB,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_CT = 0 ; _i_CT < noAB ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctAB(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTAB(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTAB(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
    
    ctBA.resize(noBA,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_CT = 0 ; _i_CT < noBA ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctBA(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTBA(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTBA(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
      
    // some more convenient storage
    
    
    _kap.resize(_bseA_size,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) );
            
        }
    }

    
    

    _kbp.resize(_bseB_size,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) );
        }
    }
    
  
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()   << "   construct projection of product functions " << flush; 
 
    _psi_AxB_dimer_basis.resize(0,0);
    _combAB.resize(0,0);
    _combA.resize(0,0);
    _combB.resize(0,0);
    // now the different spin types
            if (_doSinglets) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Evaluating singlets" << flush;
                Eigen::MatrixXd Hamiltonian_AB = _eh_s.cast<double>();
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Setup Hamiltonian" << flush;
                const Eigen::MatrixXd bseA_T = orbitalsA.BSESingletCoefficients().block(0,0,orbitalsA.BSESingletCoefficients().rows(),_levA).transpose().cast<double>();
                const Eigen::MatrixXd bseB_T =orbitalsB.BSESingletCoefficients().block(0,0,orbitalsB.BSESingletCoefficients().rows(),_levB).transpose().cast<double>();
                
                JAB_singlet = ProjectExcitons(bseA_T, bseB_T, Hamiltonian_AB);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   calculated singlet couplings " << flush;
            }



            if (_doTriplets) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Evaluating triplets" << flush;
                Eigen::MatrixXd Hamiltonian_AB = _eh_s.cast<double>();
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "  Converted Hamiltonian to double" << flush;         
                const Eigen::MatrixXd bseA_T = orbitalsA.BSETripletCoefficients().block(0,0,orbitalsA.BSETripletCoefficients().rows(),_levA).transpose().cast<double>();
                const Eigen::MatrixXd bseB_T =orbitalsB.BSETripletCoefficients().block(0,0,orbitalsB.BSETripletCoefficients().rows(),_levB).transpose().cast<double>();
                JAB_triplet = ProjectExcitons(bseA_T, bseB_T, Hamiltonian_AB);
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   calculated triplet couplings " << flush;
            }
    
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "  Done with exciton couplings" << flush;
    return true;   
};


std::vector< Eigen::MatrixXd > BSECoupling::ProjectExcitons(const Eigen::MatrixXd& bseA_T, const Eigen::MatrixXd& bseB_T, 
                                  Eigen::MatrixXd& H){
   
     // get projection of monomer excitons on dimer product functions
     Eigen::MatrixXd _proj_excA = bseA_T* _kap;
     Eigen::MatrixXd _proj_excB =  bseB_T* _kbp;
     
     _bse_exc=_levA+_levB
     int _ctAB=ctAB.rows();
     int _ctBA=ctBA.rows();
     _ct=_ctAB+_ctBA;
     int nobasisfunc=H.rows();
 
     Eigen::MatrixXd fe_states=Eigen::MatrixXd::Zero(_bse_exc,nobasisfunc);
     fe_states.block(0,0, _levA, nobasisfunc )=_proj_excA;
     fe_states.block( _levA,0,_levB,nobasisfunc )=_proj_excB;
      
     Eigen::MatrixXd ct_states=Eigen::MatrixXd::Zero(_ct,nobasisfunc);

      if (_ct > 0) {
        //orthogonalize ct-states with respect to the FE states. 
        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Orthogonalizing CT-states with respect to FE-states" << flush;

        if (_ctAB > 0) {
          ct_states.block(0, 0, _ctAB, nobasisfunc) = ctAB;
        }
        if (_ctBA > 0) {
          ct_states.block(_ctAB, 0, _ctBA, nobasisfunc) = ctBA;
        }

        //orthogonalize ct-states with respect to FE states
        Eigen::MatrixXd correction = ct_states * fe_states.transpose() * fe_states;
        ct_states -= correction;
        correction.resize(0, 0);
        //normalize

        for (int i = 0; i < _ct; i++) {
          double norm = ct_states.row(i).norm();
          if (norm < 0.95) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " WARNING: CT-state " << i << " norm is only" << norm << flush;
          }
          ct_states.row(i) /= norm;
        }
      }
      
     Eigen::MatrixXd projection(_bse_exc+_ct,nobasisfunc);
     CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << " merging projections into one vector  " << flush;
  projection.block(0 ,0, _bse_exc,nobasisfunc)=fe_states;
   
     if(_ct>0){
    projection.block(_bse_exc,0,_ct ,nobasisfunc )=ct_states;
     }
      CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Setting up coupling matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // matrix _J 
    //  E_A         J_AB        J_A_ABCT        J_A_BACT
    //  J_BA        E_B         J_B_ABCT        J_B_BACT
    //  J_ABCT_A    J_ABCT_B    E_ABCT          J_ABCT_BACT
    //  J_BACT_A   J_BACT_B    J_BACT_ABCT     E_BACT
     
     // this only works for hermitian/symmetric H so only in TDA
    
     Eigen::MatrixXd J_dimer=projection*H*projection.transpose();
     H.resize(0,0);
    
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Setting up overlap matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;    
    Eigen::MatrixXd S_dimer=projection*projection.transpose();
    
    projection.resize(0,0);
    if(tools::globals::verbose &&  _bse_exc+_ct<100){
         CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "_J_dimer[Ryd]"<<flush;
     
     CTP_LOG(ctp::logDEBUG, *_pLog) << J_dimer<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "_S_dimer"<<flush;
     
     CTP_LOG(ctp::logDEBUG, *_pLog) << S_dimer<<flush;
      CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    }
   
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_dimer);
    Eigen::MatrixXd Sm1=es.operatorInverseSqrt();
    J_dimer=Sm1*J_dimer*Sm1;
   
    
    if(tools::globals::verbose && _bse_exc+_ct<100){
         CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << "_J_ortho[Ryd]"<<flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << J_dimer<<flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << "_S-1/2"<<flush;
    CTP_LOG(ctp::logDEBUG, *_pLog) << Sm1<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    }
     CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << 
             "   Smallest value of dimer overlapmatrix is "<<es.eigenvalues()(0)<< flush;
     
    std::vector< Eigen::MatrixXd >J;
     
     CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "   Running Perturbation algorithm"<< flush;
    J.push_back( Perturbation(J_dimer));
    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp()  << "    Running Projection algorithm"<< flush;
    J.push_back( Fulldiag(J_dimer));
    
    
       if(tools::globals::verbose){
     CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "Jeff_pert[Hrt]"<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << J[0]<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "Jeff_diag[Hrt]"<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << J[1]<<flush;
     CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     }
      
     return J;
}

Eigen::MatrixXd BSECoupling::Perturbation(const Eigen::MatrixXd& J_dimer){
    
    Eigen::MatrixXd Jmatrix =Eigen::MatrixXd::Zero(_bse_exc, _bse_exc);
    bool diag_ct = true;
    Eigen::MatrixXd J_result=J_dimer;
    if (_ct > 0 && diag_ct) {
        Eigen::MatrixXd transformation = Eigen::MatrixXd::Identity(_bse_exc + _ct, _bse_exc + _ct);
        Eigen::MatrixXd Ct = J_dimer.block(_bse_exc,_bse_exc,_ct, _ct);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ct);
        transformation.block(_bse_exc, _bse_exc ,_ct,_ct) = es.eigenvectors();
        Ct.resize(0, 0);

        if (tools::globals::verbose) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << "FE state hamiltonian" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << J_dimer.block(0,0, _bse_exc,_bse_exc) << flush;
            if (_ct > 0) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "eigenvalues of CT states" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << es.eigenvalues() << flush;
            }
        }

     J_result = transformation.transpose()*J_dimer*transformation;
        if (tools::globals::verbose && _bse_exc + _ct < 100) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "_J_ortho[Hrt] CT-state diag" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << J_result << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
        }
    }
    for (int stateA = 0; stateA < _levA; stateA++) {
        double Ea = J_result(stateA, stateA);
        for (int stateB = 0; stateB < _levB; stateB++) {
            int stateBd = stateB + _levA;
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Calculating coupling between exciton A" 
                    << stateA + 1 << " and exciton B" << stateB + 1 << flush;
            double J = J_result(stateA, stateBd);

            double Eb = J_result(stateBd, stateBd);
            for (int k = _bse_exc; k < (_bse_exc + _ct); k++) {
                double Eab = J_result(k, k);
                if (std::abs(Eab - Ea) < 0.001) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "Energydifference between state A " 
                            << stateA + 1 << "and CT state " << k + 1 << " is " << Eab - Ea << "[Hrt]" << flush;
                }
                if (std::abs(Eab - Eb) < 0.001) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "Energydifference between state B "
                            << stateB + 1 << "and CT state " << k + 1 << " is " << Eab - Eb << "[Hrt]" << flush;

                }
                J += 0.5 * J_result(k, stateA) * J_result(k, stateBd)*(1 / (Ea - Eab) + 1 / (Eb - Eab)); // Have no clue why 0.5
            }
            Jmatrix(stateA, stateBd) = J;
            Jmatrix(stateBd, stateA) = J;


        }
    }           
    return Jmatrix;
}


Eigen::MatrixXd BSECoupling::Fulldiag(const Eigen::MatrixXd& J_dimer){
    Eigen::MatrixXd Jmat = Eigen::MatrixXd::Zero(_bse_exc, _bse_exc);
   
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J_dimer);
    if (tools::globals::verbose && _bse_exc + _ct < 10) {
        CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << "Eigenvectors of J" << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << es.eigenvectors() << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << "J_eigenvalues[Hrt]" << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << es.eigenvalues() << flush;
        CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
    }
    //Calculate projection on subspace for every pair of excitons separately
    for (int stateA = 0; stateA < _levA; stateA++) {
        for (int stateB = 0; stateB < _levB; stateB++) {
            int stateBd = stateB + _levA;
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Calculating coupling between exciton A" << stateA + 1 << " and exciton B" << stateB + 1 << flush;
            std::vector<int> index;
            std::vector<int> signvec;
            for (int i = 0; i < _bse_exc + _ct; i++) {
                if (i == int(stateA) || i == int(stateBd)) {
                    double close = 0.0;
                    int ind = 0;
                    int sign = 0;
                    //row
                    for (int j = 0; j < _bse_exc + _ct; j++) {
                        bool check = true;
                        // if index i is already in index
                        // should not happen but if one vector was similar to two others.
                        for (unsigned l = 0; l < index.size(); l++) {
                            if (j == index[l]) {
                                check = false;
                                break;
                            }
                        }
                        if (check && std::abs(es.eigenvalues()(i, j)) > close) {
                            ind = j;
                            close = std::abs(es.eigenvalues()(i, j));
                            if (es.eigenvalues()(i, j) >= 0) {
                                sign = 1;
                            } else {
                                sign = -1;
                            }
                        }
                    }
                    index.push_back(ind);
                    signvec.push_back(sign);
                }
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Order is: [Initial state n->nth eigenvalue]" << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "    A" << stateA + 1 << ":" << stateA + 1 << "->" << index[0] + 1 << " ";
            CTP_LOG(ctp::logDEBUG, *_pLog) << "    B" << stateB + 1 << ":" << stateBd + 1 << "->" << index[1] + 1 << " " << flush;

            //setting up transformation matrix Tmat and diagonal matrix Emat for the eigenvalues;
            Eigen::MatrixXd Emat = Eigen::MatrixXd::Zero(2, 2);
            Eigen::MatrixXd Tmat = Eigen::MatrixXd::Zero(2, 2);
            //find the eigenvectors which are most similar to the initial states
            //row 
            for (int i = 0; i < 2; i++) {
                int k = index[i];
                double sign = signvec[i];
                double normr = 1 / std::sqrt(es.eigenvectors()(stateA, k) * es.eigenvectors()(stateA, k) + es.eigenvectors()(stateBd, k) * es.eigenvectors()(stateBd, k));
                Tmat(0, i) = sign * es.eigenvectors()(stateA, k) * normr;
                Tmat(1, i) = sign * es.eigenvectors()(stateBd, k) * normr;
                Emat(i, i) = es.eigenvectors()(k);
            }

            if ((Tmat(1, 1) * Tmat(0, 0) - Tmat(1, 0) * Tmat(0, 1)) < 0) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << " Reduced state matrix is not in a right handed basis, multiplying second eigenvector by -1 " << flush;
                Tmat(0, 1) = -Tmat(0, 1);
                Tmat(1, 1) = -Tmat(1, 1);
            }

            if (tools::globals::verbose) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "_T" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << Tmat << flush;

            }

            Eigen::MatrixXd S_small = Tmat*Tmat.transpose();
            if (tools::globals::verbose) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "S_small" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << S_small << flush;

            }
            //orthogonalize that matrix
            
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ss(S_small);
            Eigen::MatrixXd sm1=ss.operatorInverseSqrt();
            Emat=sm1*Emat*sm1;
            
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "   Smallest value of dimer overlapmatrix is " << ss.eigenvalues()(0) << flush;
            if (tools::globals::verbose) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "S-1/2" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << sm1 << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "E_ortho" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << Emat << flush;
            }
            Tmat = Tmat*sm1;
           
            if (tools::globals::verbose) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "T_ortho" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << Tmat << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
            }


            Eigen::MatrixXd J_small = Tmat*Emat*Tmat.transpose();
            if (tools::globals::verbose) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "T_ortho*E_ortho*T_ortho^T" << flush;
                CTP_LOG(ctp::logDEBUG, *_pLog) << J_small << flush;
            }

            Jmat(stateA, stateBd) = J_small(0, 1);
            Jmat(stateBd, stateA) = J_small(1, 0);

        }
    }
       
    return Jmat;
}


    
}}
