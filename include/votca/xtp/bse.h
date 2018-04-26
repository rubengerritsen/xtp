/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_BSE_H
#define _VOTCA_XTP_BSE_H
#include <votca/xtp/linalg.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/ppm.h>
#include <votca/xtp/threecenter.h>



namespace votca {
namespace xtp {




class BSE {
 private:
     
struct Interaction {
    Eigen::VectorXd exchange_contrib;
    Eigen::VectorXd direct_contrib;
    Eigen::VectorXd qp_contrib;
};

struct Population {
    std::vector<Eigen::VectorXd> popH;
    std::vector<Eigen::VectorXd> popE;
    std::vector<Eigen::VectorXd> Crgs;
    Eigen::VectorXd popGs;
};   
    
 public:
 
    
  BSE(Orbitals* orbitals,ctp::Logger *log,double min_print_weight):
        _log(log),
        _orbitals(orbitals),
        _eh_s(orbitals->eh_s()),
        _eh_t(orbitals->eh_t()),
        _bse_singlet_energies(orbitals->BSESingletEnergies()),
        _bse_singlet_coefficients(orbitals->BSESingletCoefficients()),
        _bse_singlet_coefficients_AR(orbitals->BSESingletCoefficientsAR()),
        _bse_triplet_energies(orbitals->BSETripletEnergies()),
        _bse_triplet_coefficients(orbitals->BSETripletCoefficients()),
        _min_print_weight(min_print_weight){};

  ~BSE(){};
  
  void setGWData(const TCMatrix_gwbse* Mmn,const PPM* ppm,const Eigen::MatrixXd* Hqp){
      _Mmn=Mmn;
      _ppm=ppm;
      _Hqp=Hqp;
      
  }
  
  void setBSEindices(unsigned homo,int vmin, int vmax, int cmin, int cmax, int nmax) {
                _homo=homo;
                _bse_vmin = vmin;
                _bse_vmax = vmax;
                _bse_cmin = cmin;
                _bse_cmax = cmax;
                _bse_nmax = nmax;
                _bse_vtotal = _bse_vmax - _bse_vmin + 1;
                _bse_ctotal = _bse_cmax - _bse_cmin + 1;
                _bse_size = _bse_vtotal * _bse_ctotal;
                // indexing info BSE vector index to occupied/virtual orbital
                for (unsigned _v = 0; _v < _bse_vtotal; _v++) {
                    for (unsigned _c = 0; _c < _bse_ctotal; _c++) {
                        _index2v.push_back(_bse_vmin + _v);
                        _index2c.push_back(_bse_cmin + _c);
                    }
                }
                return;
            }

   
  void Solve_triplets();
  void Solve_singlets();
  void Solve_singlets_BTDA();
  void Analyze_triplets(const AOBasis& dftbasis);
  void Analyze_singlets(const AOBasis& dftbasis);
   
  void FreeMatrices(){
      _eh_t.resize(0, 0);
      _eh_s.resize(0, 0);
  }
  
  void SetupHs();
  
  void SetupHt();
  
  void FreeTriplets(){
      _bse_triplet_coefficients.resize(0,0);
  }
  
  void FreeSinglets(){
      _bse_singlet_coefficients.resize(0,0);
      _bse_singlet_coefficients_AR.resize(0,0);
  }
 

 private:
 
      
ctp::Logger *_log;
  unsigned _homo;
  unsigned  _bse_vmin;
  unsigned  _bse_vmax;
  unsigned  _bse_cmin;
  unsigned  _bse_cmax;
  unsigned  _bse_size;
  unsigned  _bse_vtotal;
  unsigned  _bse_ctotal;
  int _bse_nmax;
  
  Orbitals* _orbitals;
  
  const TCMatrix_gwbse* _Mmn;
  const PPM* _ppm;
  const Eigen::MatrixXd* _Hqp;
  

  // BSE variables and functions
  MatrixXfd& _eh_s;  // only for storage in orbitals object
  MatrixXfd& _eh_t;  // only for storage in orbitals object

  VectorXfd& _bse_singlet_energies;  // stored in orbitals object
  MatrixXfd& _bse_singlet_coefficients;  // stored in orbitals
                                                      // object
  MatrixXfd& _bse_singlet_coefficients_AR;  // stored in orbitals
                                                         // object
  VectorXfd& _bse_triplet_energies;  // stored in orbitals object
  MatrixXfd& _bse_triplet_coefficients;  // stored in orbitals
                                                      // object
  
  double _min_print_weight;

   template <typename T>
  void Add_Hqp(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H);
   template <typename T>
  void Add_Hx(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H, double factor);
   template <typename T>
   void Add_Hd(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H);
   template <typename T>
  void Add_Hd2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H, double factor);
  
  
  
  std::vector<int> _index2v;
  std::vector<int> _index2c;

 void printFragInfo(Population& pop, int i);
 void printWeights(unsigned i_bse, double weight);
 
  
  Interaction Analyze_eh_interaction(const std::string& spin);
  Eigen::VectorXd Analyze_IndividualContribution(const std::string& spin, const MatrixXfd& H);


  Population FragmentPopulations(const std::string& spin, const AOBasis& dftbasis);

  std::vector<Eigen::MatrixXd > CalcFreeTransition_Dipoles(const AOBasis& dftbasis);

  std::vector<tools::vec > CalcCoupledTransition_Dipoles(const AOBasis& dftbasis);
};
}
}

#endif /* _VOTCA_XTP_BSE_H */