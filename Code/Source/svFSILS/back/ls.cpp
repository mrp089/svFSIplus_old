
#include "fsils_api.hpp"
#include "fils_struct.hpp"

namespace fsi_linear_solver {

//-----------------
// fsils_ls_create
//-----------------
// Set solver parameters.
//
// Reproduces 'SUBROUTINE FSILS_LS_CREATE(ls, LS_type, relTol, absTol, maxItr, dimKry, relTolIn, absTolIn, maxItrIn)'.
//
void fsils_ls_create(FSILS_lsType& ls, LinearSolverType LS_type, double relTol, double absTol, int maxItr, int dimKry, 
  std::array<double,2> relTolIn, std::array<double,2> absTolIn, std::array<int,2> maxItrIn)
{
  std::cout << "[fsils_ls_create] " << std::endl;
  std::cout << "[fsils_ls_create] =========== fsils_ls_create ==========" << std::endl;

  using namespace consts;
  ls.foC = true;
  ls.LS_type = LS_type;

  // Set default parameters for each solver type.
  //
  switch (LS_type) {

    case LinearSolverType::LS_TYPE_NS:
      std::cout << "[fsils_ls_create] LS_type: LS_TYPE_NS" << std::endl;
      ls.RI.relTol = 0.4;
      ls.GM.relTol = 1.E-2;
      ls.CG.relTol = 0.2;
      ls.RI.mItr = 10;
      ls.GM.mItr = 2;
      ls.CG.mItr = 500;
      ls.GM.sD   = 100;
      ls.RI.sD   = 100;
    break;

    case LinearSolverType::LS_TYPE_GMRES:
      std::cout << "[fsils_ls_create] LS_type: LS_TYPE_GMRES " << std::endl;
      ls.RI.relTol = 0.1;
      ls.RI.mItr   = 4;
      ls.RI.sD     = 250;
    break;

    case LinearSolverType::LS_TYPE_CG:
      std::cout << "[fsils_ls_create] LS_type: LS_TYPE_CG " << std::endl;
      ls.RI.relTol = 1.E-2;
      ls.RI.mItr   = 1000;
    break;

    case LinearSolverType::LS_TYPE_BICGS:
      std::cout << "[fsils_ls_create] LS_type: LS_TYPE_BICGS " << std::endl;
      ls.RI.relTol = 1.E-2;
      ls.RI.mItr   = 500;
    break;

    default:
    break;
  }

  ls.RI.absTol = 1.E-10;
  ls.GM.absTol = 1.E-10;
  ls.CG.absTol = 1.E-10;

  // Set values optionally passed in.
  //
  if (present(relTol)) {
    std::cout << "[fsils_ls_create] Set ls.RI.relTol " << std::endl;
    ls.RI.relTol = relTol;
  }

  if (present(absTol)) {
    std::cout << "[fsils_ls_create] Set ls.RI.absTol " << std::endl;
    ls.RI.absTol = absTol;
  }

  if (present(maxItr)) {
    ls.RI.mItr = maxItr;
  }

  if (present(dimKry)) {
    ls.RI.sD = dimKry;
    ls.GM.sD = dimKry;
    std::cout << "[fsils_ls_create] Set ls.RI.sD: " << ls.RI.sD << std::endl;
    std::cout << "[fsils_ls_create] Set ls.GM.sD: " << ls.GM.sD << std::endl;
  }

  if (relTolIn.empty() != 0) {
    ls.GM.relTol = relTolIn[0];
    ls.CG.relTol = relTolIn[1];
  }

  if (absTolIn.empty() != 0) {
    ls.GM.absTol = absTolIn[0];
    ls.CG.absTol = absTolIn[1];
    std::cout << "[fsils_ls_create] Set ls.GM.absTol: " << ls.GM.absTol << std::endl;
  }

  if (maxItrIn.empty() != 0) {
    ls.GM.mItr = maxItrIn[0];
    ls.CG.mItr = maxItrIn[1];
  }

  std::cout << "[fsils_ls_create] Set ls.RI.absTol: " << ls.RI.absTol << std::endl;
  std::cout << "[fsils_ls_create] Set ls.GM.absTol: " << ls.GM.absTol << std::endl;
}

};

