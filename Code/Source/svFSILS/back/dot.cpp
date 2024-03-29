
// The dot product between two scaler-vector container vectors are
// calculated here.
//
// Only the part of U and V which are owned by this processor is
// included in dot product calculation
// In order to have the correct answer it is needed that COMMU has
// been done before calling this function (or the ansesters of U and
// V are passed through COMMU)

#include "dot.h"

#include "fils_struct.hpp"

namespace dot {

//--------------
// fsils_dot_s
//--------------
// Reproduces 'FUNCTION FSILS_DOTS(nNo, commu, U, V)'. 
//
double fsils_dot_s(const int nNo, FSILS_commuType& commu, const Vector<double>& U, const Vector<double>& V)
{
  #define n_debug_fsils_dot_s
  int tid = commu.task;
  auto msg_prefix = std::string("[debug_fsils_dot_s:") + std::to_string(tid) + "] ";
  #ifdef debug_fsils_dot_s 
  std::cout << msg_prefix << "========== fsils_dot_s ==========" << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  #endif

  double result = 0.0; 

  for (int i = 0; i < nNo; i++) {
    result = result + U(i)*V(i);
  }

  if (commu.nTasks == 1) {
    return result;
  }

  double tmp{0.0};
  MPI_Allreduce(&result, &tmp, 1, cm_mod::mpreal, MPI_SUM, commu.comm);
  // CALL MPI_ALLREDUCE(FSILS_DOTS, tmp, 1, mpreal, MPI_SUM, commu%comm, ierr)

  #ifdef debug_fsils_dot_s 
  std::cout << msg_prefix << "result: " << result << " tmp: " << tmp << std::endl;
  #endif

  return tmp;
}

//-------------
// fsils_dot_v
//-------------
// Reproduces 'FUNCTION FSILS_DOTV(dof, nNo, commu, U, V)'.
//
double fsils_dot_v(const int dof, const int nNo, FSILS_commuType& commu, const Array<double>& U, const Array<double>& V)
{
  double result = 0.0; 

  switch (dof) {
    case 1:
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*V(0,i);
      }
    break;

    case 2:
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i);
      }
    break;

    case 3:
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i) +  U(2,i)*V(2,i);
      }
    break;

    case 4:
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i) + U(2,i)*V(2,i) + U(3,i)*V(3,i);
      }
    break;

    defualt: 
      for (int i = 0; i < nNo; i++) {
        double sum{0.0};
        for (int j = 0; j < U.num_rows(); j++) {
          sum += U(j,i) * V(j,i);
        }
        result = result + sum; 
        //result = result + SUM(U(:,i)*V(:,i))
      }
  }

  if (commu.nTasks == 1) {
    return result;
  }

  double tmp{0.0};
  MPI_Allreduce(&result, &tmp, 1, cm_mod::mpreal, MPI_SUM, commu.comm);
  //CALL MPI_ALLREDUCE(result, tmp, 1, mpreal, MPI_SUM, commu%comm, ierr)

  return tmp;
}

//----------------
// fsils_nc_dot_s
//----------------
// Reproduces Fortran 'FSILS_NCDOTS(nNo, , U, V)'.
//
double fsils_nc_dot_s(const int nNo, const Vector<double>& U, const Vector<double>& V)
{
  double result{0.0};

  for (int i = 0; i < nNo; i++) { 
    result = result + U(i)*V(i);
  } 

  return result;
}

//----------------
// fsils_nc_dot_v
//----------------
// Reproduces 'FUNCTION FSILS_NCDOTV(dof, nNo, U, V) RESULT(FSILS_DOTV)'.
//
double fsils_nc_dot_v(const int dof, const int nNo, const Array<double>& U, const Array<double>& V)
{
  double result = 0.0;

  switch (dof) {
    case 1: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U(0,i)*V(0,i);
      }
    } break;

    case 2: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i);
      }
    } break;

    case 3: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i) + U(2,i)*V(2,i);
      }
    } break;

    case 4: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i) +  U(2,i)*V(2,i) + U(3,i)*V(3,i);
      }
    } break;

    default: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U.col(i) * V.col(i);
      }
    } break;
  } 

  return result;
}


};


