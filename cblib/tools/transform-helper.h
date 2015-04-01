// Copyright (c) 2012 by Zuse-Institute Berlin and the Technical University of Denmark.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. Neither the name of the copyright holders nor contributors may not
//        be used to endorse or promote products derived from this software
//        without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS NOR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef CBF_TRANSFORM_HELPER_H
#define CBF_TRANSFORM_HELPER_H

#include "cbf-helper.h"
#include "cbf-format.h"
#include "transform.h"

#include <Eigen/SparseCore>
#include <Eigen/src/Core/DenseBase.h>

#include <stdio.h>
#include <set>

//
// DEFINITIONS
//

typedef void* CBFtransformmemory;

typedef CBFresponsee (*transform_varcone_handler)(CBFtransform_param param,
                                                  CBFtransformmemory* mem,
                                                  long long int& k,
                                                  long long int& varbeg,
                                                  CBFdata* data,
                                                  CBFdyndata& newdata);

typedef CBFresponsee (*transform_mapcone_handler)(CBFtransform_param param,
                                                  CBFtransformmemory* mem,
                                                  long long int& k,
                                                  long long int& rowbeg,
                                                  long long int& bbeg,
                                                  long long int& abeg,
                                                  double* mapmaxviol,
                                                  CBFdata* data,
                                                  CBFdyndata& newdata);

//
// Row comparison tools
//
typedef Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator rowiter;

int row_lexorder(const double b1, const double b2, rowiter itA1, rowiter itA2, const double scr1, const double scr2);

int is_quadcone_member_and_nonzero(rowiter it);

struct CmpRows
{
  const Eigen::SparseMatrix<double, Eigen::RowMajor>* b;
  const Eigen::SparseMatrix<double, Eigen::RowMajor>* A;
  const Eigen::VectorXd* rscal;

  CmpRows(const Eigen::SparseMatrix<double, Eigen::RowMajor>* b_,
          const Eigen::SparseMatrix<double, Eigen::RowMajor>* A_,
          const Eigen::VectorXd* rscal_)
      : b(b_)
      , A(A_)
      , rscal(rscal_)
  {
  }

  bool operator()(long long int r1, long long int r2)
  {
    return row_lexorder(
               b->coeff(r1, 0), b->coeff(r2, 0), rowiter(*A, r1), rowiter(*A, r2), (*rscal)(r1), (*rscal)(r2)) >= 0;
  }
};

typedef struct CBFmem_colmajor_struct
{

  // Refers to columns of the combined [b, A, vech(F)]
  Eigen::VectorXi colnnz;
  Eigen::VectorXi colmap;
  Eigen::VectorXi invcolmap;

  CBFmem_colmajor_struct(CBFdata* data)
      : colnnz(1 + data->varnum)
      , colmap(1 + data->varnum)
      , invcolmap(1 + data->varnum)
  {

    colnnz.setZero();
    colmap.setZero();
    invcolmap.setZero();
  }

} CBFmem_colmajor;

//
// FUNCTIONS
//

CBFresponsee transform_stackwise(CBFdata* data,
                                 CBFtransform_param& param,
                                 CBFtransformmemory* mem,
                                 transform_varcone_handler varcone_handler,
                                 transform_mapcone_handler mapcone_handler);

CBFresponsee
get_obj_entry_A(const CBFdata* data, Eigen::SparseMatrix<double, Eigen::RowMajor>& vec, long long int& objannz);

CBFresponsee get_cone_entry_A(const CBFdata* data,
                              long long int nrow,
                              long long int rowbeg,
                              long long int abeg,
                              Eigen::SparseMatrix<double, Eigen::RowMajor>& mat,
                              long long int& annz,
                              std::map<long long int, long long int>* compressedcolumnmap);

CBFresponsee get_cone_entry_b(const CBFdata* data,
                              long long int nrow,
                              long long int rowbeg,
                              long long int bbeg,
                              Eigen::SparseMatrix<double, Eigen::RowMajor>& vec,
                              long long int& bnnz);

CBFresponsee get_cone_entry_all(const CBFdata* data,
                                long long int nrow,
                                long long int rowbeg,
                                long long int bbeg,
                                long long int abeg,
                                Eigen::SparseMatrix<double, Eigen::RowMajor>& mat,
                                long long int& bnnz,
                                long long int& annz,
                                std::map<long long int, long long int>* compressedcolumnmap);

void strip_cone_entry_fixvarnnz(Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
                                Eigen::SparseMatrix<double, Eigen::RowMajor>& b,
                                std::map<long long int, long long int>* compressedcolumnmap,
                                const double* lb,
                                const double* ub,
                                const double eps);

CBFresponsee
put_obj_entry_A(const CBFdata* data, CBFdyndata* newdata, Eigen::SparseMatrix<double, Eigen::RowMajor>& vec);

CBFresponsee put_cone_entry_A(CBFdata* data,
                              CBFdyndata* newdata,
                              long long int rowbeg,
                              long long int a_beg,
                              long long int a_end,
                              const Eigen::SparseMatrix<double, Eigen::RowMajor>& mat,
                              std::map<long long int, long long int>* compressedcolumnmap);

CBFresponsee put_cone_entry_b(CBFdata* data,
                              CBFdyndata* newdata,
                              long long int rowbeg,
                              long long int b_beg,
                              long long int b_end,
                              const Eigen::SparseMatrix<double, Eigen::RowMajor>& vec);

CBFresponsee put_cone_entry_all(CBFdata* data,
                                CBFdyndata* newdata,
                                long long int rowbeg,
                                long long int b_beg,
                                long long int b_end,
                                long long int a_beg,
                                long long int a_end,
                                const Eigen::SparseMatrix<double, Eigen::RowMajor>& mat,
                                std::map<long long int, long long int>* compressedcolumnmap);

CBFresponsee put_varbound(CBFdata* data,
                          Eigen::VectorXd& lb,
                          Eigen::VectorXd& ub,
                          std::vector<bool>* stronglb,
                          std::vector<bool>* strongub);

CBFresponsee get_cone_matrix(const CBFdata* data,
                             long long int nrow,
                             long long int rowbeg,
                             long long int bbeg,
                             long long int abeg,
                             Eigen::SparseMatrix<double, Eigen::RowMajor>& mat,
                             long long int& bnnz,
                             long long int& annz);

CBFresponsee put_cone_matrix(CBFdata* data,
                             CBFdyndata* newdata,
                             long long int rowbeg,
                             long long int b_beg,
                             long long int b_end,
                             long long int a_beg,
                             long long int a_end,
                             const Eigen::SparseMatrix<double, Eigen::RowMajor>& mat);

CBFresponsee CBFintvars_incone_fetch(CBFdata* data, char* integerarray, std::set<long long int>& intvars);

int row_lexorder(long long int a, long long int b, const Eigen::SparseMatrix<double, Eigen::RowMajor>* mat);

struct CmpPairs
{
  const Eigen::SparseMatrix<double, Eigen::RowMajor>* data_;
  CmpPairs(const Eigen::SparseMatrix<double, Eigen::RowMajor>* data)
      : data_(data)
  {
  }
  bool operator()(long long int a, long long int b)
  {
    return row_lexorder(a, b, data_) >= 0;
  }
};

#endif
