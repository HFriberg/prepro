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

#include "transform-helper.h"
#include "cbf-helper.h"

#include <Eigen/src/Core/DenseBase.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/Core/util/ForwardDeclarations.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>

#include <stdio.h>
#include <iostream>

#include "cbf-data.h"
#include "programmingstyle.h"

using namespace Eigen;

// -------------------------------------
// Function definitions
// -------------------------------------

CBFresponsee CBFintvars_incone_fetch(CBFdata *data, char* integerarray, std::set<long long int> &intvars)
{
  long long int k, klen;
  long long int row, abeg;

  row = 0;
  abeg = 0;
  for(k = 0; k < data->mapstacknum; ++k) {
    klen = data->mapstackdim[k];
    
    switch(data->mapstackdomain[k]) {
      case CBF_CONE_QUAD:
      case CBF_CONE_RQUAD:
        CBF_findforward_map(data, row, NULL, &abeg, NULL);
      
        for(; (abeg < data->annz) && (data->asubi[abeg] < row+klen); ++abeg) {
          if (integerarray[data->asubj[abeg]]) {
            intvars.insert(data->asubj[abeg]);
          }
        }
        break;
        
      default:
        break;
    }
    
    row += klen;
  }

  return CBF_RES_OK;
}

int row_lexorder(const double b1, const double b2, rowiter itA1, rowiter itA2, const double scr1, const double scr2)
{
  double v1;

  // Compare b
  if(b1 / scr1 > b2 / scr2 + 1e-9)
    return 1;
  else if(b1 / scr1 < b2 / scr2 - 1e-9)
    return -1;

  // Compare A
  for(; itA2; ++itA2) {
    while(itA1 && itA1.col() < itA2.col()) {
      if(itA1.value() / scr1 > 0.0 + 1e-9)
        return 1;
      else if(itA1.value() / scr1 < 0.0 - 1e-9)
        return -1;

      ++itA1;
    }

    if(itA1 && itA1.col() == itA2.col()) {
      v1 = itA1.value() / scr1;
      ++itA1;
    } else {
      v1 = 0.0;
    }

    if(v1 > itA2.value() / scr2 + 1e-9)
      return 1;
    else if(v1 < itA2.value() / scr2 - 1e-9) {
      return -1;
    }
  }

  while(itA1) {
    if(itA1.value() / scr1 > 0.0 + 1e-9)
      return 1;
    else if(itA1.value() / scr1 < 0.0 - 1e-9)
      return -1;

    ++itA1;
  }

  return 0;
}

int is_quadcone_member_and_nonzero(rowiter it)
{
  double radius;
  double hyperball_norm = 0;

  if(it.col() == 0 && it.value() != 0.0) {
    radius = it.value();
    ++it;
  } else {
    return 0;
  }

  for(; it; ++it) {
    hyperball_norm += it.value() * it.value();
  }

  if(hyperball_norm <= radius * radius) {
    return (radius >= 0 ? 1 : -1);
  } else {
    return 0;
  }
}

CBFresponsee transform_stackwise(CBFdata* data,
                                 CBFtransform_param& param,
                                 CBFtransformmemory *mem,
                                 transform_varcone_handler varcone_handler,
                                 transform_mapcone_handler mapcone_handler)
{
  CBFresponsee res = CBF_RES_OK;
  long long int j, k, klen, rowbeg, bbeg, abeg;
  double* mapmaxviol = NULL;
  CBFdata newdata = {
    0,
  };
  CBFdyndata dyn_newdata = {
    &newdata,
  };
  CBFdyndata dyndata;

  if(data->fnnz >= 1) {
    printf("The transform is not yet support semidefinite variables");
    return CBF_RES_ERR;
  }

  // Delete coordinates by giving them value zero.
  // Delete rows by setting 'mapmaxviol' to non-positive number indicating redundancy.
  res = CBFmapmaxviol_init(data->mapnum, &mapmaxviol);

  // Sort coordinates row major style
  if(res == CBF_RES_OK)
    res = CBF_coordinatesort(data->asubi, data->asubj, data->aval, data->annz, data->mapnum, data->varnum);
  if(res == CBF_RES_OK)
    res = CBF_coordinatesort(data->bsubi, data->bval, data->bnnz, data->mapnum);

  // Main loops
  if(varcone_handler) {
    rowbeg = 0;
    for(k = 0; k < data->varstacknum && res == CBF_RES_OK; ++k) {
      klen = data->varstackdim[k];
      res = varcone_handler(param, mem, k, rowbeg, data, dyn_newdata);

      if(klen != data->varstackdim[k]) {
        printf("varstackdim[%i] was reduced", k);
        return CBF_RES_ERR;
      }

      rowbeg += klen;
    }
  }

  if(mapcone_handler) {
    rowbeg = 0;
    bbeg = 0;
    abeg = 0;
    for(k = 0; k < data->mapstacknum && res == CBF_RES_OK; ++k) {
      klen = data->mapstackdim[k];

      res = mapcone_handler(param, mem, k, rowbeg, bbeg, abeg, mapmaxviol, data, dyn_newdata);

      if(klen != data->mapstackdim[k]) {
        printf("mapstackdim[%i] was reduced", k);
        return CBF_RES_ERR;
      }

      rowbeg += klen;
    }
  }

  // Prepare mapmaxviol
  if(newdata.mapnum >= 1 && res == CBF_RES_OK) {
    mapmaxviol = (double*)realloc(mapmaxviol, (data->mapnum + newdata.mapnum) * sizeof(double));
    if(!mapmaxviol)
      res = CBF_RES_ERR;
    else {
      for(j = data->mapnum; j < data->mapnum + newdata.mapnum; ++j) {
        mapmaxviol[j] = INFINITY;
      }
    }
  }

  // Write results
  if(res == CBF_RES_OK)
    res = CBFdyn_assign(&dyndata, data);

  if(res == CBF_RES_OK)
    res = CBFdyn_append(&dyndata, &newdata);

  if(res == CBF_RES_OK)
    res = CBF_compress_maps(data, mapmaxviol);

  // Clean memory
  CBFdyn_freedynamicallocations(&dyn_newdata);
  CBFmapmaxviol_free(&mapmaxviol);

  return res;
}

CBFresponsee get_obj_entry_A(const CBFdata* data, SparseMatrix<double, RowMajor>& vec, long long int& objannz)
{

  long long int idx;
  VectorXi wi(1);

  // Preallocate
  objannz = 0;

  for(idx = 0; idx < data->objannz; ++idx) {
    if(data->objaval[idx] != 0) {
      ++objannz;
    }
  }
  wi(0) = objannz;

  vec.resize(1, data->varnum);
  vec.reserve(wi);

  // Fill
  for(idx = 0; idx < data->objannz; ++idx)
    if(data->objaval[idx] != 0)
      vec.insertBackUncompressed(0, data->objasubj[idx]) = data->objaval[idx];

  return CBF_RES_OK;
}

CBFresponsee get_cone_entry_A(const CBFdata* data,
                              long long int nrow,
                              long long int rowbeg,
                              long long int abeg,
                              SparseMatrix<double, RowMajor>& mat,
                              long long int& annz,
                              std::map<long long int, long long int>* compressedcolumnmap)
{

  std::map<long long int, long long int> columncompression;
  long long int idx, nnzcols = 0;
  VectorXi wi(nrow);

  // Preallocate
  wi.setZero();
  annz = 0;

  for(idx = abeg; (idx < data->annz) && (data->asubi[idx] < rowbeg + nrow); ++idx) {
    if(data->aval[idx] != 0) {
      ++annz;
      ++wi(data->asubi[idx] - rowbeg);

      if(compressedcolumnmap &&
         columncompression.insert(std::pair<long long int, long long int>(data->asubj[idx], nnzcols)).second) {
        compressedcolumnmap->insert(std::pair<long long int, long long int>(nnzcols, data->asubj[idx]));
        ++nnzcols;
      }
    }
  }

  // Fill
  if(!compressedcolumnmap) {
    mat.resize(nrow, data->varnum);
    mat.reserve(wi);

    for(idx = abeg; (idx < data->annz) && (data->asubi[idx] < rowbeg + nrow); ++idx)
      if(data->aval[idx] != 0)
        mat.insertBackUncompressed(data->asubi[idx] - rowbeg, data->asubj[idx]) = data->aval[idx];

  } else {
    mat.resize(nrow, nnzcols);
    mat.reserve(wi);

    for(idx = abeg; (idx < data->annz) && (data->asubi[idx] < rowbeg + nrow); ++idx)
      if(data->aval[idx] != 0)
        mat.insertBackUncompressed(data->asubi[idx] - rowbeg, columncompression[data->asubj[idx]]) = data->aval[idx];
  }

  return CBF_RES_OK;
}

CBFresponsee get_cone_entry_b(const CBFdata* data,
                              long long int nrow,
                              long long int rowbeg,
                              long long int bbeg,
                              SparseMatrix<double, RowMajor>& vec,
                              long long int& bnnz)
{

  long long int idx;
  VectorXi wi(nrow);

  // Preallocate
  wi.setZero();
  bnnz = 0;

  for(idx = bbeg; (idx < data->bnnz) && (data->bsubi[idx] < rowbeg + nrow); ++idx) {
    if(data->bval[idx] != 0) {
      ++bnnz;
      ++wi(data->bsubi[idx] - rowbeg);
    }
  }

  vec.resize(nrow, 1);
  vec.reserve(wi);

  // Fill
  for(idx = bbeg; (idx < data->bnnz) && (data->bsubi[idx] < rowbeg + nrow); ++idx) {
    if(data->bval[idx] != 0) {
      vec.insertBackUncompressed(data->bsubi[idx] - rowbeg, 0) = data->bval[idx];
    }
  }

  return CBF_RES_OK;
}

CBFresponsee get_cone_entry_all(const CBFdata* data,
                                long long int nrow,
                                long long int rowbeg,
                                long long int bbeg,
                                long long int abeg,
                                SparseMatrix<double, RowMajor>& mat,
                                long long int& bnnz,
                                long long int& annz,
                                std::map<long long int, long long int>* compressedcolumnmap)
{

  std::map<long long int, long long int> columncompression;
  long long int idx, nnzcols = 0;
  VectorXi wi(nrow);

  // Preallocate
  wi.setZero();
  bnnz = 0, annz = 0;

  for(idx = bbeg; (idx < data->bnnz) && (data->bsubi[idx] < rowbeg + nrow); ++idx) {
    if(data->bval[idx] != 0) {
      ++bnnz;
      ++wi(data->bsubi[idx] - rowbeg);
    }
  }

  for(idx = abeg; (idx < data->annz) && (data->asubi[idx] < rowbeg + nrow); ++idx) {
    if(data->aval[idx] != 0) {
      ++annz;
      ++wi(data->asubi[idx] - rowbeg);

      if(compressedcolumnmap &&
         columncompression.insert(std::pair<long long int, long long int>(data->asubj[idx], nnzcols)).second) {
        compressedcolumnmap->insert(std::pair<long long int, long long int>(nnzcols, data->asubj[idx]));
        ++nnzcols;
      }
    }
  }

  // Fill
  if(!compressedcolumnmap) {
    mat.resize(nrow, 1 + data->varnum);
    mat.reserve(wi);

    for(idx = bbeg; (idx < data->bnnz) && (data->bsubi[idx] < rowbeg + nrow); ++idx)
      if(data->bval[idx] != 0)
        mat.insertBackUncompressed(data->bsubi[idx] - rowbeg, 0) = data->bval[idx];

    for(idx = abeg; (idx < data->annz) && (data->asubi[idx] < rowbeg + nrow); ++idx)
      if(data->aval[idx] != 0)
        mat.insertBackUncompressed(data->asubi[idx] - rowbeg, 1 + data->asubj[idx]) = data->aval[idx];

  } else {
    mat.resize(nrow, 1 + nnzcols);
    mat.reserve(wi);

    for(idx = bbeg; (idx < data->bnnz) && (data->bsubi[idx] < rowbeg + nrow); ++idx)
      if(data->bval[idx] != 0)
        mat.insertBackUncompressed(data->bsubi[idx] - rowbeg, 0) = data->bval[idx];

    for(idx = abeg; (idx < data->annz) && (data->asubi[idx] < rowbeg + nrow); ++idx)
      if(data->aval[idx] != 0)
        mat.insertBackUncompressed(data->asubi[idx] - rowbeg, 1 + columncompression[data->asubj[idx]]) =
            data->aval[idx];
  }

  return CBF_RES_OK;
}

void strip_cone_entry_fixvarnnz(Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
                                Eigen::SparseMatrix<double, Eigen::RowMajor>& b,
                                std::map<long long int, long long int>* compressedcolumnmap,
                                const double* lb,
                                const double* ub,
                                const double eps)
{
  long long int r;

  if (compressedcolumnmap) {
    long long int var;
    for(r = 0; r < A.outerSize(); ++r) {
      for(SparseMatrix<double, RowMajor>::InnerIterator it(A, r); it; ++it) {
        var = (*compressedcolumnmap)[it.col()];
        if(!isinf(lb[var]) && !isinf(ub[var]) && fabs(ub[var] - lb[var]) <= eps) {
          b.coeffRef(r, 0) += it.value() * lb[var];
          it.valueRef() = 0;
        }
      }
    }
  } else {
   for(r = 0; r < A.outerSize(); ++r) {
      for(SparseMatrix<double, RowMajor>::InnerIterator it(A, r); it; ++it) {
        if(!isinf(lb[it.col()]) && !isinf(ub[it.col()]) && fabs(ub[it.col()] - lb[it.col()]) <= eps) {
          b.coeffRef(r, 0) += it.value() * lb[it.col()];
          it.valueRef() = 0;
        }
      }
    } 
  }
  
  A.prune(1e-16);
  b.prune(1e-16);
}

CBFresponsee
put_obj_entry_A(const CBFdata* data, CBFdyndata* newdata, Eigen::SparseMatrix<double, Eigen::RowMajor>& vec)
{
  CBFresponsee res = CBF_RES_OK;
  long long int r, obja_beg = 0;

  // Preallocate
  if(vec.nonZeros() > data->objannz) {
    res = CBFdyn_obja_capacitysurplus(newdata, vec.nonZeros() - data->objannz);
  }

  // Fill
  for(r = 0; r < vec.outerSize() && res == CBF_RES_OK; ++r) {
    for(SparseMatrix<double, RowMajor>::InnerIterator it(vec, r); it && res == CBF_RES_OK; ++it) {
      if(obja_beg < data->objannz) {
        // Overwrite existing nnz
        data->objasubj[obja_beg] = it.col();
        data->objaval[obja_beg] = it.value();
        ++obja_beg;

      } else {
        // Append to end
        res = CBFdyn_obja_add(newdata, it.col(), it.value());
      }
    }
  }

  // Delete remaining old nnz
  if(res == CBF_RES_OK) {
    while(obja_beg < data->objannz)
      data->objaval[obja_beg++] = 0;
  }

  return res;
}

CBFresponsee put_cone_entry_A(CBFdata* data,
                              CBFdyndata* newdata,
                              long long int rowbeg,
                              long long int a_beg,
                              long long int a_end,
                              const SparseMatrix<double, RowMajor>& mat,
                              std::map<long long int, long long int>* compressedcolumnmap)
{
  CBFresponsee res = CBF_RES_OK;
  long long int r;

  // Preallocate
  if(mat.nonZeros() > a_end - a_beg + 1) {
    res = CBFdyn_a_capacitysurplus(newdata, mat.nonZeros() - (a_end - a_beg + 1));
  }

  // Fill
  if(!compressedcolumnmap) {
    for(r = 0; r < mat.outerSize() && res == CBF_RES_OK; ++r) {
      for(SparseMatrix<double, RowMajor>::InnerIterator it(mat, r); it && res == CBF_RES_OK; ++it) {
        if(a_beg <= a_end) {
          // Overwrite existing nnz
          data->asubi[a_beg] = rowbeg + it.row();
          data->asubj[a_beg] = it.col();
          data->aval[a_beg] = it.value();
          ++a_beg;

        } else {
          // Append to end
          res = CBFdyn_a_add(newdata, rowbeg + it.row(), it.col(), it.value());
        }
      }
    }

  } else {
    for(r = 0; r < mat.outerSize() && res == CBF_RES_OK; ++r) {
      for(SparseMatrix<double, RowMajor>::InnerIterator it(mat, r); it && res == CBF_RES_OK; ++it) {
        if(a_beg <= a_end) {
          // Overwrite existing nnz
          data->asubi[a_beg] = rowbeg + it.row();
          data->asubj[a_beg] = (*compressedcolumnmap)[it.col()];
          data->aval[a_beg] = it.value();
          ++a_beg;

        } else {
          // Append to end
          res = CBFdyn_a_add(newdata, rowbeg + it.row(), (*compressedcolumnmap)[it.col()], it.value());
        }
      }
    }
  }

  // Delete remaining old nnz
  if(res == CBF_RES_OK) {
    while(a_beg <= a_end)
      data->aval[a_beg++] = 0;
  }

  return res;
}

CBFresponsee put_cone_entry_b(CBFdata* data,
                              CBFdyndata* newdata,
                              long long int rowbeg,
                              long long int b_beg,
                              long long int b_end,
                              const Eigen::SparseMatrix<double, Eigen::RowMajor>& vec)
{
  CBFresponsee res = CBF_RES_OK;
  long long int r;

  // Preallocate
  if(vec.nonZeros() > b_end - b_beg + 1) {
    res = CBFdyn_b_capacitysurplus(newdata, vec.nonZeros() - (b_end - b_beg + 1));
  }

  // Fill
  for(r = 0; r < vec.outerSize() && res == CBF_RES_OK; ++r) {
    for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(vec, r); it && res == CBF_RES_OK; ++it) {
      if(b_beg <= b_end) {
        // Overwrite existing nnz
        data->bsubi[b_beg] = rowbeg + it.row();
        data->bval[b_beg] = it.value();
        ++b_beg;

      } else {
        // Append to end
        res = CBFdyn_b_add(newdata, rowbeg + it.row(), it.value());
      }
    }
  }

  // Overwrite remaining nnz
  if(res == CBF_RES_OK) {
    while(b_beg <= b_end)
      data->bval[b_beg++] = 0;
  }

  return res;
}

CBFresponsee put_cone_entry_all(CBFdata* data,
                                CBFdyndata* newdata,
                                long long int rowbeg,
                                long long int b_beg,
                                long long int b_end,
                                long long int a_beg,
                                long long int a_end,
                                const SparseMatrix<double, RowMajor>& mat,
                                std::map<long long int, long long int>* compressedcolumnmap)
{
  CBFresponsee res = CBF_RES_OK;

  res = put_cone_entry_b(data, newdata, rowbeg, b_beg, b_end, mat.middleCols(0, 1));

  if(res == CBF_RES_OK)
    res = put_cone_entry_A(data, newdata, rowbeg, a_beg, a_end, mat.middleCols(1, mat.cols() - 1), compressedcolumnmap);

  return res;
}

CBFresponsee
put_varbound(CBFdata* data, VectorXd& lb, VectorXd& ub, std::vector<bool>* stronglb, std::vector<bool>* strongub)
{

  CBFdyndata newdata;
  CBFresponsee res = CBF_RES_OK;
  long long int i, r;

  // Eliminate fixed variables
  res = CBF_stripfixvarnnz_maps(data, lb.data(), ub.data(), 1e-9);

  // Allow CBFdyndata to overwrite 'mapstack', 'a' and 'b' in data without CBF_RES_ERR
  newdata.data = data;
  newdata.mapstackdyncap = data->mapstacknum;
  newdata.adyncap = data->annz;
  newdata.bdyncap = data->bnnz;

  // Ensure sufficient space
  if(res == CBF_RES_OK)
    res = CBFdyn_a_capacitysurplus(&newdata, 2 * data->varnum);

  if(res == CBF_RES_OK)
    res = CBFdyn_b_capacitysurplus(&newdata, 2 * data->varnum);

  if(res == CBF_RES_OK)
    res = CBFdyn_map_capacitysurplus(&newdata, 2);

  // Write lower bounds
  r = data->mapnum;
  for(i = 0; i < data->varnum && res == CBF_RES_OK; ++i) {
    if(!isinf(lb[i]) && (!stronglb || (*stronglb)[i]) && !(!isinf(ub[i]) && fabs(ub[i] - lb[i]) <= 1e-9)) {
      res = CBFdyn_map_adddomain(&newdata, CBF_CONE_POS, 1);

      if(res == CBF_RES_OK)
        res = CBFdyn_a_add(&newdata, r, i, 1);

      if(res == CBF_RES_OK)
        res = CBFdyn_b_add(&newdata, r, -lb[i]);

      ++r;
    }
  }

  // Write upper bounds
  r = data->mapnum;
  for(i = 0; i < data->varnum && res == CBF_RES_OK; ++i) {
    if(!isinf(ub[i]) && (!strongub || (*strongub)[i]) && !(!isinf(lb[i]) && fabs(ub[i] - lb[i]) <= 1e-9)) {
      res = CBFdyn_map_adddomain(&newdata, CBF_CONE_NEG, 1);

      if(res == CBF_RES_OK)
        res = CBFdyn_a_add(&newdata, r, i, 1);

      if(res == CBF_RES_OK)
        res = CBFdyn_b_add(&newdata, r, -ub[i]);

      ++r;
    }
  }

  return res;
}

int row_lexorder(long long int a, long long int b, const SparseMatrix<double, RowMajor>* mat)
{
  long long int i;

  const long long int abeg = mat->outerIndexPtr()[a];
  const long long int alen = mat->outerIndexPtr()[a + 1] - abeg;
  const int* acol = mat->innerIndexPtr() + abeg;
  const double* aval = mat->valuePtr() + abeg;

  const long long int bbeg = mat->outerIndexPtr()[b];
  const long long int blen = mat->outerIndexPtr()[b + 1] - bbeg;
  const int* bcol = mat->innerIndexPtr() + bbeg;
  const double* bval = mat->valuePtr() + bbeg;

  for(i = 0; i < alen && i < blen; ++i) {

    if(acol[i] < bcol[i]) {
      if(aval[i] < 0) {
        return -1;
      } else if(aval[i] > 0) {
        return 1;
      } else {
        printf("ERROR: row_lexorder is not good");
      }

    } else if(acol[i] > bcol[i]) {
      if(bval[i] < 0) {
        return 1;
      } else if(bval[i] > 0) {
        return -1;
      } else {
        printf("ERROR: row_lexorder is not good");
      }

    } else {
      if(aval[i] > bval[i] + 1e-9)
        return 1;
      else if(aval[i] + 1e-9 < bval[i])
        return -1;
      else {
        // Continue to next column..
      }
    }
  }

  if(alen > blen)
    return 1;
  else if(alen < blen)
    return -1;
  else
    return 0;
}
