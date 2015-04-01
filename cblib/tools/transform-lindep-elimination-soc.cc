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

#include "transform-lindep-elimination-soc.h"
#include "transform-helper.h"
#include "cbf-helper.h"
#include "cbf-format.h"

#include <Eigen/Dense>
#include <Eigen/SparseQR>
#include <algorithm>
#include <cmath>
#include <stddef.h>
#include <stdio.h>
#include <ctime>
#include <iostream>

using namespace Eigen;

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag);

static CBFresponsee transform_callback(CBFtransform_param param,
                                       CBFtransformmemory* mem,
                                       long long int& k,
                                       long long int& rowbeg,
                                       long long int& bbeg,
                                       long long int& abeg,
                                       double* mapmaxviol,
                                       CBFdata* data,
                                       CBFdyndata& newdata);

// -------------------------------------
// Global variable
// ------------------------------------

CBFtransform const transform_lindep_elimination_soc = { "lindep-elimination-soc", transform };

clock_t timer_lindep_elimination_soc = 0;

// -------------------------------------
// Function definitions
// -------------------------------------

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag)
{
  CBFresponsee res;

  if(changeflag) {
    printf("changeflag not supported by this transformation");
    return CBF_RES_ERR;
  }

  if(data->fnnz >= 1) {
    printf("The transform '%s' does not yet support semidefinite variables", transform_lindep_elimination_soc.name);
    return CBF_RES_ERR;
  }

  clock_t start = clock();
  res = transform_stackwise(data, param, NULL, NULL, transform_callback);
  timer_lindep_elimination_soc += (clock() - start);
  return res;
}

static CBFresponsee transform_callback(CBFtransform_param param,
                                       CBFtransformmemory* mem,
                                       long long int& k,
                                       long long int& rowbeg,
                                       long long int& bbeg,
                                       long long int& abeg,
                                       double* mapmaxviol,
                                       CBFdata* data,
                                       CBFdyndata& newdata)
{
  CBFresponsee res = CBF_RES_OK;
  SparseMatrix<double, RowMajor> b;
  SparseMatrix<double, RowMajor> A;
  SparseMatrix<double, RowMajor> entries;
  std::map<long long int, long long int> compressedcolumnmap;
  long long int bnnz, annz, klen, abeg_r1, bbeg_r1, lastrow;
  redtype_t redtype = none;
  bool anymatrixchange = false;

  if(data->mapstackdomain[k] == CBF_CONE_QUAD) {
    klen = data->mapstackdim[k];

    if(res == CBF_RES_OK)
      res = CBF_findforward_map(data, rowbeg, NULL, &abeg, &bbeg);

    if(param.USE_LINDEP_ELIMINATION == CBFtransform_param::FAST) {

      if(res == CBF_RES_OK)
        res = get_cone_entry_b(data, klen, rowbeg, bbeg, b, bnnz);

      if(res == CBF_RES_OK)
        res = get_cone_entry_A(data, klen, rowbeg, abeg, A, annz, NULL);

      if(res == CBF_RES_OK)
        res = lindep_elimination_simple(rowbeg, klen, b, A, mapmaxviol, redtype, &anymatrixchange);

      if(anymatrixchange && redtype != line) {
        if(res == CBF_RES_OK)
          res = put_cone_entry_b(data, &newdata, rowbeg, bbeg, bbeg + bnnz - 1, b);

        if(res == CBF_RES_OK)
          res = put_cone_entry_A(data, &newdata, rowbeg, abeg, abeg + annz - 1, A, NULL);
      }

    } else {

      if(res == CBF_RES_OK)
        res = get_cone_entry_all(data, klen, rowbeg, bbeg, abeg, entries, bnnz, annz, &compressedcolumnmap);

      if(res == CBF_RES_OK)
        res = lindep_elimination_qr(rowbeg,
                                    klen,
                                    entries,
                                    param.USE_LINDEP_ELIMINATION_FULL_NNZCONSERVATIVE,
                                    mapmaxviol,
                                    redtype,
                                    &anymatrixchange);

      if(anymatrixchange && redtype != line) {
        if(res == CBF_RES_OK)
          res = put_cone_entry_all(
              data, &newdata, rowbeg, bbeg, bbeg + bnnz - 1, abeg, abeg + annz - 1, entries, &compressedcolumnmap);
      }
    }

    if(res == CBF_RES_OK) {
      switch(redtype) {
      case point:
//        printf("FACIAL REDUCTION TO POINT\n");
        data->mapstackdomain[k] = CBF_CONE_ZERO;
        break;

      case line:
//        printf("FACIAL REDUCTION TO LINE\n");
        data->mapstackdomain[k] = CBF_CONE_ZERO;
        mapmaxviol[rowbeg] = std::min(mapmaxviol[rowbeg], 0.0);

        // Change subi of first row to new L+ cone
        // (Assuming forward sweep to move (abeg,bbeg) past nonsorted coordinates)
        abeg_r1 = abeg, bbeg_r1 = bbeg;
        res = CBF_findforward_map(data, rowbeg + 1, NULL, &abeg, &bbeg);

        if(res == CBF_RES_OK)
          res = CBFdyn_map_capacitysurplus(&newdata, 1);

        if(res == CBF_RES_OK) {
          res = CBFdyn_map_adddomain(&newdata, CBF_CONE_POS, 1);
          lastrow = data->mapnum + newdata.data->mapnum - 1;
        }

        if(param.USE_LINDEP_ELIMINATION == CBFtransform_param::FAST) {

          if(res == CBF_RES_OK)
            res = put_cone_entry_b(data, &newdata, lastrow, bbeg_r1, bbeg - 1, b.topRows(1));

          if(res == CBF_RES_OK)
            res = put_cone_entry_A(data, &newdata, lastrow, abeg_r1, abeg - 1, A.topRows(1), NULL);

          if(res == CBF_RES_OK && anymatrixchange) {
            if(res == CBF_RES_OK)
              res = put_cone_entry_b(data, &newdata, rowbeg + 1, bbeg, bbeg + bnnz - 1, b.bottomRows(klen - 1));

            if(res == CBF_RES_OK)
              res = put_cone_entry_A(data, &newdata, rowbeg + 1, abeg, abeg + annz - 1, A.bottomRows(klen - 1), NULL);
          }

        } else {

          if(res == CBF_RES_OK)
            res = put_cone_entry_all(data,
                                     &newdata,
                                     lastrow,
                                     bbeg_r1,
                                     bbeg - 1,
                                     abeg_r1,
                                     abeg - 1,
                                     entries.topRows(1),
                                     &compressedcolumnmap);

          if(res == CBF_RES_OK && anymatrixchange)
            res = put_cone_entry_all(data,
                                     &newdata,
                                     rowbeg + 1,
                                     bbeg,
                                     bbeg + bnnz - 1,
                                     abeg,
                                     abeg + annz - 1,
                                     entries.bottomRows(klen - 1),
                                     &compressedcolumnmap);
        }

      default:
        break;
      }
    }
  }

  return res;
}

CBFresponsee lindep_elimination_simple(long long int rowbeg,
                                       long long int klen,
                                       SparseMatrix<double, RowMajor>& b,
                                       SparseMatrix<double, RowMajor>& A,
                                       double* mapmaxviol,
                                       redtype_t& redtype,
                                       bool* anymatrixchange)
{
  DiagonalMatrix<double, Dynamic> W(klen);
  VectorXd rmap(klen), rscal(klen);
  long long int r, ru, rd, emptyrows = 0;
  double lambda;
  bool isreduced = false;

  // Identity scaling and row mapping
  W.setIdentity();
  for(r = 0; r < klen; ++r) {
    rmap(r) = r;
  }

  // Normalize rows by first entry
  for(r = 0; r < klen; ++r) {
    if(b.coeff(r, 0) != 0.0) {
      rscal(r) = b.coeff(r, 0);

    } else {
      if(A.outerIndexPtr()[r] != A.outerIndexPtr()[r + 1]) {
        // First row-nnz of A
        rscal(r) = *(A.valuePtr() + A.outerIndexPtr()[r]);

      } else {
        // Empty row
        mapmaxviol[rowbeg + r] = std::min(mapmaxviol[rowbeg + r], 0.0);
        W.diagonal()[r] = 0.0;
        isreduced = true;

        if(r == 0) {
          redtype = point;
          return CBF_RES_OK;

        } else {
          // Swap to bottom and ignore
          std::swap(rmap(r), rmap(klen - emptyrows - 1));
          ++emptyrows;
        }
      }
    }
  }

  // Sort scaled rows and sweep
  std::sort(rmap.data(), rmap.data() + klen - emptyrows, CmpRows(&b, &A, &rscal));

  for(r = 0; r <= klen - emptyrows - 2; ++r) {

    // Check dublicate rows
    if(row_lexorder(b.coeff(rmap[r], 0),
                    b.coeff(rmap[r + 1], 0),
                    rowiter(A, rmap[r]),
                    rowiter(A, rmap[r + 1]),
                    rscal(rmap[r]),
                    rscal(rmap[r + 1])) == 0) {

      // Merge rows upwards
      if(rmap[r] < rmap[r + 1]) {
        rd = rmap[r + 1];
        ru = rmap[r + 1] = rmap[r];
      } else {
        rd = rmap[r];
        ru = rmap[r] = rmap[r + 1];
      }

      // Move W-value from rd to ru
      lambda = rscal(rd) / rscal(ru);
      if(ru == 0) {
        W.diagonal()[ru] -= W.diagonal()[rd] * lambda * lambda;
      } else {
        W.diagonal()[ru] += W.diagonal()[rd] * lambda * lambda;
      }
      mapmaxviol[rowbeg + rd] = std::min(mapmaxviol[rowbeg + rd], 0.0);
      W.diagonal()[rd] = 0.0;

      isreduced = true;
    }
  }

  // Interpret reduction
  if(W.diagonal()[0] < -1e-8) {
    redtype = point;

  } else if(W.diagonal()[0] < 1e-8) {
    redtype = line;

  } else if(isreduced) {
    if(anymatrixchange) {
      *anymatrixchange = true;
    }

    for(r = 0; r < klen; ++r) {
      W.diagonal()[r] = std::sqrt(W.diagonal()[r]);
    }

    b = W * b;
    A = W * A;
  }

  return CBF_RES_OK;
}

CBFresponsee lindep_elimination_qr(long long int rowbeg,
                                   long long int klen,
                                   SparseMatrix<double, RowMajor>& entries,
                                   bool nnz_conservative,
                                   double* mapmaxviol,
                                   redtype_t& redtype,
                                   bool* anymatrixchange)
{
  SparseQR<SparseMatrix<double, RowMajor>, COLAMDOrdering<int> > QR;
  SparseMatrix<double, RowMajor> redhyperball;
  SparseMatrix<double, RowMajor> lambda;
  SparseMatrix<double, RowMajor> W(klen, klen);
  VectorXd r1T;
  double norm;
  long long int j;

  entries.makeCompressed();
  QR.compute(entries.bottomRows(klen - 1));

  // Hyperball reduction
  // ||A'|| = ||Q R PT|| = ||R PT||
  redhyperball = QR.matrixR() * QR.colsPermutation().transpose();

  if(redhyperball.nonZeros() < entries.bottomRows(klen - 1).nonZeros() ||
     (!nnz_conservative && QR.rank() <= klen - 2)) {
    *anymatrixchange = true;

    entries.bottomRows(klen - 1) = redhyperball;

    for(j = 1; j <= (klen - 1) - QR.rank(); ++j)
      mapmaxviol[rowbeg + klen - j] = std::min(mapmaxviol[rowbeg + klen - j], 0.0);
  }

  // Dependency between the hyperball and its radius
  // xT A' = r1
  // (QTx)T R = r1 P
  // R^T (QTx) = (r1 P)^T

  //
  //  The easy slow way:
  //  QR.compute(entries.bottomRows(klen - 1).transpose());
  //  lambda = QR.solve(entries.row(0).transpose().eval()).topRows(klen-1).transpose();
  //

  r1T = (entries.topRows(1).eval() * QR.colsPermutation()).transpose().topRows(QR.rank());

  if(*anymatrixchange) {
    lambda = (QR.matrixR().topLeftCorner(QR.rank(), QR.rank()).transpose().triangularView<Lower>().solve(r1T))
                 .sparseView()
                 .transpose();
    lambda.conservativeResize(1, klen - 1);

  } else {
    lambda = (QR.matrixQ() *
              (MatrixXd::Identity(klen - 1, QR.rank()) *
               (QR.matrixR().topLeftCorner(QR.rank(), QR.rank()).transpose().triangularView<Lower>().solve(r1T))))
                 .sparseView()
                 .transpose();
  }

  if(row_lexorder(0.0, 0.0, rowiter(lambda * entries.bottomRows(klen - 1), 0), rowiter(entries, 0), 1.0, 1.0) == 0) {

    // Interpret reduction
    norm = lambda.norm();
    if(norm < 1 - 1e-4) {
      redtype = point;

      // Delete dependency in radius
      mapmaxviol[rowbeg] = std::min(mapmaxviol[rowbeg], 0.0);

    } else {

      if(norm < 1 + 1e-4) {
        redtype = line;
      }

      if(lambda.nonZeros() == 1) {

        if(redtype == none) {
          // Update radius
          *anymatrixchange = true;
          entries.topRows(1) = entries.topRows(1) * std::sqrt(1.0 - 1.0 / (norm * norm));
        }

        // Delete dependency in hyperball
        mapmaxviol[rowbeg + 1 + rowiter(lambda, 0).col()] =
            std::min(mapmaxviol[rowbeg + 1 + rowiter(lambda, 0).col()], 0.0);

      } else {

        // Perform Householder transformation on hyperball
        lambda = lambda * (1.0 / norm);
        lambda.coeffRef(0, 0) -= 1;
        lambda = lambda * (1.0 / lambda.norm());

        lambda.conservativeResize(1, klen);
        for(j = 0; j < lambda.innerSize(); ++j)
          ++lambda.innerIndexPtr()[j];

        W.setIdentity();
        if(redtype == none) {
          // Update radius
          W.coeffRef(0, 0) = std::sqrt(1.0 - 1.0 / (norm * norm));
        }
        W.middleRows(1, klen - 2) -=
            (SparseMatrix<double, RowMajor>)((2 * lambda.middleCols(1, klen - 2).transpose()) * lambda);
        W.coeffRef(klen - 1, klen - 1) = 0;

        // Test if we want to eliminate the radius-hyperball dependency
        redhyperball = W * entries;
        if(redhyperball.nonZeros() < entries.nonZeros() || (redtype == none) || (!nnz_conservative)) {

          *anymatrixchange = true;
          entries = W * entries;

          // Delete empty hyperball entry
          mapmaxviol[rowbeg + klen - 1] = std::min(mapmaxviol[rowbeg + klen - 1], 0.0);
        }
      }
    }
  }

  return CBF_RES_OK;
}
