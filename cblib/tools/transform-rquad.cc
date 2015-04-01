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

#include "transform-rquad.h"
#include "transform-helper.h"
#include "cbf-format.h"
#include "time.h"

#include <Eigen/SparseCore>

#include "stdio.h"

using namespace Eigen;

static CBFresponsee
transform(CBFdata *data, CBFtransform_param &param, bool *changeflag);

static CBFresponsee
transform_callback(CBFtransform_param param, CBFtransformmemory* mem, long long int &k, long long int &rowbeg, long long int &bbeg, long long int &abeg, double* mapmaxviol, CBFdata *data,
    CBFdyndata& newdata);

static CBFresponsee
revert(CBFdata *data, CBFtransform_param &param, bool *changeflag);

static CBFresponsee
revert_callback(CBFtransform_param param, CBFtransformmemory* mem, long long int &k, long long int &rowbeg, long long int &bbeg, long long int &abeg, double* mapmaxviol, CBFdata *data,
    CBFdyndata& newdata);

static void
construct_rotmat();

// -------------------------------------
// Global variables
// -------------------------------------

CBFtransform const transform_rquad = { "rquad", transform, revert };

static SparseMatrix<double, RowMajor> rotmat;

// -------------------------------------
// Function definitions
// -------------------------------------

static CBFresponsee transform(CBFdata *data, CBFtransform_param &param, bool *changeflag) {
  CBFresponsee res = CBF_RES_OK;

  if (changeflag) {
    printf("changeflag not supported by this transformation");
    return CBF_RES_ERR;
  }

  // Prepare constant global variable
  construct_rotmat();

  // Transform cones from rquad to quad in mapstack
  if (res == CBF_RES_OK)
    res = transform_stackwise(data, param, NULL, NULL, transform_callback);

  return res;
}

static CBFresponsee transform_callback(CBFtransform_param param, CBFtransformmemory* mem, long long int &k, long long int &rowbeg, long long int &bbeg, long long int &abeg,
    double* mapmaxviol, CBFdata *data, CBFdyndata& newdata) {

  CBFresponsee res = CBF_RES_OK;
  SparseMatrix<double, RowMajor> b;
  SparseMatrix<double, RowMajor> A;
  long long int bnnz, annz;

  if (data->mapstackdomain[k] == CBF_CONE_RQUAD) {
    
    if (res == CBF_RES_OK)
      res = CBF_findforward_map(data, rowbeg, NULL, &abeg, &bbeg);

    if (res == CBF_RES_OK)
      res = get_cone_entry_b(data, 2, rowbeg, bbeg, b, bnnz);
    
    if (res == CBF_RES_OK)
      res = put_cone_entry_b(data, &newdata, rowbeg, bbeg, bbeg + bnnz - 1, (rotmat*b).pruned(1e-16));
    
    if (res == CBF_RES_OK)
      res = get_cone_entry_A(data, 2, rowbeg, abeg, A, annz, NULL);
      
    if (res == CBF_RES_OK)
    {
      if (annz > 10) {
        res = put_cone_entry_A(data, &newdata, rowbeg, abeg, abeg + annz - 1, (rotmat*A).pruned(1e-16), NULL);
        
      } else {
        // Using insertion sort (fast for tiny nnz)
        SparseMatrix<double, RowMajor> rotmatA;
        rotmatA.resize(2, data->varnum);
        rotmatA.reserve(VectorXi::Constant(2, annz));
        double val;
        
        for (SparseMatrix<double, RowMajor>::InnerIterator it(A, 0); it && res == CBF_RES_OK; ++it) {
          val = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105350*it.value();
          rotmatA.insert(0, it.col()) = val;
          rotmatA.insert(1, it.col()) = val;
        }
        
        for (SparseMatrix<double, RowMajor>::InnerIterator it(A, 1); it && res == CBF_RES_OK; ++it) {
          val = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105350*it.value();
          rotmatA.coeffRef(0, it.col()) += val;
          rotmatA.coeffRef(1, it.col()) -= val;
        }
        
        rotmatA.prune(1e-16);
        res = put_cone_entry_A(data, &newdata, rowbeg, abeg, abeg + annz - 1, rotmatA, NULL);
      }
    }

    if (res == CBF_RES_OK)
      data->mapstackdomain[k] = CBF_CONE_QUAD;
  }

  return res;
}

static CBFresponsee revert(CBFdata *data, CBFtransform_param &param, bool *changeflag) {

  // Prepare constant global variable
  construct_rotmat();

  return transform_stackwise(data, param, NULL, NULL, revert_callback);
}

static CBFresponsee revert_callback(CBFtransform_param param, CBFtransformmemory* mem, long long int &k, long long int &rowbeg, long long int &bbeg, long long int &abeg, double* mapmaxviol,
    CBFdata *data, CBFdyndata& newdata) {

  CBFresponsee res = CBF_RES_OK;
  SparseMatrix<double, RowMajor> b;
  SparseMatrix<double, RowMajor> A;
  SparseMatrix<double, RowMajor> rotmatA, rotmatb;
  long long int bnnz, annz;

  if (data->mapstackdomain[k] == CBF_CONE_QUAD) {

    if (res == CBF_RES_OK)
      res = CBF_findforward_map(data, rowbeg, NULL, &abeg, &bbeg);

    if (res == CBF_RES_OK)
      res = get_cone_entry_b(data, 2, rowbeg, bbeg, b, bnnz);

    if (res == CBF_RES_OK)
      rotmatb = (rotmat*b).pruned(1e-16);

    if (res == CBF_RES_OK)
      res = get_cone_entry_A(data, 2, rowbeg, abeg, A, annz, NULL);

    if (res == CBF_RES_OK)
    {
      if (annz > 10) {
        rotmatA = rotmat*A;

      } else {
        // Using insertion sort (fast for tiny nnz)
        double val;
        rotmatA.resize(2, data->varnum);
        rotmatA.reserve(VectorXi::Constant(2, annz));

        for (SparseMatrix<double, RowMajor>::InnerIterator it(A, 0); it && res == CBF_RES_OK; ++it) {
          val = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105350*it.value();
          rotmatA.insert(0, it.col()) = val;
          rotmatA.insert(1, it.col()) = val;
        }

        for (SparseMatrix<double, RowMajor>::InnerIterator it(A, 1); it && res == CBF_RES_OK; ++it) {
          val = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105350*it.value();
          rotmatA.coeffRef(0, it.col()) += val;
          rotmatA.coeffRef(1, it.col()) -= val;
        }
      }

      rotmatA.prune(1e-16);
    }

    if (res == CBF_RES_OK && (rotmatA.nonZeros() < annz)) {
      data->mapstackdomain[k] = CBF_CONE_RQUAD;

      if (res == CBF_RES_OK)
        res = put_cone_entry_b(data, &newdata, rowbeg, bbeg, bbeg + bnnz - 1, rotmatb);

      if (res == CBF_RES_OK)
        res = put_cone_entry_A(data, &newdata, rowbeg, abeg, abeg + annz - 1, rotmatA, NULL);
    }
  }

  return res;
}

static void construct_rotmat() {
  rotmat.resize(2, 2);
  rotmat.reserve(VectorXi::Constant(2, 2));
  rotmat.insert(0, 0) = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105350;
  rotmat.insert(0, 1) = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105350;
  rotmat.insert(1, 0) = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105350;
  rotmat.insert(1, 1) = -0.70710678118654752440084436210484903928483593768847403658833986899536623923105350;
  rotmat.makeCompressed();
}
