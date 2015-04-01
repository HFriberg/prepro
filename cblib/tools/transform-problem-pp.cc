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

#include "transform-problem-pp.h"
#include "transform-presolve.h"

#include "transform-varstack-to-mapstack.h"
#include "transform-rquad.h"
#include "transform-placeholders.h"
#include "transform-lindep-elimination-soc.h"
#include "transform-analytics.h"
#include "transform-helper.h"

#include "cbf-helper.h"
#include "cbf-format.h"

#include <Eigen/Core>

#include <math.h>
#include <vector>
#include <stddef.h>
#include <stdio.h>
#include <iostream>

#include <ctime>

using namespace Eigen;

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag);

// -------------------------------------
// Global variable
// -------------------------------------

CBFtransform const transform_problem_pp = { "pp-problem", transform };

// -------------------------------------
// Function definitions
// -------------------------------------

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag)
{
  CBFresponsee res = CBF_RES_OK;
  VectorXd lb, ub;
  std::vector<bool> stronglb(data->varnum, false), strongub(data->varnum, false);
  char* integerarray = NULL;
  std::set<long long int> intvars_in_cone;
  std::set<long long int>::iterator probevar;
  double* mapmaxviol = NULL;
  long long int fixvar = 0, nnz, i;
  bool infeas = false;
  CBFtransform_param fast_param = param;
  fast_param.USE_VARBOUNDING_FAST_INEQFREE = false;
  fast_param.USE_VARBOUNDING_FAST_COMPONENTANALYSIS = 0;

  long long int n;
  SparseMatrix<double, RowMajor> obj;
  SparseMatrix<double, RowMajor> knapsack;
  double contopt, intopt1, intopt2;

  if(changeflag) {
    printf("changeflag not supported by this transformation");
    return CBF_RES_ERR;
  }

  if(res == CBF_RES_OK)
    res = init_varbound(data, lb, ub);

  if(res == CBF_RES_OK)
    res = CBFintegerarray_init(data, &integerarray);

  if(res == CBF_RES_OK) {
    // Side-effect: coordinates are sorted
    res = transform_varstack_to_mapstack.transform(data, param, NULL);
  }

  if(res == CBF_RES_OK) {
    res = transform_rquad.transform(data, param, NULL);
  }

  if(res == CBF_RES_OK) {
    res = CBFmapmaxviol_init(data->mapnum, &mapmaxviol);

    if(res == CBF_RES_OK) {
      res = update_varbound_fast(
          data, fast_param, integerarray, lb, ub, &stronglb, &strongub, mapmaxviol, &fixvar, &infeas, NULL);
    }

    if(res == CBF_RES_OK && fixvar >= 1) {
      res = CBF_stripfixvarnnz_maps(data, lb.data(), ub.data(), 1e-9);
      fixvar = 0;
    }

    // Compress representation
    if(res == CBF_RES_OK) {
      res = CBF_compress_maps(data, mapmaxviol);
    }

    CBFmapmaxviol_free(&mapmaxviol);
  }

  if(res == CBF_RES_OK) {
    transform_placeholders_init(integerarray, &stronglb, &strongub);
    res = transform_placeholders.transform(data, param, NULL);
  }

  // Compress representation
  if(res == CBF_RES_OK) {
    res = CBF_compress_maps(data, NULL);
  }

  // MAJOR ASSUMPTIONS:
  // n = varstacknum - 2
  // x = varidx 0 to n-1
  // y = varidx n to 2n-1
  // knapsack at mapidx 0

  n = data->varstacknum - 2;

  if(res == CBF_RES_OK)
    get_obj_entry_A(data, obj, nnz);
    
  if(res == CBF_RES_OK)
    res = get_cone_entry_A(data, 1, 0, 0, knapsack, nnz, NULL);
  
  if(res == CBF_RES_OK) {
    // Attainability
    for (i=0; i<n; ++i) {
      if (obj.coeff(0, i) == 0.0 && knapsack.coeff(0, i) == 0.0) {
        // Fix x[i] to infinity and y[i] to 0
        printf("DETECTED: x -> INF\n");
        lb[i] = 0;
        ub[i] = 0;
        lb[n+i] = 0;
        ub[n+i] = 0;
        
      } else if (obj.coeff(0, n+i) == 0.0) {
        // Fix y[i] to infinity and x[i] to 0
        printf("DETECTED: y -> INF\n");
        lb[i] = 0;
        ub[i] = 0;
        lb[n+i] = 0;
        ub[n+i] = 0;
      }
    }
    
    // LB
    for (i=0; i<n; ++i) {
      lb[i] = std::max(1.0, lb[i]);
      stronglb[i] = true;
    }
    
    // UB
    for (i=0; i<n; ++i) {
      contopt = std::sqrt(obj.coeff(0, n+i) / obj.coeff(0, i));
      intopt1 = std::floor(contopt);
      intopt2 = std::ceil(contopt);
      
      if (obj.coeff(0, i)*intopt1 + obj.coeff(0, n+i)/intopt1 <= obj.coeff(0, i)*intopt2 + obj.coeff(0, n+i)/intopt2) {
        ub[i] = intopt1;
      } else {
        ub[i] = intopt2;
      }
      
      strongub[i] = true;
    }
  }
  
  // Revert rquad to quad?
  if (res == CBF_RES_OK) {
    if (param.USE_RQUAD_CONVERTION) {
      transform_rquad.revert(data, param, NULL);
    }
  }

  // Append relaxation-strengthening variable bounds
  if(res == CBF_RES_OK) {
    res = put_varbound(data, lb, ub, &stronglb, &strongub);
  }

  // Compress representation
  if(res == CBF_RES_OK) {
    res = CBF_compress_maps(data, NULL);
  }

  // Clean memory
  CBFintegerarray_free(&integerarray);

  return res;
}
