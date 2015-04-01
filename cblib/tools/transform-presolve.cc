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

#include "transform-presolve.h"

#include "transform-varstack-to-mapstack.h"
#include "transform-rquad.h"
#include "transform-placeholders.h"
#include "transform-lindep-elimination-soc.h"
#include "transform-probing.h"
#include "transform-analytics.h"
#include "transform-helper.h"

#include "cbf-helper.h"
#include "cbf-format.h"

#include <Eigen/Core>

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

CBFtransform const transform_presolve = { "presolve", transform };

static VectorXd lb, ub;
static char* integerarray = NULL;
static std::vector<bool> stronglb, strongub;
static double* mapmaxviol = NULL;
static long long int mapmaxviolsize = 0;
static bool infeas = false;

extern clock_t timer_lindep_elimination_soc;
extern clock_t timer_quadcone_entrywisedecoupled;
extern clock_t timer_quadcone_inequalityfree;
extern clock_t timer_propagation_graph;

// -------------------------------------
// Function definitions
// -------------------------------------

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag)
{
  CBFresponsee res = CBF_RES_OK;
  std::set<long long int> intvars_in_cone;
  std::set<long long int>::iterator probevar;
  double* mapmaxviol = NULL;
  CBFtransform_param fast_param = param;
  fast_param.USE_VARBOUNDING_FAST_INEQFREE = false;
  fast_param.USE_VARBOUNDING_FAST_COMPONENTANALYSIS = 0;

  stronglb = std::vector<bool>(data->varnum, false);
  strongub = std::vector<bool>(data->varnum, false);

  timer_lindep_elimination_soc = 0;
  timer_quadcone_entrywisedecoupled = 0;
  timer_quadcone_inequalityfree = 0;
  timer_propagation_graph = 0;

  if(changeflag) {
    printf("changeflag not supported by this transformation");
    return CBF_RES_ERR;
  }

  if(res == CBF_RES_OK)
    res = init_varbound(data, lb, ub);

  if(res == CBF_RES_OK)
    res = CBFintegerarray_init(data, &integerarray);

  if(res == CBF_RES_OK) {
    res = transform_varstack_to_mapstack.transform(data, param, NULL);
  }

  if(res == CBF_RES_OK && param.USE_RQUAD_CONVERTION) {
    res = transform_rquad.transform(data, param, NULL);
  }

  if(res == CBF_RES_OK && param.USE_VARBOUNDING == CBFtransform_param::FAST) {
    res = fast_transform(data, fast_param, NULL);
  }

  if(res == CBF_RES_OK && param.USE_PLACEHOLD_ELIM && !infeas) {
    transform_placeholders_init(integerarray, &stronglb, &strongub);
    res = transform_placeholders.transform(data, param, NULL);
  }

  if(res == CBF_RES_OK && param.USE_LINDEP_ELIMINATION != CBFtransform_param::NONE && !infeas) {
    res = transform_lindep_elimination_soc.transform(data, param, NULL);
  }

  if(res == CBF_RES_OK && param.USE_VARBOUNDING == CBFtransform_param::FAST && !infeas) {
    res = update_varbound_fast(data, param, integerarray, lb, ub, &stronglb, &strongub, NULL, NULL, &infeas, NULL);
  }

  //  // Debug update_varbound_fast procedure?
  //  // No significant reductions should be possible for integer variables
  //  if (res == CBF_RES_OK && param.USE_VARBOUNDING == CBFtransform_param::FULL)
  //      res = update_varbound_full(data, integerarray, true, lb, ub, &stronglb, &strongub, &infeas);
  //
  //      if (res == CBF_RES_OK)
  //        res = update_varbound_fast(data, param, integerarray, lb, ub, &stronglb, &strongub, NULL, &infeas);
  //    }
  //  }

  if(res == CBF_RES_OK && param.USE_PROBE && !infeas) {
    transform_probing_init(integerarray, &lb, &ub, &stronglb, &strongub, &infeas);
    res = transform_probing.transform(data, param, NULL);
  }

  // Revert rquad to quad?
  if (res == CBF_RES_OK) {
    if (param.USE_RQUAD_CONVERTION) {
      transform_rquad.revert(data, param, NULL);
    }
  }

  // Print bound status?
  //    printf("----\n");
  //    for (long long int i = 0; i < data->varnum; ++i) {
  //      printf("%g <= x[%lli] <= %g\n", lb[i], i, ub[i]);
  //    }

  // Append relaxation-strengthening variable bounds
  if(res == CBF_RES_OK) {
    res = put_varbound(data, lb, ub, &stronglb, &strongub);
  }

  // Compress representation
  if(res == CBF_RES_OK) {
    res = CBF_compress_maps(data, NULL);
  }

  //  printf("Time: %g %g %g %g\n",
  //         timer_lindep_elimination_soc / (double)CLOCKS_PER_SEC,
  //         timer_quadcone_entrywisedecoupled / (double)CLOCKS_PER_SEC,
  //         timer_quadcone_inequalityfree / (double)CLOCKS_PER_SEC,
  //         timer_propagation_graph / (double)CLOCKS_PER_SEC);

  printf(";%g",
         (timer_lindep_elimination_soc + timer_quadcone_entrywisedecoupled + timer_quadcone_inequalityfree +
          timer_propagation_graph) /
             (double)CLOCKS_PER_SEC);

  if(data->annz == 0 || infeas) {
    printf(";SOLVED");
    if(infeas) {
      printf("(INFEASIBLE)");
    }
  } else {
    printf(";OPEN");
  }

  // Clean memory
  CBFintegerarray_free(&integerarray);
  CBFmapmaxviol_free(&mapmaxviol);

  return res;
}

CBFresponsee fast_transform(CBFdata* data, CBFtransform_param& param, bool* changeflag)
{
  CBFresponsee res = CBF_RES_OK;
  long long int fixvar = 0;

  if(mapmaxviolsize == 0) {
    res = CBFmapmaxviol_init(data->mapnum, &mapmaxviol);

  } else if(mapmaxviolsize < data->mapnum) {
    CBFmapmaxviol_free(&mapmaxviol);
    res = CBFmapmaxviol_init(data->mapnum, &mapmaxviol);
  }

  if(res == CBF_RES_OK) {
    res = update_varbound_fast(
        data, param, integerarray, lb, ub, &stronglb, &strongub, mapmaxviol, &fixvar, &infeas, NULL);
  }

  if(res == CBF_RES_OK && fixvar >= 1 && !infeas) {
    res = CBF_stripfixvarnnz_maps(data, lb.data(), ub.data(), 1e-9);
    fixvar = 0;
  }

  // Compress representation
  if(res == CBF_RES_OK) {
    res = CBF_compress_maps(data, mapmaxviol);
  }

  return res;
}
