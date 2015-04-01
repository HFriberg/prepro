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

#include "transform-probing.h"
#include "transform-analytics.h"
#include "transform-presolve.h"
#include "transform-helper.h"
#include "cbf-helper.h"

#include <stdio.h>
#include <iostream>

using namespace Eigen;

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag);

static CBFresponsee transform_coefimprove(CBFtransform_param param,
                                          CBFtransformmemory* mem,
                                          long long int& k,
                                          long long int& rowbeg,
                                          long long int& bbeg,
                                          long long int& abeg,
                                          double* mapmaxviol,
                                          CBFdata* data,
                                          CBFdyndata& newdata);

typedef struct CBFtransformmemory_coefimprove_enum
{
  double* mapmaxviol;
  long long int probevar;
  double probeval;
  double fixpoint;

} CBFtransformmemory_coefimprovee;

// -------------------------------------
// Global variables
// -------------------------------------

CBFtransform const transform_probing = {
  "probing",
  transform,
};

static VectorXd* lb, *ub;
static std::vector<bool>* stronglb, *strongub;
static char* integerarray = NULL;
static bool* infeas;

// -------------------------------------
// Function definitions
// -------------------------------------

void transform_probing_init(char* integerarray__,
                            VectorXd* lb__,
                            VectorXd* ub__,
                            std::vector<bool>* stronglb__,
                            std::vector<bool>* strongub__,
                            bool* infeas__)
{
  integerarray = integerarray__;
  lb = lb__;
  ub = ub__;
  stronglb = stronglb__;
  strongub = strongub__;
  infeas = infeas__;
}

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag)
{
  CBFresponsee res = CBF_RES_OK;

  if(changeflag) {
    printf("changeflag not supported by this transformation");
    return CBF_RES_ERR;
  }

  CBFtransform_param fast_param, probing_param;
  std::set<long long int> probelist;
  std::set<long long int>::iterator probevar;
  CBFtransformmemory_coefimprovee coefimprove;
  double* mapmaxviol = NULL;
  VectorXd probelb, probeub;
  bool repeatedprobevar;
  bool condinfeas;
  bool condredund;
  long long int stepsize;
  long long int iter;
  long long int i;

  // Prepare parameters
  fast_param = param;
  fast_param.USE_VARBOUNDING_FAST_DECOUPLEDCONE = param.USE_PROBE_VARBOUNDING_FAST_DECOUPLEDCONE;
  fast_param.USE_VARBOUNDING_FAST_INEQFREE = false;
  fast_param.USE_LINDEP_ELIMINATION = CBFtransform_param::FAST;
  fast_param.USE_VARBOUNDING_FAST_COMPONENTANALYSIS = 0;

  probing_param = param;
  probing_param.USE_VARBOUNDING_FAST_DECOUPLEDCONE = param.USE_PROBE_VARBOUNDING_FAST_DECOUPLEDCONE;
  probing_param.USE_VARBOUNDING_FAST_INEQFREE = param.USE_PROBE_VARBOUNDING_FAST_INEQFREE;
  probing_param.USE_VARBOUNDING_FAST_COMPONENTANALYSIS = param.USE_PROBE_VARBOUNDING_FAST_COMPONENTANALYSIS;

  // Prepare redundancy indicator storage
  if(res == CBF_RES_OK)
    res = CBFmapmaxviol_init(data->mapnum, &mapmaxviol);

  // Select candidates for probing
  if(res == CBF_RES_OK)
    res = CBFintvars_incone_fetch(data, integerarray, probelist);

  probevar = probelist.begin();
  iter = 0;
  stepsize = 1;
  repeatedprobevar = false;
  while(probevar != probelist.end() && !(*infeas) && iter < 50 && res == CBF_RES_OK) {

    if((*ub)[*probevar] - (*lb)[*probevar] > 0.5) {
      if(!isinf((*lb)[*probevar])) {
        ++iter;

        // Reset
        probelb = *lb;
        probeub = *ub;
        condinfeas = false;
        condredund = false;
        res = CBFmapmaxviol_reset(data->mapnum, mapmaxviol);

        // Fix to lower
//        printf("Probing on UB[%lli] = %g\n", *probevar, (*lb)[*probevar]);
        probeub[*probevar] = (*lb)[*probevar];

        // Check
        if(res == CBF_RES_OK) {
          res = update_varbound_fast(
              data, fast_param, integerarray, probelb, probeub, NULL, NULL, mapmaxviol, NULL, &condinfeas, &condredund);
        }

        if(res == CBF_RES_OK && param.USE_PROBE_VARBOUNDING_FAST_INEQFREE) {
          res = update_varbound_fast(data,
                                     probing_param,
                                     integerarray,
                                     probelb,
                                     probeub,
                                     NULL,
                                     NULL,
                                     mapmaxviol,
                                     NULL,
                                     &condinfeas,
                                     &condredund);
        }

        if(condinfeas) {

          // Strengthening variable bounds
          (*lb)[*probevar] += 1;
          (*stronglb)[*probevar] = true;

        } else if(condredund && (*ub)[*probevar] - (*lb)[*probevar] <= 1.5) {

          // Improving conic constraint coefficients
          coefimprove.probevar = *probevar;
          coefimprove.probeval = (*lb)[*probevar];
          coefimprove.fixpoint = (*ub)[*probevar];
          coefimprove.mapmaxviol = mapmaxviol;

          if(res == CBF_RES_OK) {
            res = transform_stackwise(data, param, (CBFtransformmemory*)&coefimprove, NULL, transform_coefimprove);
          }
        }

        if(condinfeas) {
          // Stay with the same variable
          repeatedprobevar = true;

        } else {
          if(repeatedprobevar || condredund) {
            // Propagate changes
            res = fast_transform(data, fast_param, NULL);
            stepsize = 1;

          } else {
            // Fast forward if no progress
            stepsize *= 2;
          }

          repeatedprobevar = false;
        }

        if(!repeatedprobevar) {
          for(i = 1; i <= stepsize && probevar != probelist.end(); ++i) {
            ++probevar;
          }
        }
        continue;
      }
    } else if(repeatedprobevar) {
      // Interval has been closed
      repeatedprobevar = false;

      // Propagate changes
      res = fast_transform(data, fast_param, NULL);
      stepsize = 1;
    }

    if(!repeatedprobevar) {
      ++probevar;
    }
  }

  if(repeatedprobevar) {
    // Propagate changes
    res = fast_transform(data, fast_param, NULL);
  }

  probevar = probelist.begin();
  iter = 0;
  stepsize = 1;
  repeatedprobevar = false;
  while(probevar != probelist.end() && !(*infeas) && iter < 50 && res == CBF_RES_OK) {

    if((*ub)[*probevar] - (*lb)[*probevar] > 0.5) {
      if(!isinf((*ub)[*probevar])) {
        ++iter;

        // Reset
        probelb = *lb;
        probeub = *ub;
        condinfeas = false;
        condredund = false;
        res = CBFmapmaxviol_reset(data->mapnum, mapmaxviol);

        // Fix to upper
//        printf("Probing on LB[%lli] = %g\n", *probevar, (*ub)[*probevar]);
        probelb[*probevar] = (*ub)[*probevar];

        // Check
        if(res == CBF_RES_OK) {
          res = update_varbound_fast(
              data, fast_param, integerarray, probelb, probeub, NULL, NULL, mapmaxviol, NULL, &condinfeas, &condredund);
        }

        if(res == CBF_RES_OK && param.USE_PROBE_VARBOUNDING_FAST_INEQFREE) {
          res = update_varbound_fast(data,
                                     probing_param,
                                     integerarray,
                                     probelb,
                                     probeub,
                                     NULL,
                                     NULL,
                                     mapmaxviol,
                                     NULL,
                                     &condinfeas,
                                     &condredund);
        }

        if(condinfeas) {

          // Strengthening variable bounds
          (*ub)[*probevar] -= 1;
          (*strongub)[*probevar] = true;

        } else if(condredund && (*ub)[*probevar] - (*lb)[*probevar] <= 1.5) {

          // Improving conic constraint coefficients
          coefimprove.probevar = *probevar;
          coefimprove.probeval = (*ub)[*probevar];
          coefimprove.fixpoint = (*lb)[*probevar];
          coefimprove.mapmaxviol = mapmaxviol;

          if(res == CBF_RES_OK) {
            res = transform_stackwise(data, param, (CBFtransformmemory*)&coefimprove, NULL, transform_coefimprove);
          }
        }

        if(condinfeas) {
          // Stay with the same variable
          repeatedprobevar = true;

        } else {
          if(repeatedprobevar || condredund) {
            // Propagate changes
            res = fast_transform(data, fast_param, NULL);
            stepsize = 1;

          } else {
            // Fast forward if no progress
            stepsize *= 2;
          }

          repeatedprobevar = false;
        }

        if(!repeatedprobevar) {
          for(i = 1; i <= stepsize && probevar != probelist.end(); ++i) {
            ++probevar;
          }
        }
        continue;
      }
    } else if(repeatedprobevar) {
      // Interval has been closed
      repeatedprobevar = false;

      // Propagate changes
      res = fast_transform(data, fast_param, NULL);
      stepsize = 1;
    }

    if(!repeatedprobevar) {
      ++probevar;
    }
  }

  if(repeatedprobevar) {
    // Propagate changes
    res = fast_transform(data, fast_param, NULL);
  }

  return res;
}

static CBFresponsee transform_coefimprove(CBFtransform_param param,
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
  long long int bnnz, annz;
  double alpha, beta;

  if(data->mapstackdomain[k] == CBF_CONE_QUAD) {

    CBFtransformmemory_coefimprovee coefimprove = *((CBFtransformmemory_coefimprovee*)mem);

    if (coefimprove.mapmaxviol[rowbeg] <= -1e-6) {
      
      alpha = (-coefimprove.mapmaxviol[rowbeg]) / (coefimprove.fixpoint - coefimprove.probeval);
      beta = -alpha * coefimprove.fixpoint;
      
      if(res == CBF_RES_OK)
        res = CBF_findforward_map(data, rowbeg, NULL, &abeg, &bbeg);

      if(res == CBF_RES_OK)
        res = get_cone_entry_b(data, 1, rowbeg, bbeg, b, bnnz);

      if(res == CBF_RES_OK)
        res = get_cone_entry_A(data, 1, rowbeg, abeg, A, annz, NULL);
  
      b.coeffRef(0,0) += beta;
      A.coeffRef(0,coefimprove.probevar) += alpha;

      if(res == CBF_RES_OK)
        res = put_cone_entry_b(data, &newdata, rowbeg, bbeg, bbeg + bnnz - 1, b);

      if(res == CBF_RES_OK)
        res = put_cone_entry_A(data, &newdata, rowbeg, abeg, abeg + annz - 1, A, NULL);
    }
  }

  return res;
}
