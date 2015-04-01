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

#include "transform-analytics.h"
#include "transform-graph.h"
#include "transform-helper.h"
#include "transform-lindep-elimination-soc.h"
#include "transform.h"
#include "solver-mosek.h"
#include "cbf-helper.h"

#include <Eigen/SparseQR>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <stdio.h>
#include <math.h>
#include <deque>
#include <vector>
#include <map>
#include <iostream>

using namespace Eigen;
using std::deque;

clock_t timer_quadcone_entrywisedecoupled = 0;
clock_t timer_quadcone_inequalityfree = 0;
clock_t timer_propagation_graph = 0;

CBFresponsee compute_varbound_extreme(const CBFsolver* solver,
                                      CBFsolvermemory* mem,
                                      CBFobjsensee objsense,
                                      long long int var,
                                      double* bound);

void compute_minmapval(const CBFdata* data,
                       const long long int r,
                       long long int bend,
                       long long int aend,
                       VectorXd& lb,
                       VectorXd& ub,
                       double& mapval,
                       long long int& varexcluded);

void compute_maxmapval(const CBFdata* data,
                       const long long int r,
                       long long int bend,
                       long long int aend,
                       VectorXd& lb,
                       VectorXd& ub,
                       double& mapval,
                       long long int& varexcluded);

CBFresponsee fixvars_at_minmapval(const CBFdata* data,
                                  char* integerarray,
                                  VectorXd& lb,
                                  VectorXd& ub,
                                  std::vector<bool>* stronglb,
                                  std::vector<bool>* strongub,
                                  const long long int r,
                                  long long int aend,
                                  bool* anychange,
                                  long long int graphr,
                                  Graph* graph,
                                  long long int* fixvar,
                                  bool* infeas);

CBFresponsee fixvars_at_maxmapval(const CBFdata* data,
                                  char* integerarray,
                                  VectorXd& lb,
                                  VectorXd& ub,
                                  std::vector<bool>* stronglb,
                                  std::vector<bool>* strongub,
                                  const long long int r,
                                  long long int aend,
                                  bool* anychange,
                                  long long int graphr,
                                  Graph* graph,
                                  long long int* fixvar,
                                  bool* infeas);

CBFresponsee strongboundvars_at_minmapval(const CBFdata* data,
                                          std::vector<bool>* stronglb,
                                          std::vector<bool>* strongub,
                                          const long long int r,
                                          long long int aend);

CBFresponsee strongboundvars_at_maxmapval(const CBFdata* data,
                                          std::vector<bool>* stronglb,
                                          std::vector<bool>* strongub,
                                          const long long int r,
                                          long long int aend);

CBFresponsee update_varbounds_lincone_ineq(const CBFdata* data,
                                           char* integerarray,
                                           VectorXd& lb,
                                           VectorXd& ub,
                                           std::vector<bool>* stronglb,
                                           std::vector<bool>* strongub,
                                           const long long int r,
                                           long long int bend,
                                           long long int aend,
                                           const CBFscalarconee lincone,
                                           double extra_bval,
                                           bool* anychange,
                                           long long int graphr,
                                           Graph* graph,
                                           long long int* fixvar,
                                           bool* infeas,
                                           double* maxviol,
                                           bool strongboundvars_at_redundancy);

CBFresponsee update_varbounds_lincone_equality(const CBFdata* data,
                                               char* integerarray,
                                               VectorXd& lb,
                                               VectorXd& ub,
                                               std::vector<bool>* stronglb,
                                               std::vector<bool>* strongub,
                                               const long long int r,
                                               long long int bend,
                                               long long int aend,
                                               double extra_bval,
                                               bool* anychange,
                                               long long int graphr,
                                               Graph* graph,
                                               long long int* fixvar,
                                               bool* infeas,
                                               double* maxviol,
                                               bool strongboundvars_at_redundancy);

CBFresponsee update_varbounds_quadcone_entrywisedecoupled(const CBFdata* data,
                                                          char* integerarray,
                                                          VectorXd& lb,
                                                          VectorXd& ub,
                                                          std::vector<bool>* stronglb,
                                                          std::vector<bool>* strongub,
                                                          const long long int klen,
                                                          const long long int row,
                                                          long long int bend,
                                                          long long int aend,
                                                          bool* anychange,
                                                          long long int graphr,
                                                          Graph* graph,
                                                          long long int* fixvar,
                                                          bool* infeas,
                                                          double* maxviol,
                                                          bool strongboundvars_at_redundancy);

CBFresponsee update_varbounds_quadcone_inequalityfree(const CBFdata* data,
                                                      char* integerarray,
                                                      VectorXd& lb,
                                                      VectorXd& ub,
                                                      std::vector<bool>* stronglb,
                                                      std::vector<bool>* strongub,
                                                      const long long int klen,
                                                      const long long int row,
                                                      long long int bend,
                                                      long long int aend,
                                                      bool* anychange,
                                                      long long int graphr,
                                                      Graph* graph,
                                                      long long int* fixvar,
                                                      bool* infeas,
                                                      double* maxviol,
                                                      bool strongboundvars_at_redundancy);

CBFresponsee update_varbounds_strongcomponents(const CBFdata* data,
                                               char* integerarray,
                                               VectorXd& lb,
                                               VectorXd& ub,
                                               std::vector<bool>* stronglb,
                                               std::vector<bool>* strongub,
                                               bool* anychange,
                                               Graph* graph,
                                               long long int* fixvar,
                                               bool* infeas);

void update_lb(const CBFdata* data,
               long long int var,
               double varbound,
               bool isinteger,
               VectorXd& lb,
               VectorXd& ub,
               std::vector<bool>* stronglb,
               std::vector<bool>* strongub,
               bool* anychange,
               long long int graphr,
               Graph* graph,
               long long int* fixvar,
               bool* infeas);

void update_ub(const CBFdata* data,
               long long int var,
               double varbound,
               bool isinteger,
               VectorXd& lb,
               VectorXd& ub,
               std::vector<bool>* stronglb,
               std::vector<bool>* strongub,
               bool* anychange,
               long long int graphr,
               Graph* graph,
               long long int* fixvar,
               bool* infeas);

void update_onevar_lincone(const CBFdata* data,
                           long long int var,
                           double aval,
                           double bval,
                           const CBFscalarconee lincone,
                           bool isinteger,
                           VectorXd& lb,
                           VectorXd& ub,
                           std::vector<bool>* stronglb,
                           std::vector<bool>* strongub,
                           bool* anychange_lb,
                           bool* anychange_ub,
                           long long int graphr,
                           Graph* graph,
                           long long int* fixvar,
                           bool* infeas);

// -------------------------------------
// Function definitions
// -------------------------------------

CBFresponsee init_varbound(const CBFdata* data, VectorXd& lb, VectorXd& ub)
{
  // Set initial bound values
  lb.resize(data->varnum);
  lb.setConstant(-INFINITY);

  ub.resize(data->varnum);
  ub.setConstant(INFINITY);

  return CBF_RES_OK;
}

CBFresponsee compute_varbound_extreme(const CBFsolver* solver,
                                      CBFsolvermemory* mem,
                                      CBFobjsensee objsense,
                                      long long int var,
                                      double* bound)
{
  CBFresponsee res = CBF_RES_OK;

  res = solver->set_objsense(mem, objsense);

  if(res == CBF_RES_OK)
    res = solver->set_cj(mem, var, 1);

  if(res == CBF_RES_OK) {
    res = solver->optimize(mem);

    if(res == CBF_RES_OK) {
      res = solver->getprimalobj(mem, bound);

      // Safeguard
      if(res == CBF_RES_OK) {
        if(objsense == CBF_OBJ_MINIMIZE) {
          *bound -= 1e-6;
        } else {
          *bound += 1e-6;
        }
      }

    } else {

//      if(objsense == CBF_OBJ_MINIMIZE) {
//        printf("skipped LB[%lli]\n", var);
//      } else {
//        printf("skipped UB[%lli]\n", var);
//      }
      res = CBF_RES_OK;
    }
  }

  if(res == CBF_RES_OK)
    res = solver->set_cj(mem, var, 0);

  return res;
}

void compute_minmapval(const CBFdata* data,
                       const long long int r,
                       long long int bend,
                       long long int aend,
                       VectorXd& lb,
                       VectorXd& ub,
                       double& mapval,
                       long long int* aexcluded)
{
  mapval = 0;
  if(aexcluded) {
    *aexcluded = -1;
  }

  // ACOORD
  while(aend >= 0 && data->asubi[aend] == r) {

    if(data->aval[aend] < 0) {
      if(!isinf(ub[data->asubj[aend]])) {
        mapval += data->aval[aend] * ub[data->asubj[aend]];
      } else {
        if(aexcluded && *aexcluded == -1) {
          *aexcluded = aend;
        } else {
          mapval = -INFINITY;
          if(aexcluded) {
            *aexcluded = -1;
          }
          return;
        }
      }

    } else if(data->aval[aend] > 0) {
      if(!isinf(lb[data->asubj[aend]])) {
        mapval += data->aval[aend] * lb[data->asubj[aend]];
      } else {
        if(aexcluded && *aexcluded == -1) {
          *aexcluded = aend;
        } else {
          mapval = -INFINITY;
          if(aexcluded) {
            *aexcluded = -1;
          }
          return;
        }
      }
    }

    --aend;
  }

  // BCOORD
  if(bend >= 0 && data->bsubi[bend] == r) {
    mapval += data->bval[bend];
  }
}

void compute_maxmapval(const CBFdata* data,
                       const long long int r,
                       long long int bend,
                       long long int aend,
                       VectorXd& lb,
                       VectorXd& ub,
                       double& mapval,
                       long long int* aexcluded)
{
  mapval = 0;
  if(aexcluded) {
    *aexcluded = -1;
  }

  // ACOORD
  while(aend >= 0 && data->asubi[aend] == r) {

    if(data->aval[aend] > 0) {
      if(!isinf(ub[data->asubj[aend]])) {
        mapval += data->aval[aend] * ub[data->asubj[aend]];
      } else {
        if(aexcluded && *aexcluded == -1) {
          *aexcluded = aend;
        } else {
          mapval = INFINITY;
          if(aexcluded) {
            *aexcluded = -1;
          }
          return;
        }
      }

    } else if(data->aval[aend] < 0) {
      if(!isinf(lb[data->asubj[aend]])) {
        mapval += data->aval[aend] * lb[data->asubj[aend]];
      } else {
        if(aexcluded && *aexcluded == -1) {
          *aexcluded = aend;
        } else {
          mapval = INFINITY;
          if(aexcluded) {
            *aexcluded = -1;
          }
          return;
        }
      }
    }

    --aend;
  }

  // BCOORD
  if(bend >= 0 && data->bsubi[bend] == r) {
    mapval += data->bval[bend];
  }
}

CBFresponsee fixvars_at_minmapval(const CBFdata* data,
                                  char* integerarray,
                                  VectorXd& lb,
                                  VectorXd& ub,
                                  std::vector<bool>* stronglb,
                                  std::vector<bool>* strongub,
                                  const long long int r,
                                  long long int aend,
                                  bool* anychange,
                                  long long int graphr,
                                  Graph* graph,
                                  long long int* fixvar,
                                  bool* infeas)
{
  while(aend >= 0 && data->asubi[aend] == r) {

    if(data->aval[aend] < 0) {
      if(!isinf(ub[data->asubj[aend]])) {

        update_lb(data,
                  data->asubj[aend],
                  ub[data->asubj[aend]],
                  integerarray[data->asubj[aend]],
                  lb,
                  ub,
                  stronglb,
                  strongub,
                  anychange,
                  graphr,
                  graph,
                  fixvar,
                  infeas);

      } else {
        printf("fixvars_at_minmapval refused to fix variable at infinite bound!\n");
        return CBF_RES_ERR;
      }

    } else if(data->aval[aend] > 0) {
      if(!isinf(lb[data->asubj[aend]])) {

        update_ub(data,
                  data->asubj[aend],
                  lb[data->asubj[aend]],
                  integerarray[data->asubj[aend]],
                  lb,
                  ub,
                  stronglb,
                  strongub,
                  anychange,
                  graphr,
                  graph,
                  fixvar,
                  infeas);

      } else {
        printf("fixvars_at_minmapval refused to fix variable at infinite bound!\n");
        return CBF_RES_ERR;
      }
    }

    --aend;
  }

  return CBF_RES_OK;
}

CBFresponsee fixvars_at_maxmapval(const CBFdata* data,
                                  char* integerarray,
                                  VectorXd& lb,
                                  VectorXd& ub,
                                  std::vector<bool>* stronglb,
                                  std::vector<bool>* strongub,
                                  const long long int r,
                                  long long int aend,
                                  bool* anychange,
                                  long long int graphr,
                                  Graph* graph,
                                  long long int* fixvar,
                                  bool* infeas)
{
  while(aend >= 0 && data->asubi[aend] == r) {

    if(data->aval[aend] > 0) {
      if(!isinf(ub[data->asubj[aend]])) {

        update_lb(data,
                  data->asubj[aend],
                  ub[data->asubj[aend]],
                  integerarray[data->asubj[aend]],
                  lb,
                  ub,
                  stronglb,
                  strongub,
                  anychange,
                  graphr,
                  graph,
                  fixvar,
                  infeas);

      } else {
        printf("fixvars_at_maxmapval refused to fix variable at infinite bound!\n");
        return CBF_RES_ERR;
      }

    } else if(data->aval[aend] < 0) {
      if(!isinf(lb[data->asubj[aend]])) {

        update_ub(data,
                  data->asubj[aend],
                  lb[data->asubj[aend]],
                  integerarray[data->asubj[aend]],
                  lb,
                  ub,
                  stronglb,
                  strongub,
                  anychange,
                  graphr,
                  graph,
                  fixvar,
                  infeas);

      } else {
        printf("fixvars_at_maxmapval refused to fix variable at infinite bound!\n");
        return CBF_RES_ERR;
      }
    }

    --aend;
  }

  return CBF_RES_OK;
}

CBFresponsee strongboundvars_at_minmapval(const CBFdata* data,
                                          std::vector<bool>* stronglb,
                                          std::vector<bool>* strongub,
                                          const long long int r,
                                          long long int aend)
{
  while(aend >= 0 && data->asubi[aend] == r) {
    if(data->aval[aend] > 0) {
      (*stronglb)[data->asubj[aend]] = true;
    } else {
      (*strongub)[data->asubj[aend]] = true;
    }
    --aend;
  }
  return CBF_RES_OK;
}

CBFresponsee strongboundvars_at_maxmapval(const CBFdata* data,
                                          std::vector<bool>* stronglb,
                                          std::vector<bool>* strongub,
                                          const long long int r,
                                          long long int aend)
{
  while(aend >= 0 && data->asubi[aend] == r) {
    if(data->aval[aend] > 0) {
      (*strongub)[data->asubj[aend]] = true;
    } else {
      (*stronglb)[data->asubj[aend]] = true;
    }
    --aend;
  }
  return CBF_RES_OK;
}

CBFresponsee update_varbound_full(const CBFdata* data,
                                  char* integerarray,
                                  bool only_intvars,
                                  VectorXd& lb,
                                  VectorXd& ub,
                                  std::vector<bool>* stronglb,
                                  std::vector<bool>* strongub,
                                  long long int* fixvar,
                                  bool* infeas)
{
  CBFresponsee res = CBF_RES_OK;
  const CBFsolver solver = solver_mosek;
  CBFsolvermemory mem = {
    0,
  };
  CBFdata __pvarbound;
  CBFdyndata pvarbound = {
    &__pvarbound,
  };
  double varbound;

  bool localchange, anychange = true;
  long long int loopcnt = 0;
  long long int j;

  // Make a shallow copy
  // (dynamic allocations are used to keep data unchanged)
  *pvarbound.data = *data;

  // Relax integrality
  pvarbound.data->intvarnum = 0;
  pvarbound.data->intvar = NULL;

  // Set objective to zero
  pvarbound.data->objannz = 0;
  pvarbound.data->objbval = 0;

  if(res == CBF_RES_OK)
    res = solver.write(pvarbound.data, &mem);

  // Start from currently known bounds
  for(j = 0; j < data->varnum; ++j) {
    if(!isinf(lb[j])) {
      solver.set_lbj(&mem, j, lb[j]);
    }

    if(!isinf(ub[j])) {
      solver.set_ubj(&mem, j, ub[j]);
    }
  }

  while(anychange && !(*infeas)) {
    ++loopcnt;
    anychange = false;

    //
    // Get lower bounds
    //
    for(j = 0; j < data->varnum && !(*infeas) && res == CBF_RES_OK; ++j) {
      if(!only_intvars || integerarray[j]) {
        varbound = lb[j];
        res = compute_varbound_extreme(&solver, &mem, CBF_OBJ_MINIMIZE, j, &varbound);

        if(res == CBF_RES_OK) {
          localchange = false;
          update_lb(
              data, j, varbound, integerarray[j], lb, ub, stronglb, strongub, &localchange, -1, NULL, fixvar, infeas);

          // Can the solver learn something?
          // if (res == CBF_RES_OK && !(*infeas) && integerarray[j] && (lb[j] - varbound) >= 1e-3) {
          if(localchange && (*stronglb)[j]) {
//            printf("LEARN: varbound[%lli]=%g, lb[%lli]=%g\n", j, varbound, j, lb[j]);
            solver.set_lbj(&mem, j, lb[j]);
            anychange = true;
          }
        }
      }
    }

    //
    // Get upper bounds
    //
    for(j = 0; j < data->varnum && !(*infeas) && res == CBF_RES_OK; ++j) {
      if(!only_intvars || integerarray[j]) {
        varbound = ub[j];
        res = compute_varbound_extreme(&solver, &mem, CBF_OBJ_MAXIMIZE, j, &varbound);

        if(res == CBF_RES_OK) {
          localchange = false;
          update_ub(
              data, j, varbound, integerarray[j], lb, ub, stronglb, strongub, &localchange, -1, NULL, fixvar, infeas);

          // Can the solver learn something?
          // if (res == CBF_RES_OK && !(*infeas) && integerarray[j] && (varbound - ub[j]) >= 1e-3) {
          if(localchange && (*strongub)[j]) {
//            printf("LEARN: varbound[%lli]=%g, ub[%lli]=%g\n", j, varbound, j, ub[j]);
            solver.set_ubj(&mem, j, ub[j]);
            anychange = true;
          }
        }
      }
    }
  }

//  printf("EXTREME VALUE BOUNDING DONE IN LOOPS: %lli\n", loopcnt);

  solver.clean(&mem);
  CBFdyn_freedynamicallocations(&pvarbound);

  return res;
}

CBFresponsee update_varbound_fast(const CBFdata* data,
                                  CBFtransform_param& param,
                                  char* integerarray,
                                  VectorXd& lb,
                                  VectorXd& ub,
                                  std::vector<bool>* stronglb,
                                  std::vector<bool>* strongub,
                                  double* mapmaxviol,
                                  long long int* fixvar,
                                  bool* infeas,
                                  bool* conestrictredund)
{
  CBFresponsee res = CBF_RES_OK;
  long long int r, row, k, bend, aend, var;
  double bval, aval;
  bool localchange_lb, localchange_ub, anychange = true;
  long long int loopcnt = 1;
  const long long int maxloops = 5;
  double maxviol;

  //
  // This graph tracks changes in domain propagation, and aids identification of
  // cyclic updates
  //   map[1] -> lb[1] -> map[2] -> ub[2] -> map[1]
  // by strong component analysis. There is one node for each
  //   map,
  //   lb,
  //   ub,
  // with node ID's numbered from 0 in this order..
  //
  Graph graph(param.USE_VARBOUNDING_FAST_COMPONENTANALYSIS <= 0 ? 0 : data->mapnum + 2 * data->varnum);
  Graph* graphptr;

  while(res == CBF_RES_OK && anychange && !(*infeas) &&
        loopcnt <= maxloops + param.USE_VARBOUNDING_FAST_COMPONENTANALYSIS) {
    anychange = false;

    graphptr = (param.USE_VARBOUNDING_FAST_COMPONENTANALYSIS &&
                loopcnt == maxloops + param.USE_VARBOUNDING_FAST_COMPONENTANALYSIS) ?
                   &graph :
                   NULL;

    // CBF files often have simple bounds in bottom, so
    // to reverse sweep over rows to extract the bounds
    // should get more work done in the first loop-around.
    bend = data->bnnz - 1;
    aend = data->annz - 1;
    row = data->mapnum - 1;
    for(k = data->mapstacknum - 1; k >= 0 && !(*infeas) && res == CBF_RES_OK; --k) {
      switch(data->mapstackdomain[k]) {

      case CBF_CONE_QUAD:

        res = CBF_findbackward_map(data, row, NULL, &aend, &bend);

        if(res == CBF_RES_OK && param.USE_VARBOUNDING_FAST_INEQFREE && loopcnt == 1) {
          res = update_varbounds_quadcone_inequalityfree(data,
                                                         integerarray,
                                                         lb,
                                                         ub,
                                                         stronglb,
                                                         strongub,
                                                         data->mapstackdim[k],
                                                         row,
                                                         bend,
                                                         aend,
                                                         &anychange,
                                                         row,
                                                         graphptr,
                                                         fixvar,
                                                         infeas,
                                                         NULL,
                                                         false);
        }

        if(res == CBF_RES_OK && !(*infeas) && param.USE_VARBOUNDING_FAST_DECOUPLEDCONE) {
          update_varbounds_quadcone_entrywisedecoupled(data,
                                                       integerarray,
                                                       lb,
                                                       ub,
                                                       stronglb,
                                                       strongub,
                                                       data->mapstackdim[k],
                                                       row,
                                                       bend,
                                                       aend,
                                                       &anychange,
                                                       row,
                                                       graphptr,
                                                       fixvar,
                                                       infeas,
                                                       &maxviol,
                                                       true);

          if(mapmaxviol && maxviol <= 0) {
            for(r = row; r >= row - data->mapstackdim[k] + 1 && res == CBF_RES_OK; --r) {
              mapmaxviol[r] = std::min(mapmaxviol[r], maxviol);
            }
            
            if (conestrictredund && maxviol <= -1e-6) {
              *conestrictredund = true;
            }
          }
        }
        break;

      case CBF_CONE_POS:
      case CBF_CONE_NEG:
      case CBF_CONE_ZERO:

        // Process each row separately
        for(r = row; r >= row - data->mapstackdim[k] + 1 && !(*infeas) && res == CBF_RES_OK; --r) {
          res = CBF_findbackward_map(data, r, NULL, &aend, &bend);

          if(res == CBF_RES_OK) {

            // Special case: empty and singleton rows
            if(aend <= 0 || data->asubi[aend - 1] != r) {

              if(loopcnt == 1) {

                if(bend >= 0 && data->bsubi[bend] == r) {
                  bval = data->bval[bend];
                } else {
                  bval = 0;
                }

                if(aend >= 0 && data->asubi[aend] == r) {
                  aval = data->aval[aend];
                  var = data->asubj[aend];
                } else {
                  aval = 0;
                  var = -1;
                }

                localchange_lb = false;
                localchange_ub = false;

                update_onevar_lincone(data,
                                      var,
                                      aval,
                                      bval,
                                      data->mapstackdomain[k],
                                      (var != -1) && integerarray[var],
                                      lb,
                                      ub,
                                      stronglb,
                                      strongub,
                                      &localchange_lb,
                                      &localchange_ub,
                                      r,
                                      graphptr,
                                      fixvar,
                                      infeas);

                if(mapmaxviol) {
                  mapmaxviol[r] = std::min(mapmaxviol[r], 0.0);

                  if(localchange_lb && stronglb) {
                    (*stronglb)[data->asubj[aend]] = true;
                  }

                  if(localchange_ub && strongub) {
                    (*strongub)[data->asubj[aend]] = true;
                  }
                }
              }

            } else {

              // General case: Update variable bounds and check for redundancy
              if(data->mapstackdomain[k] == CBF_CONE_POS) {
                res = update_varbounds_lincone_ineq(data,
                                                    integerarray,
                                                    lb,
                                                    ub,
                                                    stronglb,
                                                    strongub,
                                                    r,
                                                    bend,
                                                    aend,
                                                    CBF_CONE_POS,
                                                    0.0,
                                                    &anychange,
                                                    r,
                                                    graphptr,
                                                    fixvar,
                                                    infeas,
                                                    &maxviol,
                                                    true);

              } else if(data->mapstackdomain[k] == CBF_CONE_NEG) {
                res = update_varbounds_lincone_ineq(data,
                                                    integerarray,
                                                    lb,
                                                    ub,
                                                    stronglb,
                                                    strongub,
                                                    r,
                                                    bend,
                                                    aend,
                                                    CBF_CONE_NEG,
                                                    0.0,
                                                    &anychange,
                                                    r,
                                                    graphptr,
                                                    fixvar,
                                                    infeas,
                                                    &maxviol,
                                                    true);

              } else {
                res = update_varbounds_lincone_equality(data,
                                                        integerarray,
                                                        lb,
                                                        ub,
                                                        stronglb,
                                                        strongub,
                                                        r,
                                                        bend,
                                                        aend,
                                                        0.0,
                                                        &anychange,
                                                        r,
                                                        graphptr,
                                                        fixvar,
                                                        infeas,
                                                        &maxviol,
                                                        false);
              }

              if(mapmaxviol && maxviol <= 0) {
                mapmaxviol[r] = std::min(mapmaxviol[r], maxviol);
              }
            }
          }
        }
        break;

      case CBF_CONE_FREE:
        if(mapmaxviol) {
          for(r = row; r >= row - data->mapstackdim[k] + 1; --r) {
            mapmaxviol[r] = std::min(mapmaxviol[r], 0.0);
          }
        }
        break;

      default:
        res = CBF_RES_ERR;
        break;
      }

      row -= data->mapstackdim[k];
    }

    if(res == CBF_RES_OK && !(*infeas) && anychange) {
      if(graphptr) {
//        printf("Solving subsystems\n");

        anychange = false;
        res = update_varbounds_strongcomponents(
            data, integerarray, lb, ub, stronglb, strongub, &anychange, graphptr, fixvar, infeas);

        if(anychange) {
          loopcnt = 0;
        }
      }
    }

    ++loopcnt;
  }

  return res;
}

CBFresponsee update_varbounds_lincone_ineq(const CBFdata* data,
                                           char* integerarray,
                                           VectorXd& lb,
                                           VectorXd& ub,
                                           std::vector<bool>* stronglb,
                                           std::vector<bool>* strongub,
                                           const long long int r,
                                           long long int bend,
                                           long long int aend,
                                           const CBFscalarconee lincone,
                                           double extra_bval,
                                           bool* anychange,
                                           long long int graphr,
                                           Graph* graph,
                                           long long int* fixvar,
                                           bool* infeas,
                                           double* maxviol,
                                           bool strongboundvars_at_redundancy)
{
  CBFresponsee res = CBF_RES_OK;
  bool localchange = false;
  long long int i, aexcluded;
  double maxslack, bval;

  // Validate input cone
  if(lincone != CBF_CONE_POS && lincone != CBF_CONE_NEG) {
    printf("update_varbounds_lincone_ineq does not support this cone type!");
    return CBF_RES_ERR;
  }

  // Compute maximum slack
  if(lincone == CBF_CONE_POS) {
    compute_maxmapval(data, r, bend, aend, lb, ub, maxslack, &aexcluded);
    maxslack += extra_bval;
  } else {
    compute_minmapval(data, r, bend, aend, lb, ub, maxslack, &aexcluded);
    maxslack += extra_bval;
  }

  if(!isinf(maxslack)) {
    //
    // Convert maxslack into bounds
    //

    if(aexcluded != -1) {

      // Interior intersection with variable bounding box
      // (can only tighten bounds on one variable)
      update_onevar_lincone(data,
                            data->asubj[aexcluded],
                            data->aval[aexcluded],
                            maxslack,
                            lincone,
                            integerarray[data->asubj[aexcluded]],
                            lb,
                            ub,
                            stronglb,
                            strongub,
                            &localchange,
                            &localchange,
                            graphr,
                            graph,
                            fixvar,
                            infeas);

    } else {

      if((lincone == CBF_CONE_POS && maxslack < -1e-6) || (lincone == CBF_CONE_NEG && maxslack > 1e-6)) {

        // No intersection with variable bounding box
//        printf("INFEAS:INEQ\n");
        *infeas = true;

      } else if((lincone == CBF_CONE_POS && maxslack < 1e-6) || (lincone == CBF_CONE_NEG && maxslack > -1e-6)) {

        // Boundary intersection with variable bounding box
        if(lincone == CBF_CONE_POS) {
          res = fixvars_at_maxmapval(
              data, integerarray, lb, ub, stronglb, strongub, r, aend, &localchange, graphr, graph, fixvar, infeas);

        } else {
          res = fixvars_at_minmapval(
              data, integerarray, lb, ub, stronglb, strongub, r, aend, &localchange, graphr, graph, fixvar, infeas);
        }

      } else {

        // Interior intersection with variable bounding box
        i = aend;
        while(i >= 0 && data->asubi[i] == r) {

          if(data->aval[i] != 0.0) {
            if((lincone == CBF_CONE_POS && data->aval[i] > 0) || (lincone == CBF_CONE_NEG && data->aval[i] < 0)) {
              bval = maxslack - data->aval[i] * ub[data->asubj[i]];
            } else {
              bval = maxslack - data->aval[i] * lb[data->asubj[i]];
            }

            update_onevar_lincone(data,
                                  data->asubj[i],
                                  data->aval[i],
                                  bval,
                                  lincone,
                                  integerarray[data->asubj[i]],
                                  lb,
                                  ub,
                                  stronglb,
                                  strongub,
                                  &localchange,
                                  &localchange,
                                  graphr,
                                  graph,
                                  fixvar,
                                  infeas);
          }

          --i;
        }
      }
    }

    // Update propagation graph
    if(localchange && graph) {
      i = aend;
      while(i >= 0 && data->asubi[i] == r) {
        if(i != aexcluded) {
          if((lincone == CBF_CONE_POS && data->aval[i] > 0) || (lincone == CBF_CONE_NEG && data->aval[i] < 0)) {
            // ub was used by row r
            graph->addArc(data->mapnum + data->varnum + data->asubj[i], graphr);

          } else if((lincone == CBF_CONE_POS && data->aval[i] < 0) || (lincone == CBF_CONE_NEG && data->aval[i] > 0)) {
            // lb was used by row r
            graph->addArc(data->mapnum + data->asubj[i], graphr);
          }
        }
        --i;
      }
    }

    *anychange |= localchange;
  }

  // Compute maximum violation
  if(maxviol && res == CBF_RES_OK) {
    if(lincone == CBF_CONE_POS) {
      compute_minmapval(data, r, bend, aend, lb, ub, *maxviol, NULL);
      *maxviol = -*maxviol;

      if(strongboundvars_at_redundancy && stronglb && strongub && *maxviol <= 0) {
        res = strongboundvars_at_minmapval(data, stronglb, strongub, r, aend);
      }

    } else {
      compute_maxmapval(data, r, bend, aend, lb, ub, *maxviol, NULL);

      if(strongboundvars_at_redundancy && stronglb && strongub && *maxviol <= 0)
        res = strongboundvars_at_maxmapval(data, stronglb, strongub, r, aend);
    }
  }

  return res;
}

CBFresponsee update_varbounds_lincone_equality(const CBFdata* data,
                                               char* integerarray,
                                               VectorXd& lb,
                                               VectorXd& ub,
                                               std::vector<bool>* stronglb,
                                               std::vector<bool>* strongub,
                                               const long long int r,
                                               long long int bend,
                                               long long int aend,
                                               double extra_bval,
                                               bool* anychange,
                                               long long int graphr,
                                               Graph* graph,
                                               long long int* fixvar,
                                               bool* infeas,
                                               double* maxviol,
                                               bool strongboundvars_at_redundancy)
{
  CBFresponsee res = CBF_RES_OK;
  double maxviol1, maxviol2;

  if(res == CBF_RES_OK) {
    res = update_varbounds_lincone_ineq(data,
                                        integerarray,
                                        lb,
                                        ub,
                                        stronglb,
                                        strongub,
                                        r,
                                        bend,
                                        aend,
                                        CBF_CONE_POS,
                                        extra_bval,
                                        anychange,
                                        graphr,
                                        graph,
                                        fixvar,
                                        infeas,
                                        &maxviol1,
                                        false);
  }

  if(res == CBF_RES_OK) {
    res = update_varbounds_lincone_ineq(data,
                                        integerarray,
                                        lb,
                                        ub,
                                        stronglb,
                                        strongub,
                                        r,
                                        bend,
                                        aend,
                                        CBF_CONE_NEG,
                                        extra_bval,
                                        anychange,
                                        graphr,
                                        graph,
                                        fixvar,
                                        infeas,
                                        &maxviol2,
                                        false);
  }

  // Compute maximum violation
  if(maxviol && res == CBF_RES_OK) {
    *maxviol = std::max(maxviol1, maxviol2);

    if(strongboundvars_at_redundancy && stronglb && strongub && *maxviol <= 0) {
      res = strongboundvars_at_minmapval(data, stronglb, strongub, r, aend);

      if(res == CBF_RES_OK)
        res = strongboundvars_at_maxmapval(data, stronglb, strongub, r, aend);
    }
  }

  return res;
}

CBFresponsee update_varbounds_quadcone_entrywisedecoupled(const CBFdata* data,
                                                          char* integerarray,
                                                          VectorXd& lb,
                                                          VectorXd& ub,
                                                          std::vector<bool>* stronglb,
                                                          std::vector<bool>* strongub,
                                                          const long long int klen,
                                                          const long long int row,
                                                          long long int bend,
                                                          long long int aend,
                                                          bool* anychange,
                                                          long long int graphr,
                                                          Graph* graph,
                                                          long long int* fixvar,
                                                          bool* infeas,
                                                          double* maxviol,
                                                          bool strongboundvars_at_redundancy)
{
  clock_t start = clock();
  
  CBFresponsee res = CBF_RES_OK;
  long long int r, k_aend = aend, k_bend = bend;
  std::vector<double> minmapval(klen);
  std::vector<double> maxmapval(klen);
  double mapval;
  double maxslack = 0;
  double maxslack_hyperball = 0;
  bool localchange = false;

  // Compute min/max map values
  for(r = row; klen - 1 - (row - r) >= 0; --r) {
    res = CBF_findbackward_map(data, r, NULL, &aend, &bend);

    compute_minmapval(data, r, bend, aend, lb, ub, minmapval[klen - 1 - (row - r)], NULL);
    compute_maxmapval(data, r, bend, aend, lb, ub, maxmapval[klen - 1 - (row - r)], NULL);
  }

  // Compute maximum slack
  for(r = 1; r < klen; ++r) {
    if((minmapval[r] > 0) || (maxmapval[r] < 0)) {
      mapval = std::min(std::fabs(minmapval[r]), std::fabs(maxmapval[r]));
      maxslack_hyperball += mapval * mapval;
    }
  }
  maxslack = std::sqrt(maxslack_hyperball) - maxmapval[0];

  if(isinf(maxmapval[0])) {

    // Interior intersection with variable bounding box
    // (can only tighten bounds on radius entry)

    res = update_varbounds_lincone_ineq(data,
                                        integerarray,
                                        lb,
                                        ub,
                                        stronglb,
                                        strongub,
                                        row - (klen - 1),
                                        bend,
                                        aend,
                                        CBF_CONE_POS,
                                        -std::sqrt(maxslack_hyperball),
                                        &localchange,
                                        graphr,
                                        graph,
                                        fixvar,
                                        infeas,
                                        NULL,
                                        false);

  } else {

    //
    // Convert maxslack into bounds
    //

    if(maxslack > 1e-6) {

      // No intersection with variable bounding box
//      printf("INFEAS:CONE:%lli\n", row);
      *infeas = true;

    } else if(maxslack > -1e-6) {

//      printf("FORCING:CONE:%lli\n", row);

      // Boundary intersection with variable bounding box
      aend = k_aend;
      bend = k_bend;
      for(r = row; klen - 1 - (row - r) > 0 && res == CBF_RES_OK; --r) {
        res = CBF_findbackward_map(data, r, NULL, &aend, &bend);

        if(res == CBF_RES_OK) {
          if(minmapval[klen - 1 - (row - r)] >= 0) {
            res = fixvars_at_minmapval(
                data, integerarray, lb, ub, stronglb, strongub, r, aend, &localchange, graphr, graph, fixvar, infeas);

          } else if(maxmapval[klen - 1 - (row - r)] <= 0) {
            res = fixvars_at_maxmapval(
                data, integerarray, lb, ub, stronglb, strongub, r, aend, &localchange, graphr, graph, fixvar, infeas);
          }
        }
      }

      // Radius entry: klen - 1 - (row - r) == 0
      if(res == CBF_RES_OK) {
        res = CBF_findbackward_map(data, r, NULL, &aend, &bend);

        if(res == CBF_RES_OK)
          res = fixvars_at_maxmapval(
              data, integerarray, lb, ub, stronglb, strongub, r, aend, &localchange, graphr, graph, fixvar, infeas);
      }

    } else {

      // Interior intersection with variable bounding box
      aend = k_aend;
      bend = k_bend;
      for(r = row; klen - 1 - (row - r) > 0 && res == CBF_RES_OK; --r) {
        res = CBF_findbackward_map(data, r, NULL, &aend, &bend);

        if(res == CBF_RES_OK) {
          if((minmapval[klen - 1 - (row - r)] > 0) || (maxmapval[klen - 1 - (row - r)] < 0)) {
            mapval = std::min(std::fabs(minmapval[klen - 1 - (row - r)]), std::fabs(maxmapval[klen - 1 - (row - r)]));
          } else {
            mapval = 0.0;
          }

          res = update_varbounds_lincone_ineq(
              data,
              integerarray,
              lb,
              ub,
              stronglb,
              strongub,
              r,
              bend,
              aend,
              CBF_CONE_NEG,
              -std::sqrt(maxmapval[0] * maxmapval[0] - (maxslack_hyperball - mapval * mapval)),
              &localchange,
              graphr,
              graph,
              fixvar,
              infeas,
              NULL,
              false);

          res = update_varbounds_lincone_ineq(
              data,
              integerarray,
              lb,
              ub,
              stronglb,
              strongub,
              r,
              bend,
              aend,
              CBF_CONE_POS,
              std::sqrt(maxmapval[0] * maxmapval[0] - (maxslack_hyperball - mapval * mapval)),
              &localchange,
              graphr,
              graph,
              fixvar,
              infeas,
              NULL,
              false);
        }
      }

      // Radius entry: klen - 1 - (row - r) == 0
      if(res == CBF_RES_OK) {
        res = CBF_findbackward_map(data, r, NULL, &aend, &bend);

        if(res == CBF_RES_OK) {
          res = update_varbounds_lincone_ineq(data,
                                              integerarray,
                                              lb,
                                              ub,
                                              stronglb,
                                              strongub,
                                              r,
                                              bend,
                                              aend,
                                              CBF_CONE_POS,
                                              -std::sqrt(maxslack_hyperball),
                                              &localchange,
                                              graphr,
                                              graph,
                                              fixvar,
                                              infeas,
                                              NULL,
                                              false);
        }
      }
    }
  }

  //  // Update propagation graph
  //  if(localchange && graph) {
  //    i = aend;
  //    while(i >= 0 && data->asubi[i] == r) {
  //      if(i != aexcluded) {
  //        if((lincone == CBF_CONE_POS && data->aval[i] > 0) || (lincone == CBF_CONE_NEG && data->aval[i] < 0)) {
  //          // ub was used by row r
  //          graph->addArc(data->mapnum + data->varnum + data->asubj[i], graphr);
  //
  //        } else if((lincone == CBF_CONE_POS && data->aval[i] < 0) || (lincone == CBF_CONE_NEG && data->aval[i] > 0))
  //        {
  //          // lb was used by row r
  //          graph->addArc(data->mapnum + data->asubj[i], graphr);
  //        }
  //      }
  //      --i;
  //    }
  //  }

  *anychange |= localchange;

  if(maxviol && res == CBF_RES_OK) {
    if(localchange) {
      // Recompute min/max map values
      aend = k_aend;
      bend = k_bend;
      for(r = row; klen - 1 - (row - r) >= 0; --r) {
        res = CBF_findbackward_map(data, r, NULL, &aend, &bend);

        compute_minmapval(data, r, bend, aend, lb, ub, minmapval[klen - 1 - (row - r)], NULL);
        compute_maxmapval(data, r, bend, aend, lb, ub, maxmapval[klen - 1 - (row - r)], NULL);
      }
    }

    // Compute maximum violation
    *maxviol = 0.0;
    for(r = 1; r < klen; ++r) {
      mapval = std::max(std::fabs(minmapval[r]), std::fabs(maxmapval[r]));
      *maxviol += mapval * mapval;
    }
    *maxviol = std::sqrt(*maxviol) - minmapval[0];

    if(*maxviol <= 0) {
//      printf("REDUNDANT:CONE:%lli=%g\n", row, *maxviol);

      if(strongboundvars_at_redundancy && stronglb && strongub) {
        aend = k_aend;
        for(r = row; klen - 1 - (row - r) >= 0; --r) {
          res = CBF_findbackward_map(data, r, NULL, &aend, NULL);

          if(res == CBF_RES_OK)
            res = strongboundvars_at_minmapval(data, stronglb, strongub, r, aend);

          if(res == CBF_RES_OK)
            res = strongboundvars_at_maxmapval(data, stronglb, strongub, r, aend);
        }

        // Radius entry: klen - 1 - (row - r) == 0
        if(res == CBF_RES_OK) {
          res = CBF_findbackward_map(data, r, NULL, &aend, NULL);

          if(res == CBF_RES_OK)
            res = strongboundvars_at_minmapval(data, stronglb, strongub, r, aend);
        }
      }
    }
  }

  timer_quadcone_entrywisedecoupled += (clock() - start);

  return res;
}

CBFresponsee update_varbounds_quadcone_inequalityfree(const CBFdata* data,
                                                      char* integerarray,
                                                      VectorXd& lb,
                                                      VectorXd& ub,
                                                      std::vector<bool>* stronglb,
                                                      std::vector<bool>* strongub,
                                                      const long long int klen,
                                                      const long long int row,
                                                      long long int bend,
                                                      long long int aend,
                                                      bool* anychange,
                                                      long long int graphr,
                                                      Graph* graph,
                                                      long long int* fixvar,
                                                      bool* infeas,
                                                      double* maxviol,
                                                      bool strongboundvars_at_redundancy)
{
  clock_t start = clock();
  
  CBFresponsee res = CBF_RES_OK;
  SparseQR<SparseMatrix<double, RowMajor>, COLAMDOrdering<int> > QR;
  std::map<long long int, long long int> compressedcolumnmap;
  SparseMatrix<double, RowMajor> b, brow, n, scal(1, 1), bval(1, 1);
  SparseMatrix<double, RowMajor> lambda;
  SparseMatrix<double, RowMajor> A, Arow;
  SparseMatrix<double, RowMajor> ej_sparse;
  VectorXd ej;
  double aval;
  long long int abeg, bbeg, rowbeg, annz, bnnz;
  long long int j, r;
  redtype_t redtype = none;
  double* reduced_mapmaxviol = NULL;

  CBFdata newdata = {
    0,
  };
  CBFdyndata dyn_newdata = {
    &newdata,
  };

  rowbeg = row - (klen - 1);
  res = CBF_findbackward_map(data, rowbeg - 1, NULL, &aend, &bend);
  abeg = aend + 1;
  bbeg = bend + 1;

  if(res == CBF_RES_OK)
    res = get_cone_entry_b(data, klen, rowbeg, bbeg, b, bnnz);

  if(res == CBF_RES_OK)
    res = get_cone_entry_A(data, klen, rowbeg, abeg, A, annz, &compressedcolumnmap);

  if(res == CBF_RES_OK)
    res = CBFmapmaxviol_init(klen, &reduced_mapmaxviol);

  if(res == CBF_RES_OK) {
    strip_cone_entry_fixvarnnz(A, b, &compressedcolumnmap, lb.data(), ub.data(), 1e-9);
    res = lindep_elimination_simple(0, klen, b, A, reduced_mapmaxviol, redtype, NULL);
  }

  if(res == CBF_RES_OK) {
    switch(redtype) {
    case point:
//      std::cout << "FORCING:POINT:CONE:" << row << std::endl;

      res = put_cone_entry_A(&newdata, &dyn_newdata, 0, 0, -1, A, &compressedcolumnmap);

      if(res == CBF_RES_OK)
        res = put_cone_entry_b(&newdata, &dyn_newdata, 0, 0, -1, b);

      aend = newdata.annz - 1;
      bend = newdata.bnnz - 1;

      // All entries
      for(r = klen - 1; r >= 0 && res == CBF_RES_OK; --r) {
        if(reduced_mapmaxviol[r] > 0) {
          res = CBF_findbackward_map(&newdata, r, NULL, &aend, &bend);

          if(res == CBF_RES_OK)
            res = update_varbounds_lincone_equality(&newdata,
                                                    integerarray,
                                                    lb,
                                                    ub,
                                                    stronglb,
                                                    strongub,
                                                    r,
                                                    bend,
                                                    aend,
                                                    0.0,
                                                    anychange,
                                                    row,
                                                    graph,
                                                    fixvar,
                                                    infeas,
                                                    NULL,
                                                    false);
        }
      }
      break;

    case line:
//      std::cout << "FORCING:LINE:CONE:" << row << std::endl;

      res = put_cone_entry_A(&newdata, &dyn_newdata, 0, 0, -1, A, &compressedcolumnmap);

      if(res == CBF_RES_OK)
        res = put_cone_entry_b(&newdata, &dyn_newdata, 0, 0, -1, b);

      aend = newdata.annz - 1;
      bend = newdata.bnnz - 1;

      // Hyperball entries
      for(r = klen - 1; r >= 1 && res == CBF_RES_OK; --r) {
        if(reduced_mapmaxviol[r] > 0) {
          res = CBF_findbackward_map(&newdata, r, NULL, &aend, &bend);

          if(res == CBF_RES_OK)
            res = update_varbounds_lincone_equality(&newdata,
                                                    integerarray,
                                                    lb,
                                                    ub,
                                                    stronglb,
                                                    strongub,
                                                    r,
                                                    bend,
                                                    aend,
                                                    0.0,
                                                    anychange,
                                                    row,
                                                    graph,
                                                    fixvar,
                                                    infeas,
                                                    NULL,
                                                    false);
        }
      }

      // Radius entry: r == 0
      if(res == CBF_RES_OK) {
        if(reduced_mapmaxviol[r] > 0) {
          res = CBF_findbackward_map(&newdata, r, NULL, &aend, &bend);

          if(res == CBF_RES_OK)
            res = update_varbounds_lincone_equality(&newdata,
                                                    integerarray,
                                                    lb,
                                                    ub,
                                                    stronglb,
                                                    strongub,
                                                    r,
                                                    bend,
                                                    aend,
                                                    0.0,
                                                    anychange,
                                                    row,
                                                    graph,
                                                    fixvar,
                                                    infeas,
                                                    NULL,
                                                    false);
        }
      }
      break;

    case none:
      //
      // Analyse coupling in the quadratic cone
      //

      if(res == CBF_RES_OK) {
        b = b.transpose();

        A.makeCompressed();
        QR.compute(A);

        for(j = 1; j <= QR.rank(); ++j) {
          ej.setZero(klen);
          ej[klen - j] = 1;
          n = (QR.matrixQ() * ej).sparseView();
          n = n.transpose();

          if(row_lexorder(0.0, 0.0, rowiter(n, 0), rowiter(b, 0), 1.0, 1.0) != 0) {
            scal = -(n * b.transpose());
            break;
          }
        }

        if(scal.coeff(0, 0) != 0.0) {
          if(is_quadcone_member_and_nonzero(rowiter(n, 0)) == (scal.coeff(0, 0) > 0 ? 1 : -1)) {
//            printf("INFEAS:COUPLING:CONE:%lli\n", row);
            *infeas = true;
          }
        }
      }

      if(res == CBF_RES_OK && !(*infeas)) {
        // Bound analysis
        ej_sparse.resize(1, A.cols());

        for(rowiter it(A, 0); it; ++it) {
          ej_sparse.setZero();
          ej_sparse.coeffRef(0, it.col()) = 1;
          ej = (ej_sparse * QR.colsPermutation()).transpose().topRows(QR.rank());

          lambda = (QR.matrixQ() *
                    (MatrixXd::Identity(klen, QR.rank()) *
                     (QR.matrixR().topLeftCorner(QR.rank(), QR.rank()).transpose().triangularView<Lower>().solve(ej))))
                       .sparseView()
                       .transpose();

          if(row_lexorder(0.0, 0.0, rowiter(lambda * A, 0), rowiter(ej_sparse, 0), 1.0, 1.0) == 0) {
            aval = (double)is_quadcone_member_and_nonzero(rowiter(lambda, 0));
            if(aval) {

              bval = (aval >= 1 ? lambda * b.transpose() : -lambda * b.transpose());

              update_onevar_lincone(data,
                                    compressedcolumnmap[it.col()],
                                    aval,
                                    bval.coeff(0, 0),
                                    CBF_CONE_POS,
                                    integerarray[compressedcolumnmap[it.col()]],
                                    lb,
                                    ub,
                                    stronglb,
                                    strongub,
                                    anychange,
                                    anychange,
                                    row,
                                    graph,
                                    fixvar,
                                    infeas);
            }
          }
        }
      }
      break;
    }
  }

  CBFmapmaxviol_free(&reduced_mapmaxviol);

  timer_quadcone_inequalityfree += (clock() - start);

  return res;
}

CBFresponsee update_varbounds_strongcomponents(const CBFdata* data,
                                               char* integerarray,
                                               VectorXd& lb,
                                               VectorXd& ub,
                                               std::vector<bool>* stronglb,
                                               std::vector<bool>* strongub,
                                               bool* anychange,
                                               Graph* graph,
                                               long long int* fixvar,
                                               bool* infeas)
{
  clock_t start = clock();
  
  CBFresponsee res = CBF_RES_OK;
  const CBFsolver solver = solver_mosek;
  CBFsolvermemory mem = {
    0,
  };
  long long int slen, siter, x, mapnum = 0, varnum = 0;
  deque<long long int>* scc;
  std::deque<long long int>::reverse_iterator i;
  CBFdyndata dyndata;
  CBFdata __data = {
    0,
  };
  CBFdyn_assign(&dyndata, &__data);

  std::map<long long int, long long int> mapID;
  std::map<long long int, long long int> varID;

  long long int k, row, bbeg, abeg, var;
  double varbound;

  // Loop over all scc in graph
  scc = new deque<long long int>();
  graph->SCC(data->mapnum, scc);

  while(!scc->empty() && res == CBF_RES_OK) {
    slen = -(long long int)scc->back();
    scc->pop_back();

    // Add constraints (variables added with their bounds on first visit)
    i = scc->rbegin();
    for(siter = 1; siter <= slen && res == CBF_RES_OK; ++siter) {
      if(*i < data->mapnum) {

        // Preallocate constraint
        res = CBFdyn_map_capacitysurplus(&dyndata, 1);

        if(res == CBF_RES_OK) {

          // Find cone of constraint identified by row '*i'
          row = data->mapnum - 1;
          k = data->mapstacknum - 1;
          while(k >= 0 && *i <= row - data->mapstackdim[k] && res == CBF_RES_OK) {
            row -= data->mapstackdim[k];
            --k;
          }

          switch(data->mapstackdomain[k]) {

          case CBF_CONE_QUAD:

            res = CBFdyn_map_adddomain(&dyndata, data->mapstackdomain[k], data->mapstackdim[k]);
            mapnum = dyndata.data->mapnum;

            if(res == CBF_RES_OK) {
              bbeg = data->bnnz - 1;
              abeg = data->annz - 1;
              res = CBF_findbackward_map(data, row, NULL, &abeg, &bbeg);
            }

            while(bbeg >= 0 && row - data->bsubi[bbeg] < data->mapstackdim[k] && res == CBF_RES_OK) {
              CBFdyn_b_capacitysurplus(&dyndata, 1);
              CBFdyn_b_add(&dyndata, mapnum - (row - data->bsubi[bbeg]) - 1, data->bval[bbeg]);
              bbeg--;
            }

            while(abeg >= 0 && row - data->asubi[abeg] < data->mapstackdim[k] && res == CBF_RES_OK) {

              // Add variable if it does not exists
              if(varID.find(data->asubj[abeg]) == varID.end()) {
                varID[data->asubj[abeg]] = varnum++;

                res = CBFdyn_var_capacitysurplus(&dyndata, 1);

                if(res == CBF_RES_OK)
                  res = CBFdyn_var_adddomain(&dyndata, CBF_CONE_FREE, 1);

                if(!isinf(lb[data->asubj[abeg]]) && res == CBF_RES_OK) {
                  res = CBFdyn_varbound_capacitysurplus(&dyndata, 1);

                  if(res == CBF_RES_OK) {
                    res = CBFdyn_varbound_addlower(&dyndata, varID[data->asubj[abeg]], lb[data->asubj[abeg]]);
                  }
                }

                if(!isinf(ub[data->asubj[abeg]]) && res == CBF_RES_OK) {
                  res = CBFdyn_varbound_capacitysurplus(&dyndata, 1);

                  if(res == CBF_RES_OK) {
                    res = CBFdyn_varbound_addupper(&dyndata, varID[data->asubj[abeg]], ub[data->asubj[abeg]]);
                  }
                }
              }

              // Set variable coefficient
              CBFdyn_a_capacitysurplus(&dyndata, 1);
              CBFdyn_a_add(
                  &dyndata, mapnum - (row - data->asubi[abeg]) - 1, varID[data->asubj[abeg]], data->aval[abeg]);
              abeg--;
            }

            break;

          case CBF_CONE_POS:
          case CBF_CONE_NEG:
          case CBF_CONE_ZERO:

            res = CBFdyn_map_adddomain(&dyndata, data->mapstackdomain[k], 1);
            mapnum = dyndata.data->mapnum;

            if(res == CBF_RES_OK) {
              bbeg = data->bnnz - 1;
              abeg = data->annz - 1;
              res = CBF_findbackward_map(data, *i, NULL, &abeg, &bbeg);
            }

            if(bbeg >= 0 && data->bsubi[bbeg] == *i && res == CBF_RES_OK) {
              CBFdyn_b_capacitysurplus(&dyndata, 1);
              CBFdyn_b_add(&dyndata, mapnum - 1, data->bval[bbeg]);
              bbeg--;
            }

            while(abeg >= 0 && data->asubi[abeg] == *i && res == CBF_RES_OK) {

              // Add variable if it does not exists
              if(varID.find(data->asubj[abeg]) == varID.end()) {
                varID[data->asubj[abeg]] = varnum++;

                res = CBFdyn_var_capacitysurplus(&dyndata, 1);

                if(res == CBF_RES_OK)
                  res = CBFdyn_var_adddomain(&dyndata, CBF_CONE_FREE, 1);

                if(!isinf(lb[data->asubj[abeg]]) && res == CBF_RES_OK) {
                  res = CBFdyn_varbound_capacitysurplus(&dyndata, 1);

                  if(res == CBF_RES_OK) {
                    res = CBFdyn_varbound_addlower(&dyndata, varID[data->asubj[abeg]], lb[data->asubj[abeg]]);
                  }
                }

                if(!isinf(ub[data->asubj[abeg]]) && res == CBF_RES_OK) {
                  res = CBFdyn_varbound_capacitysurplus(&dyndata, 1);

                  if(res == CBF_RES_OK) {
                    res = CBFdyn_varbound_addupper(&dyndata, varID[data->asubj[abeg]], ub[data->asubj[abeg]]);
                  }
                }
              }

              // Set variable coefficient
              CBFdyn_a_capacitysurplus(&dyndata, 1);
              CBFdyn_a_add(&dyndata, mapnum - 1, varID[data->asubj[abeg]], data->aval[abeg]);
              abeg--;
            }
            break;

          default:
            res = CBF_RES_ERR;
          }
        }
      }

      ++i;
    }

    // Write scc constraints
    if(res == CBF_RES_OK) {
      res = solver_mosek.write(dyndata.data, &mem);
      //      solver_mosek.savetodisk(&mem, "tmp.opf");
    }

    // Extremal bounding
    for(siter = 1; siter <= slen && res == CBF_RES_OK; ++siter) {
      x = (long long int)scc->back();
      scc->pop_back();

      if(x < data->mapnum) {
        // Skip constraints

      } else if(x < data->mapnum + data->varnum) {
        // Compute extremal LB
        var = x - data->mapnum;
        varbound = lb[var];
        res = compute_varbound_extreme(&solver, &mem, CBF_OBJ_MINIMIZE, varID[var], &varbound);

        if(res == CBF_RES_OK) {
          update_lb(
              data, var, varbound, integerarray[var], lb, ub, stronglb, strongub, anychange, -1, NULL, fixvar, infeas);
        }

      } else {
        // Compute extremal UB
        var = x - data->mapnum - data->varnum;
        varbound = ub[var];
        res = compute_varbound_extreme(&solver, &mem, CBF_OBJ_MAXIMIZE, varID[var], &varbound);

        if(res == CBF_RES_OK) {
          update_ub(
              data, var, varbound, integerarray[var], lb, ub, stronglb, strongub, anychange, -1, NULL, fixvar, infeas);
        }
      }
    }

    CBFdyn_freedynamicallocations(&dyndata);
  }

  // Clean up
  delete scc;

  timer_propagation_graph += (clock() - start);

  return res;
}

void update_lb(const CBFdata* data,
               long long int var,
               double varbound,
               bool isinteger,
               VectorXd& lb,
               VectorXd& ub,
               std::vector<bool>* stronglb,
               std::vector<bool>* strongub,
               bool* anychange,
               long long int graphr,
               Graph* graph,
               long long int* fixvar,
               bool* infeas)
{
  double frac = 0;

  if(isinteger) {
    frac = ceil(varbound - 1e-3);
    varbound = frac - varbound;
    std::swap(frac, varbound);
  }

  if(varbound > lb[var] + 1e-6) {
    //    printf("%lli: LB[%lli] = %g (from %g) with frac=%g\n", graphr, var, varbound, lb[var], frac);

    lb[var] = varbound;
    *anychange = true;

    if(lb[var] >= ub[var] + 1e-9) {
//      printf("INFEAS:LB\n");
      *infeas = true;

    } else {
      if(stronglb)
        (*stronglb)[var] = (frac >= 1e-3);

      if(fixvar && lb[var] >= ub[var] - 1e-9)
        ++(*fixvar);

      if(graph)
        graph->addArc(graphr, data->mapnum + var);
    }
  }
}

void update_ub(const CBFdata* data,
               long long int var,
               double varbound,
               bool isinteger,
               VectorXd& lb,
               VectorXd& ub,
               std::vector<bool>* stronglb,
               std::vector<bool>* strongub,
               bool* anychange,
               long long int graphr,
               Graph* graph,
               long long int* fixvar,
               bool* infeas)
{
  double frac = 0;

  if(isinteger) {
    frac = floor(varbound + 1e-3);
    varbound = varbound - frac;
    std::swap(frac, varbound);
  }

  if(varbound < ub[var] - 1e-6) {
    //    printf("%lli: UB[%lli] = %g (from %g) with frac=%g\n", graphr, var, varbound, ub[var], frac);

    ub[var] = varbound;
    *anychange = true;

    if(lb[var] >= ub[var] + 1e-9) {
//      printf("INFEAS:UB\n");
      *infeas = true;

    } else {
      if(strongub)
        (*strongub)[var] = (frac >= 1e-3);

      if(fixvar && lb[var] >= ub[var] - 1e-9)
        ++(*fixvar);

      if(graph)
        graph->addArc(graphr, data->mapnum + data->varnum + var);
    }
  }
}

void update_onevar_lincone(const CBFdata* data,
                           long long int var,
                           double aval,
                           double bval,
                           const CBFscalarconee lincone,
                           bool isinteger,
                           VectorXd& lb,
                           VectorXd& ub,
                           std::vector<bool>* stronglb,
                           std::vector<bool>* strongub,
                           bool* anychange_lb,
                           bool* anychange_ub,
                           long long int graphr,
                           Graph* graph,
                           long long int* fixvar,
                           bool* infeas)
{
  double bnd;

  if(aval != 0.0) {
    bnd = -bval / aval;

    if(lincone == CBF_CONE_ZERO) {
      update_lb(data, var, bnd, isinteger, lb, ub, stronglb, strongub, anychange_lb, graphr, graph, fixvar, infeas);
      update_ub(data, var, bnd, isinteger, lb, ub, stronglb, strongub, anychange_ub, graphr, graph, fixvar, infeas);

    } else if((lincone == CBF_CONE_POS && aval > 0) || (lincone == CBF_CONE_NEG && aval < 0)) {
      update_lb(data, var, bnd, isinteger, lb, ub, stronglb, strongub, anychange_lb, graphr, graph, fixvar, infeas);

    } else if((lincone == CBF_CONE_POS && aval < 0) || (lincone == CBF_CONE_NEG && aval > 0)) {
      update_ub(data, var, bnd, isinteger, lb, ub, stronglb, strongub, anychange_ub, graphr, graph, fixvar, infeas);
    }

  } else if(infeas) {

    if((lincone == CBF_CONE_ZERO && bval != 0) || (lincone == CBF_CONE_POS && bval < 0) ||
       (lincone == CBF_CONE_NEG && bval > 0)) {
//      printf("INFEAS:EMPTY LINCONE\n");
      *infeas = true;
    }
  }
}
