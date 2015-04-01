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

#include "transform-placeholders.h"
#include "cbf-helper.h"
#include "transform-helper.h"
#include "transform-analytics.h"

#include <Eigen/SparseCore>

#include <stdio.h>
#include <iostream>

using namespace Eigen;

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag);

struct coord_t
{
  long long int rbeg;
  long long int rend;
  long long int abeg;
  long long int bbeg;
};

struct full_coord_t
{
  long long int rbeg;
  long long int rend;
  long long int abeg;
  long long int aend;
  long long int bbeg;
  long long int bend;
};

// -------------------------------------
// Global variable
// -------------------------------------

CBFtransform const transform_placeholders = {
  "placeholders",
  transform,
};

char* integerarray = NULL;
std::vector<bool>* stronglb = NULL;
std::vector<bool>* strongub = NULL;

// -------------------------------------
// Function definitions
// -------------------------------------

void transform_placeholders_init(char* integerarray__, std::vector<bool>* stronglb__, std::vector<bool>* strongub__)
{
  integerarray = integerarray__;
  stronglb = stronglb__;
  strongub = strongub__;
}

void lookupvar(CBFdata* data,
               std::vector<long long int>& varmap,
               std::vector<double>& varscal,
               long long int& var,
               double& scal);

static CBFresponsee transform_nnzsweep(CBFdata* data,
                                CBFtransform_param& param,
                                bool* changeflag,
                                bool* incompletetransformation);

static CBFresponsee transform(CBFdata* data, CBFtransform_param& param, bool* changeflag)
{
  CBFresponsee res = CBF_RES_OK;
  bool incompletetransform = true;

  while(incompletetransform && res == CBF_RES_OK) {
    incompletetransform = false;
    res = transform_nnzsweep(data, param, changeflag, &incompletetransform);
  }

  return res;
}

static CBFresponsee transform_nnzsweep(CBFdata* data,
                                       CBFtransform_param& param,
                                       bool* changeflag,
                                       bool* incompletetransformation)
{
  CBFresponsee res = CBF_RES_OK;
  long long int i, k, klen, k_bbeg, k_abeg, rowbeg, r, var, bbeg, abeg, bnnz, annz;
  long long int src_aidx, dst_aidx;
  std::map<long long int, long long int> rowmap;
  std::map<long long int, long long int>::const_iterator it_rowmap;
  SparseMatrix<double, RowMajor> tmpA, srcA;
  SparseMatrix<double, RowMajor> tmpb, srcb;
  SparseMatrix<double, RowMajor> objA;
  SparseMatrix<double, RowMajor> rowsA(0, data->varnum);
  SparseMatrix<double, RowMajor> rowsb(0, 1);
  full_coord_t dstcoord;
  double srccoeff, dstcoeff, scal;
  double* mapmaxviol = NULL;
  bool anychange = false, anyrenaming = false;

  CBFdata newdata = {
    0,
  };
  CBFdyndata dyn_newdata = {
    &newdata,
  };
  CBFdyndata dyndata;

  if(res == CBF_RES_OK)
    res = CBFmapmaxviol_init(data->mapnum, &mapmaxviol);

  // Collect 'coord_t' information from all source (always L=) and destination rows
  const coord_t unset = { -1, -1, -1, -1 };
  const coord_t fail = { -2, -2, -2, -2 };

  std::vector<coord_t> src(data->varnum, unset);
  std::vector<coord_t> dst(data->varnum, unset);
  std::vector<full_coord_t> dst_compressed;
  std::vector<double> varscal(data->varnum, 1.0);
  std::vector<long long int> varmap(data->varnum);
  for(i = 0; i < data->varnum; ++i)
    varmap[i] = i;

  rowbeg = 0;
  bbeg = 0;
  abeg = 0;
  for(k = 0; k < data->mapstacknum && res == CBF_RES_OK; ++k) {

    res = CBF_findforward_map(data, rowbeg, NULL, &abeg, &bbeg);
    klen = data->mapstackdim[k];
    k_bbeg = bbeg;
    k_abeg = abeg;

    for(r = rowbeg; r <= rowbeg + klen - 1 && res == CBF_RES_OK; ++r) {
      res = CBF_findforward_map(data, r, NULL, &abeg, &bbeg);

      if(res == CBF_RES_OK) {
        i = abeg;

        if(data->mapstackdomain[k] == CBF_CONE_ZERO && (bbeg >= data->bnnz || data->bsubi[bbeg] != r) &&
           (i + 2 == data->annz || (i + 2 < data->annz && data->asubi[i + 1] == r && data->asubi[i + 2] != r)) &&
           (!integerarray || (!integerarray[data->asubj[i]] && (!stronglb || !(*stronglb)[data->asubj[i]]) &&
                              (!strongub || !(*strongub)[data->asubj[i]])) ||
            (!integerarray[data->asubj[i + 1]] && (!stronglb || !(*stronglb)[data->asubj[i + 1]]) &&
             (!strongub || !(*strongub)[data->asubj[i + 1]])))) {

          // Simple renaming
          anychange = true;
          anyrenaming = true;

          if(!integerarray || (!integerarray[data->asubj[i + 1]] && (!stronglb || !(*stronglb)[data->asubj[i + 1]]) &&
                               (!strongub || !(*strongub)[data->asubj[i + 1]]))) {
            src_aidx = i + 1;
            dst_aidx = i;
          } else {
            src_aidx = i;
            dst_aidx = i + 1;
          }

          var = data->asubj[dst_aidx];
          scal = -data->aval[dst_aidx] / data->aval[src_aidx];
          lookupvar(data, varmap, varscal, var, scal);

          varmap[data->asubj[src_aidx]] = var;
          varscal[data->asubj[src_aidx]] = scal;

          mapmaxviol[r] = std::min(mapmaxviol[r], 0.0);

          // Skip full placeholder detection for these variables
          *incompletetransformation = true;
          dst[var].rend = fail.rend;
          dst[data->asubj[src_aidx]].rend = fail.rend;

        } else {

          // Full placeholder detection
          while(i < data->annz && data->asubi[i] == r) {
            if(data->aval[i] != 0.0) {
              var = data->asubj[i];

              if(src[var].rend == unset.rend && data->mapstackdomain[k] == CBF_CONE_ZERO &&
                 (!integerarray || !integerarray[var]) && (!stronglb || !(*stronglb)[var]) &&
                 (!strongub || !(*strongub)[var])) {
                src[var].rbeg = r;
                src[var].rend = r;
                src[var].abeg = abeg;
                src[var].bbeg = bbeg;

              } else if(dst[var].rend == unset.rend) {

                if(data->mapstackdomain[k] == CBF_CONE_QUAD) {
                  dst[var].rbeg = rowbeg;
                  dst[var].rend = rowbeg + klen - 1;
                  dst[var].abeg = k_abeg;
                  dst[var].bbeg = k_bbeg;
                } else {
                  dst[var].rbeg = r;
                  dst[var].rend = r;
                  dst[var].abeg = abeg;
                  dst[var].bbeg = bbeg;
                }

              } else if(r > dst[var].rend) {
                dst[var].rend = fail.rend;
              }
            }

            ++i;
          }
        }
      }
    }

    rowbeg += klen;
  }

  // Apply eliminations
  if(res == CBF_RES_OK)
    res = get_obj_entry_A(data, objA, annz);

  for(var = data->varnum - 1; var >= 0 && res == CBF_RES_OK; --var) {
    if(dst[var].rend != fail.rend && src[var].rend != unset.rend && mapmaxviol[src[var].rbeg] > 0) {

      klen = dst[var].rend - dst[var].rbeg + 1;

      // Load source row (L=)
      if(res == CBF_RES_OK) {
        res = get_cone_entry_A(data, 1, src[var].rbeg, src[var].abeg, srcA, annz, NULL);
      }

      if(res == CBF_RES_OK) {
        res = get_cone_entry_b(data, 1, src[var].rbeg, src[var].bbeg, srcb, bnnz);
      }

      if(res == CBF_RES_OK) {
        srccoeff = srcA.coeff(0, var);

        if(srccoeff == 0.0) {
          // Other elimination has interferred
          continue;
        }
      }

      // Register destination rows
      if(dst[var].rend != unset.rend && res == CBF_RES_OK) {

        it_rowmap = rowmap.find(dst[var].rbeg);
        if(it_rowmap == rowmap.end()) {
          it_rowmap = rowmap.insert(std::pair<long long int, long long int>(dst[var].rbeg, rowsA.rows())).first;

          rowsA.conservativeResize(rowsA.rows() + klen, data->varnum);
          rowsb.conservativeResize(rowsb.rows() + klen, 1);

          res = get_cone_entry_A(data, klen, dst[var].rbeg, dst[var].abeg, tmpA, annz, NULL);
          rowsA.bottomRows(klen) = tmpA;

          if(res == CBF_RES_OK) {
            res = get_cone_entry_b(data, klen, dst[var].rbeg, dst[var].bbeg, tmpb, bnnz);
            rowsb.bottomRows(klen) = tmpb;
          }

          dstcoord.rbeg = dst[var].rbeg;
          dstcoord.rend = dst[var].rend;
          dstcoord.abeg = dst[var].abeg;
          dstcoord.aend = dst[var].abeg + annz - 1;
          dstcoord.bbeg = dst[var].bbeg;
          dstcoord.bend = dst[var].bbeg + bnnz - 1;
          dst_compressed.push_back(dstcoord);
        }
      }

      if(res == CBF_RES_OK) {

        // Update objective
        dstcoeff = objA.coeff(0, var);
        if(dstcoeff != 0) {
          objA -= (dstcoeff / srccoeff) * srcA;
          objA.coeffRef(0, var) = 0.0;

          data->objbval -= (dstcoeff / srccoeff) * srcb.coeff(0, 0);
        }

        // Update destination rows
        if(dst[var].rend != unset.rend) {
          for(r = it_rowmap->second; r <= it_rowmap->second + klen - 1 && res == CBF_RES_OK; ++r) {
            dstcoeff = rowsA.coeff(r, var);
            if(dstcoeff != 0) {
              rowsA.row(r) -= (dstcoeff / srccoeff) * srcA;
              rowsA.coeffRef(r, var) = 0.0;

              rowsb.coeffRef(r, 0) -= (dstcoeff / srccoeff) * srcb.coeff(0, 0);
            }
          }
        }
      }

      anychange = true;
      mapmaxviol[src[var].rbeg] = std::min(mapmaxviol[src[var].rbeg], 0.0);
    }
  }

  if(anychange && res == CBF_RES_OK) {
    if(changeflag) {
      *changeflag = true;
    }

    // Write objective
    res = put_obj_entry_A(data, &dyn_newdata, objA);

    // Write destination rows
    rowbeg = 0;
    for(std::vector<full_coord_t>::iterator it = dst_compressed.begin(); it != dst_compressed.end(); ++it) {
      klen = it->rend - it->rbeg + 1;

      res = put_cone_entry_A(data, &dyn_newdata, it->rbeg, it->abeg, it->aend, rowsA.middleRows(rowbeg, klen), NULL);

      if(res == CBF_RES_OK) {
        res = put_cone_entry_b(data, &dyn_newdata, it->rbeg, it->bbeg, it->bend, rowsb.middleRows(rowbeg, klen));
      }

      rowbeg += klen;
    }

    if(res == CBF_RES_OK)
      res = CBFdyn_assign(&dyndata, data);

    if(res == CBF_RES_OK)
      res = CBFdyn_append(&dyndata, &newdata);

    if(anyrenaming && res == CBF_RES_OK) {

      // Perform renaming
      if(data->objannz >= 1) {
        for(i = 0; i < data->objannz; ++i) {
          var = data->objasubj[i];

          if(var != varmap[var]) {
            scal = data->objaval[i];
            lookupvar(data, varmap, varscal, var, scal);

            data->objasubj[i] = var;
            data->objaval[i] = scal;
          }
        }
      }

      if(data->annz >= 1) {
        for(i = 0; i < data->annz; ++i) {
          var = data->asubj[i];

          if(var != varmap[var]) {
            scal = varscal[var] * data->aval[i];
            var = varmap[var];
            lookupvar(data, varmap, varscal, var, scal);

            data->asubj[i] = var;
            data->aval[i] = scal;
          }
        }
      }
    }

    if(res == CBF_RES_OK)
      res = CBF_compress_obj(data);

    if(res == CBF_RES_OK)
      res = CBF_compress_maps(data, mapmaxviol);
  }

  return res;
}

void lookupvar(CBFdata* data,
               std::vector<long long int>& varmap,
               std::vector<double>& varscal,
               long long int& var,
               double& scal)
{
  while(var != varmap[var]) {
    scal *= varscal[var];
    var = varmap[var];
  }
}