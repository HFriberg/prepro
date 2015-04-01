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

#include "solver-mosek-x.h"
#include "solver.h"
#include "cbf-data.h"
#include "programmingstyle.h"

#include <mosek.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <math.h>

static CBFresponsee
write(CBFdata *data, CBFsolvermemory *mem);

static void
clean(CBFsolvermemory *mem);

static CBFresponsee
savetodisk(CBFsolvermemory *mem, const char *file);

static CBFresponsee
optimize(CBFsolvermemory *mem);

static CBFresponsee
getprimalobj(CBFsolvermemory *mem, double *val);

static CBFresponsee
getvarslice(CBFsolvermemory *mem, long long int start, long long int len, double *val);

static CBFresponsee
set_objsense(CBFsolvermemory *mem, CBFobjsensee objsense);

static CBFresponsee
set_cj(CBFsolvermemory *mem, long long int j, double cj);

static CBFresponsee
set_lbj(CBFsolvermemory *mem, long long int j, double bound);

static CBFresponsee
set_ubj(CBFsolvermemory *mem, long long int j, double bound);

static void
mosekprint(void *handle, MSKCONST char str[]);

static CBFresponsee
MSKresponse(CBFsolvermemory *mem, MSKrescodee res);

// -------------------------------------
// Global variable
// -------------------------------------

CBFsolver const solver_mosek = { "mosek", write, clean, savetodisk, optimize, getprimalobj, getvarslice, set_objsense, set_cj, set_lbj, set_ubj };

static char response[MSK_MAX_STR_LEN] = "";

// -------------------------------------
// Function definitions
// -------------------------------------

static CBFresponsee write(CBFdata *data, CBFsolvermemory *mem) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;
  long long int i, k, km;

  // Make CBFsolvermemory_mosek if not available
  if (*mem == NULL) {
    *mem = calloc(1, sizeof(CBFsolvermemory_mosek));
    if (!*mem) {
      return CBF_RES_ERR;
    }
  }

  // Make MOSEK environment if not available
  mosek = (CBFsolvermemory_mosek*) *mem;

  if (res == MSK_RES_OK && mosek->env == NULL) {
    res = MSK_makeenv(&mosek->env, "");
  }

  // Delete old task
  if (res == MSK_RES_OK && mosek->task != NULL)
    res = MSK_deletetask(&mosek->task);

  // Create a new task
  if (res == MSK_RES_OK) {
    res = MSK_maketask(mosek->env, 0, 0, &mosek->task);

    if (res == MSK_RES_OK)
      res = MSK_linkfunctotaskstream(mosek->task, MSK_STREAM_ERR, NULL, mosekprint);

    if (res == MSK_RES_OK)
      res = MSK_linkfunctotaskstream(mosek->task, MSK_STREAM_WRN, NULL, mosekprint);

//    if (res == MSK_RES_OK)
//      res = MSK_linkfunctotaskstream(mosek->task, MSK_STREAM_MSG, NULL, mosekprint);
  }

  // Append row, columns and non-zeros
  if (res == MSK_RES_OK)
    res = MSK_appendcons(mosek->task, data->mapnum);

  if (res == MSK_RES_OK)
    res = MSK_appendvars(mosek->task, data->varnum + data->mapnum);

  if (res == MSK_RES_OK)
    res = MSK_putmaxnumanz(mosek->task, data->annz + data->mapnum);

  // Objective sense
  if (res == MSK_RES_OK) {
    if (data->objsense == CBF_OBJ_MAXIMIZE)
      res = MSK_putobjsense(mosek->task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
    else
      res = MSK_putobjsense(mosek->task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  }

  // Objective (cx + c0)
  if (res == MSK_RES_OK)
    res = MSK_putcfix(mosek->task, data->objbval);

  for (i = 0; i < data->objannz && res == MSK_RES_OK; ++i)
    res = MSK_putcj(mosek->task, data->objasubj[i], data->objaval[i]);

  // Constraints (Ax - s = -b)
  for (i = 0; i < data->annz && res == MSK_RES_OK; ++i)
    res = MSK_putaij(mosek->task, data->asubi[i], data->asubj[i], data->aval[i]);

  for (i = 0; i < data->mapnum && res == MSK_RES_OK; ++i)
    res = MSK_putaij(mosek->task, i, data->varnum + i, -1);

  for (i = 0; i < data->mapnum && res == MSK_RES_OK; ++i)
    res = MSK_putbound(mosek->task, MSK_ACC_CON, i, MSK_BK_FX, 0, 0);

  for (i = 0; i < data->bnnz && res == MSK_RES_OK; ++i)
    res = MSK_putbound(mosek->task, MSK_ACC_CON, data->bsubi[i], MSK_BK_FX, -data->bval[i], -data->bval[i]);

  // Integer variables
  for (i = 0; i < data->intvarnum && res == MSK_RES_OK; ++i)
    res = MSK_putvartype(mosek->task, data->intvar[i], MSK_VAR_TYPE_INT);

  // Variable bounds and conic domains
  static MSKboundkeye CBF2MSKlin_cone[CBF_CONE_END];
  CBF2MSKlin_cone[CBF_CONE_FREE] = MSK_BK_FR;
  CBF2MSKlin_cone[CBF_CONE_POS] = MSK_BK_LO;
  CBF2MSKlin_cone[CBF_CONE_NEG] = MSK_BK_UP;
  CBF2MSKlin_cone[CBF_CONE_ZERO] = MSK_BK_FX;

  static MSKconetypee CBF2MSKsoc_cone[CBF_CONE_END];
  CBF2MSKsoc_cone[CBF_CONE_QUAD] = MSK_CT_QUAD;
  CBF2MSKsoc_cone[CBF_CONE_RQUAD] = MSK_CT_RQUAD;

  i = 0;
  for (k = 0; k < data->varstacknum && res == MSK_RES_OK; ++k) {
    switch (data->varstackdomain[k]) {
    case CBF_CONE_FREE:
    case CBF_CONE_POS:
    case CBF_CONE_NEG:
    case CBF_CONE_ZERO:
      for (km = 0; km < data->varstackdim[k] && res == MSK_RES_OK; ++km) {
        res = MSK_putbound(mosek->task, MSK_ACC_VAR, i, CBF2MSKlin_cone[data->varstackdomain[k]], 0.0, 0.0);
        ++i;
      }
      break;

    case CBF_CONE_QUAD:
    case CBF_CONE_RQUAD:
      res = MSK_appendconeseq(mosek->task, CBF2MSKsoc_cone[data->varstackdomain[k]], 0.0, data->varstackdim[k], i);

      for (km = 0; km < data->varstackdim[k] && res == MSK_RES_OK; ++km) {
        res = MSK_putbound(mosek->task, MSK_ACC_VAR, i, MSK_BK_FR, 0.0, 0.0);
        ++i;
      }
      break;

    default:
      return CBF_RES_ERR;
    }
  }

  for (k = 0; k < data->mapstacknum && res == MSK_RES_OK; ++k) {
    switch (data->mapstackdomain[k]) {
    case CBF_CONE_FREE:
    case CBF_CONE_POS:
    case CBF_CONE_NEG:
    case CBF_CONE_ZERO:
      for (km = 0; km < data->mapstackdim[k] && res == MSK_RES_OK; ++km) {
        res = MSK_putbound(mosek->task, MSK_ACC_VAR, i, CBF2MSKlin_cone[data->mapstackdomain[k]], 0.0, 0.0);
        ++i;
      }
      break;

    case CBF_CONE_QUAD:
    case CBF_CONE_RQUAD:
      res = MSK_appendconeseq(mosek->task, CBF2MSKsoc_cone[data->mapstackdomain[k]], 0.0, data->mapstackdim[k], i);

      for (km = 0; km < data->mapstackdim[k] && res == MSK_RES_OK; ++km) {
        res = MSK_putbound(mosek->task, MSK_ACC_VAR, i, MSK_BK_FR, 0.0, 0.0);
        ++i;
      }
      break;

    default:
      return CBF_RES_ERR;
    }
  }

  return MSKresponse(mem, res);
}

static void clean(CBFsolvermemory *mem) {

  CBFsolvermemory_mosek *mosek;

  if (*mem) {
    // Free memory allocated by MOSEK
    mosek = (CBFsolvermemory_mosek*) *mem;

    if (mosek->task)
      MSK_deletetask(&mosek->task);

    if (mosek->env)
      MSK_deleteenv(&mosek->env);

    free(*mem);
  }
}

static CBFresponsee savetodisk(CBFsolvermemory *mem, const char *file) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;

  if (!*mem)
    return CBF_RES_ERR;

  mosek = (CBFsolvermemory_mosek*) *mem;
  res = MSK_writedata(mosek->task, file);

  return MSKresponse(mem, res);
}

static CBFresponsee optimize(CBFsolvermemory *mem) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;

  if (!*mem)
    return CBF_RES_ERR;

  // Call MOSEK optimizer
  mosek = (CBFsolvermemory_mosek*) *mem;

//  if (res == MSK_RES_OK)
//    res = MSK_writedata(mosek->task, "tmp.opf");

  if (res == MSK_RES_OK)
    res = MSK_optimize(mosek->task);

//  if (res == MSK_RES_OK)
//    res = MSK_solutionsummary(mosek->task, MSK_STREAM_WRN);

  return MSKresponse(mem, res);
}

static CBFresponsee getprimalobj(CBFsolvermemory *mem, double *val) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;
  MSKsoltypee soltype;
  MSKobjsensee objsense;
  MSKsolstae solsta;
  MSKbooleant isdef;

  if (!*mem)
    return CBF_RES_ERR;

  // Get primal objective from MOSEK
  mosek = (CBFsolvermemory_mosek*) *mem;

  soltype = MSK_SOL_ITG;
  res = MSK_solutiondef(mosek->task, soltype, &isdef);
  if (res == MSK_RES_OK && !isdef) {
    soltype = MSK_SOL_ITR;
    res = MSK_solutiondef(mosek->task, soltype, &isdef);
    if (res == MSK_RES_OK && !isdef)
      return CBF_RES_ERR;
  }

  if (res == MSK_RES_OK)
    res = MSK_solutiondef(mosek->task, soltype, &isdef);

  if (isdef && res == MSK_RES_OK)
    res = MSK_getsolsta(mosek->task, soltype, &solsta);
  else
    return CBF_RES_ERR;

  if (res == MSK_RES_OK) {
    switch (solsta) {
    case MSK_SOL_STA_DUAL_INFEAS_CER:
      res = MSK_getobjsense(mosek->task, &objsense);

      if (res == MSK_RES_OK) {
        if (objsense == MSK_OBJECTIVE_SENSE_MAXIMIZE)
          *val = INFINITY;
        else
          *val = -INFINITY;
      }
      break;

    case MSK_SOL_STA_PRIM_INFEAS_CER:
      res = MSK_getobjsense(mosek->task, &objsense);

      if (res == MSK_RES_OK) {
        if (objsense == MSK_OBJECTIVE_SENSE_MAXIMIZE)
          *val = -INFINITY;
        else
          *val = INFINITY;
      }
      break;

    case MSK_SOL_STA_INTEGER_OPTIMAL:
    case MSK_SOL_STA_OPTIMAL:
      res = MSK_getprimalobj(mosek->task, soltype, val);
      break;

    default:
      return CBF_RES_ERR;
      break;
    }
  }

  return MSKresponse(mem, res);
}

static CBFresponsee getvarslice(CBFsolvermemory *mem, long long int first, long long int len, double *val) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;
  MSKsoltypee soltype;
  MSKbooleant isdef;

  if (!*mem)
    return CBF_RES_ERR;

  // Get primal solution from MOSEK
  mosek = (CBFsolvermemory_mosek*) *mem;

  soltype = MSK_SOL_ITG;
  res = MSK_solutiondef(mosek->task, soltype, &isdef);
  if (res == MSK_RES_OK && !isdef) {
    soltype = MSK_SOL_ITR;
    res = MSK_solutiondef(mosek->task, soltype, &isdef);
    if (res == MSK_RES_OK && !isdef)
      return CBF_RES_ERR;
  }

  if (res == MSK_RES_OK)
    res = MSK_getxxslice(mosek->task, soltype, first, first + len, val);

  return MSKresponse(mem, res);
}

static CBFresponsee set_objsense(CBFsolvermemory *mem, CBFobjsensee objsense) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;

  if (!*mem)
    return CBF_RES_ERR;

  // Update objective sense
  mosek = (CBFsolvermemory_mosek*) *mem;
  if (objsense == CBF_OBJ_MAXIMIZE)
    res = MSK_putobjsense(mosek->task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
  else
    res = MSK_putobjsense(mosek->task, MSK_OBJECTIVE_SENSE_MINIMIZE);

  return MSKresponse(mem, res);
}

static CBFresponsee set_cj(CBFsolvermemory *mem, long long int j, double cj) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;

  if (!*mem)
    return CBF_RES_ERR;

  // Update objective coefficient
  mosek = (CBFsolvermemory_mosek*) *mem;
  res = MSK_putcj(mosek->task, j, cj);

  return MSKresponse(mem, res);
}

static CBFresponsee set_lbj(CBFsolvermemory *mem, long long int j, double bound) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;
  MSKboundkeye bk;
  MSKrealt bl, bu;

  if (!*mem)
    return CBF_RES_ERR;

  // Update bound value
  mosek = (CBFsolvermemory_mosek*) *mem;
  res = MSK_getvarbound(mosek->task, j, &bk, &bl, &bu);

  if (res == MSK_RES_OK) {
    switch (bk) {
    case MSK_BK_FR:
      bk = MSK_BK_LO;
      break;
    case MSK_BK_UP:
      bk = MSK_BK_RA;
      break;
    case MSK_BK_FX:
      return CBF_RES_ERR;
    default:
      break;
    }

    res = MSK_putvarbound(mosek->task, j, bk, bound, bu);
  }

  return MSKresponse(mem, res);
}

static CBFresponsee set_ubj(CBFsolvermemory *mem, long long int j, double bound) {
  MSKrescodee res = MSK_RES_OK;
  CBFsolvermemory_mosek *mosek;
  MSKboundkeye bk;
  MSKrealt bl, bu;

  if (!*mem)
    return CBF_RES_ERR;

  // Update bound value
  mosek = (CBFsolvermemory_mosek*) *mem;
  res = MSK_getvarbound(mosek->task, j, &bk, &bl, &bu);

  if (res == MSK_RES_OK) {
    switch (bk) {
    case MSK_BK_FR:
      bk = MSK_BK_UP;
      break;
    case MSK_BK_LO:
      bk = MSK_BK_RA;
      break;
    case MSK_BK_FX:
      return CBF_RES_ERR;
    default:
      break;
    }

    res = MSK_putvarbound(mosek->task, j, bk, bl, bound);
  }

  return MSKresponse(mem, res);
}

static void mosekprint(void *handle, MSKCONST char str[]) {
  printf("%s", str);
}

static CBFresponsee MSKresponse(CBFsolvermemory *mem, MSKrescodee res) {
  switch (res) {
  case MSK_RES_OK:
    return CBF_RES_OK;

  default:
    MSK_getcodedesc(res, NULL, response);
    printf("%s\n", response);
    return CBF_RES_ERR;
  }
}
