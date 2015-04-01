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

#include "transform-varstack-to-mapstack.h"
#include "transform-helper.h"
#include "cbf-format.h"

#include "stdio.h"

static CBFresponsee
transform(CBFdata *data, CBFtransform_param &param, bool *changeflag);

static CBFresponsee
transform_callback(CBFtransform_param param, CBFtransformmemory* mem, long long int &k, long long int &varbeg, CBFdata *data, CBFdyndata& newdata);

// -------------------------------------
// Global variables
// -------------------------------------

CBFtransform const transform_varstack_to_mapstack = { "varstack-to-mapstack", transform };

static const bool lincones = true;
static const bool socones = true;

// -------------------------------------
// Function definitions
// -------------------------------------

static CBFresponsee transform(CBFdata *data, CBFtransform_param &param, bool *changeflag) {
  CBFresponsee res = CBF_RES_OK;

  if (changeflag) {
    printf("changeflag not supported by this transformation");
    return CBF_RES_ERR;
  }

  // Move cones from varstack to mapstack
  if (res == CBF_RES_OK)
    res = transform_stackwise(data, param, NULL, transform_callback, NULL);

  return res;
}

static CBFresponsee transform_callback(CBFtransform_param param, CBFtransformmemory* mem, long long int &k, long long int &varbeg, CBFdata *data, CBFdyndata& newdata) {

  CBFresponsee res = CBF_RES_OK;
  long long int j, mapnum;

  if ((lincones && (data->varstackdomain[k] == CBF_CONE_ZERO || data->varstackdomain[k] == CBF_CONE_POS || data->varstackdomain[k] == CBF_CONE_NEG))
      || (socones && (data->varstackdomain[k] == CBF_CONE_QUAD || data->varstackdomain[k] == CBF_CONE_RQUAD))) {

    mapnum = data->mapnum + newdata.data->mapnum;

    res = CBFdyn_map_capacitysurplus(&newdata, 1);

    if (res == CBF_RES_OK) {
      res = CBFdyn_map_adddomain(&newdata, data->varstackdomain[k], data->varstackdim[k]);
    }

    if (res == CBF_RES_OK) {
      res = CBFdyn_a_capacitysurplus(&newdata, data->varstackdim[k]);
    }

    for (j = 0; j < data->varstackdim[k] && res == CBF_RES_OK; ++j) {
      res = CBFdyn_a_add(&newdata, mapnum + j, varbeg + j, 1.0);
    }

    if (res == CBF_RES_OK) {
      data->varstackdomain[k] = CBF_CONE_FREE;
    }
  }

  return res;
}
