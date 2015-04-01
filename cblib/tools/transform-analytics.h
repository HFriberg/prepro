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

#ifndef CBF_TRANSFORM_ANALYTICS_H
#define CBF_TRANSFORM_ANALYTICS_H

#include "transform.h"
#include <Eigen/Dense>
#include <vector>

CBFresponsee init_varbound(const CBFdata* data, Eigen::VectorXd& lb, Eigen::VectorXd& ub);

CBFresponsee update_varbound_full(const CBFdata* data,
                                  char* intvars,
                                  bool only_intvars,
                                  Eigen::VectorXd& lb,
                                  Eigen::VectorXd& ub,
                                  std::vector<bool>* stronglb,
                                  std::vector<bool>* strongub,
                                  long long int* fixvar,
                                  bool* infeas);

CBFresponsee update_varbound_fast(const CBFdata* data,
                                  CBFtransform_param& param,
                                  char* intvars,
                                  Eigen::VectorXd& lb,
                                  Eigen::VectorXd& ub,
                                  std::vector<bool>* stronglb,
                                  std::vector<bool>* strongub,
                                  double* mapmaxviol,
                                  long long int* fixvar,
                                  bool* infeas,
                                  bool* conestrictredund);

#endif
