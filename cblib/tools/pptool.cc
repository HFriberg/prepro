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

#include "frontend-cbf.h"
#include "backend-cbf.h"
#include "transform-none.h"
#include "transform-problem-pp.h"
#include "transform-presolve.h"
//#include "transform-soc-q-simple.h"
//#include "transform-soc-q-conservative.h"
//#include "transform-soc-q-full.h"
//#include "transform-varbound-fast.h"
//#include "transform-varbound-full.h"

#include "console.h"

#include <string>
#include <stdio.h>


// -------------------------------------
// Function definitions
// -------------------------------------

int main (int argc, char *argv[])
{
  CBFresponsee res = CBF_RES_OK;
  const CBFfrontend *default_frontend, *frontend;
  const CBFbackend  *default_backend,  *backend;
  const CBFtransform *default_transform, *transform;
  std::string ofile;
  const char *ifile;
  const char *opath;
  const char *pfix;
  int i;

  // For debugging crashes
  setbuf(stdout, NULL);

  // List of plugins
  const CBFfrontend *plugs_frontend[] = {&frontend_cbf,
                                         NULL};

  const CBFbackend  *plugs_backend[]  = {&backend_cbf,
                                         NULL};

  const CBFtransform *plugs_transform[] = {&transform_none,
                                           &transform_problem_pp,
                                           &transform_presolve,
//                                           &transform_varbound_fast,
//                                           &transform_varbound_full,
                                           NULL};

  // Default options
  frontend  = default_frontend  = &frontend_cbf;
  backend   = default_backend   = &backend_cbf;
  transform = default_transform = &transform_none;
  opath = NULL;
  pfix  = NULL;

  // User defined options
  res = getoptions(argc, argv, plugs_frontend, plugs_backend, plugs_transform,
                   &frontend,
                   &backend,
                   &transform,
                   &opath,
                   &pfix);

  if (argc <= 1 || res != CBF_RES_OK)
  {
    printf("\nBad command, syntax is:\n");
    printf(">> pptool [OPTIONS] infile1 infile2 infile3 ...\n\n");
    printoptions(plugs_frontend, plugs_backend, plugs_transform,
            default_frontend, default_backend, default_transform);
  }
  else
  {
    printf("NAME;TIME;STAT\n");

    // All non-nullified arguments are filenames
    for (i=1; i<argc && res==CBF_RES_OK; ++i) {
      if (argv[i]) {
        ifile = argv[i];
        ofile = swapfiledirandext(ifile, opath, pfix, backend->format);

        printf("%s", ifile);
        res = processfile(frontend, backend, transform, ifile, ofile.c_str());
        printf("\n");
      }
    }
  }

  return res;
}
