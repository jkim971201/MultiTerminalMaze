////////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2018, Iowa State University All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software
// without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
////////////////////////////////////////////////////////////////////////////////

#include <vector>

#pragma once

namespace flt 
{

/*****************************/
/*  User-Defined Parameters  */
/*****************************/
#define FLUTE_ACCURACY 10  // Default accuracy
#define FLUTE_ROUTING 1    // 1 to construct routing, 0 to estimate WL only
#define FLUTE_LOCAL_REFINEMENT 1  // Suggestion: Set to 1 if ACCURACY >= 5
#define FLUTE_REMOVE_DUPLICATE_PIN \
  0  // Remove dup. pin for flute_wl() & flute()

#define FLUTE_POWVFILE "POWV9.dat"  // LUT for POWV (Wirelength Vector)
#define FLUTE_POSTFILE "POST9.dat"  // LUT for POST (Steiner Tree)
#define FLUTE_D 9  // LUT is used for d <= FLUTE_D, FLUTE_D <= 9

using DTYPE = int;

struct Branch
{
  DTYPE x, y;  // starting point of the branch
  int n;       // index of neighbor
};

struct Tree
{
  int deg;                     // degree
  DTYPE length;                // total wirelength
  std::vector<Branch> branch;  // array of tree branches
  int branchCount() const { return branch.size(); }
};

struct point
{
  DTYPE x, y;
  int o; // o = order?
};

// User-Callable Functions
// Delete LUT tables for exit so they are not leaked.
void deleteLUT();
DTYPE flute_wl(int d,
               const std::vector<DTYPE>& x,
               const std::vector<DTYPE>& y,
               int acc);
Tree flute(const std::vector<DTYPE>& x, const std::vector<DTYPE>& y, int acc);
DTYPE wirelength(Tree t);
void plottree(Tree t);
void write_svg(Tree t, const char* filename);

// Other useful functions
DTYPE flutes_wl_LD(int d,
                   const std::vector<DTYPE>& xs,
                   const std::vector<DTYPE>& ys,
                   const std::vector<int>& s);
DTYPE flutes_wl_MD(int d,
                   const std::vector<DTYPE>& xs,
                   const std::vector<DTYPE>& ys,
                   const std::vector<int>& s,
                   int acc);
DTYPE flutes_wl_RDP(int d,
                    std::vector<DTYPE> xs,
                    std::vector<DTYPE> ys,
                    std::vector<int> s,
                    int acc);
Tree flutes_LD(int d,
               const std::vector<DTYPE>& xs,
               const std::vector<DTYPE>& ys,
               const std::vector<int>& s);
Tree flutes_MD(int d,
               const std::vector<DTYPE>& xs,
               const std::vector<DTYPE>& ys,
               const std::vector<int>& s,
               int acc);
Tree flutes_RDP(int d,
                std::vector<DTYPE> xs,
                std::vector<DTYPE> ys,
                std::vector<int> s,
                int acc);

inline DTYPE flutes_wl_LMD(int d,
                           const std::vector<DTYPE>& xs,
                           const std::vector<DTYPE>& ys,
                           const std::vector<int>& s,
                           int acc)
{
  if (d <= FLUTE_D) {
    return flutes_wl_LD(d, xs, ys, s);
  }
  return flutes_wl_MD(d, xs, ys, s, acc);
}

inline DTYPE flutes_wl_ALLD(int d,
                            const std::vector<DTYPE>& xs,
                            const std::vector<DTYPE>& ys,
                            const std::vector<int>& s,
                            int acc)
{
  return flutes_wl_LMD(d, xs, ys, s, acc);
}

inline DTYPE flutes_wl(int d,
                       const std::vector<DTYPE>& xs,
                       const std::vector<DTYPE>& ys,
                       const std::vector<int>& s,
                       int acc)
{
  if (FLUTE_REMOVE_DUPLICATE_PIN == 1) {
    return flutes_wl_RDP(d, xs, ys, s, acc);
  }
  return flutes_wl_ALLD(d, xs, ys, s, acc);
}

inline Tree flutes_ALLD(int d,
                        const std::vector<DTYPE>& xs,
                        const std::vector<DTYPE>& ys,
                        const std::vector<int>& s,
                        int acc)
{
  if (d <= FLUTE_D) {
    return flutes_LD(d, xs, ys, s);
  }
  return flutes_MD(d, xs, ys, s, acc);
}

inline Tree flutes(const std::vector<DTYPE>& xs,
                   const std::vector<DTYPE>& ys,
                   const std::vector<int>& s,
                   int acc)
{
  int d = xs.size();
  if (FLUTE_REMOVE_DUPLICATE_PIN == 1) {
    return flutes_RDP(d, xs, ys, s, acc);
  }
  return flutes_ALLD(d, xs, ys, s, acc);
}

inline Tree flutes_LMD(int d,
                       const std::vector<DTYPE>& xs,
                       const std::vector<DTYPE>& ys,
                       const std::vector<int>& s,
                       int acc)
{
  if (d <= FLUTE_D) {
    return flutes_LD(d, xs, ys, s);
  }
  return flutes_MD(d, xs, ys, s, acc);
}

}  // namespace flt
