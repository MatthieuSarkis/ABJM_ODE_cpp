#pragma once

#include <iostream>
#include <complex>
#include <vector>
#include <numeric>  // to use the accumulate attribute of vectors.
#include <algorithm> // to use max_element

using namespace std;

typedef complex<double> comp;

void initialize(int M, vector<int> N, comp** q)
{
  for(int i=0; i<M; i++)
  {
    for(int j=0; j<N[i]; j++)
    {
      q[i][j] = comp(0.5-j/(4*double(N[i])), j/(400*double(N[i])));
    }
  }
}
