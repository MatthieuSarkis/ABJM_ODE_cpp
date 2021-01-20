#pragma once

#include <vector>
#include <complex>
#include "mathFunction.h"
using namespace std;

typedef complex<double> comp;
const comp I(0.0,1.0);

class Equation : public MathFunction<comp**>
{
public:
  Equation(int M, vector<int> N, vector<int> K, vector<comp> P);
  ~Equation();

  void operator()(double t, comp **q, comp **qp);

private:
  int m_M;
  vector<int> m_N;
  vector<int> m_K;
  vector<comp> m_P;
};

Equation::Equation(int M, vector<int> N, vector<int> K, vector<comp> P)
: m_M(M), m_N(N), m_K(K), m_P(P)
{

}

Equation::~Equation()
{

}

int ind(int k, int M)
{
  if(k == -1)
  {
    return M-1;
  }

  else if(k == M)
  {
    return 0;
  }

  else
  {
    return k;
  }
}

void Equation::operator()(double t, comp** q, comp** qp)
{
  comp C, T1, T2;

  for(int i=0; i<m_M; i++)
  {
    for(int j=0; j<m_N[i]; j++)
    {
        C = comp(0.0, 0.0);
        T1 = comp(0.0, 0.0);
        T2 = comp(0.0, 0.0);

      for(int k=0; k<m_N[i]; k++)
      {
        if(k != j)
        {
          C += cosh((q[i][k]-q[i][j])/2.)/sinh((q[i][k]-q[i][j])/2.);
        }
      }

      for(int k=0; k<m_N[ind(i+1, m_M)]; k++)
      {
        T1 += tanh((q[ind(i+1, m_M)][k]-q[i][j]+m_P[i])/2.);
      }

      for(int k=0; k<m_N[ind(i-1, m_M)]; k++)
      {
        T2 += tanh((q[ind(i-1, m_M)][k]-q[i][j]-m_P[ind(i-1, m_M)])/2.);
      }

      qp[i][j] = (I*double(m_K[i])*q[i][j])/(2*M_PI)-C+0.5*T1+0.5*T2;
    }
  }
}
