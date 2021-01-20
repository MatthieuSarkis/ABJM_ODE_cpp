#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include "equation.h"
using namespace std;

typedef complex<double> comp;

class Euler
{
public:
  Euler(Equation &f);
  ~Euler();

  void solve(double t0, double tf, int n, comp **q, int M, vector<int> N);

private:
  Equation m_f;
};

Euler::Euler(Equation &f)
: m_f(f)
{

}

Euler::~Euler()
{

}

void Euler::solve(double t0, double tf, int n, comp **q, int M, vector<int> N)
{
  double t = t0;
  double h = (tf-t0)/(n-1);
  comp **qp = new comp* [M];
  for(int i=0; i<M; i++) qp[i] = new comp [N[i]];
  double epsilon = 0.0000001;

  //fstream evol;
  //evol.open("evol_temp.dat", ios::out);

  cout << " Computing..." << endl << endl;

  for(int i=0; i<n; i++)
  {
    //evol << t << " ";
    //for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) evol << q[j][k] << " ";
    //evol << endl;
    comp cp1 = q[0][0];
    comp cp2 = q[int(M/2.0)][int(N[int(M/2.0)]/2.0)];
    comp cp3 = q[int(2.0*M/3.0)][int(N[int(2.0*M/3.0)]/2.0)];

    m_f(t, q, qp);
    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) q[j][k] += h*qp[j][k];
    t += h;

    if((real(q[0][0])-real(cp1) < epsilon) && (imag(q[0][0])-imag(cp1) < epsilon)
     && (real(q[int(M/3.0)][int(N[int(M/3.0)]/2.0)])-real(cp2) < epsilon)
     && (imag(q[int(M/3.0)][int(N[int(M/3.0)]/2.0)])-imag(cp2) < epsilon)
     && (real(q[int(2.0*M/3.0)][int(N[int(2.0*M/3.0)]/2.0)])-real(cp3) < epsilon)
     && (imag(q[int(2.0*M/3.0)][int(N[int(2.0*M/3.0)]/2.0)])-imag(cp3) < epsilon)) break;
  }

  //evol.close();

  for(int i=0; i<M; i++) delete[] qp[i];
  delete[] qp;
}
