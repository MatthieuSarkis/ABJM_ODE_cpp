#pragma once

#include <fstream>
#include <vector>
#include <complex>
#include "equation.h"
using namespace std;

typedef complex<double> comp;

class PredictorCorrector
{
public:
  PredictorCorrector(Equation &f);
  ~PredictorCorrector();

  void solve(double t0, double tf, int n, comp **q, int M, vector<int> N);

private:
  Equation m_f;
};

PredictorCorrector::PredictorCorrector(Equation &f)
: m_f(f)
{

}

PredictorCorrector::~PredictorCorrector()
{

}

void PredictorCorrector::solve(double t0, double tf, int n, comp **q, int M, vector<int> N)
{
  double t = t0;
  double h = (tf-t0)/(n-1);
  comp **qp = new comp* [M]; for(int i=0; i<M; i++) qp[i] = new comp [N[i]];
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
    comp cp2 = q[int(M/3.0)][int(N[int(M/3.0)]/2.0)];
    comp cp3 = q[int(2.0*M/3.0)][int(N[int(2.0*M/3.0)]/2.0)];

    m_f(t, q, qp);
    comp **temp = new comp* [M]; for(int j=0; j<M; j++) temp[j] = new comp [N[j]];
    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) temp[j][k] = q[j][k]+h*qp[j][k];
    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) q[j][k] += (h/2)*qp[j][k];
    m_f(t+h, temp, qp);
    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) q[j][k] += (h/2)*qp[j][k];
    t += h;
    for(int j=0; j<M; j++) delete[] temp[j];
    delete[] temp;

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
