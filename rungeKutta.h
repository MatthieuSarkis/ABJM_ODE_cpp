#pragma once

#include <fstream>
#include <vector>
#include <complex>
#include "equation.h"
using namespace std;

typedef complex<double> comp;

class RungeKutta
{
public:
  RungeKutta(Equation &f);
  ~RungeKutta();

  void solve(double t0, double tf, int n, comp **q, int M, vector<int> N);

private:
  Equation m_f;
};


RungeKutta::RungeKutta(Equation &f)
: m_f(f)
{

}

RungeKutta::~RungeKutta()
{

}

void RungeKutta::solve(double t0, double tf, int n, comp **q, int M, vector<int> N)
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
    //for(int j=0; j<M; j++) evol << q[j] << " ";
    //evol << endl;

    comp cp1 = q[0][0];
    comp cp2 = q[int(M/2.0)][int(N[int(M/2.0)]/2.0)];
    comp cp3 = q[int(2.0*M/3.0)][int(N[int(2.0*M/3.0)]/2.0)];

    comp **k1 = new comp* [M];
    comp **k2 = new comp* [M];
    comp **k3 = new comp* [M];
    comp **k4 = new comp* [M];
    comp **temp = new comp* [M];
    for(int i=0; i<M; i++)
    {
      k1[i] = new comp [N[i]];
      k2[i] = new comp [N[i]];
      k3[i] = new comp [N[i]];
      k4[i] = new comp [N[i]];
      temp[i] = new comp [N[i]];
    }

    m_f(t, q, qp);
    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) {k1[j][k]=h*qp[j][k]; temp[j][k]=q[j][k]+k1[j][k]/2.0;}
    m_f(t+h/2, temp, qp);
    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) {k2[j][k]=h*qp[j][k]; temp[j][k]=q[j][k]+k2[j][k]/2.0;}
    m_f(t+h/2, temp, qp);
    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) {k3[j][k]=h*qp[j][k]; temp[j][k]=q[j][k]+k3[j][k];}
    m_f(t+h, temp, qp);
    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) {k4[j][k]=h*qp[j][k];}

    for(int j=0; j<M; j++) for(int k=0; k<N[j]; k++) q[j][k] += (k1[j][k]+2.0*k2[j][k]+2.0*k3[j][k]+k4[j][k])/6.0;


    for(int i=0; i<M; i++)
    {
      delete[] k1[i];
      delete[] k2[i];
      delete[] k3[i];
      delete[] k4[i];
      delete[] temp[i];
    }
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] temp;

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
