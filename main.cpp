#include "euler.h"
#include "predictorCorrector.h"
#include "rungeKutta.h"
#include "equation.h"
#include "askUser.h"
#include "gnuplot.h"
#include "initialize.h"
#include "writeInFile.h"
#include <tuple>
#include <cstdlib>
#include <complex>
#include <fstream>

typedef complex<double> comp;

int main()
{
  int M;
  vector<int> N;
  vector<int> K;
  vector<comp> P;

  tie(M, N, K, P) = askUser();

  comp **q = new comp* [M];
  for(int i=0; i<M; i++) q[i] = new comp [N[i]];
  initialize(M, N, q);
  
  Equation eq(M, N, K, P);

  //Euler solver(eq);
  PredictorCorrector solver(eq);
  //RungeKutta solver(eq);

  solver.solve(0, 1000, 100000, q, M, N);

  string fichier = writeInFile(M, N, K, P, q);

  gnuplot p;
  p("set term postscript eps");
  p("set output \"" + fichier + ".eps\" ");
  p("plot \'./" + fichier + ".dat\'");

  for(int i=0; i<M; i++) delete[] q[i];
  delete[] q;

  return 0;
}
