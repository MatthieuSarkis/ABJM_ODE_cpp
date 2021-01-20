#pragma once

#include <iostream>
#include <tuple>
#include <vector>
#include <complex>
using namespace std;

typedef complex<double> comp;

tuple<int, vector<int>, vector<int>, vector<comp>> askUser()
{
  int M(0);    // Nombre de noeuds dans le collier.
  cout << endl << "Number of nodes in the quiver: ";
  cin >> M;
  cout << endl;

  vector<int> N;   // Rang des differents noeuds.

  for(int m=0; m<M; m++)
  {
    int temp(0);
    cout << "Rank of node number " << m+1 << ": ";
    cin >> temp;
    N.push_back(temp);
  }

  cout << endl;

  vector<int> K;   // Niveau des differents noeuds.

  for(int m=0; m<M; m++)
  {
    int temp(0);
    cout << "Level of node number " << m+1 << ": ";
    cin >> temp;
    K.push_back(temp);
  }

  cout << endl;

  vector<comp> P;   // Masses.  Enter the compex masses in the form (Re m,Im m).

  for(int m=0; m<M; m++)
  {
    comp temp(0);
    cout << "Mass for bifundamental " << m+1 << ": ";
    cin >> temp;
    P.push_back(temp);
  }

  cout << endl;

  return make_tuple(M, N, K, P);
}
