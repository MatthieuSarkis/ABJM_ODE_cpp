#include <fstream>
#include <string>

using namespace std;

string writeInFile(int M, vector<int> N, vector<int> K, vector<comp> P, comp **q)
{
  fstream res;

  string fichier;
  stringstream ss;
  ss << M << "nodes";
  for(int m=0; m<M; m++)
  {
    ss << "_N" << m+1 << "=" << N[m];
  }
  for(int m=0; m<M; m++)
  {
    ss << "_k" << m+1 << "=" << K[m];
  }
  for(int m=0; m<M; m++)
  {
    ss << "_m" << m+1 << "=" << P[m];
  }

  fichier = ss.str();

  if(fichier.length() > 200)
  {
    stringstream ss_alt;
    ss_alt << M << "nodes" << "_N=" << N[0];

    fichier = ss_alt.str();
  }

  res.open((fichier+".dat").c_str(), ios::out);   // On ecrit les parties reelles et imaginaires
                                             // des valeurs propres asymptotiques dans un fichier.

  for(int i=0; i<M; i++) for(int j=0; j<N[i]; j++) res << real(q[i][j]) << " " << imag(q[i][j]) << endl;

  res.close();

  return fichier;
}
