#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace std;


using Eigen::MatrixXd;
using Eigen::EigenSolver;

int main()
{
  double testeNum;
  MatrixXd m(2,2);
  m(0,0) = -1;
  m(1,0) = 1;
  m(0,1) = -1;
  m(1,1) = 8;
  cout << "m: " << endl << m << endl;
  EigenSolver<MatrixXd> es(m);
  cout << "AutoValores de m:" << endl << es.eigenvalues() << endl;
  cout << "Matriz de autovetores de m:" << endl << es.eigenvectors() << endl;
  testeNum = es.eigenvalues()[0].real();
  cout << testeNum << endl;
  return 0;
}
