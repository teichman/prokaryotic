#include <prokaryotic.h>
#include <iostream>

using std::cout, std::endl;
using Eigen::ArrayXd, Eigen::VectorXd;

int main(int argc, char** argv)
{
  Prokaryotic pro;
  pro.initialize(YAML::LoadFile(argv[1]));
  //pro.initializeHardcoded();
  // MoleculeVals cytosol_contents(pro);
  // cytosol_contents[0] = 1e8;
  // cout << cytosol_contents.str() << endl;  
  pro.run();
  return 0;
}
