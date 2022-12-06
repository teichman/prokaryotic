#include <prokaryotic.h>
#include <iostream>

using std::cout, std::endl;
using Eigen::ArrayXd, Eigen::VectorXd;

int main(int argc, char** argv)
{
  Prokaryotic pro;
  pro.initialize(YAML::LoadFile(argv[1]));
  pro.run();
  return 0;
}
