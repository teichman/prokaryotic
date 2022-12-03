#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include <prokaryotic.h>

using std::cout, std::endl;
using Eigen::ArrayXd, Eigen::VectorXd;

TEST_CASE("aoeau") {
  Prokaryotic pro;
  pro.initializeHardcoded();

  MoleculeVals cytosol_contents(pro);
  cytosol_contents[0] = 1e8;
  cytosol_contents[1] = 1;
  double um3 = 1.0;
  MoleculeVals concentrations = Cell::cytosolConcentrations(pro, cytosol_contents, um3);
  cout << "Counts: " << endl;
  cout << cytosol_contents.str("  ") << endl;
  cout << "Concentrations: " << endl;
  cout << concentrations.str("  ") << endl;

  cout << "Round trip test" << endl;
  MoleculeVals contents2 = Cell::cytosolContents(pro, concentrations, um3);
  cout << "Counts: " << endl;
  cout << contents2.str("  ") << endl;
  
  CHECK(1 == 1);
  CHECK(cytosol_contents.size() == contents2.size());

  // VectorXd difference = (cytosol_contents.vals_ - contents2.vals_).matrix().norm();
  // CHECK(difference.norm() == 0);

  CHECK((cytosol_contents.vals_ - contents2.vals_).matrix().norm() == 0);
  
  // for (int i = 0; i < cytosol_contents.size(); ++i)
  //   CHECK(cytosol_contents.vals_[i] == contents2.vals_[i]);
}
