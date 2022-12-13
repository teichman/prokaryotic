#include <prokaryotic.h>
#include <iostream>

using std::cout, std::endl;
using Eigen::ArrayXd, Eigen::VectorXd;

int main(int argc, char** argv)
{
  Prokaryotic pro;

  const YAML::Node& yaml = YAML::LoadFile("2022-12-08-test_config.yaml");
  pro.applyConfig(yaml);
  Cell::Ptr cell = pro.cells_[0];
  
  // The cell is taking in things in the environment at some rate.  Nevermind how, for the moment.
  cell->membrane_permeabilities_["Amino acids"] = 0.1;
  cell->membrane_permeabilities_["Starch"] = 0.01;
  cell->membrane_permeabilities_["Phosphate"] = 0.01;

  // Start us off with some of each important molecule.
  cell->cytosol_contents_["Amino acids"] = 6e7;  // 100 mM sounds typical
  cell->cytosol_contents_["ATP"] = 6e6;  // 10 mM sounds typical, maybe a bit high
  cell->cytosol_contents_["Phosphate"] = 6e6;
  cell->cytosol_contents_["Starch"] = 6e6;
  cell->cytosol_contents_["Glucose"] = 6e6;
  
  // Dunno how much of these we need, we'll see where they end up at steady state.
  cell->cytosol_contents_["ATP Synthase"] = 1e7;
  cell->cytosol_contents_["Amylase"] = 1e3;
  cell->cytosol_contents_["Ribosome"] = 2e5;
  cell->cytosol_contents_["Proteasome"] = 3e4;
  
  for (const YAML::Node& dnaif : yaml["DNA"])
    cell->addDNAIf(dnaif);

  pro.run();
  return 0;
}
