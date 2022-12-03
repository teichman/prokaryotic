#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include <prokaryotic.h>

using std::cout, std::endl;
using Eigen::ArrayXd, Eigen::VectorXd;

TEST_CASE("Basics") {
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

  SUBCASE("Cytosol contents / concentration round trip test") {
    cout << "Round trip test" << endl;
    MoleculeVals cytosol_contents2 = Cell::cytosolContents(pro, concentrations, um3);
    cout << "Counts: " << endl;
    cout << cytosol_contents2.str("  ") << endl;
    CHECK(cytosol_contents.size() == cytosol_contents2.size());
    CHECK((cytosol_contents.vals_ - cytosol_contents2.vals_).matrix().norm() == 0);  // array doesn't have norm(), really?
  }
}

TEST_CASE("Half life")
{
  double half_life_hours = 1.0;
  double p_denature = probabilityPerSecond(half_life_hours);
  cout << "half_life_hours: " << half_life_hours << endl;
  cout << "p_denature per tick: " << p_denature << endl;

  double population = 1.0;
  int num_ticks = 0;
  for (int i = 0; i < int(half_life_hours * 60 * 60); ++i) {
    population *= (1.0 - p_denature);
    num_ticks += 1;
  }
  cout << "After " << num_ticks << " ticks, population is now " << population << endl;
  CHECK(population == doctest::Approx(0.5));
}

TEST_CASE("2")
{
  Prokaryotic pro;

  biomes_.push_back(Biome::Ptr(new Biome(*this, 10, "Alkaline vents")));
  biomes_[0]->concentrations_["Phosphate"] = 0.01;
  biomes_[0]->concentrations_["R"] = 0;
  biomes_[0]->concentrations_["X"] = 10;
  
}
