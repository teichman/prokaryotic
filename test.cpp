#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include <prokaryotic.h>

using std::cout, std::endl;
using Eigen::ArrayXd, Eigen::VectorXd;

TEST_CASE("Cytosol concentrations / contents round trip") {
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

TEST_CASE("Membrane permeability")
{
  Prokaryotic pro;
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Phosphate", "P", 94.97)));
  
  Cell::Ptr cell(new Cell(pro, "cell"));
  cell->membrane_permeabilities_["Phosphate"] = 0.1;
  pro.cells_.push_back(cell);

  Biome::Ptr biome(new Biome(pro, 10, "Alkaline vents"));
  biome->concentrations_["Phosphate"] = 10;
  pro.biomes_.push_back(biome);

  cout << "Before" << endl;
  cout << biome->str() << endl;
  cout << cell->str() << endl;

  for (int i = 0; i < 1000; ++i)
    pro.tick();

  cout << "------------------------------" << endl;
  cout << "After" << endl;
  cout << biome->str() << endl;
  cout << cell->str() << endl;
  cout << "Cytosol concentrations: " << endl;
  cout << cell->cytosolConcentrations().str("  ") << endl;
  
  CHECK(cell->cytosolConcentrations()["Phosphate"] == doctest::Approx(biome->concentrations_["Phosphate"]));

  // If the cell magically shrinks, the concentration goes up.
  // This maybe actually isn't the right behavior but it's by design for now.
  cell->um3_ *= 0.5;
  CHECK(cell->cytosolConcentrations()["Phosphate"] == doctest::Approx(2.0 * biome->concentrations_["Phosphate"]));

  for (int i = 0; i < 1000; ++i)
    pro.tick();

  // Confirm it re-equalized.
  CHECK(cell->cytosolConcentrations()["Phosphate"] == doctest::Approx(biome->concentrations_["Phosphate"]));
}
