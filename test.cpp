#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include <prokaryotic.h>

using namespace std;
using Eigen::ArrayXd, Eigen::VectorXd;

TEST_CASE("Cytosol concentrations / contents round trip") {
  Prokaryotic pro;
  //pro.initializeHardcoded();
  pro.initialize(YAML::LoadFile("config.yaml"));

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

TEST_CASE("Ribosome")
{
  Prokaryotic pro;
  double half_life_hours = 0.5;
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, YAML::Load("name: ATP\n"
                                                                         "symbol: bang\n"
                                                                         "daltons: 507.18"))));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, YAML::Load("name: Amino acids\n"
                                                                         "symbol: bricks\n"
                                                                         "daltons: 100"))));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, YAML::Load("name: Phosphate\n"
                                                                         "symbol: P\n"
                                                                         "daltons: 94.97"))));
  
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ATP Synthase", ":hammer:", 5e5, 1, half_life_hours)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Ribosome", ":factory:", 2e6, 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ADP", ":briefcase:", 423.17)));
  
  Cell::Ptr cell(new Cell(pro, "cell"));
  pro.cells_.push_back(cell);
  Biome::Ptr biome(new Biome(pro, 10, "Alkaline vents"));
  pro.biomes_.push_back(biome);

  cout << pro.str() << endl;
  cout << pro.molecule("ATP Synthase")->str() << endl;
  cout << pro.molecule("ATP Synthase")->pDenature() << endl;

  // Confirm that ribosomes don't leak out or degenerate as configured above.
  cell->cytosol_contents_["Ribosome"] = 1000;
  for (int i = 0; i < 1000; ++i)
    pro.tick();
  CHECK(cell->cytosol_contents_["Ribosome"] == 1000);

  // Add a reaction that generates ATP Synthase.
  {
    MoleculeVals inputs(pro);
    inputs["ADP"] = 1;  // A strange world in which we will fabricate ATP Synthase from 1 molecule of ADP.
    MoleculeVals outputs(pro);
    outputs["ATP Synthase"] = 1;
    MoleculeVals kms(pro);
    kms["ADP"] = 1e-3;
    double kcat = 0.1;
    ReactionType::ConstPtr rt(new ReactionType(pro, inputs, outputs, kms, kcat));
    cell->dna_->synthesis_reactions_[pro.moleculeIdx("ATP Synthase")] = rt;
  }
  
  // Run it.  We shouldn't get any yet because we don't have any building blocks to make the ATP Synthase from.
  cout << "----------------------------------------" << endl;
  cout << "Added reaction to make ATP Synthase, but no building blocks here. " << endl;
  for (int i = 0; i < 1000; ++i)
    pro.tick();
  cout << cell->str() << endl;
  CHECK(cell->cytosol_contents_["ATP Synthase"] == 0);

  // Now add the building blocks.  We should get some ATP Synthase.
  cout << "----------------------------------------" << endl;
  cout << "Added building blocks to make ATP Synthase out of. " << endl;
  cell->cytosol_contents_["ADP"] = 300;
  for (int i = 0; i < 100; ++i) {
    pro.tick();
    //cout << cell->str() << endl;
  }
  cout << cell->str() << endl;
  CHECK(cell->cytosol_contents_["ATP Synthase"] > 200);

  cout << "----------------------------------------" << endl;
  cout << "ATP Synthase breaks down and we run out of building blocks. " << endl;
  for (int i = 0; i < 10000; ++i)
    pro.tick();
  cout << cell->str() << endl;
  CHECK(cell->cytosol_contents_["ATP Synthase"] < 10);
}

TEST_CASE("YAML")
{
  Prokaryotic pro;

  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, YAML::Load("name: ATP\n"
                                                                         "symbol: bang\n"
                                                                         "daltons: 507.18"))));
  cout << pro.molecule("ATP")->str() << endl;
  CHECK(pro.molecule("ATP")->name_ == "ATP");
  CHECK(pro.molecule("ATP")->symbol_ == "bang");
  CHECK(pro.molecule("ATP")->daltons_ == 507.18);
  
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ADP", "", 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Phosphate", "", 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "X", "", 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "R", "", 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ATP Synthase", "", 1, 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ATP Consumer", "", 1, 1)));

  { 
    YAML::Node yaml = YAML::Load("formula: ADP + Phosphate -> ATP\n"
                                 "protein: ATP Synthase\n"
                                 "kcat: 1e-3\n"
                                 "KMs:\n"
                                 "  ADP: 1e-1\n"
                                 "  Phosphate: 1e-2");
    ReactionType rt(pro, yaml);
    cout << rt.str() << endl;

    CHECK(rt.inputs_["ADP"] == 1);
    CHECK(rt.inputs_["Phosphate"] == 1);
    CHECK(rt.inputs_.vals_.sum() == 2);
  
    CHECK(rt.outputs_["ATP"] == 1);
    CHECK(rt.outputs_.vals_.sum() == 1);
    CHECK(rt.kcat_ == 1e-3);

    CHECK(rt.kms_["ADP"] == 1e-1);
    CHECK(rt.kms_["Phosphate"] == 1e-2);
    CHECK(rt.kms_.vals_.sum() == rt.kms_["ADP"] + rt.kms_["Phosphate"]);
  }

  { 
    YAML::Node yaml = YAML::Load("formula: 2 X + ATP -> R + ADP + Phosphate\n"
                                 "protein: ATP Consumer\n"
                                 "kcat: 1e-2\n"
                                 "KMs:\n"
                                 "  X: 1e-1\n"
                                 "  ATP: 1e-3\n");
    ReactionType rt(pro, yaml);
    cout << rt.str() << endl;

    CHECK(rt.inputs_["X"] == 2);
    CHECK(rt.inputs_["ATP"] == 1);
    CHECK(rt.inputs_.vals_.sum() == 3);
  
    CHECK(rt.outputs_["R"] == 1);
    CHECK(rt.outputs_["ADP"] == 1);
    CHECK(rt.outputs_["Phosphate"] == 1);
    CHECK(rt.outputs_.vals_.sum() == 3);
    
    CHECK(rt.kcat_ == 1e-2);
    CHECK(rt.kms_["X"] == 1e-1);
    CHECK(rt.kms_["ATP"] == 1e-3);
    CHECK(rt.kms_.vals_.sum() == rt.kms_["X"] + rt.kms_["ATP"]);
  }
}

TEST_CASE("DNAIf")
{
  Prokaryotic pro;
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Ribosome", "", 1, 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ATP Synthase", "", 1, 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ADP", ":briefcase:", 423.17)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Amino acids", "", 100)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, YAML::Load("name: ATP\n"
                                                                         "symbol: bang\n"
                                                                         "daltons: 507.18"))));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, YAML::Load("name: Amino acids\n"
                                                                         "symbol: bricks\n"
                                                                         "daltons: 100"))));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, YAML::Load("name: Phosphate\n"
                                                                         "symbol: P\n"
                                                                         "daltons: 94.97"))));
  
  Cell::Ptr cell(new Cell(pro, "cell"));
  pro.cells_.push_back(cell);
  //cell->cytosol_contents_["ATP Synthase"] = 0;
  cell->cytosol_contents_["Ribosome"] = 1000;

  SUBCASE("Testing in isolation")
  {
    DNAIf dnaif(pro, "ATP Synthase < 300");
    CHECK(dnaif.molecule_name_ == "ATP Synthase");
    CHECK(dnaif.inequality_type_ == "<");
    CHECK(dnaif.threshold_ == 300);

    CHECK(dnaif.check(*cell));
    cell->cytosol_contents_["ATP Synthase"] = 1000;
    CHECK(!dnaif.check(*cell));

    DNAThen dnathen(pro, "all = 0");
    DNAThen dnathen2(pro, "ATP Synthase += 1.0");
    
    dnathen.apply(cell->dna_.get());
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 0);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 0);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 0);
    dnathen2.apply(cell->dna_.get());
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 1.0);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 0);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 1.0);
    dnathen2.apply(cell->dna_.get());
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 2.0);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 0);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 2.0);
  }
  
  SUBCASE("Nested")
  {
    DNAIf dnaif(pro, "ATP Synthase < 300");

    CHECK(dnaif.check(*cell));
    cell->cytosol_contents_["ATP Synthase"] = 1000;
    CHECK(!dnaif.check(*cell));

    DNAThen::Ptr dnathen(new DNAThen(pro, "all = 0"));
    DNAThen::Ptr dnathen2(new DNAThen(pro, "ATP Synthase += 1.0"));

    dnaif.thens_.push_back(dnathen);
    dnaif.execute(cell.get());
    
    // Transcription factors should be unchanged because ATP Synthase is high enough.
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 1);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 1);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 2);

    // Now it drops too low, but we only have the first `dnathen` attached to `dnaif`.
    cell->cytosol_contents_["ATP Synthase"] = 200;
    dnaif.execute(cell.get());
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 0);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 0);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 0);

    // Now add the second DNAThen.
    dnaif.thens_.push_back(dnathen2);
    dnaif.execute(cell.get());
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 1.0);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 0);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 1.0);
  }

  SUBCASE("YAML construction")
  {
    DNAIf dnaif(pro, YAML::Load("if: ATP Synthase < 300\n"
                                "then:\n"
                                "  - all = 0\n"
                                "  - ATP Synthase += 1.0"));
    
    // Transcription factors should be unchanged because ATP Synthase is high enough.
    cell->cytosol_contents_["ATP Synthase"] = 1000;
    dnaif.execute(cell.get());
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 1);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 1);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 2);

    // Now transcription factors should change.
    cell->cytosol_contents_["ATP Synthase"] = 200;
    dnaif.execute(cell.get());
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 1.0);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 0);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 1.0);
  }

  SUBCASE("DNA::tick()")
  {
    Biome::Ptr biome(new Biome(pro, 10, "Alkaline vents"));
    pro.biomes_.push_back(biome);

    cell->dna_->dna_ifs_.push_back(DNAIf::Ptr(new DNAIf(pro, YAML::Load("if: ATP Synthase < 300\n"
                                                                        "then:\n"
                                                                        "  - all = 0\n"
                                                                        "  - ATP Synthase += 1.0"))));

    cell->cytosol_contents_["ATP Synthase"] = 200;
    
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 1.0);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 1.0);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 2.0);
    
    cell->tick(*biome);
    
    CHECK(cell->dna_->transcription_factors_["ATP Synthase"] == 1.0);
    CHECK(cell->dna_->transcription_factors_["Ribosome"] == 0);
    CHECK(cell->dna_->transcription_factors_.vals_.sum() == 1.0);
  }
}

TEST_CASE("Protein synthesis rate as a function of num amino acids")
{
  Prokaryotic pro;
  
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Protein X", "", 0, 100)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Protein 2X", "", 0, 200)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Amino acids", "", 0)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Ribosome", "", 0, 150000)));  // Really big so they don't really get produced
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ADP", "", 0)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ATP", "", 0)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Phosphate", "", 0)));

  // Provide infinite AAs and ATP.
  Biome::Ptr biome(new Biome(pro, 10, "Alkaline vents"));
  pro.biomes_.push_back(biome);
  biome->concentrations_["Amino acids"] = 200;
  biome->concentrations_["ATP"] = 200;

  Cell::Ptr cell(new Cell(pro, "cell"));
  pro.cells_.push_back(cell);
  // AAs and ATP come in through the membrane very fast.  ADP leaves very fast.
  cell->membrane_permeabilities_["Amino acids"] = 1.0;  
  cell->membrane_permeabilities_["ATP"] = 1.0;
  cell->membrane_permeabilities_["ADP"] = 1.0;

  // Start off with some ribosomes
  cell->cytosol_contents_["Ribosome"] = 2e4;  // should get 10k each to ribosomes, X, and 2X.

  // Ok just don't make any ribosomes, it's complicating things.
  cell->dna_->dna_ifs_.push_back(DNAIf::Ptr(new DNAIf(pro, YAML::Load("if: Ribosome > 1\n"
                                                                      "then:\n"
                                                                      "  - Ribosome = 0.0"))));

  for (int i = 0; i < 10; ++i) {  
    pro.tick();
    // With equal transcription factors (the default), confirm that 2*num_amino_acids_ means 1/2 the production rate.
    // Allow 5% error.
    cout << "cytosol_contents_" << endl << cell->cytosol_contents_.str("  ") << endl;
    //cout << "transcription_factors_" << endl << cell->dna_->transcription_factors_.str("  ") << endl;
    cout << "ribosome_assignments_" << endl << cell->dna_->ribosome_assignments_.str("  ") << endl;
    CHECK(2*cell->cytosol_contents_["Protein 2X"] / cell->cytosol_contents_["Protein X"] == doctest::Approx(1).epsilon(0.05));
  }
  
  for (int i = 0; i < 100; ++i) {
    pro.tick();
    // cout << "------" << endl;
    // cout << cell->str() << endl;
    // cout << "concentrations: " << endl << cell->cytosolConcentrations().str("  ") << endl;
    //cin.ignore();
  }
  cout << cell->str() << endl;
  cout << "concentrations: " << endl << cell->cytosolConcentrations().str("  ") << endl;
  CHECK(2*cell->cytosol_contents_["Protein 2X"] / cell->cytosol_contents_["Protein X"] == doctest::Approx(1).epsilon(0.05));

  CHECK(cell->cytosol_contents_["Ribosome"] == 2e4);
}

// We expect this cell to remain static - running tick() a lot shouldn't change anything, e.g.
// bc of the biome <> cell permeability math.
// TEST_CASE("Static cell")
// {
// }

// https://www.sciencedirect.com/topics/neuroscience/adenosine-triphosphate
// Normally cellular ATP concentration is maintained in the range of 1 to 10 mmol/L, with a normal ratio of ATP/ADP of approximately 1000.
// TEST_CASE()
// {
// }

// https://www.sciencedirect.com/topics/neuroscience/adenosine-triphosphate
// Normal cellular ATP is 1-10 mM.  (ADP is typically 1000x lower.)
// TEST_CASE()
// {
// }
