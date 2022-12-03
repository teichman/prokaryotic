#include <chrono>
#include <thread>
#include <limits>
#include <Eigen/Dense>
#include <cassert>
#include <sstream>
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include <string>

using std::shared_ptr;
using std::cout, std::endl;
using std::vector;
using Eigen::ArrayXd;

class Prokaryotic;
class Cell;
class Biome;

// Just gives the prefix option for free.
class Printable
{
public:
  virtual ~Printable() {}

  // Ugh maybe I should go back to python.
  // https://stackoverflow.com/questions/37803255/how-to-call-base-classs-function-that-is-not-defined-in-derived-class
  virtual std::string _str() const = 0;
  std::string str(const std::string& prefix="") const
  {
    std::istringstream iss(_str());
    std::ostringstream oss;
    for (std::string line; std::getline(iss, line); )
      oss << prefix << line + "\n";

    // Remove trailing newline.
    std::string result = oss.str();
    result.erase(std::remove(result.end() - 2, result.end(), '\n'), result.cend());
    return result;
  }
};

class MoleculeVals : public Printable
{
public:
  const Prokaryotic& pro_;
  Eigen::ArrayXd vals_;
  
  MoleculeVals(const Prokaryotic& pro);

  int nz() const;
  std::vector<int> nz_indices() const;
  void round() { vals_ = vals_.round(); }
  void probabilisticRound();
  
  std::string _str() const;
  const double& operator[](int idx) const { return vals_.coeffRef(idx); }
  double& operator[](int idx) { return vals_.coeffRef(idx); }
  const double& operator[](const std::string& name) const;
  double& operator[](const std::string& name);
};

class ReactionType : public Printable
{
public:
  typedef std::shared_ptr<const ReactionType> ConstPtr;

  const Prokaryotic& pro_;
  MoleculeVals inputs_;
  MoleculeVals outputs_;
  MoleculeVals kms_;  // concentrations that yield half the max rate.  Each element is mM.
  double kcat_;  // max rate per tick (1 tick is 1 second?)
  
  ReactionType(const Prokaryotic& pro,
               const MoleculeVals& inputs, const MoleculeVals& outputs,
               const MoleculeVals& kms, double kcat);
  void tick(Cell& cell, int num_protein_copies) const;
  std::string _str() const;
  // See Fig 1 of http://book.bionumbers.org/how-many-reactions-do-enzymes-carry-out-each-second/
  // substrate_concentration in mM, km in mM, kcat in reactions / sec
  static double rate(double substrate_concentration, double km, double kcat);
};

double probabilityPerSecond(double half_life_hours)
{
  if (half_life_hours == std::numeric_limits<double>::max())
    return 0.0;

  double half_life_ticks = half_life_hours * 3600;
  return 1.0 - exp(log(0.5) / half_life_ticks);
}

// Superclass for proteins and small molecules, both of which we will just call "Molecules".
class MoleculeType : public Printable
{
public:
  typedef std::shared_ptr<MoleculeType> Ptr;
  typedef std::shared_ptr<const MoleculeType> ConstPtr;
  
  int idx_;
  std::string name_;
  std::string symbol_;
  double daltons_;
  double half_life_hours_;
  std::vector<MoleculeType::ConstPtr> constituents_;
  ReactionType::ConstPtr reaction_;  // If a protein, this is the reaction it catalyzes 

  MoleculeType(const std::string& name, const std::string& symbol, double daltons,
               double half_life_hours = std::numeric_limits<double>::max(),
               ReactionType::ConstPtr reaction = ReactionType::ConstPtr(nullptr));
  MoleculeType(const std::string& name, const std::string& symbol,
               const std::vector<MoleculeType::ConstPtr>& constituents,
               double half_life_hours = std::numeric_limits<double>::max(),
               ReactionType::ConstPtr reaction = ReactionType::ConstPtr(nullptr));
  
  std::string _str() const;
  static size_t numMoleculeTypes() { return num_molecule_types_; }
  double pDenature() const { return probabilityPerSecond(half_life_hours_); }
  
private:
  static size_t num_molecule_types_;
  friend class MoleculeVals;
};

size_t MoleculeType::num_molecule_types_ = 0;

class Biome : public Printable
{
public:
  typedef std::shared_ptr<Biome> Ptr;
  typedef std::shared_ptr<const Biome> ConstPtr;

  const Prokaryotic& pro_;
  double m3_;
  MoleculeVals concentrations_;  // mM
  // MoleculeVals fluxes_;
  // MoleculeVals contents_;
  std::string name_;
  
  Biome(const Prokaryotic& pro, double m3, const std::string& name);
  std::string _str() const;
  void tick();
};


// // Contains the high level scripts for the cell.
// // e.g. high level functions like divide() are here.
// // This is where cell division happens based on some conditions, and protein synthesis regulation.
// // We're going to gloss over transcriptional vs translational regulation for now and just combine them into one big thing.
// class DNA : public Printable
// {
// public:
//   typedef std::shared_ptr<const DNA> ConstPtr;
//   typedef std::shared_ptr<DNA> Ptr;

//   const Prokaryotic& pro_;
//   // transcription_factors_ adjust the distribution of ribosomes to genes.
//   // User can set write simple scripts to change the transcription_factors_ in response to various things
//   // If the sum is zero, no ribosomes will do anything.
//   // If the sum is greater than 1, the distribution will be normalized.
//   MoleculeVals transcription_factors_;
//   std::vector<ReactionType::ConstPtr> synthesis_reactions_;
  
//   DNA(const Prokaryotic& pro) : pro_(pro), transcription_factors_(pro) {}
//   void tick(Cell& cell) const;
//   std::string _str() const { return ""; }
// };

class Cell : public Printable
{
public:
  typedef std::shared_ptr<Cell> Ptr;
  typedef std::shared_ptr<const Cell> ConstPtr;

  const Prokaryotic& pro_;
  std::string name_;
  double um3_;  // 1e-18 m3, or 1e-15 L (i.e. femtoliter)
  MoleculeVals cytosol_contents_;
  MoleculeVals membrane_contents_;
  MoleculeVals membrane_permeabilities_;
  // DNA::Ptr dna_;  // non-const so Prokaryotic can make changes to it.
  
  Cell(const Prokaryotic& pro, const std::string& name);

  std::string _str() const;
  void tick(const Biome& biome);
  MoleculeVals cytosolConcentrations() const;  // mM
  static MoleculeVals cytosolConcentrations(const Prokaryotic& pro, MoleculeVals cytosol_contents, double um3);  // mM
  // concentrations is mM
  void setCytosolContentsByConcentrations(const MoleculeVals& cytosol_concentrations);
  // concentrations is mM
  static MoleculeVals cytosolContents(const Prokaryotic& pro, const MoleculeVals& cytosol_concentrations, double um3);
};

// Contains the whole simulation model
// Data protection strategy: Data is stored as (non-const) Ptr, but is only handed out as ConstPtr.
class Prokaryotic : public Printable
{
public:
  Prokaryotic();
  
  void addMoleculeType(MoleculeType::Ptr mt);
  
  MoleculeType::ConstPtr molecule(const std::string& name) const {
    assert(molecule_map_.find(name) != molecule_map_.end());
    return molecule_map_.at(name);
  }
  
  size_t moleculeIdx(const std::string& name) const { return molecule(name)->idx_; }
  const std::string& moleculeName(size_t idx) const { return molecule_types_[idx]->name_; }
  std::vector<MoleculeType::ConstPtr> moleculeTypes() const;
  
  void tick();
  std::string _str() const;
  void run();
  void runTests();
  void initializeHardcoded();

private:
  std::vector<MoleculeType::Ptr> molecule_types_;
  std::map<std::string, MoleculeType::Ptr> molecule_map_;
  std::vector<Biome::Ptr> biomes_;
  std::vector<Cell::Ptr> cells_;

  MoleculeType::Ptr _molecule(const std::string& name) const {
    assert(molecule_map_.find(name) != molecule_map_.end());
    return molecule_map_.at(name);
  }
};

ReactionType::ReactionType(const Prokaryotic& pro,
                           const MoleculeVals& inputs, const MoleculeVals& outputs,
                           const MoleculeVals& kms, double kcat) :
  pro_(pro),
  inputs_(inputs),
  outputs_(outputs),
  kms_(kms),
  kcat_(kcat)
{
  for (int i = 0; i < inputs_.vals_.size(); ++i) {
    if (inputs_.vals_[i] > 0)
      assert(kms_.vals_[i] > 0);
  }
}

void ReactionType::tick(Cell& cell, int num_protein_copies) const
{
  // Simplify: Assume that everything follows Michaelis-Menten kinetics, even for the (presumably majority) of reactions
  // that have multiple substrates.
  // We'll just take the min reaction rate across the substrates (inputs).
  // https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics
  MoleculeVals concentrations = cell.cytosolConcentrations();
  double minrate = std::numeric_limits<double>::max();
  for (int i = 0; i < inputs_.vals_.size(); ++i)
    if (inputs_[i] > 0)
      minrate = std::min(minrate, rate(concentrations[i], kms_[i], kcat_));
  double rate = minrate;  // reactions / tick (1 tick == 1 second?)
  
  // update cell.cytosol_counts_
  cell.cytosol_contents_.vals_ -= inputs_.vals_ * rate * num_protein_copies;
  cell.cytosol_contents_.vals_ += outputs_.vals_ * rate * num_protein_copies;

  // Should only need one of these at the end of Cell::tick()
  //cell.cytosol_contents_.probabilisticRound();
}

std::string ReactionType::_str() const
{
  std::ostringstream oss;
  oss << "Reaction" << endl;
  oss << "  Inputs & kms: " << endl;
  for (int i = 0; i < inputs_.vals_.size(); ++i)
    if (inputs_.vals_[i] > 0)
      oss << "    " << inputs_.vals_[i] << "x " << pro_.moleculeName(i) << " (km: " << kms_[i] << " mM)" << endl;
  oss << "  Outputs: " << endl;
  for (int i = 0; i < outputs_.vals_.size(); ++i)
    if (outputs_.vals_[i] > 0)
      oss << "    " << outputs_.vals_[i] << "x " << pro_.moleculeName(i) << endl;
  oss << "  kcat_: " << kcat_ << endl;  
  return oss.str();
}

double ReactionType::rate(double substrate_concentration, double km, double kcat)
{
  // https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics
  return kcat * substrate_concentration / (km + substrate_concentration);
}

MoleculeType::MoleculeType(const std::string& name,
                           const std::string& symbol,
                           double daltons,
                           double half_life_hours,
                           ReactionType::ConstPtr reaction) :
  idx_(num_molecule_types_),
  name_(name),
  symbol_(symbol),
  daltons_(daltons),
  half_life_hours_(half_life_hours),
  reaction_(reaction)
{
  num_molecule_types_ += 1;
}

MoleculeType::MoleculeType(const std::string& name, const std::string& symbol,
                           const std::vector<MoleculeType::ConstPtr>& constituents,
                           double half_life_hours,
                           ReactionType::ConstPtr reaction) :
  idx_(num_molecule_types_),
  name_(name),
  symbol_(symbol),
  daltons_(0),
  half_life_hours_(half_life_hours),
  constituents_(constituents),
  reaction_(reaction)
{
  num_molecule_types_ += 1;
  daltons_ = 0;
  for (auto mt : constituents_)
    daltons_ += mt->daltons_;
}

std::string MoleculeType::_str() const
{
  std::ostringstream oss;
  oss << "MoleculeType \"" << name_ << "\"" << endl;
  oss << "  idx_: " << idx_ << endl;
  oss << "  symbol_: " << symbol_ << endl;
  oss << "  daltons_: " << daltons_ << endl;
  oss << "  Constituents: " << endl;
  for (auto mt : constituents_)
    oss << "    " << mt->name_ << endl;
  if (!reaction_)
    oss << "  Reaction: None" << endl;
  else {
    oss << "  Reaction" << endl;
    oss << reaction_->str("    ") << endl;
  }
  
  return oss.str();
}

MoleculeVals::MoleculeVals(const Prokaryotic& pro) :
  pro_(pro)
{
  vals_ = ArrayXd::Zero(pro.moleculeTypes().size());
}

int MoleculeVals::nz() const
{
  int num = 0;
  for (int idx = 0; idx < vals_.size(); ++idx)
    if (vals_.coeffRef(idx) == 0.0)
      num += 1;
  return num;
}

std::vector<int> MoleculeVals::nz_indices() const
{
  vector<int> indices;
  for (int idx = 0; idx < vals_.size(); ++idx)
    if (vals_.coeffRef(idx) != 0.0)
      indices.push_back(idx);
  return indices;
}

std::string MoleculeVals::_str() const
{
  std::ostringstream oss;
  for (int idx = 0; idx < vals_.size(); ++idx)
    oss << pro_.moleculeName(idx) << ": " << vals_[idx] << endl;
  return oss.str();
}

const double& MoleculeVals::operator[](const std::string& name) const
{
  return vals_.coeffRef(pro_.moleculeIdx(name));
}

double& MoleculeVals::operator[](const std::string& name)
{
  return vals_.coeffRef(pro_.moleculeIdx(name));
}

void MoleculeVals::probabilisticRound()
{
  for (int i = 0; i < vals_.size(); ++i) {
    double thresh = (ArrayXd::Random(1)[0] + 1.0) / 2.0;  // slow but convenient.  ArrayXd::Random generates random doubles in [-1, 1].
    if (vals_[i] - int(vals_[i]) > thresh)
      vals_[i] = std::ceil(vals_[i]);
    else
      vals_[i] = std::floor(vals_[i]);
  }
}


Biome::Biome(const Prokaryotic& pro, double m3, const std::string& name) :
  pro_(pro),
  m3_(m3),
  concentrations_(pro),
  name_(name)
{
}

std::string Biome::_str() const
{
  std::ostringstream oss;
  oss << "Biome \"" << name_ << "\"" << endl;
  oss << "  m3_: " << m3_ << endl;
  oss << "  concentrations_: " << endl << concentrations_.str("    ") << endl;
  return oss.str();
}

void Biome::tick()
{  
}

Cell::Cell(const Prokaryotic& pro, const std::string& name) :
           
  pro_(pro),
  name_(name),
  um3_(1.0),
  cytosol_contents_(pro_),
  membrane_contents_(pro_),
  membrane_permeabilities_(pro_)
  // dna_(new DNA(pro))
{
  membrane_permeabilities_["R"] = 0.0000005;
  membrane_permeabilities_["Phosphate"] = 0.1;
  membrane_permeabilities_["X"] = 0.001;
}

std::string Cell::_str() const
{
  std::ostringstream oss;
  oss << "Cell \"" << name_ << "\"" << endl;
  oss << "  um3_: " << um3_ << endl;
  oss << "  cytosol_contents_: " << endl << cytosol_contents_.str("    ") << endl;
  oss << "  membrane_contents_: " << endl << membrane_contents_.str("    ") << endl;
  return oss.str();
}

MoleculeVals Cell::cytosolConcentrations(const Prokaryotic& pro, MoleculeVals cytosol_contents, double um3)
{
  // millimoles = 1e3 * (num / 6e23)
  // liters = um3 * 1e-15, bc um3 is fL.
  // mM = 1e3 * (num / 6e23) / (um3 * 1e-15)
  // mM = 1.67e-06 * num / um3
  
  MoleculeVals concentrations(pro);
  concentrations.vals_ = cytosol_contents.vals_ * (1.67e-6 / um3);
  return concentrations;
}

MoleculeVals Cell::cytosolConcentrations() const
{
  return cytosolConcentrations(pro_, cytosol_contents_, um3_);
}

MoleculeVals Cell::cytosolContents(const Prokaryotic& pro, const MoleculeVals& cytosol_concentrations, double um3)
{
  MoleculeVals contents(pro);
  contents.vals_ = cytosol_concentrations.vals_ / (1.67e-6 / um3);
  return contents;
}

void Cell::setCytosolContentsByConcentrations(const MoleculeVals& cytosol_concentrations)
{
  cytosol_contents_.vals_ = cytosolContents(pro_, cytosol_concentrations, um3_).vals_;
}

void Cell::tick(const Biome& biome)
{
  // Passive membrane permeations
  MoleculeVals cytosol_concentrations = cytosolConcentrations();
  auto deltas = (biome.concentrations_.vals_ - cytosol_concentrations.vals_);
  auto step = membrane_permeabilities_.vals_ * deltas;
  cytosol_concentrations.vals_ += step;
  setCytosolContentsByConcentrations(cytosol_concentrations);
  
  // In-cell reactions.  For now we are assuming all reactions require a protein.
  for (auto mt : pro_.moleculeTypes()) {
    if (mt->reaction_) {
      mt->reaction_->tick(*this, cytosol_contents_[mt->idx_]);
    }
  }

  // Apply degredation of molecules.
  // TODO: Have a member of MoleculeType (or, better, ProteinType) that indicates if the protein is denatured or not.
  // Then add proteasomes that identify denatured proteins (through ubiquitin, though maybe we don't want that ... or maybe we do??)
  // Yeah just do proteasomes for now.
  for (auto mt : pro_.moleculeTypes())
    if (mt->pDenature() > 0)
      cytosol_contents_[mt->idx_] *= (1.0 - mt->pDenature());

  // // Apply "DNA programming"
  // dna_->tick(*this);
  
  cytosol_contents_.probabilisticRound();
}

// void DNA::tick(Cell& cell) const
// {
//   // Eventually this code will be programmable by the player.  For now we're just hardcoding it.

//   transcription_factors_.vals_ = 0;
  
//   1// Probably better would be something that adjusts the rate of ribosomal construction of ATP Synthase based on the
//   // concentration of ATP Synthase.
//   if (cell.cytosol_contents_["ATP Synthase"] < 200) {
//     // Run the reaction that generates ATP Synthase.
//     // It's a reaction run by Ribosomes (just like other proteins) that depends on concentration of substrates (building blocks,
//     // Approx 5 ATP for each amino acid in the sequence.
//     // 200 amino acids per minute made by ribosomes.  This is ~3 / second.  Maybe this is how we set kcat for the ribosome, and then
//     // the synthesis rate drops if there aren't enough building blocks around.  Not sure how to set the KM values though.
//     // Typically 100-600 amino acids per protein.  (Some are crazy long, but maybe that is just eukaryotes..)
//     // Probably we want to ignore the fact that it's 3aa / second and just absorb that into kcat though.
//     // ReactionType for protein synthesis is just a ReactionType stored in DNA class (which includes ribosome count) with kcat
//     // set by num_aa of the protein, and num_protein_copies is determined by the number of ribosomes assigned to this protein.
//     // That in turn is set by a distribution over genes that the user can set, also in the DNA programming - you can set promotion
//     // factors in R^+, one for each gene, and ribosomes get assigned to genes according to a normalized distribution over genes.
//     transcription_factors_["ATP Synthase"] = 1.0;
//   }

//   // This part the user cannot touch.

//   // If the transcription factors sum to greater than one, make them sum to one.
//   // Otherwise leave them alone.
//   MoleculeVals normalized_transcription_factors(pro_);
//   if (transcription_factors_.vals_.sum() > 1.0)
//     normalized_transcription_factors.vals_ = transcription_factors_.vals_ / transcription_factors_.vals_.sum();
//   else
//     normalized_transcription_factors.vals_ = transcription_factors_.vals_;
  
//   for (int i = 0; i < transcription_factors_.vals_.size(); ++i) {
//     if (synthesis_reactions_[i] && normalized_transcription_factors[i] > 0)
//       synthesis_reactions_[i]->tick(*this, cell.cytosol_contents_["Ribosome"] * normalized_transcription_factors[i]);
// }

Prokaryotic::Prokaryotic()
{}

void Prokaryotic::initializeHardcoded()
{
  // Create all MoleculeTypes first so we avoid the complexity in updating all the MoleculeVals when we add new MoleculeTypes.
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("ADP", ":briefcase:", 423.17)));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("Phosphate", "P", 94.97)));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("X", "X", 500)));
  {
    std::vector<MoleculeType::ConstPtr> constituents;
    constituents.push_back(molecule("X"));
    constituents.push_back(molecule("X"));
    addMoleculeType(MoleculeType::Ptr(new MoleculeType("R", "R", constituents)));
  }
  
  {
    std::vector<MoleculeType::ConstPtr> constituents;
    constituents.push_back(molecule("ADP"));
    constituents.push_back(molecule("Phosphate"));
    addMoleculeType(MoleculeType::Ptr(new MoleculeType("ATP", ":bang:", constituents)));
  }

  { 
    std::vector<MoleculeType::ConstPtr> constituents;
    constituents.push_back(molecule("X"));
    constituents.push_back(molecule("X"));
    constituents.push_back(molecule("R"));
    constituents.push_back(molecule("R"));
    constituents.push_back(molecule("Phosphate"));
    addMoleculeType(MoleculeType::Ptr(new MoleculeType("ATP Synthase", ":hammer:", constituents, 1.0)));
  }

  { 
    std::vector<MoleculeType::ConstPtr> constituents;
    constituents.push_back(molecule("X"));
    constituents.push_back(molecule("X"));
    constituents.push_back(molecule("R"));
    constituents.push_back(molecule("R"));
    constituents.push_back(molecule("Phosphate"));
    addMoleculeType(MoleculeType::Ptr(new MoleculeType("ATP Consumer", ":gear:", constituents)));
  }
  // { 
  //   std::vector<MoleculeType::ConstPtr> constituents;
  //   constituents.push_back(molecule("X"));
  //   constituents.push_back(molecule("X"));
  //   constituents.push_back(molecule("R"));
  //   constituents.push_back(molecule("R"));
  //   constituents.push_back(molecule("Phosphate"));
  //   addMoleculeType(MoleculeType::Ptr(new MoleculeType("Ribosome", ":factory:", constituents)));
  // }

  
  // Now add Reactions to MoleculeTypes.
  {
    MoleculeVals inputs(*this);
    inputs["ADP"] = 1;
    inputs["Phosphate"] = 1;
    MoleculeVals outputs(*this);
    outputs["ATP"] = 1;
    MoleculeVals kms(*this);
    kms["ADP"] = 1e-1;
    kms["Phosphate"] = 1e-2;
    double kcat = 0.001;
    _molecule("ATP Synthase")->reaction_ = ReactionType::ConstPtr(new ReactionType(*this, inputs, outputs, kms, kcat));
  }

  {
    MoleculeVals inputs(*this);
    inputs["X"] = 2;
    inputs["ATP"] = 1;
    MoleculeVals outputs(*this);
    outputs["R"] = 1;
    outputs["ADP"] = 1;
    outputs["Phosphate"] = 1;
    MoleculeVals kms(*this);
    kms["X"] = 1e-1;
    kms["ATP"] = 1e-3;
    double kcat = 0.01;
    _molecule("ATP Consumer")->reaction_ = ReactionType::ConstPtr(new ReactionType(*this, inputs, outputs, kms, kcat));
  }
  
  for (auto mt : molecule_types_)
    cout << mt->str() << endl;

  cells_.push_back(Cell::Ptr(new Cell(*this, "aoeu")));
  
  biomes_.push_back(Biome::Ptr(new Biome(*this, 10, "Alkaline vents")));
  biomes_[0]->concentrations_["Phosphate"] = 0.01;
  biomes_[0]->concentrations_["R"] = 0;
  biomes_[0]->concentrations_["X"] = 10;
  
  cells_[0]->setCytosolContentsByConcentrations(biomes_[0]->concentrations_);
  cells_[0]->cytosol_contents_["ADP"] = 1e6;
  cells_[0]->cytosol_contents_["ATP Synthase"] = 1e4;
  cells_[0]->cytosol_contents_["ATP Consumer"] = 1e3;

  cout << str() << endl;
}

std::string Prokaryotic::_str() const
{
  std::ostringstream oss;
  oss << "Simulation status" << endl;
  // oss << "  Reactions" << endl;  
  // for (auto rt : reaction_types_)
  //   oss << rt->str("    ") << endl;
  oss << "  Cells" << endl;  
  for (auto cell : cells_)
    oss << cell->str("    ") << endl;
  oss << "  Biomes" << endl;
  for (auto biome : biomes_)
    oss << biomes_[0]->str("    ") << endl;
  return oss.str();
}

void Prokaryotic::run()
{
  typedef std::chrono::high_resolution_clock HRC;
  
  cout << "Running." << endl;
  while(true)
  {
    auto start = HRC::now();
    tick();
    double seconds = double(std::chrono::duration_cast<std::chrono::nanoseconds>(HRC::now() - start).count()) * 1e-9;
    cout << str() << endl;
    cout << "seconds / tick: " << seconds << endl;
    cout << "ticks / second: " << 1.0 / seconds << endl;
    
    // cout << "Press RET to continue." << endl;
    // std::cin.get();
    //std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
}

void Prokaryotic::addMoleculeType(MoleculeType::Ptr mt)
{
  molecule_types_.push_back(mt);
  molecule_map_[mt->name_] = mt;
}


std::vector<MoleculeType::ConstPtr> Prokaryotic::moleculeTypes() const
{
  return std::vector<MoleculeType::ConstPtr>(molecule_types_.begin(), molecule_types_.end());
}

void Prokaryotic::tick()
{
  assert(biomes_.size() == 1);
  assert(cells_.size() == 1);
  for (auto biome : biomes_)
    biome->tick();
  for (auto cell : cells_)
    cell->tick(*biomes_[0]);
}

void Prokaryotic::runTests()
{
  // http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/
  // "rule of thumb that a concentration of 1 nM corresponds to roughly one copy of the molecule of interest per E. coli cell.
  // Hence, 100 mM means that there are roughly 10^8 copies of glutamate in each bacterium."
  MoleculeVals cytosol_contents(*this);
  cytosol_contents[0] = 1e8;
  cytosol_contents[1] = 1;
  double um3 = 1.0;
  MoleculeVals concentrations = Cell::cytosolConcentrations(*this, cytosol_contents, um3);
  cout << "Counts: " << endl;
  cout << cytosol_contents.str("  ") << endl;
  cout << "Concentrations: " << endl;
  cout << concentrations.str("  ") << endl;

  cout << "Round trip test" << endl;
  MoleculeVals contents2 = Cell::cytosolContents(*this, concentrations, um3);
  cout << "Counts: " << endl;
  cout << contents2.str("  ") << endl;

  //double half_life_hours = std::numeric_limits<double>::max();
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
}


int main(int argc, char** argv)
{
  Prokaryotic prokaryotic;
  prokaryotic.initializeHardcoded();
  //prokaryotic.runTests();
  prokaryotic.run();
  return 0;
}
