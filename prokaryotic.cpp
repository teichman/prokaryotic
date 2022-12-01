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
  
  std::string _str() const;
  const double& operator[](int idx) const { return vals_.coeffRef(idx); }
  double& operator[](int idx) { return vals_.coeffRef(idx); }
  const double& operator[](const std::string& name) const;
  double& operator[](const std::string& name);
};

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
  std::vector<MoleculeType::ConstPtr> constituents_;

  MoleculeType(const std::string& name, const std::string& symbol, double daltons);
  MoleculeType(const std::string& name, const std::string& symbol, const std::vector<MoleculeType::ConstPtr>& constituents);
  std::string _str() const;
  static size_t numMoleculeTypes() { return num_molecule_types_; }
  
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

class Cell : public Printable
{
public:
  typedef std::shared_ptr<Cell> Ptr;
  typedef std::shared_ptr<const Cell> ConstPtr;

  const Prokaryotic& pro_;
  std::string name_;
  double um3_;  // 1e-18 m3, or 1e-15 L (i.e. femtoliter)
//  ReactionTable reaction_table_;
  MoleculeVals cytosol_contents_;
  MoleculeVals membrane_contents_;
  MoleculeVals membrane_permeabilities_;
  
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

class ReactionType : public Printable
{
public:
  typedef std::shared_ptr<const ReactionType> ConstPtr;

  const Prokaryotic& pro_;
  MoleculeVals inputs_;
  MoleculeVals outputs_;
  double p_spon_;  // probability of spontaneous reaction if inputs collide
  
  ReactionType(const Prokaryotic& pro, const MoleculeVals& inputs, const MoleculeVals& outputs, double p_spon);
  void tick(const Biome& biome, Cell& cell) const;
  std::string _str() const;
};

// Contains the whole simulation model
class Prokaryotic : public Printable
{
public:
  Prokaryotic();
  std::vector<MoleculeType::ConstPtr> molecule_types_;
  std::map<std::string, MoleculeType::ConstPtr> molecule_map_;
  std::vector<ReactionType::ConstPtr> reaction_types_;

  std::vector<Biome::Ptr> biomes_;
  std::vector<Cell::Ptr> cells_;
  
  void addMoleculeType(MoleculeType::ConstPtr mt);
  
  MoleculeType::ConstPtr molecule(const std::string& name) const {
    assert(molecule_map_.find(name) != molecule_map_.end());
    return molecule_map_.at(name);
  }
  
  size_t moleculeIdx(const std::string& name) const { return molecule(name)->idx_; }
  const std::string& moleculeName(size_t idx) const { return molecule_types_[idx]->name_; }
  void tick();
  std::string _str() const;
  void run();
  void runTests();
};


ReactionType::ReactionType(const Prokaryotic& pro,
                           const MoleculeVals& inputs, const MoleculeVals& outputs,
                           double p_spon) :
  pro_(pro),
  inputs_(pro),
  outputs_(pro),
  p_spon_(p_spon)
{
}

void ReactionType::tick(const Biome& biome, Cell& cell) const
{
  int num_collisions = 0;
  
  // compute number of collisions for this type
  // If there is only one input, then it's easy: it's just that number.
  if (inputs_.nz() == 1) {
    int idx = inputs_.nz_indices()[0];
    num_collisions = cell.cytosol_contents_[idx];
  }
  // otherwise do something smart
  else {
    num_collisions = 0;
  }
  
  // apply probability of reaction occurring
  int num_reactions = p_spon_ * num_collisions;
  
  // update cell.cytosol_counts_
  cell.cytosol_contents_.vals_ -= inputs_.vals_ * num_reactions;
  cell.cytosol_contents_.vals_ += outputs_.vals_ * num_reactions;
}

std::string ReactionType::_str() const
{
  std::ostringstream oss;
  oss << "Reaction" << endl;
  oss << "  Inputs: " << endl;
  for (int i = 0; i < inputs_.vals_.size(); ++i)
    if (inputs_.vals_[i] > 0)
      oss << "    " << inputs_.vals_[i] << "x " << pro_.molecule_types_[i]->name_ << endl;
  for (int i = 0; i < outputs_.vals_.size(); ++i)
    if (outputs_.vals_[i] > 0)
      oss << "    " << outputs_.vals_[i] << "x " << pro_.molecule_types_[i]->name_ << endl;  
  return oss.str();
}

MoleculeType::MoleculeType(const std::string& name,
                           const std::string& symbol,
                           double daltons) :
  idx_(num_molecule_types_),
  name_(name),
  symbol_(symbol),
  daltons_(daltons)
{
  num_molecule_types_ += 1;
}

MoleculeType::MoleculeType(const std::string& name, const std::string& symbol,
                           const std::vector<MoleculeType::ConstPtr>& constituents) :
  idx_(num_molecule_types_),
  name_(name),
  symbol_(symbol),
  daltons_(0),
  constituents_(constituents)
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
  return oss.str();
}

MoleculeVals::MoleculeVals(const Prokaryotic& pro) :
  pro_(pro)
{
  vals_ = ArrayXd::Zero(MoleculeType::num_molecule_types_);
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
{
  membrane_permeabilities_["R"] = 0.005;
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
  
  // In-cell reactions
  // cytosol_concentrations_["R"] = std::max(0.0, cytosol_concentrations_["R"] - 0.1);
  // cytosol_concentrations_["Phosphate"] -= 0.5;
  // cytosol_concentrations_["X"] += 0.1;

  for (auto reaction : pro_.reaction_types_) {
    reaction->tick(biome, *this);
  }
}

Prokaryotic::Prokaryotic()
{
  // addMoleculeType(MoleculeType::Ptr(new MoleculeType("ATP", ":bang:", 503.15)));
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
    MoleculeVals inputs(*this);
    inputs["R"] = 1;
    MoleculeVals outputs(*this);
    outputs["X"] = 2;
    reaction_types_.push_back(ReactionType::ConstPtr(new ReactionType(*this, inputs, outputs, 0.01)));
  }
    
  for (auto mt : molecule_types_)
    cout << mt->str() << endl;

  cells_.push_back(Cell::Ptr(new Cell(*this, "aoeu")));
  
  biomes_.push_back(Biome::Ptr(new Biome(*this, 10, "Alkaline vents")));
  biomes_[0]->concentrations_["Phosphate"] = 24;
  biomes_[0]->concentrations_["R"] = 10;
  biomes_[0]->concentrations_["X"] = 0;
  
  cells_[0]->setCytosolContentsByConcentrations(biomes_[0]->concentrations_);

  cout << str() << endl;
}

std::string Prokaryotic::_str() const
{
  std::ostringstream oss;
  oss << "Simulation status" << endl;
  oss << "  Reactions" << endl;  
  for (auto rt : reaction_types_)
    oss << rt->str("    ") << endl;
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
  cout << "Running." << endl;
  while(true)
  {
    tick();
    cout << str() << endl;
    cout << "Press RET to continue." << endl;
    std::cin.get();
  }
}

void Prokaryotic::addMoleculeType(MoleculeType::ConstPtr mt)
{
  molecule_types_.push_back(mt);
  molecule_map_[mt->name_] = mt;
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
}


int main(int argc, char** argv)
{
  Prokaryotic prokaryotic;
  prokaryotic.runTests();
  //prokaryotic.run();
  return 0;
}
