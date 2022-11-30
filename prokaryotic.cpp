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
using Eigen::VectorXd;

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
  Eigen::VectorXd vals_;
  
  MoleculeVals(const Prokaryotic& pro);
  
  std::string _str() const;
  const double& operator[](int idx) const { return vals_.coeffRef(idx); }
  double operator[](int idx) { return vals_.coeffRef(idx); }
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
  double um3_;
//  ReactionTable reaction_table_;
  MoleculeVals cytosol_contents_;
  MoleculeVals membrane_contents_;
  
  Cell(const Prokaryotic& pro, const std::string& name);

  std::string _str() const;
  void tick(const Biome& biome);
};

// Contains the whole simulation model
class Prokaryotic : public Printable
{
public:
  Prokaryotic();
  std::vector<MoleculeType::ConstPtr> molecule_types_;
  std::map<std::string, MoleculeType::ConstPtr> molecule_map_;

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
};


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
  vals_ = VectorXd::Zero(MoleculeType::num_molecule_types_);
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
  um3_(0.6),
  cytosol_contents_(pro_),
  membrane_contents_(pro_)
{}

std::string Cell::_str() const
{
  std::ostringstream oss;
  oss << "Cell \"" << name_ << "\"" << endl;
  oss << "  um3_: " << um3_ << endl;
  oss << "  cytosol_contents_: " << endl << cytosol_contents_.str("    ") << endl;
  oss << "  membrane_contents_: " << endl << membrane_contents_.str("    ") << endl;
  return oss.str();
}

void Cell::tick(const Biome& biome)
{
  // Determine new cell contents as things permeate the membrane (in both directions)
  
  
  // Run reactions in cell
}

Prokaryotic::Prokaryotic()
{
  // addMoleculeType(MoleculeType::Ptr(new MoleculeType("ATP", ":bang:", 503.15)));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("ADP", ":briefcase:", 423.17)));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("Phosphate", "P", 94.97)));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("R", "R", 1000)));
  std::vector<MoleculeType::ConstPtr> constituents;
  constituents.push_back(molecule("ADP"));
  constituents.push_back(molecule("Phosphate"));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("ATP", ":bang:", constituents)));

  for (auto mt : molecule_types_)
    cout << mt->str() << endl;

  cells_.push_back(Cell::Ptr(new Cell(*this, "aoeu")));
  cout << "Cell: " << endl;
  
  biomes_.push_back(Biome::Ptr(new Biome(*this, 10, "Alkaline vents")));
  biomes_[0]->concentrations_["Phosphate"] = 24;
  biomes_[0]->concentrations_["R"] = 10;

  cout << str() << endl;
}

std::string Prokaryotic::_str() const
{
  std::ostringstream oss;
  oss << "Simulation status" << endl;
  for (auto cell : cells_)
    oss << cell->str("  ") << endl;
  for (auto biome : biomes_)
    oss << biomes_[0]->str("  ") << endl;
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

int main(int argc, char** argv)
{
  Prokaryotic prokaryotic;
  prokaryotic.run();
  return 0;
}
