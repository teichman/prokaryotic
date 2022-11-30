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

// Superclass for proteins and small molecules
class MoleculeType : public Printable
{
public:
  typedef std::shared_ptr<MoleculeType> Ptr;
  typedef std::shared_ptr<const MoleculeType> ConstPtr;
  
  int id_;
  std::string name_;
  std::string symbol_;
  double daltons_;
  std::vector<MoleculeType::Ptr> constituents_;

  MoleculeType(const std::string& name, const std::string& symbol, double daltons);
  MoleculeType(const std::string& name, const std::string& symbol, const std::vector<MoleculeType::Ptr>& constituents);
  virtual std::string _str() const;
  static size_t numMoleculeTypes() { return num_molecule_types_; }
  
private:
  static size_t num_molecule_types_;
};

size_t MoleculeType::num_molecule_types_ = 0;

MoleculeType::MoleculeType(const std::string& name,
                           const std::string& symbol,
                           double daltons) :
  id_(num_molecule_types_),
  name_(name),
  symbol_(symbol),
  daltons_(daltons)
{
  num_molecule_types_ += 1;
}

MoleculeType::MoleculeType(const std::string& name, const std::string& symbol,
                           const std::vector<MoleculeType::Ptr>& constituents) :
  id_(num_molecule_types_),
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
  oss << "  id_: " << id_ << endl;
  oss << "  symbol_: " << symbol_ << endl;
  oss << "  daltons_: " << daltons_ << endl;
  return oss.str();
}

class Biome : public Printable
{
public:
  typedef std::shared_ptr<Biome> Ptr;
  typedef std::shared_ptr<const Biome> ConstPtr;

  std::vector<size_t> flux_;  // Aligned with MoleculeType id_ values.
  std::string name_;
  
  Biome(const std::string& name);
  std::string _str() const;
};

Biome::Biome(const std::string& name): name_(name)
{
  flux_ = vector<size_t>(MoleculeType::numMoleculeTypes(), 0);
}

std::string Biome::_str() const
{
  std::ostringstream oss;
  oss << "Biome \"" << name_ << "\"" << endl;
  return oss.str();
}


// Contains the whole simulation model
class Prokaryotic
{
public:
  Prokaryotic();
  std::vector<MoleculeType::Ptr> molecule_types_;
  std::map<std::string, MoleculeType::Ptr> molecule_map_;

  void addMoleculeType(MoleculeType::Ptr mt);
  MoleculeType::Ptr molecule(const std::string& name) {
    assert(molecule_map_.find(name) != molecule_map_.end());
    return molecule_map_[name];
  } 
};

Prokaryotic::Prokaryotic()
{
  // addMoleculeType(MoleculeType::Ptr(new MoleculeType("ATP", ":bang:", 503.15)));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("ADP", ":briefcase:", 423.17)));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("Phosphate", "P", 94.97)));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("", "R", 1000)));
  std::vector<MoleculeType::Ptr> constituents;
  constituents.push_back(molecule("ADP"));
  constituents.push_back(molecule("Phosphate"));
  addMoleculeType(MoleculeType::Ptr(new MoleculeType("ATP", ":bang:", constituents)));

  for (auto mt : molecule_types_)
    cout << mt->str() << endl;
}

void Prokaryotic::addMoleculeType(MoleculeType::Ptr mt)
{
  molecule_types_.push_back(mt);
  molecule_map_[mt->name_] = mt;
}


int main(int argc, char** argv)
{
  // MoleculeType e0("ATP", "bang", 503.15);
  // MoleculeType e1("ADP", "wallet", 423.17);
  // MoleculeType e2("Phosphate", "", 94.97);
  // //cout << e1.id() << endl;

  // cout << e0.str() << endl;

  // Biome b0("Alkaline vents");
  // //cout << b0.flux_ << endl;
  // cout << b0.str() << endl;

  Prokaryotic prokaryotic;
  
  return 0;
}
