#pragma once


#include <limits>
#include <Eigen/Dense>
#include <cassert>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <string>

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
  int size() const { return vals_.size(); }
  
  std::string _str() const;

  // Probably should just always use .vals_ instead of these.
  const double& operator[](int idx) const { return vals_.coeffRef(idx); }
  double& operator[](int idx) { return vals_.coeffRef(idx); }
  
  // Provide string access though, for easy setup tasks.
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
