#pragma once

#include <limits>
#include <Eigen/Dense>
#include <cassert>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <yaml-cpp/yaml.h>

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
  MoleculeVals normalized() const;
  
  std::string _str() const;


  // Probably should just always use .vals_ instead of these.
  const double& operator[](int idx) const { return vals_.coeffRef(idx); }
  double& operator[](int idx) { return vals_.coeffRef(idx); }
  
  // Provide string access though, for easy setup tasks.
  const double& operator[](const std::string& name) const;
  double& operator[](const std::string& name);
  bool hasMolecule(const std::string& name) const;
};


class MoleculeMat : public Printable
{
public:
  const Prokaryotic& pro_;
  Eigen::ArrayXXd vals_;
  
  MoleculeMat(const Prokaryotic& pro);  
  std::string _str() const;
  MoleculeVals row(const std::string& name) const;
  
  // // Provide string access for easy setup tasks.
  // // (Don't use where performance is sensitive.)
  // double& operator[](const std::string& name0, const std::string& name1) const;
};


// See Fig 1 of http://book.bionumbers.org/how-many-reactions-do-enzymes-carry-out-each-second/
// substrate_concentration in mM, km in mM, kcat in reactions / sec
// Also https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics
double rateMM(double substrate_concentration, double km, double kcat);

double probabilityPerSecond(double half_life_hours)
{
  if (half_life_hours == std::numeric_limits<double>::max())
    return 0.0;

  double half_life_ticks = half_life_hours * 3600;
  return 1.0 - exp(log(0.5) / half_life_ticks);
}

MoleculeVals countsToConcentrations(const MoleculeVals& counts, double um3);
double countToConcentration(double count, double um3);

// Superclass for proteins and small molecules, both of which we will just call "Molecules".
class MoleculeType : public Printable
{
public:
  typedef std::shared_ptr<MoleculeType> Ptr;
  typedef std::shared_ptr<const MoleculeType> ConstPtr;
  
  int idx_;  // probably this should be managed strictly in MoleculeTable.
  std::string name_;
  std::string symbol_;
  double daltons_;
  double half_life_hours_;
  double num_amino_acids_;

  MoleculeType(const Prokaryotic& pro, const std::string& name, const std::string& symbol,
               double daltons,
               double num_amino_acids = 0,
               double half_life_hours = std::numeric_limits<double>::max());
  MoleculeType(const Prokaryotic& pro, const YAML::Node& yaml);
  
  std::string _str() const;
  double pDenature() const { return probabilityPerSecond(half_life_hours_); }
};

class ReactionType : public Printable
{
public:
  typedef std::shared_ptr<const ReactionType> ConstPtr;
  typedef std::shared_ptr<ReactionType> Ptr;

  const Prokaryotic& pro_;
  MoleculeVals inputs_;
  MoleculeVals outputs_;
  MoleculeVals kms_;  // concentrations that yield half the max rate.  Each element is mM.
  double kcat_;  // max rate per tick (1 tick is 1 second?)
  std::string protein_name_;
  size_t protein_idx_;

  ReactionType(const Prokaryotic& pro);
  ReactionType(const Prokaryotic& pro, const YAML::Node& yaml);
  ReactionType(const Prokaryotic& pro,
               const MoleculeVals& inputs, const MoleculeVals& outputs,
               const MoleculeVals& kms, double kcat);
  virtual ~ReactionType() {}
  // For ribosomal protein synthesis reactions, we tell it what protein we are synthesizing, and it
  // assigns the inputs_, outputs_, kms_, kcat_, protein_name_, protein_idx_.
  ReactionType(const Prokaryotic& pro, MoleculeType::ConstPtr protein);  
  virtual void tick(Cell& cell, int num_protein_copies) const;
  std::string _str() const;

private:
  void parseHalfFormula(const std::vector<std::string>& tokens, size_t startidx, size_t endidx, MoleculeVals* mvals);
  void parseFormula(const std::string& formula);
};

// Unfortunately needs special code because they target denatured proteins.
class ProteasomeReactionType : public ReactionType
{
public:
  typedef std::shared_ptr<const ProteasomeReactionType> ConstPtr;
  typedef std::shared_ptr<ProteasomeReactionType> Ptr;

  size_t target_idx_;
  size_t atp_idx_;
  double input_denatured_;
  
  ProteasomeReactionType(const Prokaryotic& pro, size_t target_idx);
  void tick(Cell& cell, int num_protein_copies) const;  
};

class DNA;

class DNAThen
{
public:
  typedef std::shared_ptr<DNAThen> Ptr;
  typedef std::shared_ptr<const DNAThen> ConstPtr;
  
  const Prokaryotic& pro_;
  Cell& cell_;
  std::string molecule_name_;
  int molecule_idx_;
  std::string operation_type_;
  double value_;
  
  DNAThen(const Prokaryotic& pro, Cell& cell, const std::string& thenstr);
  void apply(DNA* dna) const;
};

class DNAIf
{
public:
  typedef std::shared_ptr<DNAIf> Ptr;
  typedef std::shared_ptr<const DNAIf> ConstPtr;

  const Prokaryotic& pro_;
  Cell& cell_;
  std::string molecule_name_;
  int molecule_idx_;
  std::string inequality_type_;
  double threshold_;
  
  std::vector<DNAIf::Ptr> subifs_;
  std::vector<DNAThen::Ptr> thens_;
  
  DNAIf(const Prokaryotic& pro, Cell& cell, const std::string& ifstr);
  DNAIf(const Prokaryotic& pro, Cell& cell, const YAML::Node& yaml);
  void initializeFromString(const std::string& ifstr);

  void execute(Cell* cell) const;
  bool check(const Cell& cell) const;
};


// Contains the high level scripts for the cell.
// e.g. high level functions like divide() are here.
// This is where cell division happens based on some conditions, and protein synthesis regulation.
// We're going to gloss over transcriptional vs translational regulation for now and just combine them into one big thing.
class DNA : public Printable
{
public:
  typedef std::shared_ptr<const DNA> ConstPtr;
  typedef std::shared_ptr<DNA> Ptr;

  const Prokaryotic& pro_;
  // transcription_factors_ adjust the distribution of ribosomes to genes.
  // User can set write simple scripts to change the transcription_factors_ in response to various things
  // If the sum is zero, no ribosomes will do anything.
  // If the sum is greater than 1, the distribution will be normalized.
  MoleculeVals transcription_factors_;
  // This array contains the number of ribosomes assigned to each protein synthesis task.
  // ribosome_assignments_.vals_.sum() <= cytosol_contents_["Ribosome"].
  MoleculeVals ribosome_assignments_;
  std::vector<ReactionType::ConstPtr> synthesis_reactions_;
  std::vector<DNAIf::ConstPtr> dna_ifs_;
  // How many copies of each gene you have.  This affects ribosomal parallelization.
  //MoleculeVals gene_copy_numbers_;
  
  DNA(const Prokaryotic& pro);
  void tick(Cell& cell);
  std::string _str() const { return ""; }
  void setDefaultTranscriptionFactors();
  void assignRibosomes(Cell* cell);
};

// Collects stats on the Cell as it does its thing.
class CellObserver
{
public:
  const Prokaryotic& pro_;
  // On average, when X is consumed, what are the products of that reaction?
  MoleculeMat transformation_flux_;
  MoleculeMat transformation_flux_nr_;  // num_reactions
  MoleculeMat protein_io_flux_;
  MoleculeVals protein_synth_;
  MoleculeVals protein_den_;  // denaturing
  MoleculeVals proteasome_action_;  // denatured protein -> amino acids
  std::vector<double> num_ticks_per_division_period_;  // in ticks
  int num_ticks_since_last_division_;
  std::vector<Eigen::ArrayXd> cytosol_contents_history_;
  
  CellObserver(const Prokaryotic& pro);
  void recordReactionFlux(const MoleculeVals& flux, int protein_idx);
  void recordProteinSynthAndDen(const MoleculeVals& flux);
  void recordProteasomeAction(int target_idx_, double num_to_remove);
  void recordDivision();
  void recordCytosolContents(const MoleculeVals& cytosol_contents);
  void tick();
  double averageDivisionHours() const;
  std::string formatTransformationFlux(const std::string& prefix = "") const;
  std::string formatProteinIOFlux(const std::string& prefix = "") const;
  std::string formatProteinStateChanges(const std::string& prefix = "") const;
  std::string formatDivisionHours(const std::string& prefix = "") const;
  std::string formatCytosolContentsHistoryAvg(const std::string& prefix = "") const;

  MoleculeVals cytosolContentsHistoryAvg() const;
};

class Cell : public Printable
{
public:
  typedef std::shared_ptr<Cell> Ptr;
  typedef std::shared_ptr<const Cell> ConstPtr;

  const Prokaryotic& pro_;
  std::string name_;
  double um3_;  // 1e-18 m3, or 1e-15 L (i.e. femtoliter)
  MoleculeVals cytosol_contents_;
  MoleculeVals cytosol_contents_denatured_;
  MoleculeVals membrane_contents_;
  MoleculeVals membrane_permeabilities_;
  DNA::Ptr dna_;  // non-const so Prokaryotic can make changes to it.
  std::vector<ProteasomeReactionType::ConstPtr> proteasome_reactions_;
  CellObserver obs_;
  
  Cell(const Prokaryotic& pro, const std::string& name);

  std::string _str() const;
  void tick(const Biome& biome);
  MoleculeVals cytosolConcentrations() const;  // mM
  static MoleculeVals cytosolConcentrations(const Prokaryotic& pro, MoleculeVals cytosol_contents, double um3);  // mM
  // concentrations is mM
  void setCytosolContentsByConcentrations(const MoleculeVals& cytosol_concentrations);
  // concentrations is mM
  static MoleculeVals cytosolContents(const Prokaryotic& pro, const MoleculeVals& cytosol_concentrations, double um3);
  void addDNAIf(const YAML::Node& yaml);
  void applyReactionResult(const MoleculeVals& flux, int protein_idx);
  void divide();
};

class Biome : public Printable
{
public:
  typedef std::shared_ptr<Biome> Ptr;
  typedef std::shared_ptr<const Biome> ConstPtr;

  const Prokaryotic& pro_;
  double m3_;
  MoleculeVals concentrations_;  // mM
  std::string name_;
  // Once we have multiple biomes, these will need to be deep copies.
  // For now we'll just use the cells_ in pro_.
  std::vector<Cell::Ptr> cells_;
  Eigen::VectorXd populations_;
  
  Biome(const Prokaryotic& pro, double m3, const std::string& name);
  Biome(const Prokaryotic& pro, const YAML::Node& yaml);
  std::string _str() const;
  void tick();
  void step();  
};

// Contains the whole simulation model
// Data protection strategy: Data is stored as (non-const) Ptr, but is only handed out as ConstPtr.
class Prokaryotic : public Printable
{
public:
  Prokaryotic();

  void addReactionType(ReactionType::Ptr rt) { reaction_types_.push_back(rt); }
  void addReactionType(const YAML::Node& yaml);
  void addMoleculeType(MoleculeType::Ptr mt);
  void addMoleculeType(const YAML::Node& yaml);
  void addBiome(Biome::Ptr biome) { biomes_.push_back(biome); }
  void addBiome(const YAML::Node& yaml);
  void applyConfig(const std::string& path);
  void applyConfig(const YAML::Node& yaml2);
  
  MoleculeType::ConstPtr molecule(const std::string& name) const {
    assert(molecule_map_.find(name) != molecule_map_.end());
    return molecule_map_.at(name);
  }
  MoleculeType::ConstPtr molecule(size_t idx) const { return molecule_types_[idx]; }
  
  size_t moleculeIdx(const std::string& name) const;
  const std::string& moleculeName(size_t idx) const { return molecule_types_[idx]->name_; }
  std::vector<MoleculeType::ConstPtr> moleculeTypes() const;
  size_t numMoleculeTypes() const { return molecule_types_.size(); }
  bool hasMolecule(const std::string& name) const { return molecule_map_.find(name) != molecule_map_.end(); }
  MoleculeVals numAA() const;
  
  void tick();
  void step();
  std::string _str() const;
  void run();

//private:
  std::vector<MoleculeType::Ptr> molecule_types_;
  std::map<std::string, MoleculeType::Ptr> molecule_map_;
  std::vector<Biome::Ptr> biomes_;
  std::vector<Cell::Ptr> cells_;

  std::vector<ReactionType::Ptr> reaction_types_;
  
private:
  MoleculeType::Ptr _molecule(const std::string& name) const {
    assert(hasMolecule(name));
    return molecule_map_.at(name);
  }
  MoleculeType::Ptr _molecule(size_t idx) const {
    return molecule_types_[idx];
  }
};
