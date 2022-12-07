#include <cstdint>
#include <random>
#include <fmt/ranges.h>
#include <fmt/core.h>
#include <prokaryotic.h>
#include <chrono>
#include <thread>
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/Dense>


using std::shared_ptr;
using std::cout, std::endl;
using std::vector;
using std::string;
using Eigen::ArrayXd;

void PASSERT(bool flag, const std::string& msg)
{
  if (!flag) {
    cout << "[PASSERT] " << msg << endl;
    assert(false);
  }
}


vector<string> tokenizeSimple(const std::string& s)
{
  vector<string> tokens;
  std::stringstream ss(s);
  string word;
  while (ss >> word)
    tokens.push_back(word);
  return tokens;
}

bool is_number(const std::string& s)
{
  std::string::const_iterator it = s.begin();
  while (it != s.end() && std::isdigit(*it)) ++it;
  return !s.empty() && it == s.end();
}

void ReactionType::parseHalfFormula(const std::vector<std::string>& tokens, size_t startidx, size_t endidx, MoleculeVals* mvals)
{
  double coeff = 1.0;
  string molecule_name = "";
  for (size_t i = startidx; i < endidx; ++i) {
    if (is_number(tokens[i]))
      coeff = atof(tokens[i].c_str());
    else if (tokens[i] == "+")
      continue;
    else {
      (*mvals)[tokens[i]] = coeff;
      coeff = 1.0;
      molecule_name = "";
    }
  }
}

void ReactionType::parseFormula(const std::string& formula)
{
  vector<string> tokens = tokenizeSimple(formula);
  int sep = -1;
  for (size_t i = 0; i < tokens.size(); ++i) {
    if (tokens[i] == "->") {
      parseHalfFormula(tokens, 0, i, &inputs_);
      sep = i;
    }
  }
  parseHalfFormula(tokens, sep+1, tokens.size(), &outputs_);
}

ReactionType::ReactionType(const Prokaryotic& pro, const YAML::Node& yaml) :
  pro_(pro),
  inputs_(pro),
  outputs_(pro),
  kms_(pro)
{
  string formula = yaml["formula"].as<string>();

  // Eventually: a ReactionType will not have an associated protein_idx_, instead there will be a matrix
  // of factors that gets applied, and the catalyzing protein will have the dominating effect.
  // But this will do for now.
  if (yaml["protein"]) {
    protein_name_ = yaml["protein"].as<string>();
    protein_idx_ = pro_.moleculeIdx(protein_name_);
  }
  
  parseFormula(formula);
  
  kcat_ = yaml["kcat"].as<double>();
  
  for (YAML::const_iterator kmit = yaml["KMs"].begin(); kmit != yaml["KMs"].end(); ++kmit) {
    string molecule_name = kmit->first.as<string>();
    double km = kmit->second.as<double>();
    kms_[molecule_name] = km;
  }    
}

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
  assert(num_protein_copies > 0);
  
  // Simplify: Assume that everything follows Michaelis-Menten kinetics, even for the (presumably majority) of reactions
  // that have multiple substrates.
  // We'll just take the min reaction rate across the substrates (inputs).
  // https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics
  MoleculeVals concentrations = cell.cytosolConcentrations();
  double minrate = std::numeric_limits<double>::max();
  for (int i = 0; i < inputs_.vals_.size(); ++i)
    if (inputs_[i] > 0)
      minrate = std::min(minrate, rateMM(concentrations[i], kms_[i], kcat_));
  double rate = minrate;  // reactions / tick (1 tick == 1 second?)
  assert(rate >= 0);
  
  // update cell.cytosol_counts_
  assert((inputs_.vals_ >= 0.0).all());
  assert((outputs_.vals_ >= 0.0).all());
  cell.cytosol_contents_.vals_ -= inputs_.vals_ * rate * num_protein_copies;
  cell.cytosol_contents_.vals_ += outputs_.vals_ * rate * num_protein_copies;

  // It's possible for us to accidentally overstep our bounds here and end up making something go negative.
  // If that happens, print a warning but continue.
  for (int i = 0; i < cell.cytosol_contents_.vals_.size(); ++i) {
    if (cell.cytosol_contents_[i] < 0) {
      cout << "WARNING: " << pro_.molecule(i)->name_ << " went negative, to " << cell.cytosol_contents_[i]
           << ".  Clipping it back to zero." << endl;
      cell.cytosol_contents_.vals_[i] = 0;
    }
  }
  
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

ReactionType::ReactionType(const Prokaryotic& pro, MoleculeType::ConstPtr protein_to_synthesize) :
  pro_(pro),
  inputs_(pro),
  outputs_(pro),
  kms_(pro)
{
  protein_name_ = "Ribosome";
  protein_idx_ = pro_.moleculeIdx(protein_name_);
  
  // https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=107782
  // ~5 ATP per amino acid added to a protein.
  inputs_["ATP"] = 5;
  inputs_["Amino acids"] = 1;
  // We're trying out a probabilistic average synthesis model here.
  outputs_[protein_to_synthesize->name_] = 1.0 / protein_to_synthesize->num_amino_acids_;
  outputs_["ADP"] = 5;
  outputs_["Phosphate"] = 5;

  // I think for now we are just leaving these at some fixed reasonable number?
  // We're aiming for 20 AAs / second.  https://micro.magnet.fsu.edu/cells/ribosomes/ribosomes.html
  // Normal cellular ATP is 1-10 mM.  (ADP is typically 1000x lower.)
  kms_["ATP"] = 0.5;  // Probably for typical cells with 1-10 mM ATP, the ATP concentration is not the limiting factor.
  // http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/
  // glutamate 96 mM
  // aspartate 4.2 mM, valine 4.0, glutamine 3.8, alanine 2.5, arginine 0.57, asparagine 0.51, lysine 0.4, proline 0.38, methionine 0.14, ...
  // Ok so total AA concentration is maybe like 110 mM.  This is e. coli.
  // https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=2&id=107764
  // Hm this says 150 mM in yeast.  Ok, good enough.  We'll set KM at around half of e. coli I guess.
  kms_["Amino acids"] = 60;
  kcat_ = 20;
}

double rateMM(double substrate_concentration, double km, double kcat)
{
  assert(substrate_concentration >= 0);
  assert(km > 0);
  assert(kcat > 0);
  return kcat * substrate_concentration / (km + substrate_concentration);
}

MoleculeType::MoleculeType(const Prokaryotic& pro,
                           const std::string& name,
                           const std::string& symbol,
                           double daltons,
                           double num_amino_acids, double half_life_hours) :
  idx_(pro.numMoleculeTypes()),
  name_(name),
  symbol_(symbol),
  daltons_(daltons),
  half_life_hours_(half_life_hours),
  num_amino_acids_(num_amino_acids)
{
}

MoleculeType::MoleculeType(const Prokaryotic& pro, const YAML::Node& yaml) :
  idx_(pro.numMoleculeTypes()),
  name_(yaml["name"].as<string>()),
  symbol_(yaml["symbol"].as<string>()),
  daltons_(yaml["daltons"].as<double>()),
  half_life_hours_(std::numeric_limits<double>::max()),
  num_amino_acids_(0)
{
  if (yaml["half-life-hours"])
    half_life_hours_ = yaml["half-life-hours"].as<double>();
  if (yaml["num-amino-acids"])
    num_amino_acids_ = yaml["num-amino-acids"].as<double>();
}
                           
std::string MoleculeType::_str() const
{
  std::ostringstream oss;
  oss << "MoleculeType \"" << name_ << "\"" << endl;
  oss << "  idx_: " << idx_ << endl;
  oss << "  symbol_: " << symbol_ << endl;
  oss << "  daltons_: " << daltons_ << endl;
  oss << "  half_life_hours_: " << half_life_hours_ << endl;
  oss << "  num_amino_acids_: " << num_amino_acids_ << endl;
  
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

bool MoleculeVals::hasMolecule(const std::string& name) const
{
  return pro_.hasMolecule(name);
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
    if (vals_[i] - uint64_t(vals_[i]) > thresh)
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
  cytosol_contents_denatured_(pro_),
  membrane_contents_(pro_),
  membrane_permeabilities_(pro_),
  dna_(new DNA(pro_))
{
}

std::string Cell::_str() const
{
  std::ostringstream oss;
  oss << "Cell \"" << name_ << "\"" << endl;
  oss << "  um3_: " << um3_ << endl;
  oss << "  cytosol_contents_: " << endl << cytosol_contents_.str("    ") << endl;
  oss << "  cytosol_contents_denatured_: " << endl << cytosol_contents_denatured_.str("    ") << endl;
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
  assert((concentrations.vals_ >= 0.0).all());
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
  for (auto rt : pro_.reaction_types_) {
    rt->tick(*this, cytosol_contents_[rt->protein_idx_]);
  }
  
  // Apply degredation of molecules.
  // TODO: Have a member of MoleculeType (or, better, ProteinType) that indicates if the protein is denatured or not.
  // Then add proteasomes that identify denatured proteins (through ubiquitin, though maybe we don't want that ... or maybe we do??)
  // Yeah just do proteasomes for now.
  for (auto mt : pro_.moleculeTypes()) {
    if (mt->pDenature() > 0) {
      double num_to_denature = cytosol_contents_[mt->idx_] * mt->pDenature();
      cytosol_contents_[mt->idx_] -= num_to_denature;
      cytosol_contents_denatured_[mt->idx_] += num_to_denature;
    }
  }

  // Apply "DNA programming"
  dna_->tick(*this);
  
  cytosol_contents_.probabilisticRound();
}

DNAIf::DNAIf(const Prokaryotic& pro, const YAML::Node& yaml) :
  pro_(pro)
{
  initializeFromString(yaml["if"].as<string>());
  for (const YAML::Node& thenstr : yaml["then"]) {
    thens_.push_back(DNAThen::Ptr(new DNAThen(pro_, thenstr.as<string>())));
  }
}

DNAIf::DNAIf(const Prokaryotic& pro, const std::string& ifstr) :
  pro_(pro)
{
  initializeFromString(ifstr);
}

void DNAIf::initializeFromString(const std::string& ifstr)
{
  vector<string> tokens = tokenizeSimple(ifstr);
  
  for (size_t i = 0; i < tokens.size(); ++i) {
    if (tokens[i] == ">" || tokens[i] == "<") {
      molecule_name_ = fmt::format("{}", fmt::join(tokens.begin(), tokens.begin() + i, " "));
      molecule_idx_ = pro_.moleculeIdx(molecule_name_);
      inequality_type_ = tokens[i];
      assert(tokens.size() == i + 2);
      threshold_ = atof(tokens[i+1].c_str());
      break;
    }
  }  
}

void DNAIf::execute(Cell* cell) const
{
  if (!check(*cell))
    return;
  for (auto then : thens_)
    then->apply(cell->dna_.get());
  for (auto subif : subifs_)
    subif->execute(cell);
}

bool DNAIf::check(const Cell& cell) const
{
  MoleculeVals contents(pro_);
  contents.vals_ = cell.cytosol_contents_.vals_ + cell.membrane_contents_.vals_;
  
  if (inequality_type_ == "<")
    if (contents[molecule_idx_] < threshold_)
      return true;
  
  if (inequality_type_ == ">")
    if (contents[molecule_idx_] > threshold_)
      return true;
  
  return false;
}

DNAThen::DNAThen(const Prokaryotic& pro, const std::string& thenstr) :
  pro_(pro)
{
  vector<string> tokens = tokenizeSimple(thenstr);
  
  for (size_t i = 0; i < tokens.size(); ++i) {
    if (tokens[i] == "=" || tokens[i] == "+=" || tokens[i] == "-=") {
      molecule_name_ = fmt::format("{}", fmt::join(tokens.begin(), tokens.begin() + i, " "));
      if (molecule_name_ == "all")
        molecule_idx_ = -1;
      else
        molecule_idx_ = pro_.moleculeIdx(molecule_name_);
      
      operation_type_ = tokens[i];
      assert(operation_type_ == "=" || operation_type_ == "+=" || operation_type_ == "-=");
      assert(tokens.size() == i + 2);
      value_ = atof(tokens[i+1].c_str());
      break;
    }
  }  
}

void DNAThen::apply(DNA* dna) const
{
  if (molecule_idx_ >= 0) {
    if (operation_type_ == "=")
      dna->transcription_factors_[molecule_idx_] = value_;
    else if (operation_type_ == "+=")
      dna->transcription_factors_[molecule_idx_] += value_;
    else if (operation_type_ == "-=")
      dna->transcription_factors_[molecule_idx_] -= value_;
  }
  else {
    assert(molecule_name_ == "all");
    if (operation_type_ == "=")
      dna->transcription_factors_.vals_ = value_;
    else if (operation_type_ == "+=")
      dna->transcription_factors_.vals_ += value_;
    else if (operation_type_ == "-=")
      dna->transcription_factors_.vals_ -= value_;
    return;
  }
}


DNA::DNA(const Prokaryotic& pro) :
  pro_(pro),
  transcription_factors_(pro),
  ribosome_assignments_(pro),
  synthesis_reactions_(transcription_factors_.vals_.size(), ReactionType::ConstPtr(nullptr))
{
  setDefaultTranscriptionFactors();

  for (size_t i = 0; i < synthesis_reactions_.size(); ++i)
    if (pro_.molecule(i)->num_amino_acids_ > 0)
      synthesis_reactions_[i] = ReactionType::ConstPtr(new ReactionType(pro, pro_.molecule(i)));
}

void DNA::setDefaultTranscriptionFactors()
{
  transcription_factors_.vals_.setZero();
  for (int i = 0; i < transcription_factors_.vals_.size(); ++i)
    if (pro_.molecule(i)->num_amino_acids_ > 0)
      transcription_factors_.vals_[i] = 1;
}

void DNA::tick(Cell& cell)
{
  // In testing, sometimes we don't have ribosomes.
  // Eventually: get rid of this.
  if (!cell.cytosol_contents_.hasMolecule("Ribosome"))
    return;

  // Every tick, transcription factors are set to their defaults, then user code is executed to make changes.
  setDefaultTranscriptionFactors();
  for (auto dnaif : dna_ifs_) {
    dnaif->execute(&cell);
  }
  
  // if (cell.cytosol_contents_.hasMolecule("ATP Synthase") && cell.cytosol_contents_["ATP Synthase"] < 200) {
  //   // Run the reaction that generates ATP Synthase.
  //   // It's a reaction run by Ribosomes (just like other proteins) that depends on concentration of substrates (building blocks,
  //   // Approx 5 ATP for each amino acid in the sequence.
  //   // 200 amino acids per minute made by ribosomes.  This is ~3 / second.  Maybe this is how we set kcat for the ribosome, and then
  //   // the synthesis rate drops if there aren't enough building blocks around.  Not sure how to set the KM values though.
  //   // Typically 100-600 amino acids per protein.  (Some are crazy long, but maybe that is just eukaryotes..)
  //   // Probably we want to ignore the fact that it's 3aa / second and just absorb that into kcat though.
  //   // ReactionType for protein synthesis is just a ReactionType stored in DNA class (which includes ribosome count) with kcat
  //   // set by num_aa of the protein, and num_protein_copies is determined by the number of ribosomes assigned to this protein.
  //   // That in turn is set by a distribution over genes that the user can set, also in the DNA programming - you can set promotion
  //   // factors in R^+, one for each gene, and ribosomes get assigned to genes according to a normalized distribution over genes.
  //   transcription_factors_["ATP Synthase"] = 1.0;
  // }

  // Transcription factor normalization:
  //   * If there are negative numbers, replace them with zeros.
  //   * If the transcription factors sum to greater than one, make them sum to one.
  //   * Otherwise leave them alone.
  
  MoleculeVals ntfs(pro_);  // normalized transcription factors
  ntfs.vals_ = transcription_factors_.vals_.max(ArrayXd::Zero(transcription_factors_.vals_.size()));
  if (ntfs.vals_.sum() > 1.0)
    ntfs.vals_ /= ntfs.vals_.sum();

  // Assign free ribosomes.
  assert(ntfs.vals_.size() == synthesis_reactions_.size());
  int num_free_ribosomes = cell.cytosol_contents_["Ribosome"] - ribosome_assignments_.vals_.sum();
  assert(num_free_ribosomes >= 0);
  for (int i = 0; i < ntfs.vals_.size(); ++i) {
    if (ntfs.vals_[i] > 0) {
      assert(pro_.molecule(i)->num_amino_acids_ > 0);
      ribosome_assignments_.vals_[i] += floor(ntfs.vals_[i] * num_free_ribosomes);
    }
  }
  assert((ribosome_assignments_.vals_ >= 0).all());

  // Get random ordering to use for the synthesis reactions.
  static vector<int> random_indices;
  static std::random_device rd;
  static std::mt19937 g(rd());
  if (random_indices.empty())
    for (size_t i = 0; i < synthesis_reactions_.size(); ++i)
      random_indices.push_back(i);
  std::shuffle(random_indices.begin(), random_indices.end(), g);
  
  // Run protein synthesis reactions.  These are special reactions which consume only ATP and amino acids, and produce only proteins.
  cell.cytosol_contents_.probabilisticRound();
  MoleculeVals orig_cytosol_contents = cell.cytosol_contents_;
  for (int idx : random_indices) {
    if (synthesis_reactions_[idx] && ntfs[idx] > 0)
      synthesis_reactions_[idx]->tick(cell, ribosome_assignments_.vals_[idx]);
  }
 
  cout << "[DNA::tick] cytosol contents before pround: " << endl << cell.cytosol_contents_.str("    ") << endl;
  for (int i = 0; i < cell.cytosol_contents_.vals_.size(); ++i)
    assert(cell.cytosol_contents_.vals_[i] >= -1e-6);
  
  // Collapse probability distributions over new proteins.  
  cell.cytosol_contents_.probabilisticRound();
  for (int i = 0; i < cell.cytosol_contents_.vals_.size(); ++i)
    assert(fabs(cell.cytosol_contents_.vals_[i] - uint64_t(cell.cytosol_contents_.vals_[i])) < 1e-6);
  
  // If we produced any proteins, free up the corresponding ribosomes.  
  ArrayXd new_molecules = cell.cytosol_contents_.vals_ - orig_cytosol_contents.vals_;
  for (int i = 0; i < new_molecules.size(); ++i) {
    if (new_molecules[i] > 0) {
      // New molecules should always be proteins, ADP, or Phosphate
      MoleculeType::ConstPtr nmt = pro_.molecule(i);
      assert(nmt->num_amino_acids_ > 0 || nmt->name_ == "ADP" || nmt->name_ == "Phosphate");
      
      // Just for the proteins:
      if (nmt->num_amino_acids_ > 0) {
        // Num new proteins should always be approx an int.
        assert(fabs(new_molecules[i] - uint64_t(new_molecules[i])) < 1e-6);
        // Free up the ribosomes corresponding to the new proteins.
        // Very small proteins can be generated at a rate more than 1 protein per tick.
        // In this case, the num proteins generated will be greater than the number of ribosomes assigned.
        // So, free up the minimum of num ribosomes assigned & num proteins created.
        ribosome_assignments_.vals_[i] -= std::min(ribosome_assignments_.vals_[i], new_molecules[i]);
      }
    }
  }

  assert((ribosome_assignments_.vals_ >= 0).all());
}

Prokaryotic::Prokaryotic()
{}

// void Prokaryotic::initializeCore()
// {

// }

void Prokaryotic::initialize(const YAML::Node yaml)
{
  YAML::Node mt = yaml["MoleculeTable"];
  assert(mt);
  for (size_t i = 0; i < mt.size(); ++i) {
    //for (YAML::const_iterator it = mt.begin(); it != mt.end(); ++it) {
    addMoleculeType(MoleculeType::Ptr(new MoleculeType(*this, mt[i])));
  }
    
  YAML::Node rt = yaml["ReactionTable"];
  assert(rt);
  for (size_t i = 0; i < rt.size(); ++i) {
    addReactionType(ReactionType::Ptr(new ReactionType(*this, rt[i])));
  }
  
  cells_.push_back(Cell::Ptr(new Cell(*this, "aoeu")));
  cells_[0]->membrane_permeabilities_["R"] = 0.0000005;
  cells_[0]->membrane_permeabilities_["Phosphate"] = 0.1;
  cells_[0]->membrane_permeabilities_["X"] = 0.001;
  
  biomes_.push_back(Biome::Ptr(new Biome(*this, 10, "Alkaline vents")));
  biomes_[0]->concentrations_["Phosphate"] = 0.01;
  biomes_[0]->concentrations_["R"] = 0;
  biomes_[0]->concentrations_["X"] = 10;
  
  cells_[0]->setCytosolContentsByConcentrations(biomes_[0]->concentrations_);
  cells_[0]->cytosol_contents_["ADP"] = 1e6;
  cells_[0]->cytosol_contents_["ATP Synthase"] = 1e4;
  cells_[0]->cytosol_contents_["ATP Consumer"] = 1e3;

  // cout << str() << endl;
}

size_t Prokaryotic::moleculeIdx(const std::string& name) const
{
  if (!hasMolecule(name)) {
    cout << "\"" << name << "\" not found in MoleculeTable" << endl;
    // fmt::print("\"{}\" not found in MoleculeTable.", name);
    cout << "MoleculeTable: " << endl;
    for (size_t i = 0; i < molecule_types_.size(); ++i)
      cout << molecule_types_[i]->str("  ") << endl;
    assert(hasMolecule(name));
  }
  
  return molecule(name)->idx_;
}


std::string Prokaryotic::_str() const
{
  std::ostringstream oss;
  oss << "Simulation status" << endl;
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


