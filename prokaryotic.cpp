#include <iomanip>
#include <cstdint>
#include <algorithm>
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

using namespace std;
using Eigen::ArrayXd, Eigen::ArrayXXd, Eigen::VectorXd;

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
  kms_(pro),
  protein_idx_(-1),
  grad_(pro)
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
  kcat_(kcat),
  protein_idx_(-1),
  grad_(pro)
{
  for (int i = 0; i < inputs_.vals_.size(); ++i) {
    if (inputs_.vals_[i] > 0)
      assert(kms_.vals_[i] > 0);
  }
}

void ReactionType::computeGrad(const Cell& cell, int num_protein_copies)
{
  if (num_protein_copies == 0)
    grad_.vals_.setZero();
  
  // Simplify: Assume that everything follows Michaelis-Menten kinetics, even for the (presumably majority) of reactions
  // that have multiple substrates.
  // We'll just take the min reaction rate across the substrates (inputs).
  // https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics
  MoleculeVals concentrations = cell.cytosolConcentrations();
  double minrate = std::numeric_limits<double>::max();
  for (int i = 0; i < inputs_.vals_.size(); ++i)
    if (inputs_[i] > 0) {
      double r = rateMM(concentrations[i], kms_[i], kcat_);
      minrate = std::min(minrate, r);
      // rates[i] = r;
    }
  double rate = minrate;  // reactions / tick (1 tick == 1 second?)
  assert(rate >= 0);

  assert((inputs_.vals_ >= 0.0).all());
  assert((outputs_.vals_ >= 0.0).all());
  // cout << "Reaction running with rate " << rate << " and num_protein_copies " << num_protein_copies << endl;
  // cout << "Reaction by protein " << protein_name_ << " has " << cell.cytosol_contents_["ATP"] << " ATP available." << endl;
  // if (protein_name_ == "ADP Synthase") {
  //   cout << "=== ADP Synthase reaction debugging ===" << endl;
  //   cout << "  Rate: " << rate << endl;
  //   cout << "  Rates: " << rates.transpose() << endl;
  //   cout << "  Rates: " << endl << MoleculeVals(pro_, rates).str("    ") << endl;
  //   cout << "  Cytosol: " << endl << cell.cytosol_contents_.str("    ") << endl;
  //   cin.ignore();
  // }
  
  grad_.vals_ = (outputs_.vals_ - inputs_.vals_) * rate * num_protein_copies;

  // if (protein_name_ == "Ribosome") {
  //   cout << "Ribosome reaction flux before" << endl;
  //   cout << flux.str("  ") << endl;
  //   cout << "cytosol contents before" << endl;
  //   cout << cell.cytosol_contents_.str("  ") << endl;
  // }

  // Reaction scaling is now handled on a whole-cell level.
  // // It's not unlikely that tick() is not granular enough and we end up making something negative.
  // // In that case, shrink the flux so that minimum final cytosol_contents_ result will be zero.
  // double multiplier = 1.0;
  // for (int i = 0; i < cell.cytosol_contents_.vals_.size(); ++i)
  //   if (flux.vals_[i] < 0)
  //     multiplier = std::min(multiplier, -0.999 * cell.cytosol_contents_.vals_[i] / flux.vals_[i]);
  // // cout << flux.vals_.transpose() << endl;
  // flux.vals_ *= multiplier;
  // // cout << flux.vals_.transpose() << endl;
}

void ReactionType::apply(Cell& cell) const
{
  cell.applyReactionResult(grad_, protein_idx_);
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

// This ReactionType is only used for ribosomes.
ReactionType::ReactionType(const Prokaryotic& pro, MoleculeType::ConstPtr protein_to_synthesize) :
  pro_(pro),
  inputs_(pro),
  outputs_(pro),
  kms_(pro),
  protein_idx_(-1),
  grad_(pro)
{
  protein_name_ = "Ribosome";
  protein_idx_ = pro_.moleculeIdx(protein_name_);
  
  // https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=107782
  // ~5 ATP per amino acid added to a protein.
  inputs_["ATP"] = 5;
  inputs_["Amino acids"] = 1;
  outputs_[protein_to_synthesize->name_] = 1.0 / protein_to_synthesize->num_amino_acids_;
  outputs_["ADP"] = 5;
  outputs_["Phosphate"] = 5;

  // I think for now we are just leaving these at some fixed reasonable number?
  // We're aiming for 20 AAs / second.  https://micro.magnet.fsu.edu/cells/ribosomes/ribosomes.html
  // Normal cellular ATP is 1-10 mM.  (ADP is typically 1000x lower.)
  // Probably for typical cells with 1-10 mM ATP, the ATP concentration is not the limiting factor?
  // Dunno, we'll say if ATP concentration goes to half of the lower end of the range, ribosomes slow down by half.
  kms_["ATP"] = 0.5; 
  // http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/
  // glutamate 96 mM
  // aspartate 4.2 mM, valine 4.0, glutamine 3.8, alanine 2.5, arginine 0.57, asparagine 0.51, lysine 0.4, proline 0.38, methionine 0.14, ...
  // Ok so total AA concentration is maybe like 110 mM.  This is e. coli.
  // https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=2&id=107764
  // Hm this says 150 mM in yeast.  Ok, good enough.  We'll set KM at around half of typical e. coli AA concentration I guess.
  kms_["Amino acids"] = 60;
  kcat_ = 20;  // 20 AA/s.
}

ReactionType::ReactionType(const Prokaryotic& pro) :
  pro_(pro),
  inputs_(pro),
  outputs_(pro),
  kms_(pro),
  protein_idx_(-1),
  grad_(pro)
{
}

ProteasomeReactionType::ProteasomeReactionType(const Prokaryotic& pro, size_t target_idx) :
  ReactionType(pro),
  target_idx_(target_idx),
  target_name_(pro.moleculeName(target_idx)),
  num_denatured_to_remove_(0)
{
  protein_name_ = "Proteasome";
  protein_idx_ = pro.moleculeIdx("Proteasome");
  atp_idx_ = pro.moleculeIdx("ATP");
  
  double target_num_aa = pro.molecule(target_idx_)->num_amino_acids_;
  
  // http://book.bionumbers.org/how-fast-do-proteasomes-degrade-proteins/
  // 40 aa/s.
  // As with protein synthesis reactions, we'll just take the average here.
  // e.g. if the target protein is 80aa, then on average one proteasome consumes half a target protein per second.

  kcat_ = 40.0;  // This corresponds to the 40 aa/s.

  input_denatured_ = 1.0 / target_num_aa;
  // inputs_[target_idx_] = 1.0 / target_num_aa;
  
  // How much ATP does it take to consume one aa?
  // It takes 5 ATP to append a single aa in the ribosome.  Maybe it's cheaper than that?
  inputs_["ATP"] = 2;

  // TODO: What recovery fraction is reasonable?  Dunno.
  // Surely nature optimizes this to be as high as possible.
  outputs_["Amino acids"] = 0.9;
  outputs_["ADP"] = inputs_["ATP"];
  outputs_["Phosphate"] = inputs_["ATP"];
  
  // We're aiming for 40 AAs / second.  https://micro.magnet.fsu.edu/cells/ribosomes/ribosomes.html
  // Normal cellular ATP is 1-10 mM.
  // Probably for typical cells with 1-10 mM ATP, the ATP concentration is not the limiting factor.
  // Set the KM for ATP accordingly.  (But wtf do I know?)
  kms_["ATP"] = 0.5;

  // Ok now this one I really don't know.  At what concentration of denatured protein does the proteasome slow down by half?
  // Just making a wild guess here.  Maybe it's much higher because it's a protein-protein interaction, so they are both
  // slowly bumbling around, and less likely to hit than with say ATP which is much smaller and therefore faster moving.
  kms_[target_idx_] = 0.5;
}

void ProteasomeReactionType::computeGrad(const Cell& cell, int num_protein_copies)
{
  // Outside of this, I think in Cell::tick(), there is a special bit of code that determines how many of each
  // proteasome bump into which denatured molecule types.  This is probably not correct but we'll probably just
  // not worry about it for now, or ever.
  // (It probably should be that it's the concentration of all denatured proteins that matters for reaction rate.)

  if (num_protein_copies <= 0) {
    grad_.vals_.setZero();
    num_denatured_to_remove_ = 0;
    return;
  }
  
  //double conc_den = countsToConcentrations(cell.cytosol_contents_denatured_.vals_.sum(), cell.um3_);
  double conc_den = countToConcentration(cell.cytosol_contents_denatured_[target_idx_], cell.um3_);
  double conc_atp = countToConcentration(cell.cytosol_contents_[atp_idx_], cell.um3_);

  // Our inputs are fixed to denatured proteins and ATP, so we'll just do that more explicitly
  // than the general for loop in ReactionType::tick().
  double rate_den = rateMM(conc_den, kms_[target_idx_], kcat_);
  double rate_atp = rateMM(conc_atp, kms_[atp_idx_], kcat_);
  double rate = std::min(rate_den, rate_atp);
  assert(rate >= 0);
  
  assert((inputs_.vals_ >= 0.0).all());
  assert((outputs_.vals_ >= 0.0).all());

  //cout << "ProteasomeReactionType::tick running with " << num_protein_copies << " proteasome copies, rate " << rate << endl;
  
  // For just the target protein input: Subtract from cytosol_counts_denatured_.
  // Everything else will apply to cytosol_contents_ via applyReactionResult().
  num_denatured_to_remove_ = input_denatured_ * rate * num_protein_copies;
  grad_.vals_ = (outputs_.vals_ - inputs_.vals_) * rate * num_protein_copies;

  // Competition between reactions is handled on a whole-cell level but we do a first scaling here to ensure that
  // we aren't trying to make denatured proteins go negative.
  // In general, eventually, all this should be handled as one big state vector of cell contents
  // and all transformations of cell state are just gradient vectors, and there are no special cases.
  // Making "denatured" be its own special state rather than just another dimension in MoleculeVals was probably a mistake.
  double scale = std::min(1.0, 0.999 * cell.cytosol_contents_denatured_.vals_[target_idx_] / num_denatured_to_remove_);
  if (scale < 1.0) {
    cout << "Proteasome first stage scaling: " << scale << endl;
    num_denatured_to_remove_ *= scale;
    grad_.vals_ *= scale;
  }

  // cout << "after initial scaling: " << endl;
  // cout << "  proteasome target: " << target_name_ << endl;
  // cout << "  num_denatured_to_remove_: " << num_denatured_to_remove_ << endl;
  // cout << "  cell.cytosol_contents_denatured_.vals_[target_idx_] " << cell.cytosol_contents_denatured_.vals_[target_idx_] << endl;
  
  assert(num_denatured_to_remove_ <= cell.cytosol_contents_denatured_.vals_[target_idx_]);
}

void ProteasomeReactionType::apply(Cell& cell) const
{
  // Not sure where this edge case is coming from but we're going to throw out this bullshit way of managing denatured
  // proteins shortly anyway.
  double nr = num_denatured_to_remove_;
  if (num_denatured_to_remove_ < 1.0 && num_denatured_to_remove_ > 0.0 && cell.cytosol_contents_denatured_.vals_[target_idx_] == 0) {
    nr = 0;
  }
  else if (num_denatured_to_remove_ > cell.cytosol_contents_denatured_.vals_[target_idx_]) {
    cout << "proteasome target: " << target_name_ << endl;
    cout << "num_denatured_to_remove_: " << num_denatured_to_remove_ << endl;
    cout << "cell.cytosol_contents_denatured_.vals_[target_idx_] " << cell.cytosol_contents_denatured_.vals_[target_idx_] << endl;
    assert(false);
  }
  
  cell.cytosol_contents_denatured_.vals_[target_idx_] -= nr;
  cell.obs_.recordProteasomeAction(target_idx_, nr);
  cell.applyReactionResult(grad_, protein_idx_);

  cell.assertSanity();
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
  daltons_(0),
  half_life_hours_(std::numeric_limits<double>::max()),
  num_amino_acids_(0)
{
  if (yaml["symbol"])
    symbol_ = yaml["symbol"].as<string>();
  if (yaml["half-life-hours"])
    half_life_hours_ = yaml["half-life-hours"].as<double>();
  if (yaml["num-amino-acids"]) {
    num_amino_acids_ = yaml["num-amino-acids"].as<double>();
    daltons_ = num_amino_acids_ * 110;  // avg
  }
  if (yaml["daltons"])
    daltons_ = yaml["daltons"].as<double>();
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

MoleculeVals::MoleculeVals(const Prokaryotic& pro, const ArrayXd& vals) :
  pro_(pro),
  vals_(vals)
{}

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

MoleculeVals MoleculeVals::normalized() const
{
  MoleculeVals result(*this);
  assert(result.vals_.sum() > 0);
  result.vals_ /= result.vals_.sum();
  return result;
}

MoleculeMat::MoleculeMat(const Prokaryotic& pro) :
  pro_(pro)
{
  int num = pro.moleculeTypes().size();
  vals_ = ArrayXXd::Zero(num, num);
}

std::string MoleculeMat::_str() const
{
  assert(false);
  return "";
}

MoleculeVals MoleculeMat::row(const std::string& name) const
{
  MoleculeVals row(pro_);
  row.vals_ = vals_.row(pro_.moleculeIdx(name));
  return row;
}

// double& MoleculeMat::operator[](const std::string& name0, const std::string& name1) const
// {
//   return vals_.coeffRef(pro_.moleculeIdx(name0), pro_.moleculeIdx(name1));
// }


Biome::Biome(const Prokaryotic& pro, double m3, const std::string& name) :
  pro_(pro),
  m3_(m3),
  concentrations_(pro),
  name_(name),
  cells_(pro.cells_),
  populations_(VectorXd::Zero(cells_.size()))
{
}

Biome::Biome(const Prokaryotic& pro, const YAML::Node& yaml) :
  pro_(pro),
  concentrations_(pro),
  cells_(pro.cells_),
  populations_(VectorXd::Zero(cells_.size()))
{
  name_ = yaml["name"].as<string>();
  m3_ = yaml["m3"].as<double>();
  
  for (YAML::const_iterator cit = yaml["concentrations"].begin(); cit != yaml["concentrations"].end(); ++cit) {
    string molecule_name = cit->first.as<string>();
    concentrations_[molecule_name] = cit->second.as<double>();
  }  
  cout << "Loaded biome" << endl;
  cout << str("  ") << endl;
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
  for (Cell::Ptr cell : cells_) {
    cell->tick(*this);
  }
}

void Biome::step(Comms& comms)
{
  // Run the cell for a while.  Get a division rate.
  // Maybe we run the cell for a few hours, or a day, or something.

  // We should maybe call step() on the cells individually, and this is part of
  // what Cell::step() does.
  for (Cell::Ptr cell : cells_) {
    cell->obs_.num_ticks_per_division_period_.clear();
    cell->obs_.cytosol_contents_history_.clear();
    cell->obs_.cytosol_contents_denatured_history_.clear();
  }
  
  double num_hours = 24 * 5;
  int num_ticks = num_hours * 60 * 60;
  for (int i = 0; i < num_ticks; ++i) {
    if (i % (60*60) == 0) {
      cout << "." << std::flush;
      //cout << cells_[0]->cytosol_contents_.str("  ") << endl;
    }
    if (i % (num_ticks / 20) == 0) {
      MessageWrapper msg;
      msg.addField("step_progress", (double)i / (num_hours * 60 * 60));
      comms.broadcast(msg);
    }
    tick();
  }

  // Get the average growth rates.
  VectorXd avg_division_ticks = VectorXd::Zero(cells_.size());
  for (size_t i = 0; i < cells_.size(); ++i) {
    //cout << "Cell " << i << " average time to divide (hours): " << cells_[i]->obs_.averageDivisionHours() << endl;
    cout << "Cell " << i << " division periods (hours): " << cells_[i]->obs_.formatDivisionHours() << endl;
  }
  // Initial conditions are weird though.  We probably want to consider only the most recent few division times,
  // and once a sort of stable state has been reached.  e.g. at the start of a step, a cell may have more resources
  // than it should bc the cell census was from a previous step when more resources were available, or bc a mutation
  // has just been applied, etc.  We need to reach a steady state where division happens periodically, and report that
  // steady state division rate.

  // For now we are not actually increasing population, we are just trying to make the cell minimize its division time.
  
  
  // // Take a step in population growth based on average growth rate.
  // // I'm not sure how long of a time period this should be.  E. coli can double every
  // // 17 minutes, so going for 12 hours without adjusting the biome for population changes
  // // doesn't seem like that's going to work.
  // int ticks_per_step = 60 * 60;  // one hour
  // VectorXd num_doublings = ticks_per_step / avg_division_ticks;
  // for (int i = 0; i < populations_.size(); ++i)
  //   populations_[i] *= pow(2, num_doublings[i])

  // // What if we run the simulation until we get two division periods that are about the same length,
  // // and we declare that the doubling time, double that cell's population, and apply some amount of population growth
  // // to everyone else prorated by the amount of doubling they had?
}

CellObserver::CellObserver(const Prokaryotic& pro) :
  pro_(pro),
  transformation_flux_(pro),
  transformation_flux_nr_(pro),
  protein_io_flux_(pro),
  protein_synth_(pro),
  protein_den_(pro),
  proteasome_action_(pro),
  num_ticks_per_division_period_(0)
{
}

std::string CellObserver::formatProteinStateChanges(const std::string& prefix) const
{
  ostringstream oss;
  oss << prefix << "Protein synthesis: " << protein_synth_.vals_.transpose() << endl;
  oss << prefix << "Protein denaturing: " << protein_den_.vals_.transpose() << endl;
  oss << prefix << "Proteasome action: " << proteasome_action_.vals_.transpose() << endl;
  return oss.str();
}

MoleculeVals CellObserver::cytosolContentsHistoryAvg() const
{
  MoleculeVals avg(pro_);
  for (size_t i = 0; i < cytosol_contents_history_.size(); ++i)
    avg.vals_ += cytosol_contents_history_[i];
  avg.vals_ /= cytosol_contents_history_.size();
  return avg;
}

MoleculeVals CellObserver::cytosolContentsDenaturedHistoryAvg() const
{
  MoleculeVals avg(pro_);
  for (size_t i = 0; i < cytosol_contents_denatured_history_.size(); ++i)
    avg.vals_ += cytosol_contents_denatured_history_[i];
  avg.vals_ /= cytosol_contents_denatured_history_.size();
  return avg;
}

std::string CellObserver::formatCytosolContentsHistoryAvg(const std::string& prefix) const
{
  MoleculeVals avg = cytosolContentsHistoryAvg();
  ostringstream oss;
  oss << prefix << "Cytosol contents history average over the last " << cytosol_contents_history_.size() << " ticks." << endl;
  oss << avg.str(prefix);
  return oss.str();
}

std::string CellObserver::formatDivisionHours(const std::string& prefix) const
{
  ostringstream oss;
  oss << prefix;
  for (size_t i = 0; i < num_ticks_per_division_period_.size(); ++i) {
    // if (i > 0)
    //   oss << " ";
    oss << " " << num_ticks_per_division_period_[i] / (60 * 60);
  }
  return oss.str();
}

std::string CellObserver::formatProteinIOFlux(const std::string& prefix) const
{
  ostringstream oss;
  const ArrayXXd& mat = protein_io_flux_.vals_;
  
  int width = 0;
  for (int i = 0; i < mat.rows(); ++i)
    if (pro_.molecule(i)->num_amino_acids_ > 0)
      width = std::max<double>(width, pro_.molecule(i)->name_.size());

  for (int i = 0; i < mat.rows(); ++i)
    if (pro_.molecule(i)->num_amino_acids_ > 0)
      oss << prefix << "| " << std::setw(width) << pro_.moleculeName(i) << " | " << mat.row(i) << endl;
  
  return oss.str();
}


std::string oomStr2(double num)
{
  if (num < 1e-9)
    return "   -   ";
  else if (num > 1e-3 && num < 1e3) {
    return fmt::format("{:7.2g}", num);
  }
  else
    return fmt::format("{:2.1e}", num);
}

std::string oomStr(double num)
{
  if (num < 1e-9)
    return ".";
  else
    return fmt::format("{:d}", int(log10(num)));
}

std::string CellObserver::formatTransformationFluxMat(const std::string& prefix) const
{
  ostringstream oss;

  // normalize for num reactions that contributed to this flow from molecule i -> molecule j.
  // e.g. we want to see 1 ATP -> 1 ADP + 1 P, no matter how many different reactions use ATP in this way.
  ArrayXXd mat = transformation_flux_.vals_;
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      if (fabs(mat(i, j)) > 0) {
        assert((transformation_flux_nr_.vals_(i, j) > 0));
        mat(i, j) /= transformation_flux_nr_.vals_(i, j);
      }
    }
  }
 
  int col_width = 9;
  int max_num_chars = 0;
  vector<string> molecules_to_use;
  for (int i = 0; i < mat.rows(); ++i) {
    string display_name = pro_.molecule(i)->name_;
    // std::replace(display_name.begin(), display_name.end(), ' ', '');
    int pad_amt = col_width - display_name.size();
    if (pad_amt > 0)
      display_name.append(pad_amt, ' ');
    molecules_to_use.push_back(display_name);
    max_num_chars = std::max<double>(max_num_chars, display_name.size());
  }

  // header
  vector<string> header_lines(col_width, std::string(col_width * molecules_to_use.size(), ' '));
  // for (size_t i = 0; i < molecules_to_use.size(); ++i) {
  //   for (int j = 0; j < col_width; ++j) {
  //     header_lines[j][col_width * i + j] = molecules_to_use[i][j];
  //   }
  // }

  for (size_t row = 0; row < header_lines.size(); ++row) {
    for (size_t col = 0; col < molecules_to_use.size(); ++col) {
      cout << "row: " << row << " col: " << col << " stridx: " << col * col_width + row << " letter: " << molecules_to_use[col][row] << endl;
      header_lines[row][(col * col_width) + row] = molecules_to_use[col][row];
    }
  }

  oss << prefix << std::string(max_num_chars+4, ' ');
  for (size_t i = 0; i < header_lines[0].size(); ++i)
    oss << "-";
  oss << endl;
  
  for (const string& line : header_lines)
    oss << prefix << std::string(max_num_chars+4, ' ') << line << endl;
  
  oss << prefix << std::string(max_num_chars+4, ' ');
  for (size_t i = 0; i < header_lines[0].size(); ++i)
    oss << "-";
  oss << endl;
    
  for (int i = 0; i < mat.rows(); ++i) {
    if (mat.row(i).abs().sum() > 0) {
      oss << prefix << "| " << std::setw(max_num_chars) << molecules_to_use[i] << " |";
      for (int j = 0; j < mat.cols(); ++j) {
        oss << " " << std::setw(col_width - 2) << oomStr2(mat(i,j)) << " ";
      }
      oss << endl;
    }
  }
    
      
  return oss.str();
}

std::string CellObserver::formatTransformationFlux(const std::string& prefix) const
{
  ostringstream oss;
  // normalize for num reactions that contributed to this flow from molecule i -> molecule j.
  // e.g. we want to see 1 ATP -> 1 ADP + 1 P, no matter how many different reactions use ATP in this way.
  ArrayXXd mat = transformation_flux_.vals_;
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      if (fabs(mat(i, j)) > 0) {
        assert((transformation_flux_nr_.vals_(i, j) > 0));
        mat(i, j) /= transformation_flux_nr_.vals_(i, j);
      }
    }
  }

  int width = 0;
  for (int i = 0; i < mat.rows(); ++i)
    if (mat.row(i).abs().sum() > 0)
      width = std::max<double>(width, pro_.molecule(i)->name_.size());
  
  for (int i = 0; i < mat.rows(); ++i)
    if (mat.row(i).abs().sum() > 0)
      oss << prefix << "| " << std::setw(width) << pro_.moleculeName(i) << " | " << mat.row(i) << endl;
  
  return oss.str();  
}

double CellObserver::averageDivisionHours() const
{
  double avg = 0;
  for (size_t i = 0; i < num_ticks_per_division_period_.size(); ++i)
    avg += num_ticks_per_division_period_[i] / (60 * 60.);
  avg /= num_ticks_per_division_period_.size();
  return avg;
}

void CellObserver::tick()
{
  transformation_flux_.vals_.setZero();
  transformation_flux_nr_.vals_.setZero();
  protein_io_flux_.vals_.setZero();
  protein_synth_.vals_.setZero();
  protein_den_.vals_.setZero();
  proteasome_action_.vals_.setZero();
  
  num_ticks_since_last_division_++;
}

void CellObserver::recordProteasomeAction(int target_idx_, double num_to_remove)
{
  proteasome_action_[target_idx_] += num_to_remove;
}

void CellObserver::recordProteinSynthAndDen(const MoleculeVals& flux)
{
  for (int i = 0; i < flux.vals_.size(); ++i) {
    if (flux[i] > 0)
      protein_synth_[i] += flux[i];
    if (flux[i] < 0)
      protein_den_[i] -= flux[i];
  }
}

void CellObserver::recordCytosolContents(const MoleculeVals& cytosol_contents, const MoleculeVals& denatured)
{
  cytosol_contents_history_.push_back(cytosol_contents.vals_);
  cytosol_contents_denatured_history_.push_back(denatured.vals_);
}

void CellObserver::recordDivision()
{
  num_ticks_per_division_period_.push_back(num_ticks_since_last_division_);
  num_ticks_since_last_division_ = 0;
}

void CellObserver::recordReactionFlux(const MoleculeVals& flux, int protein_idx)
{
  // If we consumed i and produced j in this reaction, record the transformation (normalized by how much was consumed)
  for (int i = 0; i < flux.vals_.size(); ++i)
    if (flux.vals_[i] < 0)
      for (int j = 0; j < flux.vals_.size(); ++j)
        if (flux.vals_[j] > 0) {
          transformation_flux_.vals_(i, j) += -flux.vals_[j] / flux.vals_[i];
          transformation_flux_nr_.vals_(i, j) += 1.0;
        }

  // If it was a protein that did the reaction, record io for that protein type.
  if (protein_idx >= 0)
    protein_io_flux_.vals_.row(protein_idx) += flux.vals_;
}

Cell::Cell(const Prokaryotic& pro, const std::string& name) :           
  pro_(pro),
  name_(name),
  um3_(1.0),
  cytosol_contents_(pro_),
  cytosol_contents_denatured_(pro_),
  membrane_contents_(pro_),
  membrane_permeabilities_(pro_),
  dna_(new DNA(pro_)),
  obs_(pro)
{
  // some tests don't have Proteasomes defined, so only make these reactions
  // if that molecule type exists.
  if (pro_.hasMolecule("Proteasome"))
    for (size_t i = 0; i < pro_.numMoleculeTypes(); ++i)
      if (pro_.molecule(i)->num_amino_acids_ > 0)
        proteasome_reactions_.push_back(ProteasomeReactionType::Ptr(new ProteasomeReactionType(pro_, i)));
}

void Cell::addDNAIf(const YAML::Node& yaml)
{
  dna_->dna_ifs_.push_back(DNAIf::Ptr(new DNAIf(pro_, *this, yaml)));
}

void Cell::applyReactionResult(const MoleculeVals& flux, int protein_idx)
{
  // double atp_adp_balance = fabs(flux["ATP"] + flux["ADP"]);
  // double atp_p_balance = fabs(flux["ATP"] + flux["Phosphate"]);
  // if (atp_adp_balance > 1e-6)
  //   cout << "flux " << endl << flux.str("  ") << endl;
  // assert(atp_adp_balance < 1e-6); atp_adp_balance--; // -Wall -Werror
  // assert(atp_p_balance < 1e-6); atp_p_balance--;  // -Wall -Werror
  cytosol_contents_.vals_ += flux.vals_;
  obs_.recordReactionFlux(flux, protein_idx);
}

void Cell::divide()
{
  
  // Eventually: This should be determined by making lots of proteins that are necessary for dividing,
  // then simulate some amount of consumption of ATP by those proteins as they do their job.
  // But for now we want the bare minimum thing that looks like simulation of population growth.
  cout << "Dividing." << endl;
  cout << "before: " << endl;
  cout << obs_.cytosolContentsHistoryAvg().str("  ") << endl;
  cytosol_contents_.vals_ /= 2;

  // // Uh I have no way to produce ATP/ADP right now, shit.
  // cytosol_contents_["ATP"] *= 2;
  // cytosol_contents_["ADP"] *= 2;
  // cytosol_contents_["Phosphate"] *= 2;
  
  cytosol_contents_denatured_.vals_ /= 2;
  membrane_contents_.vals_ /= 2;
  obs_.recordDivision();
  cout << "after: " << endl;
  cout << obs_.cytosolContentsHistoryAvg().str("  ") << endl;
}

std::string Cell::_str() const
{
  std::ostringstream oss;
  oss << "Cell \"" << name_ << "\"" << endl;
  oss << "  um3_: " << um3_ << endl;
  oss << "  cytosol_contents_: " << endl << cytosol_contents_.str("    ") << endl;
  oss << "  cytosol_contents_denatured_: " << endl << cytosol_contents_denatured_.str("    ") << endl;
//  oss << "  membrane_contents_: " << endl << membrane_contents_.str("    ") << endl;
  return oss.str();
}

MoleculeVals countsToConcentrations(const MoleculeVals& counts, double um3)
{
  // millimoles = 1e3 * (num / 6e23)
  // liters = um3 * 1e-15, bc um3 is fL.
  // mM = 1e3 * (num / 6e23) / (um3 * 1e-15)
  // mM = 1.67e-06 * num / um3

  PASSERT((counts.vals_ >= 0.0).all(), fmt::format("{}", counts.vals_));
  MoleculeVals concentrations(counts.pro_);
  concentrations.vals_ = counts.vals_ * (1.67e-6 / um3);
  assert((concentrations.vals_ >= 0.0).all());
  return concentrations;
}

double countToConcentration(double count, double um3)
{
  return count * (1.67e-6 / um3);
}

MoleculeVals Cell::cytosolConcentrations(const Prokaryotic& pro, MoleculeVals cytosol_contents, double um3)
{
  return countsToConcentrations(cytosol_contents, um3);
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
  obs_.tick();
  obs_.recordCytosolContents(cytosol_contents_, cytosol_contents_denatured_);
  
  // Passive membrane permeations
  // Should make this depend on cell volume (once the volume changes..)
  // details here: http://book.bionumbers.org/what-is-the-permeability-of-the-cell-membrane/
  // basically flux = -p (cin - cout), the same model I guessed at originally but per unit surface area.
  // Eventually: tag metabolic intermediates with a charge to make them not permeate the membrane.
  // Positively charged molecules really do not want to go through.
  //cout << "Phosphate before membrane permeation: " << cytosol_contents_["Phosphate"] << endl;
  
  MoleculeVals cytosol_concentrations = cytosolConcentrations();
  auto deltas = (biome.concentrations_.vals_ - cytosol_concentrations.vals_);
  assert((membrane_permeabilities_.vals_ < 1.0 + 1e-6).all());
  auto step = membrane_permeabilities_.vals_ * deltas;
  cytosol_concentrations.vals_ += step;
  setCytosolContentsByConcentrations(cytosol_concentrations);
  //cout << "Phosphate after membrane permeation: " << cytosol_contents_["Phosphate"] << endl;
  
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

  // In-cell reactions.  For now we are assuming all reactions require a protein.
  for (auto rt : pro_.reaction_types_)
    rt->computeGrad(*this, cytosol_contents_[rt->protein_idx_]);
  
  // Apply proteasome reactions.  I'm skeptical of this method of assigning proteasomes to the different categories of
  // denatured proteins.  They should probably run at a rate limited by KM on concentration of all denatured proteins,
  // then statistically remove copies of degraded proteins according to the distribution of cytosol_contents_denatured_.
  // I guess that's an Eventually: ...
  if (cytosol_contents_denatured_.vals_.sum() > 0) {
    MoleculeVals denatured_distribution = cytosol_contents_denatured_.normalized();
    assert(fabs(denatured_distribution.vals_.sum() - 1) < 1e-6);
    for (ProteasomeReactionType::Ptr prt : proteasome_reactions_) {
      prt->computeGrad(*this, cytosol_contents_["Proteasome"] * denatured_distribution[prt->target_idx_]);
    }
  }
  
  // Apply "DNA programming"
  dna_->computeGrad(*this);

  // At this point, all reactions should have computed their instantaneous rate and stored it in ReactionType::rate_.
  // Now we'll apply these reaction rates to the cell contents.
  // For each input molecule, if applying all reaction grads would make it go negative, scale all reactions that use that input
  // back equally so that it does not go negative.
  // Eventually there really should just be one big state vector for the cell in like R^10k, and each reaction is just a gradient
  // that will be applied to it.  Denatured, phosphorylated, tagged with charge, ribosome occupied vs free, etc should all be in there.
  // Probably should have protein states that automatically make the cell state vector.
  // Anyway, for now ...

  // Look through all ReactionTypes and proportionally scale their grad_ values so nothing goes negative
  // when we apply them to the cell.
  scaleGradients();

  // Apply grads to cell contents.
  dna_->apply(*this);
  for (auto rt : pro_.reaction_types_)
    rt->apply(*this);
  for (auto prt : proteasome_reactions_)
    prt->apply(*this);
  
  // Check for bugs.
  assertSanity();
  
  // rounding
  cytosol_contents_.probabilisticRound();
  cytosol_contents_denatured_.probabilisticRound();
}

void Cell::assertSanity() const
{
  assertPositiveCytosolVals();
  assertPositiveDenaturedCytosolVals();
}

void Cell::assertPositiveCytosolVals() const
{
  for (int i = 0; i < cytosol_contents_.vals_.size(); ++i) {
    if (cytosol_contents_.vals_[i] < -1e-6) {
      cout << "ERROR: negative cytosol_contents_ value." << endl;
      cout << cytosol_contents_.str("  ") << endl;
      assert(false);
    }
  }
}

void Cell::assertPositiveDenaturedCytosolVals() const
{
  for (int i = 0; i < cytosol_contents_denatured_.vals_.size(); ++i) {
    if (cytosol_contents_denatured_.vals_[i] < -1e-6) {
      cout << "ERROR: negative cytosol_contents_denatured_ value." << endl;
      cout << cytosol_contents_denatured_.str("  ") << endl;
      assert(false);
    }
  }
}

void accumulateGradComponents(const MoleculeVals& grad, MoleculeVals* total_grad, MoleculeVals* ngc, MoleculeVals* pgc)
{
  assert(grad.size() == total_grad->size());
  assert(grad.size() == ngc->size());
  assert(grad.size() == pgc->size());
  
  total_grad->vals_ += grad.vals_;
  for (int i = 0; i < grad.size(); ++i) {
    if (grad[i] > 0)
      pgc->vals_[i] += grad[i];
    if (grad[i] < 0)
      ngc->vals_[i] += grad[i];
  }

  // ngc + pgc == total_grad
  assert(((ngc->vals_ + pgc->vals_ - total_grad->vals_) < 1e-6).all());
}

void Cell::scaleGradients()
{
  // Compute the total gradient (and negative & positive grad components) for everything in the cell.
  MoleculeVals total_grad(pro_);
  MoleculeVals ngc(pro_);  // negative grad components
  MoleculeVals pgc(pro_);  // positive grad components

  for (auto srt : dna_->synthesis_reactions_)
    if (srt)
      accumulateGradComponents(srt->grad_, &total_grad, &ngc, &pgc);
  for (auto prt : proteasome_reactions_)
    accumulateGradComponents(prt->grad_, &total_grad, &ngc, &pgc);
  for (auto rt : pro_.reaction_types_)
    accumulateGradComponents(rt->grad_, &total_grad, &ngc, &pgc);

  // cout << "============================================================" << endl;
  // cout << "cytosol_contents_" << endl << cytosol_contents_.str("  ") << endl;
  // cout << "ngc (before scaling)" << endl << ngc.str("  ") << endl;
  // cout << "pgc (before scaling)" << endl << pgc.str("  ") << endl;

  MoleculeVals scalefactors(pro_);
  scalefactors.vals_.setOnes();
  for (int i = 0; i < total_grad.vals_.size(); ++i) {
    if (ngc[i] == 0) continue;
    // Determine how much we need to scale it back to avoid going negative.
    // Small fudge factor to ensure we don't cross zero due to numerical discretization.
    // Ideally we could include the positive gradient components here, but that leads
    // to a system of polynomial equations that will be too slow to solve inside tick().
    // So, no reaction products will exist or be assumed to exist until after everything is scaled.
    double scale = 0.999 * (-cytosol_contents_[i] / ngc[i]);
    if (scale >= 1.0) continue;
    scalefactors.vals_[i] = scale;
  }
  // For each reaction, for each input, check the scale factor for that input molecule.
  // Apply the smallest one (but not both) to reaction grad_.
  scaleNegativeGradientComponents(scalefactors);

  // Check that we actually end up with a new total gradient that won't take us negative.
  MoleculeVals total_grad2(pro_);
  MoleculeVals ngc2(pro_);  // negative grad components
  MoleculeVals pgc2(pro_);  // positive grad components
  for (auto srt : dna_->synthesis_reactions_)
    if (srt)
      accumulateGradComponents(srt->grad_, &total_grad2, &ngc2, &pgc2);
  for (auto prt : proteasome_reactions_)
    accumulateGradComponents(prt->grad_, &total_grad2, &ngc2, &pgc2);
  for (auto rt : pro_.reaction_types_)
    accumulateGradComponents(rt->grad_, &total_grad2, &ngc2, &pgc2);

  // cout << "ngc2 (after scaling)" << endl << ngc2.str("  ") << endl;
  // cout << "pgc2 (after scaling)" << endl << pgc2.str("  ") << endl;
  // cout << "total_grad2 + cytosol_contents_: " << endl << MoleculeVals(pro_, total_grad2.vals_ + cytosol_contents_.vals_).str("  ") << endl;
  
  assert(((total_grad2.vals_ + cytosol_contents_.vals_) >= 0).all());
  for (int i = 0; i < cytosol_contents_.vals_.size(); ++i)
    assert(cytosol_contents_.vals_[i] >= -1e-6);
  
}

void Cell::scaleNegativeGradientComponents(const MoleculeVals& scalefactors)
{
  // For each reaction, for each input, check the scale factor for that input molecule.
  // Apply the smallest one (but not both) to reaction grad_.

  for (auto rt : pro_.reaction_types_) {
    double scale_to_apply = 1.0;
    for (int midx = 0; midx < rt->grad_.size(); ++midx)
      if (rt->grad_[midx] < 0)
        scale_to_apply = std::min(scale_to_apply, scalefactors[midx]);
    if (scale_to_apply < 1.0) {
      rt->grad_.vals_ *= scale_to_apply;
      //cout << "Scaled back " << rt->protein_name_ << " reaction speed by " << scale_to_apply << endl;
    }
  }

  for (auto prt : proteasome_reactions_) {
    double scale_to_apply = 1.0;
    for (int midx = 0; midx < prt->grad_.size(); ++midx)
      if (prt->grad_[midx] < 0)
        scale_to_apply = std::min(scale_to_apply, scalefactors[midx]);
    assert(scale_to_apply > 0.0);
    if (scale_to_apply < 1.0) {
      prt->grad_.vals_ *= scale_to_apply;
      prt->num_denatured_to_remove_ *= scale_to_apply;
      //cout << "Scaled back " << prt->protein_name_ << "-" << prt->target_name_ << " reaction speed by " << scale_to_apply << endl;
    }
  }

  for (auto srt : dna_->synthesis_reactions_) {
    if (!srt) continue;
    
    double scale_to_apply = 1.0;
    for (int midx = 0; midx < srt->grad_.size(); ++midx)
      if (srt->grad_[midx] < 0)
        scale_to_apply = std::min(scale_to_apply, scalefactors[midx]);
    if (scale_to_apply < 1.0) {
      srt->grad_.vals_ *= scale_to_apply;
      //cout << "Scaled back " << srt->protein_name_ << " reaction speed by " << scale_to_apply << endl;
    }
  }
  
}

DNAIf::DNAIf(const Prokaryotic& pro, Cell& cell, const YAML::Node& yaml) :
  pro_(pro),
  cell_(cell)
{
  initializeFromString(yaml["if"].as<string>());
  for (const YAML::Node& y : yaml["then"]) {
    if (y.Type() == YAML::NodeType::Map) {
      subifs_.push_back(DNAIf::Ptr(new DNAIf(pro_, cell_, y)));
    }
    else
      thens_.push_back(DNAThen::Ptr(new DNAThen(pro_, cell_, y.as<string>())));
  }
}

DNAIf::DNAIf(const Prokaryotic& pro, Cell& cell, const std::string& ifstr) :
  pro_(pro),
  cell_(cell)
{
  initializeFromString(ifstr);
}

void DNAIf::initializeFromString(const std::string& ifstr)
{
  cout << "DNAIf::initializeFromString using str: " << ifstr << endl;
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

DNAThen::DNAThen(const Prokaryotic& pro, Cell& cell, const std::string& thenstr) :
  pro_(pro),
  cell_(cell)
{
  vector<string> tokens = tokenizeSimple(thenstr);

  if (tokens.size() == 1 && tokens[0] == "divide()") {
    operation_type_ = "divide";
    return;
  }
  
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
  if (operation_type_ == "divide") {
    cell_.divide();
    return;
  }
  
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
  synthesis_reactions_(transcription_factors_.vals_.size(), ReactionType::Ptr(nullptr))
{
  setDefaultTranscriptionFactors();

  for (size_t i = 0; i < synthesis_reactions_.size(); ++i)
    if (pro_.molecule(i)->num_amino_acids_ > 0)
      synthesis_reactions_[i] = ReactionType::Ptr(new ReactionType(pro, pro_.molecule(i)));
}

void DNA::setDefaultTranscriptionFactors()
{
  transcription_factors_.vals_.setZero();
  for (int i = 0; i < transcription_factors_.vals_.size(); ++i)
    if (pro_.molecule(i)->num_amino_acids_ > 0)
      transcription_factors_.vals_[i] = 1;
}

void DNA::assignRibosomes(Cell* cell)
{
  // Compute how many net ribosomes have been lost to degradation (min zero).
  double net_lost_ribosomes = 0;
  double num_rib = floor(cell->cytosol_contents_["Ribosome"]);
  
  static double last_num_rib = -1;
  if (last_num_rib > -1) {
    net_lost_ribosomes = std::max<double>(0, last_num_rib - num_rib);
  }
  last_num_rib = floor(cell->cytosol_contents_["Ribosome"]);

  // If we have lost ribosomes, decrement the ribosome assignments accordingly.
  if (net_lost_ribosomes > 0 && ribosome_assignments_.vals_.sum() > 0) {
    ArrayXd dist = ribosome_assignments_.normalized().vals_;
    dist /= dist.sum();
    dist *= net_lost_ribosomes;
    ribosome_assignments_.vals_ -= dist;
  }

  
  // Every tick, transcription factors are set to their defaults, then user code is executed to make changes.
  setDefaultTranscriptionFactors();
  for (auto dnaif : dna_ifs_) {
    dnaif->execute(cell);
  }

  // Check for funny business.
  for (int i = 0; i < transcription_factors_.vals_.size(); ++i)
    if (transcription_factors_.vals_[i] > 0)
      assert(pro_.molecule(i)->num_amino_acids_ > 0);
  
  // Transcription factor normalization:
  //   * If there are negative numbers, replace them with zeros.
  //   * If the transcription factors sum to greater than one, make them sum to one.
  //   * Otherwise leave them alone.

  MoleculeVals ntfs(transcription_factors_);
  ArrayXd& tf = ntfs.vals_;
  // User can easily specify DNA rules that cause these numbers to go negative.
  // Just clip them to zero.
  tf = tf.max(ArrayXd::Zero(tf.size()));  
  assert((ntfs.vals_ >= 0).all());

  // Many mRNA copies for one gene can be floating around in the cell, and one ribosome takes up the space of like 15-20 bases
  // on that strand.  For now we are just subsuming all this complexity into "over what distribution of proteins would you
  // like your ribosomes assigned".  Eventually it may be nice to simulate num strands of mRNA floating around, and do
  // a reaction on that between ribosomes and mRNA.
  // Anyway, what this means is: There is essentially no limit to the number of ribosomes that can be assigned to a particular gene.
  // Previously that logic was here, and now instead there is just this long note.

  ArrayXd& ra = ribosome_assignments_.vals_;
  
  // Normalize the distribution if greater than 1.
  if (tf.sum() > 1.0)
    tf /= tf.sum();

  // Assign ribosomes.
  assert(tf.size() == synthesis_reactions_.size());
  int num_free_ribosomes = num_rib - ra.floor().sum();
  ArrayXd target_ra = tf * num_rib;
  
  // cout << "============================================================" << endl;
  // cout << "Allocating " << num_free_ribosomes << " free ribosomes" << endl;
  // cout << "Ribosome half life: " << pro_.molecule("Ribosome")->half_life_hours_ << endl;
  // cout << "transcription_factors_: " << transcription_factors_.vals_.transpose() << endl;
  // cout << "tf: " << tf.transpose() << endl;
  // cout << "target_ra: " << target_ra.transpose() << endl;
  // cout << "ribosome_assignments_: " << ribosome_assignments_.vals_.transpose() << endl;

  // It's possible for the delta in ribosomes to be negative, in part because ribosomes can degrade,
  // and in part because of DNA programming that says to stop making a certain thing.
  ArrayXd delta = target_ra - ra;
  
  // cout << "delta: " << delta.transpose() << endl;
  // cout << "ra before adding delta: " << ra.transpose() << endl;
  if (delta.sum() > 0) {
    if (delta.sum() < num_free_ribosomes)
      ra += delta.floor();
    else
      ra += (delta * (num_free_ribosomes / delta.sum())).floor();
  }
  ra = ra.max(ArrayXd::Zero(ra.size())).floor();

  num_free_ribosomes = floor(num_rib - ra.sum());
  // cout << "num_rib: " << num_rib << endl;
  // cout << "cytosol_contents_[Ribosome]: " << cell->cytosol_contents_["Ribosome"] << endl;
  // cout << "ra: " << ra.transpose() << endl;
  // cout << "num_free_ribosomes: " << num_free_ribosomes << endl;
  PASSERT(num_free_ribosomes >= 0, fmt::format("num_free_ribosomes {}", num_free_ribosomes));

  assert(!ra.isNaN().any());
}

void DNA::computeGrad(Cell& cell)
{
  // In testing, sometimes we don't have ribosomes.
  // Eventually: get rid of this.
  if (!cell.cytosol_contents_.hasMolecule("Ribosome"))
    return;

  assignRibosomes(&cell);
  
  // Run protein synthesis reactions.  These are special reactions which are run by ribosomes,
  // consume only ATP and amino acids, and produce only proteins, ADP, and phosphate.
  cell.cytosol_contents_.probabilisticRound();
  assert(synthesis_reactions_.size() == ribosome_assignments_.size());
  for (size_t idx = 0; idx < synthesis_reactions_.size(); ++idx)
    if (synthesis_reactions_[idx])
      synthesis_reactions_[idx]->computeGrad(cell, ribosome_assignments_.vals_[idx]);
}

void DNA::apply(Cell& cell)
{
  MoleculeVals orig_cytosol_contents = cell.cytosol_contents_;

  for (auto srt : synthesis_reactions_)
    if (srt)
      srt->apply(cell);
  
  // cout << "[DNA::tick] cytosol contents before pround: " << endl << cell.cytosol_contents_.str("    ") << endl;
  cell.assertPositiveCytosolVals();
  
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

  ribosome_assignments_.probabilisticRound();
  assert((ribosome_assignments_.vals_ >= 0).all());
}

Prokaryotic::Prokaryotic()
{}


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
    step();
    double seconds = double(std::chrono::duration_cast<std::chrono::nanoseconds>(HRC::now() - start).count()) * 1e-9;
    cout << "seconds / step: " << seconds << endl;
    
    string response = comms_.waitForResponses();
    applyClientResponse(response);
  }
}

void Prokaryotic::applyClientResponse(const std::string& response)
{
  cells_[0]->dna_->dna_ifs_.clear();
  cout << "Cleared DNA programming." << endl;
  YAML::Node yaml = YAML::Load(response);
  for (const YAML::Node& dnaif : yaml["DNA"]) {
    cout << "Applying dnaif: " << dnaif << endl;
    cells_[0]->addDNAIf(dnaif);
  }
}

void Prokaryotic::addMoleculeType(MoleculeType::Ptr mt)
{
  molecule_types_.push_back(mt);
  molecule_map_[mt->name_] = mt;
}

void Prokaryotic::addMoleculeType(const YAML::Node& yaml)
{
  addMoleculeType(MoleculeType::Ptr(new MoleculeType(*this, yaml)));
}


std::vector<MoleculeType::ConstPtr> Prokaryotic::moleculeTypes() const
{
  return std::vector<MoleculeType::ConstPtr>(molecule_types_.begin(), molecule_types_.end());
}

std::vector<std::string> Prokaryotic::moleculeNames() const
{
  vector<string> names;
  for (auto mt : molecule_types_)
    names.push_back(mt->name_);
  return names;
}

void Prokaryotic::step()
{
  assert(biomes_.size() == 1);
  assert(cells_.size() == 1);
  for (auto biome : biomes_)
    biome->step(comms_);

  cout << cells_[0]->obs_.cytosolContentsHistoryAvg().str("  ") << endl;
  cout << "Denatured avg: " << endl;
  cout << cells_[0]->obs_.cytosolContentsDenaturedHistoryAvg().str("  ") << endl;
  cout << "Denatured instantaneous: " << endl;
  cout << cells_[0]->cytosol_contents_denatured_.str("  ") << endl;
  
  MessageWrapper msg;
  msg.addField("molecule_names", moleculeNames());
  msg.addField("division_hours", cells_[0]->obs_.averageDivisionHours());
  msg.addField("protein_io_flux", cells_[0]->obs_.protein_io_flux_.vals_);
  msg.addField("proteasome_action", cells_[0]->obs_.proteasome_action_.vals_);
  msg.addField("cytosol_contents_hist_avg", cells_[0]->obs_.cytosolContentsHistoryAvg().vals_);
  msg.addField("cytosol_contents_denatured_hist_avg", cells_[0]->obs_.cytosolContentsDenaturedHistoryAvg().vals_);
  msg.addField("cytosol_concentration_hist_avg", countsToConcentrations(cells_[0]->obs_.cytosolContentsHistoryAvg(), cells_[0]->um3_).vals_);
  comms_.broadcast(msg);
}

void Prokaryotic::tick()
{
  assert(biomes_.size() == 1);
  assert(cells_.size() == 1);
  for (auto biome : biomes_)
    biome->tick();
}

void Prokaryotic::applyConfig(const std::string& path)
{
  applyConfig(YAML::LoadFile(path));
}

void Prokaryotic::addReactionType(const YAML::Node& yaml)
{
  addReactionType(ReactionType::Ptr(new ReactionType(*this, yaml)));
}

void Prokaryotic::addBiome(const YAML::Node& yaml)
{
  addBiome(Biome::Ptr(new Biome(*this, yaml)));
}

void Prokaryotic::applyConfig(const YAML::Node& yaml)
{
  for (const YAML::Node& mol : yaml["MoleculeTable"])
    addMoleculeType(mol);
  for (const YAML::Node& rxn : yaml["ReactionTable"])
    addReactionType(rxn);

  // Eventually there will be multiple player cells here.
  cells_.push_back(Cell::Ptr(new Cell(*this, "cell")));
  
  for (const YAML::Node& biome : yaml["BiomeTable"])
    addBiome(biome);
}

MoleculeVals Prokaryotic::numAA() const
{
  MoleculeVals num_aa(*this);
  for (size_t i = 0; i < molecule_types_.size(); ++i)
    num_aa[i] = molecule_types_[i]->num_amino_acids_;
  return num_aa;
}
