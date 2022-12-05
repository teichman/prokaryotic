#include <prokaryotic.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
  Prokaryotic pro;
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ADP", "", 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "Phosphate", "", 1)));
  pro.addMoleculeType(MoleculeType::Ptr(new MoleculeType(pro, "ATP", "", 1)));
  YAML::Node yaml = YAML::Load("formula: ADP + Phosphate -> ATP\n"
                               "protein: ATP Synthase\n"
                               "kcat: 1e-3\n"
                               "KMs:\n"
                               "  ADP: 1e-1\n"
                               "  Phosphate: 1e-2");
  ReactionType rt(pro, yaml);
  cout << rt.str() << endl;
}

