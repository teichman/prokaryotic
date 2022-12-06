#include <iostream>
#include <yaml-cpp/yaml.h>
#include <yaml-cpp/node/emit.h>

using namespace std;

int main(int argc, char** argv)
{
  YAML::Node config = YAML::LoadFile(argv[1]);
  cout << "Entire file: " << endl;
  cout << YAML::Dump(config) << endl;

  const YAML::Node& dna = config["DNA"];
  for (const YAML::Node& cond : dna) {
    cout << "condition: " << cond["if"].as<string>() << endl;
    for (const YAML::Node& result : cond["then"]) {
      cout << "  result: " << result.as<string>() << endl;
    }
  }


  
  cout << "------------------------------------------------------------" << endl;
  cout << endl;
  cout << config["ReactionTable"] << endl;
  cout << "------------------------------------------------------------" << endl;
  
  YAML::Node rt = config["ReactionTable"];
  for (YAML::const_iterator it = rt.begin(); it != rt.end(); ++it) {
    const YAML::Node& rxn = *it;
    cout << "formula: " << rxn["formula"].as<string>() << endl;
    cout << "  kcat: " << rxn["kcat"].as<double>() << endl;
    cout << "  kms:" << endl;
    for (YAML::const_iterator kmit = rxn["KMs"].begin(); kmit != rxn["KMs"].end(); ++kmit) {
      const YAML::Node& key = kmit->first;
      const YAML::Node& value = kmit->second;
      cout << "    " << key << " " << value << endl;
    }
  }
  
  return 0;
}
