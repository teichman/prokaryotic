#include <iostream>
#include <yaml-cpp/yaml.h>
#include <yaml-cpp/node/emit.h>


using namespace std;

int main(int argc, char** argv)
{
  // For defining yaml directly in the cpp code, I guess this will do.
  string mlsl =
    "This is a \n"
    "multiline\n"
    "string literal.";
  cout << mlsl << endl;

  // This lets one remove the \ns, but ofc you can't indent it properly or you get extra spaces in the raw string.
  string mlsl2 = R"(This is the first line.
This is the second line.
This is the third line.)";

  cout << mlsl2 << endl;
  
  cout << "Hi" << endl;

  YAML::Node config = YAML::LoadFile(argv[1]);
  cout << "Entire file: " << endl;
  cout << YAML::Dump(config) << endl;

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
  
  // YAML::Node config = YAML::Load("username: aoeu\nmouse: snth");
  
  // //YAML::Node config = YAML::LoadFile("config.yaml");
  // if (config["username"])
  //   cout << "username: " << config["username"].as<std::string>() << endl;
  // if (config["mouse"])
  //   cout << "mouse: " << config["mouse"].as<std::string>() << endl;

  // // const std::string username = config["username"].as<std::string>();
  // // const std::string password = config["password"].as<std::string>();

  // // std::ofstream fout("config.yaml");
  // // fout << config;
  
  return 0;
}
