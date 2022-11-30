#include <sstream>
#include <vector>
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
  virtual std::string str(const std::string& prefix="") const
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

class Entity : public Printable
{
public:
  typedef std::shared_ptr<Entity> Ptr;
  typedef std::shared_ptr<const Entity> ConstPtr;
  
  int id_;
  std::string name_;
  std::string emoji_;

  Entity(const std::string& name, const std::string& emoji);
  virtual std::string _str() const;
  static size_t numEntities() { return num_entities_; }
  
private:
  static size_t num_entities_;
};

size_t Entity::num_entities_ = 0;

Entity::Entity(const std::string& name,
               const std::string& emoji) :
  id_(num_entities_),
  name_(name),
  emoji_(emoji)
{
  num_entities_ += 1;
  cout << "num_entities_: " << num_entities_ << endl;
  cout << "Created Entity with id " << id_ << endl;
}

std::string Entity::_str() const
{
  std::ostringstream oss;
  oss << "Entity \"" << name_ << "\"" << endl;
  oss << "id_: " << id_ << endl;
  oss << "emoji_: " << emoji_ << endl;
  return oss.str();
}

class Biome : public Printable
{
public:
  typedef std::shared_ptr<Biome> Ptr;
  typedef std::shared_ptr<const Biome> ConstPtr;

  std::vector<size_t> flux_;  // Aligned with Entity id_ values.
  std::string name_;
  
  Biome(const std::string& name);
  virtual std::string _str() const;
};

Biome::Biome(const std::string& name): name_(name)
{
  flux_ = vector<size_t>(Entity::numEntities(), 0);
}

std::string Biome::_str() const
{
  std::ostringstream oss;
  oss << "Biome \"" << name_ << "\"" << endl;
  return oss.str();
}
  
int main(int argc, char** argv)
{
  Entity e0("ATP", "bang");
  Entity e1("ADP", "wallet");
  //cout << e1.id() << endl;

  cout << e0.str("  ") << endl;

  Biome b0("Alkaline vents");
  //cout << b0.flux_ << endl;
  cout << b0.str() << endl;
  return 0;
}
