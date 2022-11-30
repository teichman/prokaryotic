#include <vector>
#include <iostream>
#include <memory>
#include <string>

using std::shared_ptr;
using std::cout, std::endl;
using std::vector;

class Entity
{
public:
  typedef std::shared_ptr<Entity> Ptr;
  typedef std::shared_ptr<const Entity> ConstPtr;
  
  int id_;
  std::string name_;

  Entity(const std::string& name);

  static size_t numEntities() { return num_entities_; }
  
private:
  static size_t num_entities_;
};

size_t Entity::num_entities_ = 0;

Entity::Entity(const std::string& name) : name_(name),
                                          id_(num_entities_)
{
  num_entities_ += 1;
  cout << "num_entities_: " << num_entities_ << endl;
  cout << "Created Entity with id " << id_ << endl;
}

class Biome
{
public:
  typedef std::shared_ptr<Biome> Ptr;
  typedef std::shared_ptr<const Biome> ConstPtr;

  std::vector<size_t> flux_;  // Aligned with Entity id_ values.
  std::string name_;
  
  Biome(const std::string& name);
};

Biome::Biome(const std::string& name): name_(name)
{
  flux_ = vector<size_t>(Entity::numEntities(), 0);
}
  
  
int main(int argc, char** argv)
{
  Entity e0("A");
  Entity e1("B");
  //cout << e1.id() << endl;

  Biome b0("Alkaline vents");
  //cout << b0.flux_ << endl;
  return 0;
}
