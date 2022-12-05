#include <molecule_type.pb.h>
#include <google/protobuf/util/json_util.h>
#include <iostream>

using namespace std;

std::string pb2json(const google::protobuf::Message& msg)
{
  google::protobuf::util::JsonPrintOptions options;    
  options.add_whitespace = true;
  options.always_print_primitive_fields = true;
  options.preserve_proto_field_names = true;
  
  std::string json_string;  
  MessageToJsonString(msg, &json_string, options);
  return json_string;
}


int main(int argc, char** argv)
{
  propro::MoleculeType mt_pb;
  cout << pb2json(mt_pb) << endl;
}
