#include <comms.h>
#include <iostream>

using namespace std;

void MessageWrapper::addField(const std::string& name, const Eigen::ArrayXd& arr)
{
  append(name);
  append(dtype::ArrayXd);
  append((int)arr.size());
  uint8_t const *dptr = reinterpret_cast<uint8_t const *>(arr.data());
  for (int i = 0; i < arr.size() * 8; ++i)
    data_.push_back(dptr[i]);
}

void MessageWrapper::addField(const std::string& name, const Eigen::ArrayXXd& arr)
{
  append(name);
  append(dtype::ArrayXXd);
  append((int)arr.rows());
  append((int)arr.cols());
  uint8_t const *dptr = reinterpret_cast<uint8_t const *>(arr.data());
  for (int i = 0; i < arr.size() * 8; ++i)
    data_.push_back(dptr[i]);
}

void MessageWrapper::addField(const std::string& name, const std::vector<std::string>& strings)
{
  append(name);
  append(dtype::Strings);
  append(strings);
}

void MessageWrapper::append(int val)
{
  uint8_t* ptr = reinterpret_cast<uint8_t*>(&val);
  for (int i = 0; i < 4; ++i)
    data_.push_back(ptr[i]);
}

void MessageWrapper::append(const std::string& str)
{
  append((int)str.size());
  uint8_t const* ptr = reinterpret_cast<uint8_t const*>(str.data());
  for (size_t i = 0; i < str.size(); ++i)
    data_.push_back(ptr[i]);
}

void MessageWrapper::append(const std::vector<std::string>& strings)
{
  append((int)strings.size());
  for (string str : strings)
    append(str);
}

Comms::Comms() :
  sock_pub_(ctx_, zmq::socket_type::pub),
  sock_sub_(ctx_, zmq::socket_type::sub)
{
  sock_pub_.bind("tcp://127.0.0.1:53269");
  std::cout << "Connecting sock_pub_ to " << sock_pub_.get(zmq::sockopt::last_endpoint) << std::endl;
  
  sock_sub_.bind("tcp://127.0.0.1:53270");
  std::cout << "Connecting sock_sub_ to " << sock_sub_.get(zmq::sockopt::last_endpoint) << std::endl;
  sock_sub_.set(zmq::sockopt::subscribe, "");
}

void Comms::waitForResponses()
{
  zmq::message_t msg;
  std::optional<size_t> result = sock_sub_.recv(msg, zmq::recv_flags::none);
  assert(result.has_value());
  //assert(*result == 0);
  cout << "Got response: " << *result << " with msg: " << msg << endl;
}
