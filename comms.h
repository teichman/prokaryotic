#include <vector>
#include <string>
#include <cstdint>
#include <Eigen/Dense>
#include <zmq.hpp>

// https://github.com/zeromq/cppzmq

class MessageWrapper
{
public:  
  MessageWrapper() {  data_.push_back(13); }  // magic number
  operator zmq::message_t () const { return zmq::message_t(data_); }
  size_t size() const { return data_.size(); }
  
  void addField(const std::string& name, const Eigen::ArrayXd& arr);
  void addField(const std::string& name, const Eigen::ArrayXXd& arr);
  void addField(const std::string& name, const std::vector<std::string>& strings);

private:
  enum dtype {
    String=10,
    Strings=11,
    ArrayXd=12,
    ArrayXXd=13
  };
  
  std::vector<uint8_t> data_;

  void append(dtype val) { data_.push_back(val); }
  void append(int val);
  void append(const std::string& str);
  void append(const std::vector<std::string>& strings);
};

class Comms
{
public:
  Comms();
  ~Comms() { sock_pub_.close(); sock_sub_.close(); }
  void broadcast(const MessageWrapper& msg) { sock_pub_.send(msg, zmq::send_flags::none); }
  // blocking
  void receive();
  void waitForResponses();

private:
  zmq::context_t ctx_;
  zmq::socket_t sock_pub_;
  zmq::socket_t sock_sub_;
};

