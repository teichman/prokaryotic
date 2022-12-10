#include <Eigen/Dense>
#include <chrono>
#include <thread>
#include <iostream>
#include <zmq.hpp>

using namespace std;

class MessageWrapper
{
public:  
  MessageWrapper()
  {
    data_.push_back(13);  // magic number
  }

  size_t size() const { return data_.size(); }
  
  void addField(const std::string& name, const Eigen::ArrayXd& arr)
  {
    append(name);
    append(dtype::ArrayXd);
    append((int)arr.size());
    uint8_t const *dptr = reinterpret_cast<uint8_t const *>(arr.data());
    for (int i = 0; i < arr.size() * 8; ++i)
      data_.push_back(dptr[i]);
  }

  void addField(const std::string& name, const Eigen::ArrayXXd& arr)
  {
    append(name);
    append(dtype::ArrayXXd);
    append((int)arr.rows());
    append((int)arr.cols());
    uint8_t const *dptr = reinterpret_cast<uint8_t const *>(arr.data());
    for (int i = 0; i < arr.size() * 8; ++i)
      data_.push_back(dptr[i]);
  }
  
  operator zmq::message_t () const
  {
    return zmq::message_t(data_);
  }

private:
  enum dtype {
    String=9,
    ArrayXd=10,
    ArrayXXd=11
  };
  
  vector<uint8_t> data_;

  void append(dtype val) { data_.push_back(val); }
  
  void append(int val)
  {
    uint8_t* ptr = reinterpret_cast<uint8_t*>(&val);
    for (int i = 0; i < 4; ++i)
      data_.push_back(ptr[i]);
  }
  
  void append(const std::string& str)
  {
    append((int)str.size());
    uint8_t const* ptr = reinterpret_cast<uint8_t const*>(str.data());
    for (size_t i = 0; i < str.size(); ++i)
      data_.push_back(ptr[i]);
  }
    
};

class Comms
{
public:
  zmq::context_t ctx_;
  zmq::socket_t sock_;
  
  Comms() :
    sock_(ctx_, zmq::socket_type::pub)
  {
    sock_.bind("tcp://127.0.0.1:53269");
    string last_endpoint = sock_.get(zmq::sockopt::last_endpoint);
    std::cout << "Connecting to " << last_endpoint << std::endl;
  }

  void broadcast(const MessageWrapper& msg)
  {
    sock_.send(msg, zmq::send_flags::none);
  }
};


int main()
{
  Comms comms;
            
  while (true) {
    MessageWrapper mw;
    Eigen::ArrayXXd mat = Eigen::ArrayXXd::Zero(2, 2);
    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(1, 0) = 3;
    mat(1, 1) = 4;
    mw.addField("matrix_something", mat);
    Eigen::ArrayXd vec = Eigen::ArrayXd::Zero(4);
    vec(0) = 13;
    vec(3) = 42;
    mw.addField("vector_something", vec);
      
    //zmq::message_t msg(mw);
    cout << "Publishing msg of size " << mw.size() << endl;
    comms.broadcast(mw);
      
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
    
  return 0;
}
