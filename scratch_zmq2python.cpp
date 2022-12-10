#include <Eigen/Dense>
#include <chrono>
#include <thread>
#include <iostream>
#include <zmq_addon.hpp>

using namespace std;


class MessageWrapper
{
public:  
  MessageWrapper()
  {
    data_.push_back(13);  // magic number
  }
  
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



int main()
{
    zmq::context_t ctx;
    //zmq::socket_t sock1(ctx, zmq::socket_type::push);
    zmq::socket_t sock1(ctx, zmq::socket_type::pub);
    //zmq::socket_t sock2(ctx, zmq::socket_type::pull);
    
    sock1.bind("tcp://127.0.0.1:53269");
    const std::string last_endpoint =
        sock1.get(zmq::sockopt::last_endpoint);
    std::cout << "Connecting to "
              << last_endpoint << std::endl;
    
    //zmq::message_t msg;
    
    // std::array<zmq::const_buffer, 3> send_msgs = {
    //     zmq::str_buffer("Hello"),
    //     zmq::str_buffer("over"),
    //     zmq::str_buffer("zmq")
    // };

    while (true) {
      // if (!zmq::send_multipart(sock1, send_msgs))
      //   return 1;
      // vector<int> data;
      // data.push_back(13);
      // data.push_back(42);
      // zmq::message_t msg(data);
      //const char bytes[2] = {1, 2};
      //zmq::message_t msg(bytes, 2);

      // vector<double> data;
      // data.push_back(13);
      // data.push_back(42);
      // zmq::message_t msg(data);

      // Eigen::ArrayXd data = Eigen::ArrayXd::Zero(4);
      // data(0) = 13;
      // data(3) = 42;
      // MessageWrapper mw;
      // mw.append(data);


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
      
      zmq::message_t msg(mw);
      
      // double payload = 3.14159;
      // zmq::message_t msg(&payload, 8);
      
      cout << "Publishing " << msg << " of size " << msg.size() << endl;
      sock1.send(msg, zmq::send_flags::none);
      
      
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    // std::vector<zmq::message_t> recv_msgs;
    // const auto ret = zmq::recv_multipart(
    //     sock2, std::back_inserter(recv_msgs));
    // if (!ret)
    //     return 1;
    // std::cout << "Got " << *ret
    //           << " messages" << std::endl;
    
    return 0;
}
