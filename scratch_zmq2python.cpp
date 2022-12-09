#include <chrono>
#include <thread>
#include <iostream>
#include <zmq_addon.hpp>

using namespace std;

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

    std::array<zmq::const_buffer, 3> send_msgs = {
        zmq::str_buffer("Hello"),
        zmq::str_buffer("over"),
        zmq::str_buffer("zmq")
    };

    while (true) {
      cout << "Published." << endl;
      if (!zmq::send_multipart(sock1, send_msgs))
        return 1;
      
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
