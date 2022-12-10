#include <comms.h>
#include <chrono>
#include <thread>
#include <iostream>

using namespace std;

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
