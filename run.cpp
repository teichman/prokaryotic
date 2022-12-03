#include <prokaryotic.h>

int main(int argc, char** argv)
{
  Prokaryotic prokaryotic;
  prokaryotic.initializeHardcoded();
  //prokaryotic.runTests();
  prokaryotic.run();
  return 0;
}
