#include <comms.h>

int main(int argc, char** argv)
{
  Comms comms;
  while (true)
    comms.waitForResponses();
}

