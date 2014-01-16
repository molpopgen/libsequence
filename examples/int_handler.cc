#include "int_handler.hpp"
#include <iostream>
void cntrl_c_handler(int sig)
{
  char answer;
  std::cerr << "Interrupt signal received.\n";
  std::cerr << "Do you really want to quit [y or n]?\n";
  std::cin >> answer;
  switch (answer)
    {
    case 'Y':
      exit(0);
      break;
    case 'y':
      exit(0);
      break;
    default:
      signal(SIGINT,cntrl_c_handler);
      std::cerr << "continuing"<<endl;
      break;
    }
}
