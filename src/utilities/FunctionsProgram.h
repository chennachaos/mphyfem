
#ifndef incl_FunctionsProgram_h
#define incl_FunctionsProgram_h

#include <fstream>
#include <iostream>


inline  void  prgError(int id, char* fname, char* msg)
{
  std::cout << " ERROR (" << id << ") in " << fname << ": " << msg << "\n\n";

  exit(0);
}


inline  void  prgWarning(int id, char* fname, char* msg)
{
  std::cout << " WARNING (" << id << ") in " << fname << ": " << msg << "\n\n";

  return;
}

void       prgReadTimeFunctions          (std::ifstream &);



#endif

