#include <fstream>
#include <string>

struct stats_r {
  int A;
  int C;
  int T;
  int G;
  int nL;
};

extern stats_r stats(const char* memory, const size_t size);
extern stats_r stats2(const char* memory, const size_t size);
