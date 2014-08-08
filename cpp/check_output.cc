#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <vector>

#include <fftw3.h>

#include "helpers.h"

using namespace std;

typedef pair<int, int> location;


int main(int argc, char** argv) {
  // read program options
  if (argc !=  4) {
    printf("Wrong number of arguments.\n");
    return 1;
  }
  const char* test_filename = argv[1];
  const char* reference_filename = argv[2];
  const int n = atoi(argv[3]);
  const int n2 = n * n;

  // read test and reference data
  vector<complexf> test_data(n2);
  if (!read_input(test_data.data(), test_filename, n)) {
    return 1;
  }
  vector<complexf> reference_data(n2);
  if (!read_input(reference_data.data(), reference_filename, n)) {
    return 1;
  }

  // compute l1 norm of difference
  double sum = 0.0;
  for (int ii = 0; ii < n2; ++ii) {
    sum += abs(test_data[ii] - reference_data[ii]);
  }

  printf("l1 norm of differnce: %.6e\n", sum);

  return 0;
}
