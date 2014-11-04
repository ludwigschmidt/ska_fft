#include <chrono>
#include <iostream>
#include <vector>

#include <unistd.h>

#include "helpers.h"
#include "peaks_helpers.h"

using namespace std;

struct Options {
  string data_filename;
  long long n;
  string regions_filename;
};

bool parse_options(Options* options, int argc, char** argv);

int main(int argc, char** argv) {
  Options opts;
  if (!parse_options(&opts, argc, argv)) {
    return 1;
  }

  const long long n = opts.n;
  const long long n2 = n * n;

  vector<complexf> data(n2, 0);
  
  if (!read_input(data.data(), opts.data_filename.c_str(), n)) {
    return 1;
  }

  vector<region> regions;
  vector<peak> peaks;

  if (!read_regions(opts.regions_filename.c_str(), &regions)) {
    return 1;
  }

  if (!find_peaks(data, n, regions, &peaks)) {
    return 1;
  }

  for (size_t ii = 0; ii < peaks.size(); ++ii) {
    printf("%lld %lld %lf\n", peaks[ii].x, peaks[ii].y, peaks[ii].value);
  }

  return 0;
}

bool parse_options(Options* options, int argc, char** argv) {
  options->data_filename = "";
  options->regions_filename = "";
  options->n = -1;

  int c;
  while ((c = getopt(argc, argv, "i:n:r:")) != -1) {
    if (c == 'i') {
      options->data_filename = string(optarg);
    } else if (c == 'n') {
      options->n = stoi(string(optarg));
    } else if (c == 'r') {
      options->regions_filename = string(optarg);
    } else {
      printf("Could not parse command line options.\n");
      return false;
    }
  }

  if (options->data_filename == "") {
    printf("Data file option not set.\n");
    return false;
  }

  if (options->regions_filename == "") {
    printf("Regions file option not set.\n");
    return false;
  }

  if (options->n == -1) {
    printf("Input size option not set.\n");
    return false;
  }

  return true;
}
