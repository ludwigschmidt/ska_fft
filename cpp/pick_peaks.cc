#include <chrono>
#include <iostream>
#include <vector>

#include "helpers.h"

using namespace std;

struct region {
  int x1, y1, x2, y2;
};

bool read_regions(const char* filename, vector<region>* regions);

int main(int argc, char** argv) {
  // read program options
  if (argc != 4) {
    printf("Wrong number of arguments.\n");
    return 1;
  }
  const char* data_filename = argv[1];
  const int n = atoi(argv[2]);
  const char* region_filename = argv[3];

  const int n2 = n * n;

  vector<complexf> data(n2, 0);
  
  if (!read_input(data.data(), data_filename, n)) {
    return 1;
  }

  vector<region> regions;

  if (!read_regions(region_filename, &regions)) {
    return 1;
  }
  
  for (size_t ii = 0; ii < regions.size(); ++ii) {
    const region& cur = regions[ii];
    int maxx = cur.x1;
    int maxy = cur.y1;
    double maxv = abs(data[cur.y1 * n + cur.x1]);

    for (int yy = cur.y1; yy <= cur.y2; ++yy) {
      for (int xx = cur.x1; xx <= cur.x2; ++xx) {
        double val = abs(data[yy * n + xx]);
        if (val > maxv) {
          maxx = xx;
          maxy = yy;
          maxv = val;
        }
      }
    }

    printf("%d %d %lf\n", maxx, maxy, maxv);
  }


  return 0;
}


bool read_regions(const char* filename, vector<region>* regions) {
  FILE* f = fopen(filename, "r");
  if (!f) {
    printf("Could not open file \"%s\".\n", filename);
    return false;
  }
  regions->clear();

  region next_region;

  while (true) {
    int num_read = fscanf(f, "%d %d %d %d", &next_region.x1, &next_region.y1,
                                            &next_region.x2, &next_region.y2);
    if (num_read == EOF) {
      break;
    } else if (num_read < 4) {
      printf("Read incorrect number of values: %d\n", num_read);
      fclose(f);
      return false;
    }

    regions->push_back(next_region);
  }

  fclose(f);
  return true;
}
