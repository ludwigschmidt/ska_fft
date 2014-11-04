#ifndef __PEAKS_HELPERS_H__
#define __PEAKS_HELPERS_H__

#include <vector>

#include "helpers.h"

struct region {
  int x1, y1, x2, y2;
};

struct peak {
  long long x, y;
  double value;
};

bool read_regions(const char* filename, std::vector<region>* regions) {
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

bool find_peaks(const complexf* data,
                long long n,
                const std::vector<region>& regions,
                std::vector<peak>* peaks) {
  peaks->resize(regions.size());

  for (size_t ii = 0; ii < regions.size(); ++ii) {
    const region& cur = regions[ii];
    long long maxx = cur.x1;
    long long maxy = cur.y1;
    double maxv = abs(data[cur.y1 * n + cur.x1]);

    if (cur.y1 < 0 || cur.y2 >= n || cur.x1 < 0 || cur.x2 >= n) {
      return false;
    }

    for (long long yy = cur.y1; yy <= cur.y2; ++yy) {
      for (long long xx = cur.x1; xx <= cur.x2; ++xx) {
        double val = abs(data[yy * n + xx]);
        if (val > maxv) {
          maxx = xx;
          maxy = yy;
          maxv = val;
        }
      }
    }

    (*peaks)[ii].x = maxx;
    (*peaks)[ii].y = maxy;
    (*peaks)[ii].value = maxv;
  }

  return true;
}

bool write_peaks_to_file(const std::vector<peak>& peaks,
                         const std::string& filename) {
  FILE* f = fopen(filename.c_str(), "w");
  if (!f) {
    return false;
  }
  for (size_t ii = 0; ii < peaks.size(); ++ii) {
    fprintf(f, "%lld %lld %lf\n", peaks[ii].x, peaks[ii].y, peaks[ii].value);
  }
  fclose(f);
  return true;
}

#endif
