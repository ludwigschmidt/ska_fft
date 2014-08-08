#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <iostream>
#include <utility>
#include <vector>

#include <fftw3.h>

#include "helpers.h"

using namespace std;

typedef pair<int, int> location;


int main(int argc, char** argv) {
  // read program options
  if (argc !=  5) {
    printf("Wrong number of arguments.\n");
    return 1;
  }
  const char* filename = argv[1];
  const int n = atoi(argv[2]);
  const char* output_filename = argv[3];
  const int num_trials = atoi(argv[4]);

  vector<double> running_times;

  const int n2 = n * n;
  const int center_size = n / 16;
  const int center_size_sq = center_size * center_size;
  const int slab_size = n / 256;
  const int total_slab_size = slab_size * n;

  // detection threshold
  const float threshold = 0.35;

  // number of frequencies per bin
  const int bsize_slab = n / slab_size;
  const int bsize_center = n / center_size;

  // read input data (frequency domain)
  vector<complexf> data(n2);
  if (!read_input(data.data(), filename, n)) {
    return 1;
  }

  // output
  vector<complexf> output(n2);


  // fftw setup

  // center part
  complexf* center_data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * center_size_sq));
  fftwf_plan center_plan = fftwf_plan_dft_2d(
      center_size, center_size,
      reinterpret_cast<fftwf_complex*>(center_data),
      reinterpret_cast<fftwf_complex*>(center_data),
      FFTW_FORWARD, FFTW_MEASURE);
  // horizontal and vertical slabs
  complexf* horiz_data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * total_slab_size));
  complexf* vert_data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * total_slab_size));
  fftwf_plan horiz_plan = fftwf_plan_dft_2d(
      slab_size, n,
      reinterpret_cast<fftwf_complex*>(horiz_data),
      reinterpret_cast<fftwf_complex*>(horiz_data),
      FFTW_FORWARD, FFTW_MEASURE);
  fftwf_plan vert_plan = fftwf_plan_dft_2d(
      n, slab_size,
      reinterpret_cast<fftwf_complex*>(vert_data),
      reinterpret_cast<fftwf_complex*>(vert_data),
      FFTW_FORWARD, FFTW_MEASURE);
  // diagonal slabs
  complexf* diaghi_data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * total_slab_size));
  complexf* diaglo_data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * total_slab_size));


  for (int trial = 0; trial < num_trials; ++trial) {
    auto overall_start = chrono::high_resolution_clock::now();

    // copy in data
    // center part
    for (int row = 0; row < center_size; ++row) {
      int orig_row = n / 2 - center_size / 2 + row;

      memcpy(&(center_data[row * center_size]),
             &(data[orig_row * n + n / 2 - center_size / 2]),
             center_size * sizeof(complexf));
    }
    // horizontal and vertical slabs
    for (int row = 0; row < slab_size; ++row) {
      int orig_row = n / 2 - slab_size / 2 + row;
      memcpy(&(horiz_data[row * n]),
             &(data[orig_row * n]),
             n * sizeof(complexf));
    }
    for (int row = 0; row < n; ++row) {
      memcpy(&(vert_data[row * slab_size]),
             &(data[row * n + n / 2 - slab_size / 2]),
             slab_size * sizeof(complexf));
    }
    // diagonal slabs
    for (int col = 0; col < n; ++col) {
      for (int row = 0; row < slab_size; ++row) {
        int orig_row = col - slab_size / 2 + row;
        if (orig_row < 0) {
          orig_row += n;
        }
        if (orig_row >= n) {
          orig_row -= n;
        }
        diaghi_data[row * n + col] = data[orig_row * n + col];
      }
    }
    for (int col = 0; col < n; ++col) {
      for (int row = 0; row < slab_size; ++row) {
        int orig_row = -col + slab_size / 2 - row;
        if (orig_row < 0) {
          orig_row += n;
        }
        if (orig_row < 0) {
          orig_row += n;
        }
        if (orig_row >= n) {
          orig_row -= n;
        }
        //if (col == 0 || col == 1) {
        //  cout << orig_row << endl;
        //}
        diaglo_data[row * n + col] = data[orig_row * n + col];
      }
    }
    //cout << diaglo_data[0] << " " << diaglo_data[1] << " " << diaglo_data[2]
    //    << endl;
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> copy_time = end - overall_start;
    cout << "Data copy time: " << copy_time.count() << " s" << endl;

    // execute fftw
    auto start = chrono::high_resolution_clock::now();
    // center part
    fftwf_execute(center_plan);
    //cout << center_data[0] << " " << center_data[1] << " " << center_data[2]
    //    << endl;
    // horizontal and vertical slabs
    fftwf_execute(horiz_plan);
    fftwf_execute(vert_plan);
    //cout << horiz_data[0] << " " << horiz_data[1] << " " << horiz_data[2]
    //     << endl;
    //cout << vert_data[0] << " " << vert_data[1] << " " << vert_data[2]
    //     << endl;
    // diagonal slabs
    fftwf_execute_dft(horiz_plan,
                      reinterpret_cast<fftwf_complex*>(diaghi_data),
                      reinterpret_cast<fftwf_complex*>(diaghi_data));
    fftwf_execute_dft(horiz_plan,
                      reinterpret_cast<fftwf_complex*>(diaglo_data),
                      reinterpret_cast<fftwf_complex*>(diaglo_data));
    //cout << diaghi_data[0] << " " << diaghi_data[1] << " " << diaghi_data[2]
    //    << endl;
    //cout << diaglo_data[0] << " " << diaglo_data[1] << " " << diaglo_data[2]
    //    << endl;
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> fftw_time = end - start;
    cout << "FFTW time: " << fftw_time.count() << " s" << endl;

    start = chrono::high_resolution_clock::now();
    vector<location> large_bins;
    for (int row = 0; row < center_size; ++row) {
      for (int col = 0; col < center_size; ++col) {
        if (abs(center_data[row * center_size + col]) > threshold) {
          large_bins.push_back(make_pair(row, col));
        }
      }
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> selection_time = end - start;
    cout << "Peak selection time: " << selection_time.count() << " s" << endl;

    /*cout << "Num large: " << large_bins.size() << endl;
    for (auto loc : large_bins) {
      cout << loc.first << " " << loc.second << endl;
    }*/

    int inner_iter = 0;
    start = chrono::high_resolution_clock::now();

    for (auto loc : large_bins) {
      int col_bin = loc.second;
      int row_bin = loc.first;

      int row_start = row_bin * bsize_center - bsize_center / 2 - 1;
      int row_end = row_bin * bsize_center + bsize_center / 2 - 1;
      int col_start = col_bin * bsize_center - bsize_center / 2 - 1;
      int col_end = col_bin * bsize_center + bsize_center / 2 - 1;

      for (int row = row_start; row <= row_end; ++row) {
        for (int col = col_start; col <= col_end; ++col) {
          inner_iter += 1;

          /*if (col != 503 || row != 1511) {
            continue;
          }
          cout << row << " " << col << endl;*/
          
          if (row >= n) {
            row -= n;
          }
          if (row < 0) {
            row += n;
          }

          if (col >= n) {
            col -= n;
          }
          if (col < 0) {
            col += n;
          }

          int center_bin_row = (row + 1) / bsize_center;
          if (row + 1 - center_bin_row * bsize_center >= bsize_center / 2) {
            center_bin_row += 1;
            if (center_bin_row >= n) {
              center_bin_row -= n;
            }
          }
          int center_bin_col = (col + 1) / bsize_center;
          if (col + 1 - center_bin_col * bsize_center >= bsize_center / 2) {
            center_bin_col += 1;
            if (center_bin_col >= n) {
              center_bin_col -= n;
            }
          }
          
          //cout << center_bin_row << " " << center_bin_col << endl;

          float center_value = abs(center_data[
              center_bin_row * center_size + center_bin_col]);

          int horiz_bin_row = (row + 1) / bsize_slab;
          if (row + 1 - horiz_bin_row * bsize_slab >= bsize_slab / 2) {
            horiz_bin_row += 1;
            if (horiz_bin_row >= n) {
              horiz_bin_row -= n;
            }
          }

          //cout << "horiz_bin_row: " << horiz_bin_row << endl;
          
          float horiz_value = abs(horiz_data[horiz_bin_row * n + col]);

          int vert_bin_col = (col + 1) / bsize_slab;
          if (col + 1 - vert_bin_col * bsize_slab >= bsize_slab / 2) {
            vert_bin_col += 1;
            if (vert_bin_col >= n) {
              vert_bin_col -= n;
            }
          }

          //cout << "vert_bin_col: " << vert_bin_col << endl;
          
          float vert_value = abs(vert_data[row * slab_size + vert_bin_col]);

          int diaghi_bin_row = horiz_bin_row;
          int diaghi_bin_col = row + col;
          if (diaghi_bin_col >= n) {
            diaghi_bin_col -= n;
          }

          //cout << "diaghi_bin_row: " << diaghi_bin_row
          //     << "   diaghi_bin_col: " << diaghi_bin_col << endl;

          float diaghi_value =
              abs(diaghi_data[diaghi_bin_row * n + diaghi_bin_col]);

          int diaglo_bin_row = vert_bin_col;
          int diaglo_bin_col = -row + col;
          if (diaglo_bin_col < 0) {
            diaglo_bin_col += n;
          }

          //cout << "diaglo_bin_row: " << diaglo_bin_row
          //     << "   diaglo_bin_col: " << diaglo_bin_col << endl;

          float diaglo_value =
              abs(diaglo_data[diaglo_bin_row * n + diaglo_bin_col]);

          float average_value = (center_value + horiz_value + vert_value
              + diaghi_value + diaglo_value) / 5.0f;
          output[row * n + col] = average_value;

          /*cout << center_value << " " << horiz_value << " "
               << vert_value << " " << diaghi_value << " "
               << diaglo_value << " " << average_value << endl;*/
        }
      }
    }

    end = chrono::high_resolution_clock::now();
    chrono::duration<double> interpolation_time = end - start;
    cout << "Interpolation time: " << interpolation_time.count() << " s"
         << endl;

    chrono::duration<double> overall_time = end - overall_start;
    cout << "Overall time: " << overall_time.count() << " s" << endl;

    cout << "Data copy time: "
         << copy_time.count() / overall_time.count() * 100.0 << "%" << endl;
    cout << "FFTW time: "
         << fftw_time.count() / overall_time.count() * 100.0 << "%" << endl;
    cout << "Selection time: "
         << selection_time.count() / overall_time.count() * 100.0 << "%"
         << endl;
    cout << "Interpolation time: "
         << interpolation_time.count() / overall_time.count() * 100.0 << "%"
         << endl;

    cout << "Number of selected large bins: " << large_bins.size() << endl;
    cout << "Total inner interpolation iterations: " << inner_iter << endl;
    cout << "--------------------------------------------------------" << endl;

    running_times.push_back(overall_time.count());
  }

  // clean-up
  fftwf_destroy_plan(center_plan);
  fftwf_free(center_data);
  fftwf_destroy_plan(horiz_plan);
  fftwf_free(horiz_data);
  fftwf_destroy_plan(vert_plan);
  fftwf_free(vert_data);
  fftwf_free(diaghi_data);
  fftwf_free(diaglo_data);


  if (strcmp(output_filename, "none") != 0) {
    if (!write_data(output.data(), output_filename, n)) {
      return 1;
    }
  }

  sort(running_times.begin(), running_times.end());
  int outliers = 0.1 * running_times.size();
  double mean = 0.0;
  for (size_t ii = 0; ii < running_times.size() - outliers; ++ii) {
    mean += running_times[ii];
  }
  mean /= (running_times.size() - outliers);
  cout << "Mean running time (excluding top " << outliers << " trials): "
       << mean << " s" << endl;

  /*
  // Compute statistics
  complex<float> sum = 0;
  for (int ii = 0; ii < n2; ++ii) {
    sum += data[ii];
  }
  cout << "Sum: " << sum << endl;

  complex<float> sumsq = 0;
  for (int ii = 0; ii < n2; ++ii) {
    sumsq += data[ii] * data[ii];
  }
  cout << "Sum of squares: " << sumsq << endl;

  float sumabs = 0;
  for (int ii = 0; ii < n2; ++ii) {
    sumabs += abs(data[ii]);
  }
  cout << "Sum of absolute values: " << sumabs << endl;
  */

  return 0;
}
