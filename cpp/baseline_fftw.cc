#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <vector>

#include <fftw3.h>

#include "helpers.h"

using namespace std;


int main(int argc, char** argv) {
  // read program options
  if (argc != 5) {
    printf("Wrong number of arguments.\n");
    return 1;
  }
  const char* filename = argv[1];
  const long long n = atoi(argv[2]);
  const char* output_filename = argv[3];
  const int num_trials = atoi(argv[4]);

  vector<double> running_times;
  
  const long long n2 = n * n;

  complexf* data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * n2));

  // fftw setup
  auto start = chrono::high_resolution_clock::now();
  fftwf_plan plan = fftwf_plan_dft_2d(
      n, n,
      reinterpret_cast<fftwf_complex*>(data),
      reinterpret_cast<fftwf_complex*>(data),
      FFTW_FORWARD, FFTW_MEASURE);
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<double> fftw_setup_time = end - start;
  cout << "FFTW setup time: " << fftw_setup_time.count() << " s" << endl;

  start = chrono::high_resolution_clock::now();
  // read input data (frequency domain)
  if (!read_input(data, filename, n)) {
    return 1;
  }
  end = chrono::high_resolution_clock::now();
  chrono::duration<double> reading_input_time = end - start;
  cout << "Reading input time: " << reading_input_time.count() << " s" << endl;

  for (int trial = 0; trial < num_trials; ++trial) {
    // execute fftw
    start = chrono::high_resolution_clock::now();
    fftwf_execute(plan);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> fftw_time = end - start;
    cout << "FFTW time: " << fftw_time.count() << " s" << endl;
    
    running_times.push_back(fftw_time.count());
  }

  sort(running_times.begin(), running_times.end());
  int num_outliers = num_trials * 0.1;
  double mean = 0.0;
  for (size_t ii = 0; ii < running_times.size() - num_outliers; ++ii) {
    mean += running_times[ii];
  }
  mean /= (running_times.size() - num_outliers);
  cout << "Mean running time (excluding top " << num_outliers << " trials): "
       << mean << " s" << endl;

  if (strcmp(output_filename, "none") != 0) {
    start = chrono::high_resolution_clock::now();
    if (!write_data(data, output_filename, n)) {
      return 1;
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> writing_output_time = end - start;
    cout << "Writing output time: " << writing_output_time.count() << " s"
         << endl;
  }

  // clean-up
  fftwf_destroy_plan(plan);
  fftwf_free(data);

  return 0;
}
