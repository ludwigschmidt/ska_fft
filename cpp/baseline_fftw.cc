#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <vector>

#include <fftw3.h>
#include <unistd.h>

#include "helpers.h"
#include "peaks_helpers.h"

using namespace std;

struct Options {
  string input_filename;
  string output_filename;
  string peaks_filename;
  string regions_filename;
  long long n;
  int num_threads;
  int num_trials;
};

bool parse_options(Options* options, int argc, char** argv);


int main(int argc, char** argv) {
  Options opts;
  if (!parse_options(&opts, argc, argv)) {
    return 1;
  }
  const long long n = opts.n;

  vector<double> running_times;
  
  const long long n2 = n * n;

  if (fftwf_init_threads() == 0) {
    cout << "Error while initializing FFTW for threading." << endl;
    return 1;
  }

  fftwf_plan_with_nthreads(opts.num_threads);

  complexf* data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * n2));
  complexf* data2 = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * n2));

  // fftw setup
  auto start = chrono::high_resolution_clock::now();
  fftwf_plan plan = fftwf_plan_dft_2d(
      n, n,
      reinterpret_cast<fftwf_complex*>(data),
      reinterpret_cast<fftwf_complex*>(data2),
      FFTW_FORWARD, FFTW_MEASURE);
  auto end = chrono::high_resolution_clock::now();
  chrono::duration<double> fftw_setup_time = end - start;
  cout << "FFTW setup time: " << fftw_setup_time.count() << " s" << endl;

  start = chrono::high_resolution_clock::now();
  // read input data (frequency domain)
  if (!read_input(data, opts.input_filename.c_str(), n)) {
    return 1;
  }
  end = chrono::high_resolution_clock::now();
  chrono::duration<double> reading_input_time = end - start;
  cout << "Reading input time: " << reading_input_time.count() << " s" << endl;

  for (int trial = 0; trial < opts.num_trials; ++trial) {
    // execute fftw
    start = chrono::high_resolution_clock::now();
    fftwf_execute(plan);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> fftw_time = end - start;
    cout << "FFTW time: " << fftw_time.count() << " s" << endl;
    
    running_times.push_back(fftw_time.count());
  }

  sort(running_times.begin(), running_times.end());
  int num_outliers = opts.num_trials * 0.1;
  double mean = 0.0;
  for (size_t ii = 0; ii < running_times.size() - num_outliers; ++ii) {
    mean += running_times[ii];
  }
  mean /= (running_times.size() - num_outliers);
  cout << "Mean running time (excluding top " << num_outliers << " trials): "
       << mean << " s" << endl;

  if (opts.output_filename != "") {
    start = chrono::high_resolution_clock::now();
    if (!write_data(data2, opts.output_filename.c_str(), n)) {
      return 1;
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> writing_output_time = end - start;
    cout << "Writing output time: " << writing_output_time.count() << " s"
         << endl;
  }

  if (opts.regions_filename != "") {
    vector<region> regions;
    vector<peak> peaks;

    if (!read_regions(opts.regions_filename.c_str(), &regions)) {
      cout << "Could not read regions file.";
      return 1;
    }

    if (!find_peaks(data2, n, regions, &peaks)) {
      cout << "Could not find peaks.";
      return 1;
    }
    
    cout << "--------------------------------------------------------" << endl;
    cout << "Peaks:" << endl;
    for (size_t ii = 0; ii < peaks.size(); ++ii) {
      cout << peaks[ii].x << " " << peaks[ii].y << " " << peaks[ii].value
           << endl;
    }

    if (opts.peaks_filename != "") {
      if (!write_peaks_to_file(peaks, opts.peaks_filename)) {
        cout << "Error while writing peaks to file." << endl;
      }
    }
  }

  // clean-up
  fftwf_destroy_plan(plan);
  fftwf_free(data);
  fftwf_cleanup_threads();

  return 0;
}


bool parse_options(Options* options, int argc, char** argv) {
  options->input_filename = "";
  options->output_filename = "";
  options->peaks_filename = "";
  options->regions_filename = "";
  options->n = -1;
  options->num_threads = 1;
  options->num_trials = 1;

  int c;
  while ((c = getopt(argc, argv, "c:i:n:o:p:r:t:")) != -1) {
    if (c == 'c') {
      options->num_threads = stoi(string(optarg));
    } else if (c == 'i') {
      options->input_filename = string(optarg);
    } else if (c == 'n') {
      options->n = stoi(string(optarg));
    } else if (c == 'o') {
      options->output_filename = string(optarg);
    } else if (c == 'r') {
      options->regions_filename = string(optarg);
    } else if (c == 'p') {
      options->peaks_filename = string(optarg);
    } else if (c == 't') {
      options->num_trials = stoi(string(optarg));
    } else {
      printf("Could not parse command line options.\n");
      return false;
    }
  }

  if (options->input_filename == "") {
    printf("Input file option not set.\n");
    return false;
  }

  if (options->n == -1) {
    printf("Input size option not set.\n");
    return false;
  }

  return true;
}
