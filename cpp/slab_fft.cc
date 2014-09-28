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
#include <omp.h>
#include <unistd.h>

#include "helpers.h"

using namespace std;

typedef pair<int, int> location;

struct Options {
  string input_filename;
  string output_filename;
  long long n;
  int num_threads;
  int num_trials;
  bool use_diagonal_slabs;
};

struct Job {
  int thread_start;
  int thread_end;
  long long bin_start;
  long long bin_end;
};

bool parse_options(Options* options, int argc, char** argv);

void build_job_list(vector<Job>* jobs, const vector<vector<location> >& bins);


int main(int argc, char** argv) {
  Options opts;
  if (!parse_options(&opts, argc, argv)) {
    return 1;
  }
  
  if (fftwf_init_threads() == 0) {
    cout << "Error while initializing FFTW for threading." << endl;
    return 1;
  }

  const long long n = opts.n;

  vector<double> running_times;

  const long long n2 = n * n;
  const long long center_size = n / 16;
  const long long center_size_sq = center_size * center_size;
  const long long slab_size = n / 256;
  const long long total_slab_size = slab_size * n;

  // detection threshold
  const float threshold = 0.35;

  // number of frequencies per bin
  const long long bsize_slab = n / slab_size;
  const long long bsize_center = n / center_size;

  complexf* data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * n2));
  // output
  vector<complexf> output(n2);

  // OpenMP setup
  omp_set_num_threads(opts.num_threads);

  // fftw setup

  fftwf_plan_with_nthreads(opts.num_threads);

  // center part
  complexf* center_data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * center_size_sq));
  int center_n[2] = {static_cast<int>(center_size),
                     static_cast<int>(center_size)};
  int center_inembed[2] = {static_cast<int>(center_size),
                           static_cast<int>(n)};
  int center_onembed[2] = {static_cast<int>(center_size),
                           static_cast<int>(center_size)};
  fftwf_plan center_plan = fftwf_plan_many_dft(
      2,
      center_n,
      1,
      reinterpret_cast<fftwf_complex*>(data
          + (n / 2 - center_size / 2) * n + (n / 2 - center_size / 2)),
      center_inembed,
      1,
      0,
      reinterpret_cast<fftwf_complex*>(center_data),
      center_onembed,
      1,
      0,
      FFTW_FORWARD,
      FFTW_MEASURE);

  complexf* horiz_data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * total_slab_size));
  fftwf_plan horiz_plan = fftwf_plan_dft_2d(
      slab_size, n,
      reinterpret_cast<fftwf_complex*>(data + (n / 2 - slab_size / 2) * n),
      reinterpret_cast<fftwf_complex*>(horiz_data),
      FFTW_FORWARD, FFTW_MEASURE);

  complexf* vert_data = static_cast<complexf*>(
      fftwf_malloc(sizeof(fftwf_complex) * total_slab_size));
  int vert_n[2] = {static_cast<int>(n),
                   static_cast<int>(slab_size)};
  int vert_inembed[2] = {static_cast<int>(n),
                         static_cast<int>(n)};
  int vert_onembed[2] = {static_cast<int>(n),
                         static_cast<int>(slab_size)};
  fftwf_plan vert_plan = fftwf_plan_many_dft(
      2,
      vert_n,
      1,
      reinterpret_cast<fftwf_complex*>(data + (n / 2 - slab_size / 2)),
      vert_inembed,
      1,
      0,
      reinterpret_cast<fftwf_complex*>(vert_data),
      vert_onembed,
      1,
      0,
      FFTW_FORWARD,
      FFTW_MEASURE);

  // diagonal slabs
  complexf* diaghi_data = nullptr;
  complexf* diaglo_data = nullptr;
  if (opts.use_diagonal_slabs) {
    diaghi_data = static_cast<complexf*>(
        fftwf_malloc(sizeof(fftwf_complex) * total_slab_size));
    diaglo_data = static_cast<complexf*>(
        fftwf_malloc(sizeof(fftwf_complex) * total_slab_size));
  }


  // read input data (frequency domain)
  if (!read_input(data, opts.input_filename.c_str(), n)) {
    return 1;
  }
  
  // reserve thread-specific memory once
  vector<vector<location>> large_bins(opts.num_threads);
  vector<Job> jobs(opts.num_threads);


  for (int trial = 0; trial < opts.num_trials; ++trial) {
    auto overall_start = chrono::high_resolution_clock::now();

    // copy in data
    // diagonal slabs
    if (opts.use_diagonal_slabs) {
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
          diaglo_data[row * n + col] = data[orig_row * n + col];
        }
      }
    }
    //cout << diaglo_data[0] << " " << diaglo_data[1] << " " << diaglo_data[2]
    //    << endl;
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> copy_time = end - overall_start;

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
    if (opts.use_diagonal_slabs) {
      fftwf_execute_dft(horiz_plan,
                        reinterpret_cast<fftwf_complex*>(diaghi_data),
                        reinterpret_cast<fftwf_complex*>(diaghi_data));
      fftwf_execute_dft(horiz_plan,
                        reinterpret_cast<fftwf_complex*>(diaglo_data),
                        reinterpret_cast<fftwf_complex*>(diaglo_data));
    }
    //cout << diaghi_data[0] << " " << diaghi_data[1] << " " << diaghi_data[2]
    //    << endl;
    //cout << diaglo_data[0] << " " << diaglo_data[1] << " " << diaglo_data[2]
    //    << endl;
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> fftw_time = end - start;

    start = chrono::high_resolution_clock::now();
    
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      large_bins[id].clear();
      int row_start = id * center_size / opts.num_threads;
      int row_end = row_start + center_size / opts.num_threads;
      if (id == opts.num_threads) {
        row_end = center_size;
      }

      //printf("id = %d, row_start = %d, row_end = %d\n", id, row_start,
      //    row_end);

      for (int row = row_start; row < row_end; ++row) {
        for (int col = 0; col < center_size; ++col) {
          if (abs(center_data[row * center_size + col]) > threshold) {
            large_bins[id].push_back(make_pair(row, col));
          }
        }
      }
    }

    build_job_list(&jobs, large_bins);

    end = chrono::high_resolution_clock::now();
    chrono::duration<double> selection_time = end - start;

    /*cout << "Num large: " << large_bins.size() << endl;
    for (auto loc : large_bins) {
      cout << loc.first << " " << loc.second << endl;
    }*/

    int inner_iter = 0;
    start = chrono::high_resolution_clock::now();
    
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      Job& job = jobs[id];

      //printf("id = %d, job: (%d, %lld) to (%d, %lld)\n", id, job.thread_start,
      //    job.bin_start, job.thread_end, job.bin_end);

      for (int ithread = job.thread_start; ithread <= job.thread_end;
          ++ithread) {
        long long cur_bin_start = 0;
        if (ithread == job.thread_start) {
          cur_bin_start = job.bin_start;
        }
        long long cur_bin_end = static_cast<long long>(
            large_bins[ithread].size()) - 1;
        if (ithread == job.thread_end) {
          cur_bin_end = job.bin_end;
        }

        //printf("id = %d, thread = %d, bin_start = %lld, bin_end = %lld\n",
        //    id, ithread, cur_bin_start, cur_bin_end);
        //fflush(stdout);

        for (int ibin = cur_bin_start; ibin <= cur_bin_end; ++ibin) {
          const location& loc = large_bins[ithread][ibin];
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

              float average_value;
              if (opts.use_diagonal_slabs) {
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
                average_value = (center_value + horiz_value + vert_value
                    + diaghi_value + diaglo_value) / 5.0f;
              } else {
                average_value = (center_value + horiz_value + vert_value)
                    / 3.0f;
              }

              output[row * n + col] = average_value;

              /*cout << center_value << " " << horiz_value << " "
                   << vert_value << " " << diaghi_value << " "
                   << diaglo_value << " " << average_value << endl;*/
            }
          }
        }
      }
    }

    end = chrono::high_resolution_clock::now();
    chrono::duration<double> interpolation_time = end - start;
    chrono::duration<double> overall_time = end - overall_start;

    long long num_large_bins = 0;
    for (size_t ii = 0; ii < large_bins.size(); ++ii) {
      num_large_bins += large_bins[ii].size();
    }

    cout << "Data copy time:     " << copy_time.count() << " s   ("
         << copy_time.count() / overall_time.count() * 100.0 << "%)" << endl;
    cout << "FFTW time:          " << fftw_time.count() << " s   ("
         << fftw_time.count() / overall_time.count() * 100.0 << "%)" << endl;
    cout << "Selection time:     " << selection_time.count() << " s   (" 
         << selection_time.count() / overall_time.count() * 100.0
         << "%)" << endl;
    cout << "Interpolation time: " << interpolation_time.count() << " s   ("
         << interpolation_time.count() / overall_time.count() * 100.0 << "%)"
         << endl;
    cout << "Overall time:       " << overall_time.count() << " s" << endl;
    cout << "Number of selected large bins: " << num_large_bins << endl;
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
  if (opts.use_diagonal_slabs) {
    fftwf_free(diaghi_data);
    fftwf_free(diaglo_data);
  }
  fftwf_cleanup_threads();


  if (opts.output_filename != "") {
    if (!write_data(output.data(), opts.output_filename.c_str(), n)) {
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


bool parse_options(Options* options, int argc, char** argv) {
  options->input_filename = "";
  options->output_filename = "";
  options->n = -1;
  options->num_threads = 1;
  options->num_trials = 1;
  options->use_diagonal_slabs = false;

  int c;
  while ((c = getopt(argc, argv, "c:i:n:o:t:x")) != -1) {
    if (c == 'c') {
      options->num_threads = stoi(string(optarg));
    } else if (c == 'i') {
      options->input_filename = string(optarg);
    } else if (c == 'n') {
      options->n = stoi(string(optarg));
    } else if (c == 'o') {
      options->output_filename = string(optarg);
    } else if (c == 't') {
      options->num_trials = stoi(string(optarg));
    } else if (c == 'x') {
      options->use_diagonal_slabs = true;
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

void build_job_list(vector<Job>* jobs, const vector<vector<location> >& bins) {
  vector<Job>& jobref = *jobs;
  int num_threads = bins.size();
  long long num_large_pixels = 0;
  for (int ii = 0; ii < num_threads; ++ii) {
    num_large_pixels += bins[ii].size();
  }
  long long job_size = num_large_pixels / num_threads;
  
  int cur_thread = 0;
  long long cur_bin = 0;

  for (int ijob = 0; ijob < num_threads; ++ijob) {
    jobref[ijob].thread_start = cur_thread;
    jobref[ijob].bin_start = cur_bin;

    long long remaining_in_job = job_size;
    if (ijob == num_threads - 1) {
      jobref[ijob].thread_end = num_threads - 1;
      jobref[ijob].bin_end =
          static_cast<long long>(bins[num_threads - 1].size() - 1);
      break;
    }
    
    while (true) {
      long long remaining_in_thread = bins[cur_thread].size() - cur_bin;
      if (remaining_in_job >= remaining_in_thread) {
        cur_bin = 0;
        cur_thread += 1;
        remaining_in_job -= remaining_in_thread;
        if (remaining_in_job == 0) {
          jobref[ijob].thread_end = cur_thread;
          jobref[ijob].bin_end =
              static_cast<long long>(bins[cur_thread - 1].size()) - 1;
          break;
        }
      } else {
        jobref[ijob].thread_end = cur_thread;
        jobref[ijob].bin_end = cur_bin + remaining_in_job - 1;
        cur_bin += remaining_in_job;
        break;
      }
    }
  }
}
