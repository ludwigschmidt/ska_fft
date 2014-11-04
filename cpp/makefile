slab_fft: slab_fft.cc helpers.h
	g++ -std=c++11 -Wall -Wextra -fopenmp -lfftw3f -lfftw3f_threads -O3 -o slab_fft slab_fft.cc

baseline_fftw: baseline_fftw.cc helpers.h
	g++ -std=c++11 -Wall -Wextra -lfftw3f -lfftw3f_threads -O3 -o baseline_fftw baseline_fftw.cc

check_output: check_output.cc helpers.h
	g++ -std=c++11 -Wall -Wextra -O3 -o check_output check_output.cc

pick_peaks: pick_peaks.cc helpers.h peaks_helpers.h
	g++ -std=c++11 -Wall -Wextra -O3 -o pick_peaks pick_peaks.cc
