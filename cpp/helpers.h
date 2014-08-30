#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <cstdio>
#include <complex>
#include <vector>

typedef std::complex<float> complexf;

bool read_input(complexf* data, const char* filename, long long n) {
  long long n2 = n * n;
  std::vector<float> buffer(n2);
  FILE* f = fopen(filename, "r");
  if (!f) {
    printf("Could not open file \"%s\".\n", filename);
  }
  size_t num_read;
  num_read = fread(&(buffer[0]), sizeof buffer[0], buffer.size(), f);
  if (num_read != buffer.size()) {
    printf("Error: read only %lu real values.\n", num_read);
    return false;
  }
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < n; ++col) {
      data[row * n + col].real(buffer[col * n + row]);
    }
  }
  num_read = fread(&(buffer[0]), sizeof buffer[0], buffer.size(), f);
  if (num_read != buffer.size()) {
    printf("Error: read only %lu imaginary values.\n", num_read);
    return false;
  }
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < n; ++col) {
      data[row * n + col].imag(buffer[col * n + row]);
    }
  }
  fclose(f);
  return true;
}

bool write_data(complexf* data, const char* filename, long long n) {
  long long n2 = n * n;
  std::vector<float> buffer(n2);
  FILE* f = fopen(filename, "w");
  if (!f) {
    printf("Could not open file \"%s\".\n", filename);
  }
  
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      buffer[col * n + row] = data[row * n + col].real();
    }
  }
  size_t num_written;
  num_written = fwrite(buffer.data(), sizeof buffer[0], buffer.size(), f);
  if (num_written != buffer.size()) {
    printf("Error: wrote only %lu real values.\n", num_written);
    return false;
  }

  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      buffer[col * n + row] = data[row * n + col].imag();
    }
  }
  num_written = fwrite(buffer.data(), sizeof buffer[0], buffer.size(), f);
  if (num_written != buffer.size()) {
    printf("Error: wrote only %lu imaginary values.\n", num_written);
    return false;
  }

  fclose(f);
  return true;
}

#endif
