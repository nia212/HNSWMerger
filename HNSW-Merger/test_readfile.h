
// #include "hnswlib.h"
#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <cstdint> // For int64_t
#include <sstream> // For stringstream

using namespace std;


float *fvecs_read(const char *fname, size_t *d_out, size_t *n_out)
{
  FILE *f = fopen(fname, "r");
  printf("fvecs_read: %s\n", fname);
  if (!f)
  {
    fprintf(stderr, "could not open %s\n", fname);
    perror("");
    abort();
  }
  int d = *d_out;
  size_t n = *n_out;
  float *x = new float[n * (d + 1)];
  size_t nr = fread(x, sizeof(float), n * (d + 1), f);
  assert(nr == n * (d + 1) || !"could not read whole file");
  for (size_t i = 0; i < n; i++)
    memmove(x + i * d, x + 1 + i * (d + 1), d * sizeof(*x));

  fclose(f);
  return x;
}
int *ivecs_read(const char *fname, size_t *d_out, size_t *n_out)
{
  return (int *)fvecs_read(fname, d_out, n_out);
}
float *bvecs_read(const char *input_file, int num_vectors_to_read, size_t *d_out, size_t *n_out)
{
  float *result = nullptr;
  unsigned char *tmp = nullptr;
  FILE *file = fopen(input_file, "r"); // Open the file in binary read mode
  if (!file)
  {
    std::cerr << "Error: Unable to open file " << input_file << std::endl;
    return result;
  }
  int32_t d;
  int32_t n;
  fread(&d, sizeof(int32_t), 1, file);
  *d_out = d;
  fseek(file, 0, SEEK_END);
  long file_size = ftell(file);
  long total_n = file_size / (d + 4);
  num_vectors_to_read = std::min(num_vectors_to_read, (int)total_n);
  *n_out = total_n;
  n = num_vectors_to_read;
  int64_t tmp_elements = (int64_t)n * (d + 4);
  int64_t num_elements = (int64_t)num_vectors_to_read * d;
  result = new float[num_elements];
  tmp = new unsigned char[tmp_elements];

  fseek(file, 0, SEEK_SET);
  fread(tmp, tmp_elements, 1, file);
  for (int64_t i = 0; i < num_vectors_to_read; ++i)
  {
    for (int j = 0; j < d; ++j)
    {
      result[i * d + j] = static_cast<float>(tmp[i * (d + 4) + 4 + j]);
    }
  }
  fclose(file);
  delete[] tmp;
  return result;
}

float* load_and_convert_to_float(const char* fname, int num_vectors_to_read, size_t* d_out, size_t* n_out, size_t batch_size = 10000000) {
  FILE* f = fopen(fname, "rb");
  if (!f) {
      fprintf(stderr, "Could not open %s\n", fname);
      perror("");
      abort();
  }

  int d;
  fread(&d, sizeof(int), 1, f);
  assert((d > 0 && d < 1000000) || !"Unreasonable dimension");
  fseek(f, 0, SEEK_SET);
  struct stat st;
  fstat(fileno(f), &st);
  size_t sz = st.st_size;

  size_t per_vector_size = sizeof(int) + d * sizeof(uint8_t);
  assert(sz % per_vector_size == 0 || !"Weird file size");

  size_t n = sz / per_vector_size; 
  n = num_vectors_to_read > 0 ? std::min(num_vectors_to_read, (int)n) : n; 
  *d_out = d;
  *n_out = n;

  std::cout << "d: " << d << "  n:" << n << std::endl;
  float* result = new float[n * d];

  size_t vectors_left = n;
  size_t offset = 0;

  while (vectors_left > 0) {
      size_t current_batch_size = std::min(batch_size, vectors_left);

      std::vector<uint8_t> temp(per_vector_size * current_batch_size);
      size_t nr = fread(temp.data(), sizeof(uint8_t), temp.size(), f);
      std::cout << "vector size:"<<  current_batch_size <<" nr:" << nr << "  temp size:" << temp.size() << std::endl;
      assert(nr == temp.size() || !"Could not read batch");

      for (size_t i = 0; i < current_batch_size; i++) {
          const uint8_t* data_ptr = temp.data() + i * per_vector_size + sizeof(int);

          for (size_t j = 0; j < d; j++) {
              result[offset + i * d + j] = static_cast<float>(data_ptr[j]);
          }
      }

      offset += current_batch_size * d;
      vectors_left -= current_batch_size;
  }

  fclose(f);
  return result;
}

float* load_and_convert_to_float_range(const char* fname, size_t start_index, size_t end_index, size_t* d_out, size_t* n_out, size_t batch_size = 10000000) {
  FILE* f = fopen(fname, "rb");
  if (!f) {
      fprintf(stderr, "Could not open %s\n", fname);
      perror("");
      abort();
  }

  int d;
  fread(&d, sizeof(int), 1, f);
  assert((d > 0 && d < 1000000) || !"Unreasonable dimension");
  size_t per_vector_size = sizeof(int) + d * sizeof(uint8_t);
  fseek(f, 0, SEEK_SET);

  struct stat st;
  fstat(fileno(f), &st);
  size_t sz = st.st_size;
  size_t total_vectors = sz / per_vector_size;

  assert(start_index < end_index && end_index <= total_vectors && "Invalid range");

  size_t n = end_index - start_index;
  *d_out = d;
  *n_out = n;

  float* result = new float[n * d];

  fseek(f, start_index * per_vector_size, SEEK_SET);

  size_t vectors_left = n;
  size_t offset = 0;

  while (vectors_left > 0) {
      size_t current_batch_size = std::min(batch_size, vectors_left);
      std::vector<uint8_t> temp(per_vector_size * current_batch_size);

      size_t nr = fread(temp.data(), sizeof(uint8_t), temp.size(), f);
      assert(nr == temp.size() || !"Could not read batch");

      for (size_t i = 0; i < current_batch_size; i++) {
          const uint8_t* data_ptr = temp.data() + i * per_vector_size + sizeof(int);
          for (size_t j = 0; j < d; j++) {
              result[offset + i * d + j] = static_cast<float>(data_ptr[j]);
          }
      }

      offset += current_batch_size * d;
      vectors_left -= current_batch_size;
  }

  fclose(f);
  return result;
}

float *read_vectors(const std::string &filepath, int num, size_t *d_out, size_t *n_out) {
    if (filepath.size() >= 6) {
        std::string suffix2 = filepath.substr(filepath.size() - 6); 
        if (suffix2 == ".fvecs") {
            return fvecs_read(filepath.c_str(), d_out, n_out);
        }
        if (suffix2 == ".bvecs") {
            return bvecs_read(filepath.c_str(), num, d_out, n_out);
        }
    }
    std::cerr << "Unsupported vector file format: " << filepath << std::endl;
    std::exit(1);
}