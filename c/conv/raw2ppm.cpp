#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>

bool ppm_write(char *path, int w, int h, double *data) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) {
    return false;
  }
  fprintf(f, "P6 %d %d %d ", w, h, 255);
  char *buf = new char[3 * w * h];
  for (int i = 0; i < 3 * w * h; ++i) {
    buf[i] = (unsigned char)(255.0 * fmax(0.0, fmin(1.0, data[i])) + 0.5);
  }
  fwrite(buf, 1, 3 * w * h, f);
  fclose(f);
  return true;
}

double *raw_read(char *path, int *w, int *h) {
  FILE *f = fopen(path, "rb");
  if (f == NULL) {
    return NULL;
  }
  fseek(f, 0, SEEK_END);
  size_t size = ftell(f);
  fseek(f, 0, SEEK_SET);
  char *bytes = new char[size];
  fread(bytes, 1, size, f);
  fclose(f);
  
  memcpy(w, bytes + 0, 4);
  memcpy(h, bytes + 4, 4);
  
  double *data = new double[3 * *w * *h];
  memcpy(data, bytes + 8, 8 * 3 * *w * *h);
  
  delete[] bytes;
  return data;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    puts("too few argument");
    return EXIT_FAILURE;
  }
  int w, h;
  double *data;
  data = raw_read(argv[1], &w, &h);
  if (data == NULL) {
    puts("failed reading");
    return EXIT_FAILURE;
  }
  if (!ppm_write(argv[2], w, h, data)) {
    puts("failed writing");
    return EXIT_FAILURE;
  }
  delete[] data;
  puts("success!");
  return EXIT_SUCCESS;
}