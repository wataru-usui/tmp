#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>

void RGBE2Float(unsigned char *colRGBE, float *colFloat)
{

    if((*(colRGBE) == 0) && (*(colRGBE + 1) == 0) &&
       (*(colRGBE + 2) == 0)) { //if it is small
        *(colFloat) = 0;
        *(colFloat + 1) = 0;
        *(colFloat + 2) = 0;
        return;
    }

    int E;
    float f;

    E = *(colRGBE + 3) - 128 - 8;
    f = ldexpf(1.0f, E);

    *(colFloat) = (float(*(colRGBE)) + 0.5f) * f;
    *(colFloat + 1) = (float(*(colRGBE + 1)) + 0.5f) * f;
    *(colFloat + 2) = (float(*(colRGBE + 2)) + 0.5f) * f;
}

float *ReadHDR(std::string nameFile, float *data, int &width,
                          int &height)
{
    FILE *file = fopen(nameFile.c_str(), "rb");

    if(file == NULL) {
        return NULL;
    }

    char tmp[512];

    //Is it a Radiance file?
    fscanf(file, "%s\n", tmp);

    if(strcmp(tmp, "#?RADIANCE") != 0 && strcmp(tmp, "#?RGBE") != 0) {
        return NULL;
    }

    while(true) { //Reading Radiance Header
        std::string line = "";

        while(true) { //read property line
            char *tmp2 = fgets(tmp, 512, file);

            if(tmp2 == NULL) {
                return NULL;
            }

            line += tmp2;
            size_t pos = line.find("\n");

            if(pos != std::string::npos) {
                break;
            }
        }

        if(line.compare("\n") == 0) {
            break;
        }

        //Properties:
        if(line.find("FORMAT") != std::string::npos) { //Format
            if(line.find("32-bit_rle_rgbe") == std::string::npos) {
                return NULL;
            }
        }

        if(line.find("EXPOSURE=") != std::string::npos) { //Exposure
            //TODO: ...
        }
    }

    //width and height
    fscanf(file, "-Y %d +X %d", &height, &width);
    fgetc(file);

    if(data == NULL) {
        data = new float[width * height * 3];
    }

    //File size
    long int s_cur = ftell(file);
    fseek(file, 0 , SEEK_END);
    long int s_end = ftell(file);
    fseek(file, s_cur, SEEK_SET);
    int total = s_end - s_cur;

#ifdef PIC_DEBUG
    printf("%d %d\n", total, width * height * 4);
#endif

    //Compressed?
    if(total == (width * height * 4)) { //uncompressed
        unsigned char colRGBE[4];

        int c = 0;

        for(int i = 0; i < width; i++) {
            for(int j = 0; j < height; j++) {
                fread(colRGBE, 1, 4, file);
                RGBE2Float(colRGBE, &data[c]);
                c += 3;
            }
        }
    } else { //RLE compressed
        unsigned char *buffer = new unsigned char[total];
        fread(buffer, sizeof(unsigned char)*total, 1, file);

        int line_width3 = width * 3;
        int line_width4 = width * 4;

        unsigned char *buffer_line_start;
        unsigned char *buffer_line = new unsigned char[line_width4];
        int c = 4;
        int c_buffer_line = 0;

        //for each line
        for(int i = 0; i < height; i++) {
            buffer_line_start = &buffer[c - 4];

            int width_check  = buffer_line_start[2];
            int width_check2 = buffer_line_start[3];

            bool b1 = buffer_line_start[0] != 2;
            bool b2 = buffer_line_start[1] != 2;
            bool b3 = width_check  != (width >> 8); 
            bool b4 = width_check2 != (width & 0xFF);

            if(b1 || b2 || b3 || b4) {
                #ifdef PIC_DEBUG
                    printf("ReadHDR ERROR: the file is not a RLE encoded .hdr file.\n");
                #endif

                fclose(file);

                return NULL;
            }

            for(int j = 0; j < 4; j++) {
                int k = 0;

                //decompression of a single channel line
                while(k < width) {
                    int num = buffer[c];

                    if(num > 128) {
                        num -= 128;

                        for(int l = k; l < (k + num); l++) {
                            buffer_line[l * 4 + j] = buffer[c + 1];
                        }

                        c += 2;
                        k += num;
                    } else {
                        for(int l = 0; l < num; l++) {
                            buffer_line[(l + k) * 4 + j] = buffer[c + 1 + l];
                        }

                        c += num + 1;
                        k += num;
                    }
                }
            }

            //From RGBE to Float
            for(int j = 0; j < width; j++) {
                RGBE2Float(&buffer_line[j * 4], &data[c_buffer_line + j * 3]);
            }

            c += 4;
            c_buffer_line += line_width3;
        }

        delete[] buffer_line;
        delete[] buffer;
    }

    fclose(file);
    return data;
}

bool write(char *path, float *data, int w, int h) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) {
    return false;
  }
  fwrite(&w, 4, 1, f);
  fwrite(&h, 4, 1, f);
  double *data1 = new double[3 * w * h];
  for (int i = 0; i < 3 * w * h; ++i) {
    data1[i] = data[i];
  }
  fwrite(data1, 8, 3 * w * h, f);
  delete[] data1;
  fclose(f);
  return true;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    puts("too few argument");
    return EXIT_FAILURE;
  }
  int w, h;
  float *data = ReadHDR(argv[1], NULL, w, h);
  if (data == NULL) {
    puts("failed reading");
    return EXIT_FAILURE;
  }
  if (!write(argv[2], data, w, h)) {
    puts("failed writing");
    return EXIT_FAILURE;
  }
  puts("success!");
  delete[] data;
  return EXIT_SUCCESS;
}