gcc.exe -std=c11 -Wall -Wno-uninitialized -Ofast -march=native -ffast-math -mfpmath=sse -fopenmp -o a.exe src.c
REM gcc.exe -std=c11 -Wall -Wno-uninitialized -Ofast -march=native -ffast-math -mfpmath=sse -fopenmp -S -masm=intel -o a.asm src.c