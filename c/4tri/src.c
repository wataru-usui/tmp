#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <immintrin.h>

#define CAST(t, x) (t)(x)
#define ALLOC(t, n) malloc(sizeof(t) * (n))
#define PUSH(arr, ptr, val) arr[(ptr)++] = val
#define POP(arr, ptr) arr[--(ptr)]
#define DEQ(arr, ptr) arr[(ptr)++]

#define PI 3.141592653589793
#define EPS 1.0e-9

#define NUM_THREADS 1
#define SEED 0

typedef struct vec VEC;
typedef struct aabb AABB;
typedef struct bvh_node BVH_NODE;
typedef struct sort_entry SORT_ENTRY;

typedef struct {
  uint64_t s0, s1;
} PRNG_STATE;

struct vec {
	double x, y, z;
};

struct aabb {
	VEC min, max;
};

struct bvh_node {
	AABB b;
	int n, i;
	int c[2];
};

struct sort_entry {
	int *fs;
	double *vs;
	int i;
};

int max(int a, int b) {
	return a > b ? a : b;
}

int min(int a, int b) {
	return a < b ? a : b;
}
int bitcount(int x) {
	x = (x & 0x55555555) + (x >> 0x01 & 0x55555555);
	x = (x & 0x33333333) + (x >> 0x02 & 0x33333333);
	x = (x & 0x0f0f0f0f) + (x >> 0x04 & 0x0f0f0f0f);
	x = (x & 0x00ff00ff) + (x >> 0x08 & 0x00ff00ff);
	x = (x & 0x0000ffff) + (x >> 0x10 & 0x0000ffff);
	return x;
}

uint64_t sm64_next(uint64_t *s) {
  uint64_t x = *s += 0x9e3779b97f4a7c15;
  x = 0xbf58476d1ce4e5b9 * (x ^ (x >> 30));
  x = 0x94d049bb133111eb * (x ^ (x >> 27));
  x = x ^ (x >> 31);
  return x;
}
//xoroshiro128+, reference http://prng.di.unimi.it/#intro
uint64_t rotl(uint64_t a, uint64_t b) {
  return (a << b) | (a >> (64 - b));
}
uint64_t prng_next(PRNG_STATE *s) {
  PRNG_STATE a = *s;
  uint64_t ret = a.s0 + a.s1;
  a.s1 ^= a.s0;
  s->s0 = rotl(a.s0, 55) ^ a.s1 ^ (a.s1 << 14);
  s->s1 = rotl(a.s1, 36);
  return ret;
}
double prng_db(PRNG_STATE *s) {
  return 0x1.0p-53 * (prng_next(s) >> 11);
}
PRNG_STATE prng_seed(uint64_t seed) {
  PRNG_STATE s = {sm64_next(&seed), sm64_next(&seed)};
  return s;
}
PRNG_STATE prng_jump(PRNG_STATE s) {
  PRNG_STATE jump = {0xbeac0467eba5facb, 0xd86b048b86aa9922};
  PRNG_STATE next = {0, 0};
  for (int i = 0; i < 64; ++i) {
    if (jump.s0 & 1 << i) {
      next.s0 ^= s.s0;
      next.s1 ^= s.s1;
    }
    prng_next(&s);
  }
  for (int i = 0; i < 64; ++i) {
    if (jump.s1 & 1 << i) {
      next.s0 ^= s.s0;
      next.s1 ^= s.s1;
    }
    prng_next(&s);
  }
  return next;
}

int Compact1By1(int x) {
  x &= 0x55555555;					// x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x = (x ^ (x >> 1)) & 0x33333333;	// x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >> 2)) & 0x0f0f0f0f;	// x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >> 4)) & 0x00ff00ff;	// x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >> 8)) & 0x0000ffff;	// x = ---- ---- ---- ---- fedc ba98 7654 3210
  return x;
}
int DecodeMorton2X(int code) {
  return Compact1By1(code >> 0);
}
int DecodeMorton2Y(int code) {
  return Compact1By1(code >> 1);
}

double time_get() {
  return clock() / CAST(double, CLOCKS_PER_SEC);
}

VEC vec_init(double x, double y, double z) {
	VEC v = {x, y, z};
	return v;
}

VEC vec_init1(double x) {
	return vec_init(x, x, x);
}
VEC vec_add(VEC a, VEC b) {
	return vec_init(a.x + b.x, a.y + b.y, a.z + b.z);
}

VEC vec_sub(VEC a, VEC b) {
	return vec_init(a.x - b.x, a.y - b.y, a.z - b.z);
}

VEC vec_scale(double a, VEC b) {
	return vec_init(a * b.x, a * b.y, a * b.z);
}

double vec_dot(VEC a, VEC b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

double vec_len(VEC x) {
	return sqrt(vec_dot(x, x));
}

VEC vec_norm(VEC x) {
	return vec_scale(1.0 / vec_len(x), x);
}

VEC vec_cross(VEC a, VEC b) {
	return vec_init(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
VEC vec_min(VEC a, VEC b) {
	return vec_init(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z));
}
VEC vec_max(VEC a, VEC b) {
	return vec_init(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z));
}

void vec_pr(VEC a, FILE *o) {
	fprintf(o, "%f %f %f\n", a.x, a.y, a.z);
}

void read_scene(char *path, int *nv, double **vs, int *nf, int **fs) {
  FILE *f = fopen(path, "r");
  if (f == NULL) {
    printf("cannot open file: %s\n", path);
    exit(0);
  }
  
  int buf_size = 256, id_size = 256;
  char buf[buf_size], id[id_size];
  
  double t[3], m[9];
  
  while (true) {
    if (fgets(buf, buf_size, f) == NULL) {
      break;
    }
    sscanf(buf, "%s ", id);
    if (strcmp(id, "translate") == 0) {
      sscanf(buf, "%*s %lf %lf %lf ", t + 0, t + 1, t + 2);
    } else if (strcmp(id, "transform") == 0) {
      sscanf(buf, "%*s %lf %lf %lf %lf %lf %lf %lf %lf %lf ", m + 0, m + 1, m + 2, m + 3, m + 4, m + 5, m + 6, m + 7, m + 8);
    } else if (strcmp(id, "verticies") == 0) {
      sscanf(buf, "%*s %i ", nv);
      *vs = ALLOC(double, 3 * *nv);
      for (int i = 0; i < *nv; ++i) {
		  double x[3], y[3];
        fgets(buf, buf_size, f);
        sscanf(buf, "%lf %lf %lf ", x + 0, x + 1, x + 2);
		for (int j = 0; j < 3; ++j) {
			double s = 0.0;
			for (int k = 0; k < 3; ++k) {
				s += m[3 * j + k] * x[k];
			}
			y[j] = s + t[j];
		}
		for (int j = 0; j < 3; ++j) {
			(*vs)[3 * i + j] = y[j];
		}
      }
    } else if (strcmp(id, "faces") == 0) {
      sscanf(buf, "%*s %i ", nf);
      *fs = ALLOC(int, 3 * *nf);
      for (int i = 0; i < *nf; ++i) {
        fgets(buf, buf_size, f);
        sscanf(buf, "%*s %*s %*s %*s %i %i %i ", *fs + 3 * i + 0, *fs + 3 * i + 1, *fs + 3 * i + 2);
      }
    }
  }
  
  fclose(f);
}

void img_save(char *path, int w, int h, double *img) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) {
    printf("cannot open file: %s\n", path);
    exit(0);
  }
  fprintf(f, "Pf %d %d -1 ", w, h);
  float *data = ALLOC(float, w * h);
  for (int i = 0; i < w * h; ++i) {
	  data[i] = CAST(float, img[i]);
  }
  fwrite(data, sizeof(float), w * h, f);
  free(data);
  fclose(f);
}

void face_ld_vs(int *fs, double *vs, int i, VEC *a, VEC *b, VEC *c) {
	*a = vec_init(vs[3 * fs[3 * i + 0] + 0], vs[3 * fs[3 * i + 0] + 1], vs[3 * fs[3 * i + 0] + 2]);
	*b = vec_init(vs[3 * fs[3 * i + 1] + 0], vs[3 * fs[3 * i + 1] + 1], vs[3 * fs[3 * i + 1] + 2]);
	*c = vec_init(vs[3 * fs[3 * i + 2] + 0], vs[3 * fs[3 * i + 2] + 1], vs[3 * fs[3 * i + 2] + 2]);
}

bool tri_isect(VEC v0, VEC v1, VEC v2, VEC o, VEC d, double *t) {
  double eps = EPS;
  VEC e0 = vec_sub(v1, v0);
  VEC e1 = vec_sub(v2, v0);
  VEC pv = vec_cross(d, e1);
  double deti = vec_dot(e0, pv);
  if (fabs(deti) < eps) {
    return false;
  }
  deti = 1.0 / deti;
  VEC tv = vec_sub(o, v0);
  double u = deti * vec_dot(tv, pv);
  if (u < 0.0 || u > 1.0) {
    return false;
  }
  VEC qv = vec_cross(tv, e0);
  double v = deti * vec_dot(d, qv);
  if (v < 0.0 || u + v > 1.0) {
    return false;
  }
  double k = deti * vec_dot(e1, qv);
  if (k < eps || k > *t) {
	  return false;
  }
  *t = k;
  return true;
}

VEC tri_norm(VEC v0, VEC v1, VEC v2) {
	return vec_norm(vec_cross(vec_sub(v1, v0), vec_sub(v2, v0)));
}

int tri_isect_simd(__m256d v0[3], __m256d v1[3], __m256d v2[3], __m256d o[3], __m256d d[3], __m256d *td) {
	__m256d e0[3];
	e0[0] = _mm256_sub_pd(v1[0], v0[0]);
	e0[1] = _mm256_sub_pd(v1[1], v0[1]);
	e0[2] = _mm256_sub_pd(v1[2], v0[2]);
	
	__m256d e1[3];
	e1[0] = _mm256_sub_pd(v2[0], v0[0]);
	e1[1] = _mm256_sub_pd(v2[1], v0[1]);
	e1[2] = _mm256_sub_pd(v2[2], v0[2]);
	
	__m256d pv[3];
	pv[0] = _mm256_sub_pd(_mm256_mul_pd(d[1], e1[2]), _mm256_mul_pd(d[2], e1[1]));
	pv[1] = _mm256_sub_pd(_mm256_mul_pd(d[2], e1[0]), _mm256_mul_pd(d[0], e1[2]));
	pv[2] = _mm256_sub_pd(_mm256_mul_pd(d[0], e1[1]), _mm256_mul_pd(d[1], e1[0]));
	
	__m256d deti = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(e0[0], pv[0]), _mm256_mul_pd(e0[1], pv[1])), _mm256_mul_pd(e0[2], pv[2]));
	
	__m256d flag = _mm256_cmp_pd(_mm256_mul_pd(deti, deti), _mm256_set1_pd(EPS), _CMP_GT_OQ);
	
	deti = _mm256_div_pd(_mm256_set1_pd(1.0), deti);
	
	__m256d tv[3];
	tv[0] = _mm256_sub_pd(o[0], v0[0]);
	tv[1] = _mm256_sub_pd(o[1], v0[1]);
	tv[2] = _mm256_sub_pd(o[2], v0[2]);
	
	__m256d u = _mm256_mul_pd(deti, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(tv[0], pv[0]), _mm256_mul_pd(tv[1], pv[1])), _mm256_mul_pd(tv[2], pv[2])));
	
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(u, _mm256_set1_pd(0.0), _CMP_GT_OQ));
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(u, _mm256_set1_pd(1.0), _CMP_LT_OQ));
	
	__m256d qv[3];
	qv[0] = _mm256_sub_pd(_mm256_mul_pd(tv[1], e0[2]), _mm256_mul_pd(tv[2], e0[1]));
	qv[1] = _mm256_sub_pd(_mm256_mul_pd(tv[2], e0[0]), _mm256_mul_pd(tv[0], e0[2]));
	qv[2] = _mm256_sub_pd(_mm256_mul_pd(tv[0], e0[1]), _mm256_mul_pd(tv[1], e0[0]));
	
	__m256d v = _mm256_mul_pd(deti, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(d[0], qv[0]), _mm256_mul_pd(d[1], qv[1])), _mm256_mul_pd(d[2], qv[2])));
	
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(v, _mm256_set1_pd(0.0), _CMP_GT_OQ));
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(_mm256_add_pd(u, v), _mm256_set1_pd(1.0), _CMP_LT_OQ));
	
	*td = _mm256_mul_pd(deti, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(e1[0], qv[0]), _mm256_mul_pd(e1[1], qv[1])), _mm256_mul_pd(e1[2], qv[2])));
	
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(*td, _mm256_set1_pd(EPS), _CMP_GT_OQ));
	
	return _mm256_movemask_pd(flag);
}

int tri_isect_simd2(__m256d A[3], __m256d B[3], __m256d C[3], __m256d O[3], __m256d D[3], __m256d *t) {
	__m256d T[3];
	T[0] = _mm256_sub_pd(O[0], A[0]);
	T[1] = _mm256_sub_pd(A[1], O[1]);
	T[2] = _mm256_sub_pd(O[2], A[2]);
	
	__m256d E[3];
	E[0] = _mm256_sub_pd(B[0], A[0]);
	E[1] = _mm256_sub_pd(B[1], A[1]);
	E[2] = _mm256_sub_pd(B[2], A[2]);
	
	__m256d F[3];
	F[0] = _mm256_sub_pd(C[0], A[0]);
	F[1] = _mm256_sub_pd(C[1], A[1]);
	F[2] = _mm256_sub_pd(C[2], A[2]);
	
	__m256d G[3];
	G[0] = _mm256_sub_pd(_mm256_mul_pd(E[1], F[2]), _mm256_mul_pd(E[2], F[1]));
	G[1] = _mm256_sub_pd(_mm256_mul_pd(D[1], F[2]), _mm256_mul_pd(D[2], F[1]));
	G[2] = _mm256_sub_pd(_mm256_mul_pd(D[1], E[2]), _mm256_mul_pd(D[2], E[1]));
	
	__m256d deti = _mm256_add_pd(_mm256_sub_pd(
	_mm256_mul_pd(D[0], G[0]),
	_mm256_mul_pd(E[0], G[1])),
	_mm256_mul_pd(F[0], G[2]));
	
	__m256d flag = _mm256_cmp_pd(_mm256_mul_pd(deti, deti), _mm256_set1_pd(EPS), _CMP_GT_OQ);
	
	deti = _mm256_div_pd(_mm256_set1_pd(1.0), deti);
	
	*t = _mm256_sub_pd(_mm256_set1_pd(0.0), _mm256_mul_pd(deti, _mm256_add_pd(_mm256_add_pd(
	_mm256_mul_pd(T[0], G[0]),
	_mm256_mul_pd(T[1], _mm256_sub_pd(_mm256_mul_pd(E[0], F[2]), _mm256_mul_pd(E[2], F[0])))),
	_mm256_mul_pd(T[2], _mm256_sub_pd(_mm256_mul_pd(E[0], F[1]), _mm256_mul_pd(E[1], F[0]))))));
	
	__m256d u = _mm256_sub_pd(_mm256_set1_pd(0.0), _mm256_mul_pd(deti, _mm256_add_pd(_mm256_add_pd(
	_mm256_mul_pd(T[0], G[1]),
	_mm256_mul_pd(T[1], _mm256_sub_pd(_mm256_mul_pd(D[0], F[2]), _mm256_mul_pd(D[2], F[0])))),
	_mm256_mul_pd(T[2], _mm256_sub_pd(_mm256_mul_pd(D[0], F[1]), _mm256_mul_pd(D[1], F[0]))))));
	
	__m256d v = _mm256_mul_pd(deti, _mm256_add_pd(_mm256_add_pd(
	_mm256_mul_pd(T[0], G[2]),
	_mm256_mul_pd(T[1], _mm256_sub_pd(_mm256_mul_pd(D[0], E[2]), _mm256_mul_pd(D[2], E[0])))),
	_mm256_mul_pd(T[2], _mm256_sub_pd(_mm256_mul_pd(D[0], E[1]), _mm256_mul_pd(D[1], E[0])))));
	
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(*t, _mm256_set1_pd(EPS), _CMP_GT_OQ));
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(u, _mm256_set1_pd(0.0), _CMP_GT_OQ));
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(u, _mm256_set1_pd(1.0), _CMP_LT_OQ));
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(v, _mm256_set1_pd(0.0), _CMP_GT_OQ));
	flag = _mm256_and_pd(flag, _mm256_cmp_pd(_mm256_add_pd(u, v), _mm256_set1_pd(1.0), _CMP_LT_OQ));
	
	return _mm256_movemask_pd(flag);
}

int tri_isect_simd_wrap(double *vs, int nf, int *fs, VEC o, VEC d) {
	__m256d _o[3];
	_o[0] = _mm256_set1_pd(o.x);
	_o[1] = _mm256_set1_pd(o.y);
	_o[2] = _mm256_set1_pd(o.z);
	
	__m256d _d[3];
	_d[0] = _mm256_set1_pd(d.x);
	_d[1] = _mm256_set1_pd(d.y);
	_d[2] = _mm256_set1_pd(d.z);
	
	int id = -1;
	double td = INFINITY;
	for (int j = 0; j < nf; j += 4) {
		__m256d v[3][3];
		v[0][0] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 0] + 0], vs[3 * fs[3 * (j + 2) + 0] + 0], vs[3 * fs[3 * (j + 1) + 0] + 0], vs[3 * fs[3 * (j + 0) + 0] + 0]);
		v[0][1] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 0] + 1], vs[3 * fs[3 * (j + 2) + 0] + 1], vs[3 * fs[3 * (j + 1) + 0] + 1], vs[3 * fs[3 * (j + 0) + 0] + 1]);
		v[0][2] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 0] + 2], vs[3 * fs[3 * (j + 2) + 0] + 2], vs[3 * fs[3 * (j + 1) + 0] + 2], vs[3 * fs[3 * (j + 0) + 0] + 2]);
		v[1][0] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 1] + 0], vs[3 * fs[3 * (j + 2) + 1] + 0], vs[3 * fs[3 * (j + 1) + 1] + 0], vs[3 * fs[3 * (j + 0) + 1] + 0]);
		v[1][1] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 1] + 1], vs[3 * fs[3 * (j + 2) + 1] + 1], vs[3 * fs[3 * (j + 1) + 1] + 1], vs[3 * fs[3 * (j + 0) + 1] + 1]);
		v[1][2] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 1] + 2], vs[3 * fs[3 * (j + 2) + 1] + 2], vs[3 * fs[3 * (j + 1) + 1] + 2], vs[3 * fs[3 * (j + 0) + 1] + 2]);
		v[2][0] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 2] + 0], vs[3 * fs[3 * (j + 2) + 2] + 0], vs[3 * fs[3 * (j + 1) + 2] + 0], vs[3 * fs[3 * (j + 0) + 2] + 0]);
		v[2][1] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 2] + 1], vs[3 * fs[3 * (j + 2) + 2] + 1], vs[3 * fs[3 * (j + 1) + 2] + 1], vs[3 * fs[3 * (j + 0) + 2] + 1]);
		v[2][2] = _mm256_set_pd(vs[3 * fs[3 * (j + 3) + 2] + 2], vs[3 * fs[3 * (j + 2) + 2] + 2], vs[3 * fs[3 * (j + 1) + 2] + 2], vs[3 * fs[3 * (j + 0) + 2] + 2]);
		
		__m256d t;
		int flag = tri_isect_simd(v[0], v[1], v[2], _o, _d, &t);
		
		id = flag & 1 && t[0] < td ? j + 0 : id;
		td = flag & 1 && t[0] < td ? t[0] : td;
		id = flag & 2 && t[1] < td ? j + 1 : id;
		td = flag & 2 && t[1] < td ? t[1] : td;
		id = flag & 4 && t[2] < td ? j + 2 : id;
		td = flag & 4 && t[2] < td ? t[2] : td;
		id = flag & 8 && t[3] < td ? j + 3 : id;
		td = flag & 8 && t[3] < td ? t[3] : td;
	}
	
	return id;
}

AABB aabb_default(void) {
	AABB box = {vec_init1(INFINITY), vec_init1(-INFINITY)};
	return box;
}
double aabb_area(AABB box) {
	VEC d = vec_sub(box.max, box.min);
	return 2.0 * (d.x * d.y + d.y * d.z + d.z * d.x);
}
VEC aabb_center(AABB box) {
	return vec_scale(0.5, vec_add(box.min, box.max));
}
AABB aabb_merge(AABB a, AABB b) {
	AABB box = {vec_min(a.min, b.min), vec_max(a.max, b.max)};
	return box;
}
bool aabb_isect(AABB box, VEC o, VEC d, double *t) {
  double eps = EPS;
  double tmin = -INFINITY, tmax = INFINITY;
  if (fabs(d.x) > eps) {
	double tx0 = (box.min.x - o.x) / d.x;
	double tx1 = (box.max.x - o.x) / d.x;
	tmin = fmax(tmin, fmin(tx0, tx1));
	tmax = fmin(tmax, fmax(tx0, tx1));
	if (tmax < tmin) {
		return false;
	}
  }
  if (fabs(d.y) > eps) {
	double ty0 = (box.min.y - o.y) / d.y;
	double ty1 = (box.max.y - o.y) / d.y;
	tmin = fmax(tmin, fmin(ty0, ty1));
	tmax = fmin(tmax, fmax(ty0, ty1));
	if (tmax < tmin) {
		return false;
	}
  }
  if (fabs(d.z) > eps) {
	double tz0 = (box.min.z - o.z) / d.z;
	double tz1 = (box.max.z - o.z) / d.z;
	tmin = fmax(tmin, fmin(tz0, tz1));
	tmax = fmin(tmax, fmax(tz0, tz1));
	if (tmax < tmin) {
		return false;
	}
  }
  *t = fmax(0.0, tmin);
  return tmax > eps;
}
AABB tri_aabb(VEC v0, VEC v1, VEC v2) {
  AABB box = aabb_default();
  box.min = vec_min(box.min, v0);
  box.min = vec_min(box.min, v1);
  box.min = vec_min(box.min, v2);
  box.max = vec_max(box.max, v0);
  box.max = vec_max(box.max, v1);
  box.max = vec_max(box.max, v2);
  return box;
}
int bvh_cmp_x(void const *_a, void const *_b) {
	SORT_ENTRY const *a = _a;
	SORT_ENTRY const *b = _b;
	VEC av0, av1, av2;
	VEC bv0, bv1, bv2;
	face_ld_vs(a->fs, a->vs, a->i, &av0, &av1, &av2);
	face_ld_vs(b->fs, b->vs, b->i, &bv0, &bv1, &bv2);
	AABB ab = tri_aabb(av0, av1, av2);
	AABB bb = tri_aabb(bv0, bv1, bv2);
	VEC ac = aabb_center(ab);
	VEC bc = aabb_center(bb);
	return ac.x < bc.x ? -1 : ac.x > bc.x ? 1 : 0;
}
int bvh_cmp_y(void const *_a, void const *_b) {
	SORT_ENTRY const *a = _a;
	SORT_ENTRY const *b = _b;
	VEC av0, av1, av2;
	VEC bv0, bv1, bv2;
	face_ld_vs(a->fs, a->vs, a->i, &av0, &av1, &av2);
	face_ld_vs(b->fs, b->vs, b->i, &bv0, &bv1, &bv2);
	AABB ab = tri_aabb(av0, av1, av2);
	AABB bb = tri_aabb(bv0, bv1, bv2);
	VEC ac = aabb_center(ab);
	VEC bc = aabb_center(bb);
	return ac.y < bc.y ? -1 : ac.y > bc.y ? 1 : 0;
}
int bvh_cmp_z(void const *_a, void const *_b) {
	SORT_ENTRY const *a = _a;
	SORT_ENTRY const *b = _b;
	VEC av0, av1, av2;
	VEC bv0, bv1, bv2;
	face_ld_vs(a->fs, a->vs, a->i, &av0, &av1, &av2);
	face_ld_vs(b->fs, b->vs, b->i, &bv0, &bv1, &bv2);
	AABB ab = tri_aabb(av0, av1, av2);
	AABB bb = tri_aabb(bv0, bv1, bv2);
	VEC ac = aabb_center(ab);
	VEC bc = aabb_center(bb);
	return ac.z < bc.z ? -1 : ac.z > bc.z ? 1 : 0;
}

void sah_build_bvh_recur(SORT_ENTRY *ifs, double *la, double *ra, BVH_NODE *nodes, int s, int e, int *ptr) {
	double time_tri = 1.0, time_aabb = 1.0;
	AABB box = aabb_default();
	for (int i = s; i < e; ++i) {
		VEC v0, v1, v2;
		face_ld_vs(ifs[i].fs, ifs[i].vs, ifs[i].i, &v0, &v1, &v2);
		box = aabb_merge(box, tri_aabb(v0, v1, v2));
	}
	double area = aabb_area(box);
	double cost = time_tri * (e - s);
	int axis = -1, idx = -1;
	for (int i = 0; i < 3; ++i) {
		switch (i) {
			case 0 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_x);
			} break;
			case 1 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_y);
			} break;
			case 2 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_z);
			} break;
		}
		AABB lb = aabb_default();
		for (int j = 0; j < e - s - 1; ++j) {
			VEC v0, v1, v2;
			face_ld_vs(ifs[s + j].fs, ifs[s + j].vs, ifs[s + j].i, &v0, &v1, &v2);
			lb = aabb_merge(lb, tri_aabb(v0, v1, v2));
			la[j] = aabb_area(lb);
		}
		AABB rb = aabb_default();
		for (int j = e - s - 1; j > 0; --j) {
			VEC v0, v1, v2;
			face_ld_vs(ifs[s + j].fs, ifs[s + j].vs, ifs[s + j].i, &v0, &v1, &v2);
			rb = aabb_merge(rb, tri_aabb(v0, v1, v2));
			ra[j] = aabb_area(rb);
		}
		for (int j = 0; j < e - s - 1; ++j) {
			double cost1 = 2.0 * time_aabb + (la[j] * (j + 1) + ra[j + 1] * (e - s - j - 1)) * time_tri / area;
			if (cost1 < cost) {
				cost = cost1;
				axis = i;
				idx = s + j + 1;
			}
		}
	}
	if (axis == -1) {
		BVH_NODE node;
		node.b = box;
		node.n = s - e;
		node.i = s;
		nodes[*ptr] = node;
		++*ptr;
	} else {
		switch (axis) {
			case 0 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_x);
			} break;
			case 1 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_y);
			} break;
			case 2 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_z);
			} break;
		}
		int tmp = *ptr;
		BVH_NODE node;
		node.b = box;
		++*ptr;
		node.c[0] = *ptr;
		sah_build_bvh_recur(ifs, la, ra, nodes, s, idx, ptr);
		node.c[1] = *ptr;
		sah_build_bvh_recur(ifs, la, ra, nodes, idx, e, ptr);
		nodes[tmp] = node;
	}
}

int sah_build_bvh(int nf, int **fs, double *vs, BVH_NODE *nodes) {
	SORT_ENTRY *ifs = ALLOC(SORT_ENTRY, nf);
	for (int i = 0; i < nf; ++i) {
		SORT_ENTRY a = {*fs, vs, i};
		ifs[i] = a;
	}
	
	double *la = ALLOC(double, nf);
	double *ra = ALLOC(double, nf);
	
	int ptr = 0;
	sah_build_bvh_recur(ifs, la, ra, nodes, 0, nf, &ptr);
	
	free(la);
	free(ra);
	
	int *fsd = ALLOC(int, 3 * nf);
	for (int i = 0; i < nf; ++i) {
		for (int j = 0; j < 3; ++j) {
			fsd[3 * i + j] = (*fs)[3 * ifs[i].i + j];
		}
	}
	
	free(ifs);
	free(*fs);
	
	*fs = fsd;
	
	return ptr;
}

void bvh_traverse_recurse(int *fs, double *vs, BVH_NODE *nodes, BVH_NODE node, VEC o, VEC d, int *id, double *t) {
	if (node.n < 0) {
		for (int i = 0; i < -node.n; ++i) {
			VEC v0, v1, v2;
			face_ld_vs(fs, vs, node.i + i, &v0, &v1, &v2);
			
			if (tri_isect(v0, v1, v2, o, d, t)) {
				*id = node.i + i;
			}
		}
	} else {
		BVH_NODE l = nodes[node.c[0]], r = nodes[node.c[1]];
		bool lh, rh;
		double lt, rt;
		lh = aabb_isect(l.b, o, d, &lt);
		rh = aabb_isect(r.b, o, d, &rt);
		
		if (lh && rh) {
			if (lt < rt) {
				bvh_traverse_recurse(fs, vs, nodes, l, o, d, id, t);
				if (*t < rt) {
					return;
				} else {
					bvh_traverse_recurse(fs, vs, nodes, r, o, d, id, t);
				}
			} else {
				bvh_traverse_recurse(fs, vs, nodes, r, o, d, id, t);
				if (*t < lt) {
					return;
				} else {
					bvh_traverse_recurse(fs, vs, nodes, l, o, d, id, t);
				}
			}
		} else if (lh) {
			bvh_traverse_recurse(fs, vs, nodes, l, o, d, id, t);
		} else if (rh) {
			bvh_traverse_recurse(fs, vs, nodes, r, o, d, id, t);
		} else {
			
		}
	}
}

int bvh_traverse(int *fs, double *vs, BVH_NODE *nodes, VEC o, VEC d) {
	int id = -1;
	double t = INFINITY;
	
	BVH_NODE node = nodes[0];
	double tmp;
	if (aabb_isect(node.b, o, d, &tmp)) {
		bvh_traverse_recurse(fs, vs, nodes, node, o, d, &id, &t);
	}
	
	return id;
}

typedef struct bvh_simd_build_entry BVH_SIMD_BUILD_ENTRY;

struct bvh_simd_build_entry {
	int s, e;
};

typedef struct bvh_node_simd1 BVH_NODE_SIMD1;
typedef struct bvh_simd_build_entry1 BVH_SIMD_BUILD_ENTRY1;

struct bvh_node_simd1 {
	double b[2][3][4];
	int c[4], nb, nt, it;
};

struct bvh_simd_build_entry1 {
	int s, e;
};

void bvh_simd_build1_inter(SORT_ENTRY *ifs, double *as[2], int s, int e, int d, int *axis, int *pivot) {
	double cost = INFINITY;
	for (int i = 0; i < 3; ++i) {
		switch (i) {
			case 0 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_x);
			} break;
			case 1 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_y);
			} break;
			case 2 : {
				qsort(ifs + s, e - s, sizeof(SORT_ENTRY), bvh_cmp_z);
			} break;
		}
		AABB box;
		box = aabb_default();
		for (int j = s; j < e; j += 4) {
			for (int k = 0; k < 4; ++k) {
				VEC v0, v1, v2;
				face_ld_vs(ifs[j + k].fs, ifs[j + k].vs, ifs[j + k].i, &v0, &v1, &v2);
				box = aabb_merge(box, tri_aabb(v0, v1, v2));
			}
			as[0][(j - s) / 4] = aabb_area(box);
		}
		box = aabb_default();
		for (int j = e - 4; j >= s; j -= 4) {
			for (int k = 3; k >= 0; --k) {
				VEC v0, v1, v2;
				face_ld_vs(ifs[j + k].fs, ifs[j + k].vs, ifs[j + k].i, &v0, &v1, &v2);
				box = aabb_merge(box, tri_aabb(v0, v1, v2));
			}
			as[1][(j - s) / 4] = aabb_area(box);
		}
		for (int j = s + d; j <= e - d; j += 4) {
			double costd = as[0][(j - s) / 4 - 1] * (j - s) + as[1][(j - s) / 4] * (e - j);
			if (costd < cost) {
				cost = costd;
				*axis = i;
				*pivot = j;
			}
		}
	}
}

void bvh_simd_build1(SORT_ENTRY *ifs, SORT_ENTRY *ifsd, double *as[2], BVH_NODE_SIMD1 *nodes, BVH_SIMD_BUILD_ENTRY *stack, int s, int e, int *ptr_deq, int *ptr_enq) {
	int ranges[5] = {s, s + 4, s + 8, s + 12, e};
	if (e - s > 16) {
		int axes[3];
		bvh_simd_build1_inter(ifs, as, ranges[0], ranges[4], 8, axes + 0, ranges + 2);
		switch (axes[0]) {
			case 0 : {
				qsort(ifs + ranges[0], ranges[4] - ranges[0], sizeof(SORT_ENTRY), bvh_cmp_x);
			} break;
			case 1 : {
				qsort(ifs + ranges[0], ranges[4] - ranges[0], sizeof(SORT_ENTRY), bvh_cmp_y);
			} break;
			case 2 : {
				qsort(ifs + ranges[0], ranges[4] - ranges[0], sizeof(SORT_ENTRY), bvh_cmp_z);
			} break;
		}
		ranges[1] = ranges[0] + 4;
		if (ranges[2] - ranges[0] > 8) {
			bvh_simd_build1_inter(ifs, as, ranges[0], ranges[2], 4, axes + 1, ranges + 1);
			switch (axes[1]) {
				case 0 : {
					qsort(ifs + ranges[0], ranges[2] - ranges[0], sizeof(SORT_ENTRY), bvh_cmp_x);
				} break;
				case 1 : {
					qsort(ifs + ranges[0], ranges[2] - ranges[0], sizeof(SORT_ENTRY), bvh_cmp_y);
				} break;
				case 2 : {
					qsort(ifs + ranges[0], ranges[2] - ranges[0], sizeof(SORT_ENTRY), bvh_cmp_z);
				} break;
			}
		}
		ranges[3] = ranges[2] + 4;
		if (ranges[4] - ranges[2] > 8) {
			bvh_simd_build1_inter(ifs, as, ranges[2], ranges[4], 4, axes + 2, ranges + 3);
			switch (axes[2]) {
				case 0 : {
					qsort(ifs + ranges[2], ranges[4] - ranges[2], sizeof(SORT_ENTRY), bvh_cmp_x);
				} break;
				case 1 : {
					qsort(ifs + ranges[2], ranges[4] - ranges[2], sizeof(SORT_ENTRY), bvh_cmp_y);
				} break;
				case 2 : {
					qsort(ifs + ranges[2], ranges[4] - ranges[2], sizeof(SORT_ENTRY), bvh_cmp_z);
				} break;
			}
		}
	}
	int num = min(4, (e - s) / 4);
	int acc = 0;
	int cnt = 0;
	int rangesd[5];
	rangesd[4] = e;
	int numbox = 0;
	for (int i = 0; i < num; ++i) {
		int size = ranges[i + 1] - ranges[i];
		if (size > 4) {
			memcpy(ifsd + acc, ifs + ranges[i], sizeof(SORT_ENTRY) * size);
			rangesd[cnt] = s + acc;
			++cnt;
			acc += size;
			++numbox;
		}
	}
	int numtri = 0;
	for (int i = 0; i < num; ++i) {
		int size = ranges[i + 1] - ranges[i];
		if (size == 4) {
			memcpy(ifsd + acc, ifs + ranges[i], sizeof(SORT_ENTRY) * size);
			rangesd[cnt] = s + acc;
			++cnt;
			acc += size;
			++numtri;
		}
	}
	
	memcpy(ifs + s, ifsd, sizeof(SORT_ENTRY) * (e - s));
	
	BVH_NODE_SIMD1 node;
	
	for (int i = 0; i < numbox; ++i) {
		AABB bc = aabb_default();
		int m = rangesd[i], n = rangesd[i + 1];
		for (int j = m; j < n; ++j) {
			VEC v0, v1, v2;
			face_ld_vs(ifs[j].fs, ifs[j].vs, ifs[j].i, &v0, &v1, &v2);
			bc = aabb_merge(bc, tri_aabb(v0, v1, v2));
		}
		node.b[0][0][i] = bc.min.x;
		node.b[0][1][i] = bc.min.y;
		node.b[0][2][i] = bc.min.z;
		node.b[1][0][i] = bc.max.x;
		node.b[1][1][i] = bc.max.y;
		node.b[1][2][i] = bc.max.z;
	}
	
	node.it = rangesd[numbox];
	node.nt = numtri;
	node.nb = numbox;
	
	for (int i = 0; i < numbox; ++i) {
		node.c[i] = *ptr_enq + i + 1;
	}
	
	nodes[*ptr_deq] = node;
	
	for (int i = 0; i < numbox; ++i) {
		BVH_SIMD_BUILD_ENTRY e = {rangesd[i], rangesd[i + 1]};
		PUSH(stack, *ptr_enq, e);
	}
	//printf("r: %d, %d %d %d %d %d\n", *ptr_deq, ranges[0], ranges[1], ranges[2], ranges[3], ranges[4]);
	//printf("rd: %d, %d %d %d %d %d\n", *ptr_deq, rangesd[0], rangesd[1], rangesd[2], rangesd[3], rangesd[4]);
	//printf("node: %d, children: %d %d %d %d\n", *ptr_deq, node.c[0], node.c[1], node.c[2], node.c[3]);
	//printf("nb nt ti: %d %d %d\n", node.nb, node.nt, node.it);
	
}

int bvh_simd_build_wrap1(int nf, int **fs, double *vs, int num_nodes, BVH_NODE_SIMD1 *nodes) {
	SORT_ENTRY *ifs = ALLOC(SORT_ENTRY, nf);
	for (int i = 0; i < nf; ++i) {
		SORT_ENTRY a = {*fs, vs, i};
		ifs[i] = a;
	}
	
	double *areas[2];
	for (int i = 0; i < 2; ++i) {
		areas[i] = ALLOC(double, nf / 4);
	}
	
	BVH_SIMD_BUILD_ENTRY *stack = ALLOC(BVH_SIMD_BUILD_ENTRY, num_nodes);
	SORT_ENTRY *ifsd = ALLOC(SORT_ENTRY, nf);
	
	BVH_SIMD_BUILD_ENTRY query = {0, nf};
	int ptr_deq = 0, ptr_enq = 0;
	
	while (true) {
		bvh_simd_build1(ifs, ifsd, areas, nodes, stack, query.s, query.e, &ptr_deq, &ptr_enq);
		if (ptr_deq == ptr_enq) {
			break;
		}
		query = stack[ptr_deq++];
	}
	free(ifsd);
	free(stack);
	for (int i = 0; i < 2; ++i) {
		free(areas[i]);
	}
	
	int *fsd = ALLOC(int, 3 * nf);
	for (int i = 0; i < nf; ++i) {
		for (int j = 0; j < 3; ++j) {
			fsd[3 * i + j] = (*fs)[3 * ifs[i].i + j];
		}
	}
	
	free(ifs);
	free(*fs);
	
	*fs = fsd;
	
	int numtri = 0;
	for (int i = 0; i <= ptr_deq; ++i) {
		numtri += nodes[i].nt;
	}
	printf("num tris : %d\n", 4 * numtri);
	
	return ptr_deq + 1;
}

int aabb_isect_simd(__m256d box[2][3], __m256d o[3], __m256d di[3], int s[3], int f[3], __m256d *t) {
	//printf("align %d\n", (int)(__alignof__(box)));
	//printf("ad %d\n", &box[0][0]);
	__m256d tmin = _mm256_set1_pd(-INFINITY);
	__m256d tmax = _mm256_set1_pd(INFINITY);
	for (int i = 0; i < 3; ++i) {
		if (f[i]) {
			tmin = _mm256_max_pd(tmin, _mm256_mul_pd(di[i], _mm256_sub_pd(box[s[i]][i], o[i])));
			tmax = _mm256_min_pd(tmax, _mm256_mul_pd(di[i], _mm256_sub_pd(box[1 - s[i]][i], o[i])));
			//puts("PASSED");
		}
	}
	*t = _mm256_max_pd(_mm256_set1_pd(0.0), tmin);
	return _mm256_movemask_pd(_mm256_and_pd(_mm256_cmp_pd(tmin, tmax, _CMP_LE_OQ), _mm256_cmp_pd(tmax, _mm256_set1_pd(0.0), _CMP_GT_OQ)));
}

void bvh_simd_traverse(int *fs, double *vs, BVH_NODE_SIMD1 node, int *stack, int *ptr, __m256d o[3], __m256d d[3], __m256d di[3], int s[3], int f[3], double *t, int *id) {
	//printf("id: %d\n", *ptr);
	//printf("child id: %d %d %d %d\n", node.c[0], node.c[1], node.c[2], node.c[3]);
	//printf("nb nt ti: %d %d %d\n", node.nb, node.nt, node.it);
	for (int i = 0; i < node.nt; ++i) {
		__m256d vst[3][3];
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				vst[j][k] = _mm256_set_pd(
				vs[3 * fs[3 * (node.it + 4 * i + 3) + j] + k],
				vs[3 * fs[3 * (node.it + 4 * i + 2) + j] + k],
				vs[3 * fs[3 * (node.it + 4 * i + 1) + j] + k],
				vs[3 * fs[3 * (node.it + 4 * i + 0) + j] + k]);
			}
		}
		__m256d td;
		int hit = tri_isect_simd(vst[0], vst[1], vst[2], o, d, &td);
		for (int j = 0; j < 4; ++j) {
			if (((hit >> j) & 1) && (td[j] < *t)) {
				*t = td[j];
				*id = node.it + 4 * i + j;
			}
		}
	}
	__m256d b[2][3];
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 3; ++j) {
			b[i][j] = _mm256_loadu_pd(node.b[i][j]);
		}
	}
	__m256d td;
	int hit = aabb_isect_simd(b, o, di, s, f, &td);
	for (int i = 0; i < node.nb; ++i) {
		if (((hit >> i) & 1) && (td[i] < *t)) {
			PUSH(stack, *ptr, node.c[i]);
		}
	}
}

int bvh_simd_traverse_wrap(int *fs, double *vs, BVH_NODE_SIMD1 *nodes, int *stack, VEC od, VEC dd) {
	BVH_NODE_SIMD1 node = nodes[0];
	int ptr = 0;
	__m256d o[3] = {_mm256_set1_pd(od.x), _mm256_set1_pd(od.y), _mm256_set1_pd(od.z)};
	__m256d d[3] = {_mm256_set1_pd(dd.x), _mm256_set1_pd(dd.y), _mm256_set1_pd(dd.z)};
	__m256d di[3] = {_mm256_set1_pd(1.0 / dd.x), _mm256_set1_pd(1.0 / dd.y), _mm256_set1_pd(1.0 / dd.z)};
	int s[3] = {dd.x > 0 ? 0 : 1, dd.y > 0 ? 0 : 1, dd.z > 0 ? 0 : 1};
	int f[3] = {fabs(dd.x) < EPS ? 0 : 1, fabs(dd.y) < EPS ? 0 : 1, fabs(dd.z) < EPS ? 0 : 1};
	double t = INFINITY;
	int id = -1;
	while (true) {
		bvh_simd_traverse(fs, vs, node, stack, &ptr, o, d, di, s, f, &t, &id);
		if (ptr == 0) {
			break;
		}
		int c = POP(stack, ptr);
		node = nodes[c];
	}
	return id;
}

int main(int argc, char **argv) {
	char *fn_scene = "dragon.txt";
	
	int n = 1, w = 512 * 1, h = 512 * 1;
	double fovh = 60.0;
	double e = 1.0;
	omp_set_num_threads(NUM_THREADS);
	
	double fx = tan(fovh * PI / 360.0);
	double fy = h * fx / w;
	
	int nv, nf;
	double *vs;
	int *fs;
	read_scene(fn_scene, &nv, &vs, &nf, &fs);
	
	double *img = ALLOC(double, w * h);
	
	double t;
	
	t = time_get();
	
	int bvh_simd_num_alloc = 1000000;
	BVH_NODE_SIMD1 *nodes_simd = ALLOC(BVH_NODE_SIMD1, bvh_simd_num_alloc);
	int bvh_simd_num = bvh_simd_build_wrap1(nf, &fs, vs, bvh_simd_num_alloc, nodes_simd);
	
	t = time_get() - t;
	
	printf("bvh simd build: %f\n", t);
	
	printf("bvh simd usage: %d/%d\n", bvh_simd_num, bvh_simd_num_alloc);
	
	t = time_get();
	
	{
		for (int i = 0; i < w * h; ++i) {
			img[i] = 0.0;
		}
		for (int i = 0; i < bvh_simd_num; ++i) {
			BVH_NODE_SIMD1 node = nodes_simd[i];
			for (int j = 0; j < node.nb; ++j) {
				AABB box;
				box.min.x = node.b[0][0][j];
				box.min.y = node.b[0][1][j];
				box.min.z = node.b[0][2][j];
				box.max.x = node.b[1][0][j];
				box.max.y = node.b[1][1][j];
				box.max.z = node.b[1][2][j];
				VEC p0[12], p1[12];
				p0[0] = vec_init(box.min.x, box.min.y, box.min.z);
				p1[0] = vec_init(box.max.x, box.min.y, box.min.z);
				p0[1] = vec_init(box.min.x, box.min.y, box.min.z);
				p1[1] = vec_init(box.min.x, box.max.y, box.min.z);
				p0[2] = vec_init(box.min.x, box.min.y, box.min.z);
				p1[2] = vec_init(box.min.x, box.min.y, box.max.z);
				p0[3] = vec_init(box.max.x, box.min.y, box.max.z);
				p1[3] = vec_init(box.max.x, box.max.y, box.max.z);
				p0[4] = vec_init(box.max.x, box.min.y, box.max.z);
				p1[4] = vec_init(box.max.x, box.min.y, box.min.z);
				p0[5] = vec_init(box.max.x, box.min.y, box.max.z);
				p1[5] = vec_init(box.min.x, box.min.y, box.max.z);
				p0[6] = vec_init(box.max.x, box.max.y, box.min.z);
				p1[6] = vec_init(box.min.x, box.max.y, box.min.z);
				p0[7] = vec_init(box.max.x, box.max.y, box.min.z);
				p1[7] = vec_init(box.max.x, box.min.y, box.min.z);
				p0[8] = vec_init(box.max.x, box.max.y, box.min.z);
				p1[8] = vec_init(box.max.x, box.max.y, box.max.z);
				p0[9] = vec_init(box.min.x, box.max.y, box.max.z);
				p1[9] = vec_init(box.min.x, box.min.y, box.max.z);
				p0[10] = vec_init(box.min.x, box.max.y, box.max.z);
				p1[10] = vec_init(box.min.x, box.max.y, box.min.z);
				p0[11] = vec_init(box.min.x, box.max.y, box.max.z);
				p1[11] = vec_init(box.max.x, box.max.y, box.max.z);
				for (int j = 0; j < 12; ++j) {
					VEC p;
					VEC a = p0[j];
					p = vec_init(a.x / (fx * a.z), a.y / (fy * a.z), a.z);
					if (p.z > 0.0 || p.x < -1.0 || p.x > 1.0 || p.y < -1.0 || p.y > 1.0) {
						continue;
					}
					a = vec_scale(-1.0 / a.z, a);
					a = vec_init(0.5 * w * (a.x / fx + 1.0), 0.5 * h * (a.y / fy + 1.0), a.z);
					VEC b = p1[j];
					p = vec_init(b.x / (fx * b.z), b.y / (fy * b.z), b.z);
					if (p.z > 0.0 || p.x < -1.0 || p.x > 1.0 || p.y < -1.0 || p.y > 1.0) {
						continue;
					}
					b = vec_scale(-1.0 / b.z, b);
					b = vec_init(0.5 * w * (b.x / fx + 1.0), 0.5 * h * (b.y / fy + 1.0), b.z);
					a.x = fmax(0, fmin(w, a.x));
					a.y = fmax(0, fmin(h, a.y));
					b.x = fmax(0, fmin(w, b.x));
					b.y = fmax(0, fmin(h, b.y));
					VEC d = vec_sub(b, a);
					int px = fabs(d.x) > fabs(d.y) ? 1 : 0;
					int py = fabs(d.x) > fabs(d.y) ? 0 : 1;
					int sx = d.x > 0.0 ? 1 : -1;
					int sy = d.y > 0.0 ? 1 : -1;
					int ix = min(w - 1, a.x);
					int iy = min(h - 1, a.y);
					while (true) {
						double tmin = -INFINITY, tmax = INFINITY;
						if (fabs(d.x) > EPS) {
							double tx0 = (ix - a.x) / d.x;
							double tx1 = (ix + 1 - a.x) / d.x;
							tmin = fmax(tmin, fmin(tx0, tx1));
							tmax = fmin(tmax, fmax(tx0, tx1));
							if (tmax < tmin) {
								ix = ix + sx * (py - px);
								iy = iy + sy * (px - py);
								continue;
							}
						}
						if (fabs(d.y) > EPS) {
							double ty0 = (iy - a.y) / d.y;
							double ty1 = (iy + 1 - a.y) / d.y;
							tmin = fmax(tmin, fmin(ty0, ty1));
							tmax = fmin(tmax, fmax(ty0, ty1));
							if (tmax < tmin) {
								ix = ix + sx * (py - px);
								iy = iy + sy * (px - py);
								continue;
							}
						}
						tmin = fmax(0.0, tmin);
						tmax = fmin(1.0, tmax);
						img[ix + w * iy] += e * (tmax - tmin) * vec_len(d);
						if (1.0 - tmax < EPS) {
							break;
						}
						ix += sx * px;
						iy += sy * py;
					}
				}
			}
		}
	}
	
	t = time_get() - t;
	
	printf("bvhsimd: %f\n", t);
	
	img_save("bvhsimd.pfm", w, h, img);
	
	t = time_get();
	
	int bvh_num_alloc = 1000000;
	BVH_NODE *nodes = ALLOC(BVH_NODE, bvh_num_alloc);
	int bvh_num = sah_build_bvh(nf, &fs, vs, nodes);
	
	t = time_get() - t;
	
	printf("bvh build: %f\n", t);
	
	printf("bvh usage: %d/%d\n", bvh_num, bvh_num_alloc);
	
	t = time_get();
	
	{
		for (int i = 0; i < w * h; ++i) {
			img[i] = 0.0;
		}
		for (int i = 0; i < bvh_num; ++i) {
			AABB box = nodes[i].b;
			VEC p0[12], p1[12];
			p0[0] = vec_init(box.min.x, box.min.y, box.min.z);
			p1[0] = vec_init(box.max.x, box.min.y, box.min.z);
			p0[1] = vec_init(box.min.x, box.min.y, box.min.z);
			p1[1] = vec_init(box.min.x, box.max.y, box.min.z);
			p0[2] = vec_init(box.min.x, box.min.y, box.min.z);
			p1[2] = vec_init(box.min.x, box.min.y, box.max.z);
			p0[3] = vec_init(box.max.x, box.min.y, box.max.z);
			p1[3] = vec_init(box.max.x, box.max.y, box.max.z);
			p0[4] = vec_init(box.max.x, box.min.y, box.max.z);
			p1[4] = vec_init(box.max.x, box.min.y, box.min.z);
			p0[5] = vec_init(box.max.x, box.min.y, box.max.z);
			p1[5] = vec_init(box.min.x, box.min.y, box.max.z);
			p0[6] = vec_init(box.max.x, box.max.y, box.min.z);
			p1[6] = vec_init(box.min.x, box.max.y, box.min.z);
			p0[7] = vec_init(box.max.x, box.max.y, box.min.z);
			p1[7] = vec_init(box.max.x, box.min.y, box.min.z);
			p0[8] = vec_init(box.max.x, box.max.y, box.min.z);
			p1[8] = vec_init(box.max.x, box.max.y, box.max.z);
			p0[9] = vec_init(box.min.x, box.max.y, box.max.z);
			p1[9] = vec_init(box.min.x, box.min.y, box.max.z);
			p0[10] = vec_init(box.min.x, box.max.y, box.max.z);
			p1[10] = vec_init(box.min.x, box.max.y, box.min.z);
			p0[11] = vec_init(box.min.x, box.max.y, box.max.z);
			p1[11] = vec_init(box.max.x, box.max.y, box.max.z);
			for (int j = 0; j < 12; ++j) {
				VEC p;
				VEC a = p0[j];
				p = vec_init(a.x / (fx * a.z), a.y / (fy * a.z), a.z);
				if (p.z > 0.0 || p.x < -1.0 || p.x > 1.0 || p.y < -1.0 || p.y > 1.0) {
					continue;
				}
				a = vec_scale(-1.0 / a.z, a);
				a = vec_init(0.5 * w * (a.x / fx + 1.0), 0.5 * h * (a.y / fy + 1.0), a.z);
				VEC b = p1[j];
				p = vec_init(b.x / (fx * b.z), b.y / (fy * b.z), b.z);
				if (p.z > 0.0 || p.x < -1.0 || p.x > 1.0 || p.y < -1.0 || p.y > 1.0) {
					continue;
				}
				b = vec_scale(-1.0 / b.z, b);
				b = vec_init(0.5 * w * (b.x / fx + 1.0), 0.5 * h * (b.y / fy + 1.0), b.z);
				a.x = fmax(0, fmin(w, a.x));
				a.y = fmax(0, fmin(h, a.y));
				b.x = fmax(0, fmin(w, b.x));
				b.y = fmax(0, fmin(h, b.y));
				VEC d = vec_sub(b, a);
				int px = fabs(d.x) > fabs(d.y) ? 1 : 0;
				int py = fabs(d.x) > fabs(d.y) ? 0 : 1;
				int sx = d.x > 0.0 ? 1 : -1;
				int sy = d.y > 0.0 ? 1 : -1;
				int ix = min(w - 1, a.x);
				int iy = min(h - 1, a.y);
				while (true) {
					double tmin = -INFINITY, tmax = INFINITY;
					if (fabs(d.x) > EPS) {
						double tx0 = (ix - a.x) / d.x;
						double tx1 = (ix + 1 - a.x) / d.x;
						tmin = fmax(tmin, fmin(tx0, tx1));
						tmax = fmin(tmax, fmax(tx0, tx1));
						if (tmax < tmin) {
							ix = ix + sx * (py - px);
							iy = iy + sy * (px - py);
							continue;
						}
					}
					if (fabs(d.y) > EPS) {
						double ty0 = (iy - a.y) / d.y;
						double ty1 = (iy + 1 - a.y) / d.y;
						tmin = fmax(tmin, fmin(ty0, ty1));
						tmax = fmin(tmax, fmax(ty0, ty1));
						if (tmax < tmin) {
							ix = ix + sx * (py - px);
							iy = iy + sy * (px - py);
							continue;
						}
					}
					tmin = fmax(0.0, tmin);
					tmax = fmin(1.0, tmax);
					img[ix + w * iy] += e * (tmax - tmin) * vec_len(d);
					if (1.0 - tmax < EPS) {
						break;
					}
					ix += sx * px;
					iy += sy * py;
				}
			}
		}
	}
	
	t = time_get() - t;
	
	printf("bvh: %f\n", t);
	
	img_save("bvh.pfm", w, h, img);
	
	t = time_get();
	
	{
		for (int i = 0; i < w * h; ++i) {
			img[i] = 0.0;
		}
		for (int i = 0; i < nf; ++i) {
			for (int j = 0; j < 3; ++j) {
				VEC p;
				VEC a = vec_init(vs[3 * fs[3 * i + (j + 0) % 3] + 0], vs[3 * fs[3 * i + (j + 0) % 3] + 1], vs[3 * fs[3 * i + (j + 0) % 3] + 2]);
				p = vec_init(a.x / (fx * a.z), a.y / (fy * a.z), a.z);
				if (p.z > 0.0 || p.x < -1.0 || p.x > 1.0 || p.y < -1.0 || p.y > 1.0) {
					continue;
				}
				a = vec_scale(-1.0 / a.z, a);
				a = vec_init(0.5 * w * (a.x / fx + 1.0), 0.5 * h * (a.y / fy + 1.0), a.z);
				VEC b = vec_init(vs[3 * fs[3 * i + (j + 1) % 3] + 0], vs[3 * fs[3 * i + (j + 1) % 3] + 1], vs[3 * fs[3 * i + (j + 1) % 3] + 2]);
				p = vec_init(b.x / (fx * b.z), b.y / (fy * b.z), b.z);
				if (p.z > 0.0 || p.x < -1.0 || p.x > 1.0 || p.y < -1.0 || p.y > 1.0) {
					continue;
				}
				b = vec_scale(-1.0 / b.z, b);
				b = vec_init(0.5 * w * (b.x / fx + 1.0), 0.5 * h * (b.y / fy + 1.0), b.z);
				a.x = fmax(0, fmin(w, a.x));
				a.y = fmax(0, fmin(h, a.y));
				b.x = fmax(0, fmin(w, b.x));
				b.y = fmax(0, fmin(h, b.y));
				VEC d = vec_sub(b, a);
				int px = fabs(d.x) > fabs(d.y) ? 1 : 0;
				int py = fabs(d.x) > fabs(d.y) ? 0 : 1;
				int sx = d.x > 0.0 ? 1 : -1;
				int sy = d.y > 0.0 ? 1 : -1;
				int ix = min(w - 1, a.x);
				int iy = min(h - 1, a.y);
				while (true) {
					double tmin = -INFINITY, tmax = INFINITY;
					if (fabs(d.x) > EPS) {
						double tx0 = (ix - a.x) / d.x;
						double tx1 = (ix + 1 - a.x) / d.x;
						tmin = fmax(tmin, fmin(tx0, tx1));
						tmax = fmin(tmax, fmax(tx0, tx1));
						if (tmax < tmin) {
							ix = ix + sx * (py - px);
							iy = iy + sy * (px - py);
							continue;
						}
					}
					if (fabs(d.y) > EPS) {
						double ty0 = (iy - a.y) / d.y;
						double ty1 = (iy + 1 - a.y) / d.y;
						tmin = fmax(tmin, fmin(ty0, ty1));
						tmax = fmin(tmax, fmax(ty0, ty1));
						if (tmax < tmin) {
							ix = ix + sx * (py - px);
							iy = iy + sy * (px - py);
							continue;
						}
					}
					tmin = fmax(0.0, tmin);
					tmax = fmin(1.0, tmax);
					img[ix + w * iy] += e * (tmax - tmin) * vec_len(d);
					if (1.0 - tmax < EPS) {
						break;
					}
					ix += sx * px;
					iy += sy * py;
				}
			}
		}
	}
	
	t = time_get() - t;
	
	printf("mesh: %f\n", t);
	
	img_save("mesh.pfm", w, h, img);
	
	t = time_get();
	
	#pragma omp parallel
	{
		for (int i = 0; i < w * h; ++i) {
			img[i] = 0.0;
		}
		
		PRNG_STATE s = prng_seed(SEED);
		for (int i = 0; i < omp_get_thread_num(); ++i) {
			s = prng_jump(s);
		}
		
		for (int i = omp_get_thread_num(); i < w * h; i += omp_get_max_threads()) {
			int ix = 0, iy = 0;
			ix = DecodeMorton2X(i);
			iy = DecodeMorton2Y(i);
			double sum = 0.0;
			for (int j = 0; j < n; ++j) {
				VEC o = vec_init(0.0, 0.0, 0.0);
				VEC d = vec_norm(vec_init(fx * (2.0 * (prng_db(&s) + ix) / w - 1.0), fy * (2.0 * (prng_db(&s) + iy) / h - 1.0), -1.0));
				
				int id = bvh_traverse(fs, vs, nodes, o, d);
				if (id >= 0) {
					VEC v0, v1, v2;
					face_ld_vs(fs, vs, id, &v0, &v1, &v2);
					
					VEC normal = tri_norm(v0, v1, v2);
					
					sum += e * fabs(vec_dot(normal, d)) / PI;
				}
			}
			img[ix + w * iy] = sum / n;
		}
	}
	
	t = time_get() - t;
	
	printf("nosimd: %f\n", t);
	
	img_save("nosimd.pfm", w, h, img);
	
	t = time_get();
	
	#pragma omp parallel
	{
		PRNG_STATE s = prng_seed(SEED);
		for (int i = 0; i < omp_get_thread_num(); ++i) {
			s = prng_jump(s);
		}
		
		int *stack = ALLOC(int, bvh_simd_num_alloc);
		
		for (int i = omp_get_thread_num(); i < w * h; i += omp_get_max_threads()) {
			int ix = 0, iy = 0;
			ix = DecodeMorton2X(i);
			iy = DecodeMorton2Y(i);
			double sum = 0.0;
			for (int j = 0; j < n; ++j) {
				VEC o = vec_init(0.0, 0.0, 0.0);
				VEC d = vec_norm(vec_init(fx * (2.0 * (prng_db(&s) + ix) / w - 1.0), fy * (2.0 * (prng_db(&s) + iy) / h - 1.0), -1.0));
				
				int id = bvh_simd_traverse_wrap(fs, vs, nodes_simd, stack, o, d);
				//int id = tri_isect_simd_wrap(vs, nf, fs, o, d);
				if (id >= 0) {
					VEC v0 = vec_init(vs[3 * fs[3 * id + 0] + 0], vs[3 * fs[3 * id + 0] + 1], vs[3 * fs[3 * id + 0] + 2]);
					VEC v1 = vec_init(vs[3 * fs[3 * id + 1] + 0], vs[3 * fs[3 * id + 1] + 1], vs[3 * fs[3 * id + 1] + 2]);
					VEC v2 = vec_init(vs[3 * fs[3 * id + 2] + 0], vs[3 * fs[3 * id + 2] + 1], vs[3 * fs[3 * id + 2] + 2]);
					VEC normal = tri_norm(v0, v1, v2);
					
					sum += e * fmax(0.0, -vec_dot(normal, d)) / PI;
				}
			}
			img[ix + w * iy] = sum / n;
		}
	}
	
	t = time_get() - t;
	
	printf("simd: %f\n", t);
	
	img_save("simd.pfm", w, h, img);
	
	puts("\a");
	
	return 0;
}