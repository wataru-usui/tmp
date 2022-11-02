#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define CAST(t, x) (t)(x)
#define ALLOC(t, n) malloc(sizeof(t) * (n))
#define PI 3.141592653589793
#define EPS 1.0e-6
#define SEED 0
#define DEST_PATH "o.pfm"
#define BETA 2.0

typedef enum {
  EDF_UNI, EDF_COS, EDF_DELTA
} EDF;
typedef enum {
  BSDF_NULL, BSDF_LAMBERT, BSDF_CONDUCTOR, BSDF_DIELECTRIC
} BSDF;
typedef struct {
  uint64_t s0, s1;
} PRNG_STATE;
typedef struct {
  double x, y, z;
} DOUBLE3;
typedef struct {
  int n, w, h;
  double *data;
} IMAGE;
typedef struct {
  double n;
  DOUBLE3 a;
} MEDIUM;
typedef struct {
  DOUBLE3 e, s;
  EDF edf;
  BSDF bsdf;
  MEDIUM mp, mn;
} MATERIAL;
typedef struct {
  double r;
  DOUBLE3 c;
} SPHERE;
typedef struct {
  MATERIAL m;
  SPHERE s;
} OBJECT;
//wi is not used for y0 and z0
typedef struct {
  int id;
  bool c;
  double pf, pb;
  DOUBLE3 a, x, n, wi;
} VERTEX;
typedef struct {
  int n;
  VERTEX *vs;
} PATH;
typedef struct {
  PATH *paths;
  int i, j;
} NODE;
typedef struct {
  int num_threads, num_samples, num_objs, num_verts, num_phots, num_nodes;
  int wi, hi;
  double ef, wf, hf, df;
  double r, f;
  OBJECT *objs;
  double *pdf, *cdf;
  int num_iters;
  double init_r;
} SCENE;

double sign(double x) {
  return x < 0.0 ? -1.0 : x > 0.0 ? 1.0 : 0.0;
}
double sqr(double x) {
  return x * x;
}
double cube(double x) {
  return x * x * x;
}
double opst(double x) {
  return sqrt(fmax(0.0, 1.0 - sqr(x)));
}

double time_get() {
  return clock() / CAST(double, CLOCKS_PER_SEC);
}

uint64_t sm64_next(uint64_t *s) {
  uint64_t x = *s += 0x9E3779B97F4A7C15;
  x = 0xBF58476D1CE4E5B9 * (x ^ (x >> 30));
  x = 0x94D049BB133111EB * (x ^ (x >> 27));
  x = x ^ (x >> 31);
  return x;
}
uint64_t rotl(uint64_t a, uint64_t b) {
  return (a << b) | (a >> (64 - b));
}
uint64_t prng_next_u64(PRNG_STATE *s) {
  PRNG_STATE a = *s;
  uint64_t ret = a.s0 + a.s1;
  a.s1 ^= a.s0;
  s->s0 = rotl(a.s0, 55) ^ a.s1 ^ (a.s1 << 14);
  s->s1 = rotl(a.s1, 36);
  return ret;
}
double prng_next_double(PRNG_STATE *s) {
  return 0x1.0p-53 * (prng_next_u64(s) >> 11);
}
PRNG_STATE prng_seed(uint64_t seed) {
  PRNG_STATE s = {sm64_next(&seed), sm64_next(&seed)};
  return s;
}
PRNG_STATE prng_jump(PRNG_STATE s) {
  PRNG_STATE jump = {0xBEAC0467EBA5FACB, 0xD86B048B86AA9922};
  PRNG_STATE next = {0, 0};
  for (int i = 0; i < 64; ++i) {
    if (jump.s0 & 1 << i) {
      next.s0 ^= s.s0;
      next.s1 ^= s.s1;
    }
    prng_next_u64(&s);
  }
  for (int i = 0; i < 64; ++i) {
    if (jump.s1 & 1 << i) {
      next.s0 ^= s.s0;
      next.s1 ^= s.s1;
    }
    prng_next_u64(&s);
  }
  return next;
}

DOUBLE3 double3_init3(double x, double y, double z) {
  DOUBLE3 v = {x, y, z};
  return v;
}
DOUBLE3 double3_init(double a) {
  return double3_init3(a, a, a);
}
DOUBLE3 double3_add(DOUBLE3 a, DOUBLE3 b) {
  return double3_init3(a.x + b.x, a.y + b.y, a.z + b.z);
}
DOUBLE3 double3_sub(DOUBLE3 a, DOUBLE3 b) {
  return double3_init3(a.x - b.x, a.y - b.y, a.z - b.z);
}
DOUBLE3 double3_mul(DOUBLE3 a, DOUBLE3 b) {
  return double3_init3(a.x * b.x, a.y * b.y, a.z * b.z);
}
DOUBLE3 double3_div(DOUBLE3 a, DOUBLE3 b) {
  return double3_init3(a.x / b.x, a.y / b.y, a.z / b.z);
}
DOUBLE3 double3_scale(double a, DOUBLE3 b) {
  return double3_init3(a * b.x, a * b.y, a * b.z);
}
DOUBLE3 double3_cross(DOUBLE3 a, DOUBLE3 b) {
  return double3_init3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
double double3_dot(DOUBLE3 a, DOUBLE3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
double double3_len(DOUBLE3 a) {
  return sqrt(double3_dot(a, a));
}
DOUBLE3 double3_norm(DOUBLE3 a) {
  return double3_scale(1.0 / double3_len(a), a);
}
DOUBLE3 double3_trafo(DOUBLE3 x, DOUBLE3 y, DOUBLE3 z, DOUBLE3 a) {
  return double3_add(double3_add(
  double3_scale(a.x, x), double3_scale(a.y, y)), double3_scale(a.z, z));
}
DOUBLE3 double3_pow(DOUBLE3 a, DOUBLE3 b) {
  return double3_init3(pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z));
}
DOUBLE3 double3_exp(DOUBLE3 a) {
  return double3_init3(exp(a.x), exp(a.y), exp(a.z));
}
DOUBLE3 double3_min(DOUBLE3 a, DOUBLE3 b) {
  return double3_init3(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z));
}
DOUBLE3 double3_max(DOUBLE3 a, DOUBLE3 b) {
  return double3_init3(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z));
}
bool double3_equal(DOUBLE3 a, DOUBLE3 b) {
  double eps = EPS;
  DOUBLE3 d = double3_sub(a, b);
  return fabs(d.x) <= eps && fabs(d.y) <= eps && fabs(d.z) <= eps;
}
void double3_print(DOUBLE3 a, FILE *o) {
  fprintf(o, "%f %f %f\n", a.x, a.y, a.z);
}
void orthogonal_basis(DOUBLE3 *x, DOUBLE3 *y, DOUBLE3 z) {
  if (fabs(z.x) > fabs(z.y)) {
    *x = double3_scale(1.0 / sqrt(sqr(z.x) + sqr(z.z)), double3_init3(-z.z, 0.0, z.x));
  } else {
    *x = double3_scale(1.0 / sqrt(sqr(z.y) + sqr(z.z)), double3_init3(0.0, z.z, -z.y));
  }
  *y = double3_cross(z, *x);
}
DOUBLE3 spherical_coord(double cost, double sint, double cosp, double sinp) {
  return double3_init3(sint * cosp, sint * sinp, cost);
}
DOUBLE3 refl(DOUBLE3 n, DOUBLE3 d, double cos) {
  return double3_sub(d, double3_scale(2.0 * cos, n));
}
DOUBLE3 refr(DOUBLE3 n, DOUBLE3 d, double eta, double cosi, double cost) {
  return double3_sub(double3_scale(eta, d), double3_scale(eta * cosi - cost, n));
}

double luminance(DOUBLE3 a) {
  return double3_dot(double3_init3(0.2126, 0.7152, 0.0722), a);
}
DOUBLE3 tonemap_linear(DOUBLE3 a) {
  return double3_max(double3_init(0.0), double3_min(double3_init(1.0), a));
}
DOUBLE3 tonemap_exp(DOUBLE3 a) {
  return double3_sub(double3_init(1.0), double3_exp(double3_mul(double3_init3(-1.0, -1.0, -1.0), a)));
}
DOUBLE3 tonemap_filmic(DOUBLE3 a) {
  a = double3_max(double3_init(0.0), double3_sub(a, double3_init(0.004)));
  DOUBLE3 m, n;
  m = double3_mul(a, double3_add(double3_init(0.5),
  double3_mul(double3_init(6.2), a)));
  n = double3_add(double3_init(0.06),
  double3_mul(a, double3_add(double3_init(1.7),
  double3_mul(double3_init(6.2), a))));
  return double3_div(m, n);
}

IMAGE image_init(int n, int w, int h) {
  IMAGE img = {n, w, h, ALLOC(double, n * w * h)};
  return img;
}
void image_del(IMAGE img) {
  free(img.data);
}
void image_fill(IMAGE img, double x) {
  for (int i = 0; i < img.n * img.w * img.h; ++i) {
    img.data[i] = x;
  }
}
bool image_save_ppm(IMAGE img, char *path) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) {
    return false;
  }
  fprintf(f, "P6 %d %d %d ", img.w, img.h, 255);
  char *data = ALLOC(char, 3 * img.w * img.h);
  for (int i = 0; i < 3 * img.w * img.h; ++i) {
    data[i] = floor(255.0 * img.data[i] + 0.5);
  }
  fwrite(data, 1, 3 * img.w * img.h, f);
  free(data);
  return true;
}
bool image_save_pfm(IMAGE img, char *path) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) {
    return false;
  }
  fprintf(f, "PF %d %d -1 ", img.w, img.h);
  float *data = ALLOC(float, 3 * img.w * img.h);
  for (int iy = 0; iy < img.h; ++iy) {
    for (int ix = 0; ix < img.w; ++ix) {
      data[3 * (ix + img.w * iy) + 0] = CAST(float, img.data[3 * (ix + img.w * (img.h - iy - 1)) + 0]);
      data[3 * (ix + img.w * iy) + 1] = CAST(float, img.data[3 * (ix + img.w * (img.h - iy - 1)) + 1]);
      data[3 * (ix + img.w * iy) + 2] = CAST(float, img.data[3 * (ix + img.w * (img.h - iy - 1)) + 2]);
    }
  }
  fwrite(data, sizeof(float), 3 * img.w * img.h, f);
  free(data);
  return true;
}

double fresnel_dielectric(double eta, double cosi, double cost) {
  double perp = (eta * cosi - cost) / (eta * cosi + cost);
  double paral = (eta * cost - cosi) / (eta * cost + cosi);
  return (sqr(perp) + sqr(paral)) / 2.0;
}
double fresnel_conductor_approx(double eta, double k, double cosi) {
  cosi = fabs(cosi);
  double a = sqr(eta) + sqr(k);
  double perp = (a * sqr(cosi) - 2.0 * eta * cosi + 1.0) / (a * sqr(cosi) + 2.0 * eta * cosi + 1.0);
  double paral = (a - 2.0 * eta * cosi + sqr(cosi)) / (a + 2.0 * eta * cosi + sqr(cosi));
  return (perp + paral) / 2.0;
}
double fresnel_conductor(double eta, double k, double cosi) {
  cosi = fabs(cosi);
  double s = 1.0 - sqr(cosi);
  double i = sqr(eta) - sqr(k) - s;
  double w = sqrt(sqr(i) + sqr(2.0 * eta * k));
  double a = sqrt(fmax(0.0, (w + i) / 2.0));
  double perp = (w + sqr(cosi) - 2.0 * a * cosi) / (w + sqr(cosi) + 2.0 * a * cosi);
  double paral = (sqr(cosi) * w + sqr(s) - 2.0 * a * cosi * s) / (sqr(cosi) * w + sqr(s) + 2.0 * a * cosi * s);
  return (perp + perp * paral) / 2.0;
}

double hemisphere_cosine_weighted_pdf(DOUBLE3 n, DOUBLE3 wo) {
  double cos = double3_dot(n, wo);
  return cos > 0.0 ? cos / PI : 0.0;
}
DOUBLE3 hemisphere_cosine_weighted_sample(DOUBLE3 n, PRNG_STATE *state) {
  double u0, u1;
  u0 = prng_next_double(state);
  u1 = prng_next_double(state);
  double phi = 2.0 * PI * u1;
  double cost, sint, cosp, sinp;
  cost = sqrt(u0);
  sint = opst(cost);
  cosp = cos(phi);
  sinp = sin(phi);
  DOUBLE3 x, y;
  orthogonal_basis(&x, &y, n);
  return double3_trafo(x, y, n, spherical_coord(cost, sint, cosp, sinp));
}
double cos_pow_pdf(DOUBLE3 n, DOUBLE3 wo) {
  double a = 256.0;
  double cos = double3_dot(n, wo);
  return cos > 0.0 ? (2.0 + a) * pow(cos, 1.0 + a) / (2.0 * PI) : 0.0;
}
DOUBLE3 cos_pow_sample(DOUBLE3 n, PRNG_STATE *state) {
  double a = 256.0;
  double u0, u1;
  u0 = prng_next_double(state);
  u1 = prng_next_double(state);
  double phi = 2.0 * PI * u1;
  double cost, sint, cosp, sinp;
  cost = pow(1.0 - u0, 1.0 / (2.0 + a));
  sint = opst(cost);
  cosp = cos(phi);
  sinp = sin(phi);
  DOUBLE3 x, y;
  orthogonal_basis(&x, &y, n);
  return double3_trafo(x, y, n, spherical_coord(cost, sint, cosp, sinp));
}
double sphere_area(SPHERE s) {
  return 4.0 * PI * sqr(s.r);
}
double sphere_uniform_pdf(SPHERE s) {
  return 1.0 / sphere_area(s);
}
DOUBLE3 sphere_uniform_sample(SPHERE s, PRNG_STATE *state) {
  double u0, u1;
  u0 = prng_next_double(state);
  u1 = prng_next_double(state);
  double phi = 2.0 * PI * u1;
  double cost, sint, cosp, sinp;
  cost = 2.0 * u0 - 1.0;
  sint = opst(cost);
  cosp = cos(phi);
  sinp = sin(phi);
  return double3_add(s.c, double3_scale(s.r, spherical_coord(cost, sint, cosp, sinp)));
}
int discrete_sample(int n, double *cdf, PRNG_STATE *state) {
  double x = prng_next_double(state);
  for(int i = 0; i < n; ++i) {
    if (x < cdf[i]) {
      return i;
    }
  }
  return -1;
}

DOUBLE3 edf_eval(EDF edf, DOUBLE3 n, DOUBLE3 wo) {
  switch (edf) {
    case EDF_UNI : {
      double cos = double3_dot(n, wo);
      return double3_init(cos > 0.0 ? 1.0 / PI : 0.0);
    }
    case EDF_COS : {
      double a = 256.0;
      double cos = double3_dot(n, wo);
      return double3_init(cos > 0.0 ? (2.0 + a) * pow(cos, 1.0 + a) / (2.0 * PI) : 0.0);
    }
    case EDF_DELTA : {
      return double3_init(double3_equal(wo, n) ? 1.0 : 0.0);
    }
    default : {
      return double3_init(0.0);
    }
  }
}
double edf_pdf(EDF edf, DOUBLE3 n, DOUBLE3 wo) {
  switch (edf) {
    case EDF_UNI : {
      return hemisphere_cosine_weighted_pdf(n, wo);
    }
    case EDF_COS : {
      return cos_pow_pdf(n, wo);
    }
    case EDF_DELTA : {
      return double3_equal(wo, n) ? 1.0 : 0.0;
    }
    default : {
      return 0.0;
    }
  }
}
DOUBLE3 edf_sample(EDF edf, DOUBLE3 n, PRNG_STATE *state) {
  switch (edf) {
    case EDF_UNI : {
      return hemisphere_cosine_weighted_sample(n, state);
    }
    case EDF_COS : {
      return cos_pow_sample(n, state);
    }
    case EDF_DELTA : {
      return n;
    }
    default : {
      return double3_init(0.0);
    }
  }
}
bool edf_delta(EDF edf) {
  switch (edf) {
    case EDF_UNI : {
      return false;
    }
    case EDF_COS : {
      return false;
    }
    case EDF_DELTA : {
      return true;
    }
    default : {
      return false;
    }
  }
}

DOUBLE3 bsdf_eval(BSDF bsdf, double ni, double nt, DOUBLE3 n, DOUBLE3 wi, DOUBLE3 wo) {
  switch (bsdf) {
    case BSDF_NULL : {
      double cosi = double3_dot(n, wi);
      return double3_init(double3_equal(wi, wo) ? 1.0 / fabs(cosi) : 0.0);
    }
    case BSDF_LAMBERT : {
      double cosi = double3_dot(n, wi), coso = double3_dot(n, wo);
      return double3_init(cosi * coso < 0.0 ? 1.0 / PI : 0.0);
    }
    case BSDF_CONDUCTOR : {
      double m = 2.7732, k = 2.9278;
      double cosi = double3_dot(n, wi);
      DOUBLE3 r = refl(n, wi, cosi);
      return double3_init(double3_equal(wo, r) ? fresnel_conductor(ni / m, k, cosi) / fabs(cosi) : 0.0);
    }
    case BSDF_DIELECTRIC : {
      double cosi = double3_dot(n, wi);
      double eta = ni / nt;
      double eta2 = sqr(eta);
      double det = 1.0 - eta2 * (1.0 - sqr(cosi));
      if (det >= 0.0) {
        double cost = sign(cosi) * sqrt(det);
        double fr = fresnel_dielectric(eta, cosi, cost);
        DOUBLE3 r = refl(n, wi, cosi);
        if (double3_equal(r, wo)) {
          return double3_init(fr / fabs(cosi));
        }
        DOUBLE3 t = refr(n, wi, eta, cosi, cost);
        if (double3_equal(t, wo)) {
          return double3_init((1.0 - fr) / (eta2 * fabs(cosi)));
        }
      }
      DOUBLE3 r = refl(n, wi, cosi);
      return double3_init(double3_equal(r, wo) ? 1.0 / fabs(cosi) : 0.0);
    }
    default : {
      double cosi = double3_dot(n, wi);
      return double3_init(double3_equal(wi, wo) ? 1.0 / fabs(cosi) : 0.0);
    }
  }
}
double bsdf_pdf(BSDF bsdf, double ni, double nt, DOUBLE3 n, DOUBLE3 wi, DOUBLE3 wo) {
  switch (bsdf) {
    case BSDF_NULL : {
      return double3_equal(wi, wo) ? 1.0 : 0.0;
    }
    case BSDF_LAMBERT : {
      double cosi = double3_dot(n, wi);
      return hemisphere_cosine_weighted_pdf(double3_scale(-sign(cosi), n), wo);
    }
    case BSDF_CONDUCTOR : {
      double cosi = double3_dot(n, wi);
      DOUBLE3 r = refl(n, wi, cosi);
      return double3_equal(wo, r) ? 1.0 : 0.0;
    }
    case BSDF_DIELECTRIC : {
      double cosi = double3_dot(n, wi);
      double eta = ni / nt;
      double eta2 = sqr(eta);
      double det = 1.0 - eta2 * (1.0 - sqr(cosi));
      if (det >= 0.0) {
        double cost = sign(cosi) * sqrt(det);
        double fr = fresnel_dielectric(eta, cosi, cost);
        DOUBLE3 r = refl(n, wi, cosi);
        if (double3_equal(r, wo)) {
          return fr;
        }
        DOUBLE3 t = refr(n, wi, eta, cosi, cost);
        if (double3_equal(t, wo)) {
          return 1.0 - fr;
        }
      }
      DOUBLE3 r = refl(n, wi, cosi);
      return double3_equal(r, wo) ? 1.0 : 0.0;
    }
    default : {
      return double3_equal(wi, wo) ? 1.0 : 0.0;
    }
  }
}
DOUBLE3 bsdf_sample(BSDF bsdf, double ni, double nt, DOUBLE3 n, DOUBLE3 wi, PRNG_STATE *state) {
  switch (bsdf) {
    case BSDF_NULL : {
      return wi;
    }
    case BSDF_LAMBERT : {
      double cosi = double3_dot(n, wi);
      return
      hemisphere_cosine_weighted_sample(double3_scale(-sign(cosi), n), state);
    }
    case BSDF_CONDUCTOR : {
      double cosi = double3_dot(n, wi);
      DOUBLE3 r = refl(n, wi, cosi);
      return r;
    }
    case BSDF_DIELECTRIC : {
      double cosi = double3_dot(n, wi);
      double eta = ni / nt;
      double eta2 = sqr(eta);
      double det = 1.0 - eta2 * (1.0 - sqr(cosi));
      if (det >= 0.0) {
        double cost = sign(cosi) * sqrt(det);
        double fr = fresnel_dielectric(eta, cosi, cost);
        if (prng_next_double(state) >= fr) {
          return refr(n, wi, eta, cosi, cost);
        }
      }
      return refl(n, wi, cosi);
    }
    default : {
      return wi;
    }
  }
}
bool bsdf_delta(BSDF bsdf) {
  switch (bsdf) {
    case BSDF_NULL : {
      return true;
    }
    case BSDF_LAMBERT : {
      return false;
    }
    case BSDF_CONDUCTOR : {
      return true;
    }
    case BSDF_DIELECTRIC : {
      return true;
    }
    default : {
      return false;
    }
  }
}

MEDIUM medium_init(double n, DOUBLE3 a) {
  MEDIUM m = {n, a};
  return m;
}
DOUBLE3 transmittance(DOUBLE3 tr, double t) {
  return double3_exp(double3_scale(-t, tr));
}
void media_sides(double cos, MEDIUM mp, MEDIUM mn, MEDIUM *mr, MEDIUM *mt) {
  bool into = cos < 0.0;
  *mr = into ? mp : mn;
  *mt = into ? mn : mp;
}

MATERIAL material_init(DOUBLE3 e, DOUBLE3 s, EDF edf, BSDF bsdf, MEDIUM mp, MEDIUM mn) {
  MATERIAL m = {e, s, edf, bsdf, mp, mn};
  return m;
}

SPHERE sphere_init(double r, DOUBLE3 c) {
  SPHERE s = {r, c};
  return s;
}
DOUBLE3 sphere_normal(SPHERE s, DOUBLE3 x) {
  return double3_scale(1.0 / s.r, double3_sub(x, s.c));
}
bool sphere_intersect(SPHERE s, DOUBLE3 o, DOUBLE3 d, double *t) {
  double eps = EPS;
  DOUBLE3 co = double3_sub(o, s.c);
  double a = double3_dot(d, co);
  double det = sqr(s.r) + sqr(a) - double3_dot(co, co);
  if (det < 0.0) {
    return false;
  }
  det = sqrt(det);
  return (*t = -a - det) > eps ? true : (*t = -a + det) > eps ? true : false;
}

OBJECT object_init(MATERIAL m, SPHERE s) {
  OBJECT o = {m, s};
  return o;
}

bool intersect(int num_objs, OBJECT *objs, DOUBLE3 o, DOUBLE3 d, int *id, double *t) {
  *t = INFINITY;
  for (int i = 0; i < num_objs; ++i) {
    OBJECT obj_cur = objs[i];
    double t_cur;
    if (sphere_intersect(obj_cur.s, o, d, &t_cur) && t_cur < *t) {
      *id = i;
      *t = t_cur;
    }
  }
  return *t < INFINITY;
}
bool visible(int num_objs, OBJECT *objs, DOUBLE3 o, DOUBLE3 d, double t) {
  double eps = EPS;
  for (int i = 0; i < num_objs; ++i) {
    OBJECT obj_cur = objs[i];
    double t_cur;
    if (sphere_intersect(obj_cur.s, o, d, &t_cur) && t_cur + eps < t) {
      return false;
    }
  }
  return true;
}

double circ_area(double r) {
  return PI * sqr(r);
}
int cmp_x(void const *a, void const *b) {
  NODE const *A = a, *B = b;
  PATH *paths = A->paths;
  double d =  paths[A->i].vs[A->j].x.x - paths[B->i].vs[B->j].x.x;
  return d > 0.0 ? 1 : d < 0.0 ? -1 : 0;
}
int cmp_y(void const *a, void const *b) {
  NODE const *A = a, *B = b;
  PATH *paths = A->paths;
  double d =  paths[A->i].vs[A->j].x.y - paths[B->i].vs[B->j].x.y;
  return d > 0.0 ? 1 : d < 0.0 ? -1 : 0;
}
int cmp_z(void const *a, void const *b) {
  NODE const *A = a, *B = b;
  PATH *paths = A->paths;
  double d =  paths[A->i].vs[A->j].x.z - paths[B->i].vs[B->j].x.z;
  return d > 0.0 ? 1 : d < 0.0 ? -1 : 0;
}
void kdtree_construct(int a, int b, NODE *nodes, int depth) {
  if (b - a <= 1) {
    return;
  }
  int axis = depth % 3;
  switch (axis) {
    case 0 : {
      qsort(nodes + a, b - a, sizeof(NODE), cmp_x);
    } break;
    case 1 : {
      qsort(nodes + a, b - a, sizeof(NODE), cmp_y);
    } break;
    case 2 : {
      qsort(nodes + a, b - a, sizeof(NODE), cmp_z);
    } break;
    default : {
    } break;
  }
  int im = (a + b) / 2;
  kdtree_construct(a, im, nodes, depth + 1);
  kdtree_construct(im + 1, b, nodes, depth + 1);
}
void kdtree_gather(int a, int b, NODE *nodes, int depth,
OBJECT *objs, int N, double r, PATH pe, VERTEX ve, int je, DOUBLE3 *sum) {
  if (b - a <= 0) {
    return;
  }
  int im = (a + b) / 2;
  NODE node = nodes[im];
  int il = node.i, jl = node.j;
  PATH pl = node.paths[il];
  VERTEX vl = pl.vs[jl];
  DOUBLE3 d = double3_sub(vl.x, ve.x);
  double dist = double3_len(d);
  if (dist < r) {
    OBJECT obj = objs[ve.id];
    
    VERTEX vlm1 = pl.vs[jl - 1];
    VERTEX vem1 = pe.vs[je - 1];
    
    double tl = double3_len(double3_sub(vl.x, vlm1.x));
    double te = double3_len(double3_sub(ve.x, vem1.x));
    
    double gl = fabs(double3_dot(vlm1.n, vl.wi)) / sqr(tl);
    double ge = fabs(double3_dot(vem1.n, ve.wi)) / sqr(te);
    
    MEDIUM mlr, mlt, mer, met;
    media_sides(double3_dot(vl.n, vl.wi), obj.m.mp, obj.m.mn, &mlr, &mlt);
    media_sides(double3_dot(ve.n, ve.wi), obj.m.mp, obj.m.mn, &mer, &met);
    
    double psl = bsdf_pdf(obj.m.bsdf, mlr.n, mlt.n, vl.n, vl.wi, double3_scale(-1.0, ve.wi));
    double pse = bsdf_pdf(obj.m.bsdf, mer.n, met.n, ve.n, ve.wi, double3_scale(-1.0, vl.wi));
    
    double pblr = luminance(double3_mul(transmittance(mer.a, te), obj.m.s));
    double pbla = gl * pse;
    double pbl = pblr * pbla;
    
    double pber = luminance(double3_mul(transmittance(mlr.a, tl), obj.m.s));
    double pbea = ge * psl;
    double pbe = pber * pbea;
    
    double w = 1.0;
    
    if (jl >= 2) {
      double p = 1.0;
      p *= pbl / vl.pf;
      if (vlm1.c) {
        w += pow(p, BETA);
      }
      for (int kl = jl - 2; kl >= 1; --kl) {
        VERTEX vlp1 = pl.vs[kl + 1], vl = pl.vs[kl];
        p *= vl.pb / vlp1.pf;
        if (vl.c) {
          w += pow(p, BETA);
        }
      }
    }
    
    if (jl == 1) {
      double p = 1.0 / (N * circ_area(r));
      p *= 1.0 / vl.pf;
      if (vlm1.c) {
        w += pow(p, BETA);
      }
      p *= pbl / vlm1.pf;
      if (vlm1.c) {
        w += pow(p, BETA);
      }
    } else if (jl == 2) {
      VERTEX vlm2 = pl.vs[0];
      double p = 1.0 / (N * circ_area(r));
      p *= 1.0 / vl.pf;
      if (vlm1.c) {
        w += pow(p, BETA);
      }
      p *= pbl / vlm1.pf;
      if (vlm1.c && vlm2.c) {
        w += pow(p, BETA);
      }
      p *= vlm2.pb / vlm2.pf;
      if (vlm2.c) {
        w += pow(p, BETA);
      }
    } else {
      VERTEX vlm2 = pl.vs[jl - 2];
      double p = 1.0 / (N * circ_area(r));
      p *= 1.0 / vl.pf;
      if (vlm1.c) {
        w += pow(p, BETA);
      }
      p *= pbl / vlm1.pf;
      if (vlm1.c && vlm2.c) {
        w += pow(p, BETA);
      }
      for (int kl = jl - 2; kl >= 2; --kl) {
        VERTEX vl = pl.vs[kl], vlm1 = pl.vs[kl - 1];
        p *= vl.pb / vl.pf;
        if (vl.c && vlm1.c) {
          w += pow(p, BETA);
        }
      }
      VERTEX vl1 = pl.vs[1], vl0 = pl.vs[0];
      p *= vl1.pb / vl1.pf;
      if (vl1.c && vl0.c) {
        w += pow(p, BETA);
      }
      p *= vl0.pb / vl0.pf;
      if (vl0.c) {
        w += pow(p, BETA);
      }
    }
    
    if (je >= 2) {
      double p = 1.0;
      p *= pbe / ve.pf;
      if (vem1.c) {
        w += pow(p, BETA);
      }
      for (int ke = je - 2; ke >= 1; --ke) {
        VERTEX vep1 = pe.vs[ke + 1], ve = pe.vs[ke];
        p *= ve.pb / vep1.pf;
        if (ve.c) {
          w += pow(p, BETA);
        }
      }
    }
    
    if (je == 1) {
    } else if (je == 2) {
      double p = 1.0 / (N * circ_area(r));
      p *= 1.0 / ve.pf;
      if (vem1.c) {
        w += pow(p, BETA);
      }
    } else {
      VERTEX vem2 = pe.vs[je - 2];
      double p = 1.0 / (N * circ_area(r));
      p *= 1.0 / ve.pf;
      if (vem1.c) {
        w += pow(p, BETA);
      }
      p *= pbe / vem1.pf;
      if (vem1.c && vem2.c) {
        w += pow(p, BETA);
      }
      for (int ke = je - 2; ke >= 2; --ke) {
        VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
        p *= ve.pb / ve.pf;
        if (ve.c && vem1.c) {
          w += pow(p, BETA);
        }
      }
    }
    
    w = 1.0 / w;
    
    double k = 1.0 / circ_area(r);
    DOUBLE3 fs = double3_mul(obj.m.s, bsdf_eval(obj.m.bsdf, mlr.n, mlt.n, ve.n, vl.wi, double3_scale(-1.0, ve.wi)));
    DOUBLE3 I = double3_scale(w * k, double3_mul(double3_mul(ve.a, vl.a), fs));
    *sum = double3_add(*sum, I);
  }
  double t;
  int axis = depth % 3;
  switch (axis) {
    case 0 : {
      t = ve.x.x - vl.x.x;
    } break;
    case 1 : {
      t = ve.x.y - vl.x.y;
    } break;
    case 2 : {
      t = ve.x.z - vl.x.z;
    } break;
    default : {
      t = 0.0;
    } break;
  }
  if (t < r) {
    kdtree_gather(a, im, nodes, depth + 1, objs, N, r, pe, ve, je, sum);
  }
  if (t > -r) {
    kdtree_gather(im + 1, b, nodes, depth + 1, objs, N, r, pe, ve, je, sum);
  }
}
DOUBLE3 vm(SCENE scene, int N, double r, int num_nodes, NODE *nodes, PATH pe) {
  DOUBLE3 sum0 = double3_init(0.0);
  for (int je = 1; je < pe.n; ++je) {
    VERTEX ve = pe.vs[je];
    OBJECT obje = scene.objs[ve.id];
    if (bsdf_delta(obje.m.bsdf)) {
      continue;
    }
    DOUBLE3 sum1 = double3_init(0.0);
    kdtree_gather(0, num_nodes, nodes, 0, scene.objs, N, r, pe, ve, je, &sum1);
    sum0 = double3_add(sum0, sum1);
  }
  sum0 = double3_scale(1.0 / N, sum0);
  return sum0;
}
DOUBLE3 vc(SCENE scene, int N, double r, PATH pl, PATH pe) {
  DOUBLE3 sum0 = double3_init(0.0);
  //s = 0, t >= 2
  for (int je = 1; je < pe.n; ++je) {
    VERTEX ve = pe.vs[je];
    OBJECT obje = scene.objs[ve.id];
    bool ce = !edf_delta(obje.m.edf);
    if (ce) {
      VERTEX vem1 = pe.vs[je - 1];
      
      double te = double3_len(double3_sub(ve.x, vem1.x));
      
      double ge1 = fabs(double3_dot(vem1.n, ve.wi)) / sqr(te);
      
      double pse = edf_pdf(obje.m.edf, ve.n, double3_scale(-1.0, ve.wi));
      
      double pbea0 = scene.pdf[ve.id] * sphere_uniform_pdf(obje.s);
      double pbe0 = pbea0;
      
      double pbea1 = ge1 * pse;
      double pbe1 = pbea1;
      
      double w = 1.0;
      
      if (je == 1) {
      } else if (je == 2) {
        double p = 1.0;
        p *= pbe0 / ve.pf;
        if (vem1.c) {
          w += pow(p, BETA);
        }
      } else {
        VERTEX vem2 = pe.vs[je - 2];
        double p = 1.0;
        p *= pbe0 / ve.pf;
        if (vem1.c) {
          w += pow(p, BETA);
        }
        p *= pbe1 / vem1.pf;
        if (vem1.c && vem2.c) {
          w += pow(p, BETA);
        }
        for (int ke = je - 2; ke >= 2; --ke) {
          VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
          p *= ve.pb / ve.pf;
          if (ve.c && vem1.c) {
            w += pow(p, BETA);
          }
        }
      }
      
      if (je == 1) {
      } else {
        double p = N * circ_area(r);
        p *= pbe0;
        p *= pbe1 / ve.pf;
        if (vem1.c) {
          w += pow(p, BETA);
        }
        for (int ke = je - 2; ke >= 1; --ke) {
          VERTEX vep1 = pe.vs[ke + 1], ve = pe.vs[ke];
          p *= ve.pb / vep1.pf;
          if (ve.c) {
            w += pow(p, BETA);
          }
        }
      }
      
      w = 1.0 / w;
      
      DOUBLE3 fsl = double3_mul(obje.m.e, edf_eval(obje.m.edf, ve.n, double3_scale(-1.0, ve.wi)));
      DOUBLE3 I = double3_scale(w, double3_mul(ve.a, fsl));
      sum0 = double3_add(sum0, I);
    }
  }
  //s = 1, t >= 2
  {
    VERTEX vl = pl.vs[0];
    OBJECT objl = scene.objs[vl.id];
    bool cl = vl.c;
    if (cl) {
      for (int je = 1; je < pe.n; ++je) {
        VERTEX ve = pe.vs[je];
        OBJECT obje = scene.objs[ve.id];
        bool ce = ve.c;
        if (ce) {
          DOUBLE3 d = double3_sub(vl.x, ve.x);
          double t = double3_len(d);
          DOUBLE3 dl = double3_scale(-1.0 / t, d);
          DOUBLE3 de = double3_scale(1.0 / t, d);
          if (visible(scene.num_objs, scene.objs, ve.x, de, t)) {
            VERTEX vem1 = pe.vs[je - 1];
            
            double te = double3_len(double3_sub(ve.x, vem1.x));
            
            double g = fabs(double3_dot(vl.n, dl) * double3_dot(ve.n, de)) / sqr(t);
            double gl0 = fabs(double3_dot(vl.n, dl)) / sqr(t);
            double ge0 = fabs(double3_dot(ve.n, de)) / sqr(t);
            double ge1 = fabs(double3_dot(vem1.n, ve.wi)) / sqr(te);
            
            MEDIUM meir, meit, meor, meot;
            media_sides(double3_dot(ve.n, ve.wi), obje.m.mp, obje.m.mn, &meir, &meit);
            media_sides(double3_dot(ve.n, de), obje.m.mp, obje.m.mn, &meor, &meot);
            
            DOUBLE3 tr = transmittance(meot.a, t);
            
            double psl = edf_pdf(objl.m.edf, vl.n, dl);
            double pse = bsdf_pdf(obje.m.bsdf, meir.n, meit.n, ve.n, ve.wi, de);
            
            double pblr0 = luminance(double3_mul(transmittance(meir.a, te), obje.m.s));
            double pbla0 = gl0 * pse;
            double pbl0 = pblr0 * pbla0;
            
            double pbea0 = ge0 * psl;
            double pbe0 = pbea0;
            
            double pber1 = luminance(double3_mul(tr, obje.m.s));
            double pbea1 = ge1 * pse;
            double pbe1 = pber1 * pbea1;
            
            double w = 1.0;
            
            {
              double p = 1.0;
              p *= pbl0 / vl.pf;
              w += pow(p, BETA);
            }
            
            if (je == 1) {
            } else if (je == 2) {
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += pow(p, BETA);
              }
            } else {
              VERTEX vem2 = pe.vs[je - 2];
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += pow(p, BETA);
              }
              p *= pbe1 / vem1.pf;
              if (vem1.c && vem2.c) {
                w += pow(p, BETA);
              }
              for (int ke = je - 2; ke >= 2; --ke) {
                VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
                p *= ve.pb / ve.pf;
                if (ve.c && vem1.c) {
                  w += pow(p, BETA);
                }
              }
            }
            
            if (je == 1) {
              double p = N * circ_area(r);
              p *= pbe0;
              w += pow(p, BETA);
            } else {
              double p = N * circ_area(r);
              p *= pbe0;
              w += pow(p, BETA);
              p *= pbe1 / ve.pf;
              if (vem1.c) {
                w += pow(p, BETA);
              }
              for (int ke = je - 2; ke >= 1; --ke) {
                VERTEX vep1 = pe.vs[ke + 1], ve = pe.vs[ke];
                p *= ve.pb / vep1.pf;
                if (ve.c) {
                  w += pow(p, BETA);
                }
              }
            }
            
            w = 1.0 / w;
            
            DOUBLE3 fsl = double3_mul(objl.m.e, edf_eval(objl.m.edf, vl.n, dl));
            DOUBLE3 fse = double3_mul(obje.m.s, bsdf_eval(obje.m.bsdf, meot.n, meor.n, ve.n, dl, double3_scale(-1.0, ve.wi)));
            DOUBLE3 I = double3_scale(w * g, double3_mul(double3_mul(double3_mul(double3_mul(vl.a, ve.a), fsl), fse), tr));
            sum0 = double3_add(sum0, I);
          }
        }
      }
    }
  }
  //s >= 2, t >= 2
  for (int jl = 1; jl < pl.n; ++jl) {
    VERTEX vl = pl.vs[jl];
    OBJECT objl = scene.objs[vl.id];
    bool cl = vl.c;
    if (cl) {
      for (int je = 1; je < pe.n; ++je) {
        VERTEX ve = pe.vs[je];
        OBJECT obje = scene.objs[ve.id];
        bool ce = ve.c;
        if (ce) {
          DOUBLE3 d = double3_sub(vl.x, ve.x);
          double t = double3_len(d);
          DOUBLE3 dl = double3_scale(-1.0 / t, d);
          DOUBLE3 de = double3_scale(1.0 / t, d);
          if (visible(scene.num_objs, scene.objs, ve.x, de, t)) {
            VERTEX vlm1 = pl.vs[jl - 1];
            VERTEX vem1 = pe.vs[je - 1];
            
            double tl = double3_len(double3_sub(vl.x, vlm1.x));
            double te = double3_len(double3_sub(ve.x, vem1.x));
            
            double g = fabs(double3_dot(vl.n, dl) * double3_dot(ve.n, de)) / sqr(t);
            double gl0 = fabs(double3_dot(vl.n, dl)) / sqr(t);
            double gl1 = fabs(double3_dot(vlm1.n, vl.wi)) / sqr(tl);
            double ge0 = fabs(double3_dot(ve.n, de)) / sqr(t);
            double ge1 = fabs(double3_dot(vem1.n, ve.wi)) / sqr(te);
            
            MEDIUM mlir, mlit, mlor, mlot, meir, meit, meor, meot;
            media_sides(double3_dot(vl.n, vl.wi), objl.m.mp, objl.m.mn, &mlir, &mlit);
            media_sides(double3_dot(vl.n, dl), objl.m.mp, objl.m.mn, &mlor, &mlot);
            media_sides(double3_dot(ve.n, ve.wi), obje.m.mp, obje.m.mn, &meir, &meit);
            media_sides(double3_dot(ve.n, de), obje.m.mp, obje.m.mn, &meor, &meot);
            
            DOUBLE3 tr = transmittance(meot.a, t);
            
            double psl = bsdf_pdf(objl.m.bsdf, mlir.n, mlit.n, vl.n, vl.wi, dl);
            double pse = bsdf_pdf(obje.m.bsdf, meir.n, meit.n, ve.n, ve.wi, de);
            
            double pblr0 = luminance(double3_mul(transmittance(meir.a, te), obje.m.s));
            double pbla0 = gl0 * pse;
            double pbl0 = pblr0 * pbla0;
            
            double pblr1 = luminance(double3_mul(tr, objl.m.s));
            double pbla1 = gl1 * psl;
            double pbl1 = pblr1 * pbla1;
            
            double pber0 = luminance(double3_mul(transmittance(mlir.a, tl), objl.m.s));
            double pbea0 = ge0 * psl;
            double pbe0 = pber0 * pbea0;
            
            double pber1 = luminance(double3_mul(tr, obje.m.s));
            double pbea1 = ge1 * pse;
            double pbe1 = pber1 * pbea1;
            
            double w = 1.0;
            
            if (jl == 1) {
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += pow(p, BETA);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c) {
                w += pow(p, BETA);
              }
            } else if (jl == 2) {
              VERTEX vlm2 = pl.vs[jl - 2];
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += pow(p, BETA);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c && vlm2.c) {
                w += pow(p, BETA);
              }
              p *= vlm2.pb / vlm2.pf;
              if (vlm2.c) {
                w += pow(p, BETA);
              }
            } else {
              VERTEX vlm2 = pl.vs[jl - 2];
              double p = 1.0;
              p *= pbl0 / vl.pf;
              if (vlm1.c) {
                w += pow(p, BETA);
              }
              p *= pbl1 / vlm1.pf;
              if (vlm1.c && vlm2.c) {
                w += pow(p, BETA);
              }
              for (int kl = jl - 2; kl >= 2; --kl) {
                VERTEX vl = pl.vs[kl], vlm1 = pl.vs[kl - 1];
                p *= vl.pb / vl.pf;
                if (vl.c && vlm1.c) {
                  w += pow(p, BETA);
                }
              }
              VERTEX vl1 = pl.vs[1], vl0 = pl.vs[0];
              p *= vl1.pb / vl1.pf;
              if (vl1.c && vl0.c) {
                w += pow(p, BETA);
              }
              p *= vl0.pb / vl0.pf;
              if (vl0.c) {
                w += pow(p, BETA);
              }
            }
            
            if (jl == 1) {
              double p = N * circ_area(r);
              p *= pbl0;
              w += pow(p, BETA);
            } else {
              double p = N * circ_area(r);
              p *= pbl0;
              w += pow(p, BETA);
              p *= pbl1 / vl.pf;
              if (vlm1.c) {
                w += pow(p, BETA);
              }
              for (int kl = jl - 2; kl >= 1; --kl) {
                VERTEX vlp1 = pl.vs[kl + 1], vl = pl.vs[kl];
                p *= vl.pb / vlp1.pf;
                if (vl.c) {
                  w += pow(p, BETA);
                }
              }
            }
            
            if (je == 1) {
            } else if (je == 2) {
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += pow(p, BETA);
              }
            } else {
              VERTEX vem2 = pe.vs[je - 2];
              double p = 1.0;
              p *= pbe0 / ve.pf;
              if (vem1.c) {
                w += pow(p, BETA);
              }
              p *= pbe1 / vem1.pf;
              if (vem1.c && vem2.c) {
                w += pow(p, BETA);
              }
              for (int ke = je - 2; ke >= 2; --ke) {
                VERTEX ve = pe.vs[ke], vem1 = pe.vs[ke - 1];
                p *= ve.pb / ve.pf;
                if (ve.c && vem1.c) {
                  w += pow(p, BETA);
                }
              }
            }
            
            if (je == 1) {
              double p = N * circ_area(r);
              p *= pbe0;
              w += pow(p, BETA);
            } else {
              double p = N * circ_area(r);
              p *= pbe0;
              w += pow(p, BETA);
              p *= pbe1 / ve.pf;
              if (vem1.c) {
                w += pow(p, BETA);
              }
              for (int ke = je - 2; ke >= 1; --ke) {
                VERTEX vep1 = pe.vs[ke + 1], ve = pe.vs[ke];
                p *= ve.pb / vep1.pf;
                if (ve.c) {
                  w += pow(p, BETA);
                }
              }
            }
            
            w = 1.0 / w;
            
            DOUBLE3 fsl = double3_mul(objl.m.s, bsdf_eval(objl.m.bsdf, mlir.n, mlit.n, vl.n, vl.wi, dl));
            DOUBLE3 fse = double3_mul(obje.m.s, bsdf_eval(obje.m.bsdf, meot.n, meor.n, ve.n, dl, double3_scale(-1.0, ve.wi)));
            DOUBLE3 I = double3_scale(w * g, double3_mul(double3_mul(double3_mul(double3_mul(vl.a, ve.a), fsl), fse), tr));
            sum0 = double3_add(sum0, I);
          }
        }
      }
    }
  }
  return sum0;
}
void trace(SCENE scene, PRNG_STATE *s, PATH *path, int id, bool c, double pa, double ps,
DOUBLE3 n, DOUBLE3 o, DOUBLE3 d, DOUBLE3 fs, bool light) {
  DOUBLE3 a = double3_init(1.0 / pa);
  path->vs[0].id = id;
  path->vs[0].c = c;
  path->vs[0].pf = pa;
  path->vs[0].a = a;
  path->vs[0].x = o;
  path->vs[0].n = n;
  a = double3_scale(fabs(double3_dot(n, d)) / ps, double3_mul(a, fs));
  double prf = 1.0;
  for (int i = 1; i < scene.num_verts; ++i) {
    if (prng_next_double(s) >= prf) {
      path->n = i;
      return;
    }
    int id;
    double t;
    if (!intersect(scene.num_objs, scene.objs, o, d, &id, &t)) {
      path->n = i;
      return;
    }
    OBJECT obj = scene.objs[id];
    o = double3_add(o, double3_scale(t, d));
    n = sphere_normal(obj.s, o);
    
    double cosi = double3_dot(n, d);
    MEDIUM mir, mit;
    media_sides(cosi, obj.m.mp, obj.m.mn, &mir, &mit);
    
    DOUBLE3 wo = bsdf_sample(obj.m.bsdf, mir.n, mit.n, n, d, s);
    
    double coso = double3_dot(n, wo);
    MEDIUM mor, mot;
    media_sides(coso, obj.m.mp, obj.m.mn, &mor, &mot);
    
    DOUBLE3 tr = transmittance(mir.a, t);
    
    a = double3_scale(1.0 / prf, double3_mul(tr, a));
    
    VERTEX v;
    v.id = id;
    v.c = !bsdf_delta(obj.m.bsdf);
    v.pf = prf * ps * fabs(cosi) / sqr(t);
    if (i >= 2) {
      VERTEX v_m1 = path->vs[i - 1], *v_m2 = path->vs + i - 2;
      OBJECT obj_m1 = scene.objs[v_m1.id];
      double prb_m2 = luminance(double3_mul(tr, obj_m1.m.s));
      double coso_m2 = double3_dot(v_m2->n, v_m1.wi);
      double t_m1 = double3_len(double3_sub(v_m1.x, v_m2->x));
      v_m2->pb = prb_m2 * ps * fabs(coso_m2) / sqr(t_m1);
    }
    v.a = a;
    v.x = o;
    v.n = n;
    v.wi = d;
    path->vs[i] = v;
    
    fs = light ?
    bsdf_eval(obj.m.bsdf, mir.n, mit.n, n,
    double3_scale(1.0, d), double3_scale(1.0, wo)) :
    bsdf_eval(obj.m.bsdf, mot.n, mor.n, n,
    double3_scale(-1.0, wo), double3_scale(-1.0, d));
    
    ps = bsdf_pdf(obj.m.bsdf, mir.n, mit.n, n, d, wo);
    
    a = double3_scale(fabs(coso) / ps, double3_mul(double3_mul(a, obj.m.s), fs));
    
    prf = luminance(double3_mul(tr, obj.m.s));
    d = wo;
  }
  path->n = scene.num_verts;;
}
void sample_eye(SCENE scene, PRNG_STATE *s, int ix, int iy,
DOUBLE3 *n, DOUBLE3 *o, DOUBLE3 *d, DOUBLE3 *We, double *pa, double *ps, bool *c) {
  *n = double3_init3(0.0, 0.0, 1.0);
  double r = scene.r * sqrt(prng_next_double(s));
  double t = 2.0 * PI * prng_next_double(s);
  *o = double3_scale(r, double3_init3(cos(t), sin(t), 0.0));
  *d = double3_init3(
  scene.wf * (0.5 - (ix + prng_next_double(s)) / scene.wi),
  scene.hf * ((iy + prng_next_double(s)) / scene.hi - 0.5),
  scene.df);
  *d = double3_norm(double3_sub(double3_scale(
  scene.f / (scene.f - scene.df), *d), *o));
  *We = double3_init(
  scene.ef *
  (1.0 / (PI * sqr(scene.r))) *
  (sqr(scene.df) / (cube(fabs(d->z)) * scene.wf * scene.hf)));
  *pa = 1.0 / (PI * sqr(scene.r));
  *ps = sqr(scene.df) / (cube(fabs(d->z)) * scene.wf * scene.hf);
  *c = true;
}
void sample_light(SCENE scene, PRNG_STATE *s,
DOUBLE3 *n, DOUBLE3 *o, DOUBLE3 *d, DOUBLE3 *Le, double *pa, double *ps, int *id, bool *c) {
  *id = discrete_sample(scene.num_objs, scene.cdf, s);
  OBJECT obj = scene.objs[*id];
  *o = sphere_uniform_sample(obj.s, s);
  *n = sphere_normal(obj.s, *o);
  *d = edf_sample(obj.m.edf, *n, s);
  *Le = double3_mul(obj.m.e, edf_eval(obj.m.edf, *n, *d));
  *pa = scene.pdf[*id] * sphere_uniform_pdf(obj.s);
  *ps = edf_pdf(obj.m.edf, *n, *d);
  *c = !edf_delta(obj.m.edf);
}
void pt(SCENE scene, IMAGE img, PRNG_STATE *s) {
  PATH path;
  path.vs = ALLOC(VERTEX, scene.num_verts);
  for (int i = omp_get_thread_num(); i < img.w * img.h;
  i += omp_get_max_threads()) {
    int ix = i % img.w, iy = i / img.w;
    DOUBLE3 sum = double3_init(0.0);
    for (int i = 0; i < scene.num_samples; ++i) {
      DOUBLE3 n, o, d, We;
      double pa, ps;
      bool c;
      sample_eye(scene, s, ix, iy, &n, &o, &d, &We, &pa, &ps, &c);
      trace(scene, s, &path, -1, c, pa, ps, n, o, d, We, false);
      DOUBLE3 sum1 = double3_init(0.0);
      for (int i = 1; i < path.n; ++i) {
        VERTEX v = path.vs[i];
        OBJECT obj = scene.objs[v.id];
        DOUBLE3 Le = double3_mul(obj.m.e, edf_eval(obj.m.edf, v.n, double3_scale(-1.0, v.wi)));
        sum1 = double3_add(sum1, double3_mul(v.a, Le));
      }
      sum = double3_add(sum, double3_scale(1.0 / scene.num_samples, sum1));
    }
    sum = double3_scale(1.0 / (img.w * img.h), sum);
    img.data[3 * i + 0] = sum.x;
    img.data[3 * i + 1] = sum.y;
    img.data[3 * i + 2] = sum.z;
  }
}
IMAGE render_pt(SCENE scene) {
  IMAGE img = image_init(3, scene.wi, scene.hi);
  omp_set_num_threads(scene.num_threads);
  #pragma omp parallel
  {
    PRNG_STATE s = prng_seed(SEED);
    for (int i = 0; i < omp_get_thread_num(); ++i) {
      s = prng_jump(s);
    }
    pt(scene, img, &s);
  }
  return img;
}
void vcm(SCENE scene, IMAGE img, PRNG_STATE *states) {
  image_fill(img, 0.0);
  int N = scene.num_samples * scene.wi * scene.hi;
  VERTEX *vss[omp_get_max_threads()];
  int num_phots_per_thread = scene.num_phots / omp_get_max_threads();
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    vss[i] = ALLOC(VERTEX, num_phots_per_thread);
  }
  PATH *pls = ALLOC(PATH, N);
  PATH pes[omp_get_max_threads()];
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    pes[i].vs = ALLOC(VERTEX, scene.num_verts);
  }
  NODE *nodes = ALLOC(NODE, scene.num_nodes);
  for (int iter_num = 0; iter_num < scene.num_iters; ++iter_num) {
    double t;
    printf("number of iterations: %d\n", iter_num);
    double r = scene.init_r * pow(iter_num + 1, ((2.0 * BETA + 1.0) / (2.0 * BETA + 2.0) - 1.0) / 2.0);
    t = time_get();
    #pragma omp parallel
    {
      PRNG_STATE *s = states + omp_get_thread_num();
      int offset = 0;
      for (int i = omp_get_thread_num(); i < N; i += omp_get_max_threads()) {
        PATH *pl = pls + i;
        pl->vs = vss[omp_get_thread_num()] + offset;
        DOUBLE3 n, o, d, Le;
        double pa, ps;
        int id;
        bool c;
        sample_light(scene, s, &n, &o, &d, &Le, &pa, &ps, &id, &c);
        trace(scene, s, pl, id, c, pa, ps, n, o, d, Le, true);
        offset += pl->n;
      }
    }
    t = time_get() - t;
    printf("light paths trace time: %f\n", t);
    int num_phots = 0;
    for (int i = 0; i < N; ++i) {
      num_phots += pls[i].n;
    }
    printf("number of photons: %d/%d\n", num_phots, scene.num_phots);
    int num_nodes = 0;
    for (int i = 0; i < N; ++i) {
      PATH path = pls[i];
      for (int j = 1; j < path.n; ++j) {
        if (!bsdf_delta(scene.objs[pls[i].vs[j].id].m.bsdf)) {
          NODE node = {pls, i, j};
          nodes[num_nodes] = node;
          ++num_nodes;
        }
      }
    }
    printf("number of nodes: %d/%d\n", num_nodes, scene.num_nodes);
    printf("average number of verticies: %f\n", num_phots / CAST(double, N));
    t = time_get();
    kdtree_construct(0, num_nodes, nodes, 0);
    t = time_get() - t;
    printf("kd tree build time: %f\n", t);
    t = time_get();
    #pragma omp parallel
    {
      PRNG_STATE *s = states + omp_get_thread_num();
      for (int i = omp_get_thread_num(); i < img.w * img.h;
      i += omp_get_max_threads()) {
        int ix = i % img.w, iy = i / img.w;
        DOUBLE3 sum0 = double3_init(0.0);
        for (int j = 0; j < scene.num_samples; ++j) {
          PATH *pl = pls + scene.num_samples * i + j;
          PATH *pe = pes + omp_get_thread_num();
          DOUBLE3 n, o, d, We;
          double pa, ps;
          bool c;
          sample_eye(scene, s, ix, iy, &n, &o, &d, &We, &pa, &ps, &c);
          trace(scene, s, pe, -1, c, pa, ps, n, o, d, We, false);
          sum0 = double3_add(sum0, vm(scene, N, r, num_nodes, nodes, *pe));
          sum0 = double3_add(sum0, vc(scene, N, r, *pl, *pe));
        }
        sum0 = double3_scale(1.0 / N, sum0);
        img.data[3 * i + 0] += sum0.x;
        img.data[3 * i + 1] += sum0.y;
        img.data[3 * i + 2] += sum0.z;
      }
    }
    t = time_get() - t;
    printf("evaluation time: %f\n", t);
  }
  for (int i = 0; i < img.n * img.w * img.h; ++i) {
    img.data[i] /= scene.num_iters;
  }
}
IMAGE render_vcm(SCENE scene) {
  IMAGE img = image_init(3, scene.wi, scene.hi);
  omp_set_num_threads(scene.num_threads);
  PRNG_STATE states[omp_get_max_threads()];
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    states[i] = prng_seed(SEED);
    for (int j = 0; j < i; ++j) {
      states[i] = prng_jump(states[i]);
    }
  }
  vcm(scene, img, states);
  return img;
}

SCENE scene_init0(void) {
  int wi, hi;
  double fovh, fovv, v, u;
  
  wi = 960 * 2;
  hi = 540 * 2;
  fovh = PI * 60.0 / 180.0;
  fovv = 2.0 * atan(tan(fovh / 2.0) * hi / wi);
  v = 1.0;
  u = 1.8;
  
  SCENE scene;
  scene.num_threads = 4;
  scene.num_samples = 1;
  scene.num_verts = 1024;
  scene.num_phots = 8000000;
  scene.num_nodes = 4000000;
  scene.wi = wi;
  scene.hi = hi;
  scene.ef = 0.5 * scene.wi * scene.hi;
  scene.wf = 2.0 * v * tan(fovh / 2.0);
  scene.hf = 2.0 * v * tan(fovv / 2.0);
  scene.df = v;
  scene.r = 0.05;
  scene.f = 1.0 / (1.0 / u + 1.0 / v);
  scene.num_iters = 1;
  scene.init_r = 1.0e-2;
  
  scene.num_objs = 10;
  scene.objs = ALLOC(OBJECT, scene.num_objs);
  
  MEDIUM m = medium_init(1.0, double3_init3(0.0, 0.0, 0.0));
  double r = 1.0e+3;
  double w = 1.0, h = 9.0 / 16.0, d = 3.0;
  
  scene.objs[0] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.5, 0.25, 0.25),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(-r - w, 0.0, 0.0)));
  scene.objs[1] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.25, 0.25, 0.5),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(r + w, 0.0, 0.0)));
  scene.objs[2] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.5, 0.5, 0.5),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, -r - h, 0.0)));
  scene.objs[3] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.5, 0.5, 0.5),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, r + h, 0.0)));
  scene.objs[4] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.25, 0.5, 0.25),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, 0.0, -r - d)));
  scene.objs[5] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.5, 0.5, 0.5),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, 0.0, r + d)));
  scene.objs[6] =
  object_init(material_init(
  double3_scale(50.0 / (4.0 * PI * sqr(0.1)), double3_init3(255 / 255.0, 64 / 255.0, 64 / 255.0)),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_NULL,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(0.1, double3_init3((w / 2.0 - 0.1) / 2.0 - 0.5, h - 0.1, -2.0)));
  scene.objs[7] =
  object_init(material_init(
  double3_scale(50.0 / (4.0 * PI * sqr(0.1)), double3_init3(64 / 255.0, 64 / 255.0, 255 / 255.0)),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_NULL,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(0.1, double3_init3((0.1 - w / 2.0) / 2.0 + 0.5, h - 0.1, -2.0)));
  scene.objs[8] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(255 / 255.0, 215 / 255.0, 50 / 255.0),
  EDF_UNI, BSDF_CONDUCTOR,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(0.3, double3_init3((w / 2.0 - 0.3) / 2.0 - 0.5, 0.3 - h, -2.25)));
  scene.objs[9] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_DIELECTRIC,
  m,
  medium_init(1.5, double3_init3(-log(0.99), -log(0.99), -log(0.99)))),
  sphere_init(0.3, double3_init3((0.3 - w / 2.0) / 2.0 + 0.5, 0.3 - h, -2.0)));
  return scene;
}
SCENE scene_init1(void) {
  int wi, hi;
  double fovh, fovv, v, u;
  
  wi = 960 * 4;
  hi = 540 * 4;
  fovh = PI * 60.0 / 180.0;
  fovv = 2.0 * atan(tan(fovh / 2.0) * hi / wi);
  v = 1.0;
  u = 10.0;
  
  SCENE scene;
  scene.num_threads = 8;
  scene.num_samples = 4;
  scene.num_verts = 1024;
  scene.num_phots = 80000000;
  scene.num_nodes = 20000000;
  scene.wi = wi;
  scene.hi = hi;
  scene.ef = 1.0 * scene.wi * scene.hi;
  scene.wf = 2.0 * v * tan(fovh / 2.0);
  scene.hf = 2.0 * v * tan(fovv / 2.0);
  scene.df = v;
  scene.r = 0.001;
  scene.f = 1.0 / (1.0 / u + 1.0 / v);
  scene.num_iters = 1;
  scene.init_r = 1.0e-3;
  
  int nt = 32, nh = 8;
  scene.num_objs = 2 + nt * nh;
  scene.objs = ALLOC(OBJECT, scene.num_objs);
  
  MEDIUM m = medium_init(1.0, double3_init3(0.0, 0.0, 0.0));
  
  PRNG_STATE state = prng_seed(1);
  
  double r0 = 1.0;
  double r1 = r0 * sqrt(sqr(cos(2.0 * PI / nt) - 1.0) + sqr(sin(2.0 * PI / nt))) / 2.0;
  double r2 = 2.0;
  double a = PI * 0.0 / 180.0;
  
  double r = 10.0, theta = 60.0 * PI / 180.0, phi = 180.0 * PI / 180.0;
  double cost = cos(theta), sint = sin(theta), cosp = cos(phi), sinp = sin(phi);
  
  DOUBLE3 x, y, z, t;
  x = double3_init3(cosp, -cost * sinp, sint * sinp);
  y = double3_init3(0.0, sint, cost);
  z = double3_init3(-sinp, -cost * cosp, sint * cosp);
  t = double3_init3(0.0, 0.0, -r);
  
  scene.objs[0] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_LAMBERT,
  m,
  m),
  sphere_init(1.0e+3, double3_init3(0.0, -1.0e+3 - r1, 0.0)));
  scene.objs[1] =
  object_init(material_init(
  double3_scale(1.0e+2 / (4.0 * PI * sqr(1.0e-6)), double3_init3(255 / 255.0, 255 / 255.0, 255 / 255.0)),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_NULL,
  m,
  m),
  sphere_init(1.0e-6, double3_init3(0.0, r2 * cos(a), r2 * sin(a))));
  for (int i = 0; i < nh; ++i) {
    for (int j = 0; j < nt; ++j) {
      double theta = 2.0 * PI * j / nt + PI * i / nt;;
      double h = sqrt(3.0) * r1 * i;
      double cost = cos(theta), sint = sin(theta);
      MEDIUM mn;
      mn.n = 1.1 + 0.9 * prng_next_double(&state);
      mn.a.x = -1.0 * log(1.0e-3 + 0.1 * prng_next_double(&state));
      mn.a.y = -1.0 * log(1.0e-3 + 0.1 * prng_next_double(&state));
      mn.a.z = -1.0 * log(1.0e-3 + 0.1 * prng_next_double(&state));
      OBJECT obj;
      obj.m.e = double3_init(0.0);
      obj.m.s = double3_init(1.0);
      obj.m.edf = EDF_UNI;
      obj.m.bsdf = BSDF_DIELECTRIC;
      obj.m.mp = m;
      obj.m.mn = mn;
      obj.s.r = r1;
      obj.s.c = double3_init3(r0 * sint, h, r0 * cost);
      scene.objs[2 + j + i * nt] = obj;
    }
  }
  for (int i = 0; i < scene.num_objs; ++i) {
    scene.objs[i].s.c = double3_add(t, double3_trafo(x, y, z, scene.objs[i].s.c));
  }
  return scene;
}
SCENE scene_init2(void) {
  int wi, hi;
  double fovh, fovv, v, u;
  
  wi = 960 * 8;
  hi = 540 * 8;
  fovh = PI * 60.0 / 180.0;
  fovv = 2.0 * atan(tan(fovh / 2.0) * hi / wi);
  v = 1.0;
  u = 15.0;
  
  SCENE scene;
  scene.num_threads = 9;
  scene.num_samples = 1;
  scene.num_verts = 1024;
  scene.num_phots = 64000000;
  scene.num_nodes = 12800000;
  scene.wi = wi;
  scene.hi = hi;
  scene.ef = 1.0 * scene.wi * scene.hi;
  scene.wf = 2.0 * v * tan(fovh / 2.0);
  scene.hf = 2.0 * v * tan(fovv / 2.0);
  scene.df = v;
  scene.r = 0.001;
  scene.f = 1.0 / (1.0 / u + 1.0 / v);
  scene.num_iters = 1;
  scene.init_r = 1.0e-2;
  
  int n = 256;
  scene.num_objs = 3 + n;
  scene.objs = ALLOC(OBJECT, scene.num_objs);
  
  MEDIUM m = medium_init(1.0, double3_init3(0.0, 0.0, 0.0));
  MEDIUM m0 = medium_init(1.35, double3_init(0.0));
  PRNG_STATE state = prng_seed(2);
  
  double r0 = 2.0;
  double r1 = 3.0;
  double a = PI * 300.0 / 180.0;
  
  double r = 15.0, theta = 75.0 * PI / 180.0, phi = 45.0 * PI / 180.0;
  double cost = cos(theta), sint = sin(theta), cosp = cos(phi), sinp = sin(phi);
  
  DOUBLE3 x, y, z, t;
  x = double3_init3(cosp, -cost * sinp, sint * sinp);
  y = double3_init3(0.0, sint, cost);
  z = double3_init3(-sinp, -cost * cosp, sint * cosp);
  t = double3_init3(0.0, 0.0, -r);
  
  scene.objs[0] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_LAMBERT,
  m,
  m),
  sphere_init(1.0e+3, double3_init3(0.0, -1.0e+3 - r0, 0.0)));
  scene.objs[1] =
  object_init(material_init(
  double3_scale(1.0e+3 / (4.0 * PI * sqr(0.001)), double3_init3(255 / 255.0, 255 / 255.0, 255 / 255.0)),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_NULL,
  m,
  m),
  sphere_init(0.001, double3_init3(0.0, r1 * cos(a), r1 * sin(a))));
  scene.objs[2] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_DIELECTRIC,
  m,
  m0),
  sphere_init(r0, double3_init3(0.0, 0.0, 0.0)));
  for (int i = 0; i < n;) {
    double u0 = prng_next_double(&state);
    double u1 = prng_next_double(&state);
    double u2 = prng_next_double(&state);
    double r, cost, sint, cosp, sinp;
    r = r0 * pow(u0, 1.0 / 3.0);
    cost = 2.0 * u1 - 1.0;
    sint = opst(cost);
    double phi = 2.0 * PI * u2;
    cosp = cos(phi);
    sinp = sin(phi);
    SPHERE s;
    s.r = r0 * prng_next_double(&state);
    s.c = double3_init3(r * sint * cosp, r * sint * sinp, r * cost);
    if (double3_len(s.c) + s.r >= r0) {
      continue;
    }
    bool cont = false;
    for (int j = 0; j < i; ++j) {
      SPHERE sj = scene.objs[3 + j].s;
      double d = double3_len(double3_sub(s.c, sj.c));
      if (s.r + sj.r >= d) {
        cont = true;
        break;
      }
    }
    if (cont) {
      continue;
    }
    MEDIUM mn;
    mn.n = 1.5 + prng_next_double(&state);
    mn.a.x = -1.0 / s.r * log(1.0e-3 + prng_next_double(&state));
    mn.a.y = -1.0 / s.r * log(1.0e-3 + prng_next_double(&state));
    mn.a.z = -1.0 / s.r * log(1.0e-3 + prng_next_double(&state));
    OBJECT obj;
    obj.m.e = double3_init(0.0);
    obj.m.s = double3_init(1.0);
    obj.m.edf = EDF_UNI;
    obj.m.bsdf = BSDF_DIELECTRIC;
    obj.m.mp = m0;
    obj.m.mn = mn;
    obj.s = s;
    scene.objs[3 + i] = obj;
    ++i;
  }
  for (int i = 0; i < scene.num_objs; ++i) {
    scene.objs[i].s.c = double3_add(t, double3_trafo(x, y, z, scene.objs[i].s.c));
  }
  return scene;
}
SCENE scene_init3(void) {
  int wi, hi;
  double fovh, fovv, v, u;
  
  wi = 960 * 2;
  hi = 540 * 2;
  fovh = PI * 60.0 / 180.0;
  fovv = 2.0 * atan(tan(fovh / 2.0) * hi / wi);
  v = 1.0;
  u = 1.8;
  
  SCENE scene;
  scene.num_threads = 8;
  scene.num_samples = 1;
  scene.num_verts = 1024;
  scene.num_phots = 8000000;
  scene.num_nodes = 4000000;
  scene.wi = wi;
  scene.hi = hi;
  scene.ef = 3.0 * scene.wi * scene.hi;
  scene.wf = 2.0 * v * tan(fovh / 2.0);
  scene.hf = 2.0 * v * tan(fovv / 2.0);
  scene.df = v;
  scene.r = 0.01;
  scene.f = 1.0 / (1.0 / u + 1.0 / v);
  scene.num_iters = 1;
  scene.init_r = 1.0e-3;
  
  scene.num_objs = 11;
  scene.objs = ALLOC(OBJECT, scene.num_objs);
  
  MEDIUM m = medium_init(1.0, double3_init3(0.0, 0.0, 0.0));
  double r = 1.0e+3;
  double w = 1.0, h = 9.0 / 16.0, d = 3.0;
  
  scene.objs[0] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.9, 0.3, 0.3),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(-r - w, 0.0, 0.0)));
  scene.objs[1] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.3, 0.3, 0.9),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(r + w, 0.0, 0.0)));
  scene.objs[2] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.6, 0.6, 0.6),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, -r - h, 0.0)));
  scene.objs[3] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.6, 0.6, 0.6),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, r + h, 0.0)));
  scene.objs[4] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.3, 0.9, 0.3),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, 0.0, -r - d)));
  scene.objs[5] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.6, 0.6, 0.6),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, 0.0, r + d)));
  scene.objs[6] =
  object_init(material_init(
  double3_scale(5.0 / (4.0 * PI * sqr(0.001)), double3_init3(0.9, 0.9, 0.9)),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_NULL,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(0.001, double3_init3(0.0, h - 0.001, -2.0)));
  scene.objs[7] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_DIELECTRIC,
  m,
  medium_init(1.5, double3_init3(3.0, 1.0, 1.0))),
  sphere_init(0.3, double3_init3((w / 2.0 - 0.3) / 2.0 - 0.5, 0.3 - h, -2.25)));
  scene.objs[8] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.9, 0.9, 0.9),
  EDF_UNI, BSDF_LAMBERT,
  medium_init(1.5, double3_init3(3.0, 1.0, 1.0)),
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(0.1, double3_init3((w / 2.0 - 0.3) / 2.0 - 0.5, 0.3 - h, -2.25)));
  scene.objs[9] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_DIELECTRIC,
  m,
  medium_init(1.5, double3_init3(1.0, 1.0, 3.0))),
  sphere_init(0.3, double3_init3((0.3 - w / 2.0) / 2.0 + 0.5, 0.3 - h, -2.0)));
  scene.objs[10] =
  object_init(material_init(
  double3_scale(5.0 / (4.0 * PI * sqr(0.1)), double3_init3(0.9, 0.9, 0.9)),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_NULL,
  medium_init(1.5, double3_init3(1.0, 1.0, 3.0)),
  medium_init(1.5, double3_init3(1.0, 1.0, 3.0))),
  sphere_init(0.1, double3_init3((0.3 - w / 2.0) / 2.0 + 0.5, 0.3 - h, -2.0)));
  return scene;
}
SCENE scene_init4(void) {
  int wi, hi;
  double fovh, fovv, v, u;
  
  wi = 960 * 2;
  hi = 540 * 2;
  fovh = PI * 60.0 / 180.0;
  fovv = 2.0 * atan(tan(fovh / 2.0) * hi / wi);
  v = 1.0;
  u = 1.8;
  
  SCENE scene;
  scene.num_threads = 8;
  scene.num_samples = 1;
  scene.num_verts = 1024;
  scene.num_phots = 14000000;
  scene.num_nodes = 10000000;
  scene.wi = wi;
  scene.hi = hi;
  scene.ef = 15.0 * scene.wi * scene.hi;
  scene.wf = 2.0 * v * tan(fovh / 2.0);
  scene.hf = 2.0 * v * tan(fovv / 2.0);
  scene.df = v;
  scene.r = 0.01;
  scene.f = 1.0 / (1.0 / u + 1.0 / v);
  scene.num_iters = 1;
  scene.init_r = 1.0e-3;
  
  scene.num_objs = 10;
  scene.objs = ALLOC(OBJECT, scene.num_objs);
  
  MEDIUM m = medium_init(1.0, double3_init3(0.0, 0.0, 0.0));
  double r = 1.0e+3;
  double w = 1.0, h = 9.0 / 16.0, d = 3.0;
  
  scene.objs[0] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.9, 0.9, 0.3),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(-r - w, 0.0, 0.0)));
  scene.objs[1] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.9, 0.3, 0.9),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(r + w, 0.0, 0.0)));
  scene.objs[2] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.9, 0.9, 0.9),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, -r - h, 0.0)));
  scene.objs[3] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.9, 0.9, 0.9),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, r + h, 0.0)));
  scene.objs[4] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.3, 0.9, 0.9),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, 0.0, -r - d)));
  scene.objs[5] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(0.9, 0.9, 0.9),
  EDF_UNI, BSDF_LAMBERT,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(r, double3_init3(0.0, 0.0, r + d)));
  scene.objs[6] =
  object_init(material_init(
  double3_scale(1.0 / (4.0 * PI * sqr(0.001)), double3_init3(0.9, 0.9, 0.9)),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_NULL,
  m,
  medium_init(1.0, double3_init3(0.0, 0.0, 0.0))),
  sphere_init(0.001, double3_init3(0.0, h - 0.001 - 1.0e-3, -2.0)));
  scene.objs[7] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_DIELECTRIC,
  m,
  medium_init(1.5, double3_init3(2.0, 2.0, 0.0))),
  sphere_init(0.25, double3_init3(-0.625, 0.0, -2.0)));
  scene.objs[8] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_DIELECTRIC,
  m,
  medium_init(1.5, double3_init3(2.0, 0.0, 2.0))),
  sphere_init(0.25, double3_init3(0.0, 0.0, -2.0)));
  scene.objs[9] =
  object_init(material_init(
  double3_init3(0.0, 0.0, 0.0),
  double3_init3(1.0, 1.0, 1.0),
  EDF_UNI, BSDF_DIELECTRIC,
  m,
  medium_init(1.5, double3_init3(0.0, 2.0, 2.0))),
  sphere_init(0.25, double3_init3(0.625, 0.0, -2.0)));
  return scene;
}
void scene_prep(SCENE *scene) {
  scene->pdf = ALLOC(double, scene->num_objs);
  scene->cdf = ALLOC(double, scene->num_objs);
  
  double sum = 0.0;
  for (int i = 0; i < scene->num_objs; ++i) {
    OBJECT obj = scene->objs[i];
    double A = sphere_area(obj.s);
    double I = double3_dot(double3_init3(0.2126, 0.7152, 0.0722), obj.m.e);
    double W = A * I;
    scene->pdf[i] = W;
    sum += W;
    scene->cdf[i] = sum;
  }
  if (sum == 0.0) {
    scene->pdf[0] = 1.0;
    scene->cdf[0] = 1.0;
    for (int i = 1; i < scene->num_objs; ++i) {
      scene->pdf[i] = 0.0;
      scene->cdf[i] = 0.0;
    }
  } else {
    for (int i = 0; i < scene->num_objs; ++i) {
      scene->pdf[i] /= sum;
      scene->cdf[i] /= sum;
    }
  }
  for (int i = 0; i < scene->num_objs; ++i) {
    fprintf(stdout, "%d : pdf = %f, cdf = %f\n", i, scene->pdf[i], scene->cdf[i]);
  }
}
void scene_del(SCENE scene) {
  free(scene.objs);
  free(scene.pdf);
  free(scene.cdf);
}

int main(int argc, char **argv) {
  SCENE scene = scene_init2();
  scene_prep(&scene);
  fprintf(stdout, "rendering...\n");
  double t = time_get();
  IMAGE img = render_vcm(scene);
  t = time_get() - t;
  fprintf(stdout, "render time = %f sec\n", t);
  scene_del(scene);
  if (image_save_pfm(img, DEST_PATH)) {
    fprintf(stdout, "file saved !\n\a");
  } else {
    fprintf(stdout, "file unsaved !\n\a");
  }
  image_del(img);
  return 0;
}