#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define pi 3.141592653589793
#define eps 1.0e-12
#define inf INFINITY
typedef char i1;
typedef int8_t i1s;
typedef int16_t i2s;
typedef int32_t i4s;
typedef int64_t i8s;
typedef uint8_t i1u;
typedef uint16_t i2u;
typedef uint32_t i4u;
typedef uint64_t i8u;
typedef float f4;
typedef double f8;
typedef struct i4s2 i4s2;
typedef struct i4s3 i4s3;
typedef struct i8u2 i8u2;
typedef struct f82 f82;
typedef struct f83 f83;
typedef struct f832 f832;
typedef struct f833 f833;
typedef struct med med;
typedef struct sdf sdf;
typedef struct face face;
typedef struct scene scene;
struct i4s2 {
	i4s x, y;
};
struct i4s3 {
	i4s x, y, z;
};
struct i8u2 {
	i8u x, y;
};
struct f82 {
	f8 x, y;
};
struct f83 {
	f8 x, y, z;
};
struct f832 {
	f83 x, y;
};
struct f833 {
	f83 x, y, z;
};
struct med {
	f8 n;
	f83 a;
};
struct sdf {
	i4s t;
	f8 ps[8];
};
struct face {
	i4s ei, bi;
	i4s2 mis;
	i4s3 vis, nis;
};
struct scene {
	i4s esn, bsn, msn, vsn, nsn, fsn;
	sdf *es, *bs;
	med *ms;
	f83 *vs, *ns;
	face *fs;
};
f8 timer(void) {
	return clock() / (f8)(CLOCKS_PER_SEC);
}
i4s i4s_min(i4s a, i4s b) {
	return a < b ? a : b;
}
i4s i4s_max(i4s a, i4s b) {
	return a > b ? a : b;
}
i4s i4s_cl(i4s a, i4s b, i4s x) {
	return i4s_min(i4s_max(x, a), b);
}
f8 f8_min(f8 a, f8 b) {
	return fmin(a, b);
}
f8 f8_max(f8 a, f8 b) {
	return fmax(a, b);
}
f8 f8_abs(f8 x) {
	return fabs(x);
}
f8 f8_cl(f8 a, f8 b, f8 x) {
	return f8_min(f8_max(x, a), b);
}
f8 f8_sign(f8 x) {
	return copysign(1.0, x);
}
f8 f8_p2(f8 x) {
	return x * x;
}
f8 f8_p3(f8 x) {
	return x * x * x;
}
f8 f8_p4(f8 x) {
	return x * x * x * x;
}
f8 f8_sqrt(f8 x) {
	return sqrt(x);
}
f8 f8_cos(f8 x) {
	return cos(x);
}
f8 f8_sin(f8 x) {
	return sin(x);
}
f8 f8_exp(f8 x) {
	return exp(x);
}
f8 f8_op2u(f8 x) {
	return f8_sqrt(1.0 - f8_p2(x));
}
f8 f8_op2us(f8 x) {
	return f8_sqrt(f8_max(1.0 - f8_p2(x), 0.0));
}
f8 f8_op3u(f8 x, f8 y) {
	return f8_sqrt(1.0 - f8_p2(x) - f8_p2(y));
}
f8 f8_op3us(f8 x, f8 y) {
	return f8_sqrt(f8_max(1.0 - f8_p2(x) - f8_p2(y), 0.0));
}
i4s2 i4s2_set(i4s x, i4s y) {
	i4s2 a = {x, y};
	return a;
}
i8u2 i8u2_set(i8u x, i8u y) {
	i8u2 a = {x, y};
	return a;
}
f82 f82_set(f8 x, f8 y) {
	f82 a = {x, y};
	return a;
}
f83 f83_set(f8 x, f8 y, f8 z) {
	f83 a = {x, y, z};
	return a;
}
f833 f833_set(f83 x, f83 y, f83 z) {
	f833 a = {x, y, z};
	return a;
}
f83 f83_set1(f8 a) {
	return f83_set(a, a, a);
}
f83 f83_add(f83 a, f83 b) {
	return f83_set(a.x + b.x, a.y + b.y, a.z + b.z);
}
f83 f83_sub(f83 a, f83 b) {
	return f83_set(a.x - b.x, a.y - b.y, a.z - b.z);
}
f83 f83_mul(f83 a, f83 b) {
	return f83_set(a.x * b.x, a.y * b.y, a.z * b.z);
}
f83 f83_div(f83 a, f83 b) {
	return f83_set(a.x / b.x, a.y / b.y, a.z / b.z);
}
f83 f83_scl(f8 a, f83 b) {
	return f83_set(a * b.x, a * b.y, a * b.z);
}
f83 f83_crs(f83 a, f83 b) {
	return f83_set(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
f8 f83_dot(f83 a, f83 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
f8 f83_len2(f83 a) {
	return f83_dot(a, a);
}
f8 f83_len(f83 a) {
	return f8_sqrt(f83_len2(a));
}
f83 f83_norm(f83 a) {
	return f83_scl(1.0 / f83_len(a), a);
}
f83 f83_exp(f83 a) {
	return f83_set(exp(a.x), exp(a.y), exp(a.z));
}
f83 f83_min(f83 a, f83 b) {
	return f83_set(f8_min(a.x, b.x), f8_min(a.y, b.y), f8_min(a.z, b.z));
}
f83 f83_max(f83 a, f83 b) {
	return f83_set(f8_max(a.x, b.x), f8_max(a.y, b.y), f8_max(a.z, b.z));
}
f8 f83_minc(f83 a) {
	return f8_min(f8_min(a.x, a.y), a.z);
}
f8 f83_maxc(f83 a) {
	return f8_max(f8_max(a.x, a.y), a.z);
}
f833 f83_onb(f83 n) {
	f8 s = f8_sign(n.z);
	f8 a = -1.0 / (s + n.z);
	f8 b = n.x * n.y * a;
	return f833_set(f83_set(1.0 + s * n.x * n.x * a, s * b, -s * n.x), f83_set(b, s + n.y * n.y * a, -n.y), n);
}
void f83_print(f83 a) {
	printf("%.*f %.*f %.*f\n", DECIMAL_DIG, a.x, DECIMAL_DIG, a.y, DECIMAL_DIG, a.z);
	return;
}
f83 f833_trafo(f833 a, f83 b) {
	return f83_set(a.x.x * b.x + a.y.x * b.y + a.z.x * b.z, a.x.y * b.x + a.y.y * b.y + a.z.y * b.z, a.x.z * b.x + a.y.z * b.y + a.z.z * b.z);
}
f833 f833_trapo(f833 a) {
	return f833_set(f83_set(a.x.x, a.y.x, a.z.x), f83_set(a.x.y, a.y.y, a.z.y), f83_set(a.x.z, a.y.z, a.z.z));
}
f8 rand_i8u_f8(i8u x) {
	return 0x1.0p-53 * (x >> 11);
}
i8u sm64_i8u(i8u *s) {
	i8u x = *s += 0x9E3779B97F4A7C15;
	x = 0xBF58476D1CE4E5B9 * (x ^ (x >> 30));
	x = 0x94D049BB133111EB * (x ^ (x >> 27));
	return x ^ (x >> 31);
}
i8u x128p_rotl(i8u a, i8u b) {
	return (a << b) | (a >> (64 - b));
}
i8u x128p_i8u(i8u2 *s) {
	i8u2 a = *s;
	i8u r = a.x + a.y;
	a.y ^= a.x;
	s->x = x128p_rotl(a.x, 55) ^ a.y ^ (a.y << 14);
	s->y = x128p_rotl(a.y, 36);
	return r;
}
f8 x128p_f8(i8u2 *s) {
	return rand_i8u_f8(x128p_i8u(s));
}
i8u2 x128p_seed(i8u x) {
	return i8u2_set(sm64_i8u(&x), sm64_i8u(&x));
}
i8u2 x128p_jump(i8u2 s) {
	i8u2 a = i8u2_set(0xBEAC0467EBA5FACB, 0xD86B048B86AA9922);
	i8u2 b = i8u2_set(0, 0);
	for (i4s i = 0; i < 64; ++i) {
		if (a.x & 1 << i) {
			b.x ^= s.x;
			b.y ^= s.y;
		}
		x128p_i8u(&s);
	}
	for (i4s i = 0; i < 64; ++i) {
		if (a.y & 1 << i) {
			b.x ^= s.x;
			b.y ^= s.y;
		}
		x128p_i8u(&s);
	}
	return b;
}
i8u2 x128p_jumpn(i8u2 s, i4s n) {
	for (i4s i = 0; i < n; ++i) {
		s = x128p_jump(s);
	}
	return s;
}
void scene_free(scene s) {
	free(s.es);
	free(s.bs);
	free(s.ms);
	free(s.vs);
	free(s.ns);
	free(s.fs);
	return;
}
scene scene_read(i1 *p) {
	FILE *f = fopen(p, "rb");
	if (f == NULL) {
		puts("i/o err");
		exit(1);
	}
	i4s ln = 256;
	i1 l[ln];
	i4s kn = 256;
	i1 k[kn];
	scene s;
	s.esn = 0;
	s.bsn = 0;
	s.msn = 0;
	s.vsn = 0;
	s.nsn = 0;
	s.fsn = 0;
	while (fgets(l, ln, f)) {
		sscanf(l, "%s", k);
		if (strcmp(k, "e") == 0) {
			++s.esn;
			continue;
		}
		if (strcmp(k, "b") == 0) {
			++s.bsn;
			continue;
		}
		if (strcmp(k, "m") == 0) {
			++s.msn;
			continue;
		}
		if (strcmp(k, "v") == 0) {
			++s.vsn;
			continue;
		}
		if (strcmp(k, "n") == 0) {
			++s.nsn;
			continue;
		}
		if (strcmp(k, "f") == 0) {
			++s.fsn;
			continue;
		}
	}
	s.es = malloc(s.esn * sizeof(sdf));
	s.bs = malloc(s.bsn * sizeof(sdf));
	s.ms = malloc(s.msn * sizeof(med));
	s.vs = malloc(s.vsn * sizeof(f83));
	s.ns = malloc(s.nsn * sizeof(f83));
	s.fs = malloc(s.fsn * sizeof(face));
	s.esn = 0;
	s.bsn = 0;
	s.msn = 0;
	s.vsn = 0;
	s.nsn = 0;
	s.fsn = 0;
	rewind(f);
	while (fgets(l, ln, f)) {
		sscanf(l, "%s", k);
		if (strcmp(k, "e") == 0) {
			sdf a;
			sscanf(l, "%*s%d%lf%lf%lf%lf%lf%lf%lf%lf", &a.t, a.ps + 0, a.ps + 1, a.ps + 2, a.ps + 3, a.ps + 4, a.ps + 5, a.ps + 6, a.ps + 7);
			s.es[s.esn++] = a;
		}
		if (strcmp(k, "b") == 0) {
			sdf a;
			sscanf(l, "%*s%d%lf%lf%lf%lf%lf%lf%lf%lf", &a.t, a.ps + 0, a.ps + 1, a.ps + 2, a.ps + 3, a.ps + 4, a.ps + 5, a.ps + 6, a.ps + 7);
			s.bs[s.bsn++] = a;
		}
		if (strcmp(k, "m") == 0) {
			med a;
			sscanf(l, "%*s%lf%lf%lf%lf", &a.n, &a.a.x, &a.a.y, &a.a.z);
			s.ms[s.msn++] = a;
		}
		if (strcmp(k, "v") == 0) {
			f83 a;
			sscanf(l, "%*s%lf%lf%lf", &a.x, &a.y, &a.z);
			s.vs[s.vsn++] = a;
		}
		if (strcmp(k, "n") == 0) {
			f83 a;
			sscanf(l, "%*s%lf%lf%lf", &a.x, &a.y, &a.z);
			s.ns[s.nsn++] = a;
		}
		if (strcmp(k, "f") == 0) {
			face a;
			sscanf(l, "%*s%d%d%d%d%d%d%d%d%d%d", &a.ei, &a.bi, &a.mis.x, &a.mis.y, &a.vis.x, &a.vis.y, &a.vis.z, &a.nis.x, &a.nis.y, &a.nis.z);
			s.fs[s.fsn++] = a;
		}
	}
	fclose(f);
	return s;
}
typedef struct conf conf;
struct conf {
	i1 ps[16];
	i1 pr[16];
	i1 t[16];
	i4s s;
	i4s tn;
	i4s sn;
	i4s en;
	i4s bnsna;
	i4s2 rsi;
	f82 rsfh;
	f8 lr;
	f8 lf;
	i1 pe[16];
	f833 m;
};
conf conf_read(i1 *p) {
	FILE *f = fopen(p, "rb");
	if (f == NULL) {
		puts("io err");
		exit(1);
	}
	i4s ln = 256;
	i1 l[ln];
	i4s kn = 256;
	i1 k[kn];
	conf a;
	while (fgets(l, ln, f)) {
		sscanf(l, "%s", k);
		if (strcmp(k, "ps") == 0) {
			sscanf(l, "%*s%s", a.ps);
			continue;
		}
		if (strcmp(k, "pr") == 0) {
			sscanf(l, "%*s%s", a.pr);
			continue;
		}
		if (strcmp(k, "t") == 0) {
			sscanf(l, "%*s%s", a.t);
			continue;
		}
		if (strcmp(k, "s") == 0) {
			sscanf(l, "%*s%d", &a.s);
			continue;
		}
		if (strcmp(k, "tn") == 0) {
			sscanf(l, "%*s%d", &a.tn);
			continue;
		}
		if (strcmp(k, "sn") == 0) {
			sscanf(l, "%*s%d", &a.sn);
			continue;
		}
		if (strcmp(k, "en") == 0) {
			sscanf(l, "%*s%d", &a.en);
			continue;
		}
		if (strcmp(k, "bnsna") == 0) {
			sscanf(l, "%*s%d", &a.bnsna);
			continue;
		}
		if (strcmp(k, "rsi") == 0) {
			sscanf(l, "%*s%d%d", &a.rsi.x, &a.rsi.y);
			continue;
		}
		if (strcmp(k, "rsfh") == 0) {
			sscanf(l, "%*s%lf%lf", &a.rsfh.x, &a.rsfh.y);
			continue;
		}
		if (strcmp(k, "lr") == 0) {
			sscanf(l, "%*s%lf", &a.lr);
			continue;
		}
		if (strcmp(k, "lf") == 0) {
			sscanf(l, "%*s%lf", &a.lf);
			continue;
		}
		if (strcmp(k, "pe") == 0) {
			sscanf(l, "%*s%s", a.pe);
			continue;
		}
		if (strcmp(k, "m") == 0) {
			sscanf(l, "%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf", &a.m.x.x, &a.m.x.y, &a.m.x.z, &a.m.y.x, &a.m.y.y, &a.m.y.z, &a.m.z.x, &a.m.z.y, &a.m.z.z);
			continue;
		}
	}
	fclose(f);
	return a;
}
typedef struct bvh_node bvh_node;
struct bvh_node {
	f832 b;
	i4s2 c;
};
f832 box_default(void) {
	f832 r = {f83_set1(inf), f83_set1(-inf)};
	return r;
}
f8 box_area(f832 b) {
	f83 d = f83_sub(b.y, b.x);
	return 2.0 * (d.x * d.y + d.y * d.z + d.z * d.x);
}
f832 box_merge(f832 a, f832 b) {
	f832 r = {f83_min(a.x, b.x), f83_max(a.y, b.y)};
	return r;
}
i4s box_isect(f832 b, f83 o, f83 di, f8 *t, f8 k) {
	f83 t0s = f83_mul(f83_sub(b.x, o), di);
	f83 t1s = f83_mul(f83_sub(b.y, o), di);
	f83 tmins = f83_min(t0s, t1s);
	f83 tmaxs = f83_max(t0s, t1s);
	f8 tmin = f83_maxc(tmins);
	f8 tmax = f83_minc(tmaxs);
	*t = f8_max(tmin, 0.0);
	return tmin < tmax + eps && tmax > eps && tmin < k - eps;
}
f8 tri_area(f833 a) {
	return 0.5 * f83_len(f83_crs(f83_sub(a.y, a.x), f83_sub(a.z, a.x)));
}
f8 tri_pdf(f833 a) {
	return 1.0 / tri_area(a);
}
f82 tri_sample(i8u2 *s) {
	f8 a = f8_sqrt(x128p_f8(s));
	return f82_set(1.0 - a, x128p_f8(s) * a);
}
f832 tri_box(f833 a) {
	f832 b = box_default();
	b.x = f83_min(b.x, a.x);
	b.x = f83_min(b.x, a.y);
	b.x = f83_min(b.x, a.z);
	b.y = f83_max(b.y, a.x);
	b.y = f83_max(b.y, a.y);
	b.y = f83_max(b.y, a.z);
	return b;
}
f83 tri_cent(f833 a) {
	return f83_scl(1.0 / 3.0, f83_add(f83_add(a.x, a.y), a.z));
}
f83 tri_interp(f833 a, f82 p) {
	return f83_add(f83_add(f83_scl(1.0 - p.x - p.y, a.x), f83_scl(p.x, a.y)), f83_scl(p.y, a.z));
}
f83 tri_ng(f833 a) {
	return f83_norm(f83_crs(f83_sub(a.y, a.x), f83_sub(a.z, a.x)));
}
f83 tri_ns(f833 a, f82 p) {
	return f83_norm(tri_interp(a, p));
}
i4s tri_isect(f833 vs, f83 o, f83 d, f8 *t, f82 *uv) {
	f83 e0 = f83_sub(vs.y, vs.x);
	f83 e1 = f83_sub(vs.z, vs.x);
	f83 pv = f83_crs(d, e1);
	f8 deti = f83_dot(e0, pv);
	if (f8_abs(deti) < eps) {
		return 0;
	}
	deti = 1.0 / deti;
	f83 tv = f83_sub(o, vs.x);
	f8 u = deti * f83_dot(tv, pv);
	if (u < 0.0 || u > 1.0) {
		return 0;
	}
	f83 qv = f83_crs(tv, e0);
	f8 v = deti * f83_dot(d, qv);
	if (v < 0.0 || u + v > 1.0) {
		return 0;
	}
	f8 td = deti * f83_dot(e1, qv);
	if (td < eps || td > *t - eps) {
		return 0;
	}
	*t = td;
	uv->x = u;
	uv->y = v;
	return 1;
}
f833 get_vs(f83 *vs, i4s3 vis) {
	return f833_set(vs[vis.x], vs[vis.y], vs[vis.z]);
}
typedef struct bvh_se bvh_se;
struct bvh_se {
	i4s i;
	f83 p;
};
i4s bvh_cmp_x(void const *n, void const *m) {
	bvh_se const *a = n;
	bvh_se const *b = m;
	return a->p.x < b->p.x ? -1 : 1;
}
i4s bvh_cmp_y(void const *n, void const *m) {
	bvh_se const *a = n;
	bvh_se const *b = m;
	return a->p.y < b->p.y ? -1 : 1;
}
i4s bvh_cmp_z(void const *n, void const *m) {
	bvh_se const *a = n;
	bvh_se const *b = m;
	return a->p.z < b->p.z ? -1 : 1;
}
void bvh_make(f83 *vs, face *fs, bvh_se *fis, f82 *ass, bvh_node *ns, i4s2 *st, i4s2 q, i4s pd, i4s *pe) {
	f8 tt = 1.0, tb = 1.0;
	f832 b = box_default();
	for (i4s i = q.x; i < q.y; ++i) {
		b = box_merge(b, tri_box(get_vs(vs, fs[fis[i].i].vis)));
	}
	f8 ai = 1.0 / box_area(b);
	f8 c = tt * (q.y - q.x);
	i4s axis = -1, piv;
	for (i4s i = 0; i < 3; ++i) {
		switch (i) {
			case 0 : {
				qsort(fis + q.x, q.y - q.x, sizeof(bvh_se), bvh_cmp_x);
				break;
			}
			case 1 : {
				qsort(fis + q.x, q.y - q.x, sizeof(bvh_se), bvh_cmp_y);
				break;
			}
			case 2 : {
				qsort(fis + q.x, q.y - q.x, sizeof(bvh_se), bvh_cmp_z);
				break;
			}
		}
		f832 bc;
		bc = box_default();
		for (i4s j = q.x; j < q.y - 1; ++j) {
			bc = box_merge(bc, tri_box(get_vs(vs, fs[fis[j].i].vis)));
			ass[j - q.x].x = box_area(bc);
		}
		bc = box_default();
		for (i4s j = q.y - 1; j > q.x; --j) {
			bc = box_merge(bc, tri_box(get_vs(vs, fs[fis[j].i].vis)));
			ass[j - q.x - 1].y = box_area(bc);
		}
		for (i4s j = q.x; j < q.y - 1; ++j) {
			f8 cd = tb + tt * ai * (ass[j - q.x].x * (j - q.x + 1) + ass[j - q.x].y * (q.y - j - 1));
			if (cd < c) {
				c = cd;
				axis = i;
				piv = j + 1;
			}
		}
	}
	if (axis == -1) {
		bvh_node n = {b, {q.x, -q.y}};
		ns[pd] = n;
		return;
	}
	switch (axis) {
		case 0 : {
			qsort(fis + q.x, q.y - q.x, sizeof(bvh_se), bvh_cmp_x);
			break;
		}
		case 1 : {
			qsort(fis + q.x, q.y - q.x, sizeof(bvh_se), bvh_cmp_y);
			break;
		}
		case 2 : {
			qsort(fis + q.x, q.y - q.x, sizeof(bvh_se), bvh_cmp_z);
			break;
		}
	}
	bvh_node n = {b, {*pe + 1, *pe + 2}};
	ns[pd] = n;
	st[(*pe)++] = i4s2_set(q.x, piv);
	st[(*pe)++] = i4s2_set(piv, q.y);
	return;
}
i4s bvh_makert(scene *s, i4s nsna, bvh_node *ns) {
	bvh_se *fis = malloc(sizeof(bvh_se) * s->fsn);
	for (i4s i = 0; i < s->fsn; ++i) {
		bvh_se a = {i, tri_cent(get_vs(s->vs, s->fs[i].vis))};
		fis[i] = a;
	}
	f82 *as = malloc(sizeof(f82) * (s->fsn - 1));
	i4s2 *st = malloc(sizeof(i4s2) * nsna);
	i4s pd = 0, pe = 0;
	i4s2 q = i4s2_set(0, s->fsn);
	while (1) {
		bvh_make(s->vs, s->fs, fis, as, ns, st, q, pd, &pe);
		if (pd == pe) {
			break;
		}
		q = st[pd++];
	}
	free(st);
	free(as);
	face *fs = malloc(sizeof(face) * s->fsn);
	for (i4s i = 0; i < s->fsn; ++i) {
		fs[i] = s->fs[fis[i].i];
	}
	free(fis);
	free(s->fs);
	s->fs = fs;
	return pd + 1;
}
void bvh_trav(f83 *vs, face *fs, bvh_node *ns, bvh_node n, i4s *st, i4s *p, f83 o, f83 d, f83 di, f8 *t, f82 *uv, i4s *id) {
	if (n.c.y < 0) {
		for (i4s i = n.c.x; i < -n.c.y; ++i) {
			if(tri_isect(get_vs(vs, fs[i].vis), o, d, t, uv)) {
				*id = i;
			}
		}
	} else {
		f8 tl;
		i4s hl = box_isect(ns[n.c.x].b, o, di, &tl, *t);
		f8 tr;
		i4s hr = box_isect(ns[n.c.y].b, o, di, &tr, *t);
		if (hl && hr) {
			if (tl < tr) {
				st[(*p)++] = n.c.y;
				st[(*p)++] = n.c.x;
			} else {
				st[(*p)++] = n.c.x;
				st[(*p)++] = n.c.y;
			}
		} else if (hl) {
			st[(*p)++] = n.c.x;
		} else if (hr) {
			st[(*p)++] = n.c.y;
		}
	}
	return;
}
i4s bvh_travrt(f83 *vs, face *fs, bvh_node *ns, i4s *st, f83 o, f83 d, f8 *t, f82 *uv, i4s *id) {
	*t = inf;
	*id = -1;
	f83 di = f83_div(f83_set1(1.0), d);
	bvh_node n = ns[0];
	f8 td;
	if (box_isect(n.b, o, di, &td, *t)) {
		i4s p = 0;
		while (1) {
			bvh_trav(vs, fs, ns, n, st, &p, o, d, di, t, uv, id);
			if (p == 0) {
				break;
			}
			n = ns[st[--p]];
		}
	}
	return *id >= 0;
}
i4s bvh_vis(f83 *vs, face *fs, bvh_node *ns, bvh_node n, i4s *st, i4s *p, f83 o, f83 d, f83 di, f8 *t) {
	if (n.c.y < 0) {
		for (i4s i = n.c.x; i < -n.c.y; ++i) {
			f82 uv;
			if(tri_isect(get_vs(vs, fs[i].vis), o, d, t, &uv)) {
				return 0;
			}
		}
	} else {
		f8 tl;
		i4s hl = box_isect(ns[n.c.x].b, o, di, &tl, *t);
		f8 tr;
		i4s hr = box_isect(ns[n.c.y].b, o, di, &tr, *t);
		if (hl && hr) {
			if (tl < tr) {
				st[(*p)++] = n.c.y;
				st[(*p)++] = n.c.x;
			} else {
				st[(*p)++] = n.c.x;
				st[(*p)++] = n.c.y;
			}
		} else if (hl) {
			st[(*p)++] = n.c.x;
		} else if (hr) {
			st[(*p)++] = n.c.y;
		}
	}
	return 1;
}
i4s bvh_visrt(f83 *vs, face *fs, bvh_node *ns, i4s *st, f83 o, f83 d, f8 t) {
	f83 di = f83_div(f83_set1(1.0), d);
	bvh_node n = ns[0];
	f8 td;
	if (box_isect(n.b, o, di, &td, t)) {
		i4s p = 0;
		while (1) {
			if (!bvh_vis(vs, fs, ns, n, st, &p, o, d, di, &t)) {
				return 0;
			}
			if (p == 0) {
				break;
			}
			n = ns[st[--p]];
		}
	}
	return 1;
}
f8 disc_area(f8 r) {
	return pi * f8_p2(r);
}
f8 disc_pdf(f8 r) {
	return 1.0 / disc_area(r);
}
f83 disc_sample(f8 rd, i8u2 *s) {
	f8 r = rd * f8_sqrt(x128p_f8(s)), t = 2.0 * pi * x128p_f8(s);
	return f83_scl(r, f83_set(f8_cos(t), f8_sin(t), 0.0));
}
f83 disc_ng(void) {
	return f83_set(0.0, 0.0, 1.0);
}
i4s disc_isect(f8 r, f83 o, f83 d, f8 *t) {
	if (f8_abs(d.z) < eps) {
		return 0;
	}
	f8 td = -o.z / d.z;
	if (td < 0.0 || td > *t) {
		return 0;
	}
	f83 p = f83_add(o, f83_scl(td, d));
	if (f83_len2(p) > f8_p2(r)) {
		return 0;
	}
	*t = td;
	return 1;
}
void image_read(i1 *p, i4s2 *rs, f83 **cs) {
	FILE *f = fopen(p, "rb");
	if (f == NULL) {
		puts("i/o err");
		exit(1);
	}
	i4s ln = 256;
	i1 l[ln];
	fgets(l, ln, f);
	fgets(l, ln, f);
	sscanf(l, "%d%d", &rs->x, &rs->y);
	fgets(l, ln, f);
	float *ds = malloc(sizeof(float) * 3 * rs->x * rs->y);
	fread(ds, sizeof(float), 3 * rs->x * rs->y, f);
	*cs = malloc(sizeof(f83) * rs->x * rs->y);
	for (i4s i = 0; i < rs->x * rs->y; ++i) {
		(*cs)[i] = f83_set(ds[3 * i + 0], ds[3 * i + 1], ds[3 * i + 2]);
	}
	free(ds);
	fclose(f);
	return;
}
void image_write(i1 *p, i4s2 rs, f83 *cs) {
	FILE *f = fopen(p, "wb");
	if (f == NULL) {
		puts("i/o err");
		exit(1);
	}
	fprintf(f, "PF\n%d %d\n-1.0\n", rs.x, rs.y);
	f4 *ds = malloc(sizeof(f4) * 3 * rs.x * rs.y);
	for (int i = 0; i < rs.x * rs.y; ++i) {
		ds[3 * i + 0] = cs[i].x;
		ds[3 * i + 1] = cs[i].y;
		ds[3 * i + 2] = cs[i].z;
	}
	fwrite(ds, sizeof(f4), 3 * rs.x * rs.y, f);
	free(ds);
	fclose(f);
	return;
}
f83 image_interp(i4s2 rs, f83 *cs, f82 pif) {
	i4s2 pii = i4s2_set(i4s_cl(0, rs.x - 2, pif.x), i4s_cl(0, rs.y - 2, pif.y));
	f82 d = f82_set(pif.x - pii.x, pif.y - pii.y);
	f83 c = f83_set1(0.0);
	c = f83_add(c, f83_scl((1.0 - d.x) * (1.0 - d.y), cs[rs.x * pii.y + pii.x]));
	c = f83_add(c, f83_scl(d.x * (1.0 - d.y), cs[rs.x * pii.y + pii.x + 1]));
	c = f83_add(c, f83_scl((1.0 - d.x) * d.y, cs[rs.x * (pii.y + 1) + pii.x]));
	c = f83_add(c, f83_scl(d.x * d.y, cs[rs.x * (pii.y + 1) + pii.x + 1]));
	return c;
}
f83 imgf2film(f82 rsif, f82 rsfh, f82 pif) {
	return f83_set(rsfh.x * (1.0 - 2.0 * pif.x / rsif.x), rsfh.y * (1.0 - 2.0 * pif.y / rsif.y), 1.0);
}
i4s2 film2imgi(i4s2 rsi, f82 rsif, f82 rsfh, f82 pf) {
	return i4s2_set(i4s_cl(0, rsi.x - 1, 0.5 * rsif.x * (1.0 - pf.x / rsfh.x)), i4s_cl(0, rsi.y - 1, 0.5 * rsif.y * (1.0 - pf.y / rsfh.y)));
}
f83 pol2cartzu0(f8 cp, f8 sp, f8 ct, f8 st) {
	return f83_set(cp * st, sp * st, ct);
}
f83 pol2cartzu1(f82 p) {
	return pol2cartzu0(f8_cos(p.x), f8_sin(p.x), f8_cos(p.y), f8_sin(p.y));
}
f82 cart2polzu(f83 a) {
	return f82_set(atan2(a.y, a.x) + (a.y < 0.0 ? 2.0 * pi : 0.0), acos(f8_cl(-1.0, 1.0, a.z)));
}
f83 refl(f83 n, f83 d, f8 c) {
	return f83_sub(d, f83_scl(2.0 * c, n));
}
f83 refr(f83 n, f83 d, f8 e, f8 ci, f8 ct) {
  return f83_sub(f83_scl(e, d), f83_scl(e * ci - ct, n));
}
f8 rgb2lum(f83 a) {
  return f83_dot(f83_set(0.2126, 0.7152, 0.0722), a);
}
f8 fresd(f8 e, f8 ci, f8 ct) {
	f8 perp = (e * ci - ct) / (e * ci + ct);
	f8 para = (e * ct - ci) / (e * ct + ci);
	return 0.5 * (f8_p2(perp) + f8_p2(para));
}
f8 fresc(f8 ei, f8 k, f8 ci) {
	ci = f8_abs(ci);
	f8 s = 1.0 - f8_p2(ci);
	f8 i = f8_p2(ei) - f8_p2(k) - s;
	f8 w = f8_sqrt(f8_p2(i) + f8_p2(2.0 * ei * k));
	f8 a = f8_sqrt(f8_max(0.5 * (w + i), 0.0));
	f8 perp = (w + f8_p2(ci) - 2.0 * a * ci) / (w + f8_p2(ci) + 2.0 * a * ci);
	f8 para = (f8_p2(ci) * w + f8_p2(s) - 2.0 * a * ci * s) / (f8_p2(ci) * w + f8_p2(s) + 2.0 * a * ci * s);
	return 0.5 * (perp + perp * para);
}
f83 fresc3(f83 ei, f83 k, f8 ci) {
	return f83_set(fresc(ei.x, k.x, ci), fresc(ei.y, k.y, ci), fresc(ei.z, k.z, ci));
}
f83 ggx_sample_ndf(f83 n, f8 r, i8u2 *s) {
	f8 u0 = x128p_f8(s), u1 = x128p_f8(s);
	return f833_trafo(f83_onb(n), pol2cartzu1(f82_set(2.0 * pi * u1, atan(f8_sqrt(f8_p2(r) * u0 / (1.0 - u0))))));
}
f83 ggx_sample_vndf(f83 n, f8 r, f83 wi, i8u2 *s) {
	f833 m = f83_onb(n);
	f833 mi = f833_trapo(m);
	wi = f833_trafo(mi, wi);
	f83 vh = f83_norm(f83_set(r * wi.x, r * wi.y, wi.z));
	f833 a = f83_onb(vh);
	f8 rad = f8_sqrt(x128p_f8(s));
	f8 p = 2.0 * pi * x128p_f8(s);
	f83 t;
	t.x = rad * f8_cos(p);
	t.y = rad * f8_sin(p);
	f8 k = 0.5 * (1.0 - vh.z);
	t.y = (1.0 - k) * f8_op2us(t.x) + k * t.y;
	t.z = f8_op3us(t.x, t.y);
	f83 nh = f833_trafo(a, t);
	f83 ne = f83_norm(f83_set(r * nh.x, r * nh.y, f8_max(nh.z, 0.0)));
	return f833_trafo(m, ne);
}
f8 ggx_d(f8 r, f8 cnm) {
	return 1.0 / (pi * f8_p2(r * f8_p2(cnm) + f8_abs(1.0 - f8_p2(cnm)) / r));
}
f8 ggx_lam(f8 r, f8 ct) {
	if (ct > 1.0 - eps) {
		return 0.0;
	}
	if (ct < eps - 1.0) {
		return -1.0;
	}
	return (f8_sign(ct) * sqrt(1.0 + (1.0 - f8_p2(ct)) * f8_p2(r / ct)) - 1.0) / 2.0;
}
f8 ggx_g1(f8 r, f8 ct) {
	if (f8_abs(ct) < eps) {
		return 0.0;
	}
	return 1.0 / (1.0 + ggx_lam(r, ct));
}
f8 ggx_g2(f8 r, f8 cni, f8 cno) {
	return ggx_g1(r, cni) * ggx_g1(r, cno);
}
i4s edf_con(sdf e) {
	switch (e.t) {
		case 0 : {
			return 0;
		}
		case 1 : {
			return 1;
		}
		case 2 : {
			return 1;
		}
		case 3 : {
			return 0;
		}
		default : {
			return 0;
		}
	}
}
i4s edf_eval_b(sdf f, f8 cno, f83 *fs, f8 *psb) {
	switch (f.t) {
		case 1 : {
			if (cno < 0.0) {
				return 0;
			}
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			*fs = f83_scl(1.0 / pi, a);
			*psb = 1.0 / pi;
			return 1;
		}
		case 2 : {
			if (cno < 0.0) {
				return 0;
			}
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 r = f.ps[3];
			f8 d = ggx_d(r, cno);
			*fs = f83_scl(d, a);
			*psb = d;
			return 1;
		}
		default : {
			return 0;
		}
	}
}
i4s edf_eval(sdf f, f8 cno, f83 *fs) {
	switch (f.t) {
		case 1 : {
			if (cno < 0.0) {
				return 0;
			}
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			*fs = f83_scl(1.0 / pi, a);
			return 1;
		}
		case 2 : {
			if (cno < 0.0) {
				return 0;
			}
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 r = f.ps[3];
			f8 d = ggx_d(r, cno);
			*fs = f83_scl(d, a);
			return 1;
		}
		default : {
			return 0;
		}
	}
}
i4s edf_sample_f(sdf f, f83 n, i8u2 *s, f83 *wo, f83 *tp, f8 *psf) {
	switch (f.t) {
		case 1 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 p = 2.0 * pi * x128p_f8(s);
			f8 cp = f8_cos(p), sp = f8_sin(p), ct = f8_sqrt(x128p_f8(s)), st = f8_op2u(ct);
			*wo = f833_trafo(f83_onb(n), pol2cartzu0(cp, sp, ct, st));
			*tp = a;
			*psf = 1.0 / pi;
			return 1;
		}
		case 2 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 r = f.ps[3];
			*wo = ggx_sample_ndf(n, r, s);
			*tp = a;
			*psf = ggx_d(r, f83_dot(n, *wo));
			return 1;
		}
		case 3 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			*wo = n;
			*tp = a;
			*psf = 1.0;
			return 1;
		}
		default : {
			return 0;
		}
	}
}
i4s edf_sample(sdf f, f83 n, i8u2 *s, f83 *wo, f83 *tp) {
	switch (f.t) {
		case 1 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 p = 2.0 * pi * x128p_f8(s);
			f8 cp = f8_cos(p), sp = f8_sin(p), ct = f8_sqrt(x128p_f8(s)), st = f8_op2u(ct);
			*wo = f833_trafo(f83_onb(n), pol2cartzu0(cp, sp, ct, st));
			*tp = a;
			return 1;
		}
		case 2 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 r = f.ps[3];
			*wo = ggx_sample_ndf(n, r, s);
			*tp = a;
			return 1;
		}
		case 3 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			*wo = n;
			*tp = a;
			return 1;
		}
		default : {
			return 0;
		}
	}
}
i4s bsdf_con(sdf f) {
	switch (f.t) {
		case 0 : {
			return 0;
		}
		case 1 : {
			return 0;
		}
		case 2 : {
			return 0;
		}
		case 3 : {
			return 1;
		}
		case 4 : {
			return 1;
		}
		case 5 : {
			return 1;
		}
		case 6 : {
			return 1;
		}
		default : {
			return 0;
		}
	}
}
i4s bsdf_eval_fb(sdf f, f8 nr, f8 nt, f83 n, f83 wi, f83 wo, f8 cni, f8 cno, i4s rad, f83 *fs, f83 *rf, f83 *rb, f8 *psf, f8 *psb) {
	switch (f.t) {
		case 3 : {
			if (cni * cno > 0.0) {
				return 0;
			}
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			*fs = f83_scl(1.0 / pi, a);
			*rf = a;
			*rb = a;
			*psf = 1.0 / pi;
			*psb = 1.0 / pi;
			return 1;
		}
		case 4 : {
			if (cni * cno < 0.0) {
				return 0;
			}
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			*fs = f83_scl(1.0 / pi, a);
			*rf = a;
			*rb = a;
			*psf = 1.0 / pi;
			*psb = 1.0 / pi;
			return 1;
		}
		case 5 : {
			if (cni * cno > 0.0) {
				return 0;
			}
			f83 e = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f83 k = f83_set(f.ps[3], f.ps[4], f.ps[5]);
			f8 r = f.ps[6];
			f83 m = f83_norm(f83_sub(wo, wi));
			f8 cnm = f83_dot(n, m), cmi = f83_dot(m, wi);
			f83 fr = fresc3(f83_set(e.x / nr, e.y / nr, e.z / nr), k, cmi);
			f8 g1i = ggx_g1(r, cni), g1o = ggx_g1(r, cno);
			f8 d = ggx_d(r, cnm);
			*fs = f83_scl(g1i * g1o * d / f8_abs(4.0 * cni * cno), fr);
			*rf = f83_scl(g1o, fr);
			*rb = f83_scl(g1i, fr);
			*psf = g1i * d / f8_abs(4.0 * cni * cno);
			*psb = g1o * d / f8_abs(4.0 * cni * cno);
			return 1;
		}
		case 6 : {
			return 0;
		}
		default : {
			return 0;
		}
	}
}
i4s bsdf_eval_f(sdf f, f8 nr, f8 nt, f83 n, f83 wi, f83 wo, f8 cni, f8 cno, i4s rad, f83 *fs, f83 *rf, f8 *psf) {
	switch (f.t) {
		case 3 : {
			if (cni * cno > 0.0) {
				return 0;
			}
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			*fs = f83_scl(1.0 / pi, a);
			*rf = a;
			*psf = 1.0 / pi;
			return 1;
		}
		case 4 : {
			if (cni * cno < 0.0) {
				return 0;
			}
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			*fs = f83_scl(1.0 / pi, a);
			*rf = a;
			*psf = 1.0 / pi;
			return 1;
		}
		case 5 : {
			if (cni * cno > 0.0) {
				return 0;
			}
			f83 e = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f83 k = f83_set(f.ps[3], f.ps[4], f.ps[5]);
			f8 r = f.ps[6];
			f83 m = f83_norm(f83_sub(wo, wi));
			f8 cnm = f83_dot(n, m), cmi = f83_dot(m, wi);
			f83 fr = fresc3(f83_set(e.x / nr, e.y / nr, e.z / nr), k, cmi);
			f8 g1i = ggx_g1(r, cni), g1o = ggx_g1(r, cno);
			f8 d = ggx_d(r, cnm);
			*fs = f83_scl(g1i * g1o * d / f8_abs(4.0 * cni * cno), fr);
			*rf = f83_scl(g1o, fr);
			*psf = g1i * d / f8_abs(4.0 * cni * cno);
			return 1;
		}
		case 6 : {
			return 0;
		}
		default : {
			return 0;
		}
	}
}
i4s bsdf_sample_fb(sdf f, f8 nr, f8 nt, f83 n, f83 wi, i8u2 *s, f8 cni, i4s rad, f83 *wo, f83 *tp, f83 *rf, f83 *rb, f8 *psf, f8 *psb) {
	switch (f.t) {
		case 0 : {
			*wo = wi;
			*tp = f83_set1(1.0);
			*rf = f83_set1(1.0);
			*rb = f83_set1(1.0);
			*psf = 1.0 / f8_abs(cni);
			*psb = 1.0 / f8_abs(cni);
			return 1;
		}
		case 1 : {
			f83 e = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f83 k = f83_set(f.ps[3], f.ps[4], f.ps[5]);
			f83 fr = fresc3(f83_set(e.x / nr, e.y / nr, e.z / nr), k, cni);
			*wo = refl(n, wi, cni);
			*tp = fr;
			*rf = fr;
			*rb = fr;
			*psf = f8_abs(cni);
			*psb = f8_abs(cni);
			return 1;
		}
		case 2 : {
			f8 e = nr / nt;
			f8 e2 = f8_p2(e);
			f8 d = 1.0 - e2 * (1.0 - f8_p2(cni));
			f8 fr = 1.0;
			if (d > 0.0) {
				f8 cnt = f8_sign(cni) * f8_sqrt(d);
				fr = fresd(e, cni, cnt);
				if (x128p_f8(s) > fr) {
					*wo = refr(n, wi, e, cni, cnt);
					*tp = f83_set1(rad ? e2 : 1.0);
					*rf = f83_set1(1.0);
					*rb = f83_set1(1.0);
					*psf = (1.0 - fr) / f8_abs(cnt);
					*psb = (1.0 - fr) / f8_abs(cni);
					return 1;
				}
			}
			*wo = refl(n, wi, cni);
			*tp = f83_set1(1.0);
			*rf = f83_set1(1.0);
			*rb = f83_set1(1.0);
			*psf = fr / f8_abs(cni);
			*psb = fr / f8_abs(cni);
			return 1;
		}
		case 3 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 p = 2.0 * pi * x128p_f8(s);
			f8 cp = f8_cos(p), sp = f8_sin(p), ct = f8_sqrt(x128p_f8(s)), st = f8_op2u(ct);
			*wo = f833_trafo(f83_onb(f83_scl(-f8_sign(cni), n)), pol2cartzu0(cp, sp, ct, st));
			*tp = a;
			*rf = a;
			*rb = a;
			*psf = 1.0 / pi;
			*psb = 1.0 / pi;
			return 1;
		}
		case 4 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 p = 2.0 * pi * x128p_f8(s);
			f8 cp = f8_cos(p), sp = f8_sin(p), ct = f8_sqrt(x128p_f8(s)), st = f8_op2u(ct);
			*wo = f833_trafo(f83_onb(f83_scl(f8_sign(cni), n)), pol2cartzu0(cp, sp, ct, st));
			*tp = a;
			*rf = a;
			*rb = a;
			*psf = 1.0 / pi;
			*psb = 1.0 / pi;
			return 1;
		}
		case 5 : {
			f83 e = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f83 k = f83_set(f.ps[3], f.ps[4], f.ps[5]);
			f8 r = f.ps[6];
			f83 m = ggx_sample_vndf(n, r, f83_scl(f8_sign(cni), wi), s);
			f8 cnm = f83_dot(n, m), cmi = f83_dot(m, wi);
			*wo = refl(m, wi, cmi);
			f8 cno = f83_dot(n, *wo);
			if (cni * cno > 0.0) {
				return 0;
			}
			f83 fr = fresc3(f83_set(e.x / nr, e.y / nr, e.z / nr), k, cmi);
			f8 g1i = ggx_g1(r, cni), g1o = ggx_g1(r, cno);
			f8 d = ggx_d(r, cnm);
			*tp = f83_scl(g1o, fr);
			*rf = f83_scl(g1o, fr);
			*rb = f83_scl(g1i, fr);
			*psf = g1i * d / f8_abs(4.0 * cni * cno);
			*psb = g1o * d / f8_abs(4.0 * cni * cno);
			return 1;
		}
		case 6 : {
			return 0;
		}
		default : {
			return 0;
		}
	}
}
i4s bsdf_sample_f(sdf f, f8 nr, f8 nt, f83 n, f83 wi, i8u2 *s, f8 cni, i4s rad, f83 *wo, f83 *tp, f83 *rf, f8 *psf) {
	switch (f.t) {
		case 0 : {
			*wo = wi;
			*tp = f83_set1(1.0);
			*rf = f83_set1(1.0);
			*psf = 1.0 / f8_abs(cni);
			return 1;
		}
		case 1 : {
			f83 e = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f83 k = f83_set(f.ps[3], f.ps[4], f.ps[5]);
			f83 fr = fresc3(f83_set(e.x / nr, e.y / nr, e.z / nr), k, cni);
			*wo = refl(n, wi, cni);
			*tp = fr;
			*rf = fr;
			*psf = f8_abs(cni);
			return 1;
		}
		case 2 : {
			f8 e = nr / nt;
			f8 e2 = f8_p2(e);
			f8 d = 1.0 - e2 * (1.0 - f8_p2(cni));
			f8 fr = 1.0;
			if (d > 0.0) {
				f8 cnt = f8_sign(cni) * f8_sqrt(d);
				fr = fresd(e, cni, cnt);
				if (x128p_f8(s) > fr) {
					*wo = refr(n, wi, e, cni, cnt);
					*tp = f83_set1(rad ? e2 : 1.0);
					*rf = f83_set1(1.0);
					*psf = (1.0 - fr) / f8_abs(cnt);
					return 1;
				}
			}
			*wo = refl(n, wi, cni);
			*tp = f83_set1(1.0);
			*rf = f83_set1(1.0);
			*psf = fr / f8_abs(cni);
			return 1;
		}
		case 3 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 p = 2.0 * pi * x128p_f8(s);
			f8 cp = f8_cos(p), sp = f8_sin(p), ct = f8_sqrt(x128p_f8(s)), st = f8_op2u(ct);
			*wo = f833_trafo(f83_onb(f83_scl(-f8_sign(cni), n)), pol2cartzu0(cp, sp, ct, st));
			*tp = a;
			*rf = a;
			*psf = 1.0 / pi;
			return 1;
		}
		case 4 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 p = 2.0 * pi * x128p_f8(s);
			f8 cp = f8_cos(p), sp = f8_sin(p), ct = f8_sqrt(x128p_f8(s)), st = f8_op2u(ct);
			*wo = f833_trafo(f83_onb(f83_scl(f8_sign(cni), n)), pol2cartzu0(cp, sp, ct, st));
			*tp = a;
			*rf = a;
			*psf = 1.0 / pi;
			return 1;
		}
		case 5 : {
			f83 e = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f83 k = f83_set(f.ps[3], f.ps[4], f.ps[5]);
			f8 r = f.ps[6];
			f83 m = ggx_sample_vndf(n, r, f83_scl(f8_sign(cni), wi), s);
			f8 cnm = f83_dot(n, m), cmi = f83_dot(m, wi);
			*wo = refl(m, wi, cmi);
			f8 cno = f83_dot(n, *wo);
			if (cni * cno > 0.0) {
				return 0;
			}
			f83 fr = fresc3(f83_set(e.x / nr, e.y / nr, e.z / nr), k, cmi);
			f8 g1i = ggx_g1(r, cni), g1o = ggx_g1(r, cno);
			f8 d = ggx_d(r, cnm);
			*tp = f83_scl(g1o, fr);
			*rf = f83_scl(g1o, fr);
			*psf = g1i * d / f8_abs(4.0 * cni * cno);
			return 1;
		}
		case 6 : {
			return 0;
		}
		default : {
			return 0;
		}
	}
}
i4s bsdf_sample(sdf f, f8 nr, f8 nt, f83 n, f83 wi, i8u2 *s, f8 cni, i4s rad, f83 *wo, f83 *tp, f83 *rf) {
	switch (f.t) {
		case 0 : {
			*wo = wi;
			*tp = f83_set1(1.0);
			*rf = f83_set1(1.0);
			return 1;
		}
		case 1 : {
			f83 e = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f83 k = f83_set(f.ps[3], f.ps[4], f.ps[5]);
			f83 fr = fresc3(f83_set(e.x / nr, e.y / nr, e.z / nr), k, cni);
			*wo = refl(n, wi, cni);
			*tp = fr;
			*rf = fr;
			return 1;
		}
		case 2 : {
			f8 e = nr / nt;
			f8 e2 = f8_p2(e);
			f8 d = 1.0 - e2 * (1.0 - f8_p2(cni));
			f8 fr = 1.0;
			if (d > 0.0) {
				f8 cnt = f8_sign(cni) * f8_sqrt(d);
				fr = fresd(e, cni, cnt);
				if (x128p_f8(s) > fr) {
					*wo = refr(n, wi, e, cni, cnt);
					*tp = f83_set1(rad ? e2 : 1.0);
					*rf = f83_set1(1.0);
					return 1;
				}
			}
			*wo = refl(n, wi, cni);
			*tp = f83_set1(1.0);
			*rf = f83_set1(1.0);
			return 1;
		}
		case 3 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 p = 2.0 * pi * x128p_f8(s);
			f8 cp = f8_cos(p), sp = f8_sin(p), ct = f8_sqrt(x128p_f8(s)), st = f8_op2u(ct);
			*wo = f833_trafo(f83_onb(f83_scl(-f8_sign(cni), n)), pol2cartzu0(cp, sp, ct, st));
			*tp = a;
			*rf = a;
			return 1;
		}
		case 4 : {
			f83 a = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f8 p = 2.0 * pi * x128p_f8(s);
			f8 cp = f8_cos(p), sp = f8_sin(p), ct = f8_sqrt(x128p_f8(s)), st = f8_op2u(ct);
			*wo = f833_trafo(f83_onb(f83_scl(f8_sign(cni), n)), pol2cartzu0(cp, sp, ct, st));
			*tp = a;
			*rf = a;
			return 1;
		}
		case 5 : {
			f83 e = f83_set(f.ps[0], f.ps[1], f.ps[2]);
			f83 k = f83_set(f.ps[3], f.ps[4], f.ps[5]);
			f8 r = f.ps[6];
			f83 m = ggx_sample_vndf(n, r, f83_scl(f8_sign(cni), wi), s);
			f8 cnm = f83_dot(n, m), cmi = f83_dot(m, wi);
			*wo = refl(m, wi, cmi);
			f8 cno = f83_dot(n, *wo);
			if (cni * cno > 0.0) {
				return 0;
			}
			f83 fr = fresc3(f83_set(e.x / nr, e.y / nr, e.z / nr), k, cmi);
			f8 g1i = ggx_g1(r, cni), g1o = ggx_g1(r, cno);
			f8 d = ggx_d(r, cnm);
			*tp = f83_scl(g1o, fr);
			*rf = f83_scl(g1o, fr);
			return 1;
		}
		case 6 : {
			return 0;
		}
		default : {
			return 0;
		}
	}
}
void medr(f8 c, i4s2 mis, med *ms, med *m) {
	*m = c < 0.0 ? ms[mis.y] : ms[mis.x];
	return;
}
void medt(f8 c, i4s2 mis, med *ms, med *m) {
	*m = c < 0.0 ? ms[mis.x] : ms[mis.y];
	return;
}
void medrt(f8 c, i4s2 mis, med *ms, med *mr, med *mt) {
	medr(c, mis, ms, mr);
	medt(c, mis, ms, mt);
	return;
}
f8 cdf_pdf(f8 *cdf, i4s i) {
	return i == 0 ? cdf[0] : cdf[i] - cdf[i - 1];
}
i4s cdf_sample(i4s cdfn, f8 *cdf, i8u2 *s) {
	i4s a = 0, b = cdfn;
	f8 u = x128p_f8(s);
	while (a < b) {
		i4s m = (a + b) / 2;
		if (u < (m == 0 ? 0.0 : cdf[m - 1])) {
			b = m;
			continue;
		}
		if (u >= cdf[m]) {
			a = m + 1;
			continue;
		}
		return m;
	}
	return -1;
}
i4s bins(i4s isn, i4s *is, i4s x) {
	i4s a = 0, b = isn;
	while (a < b) {
		i4s m = (a + b) / 2;
		if (x < is[m]) {
			b = m;
			continue;
		}
		if (x > is[m]) {
			a = m + 1;
			continue;
		}
		return m;
	}
	return -1;
}
void make_cdf_env(i4s2 rsie, f83 *cse, i4s2 rsied, f82 rsiedf, i4s *cdfne, f8 **cdfe) {
	*cdfne = rsied.x * rsied.y;
	*cdfe = malloc(sizeof(f8) * *cdfne);
	f8 s = 0.0;
	for (i4s i = 0; i < *cdfne; ++i) {
		i4s2 pii = i4s2_set(i % rsied.x, i / rsied.x);
		f8 w = 0.0;
		w += (rsiedf.y * (f8_sin(pi * pii.y / rsiedf.y) - f8_sin(pi * (pii.y + 1) / rsiedf.y)) + pi * f8_cos(pi * pii.y / rsiedf.y)) / rsiedf.x * rgb2lum(cse[rsie.x * pii.y + pii.x]);
		w += (rsiedf.y * (f8_sin(pi * pii.y / rsiedf.y) - f8_sin(pi * (pii.y + 1) / rsiedf.y)) + pi * f8_cos(pi * pii.y / rsiedf.y)) / rsiedf.x * rgb2lum(cse[rsie.x * pii.y + pii.x + 1]);
		w += (rsiedf.y * (f8_sin(pi * (pii.y + 1) / rsiedf.y) - f8_sin(pi * pii.y / rsiedf.y)) - pi * f8_cos(pi * (pii.y + 1) / rsiedf.y)) / rsiedf.x * rgb2lum(cse[rsie.x * (pii.y + 1) + pii.x]);
		w += (rsiedf.y * (f8_sin(pi * (pii.y + 1) / rsiedf.y) - f8_sin(pi * pii.y / rsiedf.y)) - pi * f8_cos(pi * (pii.y + 1) / rsiedf.y)) / rsiedf.x * rgb2lum(cse[rsie.x * (pii.y + 1) + pii.x + 1]);
		s += w;
		(*cdfe)[i] = s;
	}
	for (i4s i = 0; i < rsied.x * rsied.y; ++i) {
		(*cdfe)[i] /= s;
	}
	return;
}
void make_cdf_scene(scene sc, i4s *cdfn, f8 **cdf, i4s **fis) {
	*cdfn = 0;
	for (i4s i = 0; i < sc.fsn; ++i) {
		if (sc.es[sc.fs[i].ei].t > 0) {
			++(*cdfn);
		}
	}
	*cdf = malloc(sizeof(i4s) * *cdfn);
	*fis = malloc(sizeof(i4s) * *cdfn);
	f8 s = 0.0;
	i4s c = 0;
	for (i4s i = 0; i < sc.fsn; ++i) {
		face f = sc.fs[i];
		sdf e = sc.es[f.ei];
		if (e.t > 0) {
			f83 a = f83_set(e.ps[0], e.ps[1], e.ps[2]);
			s += rgb2lum(a) * tri_area(get_vs(sc.vs, f.vis));
			(*cdf)[c] = s;
			(*fis)[c] = i;
			++c;
		}
	}
	for (i4s i = 0; i < *cdfn; ++i) {
		(*cdf)[i] /= s;
	}
	return;
}
f8 fdist(f8 a, f8 f) {
	return a * f / (a - f);
}
f83 trmit(f83 a, f8 t) {
	return f83_exp(f83_scl(-t, a));
}
i4s surf_kern(f8 cng, f8 cns) {
	return cng * cns > 0.0 && f8_abs(cng) > eps && f8_abs(cns) > eps;
}
i4s light_eval_00(sdf e, f8 cno, f83 *fs) {
	return edf_eval(e, cno, fs);
}
i4s light_eval_01(scene sc, i4s cdfn, f8 *cdf, i4s *fis, i4s id, f8 cno, f83 *fs, f8 *pab) {
	if (!edf_eval(sc.es[sc.fs[id].ei], cno, fs)) {
		return 0;
	}
	*pab = cdf_pdf(cdf, bins(cdfn, fis, id)) * tri_pdf(get_vs(sc.vs, sc.fs[id].vis));
	return 1;
}
i4s light_eval_02(scene sc, i4s cdfn, f8 *cdf, i4s *fis, i4s id, f8 cno, f83 *fs, f8 *pab, f8 *psb) {
	if (!edf_eval_b(sc.es[sc.fs[id].ei], cno, fs, psb)) {
		return 0;
	}
	*pab = cdf_pdf(cdf, bins(cdfn, fis, id)) * tri_pdf(get_vs(sc.vs, sc.fs[id].vis));
	return 1;
}
i4s light_eval_11(scene sc, i4s id, f8 cno, f83 *fs, f8 *psb) {
	if (!edf_eval_b(sc.es[sc.fs[id].ei], cno, fs, psb)) {
		return 0;
	}
	return 1;
}
void light_sample_0_f(face *fs, f83 *vs, i4s cdfn, f8 *cdf, i4s *fis, i8u2 *s, i4s *id, f83 *tp, f82 *uv, f8 *paf) {
	i4s i = cdf_sample(cdfn, cdf, s);
	*id = fis[i];
	f833 t = get_vs(vs, fs[*id].vis);
	*tp = f83_set1(1.0 / (cdf_pdf(cdf, i) * tri_pdf(t)));
	*uv = tri_sample(s);
	*paf = cdf_pdf(cdf, i) * tri_pdf(t);
}
i4s light_sample_1(scene sc, i4s cdfn, f8 *cdf, i4s *fis, i8u2 *s, f83 *tp, f83 *o, f83 *d, med *m) {
	i4s i = cdf_sample(cdfn, cdf, s);
	face f = sc.fs[fis[i]];
	f82 uv = tri_sample(s);
	f83 ns = tri_ns(get_vs(sc.ns, f.nis), uv);
	if (!edf_sample(sc.es[f.ei], ns, s, d, tp)) {
		return 0;
	}
	f833 t = get_vs(sc.vs, f.vis);
	f83 ng = tri_ng(t);
	f8 cngo = f83_dot(ng, *d), cnso = f83_dot(ns, *d);
	if (!surf_kern(cngo, cnso)) {
		return 0;
	}
	*o = tri_interp(t, uv);
	*tp = f83_scl(1.0 / (cdf_pdf(cdf, i) * tri_pdf(t)), *tp);
	medt(cngo, f.mis, sc.ms, m);
	return 1;
}
i4s light_sample_1_f(scene sc, i4s cdfn, f8 *cdf, i4s *fis, i8u2 *s, i4s *id, f83 *tp, f83 *tpd, f82 *uv, f83 *d, med *m, f8 *pt) {
	i4s i = cdf_sample(cdfn, cdf, s);
	*id = fis[i];
	face f = sc.fs[*id];
	*uv = tri_sample(s);
	f83 ns = tri_ns(get_vs(sc.ns, f.nis), *uv);
	if (!edf_sample_f(sc.es[f.ei], ns, s, d, tpd, pt)) {
		return 0;
	}
	f833 t = get_vs(sc.vs, f.vis);
	f83 ng = tri_ng(t);
	f8 cngo = f83_dot(ng, *d), cnso = f83_dot(ns, *d);
	if (!surf_kern(cngo, cnso)) {
		return 0;
	}
	*tp = f83_set1(1.0 / (cdf_pdf(cdf, i) * tri_pdf(t)));
	*tpd = f83_mul(*tp, *tpd);
	medt(cngo, f.mis, sc.ms, m);
	*pt *= fabs(cnso);
	return 1;
}
i4s eye_eval_00(i4s2 rsi, f82 rsif, f82 rsfh, f8 lr, f8 lf, f83 o, f83 od, i4s2 *pii, f83 *fs) {
	if (od.z + eps > 0.0) {
		return 0;
	}
	f83 pf3;
	if (f8_abs(lf + od.z) > eps) {
		pf3 = f83_sub(f83_scl(fdist(-od.z, lf) / od.z, od), o);
		pf3 = f83_add(o, f83_scl(1.0 / pf3.z, pf3));
	} else {
		pf3 = f83_add(o, f83_scl(1.0 / od.z, od));
	}
	f82 pf = f82_set(pf3.x, pf3.y);
	if (f8_abs(pf.x) > rsfh.x || f8_abs(pf.y) > rsfh.y) {
		return 0;
	}
	f83 d = f83_sub(od, o);
	f8 c4 = f8_p4(f8_abs(d.z)) / f8_p2(f83_len2(d));
	*pii = film2imgi(rsi, rsif, rsfh, pf);
	*fs = f83_set1(rsif.x * rsif.y / (disc_area(lr) * 4.0 * rsfh.x * rsfh.y * c4));
	return 1;
}
i4s eye_eval_01(i4s2 rsi, f82 rsif, f82 rsfh, f8 lr, f8 lf, f83 o, f83 od, i4s2 *pii, f83 *fs, f8 *pab) {
	if (od.z + eps > 0.0) {
		return 0;
	}
	f83 pf3;
	if (f8_abs(lf + od.z) > eps) {
		pf3 = f83_sub(f83_scl(fdist(-od.z, lf) / od.z, od), o);
		pf3 = f83_add(o, f83_scl(1.0 / pf3.z, pf3));
	} else {
		pf3 = f83_add(o, f83_scl(1.0 / od.z, od));
	}
	f82 pf = f82_set(pf3.x, pf3.y);
	if (f8_abs(pf.x) > rsfh.x || f8_abs(pf.y) > rsfh.y) {
		return 0;
	}
	f83 d = f83_sub(od, o);
	f8 c4 = f8_p4(f8_abs(d.z)) / f8_p2(f83_len2(d));
	*pii = film2imgi(rsi, rsif, rsfh, pf);
	*fs = f83_set1(rsif.x * rsif.y / (disc_area(lr) * 4.0 * rsfh.x * rsfh.y * c4));
	*pab = disc_pdf(lr);
	return 1;
}
i4s eye_eval_02(i4s2 rsi, f82 rsif, f82 rsfh, f8 lr, f8 lf, f83 o, f83 od, i4s2 *pii, f83 *fs, f8 *pab, f8 *psb) {
	if (od.z + eps > 0.0) {
		return 0;
	}
	f83 pf3;
	if (f8_abs(lf + od.z) > eps) {
		pf3 = f83_sub(f83_scl(fdist(-od.z, lf) / od.z, od), o);
		pf3 = f83_add(o, f83_scl(1.0 / pf3.z, pf3));
	} else {
		pf3 = f83_add(o, f83_scl(1.0 / od.z, od));
	}
	f82 pf = f82_set(pf3.x, pf3.y);
	if (f8_abs(pf.x) > rsfh.x || f8_abs(pf.y) > rsfh.y) {
		return 0;
	}
	f83 d = f83_sub(od, o);
	f8 c4 = f8_p4(f8_abs(d.z)) / f8_p2(f83_len2(d));
	*pii = film2imgi(rsi, rsif, rsfh, pf);
	*fs = f83_set1(rsif.x * rsif.y / (disc_area(lr) * 4.0 * rsfh.x * rsfh.y * c4));
	*pab = disc_pdf(lr);
	*psb = 1.0 / (4.0 * rsfh.x * rsfh.y * c4);
	return 1;
}
i4s eye_eval_11(i4s2 rsi, f82 rsif, f82 rsfh, f8 lr, f8 lf, f83 o, f83 od, i4s2 *pii, f83 *fs, f8 *psb) {
	if (od.z + eps > 0.0) {
		return 0;
	}
	f83 pf3;
	if (f8_abs(lf + od.z) > eps) {
		pf3 = f83_sub(f83_scl(fdist(-od.z, lf) / od.z, od), o);
		pf3 = f83_add(o, f83_scl(1.0 / pf3.z, pf3));
	} else {
		pf3 = f83_add(o, f83_scl(1.0 / od.z, od));
	}
	f82 pf = f82_set(pf3.x, pf3.y);
	if (f8_abs(pf.x) > rsfh.x || f8_abs(pf.y) > rsfh.y) {
		return 0;
	}
	f83 d = f83_sub(od, o);
	f8 c4 = f8_p4(f8_abs(d.z)) / f8_p2(f83_len2(d));
	*pii = film2imgi(rsi, rsif, rsfh, pf);
	*fs = f83_set1(rsif.x * rsif.y / (disc_area(lr) * 4.0 * rsfh.x * rsfh.y * c4));
	*psb = 1.0 / (4.0 * rsfh.x * rsfh.y * c4);
	return 1;
}
void eye_sample_0_f(f8 lr, i8u2 *s, f83 *tp, f83 *o, f8 *paf) {
	*paf = disc_pdf(lr);
	*o = disc_sample(lr, s);
	*tp = f83_set1(1.0 / *paf);
}
void eye_sample_1(f82 rsif, f82 rsfh, f8 lr, f8 lf, i8u2 *s, f82 pif, f83 *tp, f83 *o, f83 *d) {
	*tp = f83_set1(rsif.x * rsif.y);
	*o = disc_sample(lr, s);
	*d = f83_norm(f83_sub(f83_scl(-fdist(1.0, lf), imgf2film(rsif, rsfh, pif)), *o));
}
void eye_sample_1_f(f82 rsif, f82 rsfh, f8 lr, f8 lf, i8u2 *s, f82 pif, f83 *tp, f83 *tpd, f83 *o, f83 *d, f8 *pt) {
	*tp = f83_set1(1.0 / disc_pdf(lr));
	*tpd = f83_set1(rsif.x * rsif.y);
	*o = disc_sample(lr, s);
	*d = f83_norm(f83_sub(f83_scl(-fdist(1.0, lf), imgf2film(rsif, rsfh, pif)), *o));
	*pt = 1.0 / (4.0 * rsfh.x * rsfh.y * f8_p3(f8_abs(d->z)));
}
typedef struct vert vert;
struct vert {
	f83 x, ng, ns, wi, tp, tr;
	f8 t, cnsi, cngo, pf, pb;
	i4s id, fl;
};
void eye_sample_vert(f8 lr, f8 lf, f82 rsif, f82 rfh, i8u2 *s, f82 pif, f83 *tp, f83 *o, f83 *d, f8 *pt, vert *v) {
	*tp = f83_set1(rsif.x * rsif.y);
	*o = disc_sample(lr, s);
	*d = f83_norm(f83_sub(f83_scl(-fdist(1.0, lf), imgf2film(rsif, rfh, pif)), *o));
	*pt = 1.0 / (4.0 * rfh.x * rfh.y * f8_p3(f8_abs(d->z)));
	v->x = *o;
	v->ng = disc_ng();
	v->ns = disc_ng();
	v->tp = f83_set1(1.0 / disc_pdf(lr));
	v->cngo = d->z;
	v->pf = disc_pdf(lr);
	return;
}
i4s light_sample_vert(scene sc, i4s cdfn, f8 *cdf, i4s *fis, i8u2 *s, f83 *tp, f83 *o, f83 *d, f8 *pt, vert *v) {
	i4s i = cdf_sample(cdfn, cdf, s);
	v->id = fis[i];
	face f = sc.fs[v->id];
	f82 uv = tri_sample(s);
	v->ns = tri_ns(get_vs(sc.ns, f.nis), uv);
	f83 tpd;
	if (!edf_sample_f(sc.es[f.ei], v->ns, s, d, &tpd, pt)) {
		return 0;
	}
	f833 t = get_vs(sc.vs, f.vis);
	v->ng = tri_ng(t);
	v->cngo = f83_dot(v->ng, *d);
	f8 cnso = f83_dot(v->ns, *d);
	if (!surf_kern(v->cngo, cnso)) {
		return 0;
	}
	v->pf = cdf_pdf(cdf, i) * tri_pdf(t);
	v->tp = f83_set1(1.0 / v->pf);
	*tp = f83_mul(v->tp, tpd);
	*o = tri_interp(t, uv);
	v->x = *o;
	*pt *= fabs(cnso);
	return 1;
}
f82 pol2imgf(f82 rsiedf, f82 p) {
	return f82_set(rsiedf.x * (1.0 - p.x / (2.0 * pi)), rsiedf.y * (1.0 - p.y / pi));
}
f82 imgf2pol(f82 rsiedf, f82 p) {
	return f82_set(2.0 * pi * (1.0 - p.x / rsiedf.x), pi * (1.0 - p.y / rsiedf.y));
}
void env_eval(i4s2 rsie, f83 *cse, f833 m, f82 rsiedf, f83 d, f83 *fs) {
	d = f833_trafo(m, d);
	f82 pif = pol2imgf(rsiedf, cart2polzu(d));
	*fs = image_interp(rsie, cse, pif);
}
void env_eval_b(i4s2 rsie, f83 *cse, f833 m, i4s2 rsied, f82 rsiedf, f8 *cdfe, f83 d, f83 *fs, f8 *psb) {
	d = f833_trafo(m, d);
	f82 pif = pol2imgf(rsiedf, cart2polzu(d));
	i4s2 pii = i4s2_set(pif.x, pif.y);
	*fs = image_interp(rsie, cse, pif);
	*psb = rsiedf.x * rsiedf.y * cdf_pdf(cdfe, rsied.x * pii.y + pii.x) / (2.0 * f8_p2(pi) * f8_max(f8_op2us(d.z), eps));
}
void env_sample_f(i4s2 rsie, f83 *cse, f833 mi, i4s2 rsied, f82 rsiedf, i4s cdfne, f8 *cdfe, i8u2 *s, f83 *tp, f83 *d, f8 *psf) {
	i4s i = cdf_sample(cdfne, cdfe, s);
	f82 pif = f82_set(i % rsied.x + x128p_f8(s), i / rsied.x + x128p_f8(s));
	*d = pol2cartzu1(imgf2pol(rsiedf, pif));
	*psf = rsiedf.x * rsiedf.y * cdf_pdf(cdfe, i) / (2.0 * f8_p2(pi) * f8_max(f8_op2us(d->z), eps));
	*tp = f83_scl(1.0 / *psf, image_interp(rsie, cse, pif));
	*d = f833_trafo(mi, *d);
}
i4s trav_bpt(scene sc, i4s en, bvh_node *bns, i4s *st, i8u2 *s, i4s rad, f83 o, f83 d, f83 tp, f8 pt, vert *vs) {
	f83 rb;
	f8 psb;
	for (i4s i = 0; i < en; ++i) {
		f8 t;
		f82 uv;
		i4s id;
		if (!bvh_travrt(sc.vs, sc.fs, bns, st, o, d, &t, &uv, &id)) {
			vert v;
			v.wi = d;
			v.pf = pt;
			v.fl = 0;
			v.fl |= 1 << 2;
			vs[i + 1] = v;
			return i + 2;
		}
		o = f83_add(o, f83_scl(t, d));
		face f = sc.fs[id];
		f83 ng = tri_ng(get_vs(sc.vs, f.vis)), ns = tri_ns(get_vs(sc.ns, f.nis), uv);
		f8 cngi = f83_dot(ng, d), cnsi = f83_dot(ns, d);
		if (!surf_kern(cngi, cnsi)) {
			return i + 1;
		}
		med mr, mt;
		medrt(cngi, f.mis, sc.ms, &mr, &mt);
		f83 tr = trmit(mr.a, t);
		tp = f83_mul(tp, f83_scl(f8_abs(cnsi / cngi), tr));
		if (i > 0) {
			vert *vm2 = vs + i - 1, *vm1 = vs + i;
			vm2->pb = rgb2lum(f83_mul(tr, rb)) * psb * f8_abs(vm1->cnsi * vm2->cngo) / f8_p2(vm1->t);
		}
		sdf e = sc.es[f.ei], b = sc.bs[f.bi];
		vert v;
		v.x = o;
		v.ng = ng;
		v.ns = ns;
		v.wi = d;
		v.tp = tp;
		v.tr = tr;
		v.t = t;
		v.cnsi = cnsi;
		v.pf = pt * f8_abs(cngi) / f8_p2(t);
		v.id = id;
		v.fl = 0;
		v.fl |= bsdf_con(b) << 0;
		v.fl |= edf_con(e) << 1;
		f83 dd, tpd, rf;
		f8 psf;
		if (!bsdf_sample_fb(b, mr.n, mt.n, ns, d, s, cnsi, rad, &dd, &tpd, &rf, &rb, &psf, &psb)) {
			vs[i + 1] = v;
			return i + 2;
		}
		f8 cngo = f83_dot(ng, dd), cnso = f83_dot(ns, dd);
		if (!surf_kern(cngo, cnso)) {
			vs[i + 1] = v;
			return i + 2;
		}
		v.cngo = cngo;
		vs[i + 1] = v;
		f8 pr = rgb2lum(f83_mul(tr, rf));
		if (x128p_f8(s) >= pr) {
			return i + 2;
		}
		tp = f83_mul(tp, f83_scl(1.0 / pr, tpd));
		pt = psf * f8_abs(cnso) * pr;
		d = dd;
	}
	return en + 1;
}
f83 trav_pt(scene sc, i4s en, bvh_node *bns, i4s *st, i8u2 *s, f83 o, f83 d, f83 tp) {
	f83 c = f83_set1(0.0);
	for (i4s i = 0; i < en; ++i) {
		f8 t;
		f82 uv;
		i4s id;
		if (!bvh_travrt(sc.vs, sc.fs, bns, st, o, d, &t, &uv, &id)) {
			return c;
		}
		o = f83_add(o, f83_scl(t, d));
		face f = sc.fs[id];
		f83 ng = tri_ng(get_vs(sc.vs, f.vis)), ns = tri_ns(get_vs(sc.ns, f.nis), uv);
		f8 cngi = f83_dot(ng, d), cnsi = f83_dot(ns, d);
		if (!surf_kern(cngi, cnsi)) {
			return c;
		}
		med mr, mt;
		medrt(cngi, f.mis, sc.ms, &mr, &mt);
		f83 tr = trmit(mr.a, t);
		tp = f83_mul(tp, f83_scl(f8_abs(cnsi / cngi), tr));
		sdf e = sc.es[f.ei];
		if (edf_con(e)) {
			f83 fs;
			if (edf_eval(e, -cnsi, &fs)) {
				c = f83_add(c, f83_mul(tp, fs));
			}
		}
		sdf b = sc.bs[f.bi];
		f83 dd, tpd, rf;
		if (!bsdf_sample(b, mr.n, mt.n, ns, d, s, cnsi, 1, &dd, &tpd, &rf)) {
			return c;
		}
		f8 cngo = f83_dot(ng, dd), cnso = f83_dot(ns, dd);
		if (!surf_kern(cngo, cnso)) {
			return c;
		}
		f8 pr = rgb2lum(f83_mul(tr, rf));
		if (x128p_f8(s) >= pr) {
			return c;
		}
		tp = f83_mul(tp, f83_scl(1.0 / pr, tpd));
		d = dd;
	}
	return c;
}
f83 trav_ptd(scene sc, i4s en, i4s cdfn, f8 *cdf, i4s *fis, bvh_node *bns, i4s *st, i8u2 *s, f83 o, f83 d, f83 tp, f8 pt) {
	i4s con = 1;
	f83 c = f83_set1(0.0);
	for (i4s i = 0; i < en; ++i) {
		f8 t;
		f82 uv;
		i4s id;
		if (!bvh_travrt(sc.vs, sc.fs, bns, st, o, d, &t, &uv, &id)) {
			return c;
		}
		o = f83_add(o, f83_scl(t, d));
		face f = sc.fs[id];
		f83 ng = tri_ng(get_vs(sc.vs, f.vis)), ns = tri_ns(get_vs(sc.ns, f.nis), uv);
		f8 cngi = f83_dot(ng, d), cnsi = f83_dot(ns, d);
		if (!surf_kern(cngi, cnsi)) {
			return c;
		}
		med mr, mt;
		medrt(cngi, f.mis, sc.ms, &mr, &mt);
		f83 tr = trmit(mr.a, t);
		tp = f83_mul(tp, f83_scl(f8_abs(cnsi / cngi), tr));
		sdf e = sc.es[f.ei];
		if (edf_con(e)) {
			f83 fsl;
			f8 pb;
			if (light_eval_01(sc, cdfn, cdf, fis, id, -cnsi, &fsl, &pb)) {
				if (con) {
					f8 pf = pt * f8_abs(cngi) / f8_p2(t);
					f8 w = 1.0 / (1.0 + f8_p2(pb / pf));
					c = f83_add(c, f83_scl(w, f83_mul(tp, fsl)));
				} else {
					c = f83_add(c, f83_mul(tp, fsl));
				}
			}
		}
		sdf b = sc.bs[f.bi];
		while (bsdf_con(b) && i + 1 < en) {
			i4s idl;
			f83 tpl;
			f82 uvl;
			f8 pf;
			light_sample_0_f(sc.fs, sc.vs, cdfn, cdf, fis, s, &idl, &tpl, &uvl, &pf);
			face fl = sc.fs[idl];
			f83 od = tri_interp(get_vs(sc.vs, fl.vis), uvl);
			f83 dd = f83_sub(od, o);
			f8 td = f83_len(dd);
			dd = f83_scl(1.0 / td, dd);
			f8 cngoe = f83_dot(ng, dd), cnsoe = f83_dot(ns, dd);
			if (!surf_kern(cngoe, cnsoe)) {
				break;
			}
			f83 fse, rfe;
			f8 psfe;
			if (!bsdf_eval_f(b, mr.n, mt.n, ns, d, dd, cnsi, cnsoe, 1, &fse, &rfe, &psfe)) {
				break;
			}
			f83 ngl = tri_ng(get_vs(sc.vs, fl.vis)), nsl = tri_ns(get_vs(sc.ns, fl.nis), uvl);
			f8 cngol = -f83_dot(ngl, dd), cnsol = -f83_dot(nsl, dd);
			if (!surf_kern(cngol, cnsol)) {
				break;
			}
			f83 fsl;
			if (!edf_eval(sc.es[fl.ei], cnsol, &fsl)) {
				break;
			}
			if (!bvh_visrt(sc.vs, sc.fs, bns, st, o, dd, td)) {
				break;
			}
			med m;
			medt(cngol, fl.mis, sc.ms, &m);
			f83 trd = trmit(m.a, td);
			f8 prb = rgb2lum(f83_mul(tr, rfe));
			f8 pab = psfe * f8_abs(cnsoe * cngol) / f8_p2(td);
			f8 pb = prb * pab;
			f8 w = 1.0 / (1.0 + f8_p2(pb / pf));
			f8 g = f8_abs(cnsoe * cnsol) / f8_p2(td);
			c = f83_add(c, f83_scl(w * g, f83_mul(trd, f83_mul(f83_mul(tp, fse), f83_mul(tpl, fsl)))));
			break;
		}
		f83 dd, tpd, rf;
		f8 psf;
		if (!bsdf_sample_f(b, mr.n, mt.n, ns, d, s, cnsi, 1, &dd, &tpd, &rf, &psf)) {
			return c;
		}
		f8 cngo = f83_dot(ng, dd), cnso = f83_dot(ns, dd);
		if (!surf_kern(cngo, cnso)) {
			return c;
		}
		f8 pr = rgb2lum(f83_mul(tr, rf));
		if (x128p_f8(s) >= pr) {
			return c;
		}
		tp = f83_mul(tp, f83_scl(1.0 / pr, tpd));
		d = dd;
		pt = pr * psf * f8_abs(cnso);
		con = bsdf_con(b);
	}
	return c;
}
void trav_lt(scene sc, i4s en, i4s2 rsi, f83 *cs, f82 rsif, f82 rsfh, f8 lr, f8 lf, f8 ln, bvh_node *bns, i4s *st, i8u2 *s, f83 o, f83 d, f83 tp, med m) {
	for (i4s i = 0; i < en; ++i) {
		f8 t;
		f82 uv;
		i4s id;
		i4s h = bvh_travrt(sc.vs, sc.fs, bns, st, o, d, &t, &uv, &id);
		f8 td = t;
		if (disc_isect(lr, o, d, &td)) {
			f83 od = f83_add(o, f83_scl(td, d));
			i4s2 pii;
			f83 fs;
			if (!eye_eval_00(rsi, rsif, rsfh, lr, lf, od, o, &pii, &fs)) {
				return;
			}
			f83 tr = trmit(m.a, td);
			cs[rsi.x * pii.y + pii.x] = f83_add(cs[rsi.x * pii.y + pii.x], f83_scl(f8_abs(f83_dot(disc_ng(), d)) / ln, f83_mul(tr, f83_mul(tp, fs))));
		}
		if (!h) {
			return;
		}
		o = f83_add(o, f83_scl(t, d));
		face f = sc.fs[id];
		f83 ng = tri_ng(get_vs(sc.vs, f.vis)), ns = tri_ns(get_vs(sc.ns, f.nis), uv);
		f8 cngi = f83_dot(ng, d), cnsi = f83_dot(ns, d);
		if (!surf_kern(cngi, cnsi)) {
			return;
		}
		med mr, mt;
		medrt(cngi, f.mis, sc.ms, &mr, &mt);
		f83 tr = trmit(mr.a, t);
		tp = f83_mul(tp, f83_scl(f8_abs(cnsi / cngi), tr));
		sdf b = sc.bs[f.bi];
		f83 dd, tpd, rf;
		if (!bsdf_sample(b, mr.n, mt.n, ns, d, s, cnsi, 0, &dd, &tpd, &rf)) {
			return;
		}
		f8 cngo = f83_dot(ng, dd), cnso = f83_dot(ns, dd);
		if (!surf_kern(cngo, cnso)) {
			return;
		}
		f8 pr = rgb2lum(f83_mul(tr, rf));
		if (x128p_f8(s) >= pr) {
			return;
		}
		tp = f83_mul(tp, f83_scl(1.0 / pr, tpd));
		d = dd;
		medt(cngo, f.mis, sc.ms, &m);
	}
	return;
}
void trav_ltd(scene sc, i4s en, i4s2 rsi, f83 *cs, f82 rsif, f82 rsfh, f8 lr, f8 lf, f8 ln, bvh_node *bns, i4s *st, i8u2 *s, f83 o, f83 d, f83 tp, med m, f8 pt) {
	i4s con = 1;
	for (i4s i = 0; i < en; ++i) {
		f8 t;
		f82 uv;
		i4s id;
		i4s h = bvh_travrt(sc.vs, sc.fs, bns, st, o, d, &t, &uv, &id);
		f8 td = t;
		if (disc_isect(lr, o, d, &td)) {
			f83 od = f83_add(o, f83_scl(td, d));
			i4s2 pii;
			f83 fs;
			f8 pb;
			if (!eye_eval_01(rsi, rsif, rsfh, lr, lf, od, o, &pii, &fs, &pb)) {
				return;
			}
			f83 nge = disc_ng();
			f8 cngie = f83_dot(nge, d);
			f83 tr = trmit(m.a, td);
			if (con) {
				f8 pf = pt * f8_abs(cngie) / f8_p2(t);
				f8 w = 1.0 / (1.0 + f8_p2(pb / pf));
				cs[rsi.x * pii.y + pii.x] = f83_add(cs[rsi.x * pii.y + pii.x], f83_scl(w * f8_abs(cngie) / ln, f83_mul(tr, f83_mul(tp, fs))));
			} else {
				cs[rsi.x * pii.y + pii.x] = f83_add(cs[rsi.x * pii.y + pii.x], f83_scl(f8_abs(cngie) / ln, f83_mul(tr, f83_mul(tp, fs))));
			}
		}
		if (!h) {
			return;
		}
		o = f83_add(o, f83_scl(t, d));
		face f = sc.fs[id];
		f83 ng = tri_ng(get_vs(sc.vs, f.vis)), ns = tri_ns(get_vs(sc.ns, f.nis), uv);
		f8 cngi = f83_dot(ng, d), cnsi = f83_dot(ns, d);
		if (!surf_kern(cngi, cnsi)) {
			return;
		}
		med mr, mt;
		medrt(cngi, f.mis, sc.ms, &mr, &mt);
		f83 tr = trmit(mr.a, t);
		tp = f83_mul(tp, f83_scl(f8_abs(cnsi / cngi), tr));
		sdf b = sc.bs[f.bi];
		while (bsdf_con(b) && i + 1 < en) {
			f83 tpe;
			f83 od;
			f8 pf;
			eye_sample_0_f(lr, s, &tpe, &od, &pf);
			f83 dd = f83_sub(o, od);
			f8 td = f83_len(dd);
			dd = f83_scl(1.0 / td, dd);
			f83 nge = disc_ng(), nse = disc_ng();
			f8 cngoe = f83_dot(nge, dd), cnsoe = f83_dot(nse, dd);
			i4s2 pii;
			f83 fse;
			if (!eye_eval_00(rsi, rsif, rsfh, lr, lf, od, o, &pii, &fse)) {
				break;
			}
			f8 cngol = -f83_dot(ng, dd), cnsol = -f83_dot(ns, dd);
			if (!surf_kern(cngol, cnsol)) {
				break;
			}
			f83 fsl, rbl;
			f8 psbe;
			if (!bsdf_eval_f(b, mr.n, mt.n, ns, d, f83_scl(-1.0, dd), cnsi, cnsol, 0, &fsl, &rbl, &psbe)) {
				break;
			}
			if (!bvh_visrt(sc.vs, sc.fs, bns, st, od, dd, td)) {
				break;
			}
			med md;
			medt(cngol, f.mis, sc.ms, &md);
			f83 trd = trmit(m.a, td);
			f8 pb = psbe * f8_abs(cnsol * cngoe) / f8_p2(td);
			f8 w = 1.0 / (1.0 + f8_p2(pb / pf));
			f8 g = f8_abs(cnsoe * cnsol) / f8_p2(td);
			cs[pii.y * rsi.x + pii.x] = f83_add(cs[pii.y * rsi.x + pii.x], f83_scl(w * g / ln, f83_mul(trd, f83_mul(f83_mul(tpe, fse), f83_mul(tp, fsl)))));
			break;
		}
		f83 dd, tpd, rf;
		f8 psf;
		if (!bsdf_sample_f(b, mr.n, mt.n, ns, d, s, cnsi, 0, &dd, &tpd, &rf, &psf)) {
			return;
		}
		f8 cngo = f83_dot(ng, dd), cnso = f83_dot(ns, dd);
		if (!surf_kern(cngo, cnso)) {
			return;
		}
		f8 pr = rgb2lum(f83_mul(tr, rf));
		if (x128p_f8(s) >= pr) {
			return;
		}
		tp = f83_mul(tp, f83_scl(1.0 / pr, tpd));
		d = dd;
		medt(cngo, f.mis, sc.ms, &m);
		pt = pr * psf * f8_abs(cnso);
		con = bsdf_con(b);
	}
	return;
}
f83 trav_pte(scene sc, i4s en, i4s2 rsie, f83 *cse, f833 m, f82 rsiedf, bvh_node *bns, i4s *st, i8u2 *s, i4s rad, f83 o, f83 d, f83 tp) {
	f83 c = f83_set1(0.0);
	for (i4s i = 0; i < en; ++i) {
		f8 t;
		f82 uv;
		i4s id;
		if (!bvh_travrt(sc.vs, sc.fs, bns, st, o, d, &t, &uv, &id)) {
			f83 fs;
			env_eval(rsie, cse, m, rsiedf, d, &fs);
			c = f83_add(c, f83_mul(tp, fs));
			return c;
		}
		o = f83_add(o, f83_scl(t, d));
		face f = sc.fs[id];
		f83 ng = tri_ng(get_vs(sc.vs, f.vis)), ns = tri_ns(get_vs(sc.ns, f.nis), uv);
		f8 cngi = f83_dot(ng, d), cnsi = f83_dot(ns, d);
		if (!surf_kern(cngi, cnsi)) {
			return c;
		}
		med mr, mt;
		medrt(cngi, f.mis, sc.ms, &mr, &mt);
		f83 tr = trmit(mr.a, t);
		tp = f83_mul(tp, f83_scl(f8_abs(cnsi / cngi), tr));
		sdf b = sc.bs[f.bi];
		f83 dd, tpd, rf;
		if (!bsdf_sample(b, mr.n, mt.n, ns, d, s, cnsi, rad, &dd, &tpd, &rf)) {
			return c;
		}
		f8 cngo = f83_dot(ng, dd), cnso = f83_dot(ns, dd);
		if (!surf_kern(cngo, cnso)) {
			return c;
		}
		f8 pr = rgb2lum(f83_mul(tr, rf));
		if (x128p_f8(s) >= pr) {
			return c;
		}
		tp = f83_mul(tp, f83_scl(1.0 / pr, tpd));
		d = dd;
	}
	return c;
}
f83 trav_pted(scene sc, i4s en, i4s2 rsie, f83 *cse, f833 m, f833 mi, i4s2 rsied, f82 rsiedf, i4s cdfne, f8 *cdfe, bvh_node *bns, i4s *st, i8u2 *s, i4s rad, f83 o, f83 d, f83 tp, f8 pt) {
	i4s con = 0;
	f83 c = f83_set1(0.0);
	for (i4s i = 0; i < en; ++i) {
		f8 t;
		f82 uv;
		i4s id;
		if (!bvh_travrt(sc.vs, sc.fs, bns, st, o, d, &t, &uv, &id)) {
			f83 fsl;
			f8 pb;
			env_eval_b(rsie, cse, m, rsied, rsiedf, cdfe, d, &fsl, &pb);
			if (con) {
				f8 pf = pt;
				f8 w = 1.0 / (1.0 + f8_p2(pb / pf));
				c = f83_add(c, f83_scl(w, f83_mul(tp, fsl)));
			} else {
				c = f83_add(c, f83_mul(tp, fsl));
			}
			return c;
		}
		o = f83_add(o, f83_scl(t, d));
		face f = sc.fs[id];
		f83 ng = tri_ng(get_vs(sc.vs, f.vis)), ns = tri_ns(get_vs(sc.ns, f.nis), uv);
		f8 cngi = f83_dot(ng, d), cnsi = f83_dot(ns, d);
		if (!surf_kern(cngi, cnsi)) {
			return c;
		}
		med mr, mt;
		medrt(cngi, f.mis, sc.ms, &mr, &mt);
		f83 tr = trmit(mr.a, t);
		tp = f83_mul(tp, f83_scl(f8_abs(cnsi / cngi), tr));
		sdf b = sc.bs[f.bi];
		while (bsdf_con(b) && i + 1 < en) {
			f83 tpl, dd;
			f8 pf;
			env_sample_f(rsie, cse, mi, rsied, rsiedf, cdfne, cdfe, s, &tpl, &dd, &pf);
			f8 cngoe = f83_dot(ng, dd), cnsoe = f83_dot(ns, dd);
			if (!surf_kern(cngoe, cnsoe)) {
				break;
			}
			f83 fse, rfe;
			f8 psfe;
			if (!bsdf_eval_f(b, mr.n, mt.n, ns, d, dd, cnsi, cnsoe, 1, &fse, &rfe, &psfe)) {
				break;
			}
			if (!bvh_visrt(sc.vs, sc.fs, bns, st, o, dd, inf)) {
				break;
			}
			f8 pb = rgb2lum(f83_mul(tr, rfe)) * psfe * f8_abs(cnsoe);
			f8 w = 1.0 / (1.0 + f8_p2(pb / pf));
			f8 g = f8_abs(cnsoe);
			c = f83_add(c, f83_scl(w * g, f83_mul(f83_mul(tp, fse), tpl)));
			break;
		}
		f83 dd, tpd, rf;
		f8 psf;
		if (!bsdf_sample_f(b, mr.n, mt.n, ns, d, s, cnsi, rad, &dd, &tpd, &rf, &psf)) {
			return c;
		}
		f8 cngo = f83_dot(ng, dd), cnso = f83_dot(ns, dd);
		if (!surf_kern(cngo, cnso)) {
			return c;
		}
		f8 pr = rgb2lum(f83_mul(tr, rf));
		if (x128p_f8(s) >= pr) {
			return c;
		}
		tp = f83_mul(tp, f83_scl(1.0 / pr, tpd));
		d = dd;
		pt = pr * psf * f8_abs(cnso);
		con = bsdf_con(b);
	}
	return c;
}
void render_pt(conf co) {
	f82 rsif = f82_set(co.rsi.x, co.rsi.y);
	f8 ln = rsif.x * rsif.y * co.sn;
	scene sc = scene_read(co.ps);
	bvh_node *bns = malloc(sizeof(bvh_node) * co.bnsna);
	i4s bnsn = bvh_makert(&sc, co.bnsna, bns);
	printf("bvh usage:%d/%d\n", bnsn, co.bnsna);
	omp_set_num_threads(co.tn);
	f83 *css[omp_get_max_threads()];
	f8 t = timer();
	#pragma omp parallel
	{
		css[omp_get_thread_num()] = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
		f83 *cs = css[omp_get_thread_num()];
		i4s *st = malloc(sizeof(i4s) * bnsn);
		i8u2 s = x128p_seed(co.s);
		s = x128p_jumpn(s, omp_get_thread_num());
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			f83 c = f83_set1(0.0);
			for (i4s j = 0; j < co.sn; ++j) {
				f82 pif = f82_set(i % co.rsi.x + x128p_f8(&s), i / co.rsi.x + x128p_f8(&s));
				f83 tp, o, d;
				eye_sample_1(rsif, co.rsfh, co.lr, co.lf, &s, pif, &tp, &o, &d);
				c = f83_add(c, trav_pt(sc, co.en, bns, st, &s, o, d, tp));
			}
			cs[i] = f83_scl(1.0 / ln, c);
		}
		free(st);
	}
	t = timer() - t;
	printf("time:%f\n", t);
	free(bns);
	scene_free(sc);
	f83 *cs = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
	for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
		cs[i] = f83_set1(0.0);
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		for (i4s j = 0; j < co.rsi.x * co.rsi.y; ++j) {
			cs[j] = f83_add(cs[j], f83_scl(1.0 / omp_get_max_threads(), css[i][j]));
		}
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		free(css[i]);
	}
	image_write(co.pr, co.rsi, cs);
	free(cs);
	return;
}
void render_ptd(conf co) {
	f82 rsif = f82_set(co.rsi.x, co.rsi.y);
	f8 ln = rsif.x * rsif.y * co.sn;
	scene sc = scene_read(co.ps);
	bvh_node *bns = malloc(sizeof(bvh_node) * co.bnsna);
	i4s bnsn = bvh_makert(&sc, co.bnsna, bns);
	printf("bvh usage:%d/%d\n", bnsn, co.bnsna);
	i4s cdfn;
	f8 *cdf;
	i4s *fis;
	make_cdf_scene(sc, &cdfn, &cdf, &fis);
	omp_set_num_threads(co.tn);
	f83 *css[omp_get_max_threads()];
	f8 t = timer();
	#pragma omp parallel
	{
		css[omp_get_thread_num()] = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
		f83 *cs = css[omp_get_thread_num()];
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			cs[i] = f83_set1(0.0);
		}
		i4s *st = malloc(sizeof(i4s) * bnsn);
		i8u2 s = x128p_seed(co.s);
		s = x128p_jumpn(s, omp_get_thread_num());
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			f83 c = f83_set1(0.0);
			for (i4s j = 0; j < co.sn; ++j) {
				f82 pif = f82_set(i % co.rsi.x + x128p_f8(&s), i / co.rsi.x + x128p_f8(&s));
				f83 tpe, tp, o, d;
				f8 pt;
				eye_sample_1_f(rsif, co.rsfh, co.lr, co.lf, &s, pif, &tpe, &tp, &o, &d, &pt);
				while (0 < co.en) {
					i4s idl;
					f83 tpl;
					f82 uvl;
					f8 pf;
					light_sample_0_f(sc.fs, sc.vs, cdfn, cdf, fis, &s, &idl, &tpl, &uvl, &pf);
					face fl = sc.fs[idl];
					f83 od = tri_interp(get_vs(sc.vs, fl.vis), uvl);
					f83 dd = f83_sub(od, o);
					f8 td = f83_len(dd);
					dd = f83_scl(1.0 / td, dd);
					f83 nge = disc_ng(), nse = disc_ng();
					f8 cngoe = f83_dot(nge, dd), cnsoe = f83_dot(nse, dd);
					i4s2 pii;
					f83 fse;
					f8 psfl;
					if (!eye_eval_11(co.rsi, rsif, co.rsfh, co.lr, co.lf, o, od, &pii, &fse, &psfl)) {
						break;
					}
					f83 ngl = tri_ng(get_vs(sc.vs, fl.vis)), nsl = tri_ns(get_vs(sc.ns, fl.nis), uvl);
					f8 cngol = -f83_dot(ngl, dd), cnsol = -f83_dot(nsl, dd);
					if (!surf_kern(cngol, cnsol)) {
						break;
					}
					f83 fsl;
					if (!edf_eval(sc.es[fl.ei], cnsol, &fsl)) {
						break;
					}
					if (!bvh_visrt(sc.vs, sc.fs, bns, st, o, dd, td)) {
						break;
					}
					med m;
					medt(cngol, fl.mis, sc.ms, &m);
					f83 trd = trmit(m.a, td);
					f8 pb = psfl * f8_abs(cnsoe * cngol) / f8_p2(td);
					f8 w = 1.0 / (1.0 + f8_p2(pb / pf));
					f8 g = f8_abs(cnsoe * cnsol) / f8_p2(td);
					cs[pii.y * co.rsi.x + pii.x] = f83_add(cs[pii.y * co.rsi.x + pii.x], f83_scl(w * g / ln, f83_mul(trd, f83_mul(f83_mul(tpe, fse), f83_mul(tpl, fsl)))));
					break;
				}
				c = f83_add(c, trav_ptd(sc, co.en, cdfn, cdf, fis, bns, st, &s, o, d, tp, pt));
			}
			cs[i] = f83_add(cs[i], f83_scl(1.0 / ln, c));
		}
		free(st);
	}
	t = timer() - t;
	printf("time:%f\n", t);
	free(cdf);
	free(fis);
	free(bns);
	scene_free(sc);
	f83 *cs = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
	for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
		cs[i] = f83_set1(0.0);
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		for (i4s j = 0; j < co.rsi.x * co.rsi.y; ++j) {
			cs[j] = f83_add(cs[j], f83_scl(1.0 / omp_get_max_threads(), css[i][j]));
		}
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		free(css[i]);
	}
	image_write(co.pr, co.rsi, cs);
	free(cs);
	return;
}
void render_lt(conf co) {
	f82 rsif = f82_set(co.rsi.x, co.rsi.y);
	f8 ln = rsif.x * rsif.y * co.sn;
	scene sc = scene_read(co.ps);
	bvh_node *bns = malloc(sizeof(bvh_node) * co.bnsna);
	i4s bnsn = bvh_makert(&sc, co.bnsna, bns);
	printf("bvh usage:%d/%d\n", bnsn, co.bnsna);
	i4s cdfn;
	f8 *cdf;
	i4s *fis;
	make_cdf_scene(sc, &cdfn, &cdf, &fis);
	omp_set_num_threads(co.tn);
	f83 *css[omp_get_max_threads()];
	f8 t = timer();
	#pragma omp parallel
	{
		css[omp_get_thread_num()] = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
		f83 *cs = css[omp_get_thread_num()];
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			cs[i] = f83_set1(0.0);
		}
		i4s *st = malloc(sizeof(i4s) * bnsn);
		i8u2 s = x128p_seed(co.s);
		s = x128p_jumpn(s, omp_get_thread_num());
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			for (i4s j = 0; j < co.sn; ++j) {
				f83 tp, o, d;
				med m;
				light_sample_1(sc, cdfn, cdf, fis, &s, &tp, &o, &d, &m);
				trav_lt(sc, co.en, co.rsi, cs, rsif, co.rsfh, co.lr, co.lf, ln, bns, st, &s, o, d, tp, m);
			}
		}
		free(st);
	}
	t = timer() - t;
	printf("time:%f\n", t);
	free(cdf);
	free(fis);
	free(bns);
	scene_free(sc);
	f83 *cs = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
	for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
		cs[i] = f83_set1(0.0);
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		for (i4s j = 0; j < co.rsi.x * co.rsi.y; ++j) {
			cs[j] = f83_add(cs[j], f83_scl(1.0 / omp_get_max_threads(), css[i][j]));
		}
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		free(css[i]);
	}
	image_write(co.pr, co.rsi, cs);
	free(cs);
	return;
}
void render_ltd(conf co) {
	f82 rsif = f82_set(co.rsi.x, co.rsi.y);
	f8 ln = rsif.x * rsif.y * co.sn;
	scene sc = scene_read(co.ps);
	bvh_node *bns = malloc(sizeof(bvh_node) * co.bnsna);
	i4s bnsn = bvh_makert(&sc, co.bnsna, bns);
	printf("bvh usage:%d/%d\n", bnsn, co.bnsna);
	i4s cdfn;
	f8 *cdf;
	i4s *fis;
	make_cdf_scene(sc, &cdfn, &cdf, &fis);
	omp_set_num_threads(co.tn);
	f83 *css[omp_get_max_threads()];
	f8 t = timer();
	#pragma omp parallel
	{
		css[omp_get_thread_num()] = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
		f83 *cs = css[omp_get_thread_num()];
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			cs[i] = f83_set1(0.0);
		}
		i4s *st = malloc(sizeof(i4s) * bnsn);
		i8u2 s = x128p_seed(co.s);
		s = x128p_jumpn(s, omp_get_thread_num());
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			for (i4s j = 0; j < co.sn; ++j) {
				i4s idl;
				f83 tpl, tp, d;
				f82 uvl;
				med m;
				f8 pt;
				light_sample_1_f(sc, cdfn, cdf, fis, &s, &idl, &tpl, &tp, &uvl, &d, &m, &pt);
				face fl = sc.fs[idl];
				f83 o = tri_interp(get_vs(sc.vs, fl.vis), uvl);
				while (0 < co.en) {
					f83 tpe;
					f83 od;
					f8 pf;
					eye_sample_0_f(co.lr, &s, &tpe, &od, &pf);
					f83 dd = f83_sub(o, od);
					f8 td = f83_len(dd);
					dd = f83_scl(1.0 / td, dd);
					f83 nge = disc_ng(), nse = disc_ng();
					f8 cngoe = f83_dot(nge, dd), cnsoe = f83_dot(nse, dd);
					i4s2 pii;
					f83 fse;
					if (!eye_eval_00(co.rsi, rsif, co.rsfh, co.lr, co.lf, od, o, &pii, &fse)) {
						break;
					}
					f83 ngl = tri_ng(get_vs(sc.vs, fl.vis)), nsl = tri_ns(get_vs(sc.ns, fl.nis), uvl);
					f8 cngol = -f83_dot(ngl, dd), cnsol = -f83_dot(nsl, dd);
					if (!surf_kern(cngol, cnsol)) {
						break;
					}
					f83 fsl;
					f8 psbe;
					if (!edf_eval_b(sc.es[fl.ei], cnsol, &fsl, &psbe)) {
						break;
					}
					if (!bvh_visrt(sc.vs, sc.fs, bns, st, od, dd, td)) {
						break;
					}
					med md;
					medt(cngol, fl.mis, sc.ms, &md);
					f83 trd = trmit(m.a, td);
					f8 pb = psbe * f8_abs(cnsol * cngoe) / f8_p2(td);
					f8 w = 1.0 / (1.0 + f8_p2(pb / pf));
					f8 g = f8_abs(cnsoe * cnsol) / f8_p2(td);
					cs[pii.y * co.rsi.x + pii.x] = f83_add(cs[pii.y * co.rsi.x + pii.x], f83_scl(w * g / ln, f83_mul(trd, f83_mul(f83_mul(tpe, fse), f83_mul(tpl, fsl)))));
					break;
				}
				trav_ltd(sc, co.en, co.rsi, cs, rsif, co.rsfh, co.lr, co.lf, ln, bns, st, &s, o, d, tp, m, pt);
			}
		}
		free(st);
	}
	t = timer() - t;
	printf("time:%f\n", t);
	free(cdf);
	free(fis);
	free(bns);
	scene_free(sc);
	f83 *cs = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
	for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
		cs[i] = f83_set1(0.0);
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		for (i4s j = 0; j < co.rsi.x * co.rsi.y; ++j) {
			cs[j] = f83_add(cs[j], f83_scl(1.0 / omp_get_max_threads(), css[i][j]));
		}
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		free(css[i]);
	}
	image_write(co.pr, co.rsi, cs);
	free(cs);
	return;
}
void render_pte(conf co) {
	f82 rsif = f82_set(co.rsi.x, co.rsi.y);
	f8 ln = rsif.x * rsif.y * co.sn;
	scene sc = scene_read(co.ps);
	bvh_node *bns = malloc(sizeof(bvh_node) * co.bnsna);
	i4s bnsn = bvh_makert(&sc, co.bnsna, bns);
	printf("bvh usage:%d/%d\n", bnsn, co.bnsna);
	i4s2 rsie;
	f83 *cse;
	image_read(co.pe, &rsie, &cse);
	i4s2 rsied = i4s2_set(rsie.x - 1, rsie.y - 1);
	f82 rsiedf = f82_set(rsied.x, rsied.y);
	omp_set_num_threads(co.tn);
	f83 *css[omp_get_max_threads()];
	f8 t = timer();
	#pragma omp parallel
	{
		css[omp_get_thread_num()] = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
		f83 *cs = css[omp_get_thread_num()];
		i4s *st = malloc(sizeof(i4s) * bnsn);
		i8u2 s = x128p_seed(co.s);
		s = x128p_jumpn(s, omp_get_thread_num());
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			f83 c = f83_set1(0.0);
			for (i4s j = 0; j < co.sn; ++j) {
				f82 pif = f82_set(i % co.rsi.x + x128p_f8(&s), i / co.rsi.x + x128p_f8(&s));
				f83 tp, o, d;
				eye_sample_1(rsif, co.rsfh, co.lr, co.lf, &s, pif, &tp, &o, &d);
				c = f83_add(c, trav_pte(sc, co.en, rsie, cse, co.m, rsiedf, bns, st, &s, 1, o, d, tp));
			}
			cs[i] = f83_scl(1.0 / ln, c);
		}
		free(st);
	}
	t = timer() - t;
	printf("time:%f\n", t);
	free(cse);
	free(bns);
	scene_free(sc);
	f83 *cs = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
	for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
		cs[i] = f83_set1(0.0);
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		for (i4s j = 0; j < co.rsi.x * co.rsi.y; ++j) {
			cs[j] = f83_add(cs[j], f83_scl(1.0 / omp_get_max_threads(), css[i][j]));
		}
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		free(css[i]);
	}
	image_write(co.pr, co.rsi, cs);
	free(cs);
	return;
}
void render_pted(conf co) {
	f82 rsif = f82_set(co.rsi.x, co.rsi.y);
	f8 ln = rsif.x * rsif.y * co.sn;
	scene sc = scene_read(co.ps);
	bvh_node *bns = malloc(sizeof(bvh_node) * co.bnsna);
	i4s bnsn = bvh_makert(&sc, co.bnsna, bns);
	printf("bvh usage:%d/%d\n", bnsn, co.bnsna);
	i4s2 rsie;
	f83 *cse;
	image_read(co.pe, &rsie, &cse);
	i4s2 rsied = i4s2_set(rsie.x - 1, rsie.y - 1);
	f82 rsiedf = f82_set(rsied.x, rsied.y);
	i4s cdfne;
	f8 *cdfe;
	make_cdf_env(rsie, cse, rsied, rsiedf, &cdfne, &cdfe);
	f833 mi = f833_trapo(co.m);
	omp_set_num_threads(co.tn);
	f83 *css[omp_get_max_threads()];
	f8 t = timer();
	#pragma omp parallel
	{
		css[omp_get_thread_num()] = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
		f83 *cs = css[omp_get_thread_num()];
		i4s *st = malloc(sizeof(i4s) * bnsn);
		i8u2 s = x128p_seed(co.s);
		s = x128p_jumpn(s, omp_get_thread_num());
		for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
			f83 c = f83_set1(0.0);
			for (i4s j = 0; j < co.sn; ++j) {
				f82 pif = f82_set(i % co.rsi.x + x128p_f8(&s), i / co.rsi.x + x128p_f8(&s));
				f83 tpe, tp, o, d;
				f8 pt;
				eye_sample_1_f(rsif, co.rsfh, co.lr, co.lf, &s, pif, &tpe, &tp, &o, &d, &pt);
				c = f83_add(c, trav_pted(sc, co.en, rsie, cse, co.m, mi, rsied, rsiedf, cdfne, cdfe, bns, st, &s, 1, o, d, tp, pt));
			}
			cs[i] = f83_scl(1.0 / ln, c);
		}
		free(st);
	}
	t = timer() - t;
	printf("time:%f\n", t);
	free(cdfe);
	free(cse);
	free(bns);
	scene_free(sc);
	f83 *cs = malloc(sizeof(f83) * co.rsi.x * co.rsi.y);
	for (i4s i = 0; i < co.rsi.x * co.rsi.y; ++i) {
		cs[i] = f83_set1(0.0);
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		for (i4s j = 0; j < co.rsi.x * co.rsi.y; ++j) {
			cs[j] = f83_add(cs[j], f83_scl(1.0 / omp_get_max_threads(), css[i][j]));
		}
	}
	for (i4s i = 0; i < omp_get_max_threads(); ++i) {
		free(css[i]);
	}
	image_write(co.pr, co.rsi, cs);
	free(cs);
	return;
}
i4s main(i4s argc, i1 **argv) {
	conf co = conf_read(argv[1]);
	if (strcmp(co.t, "pt") == 0) {
		render_pt(co);
		return 0;
	}
	if (strcmp(co.t, "ptd") == 0) {
		render_ptd(co);
		return 0;
	}
	if (strcmp(co.t, "lt") == 0) {
		render_lt(co);
		return 0;
	}
	if (strcmp(co.t, "ltd") == 0) {
		render_ltd(co);
		return 0;
	}
	if (strcmp(co.t, "pte") == 0) {
		render_pte(co);
		return 0;
	}
	if (strcmp(co.t, "pted") == 0) {
		render_pted(co);
		return 0;
	}
	return 0;
}
