#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define CAST(t, x) (t)(x)
#define ALLOC(t, n) malloc(sizeof(t) * (n))
#define PUSH(arr, ptr, val) arr[(ptr)++] = val
#define POP(arr, ptr) arr[--(ptr)]
#define DEQ(arr, ptr) arr[(ptr)++]

#define PI 3.141592653589793
#define EPS 1.0e-12
#define SEED 0

typedef struct vec VEC;
typedef struct edf EDF;
typedef struct bsdf BSDF;
typedef struct medium MEDIUM;
typedef struct mat MAT;
typedef struct face FACE;
typedef struct scene SCENE;
typedef struct prng_state PRNG_STATE;

struct vec {
	double x, y, z;
};
struct edf {
	int t;
	VEC a;
	double r;
};
struct bsdf {
	int t;
	VEC a;
	double r;
	VEC n;
	VEC k;
};
struct medium {
	double n;
	VEC a;
};
struct mat {
	int ie, ib, imp, imn;
};
struct face {
	int im, iv0, iv1, iv2, in0, in1, in2;
};
struct scene {
	int nedfs, nbsdfs, nmeds, nmats, nverts, nnorms, nfaces;
	EDF *edfs;
	BSDF *bsdfs;
	MEDIUM *meds;
	MAT *mats;
	VEC *verts;
	VEC *norms;
	FACE *faces;
};
struct prng_state {
	uint64_t s0, s1;
};

double time_get() {
	return clock() / CAST(double, CLOCKS_PER_SEC);
}
void scene_del(SCENE scene) {
	free(scene.edfs);
	free(scene.bsdfs);
	free(scene.meds);
	free(scene.mats);
	free(scene.verts);
	free(scene.norms);
	free(scene.faces);
	return;
}
SCENE read_scene(char *path) {
	printf("proc read-scene\n");
	double t = time_get();
	
	FILE *file = fopen(path, "r");
	if (file == NULL) {
		puts("file error");
		exit(0);
	}
	
	SCENE scene;
	
	int nbuf = 256, nid = 256;
	char buf[nbuf], id[nid];
	
	while (true) {
		if (fgets(buf, nbuf, file) == NULL) {
			break;
		}
		sscanf(buf, "%s", id);
		if (strcmp(id, "edfs") == 0) {
			sscanf(buf, "%*s%d", &scene.nedfs);
			scene.edfs = ALLOC(EDF, scene.nedfs);
			for (int i = 0; i < scene.nedfs; ++i) {
				fgets(buf, nbuf, file);
				EDF edf;
				sscanf(buf, "%d%lf%lf%lf%lf", &edf.t, &edf.a.x, &edf.a.y, &edf.a.z, &edf.r);
				scene.edfs[i] = edf;
			}
			continue;
		}
		if (strcmp(id, "bsdfs") == 0) {
			sscanf(buf, "%*s%d", &scene.nbsdfs);
			scene.bsdfs = ALLOC(BSDF, scene.nbsdfs);
			for (int i = 0; i < scene.nbsdfs; ++i) {
				fgets(buf, nbuf, file);
				BSDF bsdf;
				sscanf(buf, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &bsdf.t, &bsdf.a.x, &bsdf.a.y, &bsdf.a.z,
				&bsdf.r, &bsdf.n.x, &bsdf.n.y, &bsdf.n.z, &bsdf.k.x, &bsdf.k.y, &bsdf.k.z);
				scene.bsdfs[i] = bsdf;
			}
			continue;
		}
		if (strcmp(id, "media") == 0) {
			sscanf(buf, "%*s%d", &scene.nmeds);
			scene.meds = ALLOC(MEDIUM, scene.nmeds);
			for (int i = 0; i < scene.nmeds; ++i) {
				fgets(buf, nbuf, file);
				MEDIUM med;
				sscanf(buf, "%lf%lf%lf%lf", &med.n, &med.a.x, &med.a.y, &med.a.z);
				scene.meds[i] = med;
			}
			continue;
		}
		if (strcmp(id, "materials") == 0) {
			sscanf(buf, "%*s%d", &scene.nmats);
			scene.mats = ALLOC(MAT, scene.nmats);
			for (int i = 0; i < scene.nmats; ++i) {
				fgets(buf, nbuf, file);
				MAT mat;
				sscanf(buf, "%d%d%d%d", &mat.ie, &mat.ib, &mat.imp, &mat.imn);
				scene.mats[i] = mat;
			}
			continue;
		}
		if (strcmp(id, "verticies") == 0) {
			sscanf(buf, "%*s%d", &scene.nverts);
			scene.verts = ALLOC(VEC, scene.nverts);
			for (int i = 0; i < scene.nverts; ++i) {
				fgets(buf, nbuf, file);
				VEC vert;
				sscanf(buf, "%lf%lf%lf", &vert.x, &vert.y, &vert.z);
				scene.verts[i] = vert;
			}
			continue;
		}
		if (strcmp(id, "normals") == 0) {
			sscanf(buf, "%*s%d", &scene.nnorms);
			scene.norms = ALLOC(VEC, scene.nnorms);
			for (int i = 0; i < scene.nnorms; ++i) {
				fgets(buf, nbuf, file);
				VEC norm;
				sscanf(buf, "%lf%lf%lf", &norm.x, &norm.y, &norm.z);
				scene.norms[i] = norm;
			}
			continue;
		}
		if (strcmp(id, "faces") == 0) {
			sscanf(buf, "%*s%d", &scene.nfaces);
			scene.faces = ALLOC(FACE, scene.nfaces);
			for (int i = 0; i < scene.nfaces; ++i) {
				fgets(buf, nbuf, file);
				FACE face;
				sscanf(buf, "%d%d%d%d%d%d%d", &face.im, &face.iv0, &face.iv1, &face.iv2, &face.in0, &face.in1, &face.in2);
				scene.faces[i] = face;
			}
			continue;
		}
	}
	
	fclose(file);
	
	t = time_get() - t;
	
	printf("time %f\n", t);
	
	return scene;
}

int imin(int a, int b) {
	return a < b ? a : b;
}
int imax(int a, int b) {
	return a > b ? a : b;
}
int iclamp(int a, int b, int x) {
	return imin(imax(x, a), b);
}
double clamp(double a, double b, double x) {
	return fmin(fmax(x, a), b);
}
double sign(double x) {
	return x < 0.0 ? -1.0 : x > 0.0 ? 1.0 : 0.0;
}
double sqr(double x) {
	return x * x;
}
double cube(double x) {
	return x * x * x;
}
double quad(double x) {
	return x * x * x * x;
}
double opst(double x) {
	return sqrt(fmax(0.0, 1.0 - sqr(x)));
}

//splitmix64
uint64_t sm64_next(uint64_t *s) {
	uint64_t x = *s += 0x9E3779B97F4A7C15;
	x = 0xBF58476D1CE4E5B9 * (x ^ (x >> 30));
	x = 0x94D049BB133111EB * (x ^ (x >> 27));
	x = x ^ (x >> 31);
	return x;
}
//xoroshiro128+
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
	PRNG_STATE jump = {0xBEAC0467EBA5FACB, 0xD86B048B86AA9922};
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

VEC vec_init(double x, double y, double z) {
	VEC v = {x, y, z};
	return v;
}
VEC vec_init1(double a) {
	return vec_init(a, a, a);
}
VEC vec_add(VEC a, VEC b) {
	return vec_init(a.x + b.x, a.y + b.y, a.z + b.z);
}
VEC vec_sub(VEC a, VEC b) {
	return vec_init(a.x - b.x, a.y - b.y, a.z - b.z);
}
VEC vec_mul(VEC a, VEC b) {
	return vec_init(a.x * b.x, a.y * b.y, a.z * b.z);
}
VEC vec_scale(double a, VEC b) {
	return vec_init(a * b.x, a * b.y, a * b.z);
}
VEC vec_cross(VEC a, VEC b) {
	return vec_init(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
double vec_dot(VEC a, VEC b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
double vec_len2(VEC a) {
	return vec_dot(a, a);
}
double vec_len(VEC a) {
	return sqrt(vec_len2(a));
}
VEC vec_norm(VEC a) {
	return vec_scale(1.0 / vec_len(a), a);
}
VEC vec_exp(VEC a) {
	return vec_init(exp(a.x), exp(a.y), exp(a.z));
}
VEC vec_min(VEC a, VEC b) {
	return vec_init(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z));
}
VEC vec_max(VEC a, VEC b) {
	return vec_init(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z));
}
VEC vec_abs(VEC a) {
	return vec_init(fabs(a.x), fabs(a.y), fabs(a.z));
}
VEC vec_trafo(VEC x, VEC y, VEC z, VEC a) {
	return vec_add(vec_add(vec_scale(a.x, x), vec_scale(a.y, y)), vec_scale(a.z, z));
}
bool vec_eq(VEC a, VEC b) {
	VEC d = vec_sub(a, b);
	return fabs(d.x) < EPS && fabs(d.y) < EPS && fabs(d.z) < EPS;
}
void vec_pr(VEC a) {
	printf("%f %f %f\n", a.x, a.y, a.z);
	return;
}
typedef struct matrix MATRIX;
struct matrix {
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
};
MATRIX matrix_init(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) {
	MATRIX m = {xx, xy, xz, yx, yy, yz, zx, zy, zz};
	return m;
}
MATRIX matrix_initb(VEC x, VEC y, VEC z) {
	return matrix_init(x.x, y.x, z.x, x.y, y.y, z.y, x.z, y.z, z.z);
}
MATRIX matrix_rot_y(double t, double p) {
	double cost = cos(t), sint = sin(t), cosp = cos(p), sinp = sin(p);
	return matrix_init(cosp, sint * sinp, cost * sinp, 0.0, cost, -sint, -sinp, sint * cosp, cost * cosp);
}
MATRIX matrix_rot_y_inv(double t, double p) {
	double cost = cos(t), sint = sin(t), cosp = cos(p), sinp = sin(p);
	return matrix_init(cosp, 0.0, -sinp, sint * sinp, cost, sint * cosp, cost * sinp, -sint, cost * cosp);
}
VEC matrix_trafo(MATRIX m, VEC v) {
	return vec_init(m.xx * v.x + m.xy * v.y + m.xz * v.z, m.yx * v.x + m.yy * v.y + m.yz * v.z, m.zx * v.x + m.zy * v.y + m.zz * v.z);
}
MATRIX matrix_trapo(MATRIX a) {
	return matrix_init(a.xx, a.yx, a.zx, a.xy, a.yy, a.zy, a.xz, a.yz, a.zz);
}

typedef struct aabb AABB;
typedef struct bvh_node BVH_NODE;
typedef struct bvh_sort_entry BVH_SORT_ENTRY;
typedef struct bvh_build_query BVH_BUILD_QUERY;

struct aabb {
	VEC min, max;
};
struct bvh_sort_entry {
	int i;
	VEC c;
};
struct bvh_node {
	AABB b;
	int l, r;
};
struct bvh_build_query {
	int s, e;
};

AABB aabb_default(void) {
	AABB box = {vec_init1(INFINITY), vec_init1(-INFINITY)};
	return box;
}
double aabb_area(AABB box) {
	VEC d = vec_sub(box.max, box.min);
	return 2.0 * (d.x * d.y + d.y * d.z + d.z * d.x);
}
AABB aabb_merge(AABB a, AABB b) {
	AABB box = {vec_min(a.min, b.min), vec_max(a.max, b.max)};
	return box;
}
bool aabb_isect(AABB box, VEC o, VEC d, double *t, double k) {
	double tmin = -INFINITY, tmax = INFINITY;
	if (fabs(d.x) > EPS) {
		double tx0 = (box.min.x - o.x) / d.x;
		double tx1 = (box.max.x - o.x) / d.x;
		tmin = fmax(tmin, fmin(tx0, tx1));
		tmax = fmin(tmax, fmax(tx0, tx1));
		if (tmax < tmin) {
			return false;
		}
	}
	if (fabs(d.y) > EPS) {
		double ty0 = (box.min.y - o.y) / d.y;
		double ty1 = (box.max.y - o.y) / d.y;
		tmin = fmax(tmin, fmin(ty0, ty1));
		tmax = fmin(tmax, fmax(ty0, ty1));
		if (tmax < tmin) {
			return false;
		}
	}
	if (fabs(d.z) > EPS) {
		double tz0 = (box.min.z - o.z) / d.z;
		double tz1 = (box.max.z - o.z) / d.z;
		tmin = fmax(tmin, fmin(tz0, tz1));
		tmax = fmin(tmax, fmax(tz0, tz1));
		if (tmax < tmin) {
			return false;
		}
	}
	*t = fmax(0.0, tmin);
	return tmin < k - EPS && tmax > EPS;
}
double tri_area(VEC v0, VEC v1, VEC v2) {
	return 0.5 * vec_len(vec_cross(vec_sub(v1, v0), vec_sub(v2, v0)));
}
double tri_pdf(VEC v0, VEC v1, VEC v2) {
	return 1.0 / tri_area(v0, v1, v2);
}
void tri_sample(PRNG_STATE *s, double *u, double *v) {
	double a = sqrt(prng_db(s));
	*u = 1.0 - a;
	*v = prng_db(s) * a;
	return;
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
VEC tri_cent(VEC v0, VEC v1, VEC v2) {
	return vec_scale(1.0 / 3.0, vec_add(vec_add(v0, v1), v2));
}
VEC tri_normg(VEC v0, VEC v1, VEC v2) {
	return vec_norm(vec_cross(vec_sub(v1, v0), vec_sub(v2, v0)));
}
VEC tri_norms(VEC n0, VEC n1, VEC n2, double u, double v) {
	return vec_norm(vec_add(vec_add(vec_scale(1.0 - u - v, n0), vec_scale(u, n1)), vec_scale(v, n2)));
}
VEC tri_point(VEC v0, VEC v1, VEC v2, double u, double v) {
	return vec_add(vec_add(vec_scale(1.0 - u - v, v0), vec_scale(u, v1)), vec_scale(v, v2));
}
bool tri_isect(VEC v0, VEC v1, VEC v2, VEC o, VEC d, VEC *tuv) {
	VEC e0 = vec_sub(v1, v0);
	VEC e1 = vec_sub(v2, v0);
	VEC pv = vec_cross(d, e1);
	double deti = vec_dot(e0, pv);
	if (fabs(deti) < EPS) {
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
	double t = deti * vec_dot(e1, qv);
	if (t < EPS || t > tuv->x - EPS) {
		return false;
	}
	*tuv = vec_init(t, u, v);
	return true;
}

void face_get_verts(VEC *verts, FACE face, VEC *v0, VEC *v1, VEC *v2) {
	*v0 = verts[face.iv0];
	*v1 = verts[face.iv1];
	*v2 = verts[face.iv2];
	return;
}
void face_get_norms(VEC *norms, FACE face, VEC *n0, VEC *n1, VEC *n2) {
	*n0 = norms[face.in0];
	*n1 = norms[face.in1];
	*n2 = norms[face.in2];
	return;
}

int bvh_cmp_x(void const *_a, void const *_b) {
	BVH_SORT_ENTRY const *a = _a;
	BVH_SORT_ENTRY const *b = _b;
	return a->c.x < b->c.x ? -1 : a->c.x > b->c.x ? 1 : 0;
}
int bvh_cmp_y(void const *_a, void const *_b) {
	BVH_SORT_ENTRY const *a = _a;
	BVH_SORT_ENTRY const *b = _b;
	return a->c.y < b->c.y ? -1 : a->c.y > b->c.y ? 1 : 0;
}
int bvh_cmp_z(void const *_a, void const *_b) {
	BVH_SORT_ENTRY const *a = _a;
	BVH_SORT_ENTRY const *b = _b;
	return a->c.z < b->c.z ? -1 : a->c.z > b->c.z ? 1 : 0;
}
void bvh_build(VEC *verts, FACE *faces, BVH_SORT_ENTRY *ifs, double *al, double *ar, BVH_NODE *bvh_nodes, BVH_BUILD_QUERY *stack, BVH_BUILD_QUERY q, int ptr_deq, int *ptr_enq) {
	double tt = 1.0, tb = 2.0;
	AABB box = aabb_default();
	for (int i = q.s; i < q.e; ++i) {
		VEC v0, v1, v2;
		face_get_verts(verts, faces[ifs[i].i], &v0, &v1, &v2);
		box = aabb_merge(box, tri_aabb(v0, v1, v2));
	}
	double areai = 1.0 / aabb_area(box);
	double cost = tt * (q.e - q.s);
	int axis = -1, piv = -1;
	for (int i = 0; i < 3; ++i) {
		switch (i) {
			case 0 : {
				qsort(ifs + q.s, q.e - q.s, sizeof(BVH_SORT_ENTRY), bvh_cmp_x);
			} break;
			case 1 : {
				qsort(ifs + q.s, q.e - q.s, sizeof(BVH_SORT_ENTRY), bvh_cmp_y);
			} break;
			case 2 : {
				qsort(ifs + q.s, q.e - q.s, sizeof(BVH_SORT_ENTRY), bvh_cmp_z);
			} break;
		}
		AABB boxc;
		boxc = aabb_default();
		al[0] = 0.0;
		for (int j = q.s; j < q.e; ++j) {
			VEC v0, v1, v2;
			face_get_verts(verts, faces[ifs[j].i], &v0, &v1, &v2);
			boxc = aabb_merge(boxc, tri_aabb(v0, v1, v2));
			al[j - q.s + 1] = aabb_area(boxc);
		}
		boxc = aabb_default();
		ar[q.e - q.s] = 0.0;
		for (int j = q.e - 1; j >= q.s; --j) {
			VEC v0, v1, v2;
			face_get_verts(verts, faces[ifs[j].i], &v0, &v1, &v2);
			boxc = aabb_merge(boxc, tri_aabb(v0, v1, v2));
			ar[j - q.s] = aabb_area(boxc);
		}
		for (int j = q.s; j <= q.e; ++j) {
			double cost1 = tb + tt * areai * (al[j - q.s] * (j - q.s) + ar[j - q.s] * (q.e - j));
			if (cost1 < cost) {
				cost = cost1;
				axis = i;
				piv = j;
			}
		}
	}
	if (axis == -1) {
		BVH_NODE node = {box, q.s - q.e, q.s};
		bvh_nodes[ptr_deq] = node;
		return;
	}
	switch (axis) {
		case 0 : {
			qsort(ifs + q.s, q.e - q.s, sizeof(BVH_SORT_ENTRY), bvh_cmp_x);
		} break;
		case 1 : {
			qsort(ifs + q.s, q.e - q.s, sizeof(BVH_SORT_ENTRY), bvh_cmp_y);
		} break;
		case 2 : {
			qsort(ifs + q.s, q.e - q.s, sizeof(BVH_SORT_ENTRY), bvh_cmp_z);
		} break;
	}
	BVH_NODE node = {box, *ptr_enq + 1, *ptr_enq + 2};
	bvh_nodes[ptr_deq] = node;
	BVH_BUILD_QUERY ql = {q.s, piv};
	PUSH(stack, *ptr_enq, ql);
	BVH_BUILD_QUERY qr = {piv, q.e};
	PUSH(stack, *ptr_enq, qr);
	return;
}
int bvh_build_wrap(SCENE *scene, int nalbvh_nodes, BVH_NODE *bvh_nodes) {
	printf("proc bvh-build\n");
	double t = time_get();
	BVH_SORT_ENTRY *ifs = ALLOC(BVH_SORT_ENTRY, scene->nfaces);
	for (int i = 0; i < scene->nfaces; ++i) {
		VEC v0, v1, v2;
		face_get_verts(scene->verts, scene->faces[i], &v0, &v1, &v2);
		BVH_SORT_ENTRY e = {i, tri_cent(v0, v1, v2)};
		ifs[i] = e;
	}
	double *al = ALLOC(double, scene->nfaces + 1);
	double *ar = ALLOC(double, scene->nfaces + 1);
	BVH_BUILD_QUERY *stack = ALLOC(BVH_BUILD_QUERY, nalbvh_nodes);
	int ptr_deq = 0, ptr_enq = 0;
	BVH_BUILD_QUERY q = {0, scene->nfaces};
	while (true) {
		bvh_build(scene->verts, scene->faces, ifs, al, ar, bvh_nodes, stack, q, ptr_deq, &ptr_enq);
		if (ptr_deq == ptr_enq) {
			break;
		}
		q = DEQ(stack, ptr_deq);
	}
	free(stack);
	free(ar);
	free(al);
	FACE *faces = ALLOC(FACE, scene->nfaces);
	for (int i = 0; i < scene->nfaces; ++i) {
		faces[i] = scene->faces[ifs[i].i];
	}
	free(ifs);
	free(scene->faces);
	scene->faces = faces;
	t = time_get() - t;
	printf("time %f\n", t);
	printf("usage %d/%d\n", ptr_deq + 1, nalbvh_nodes);
	return ptr_deq + 1;
}
double disc_area(double r) {
	return PI * sqr(r);
}
double disc_pdf(double r) {
	return 1.0 / disc_area(r);
}
VEC disc_sample(double r, PRNG_STATE *s) {
	double u = r * sqrt(prng_db(s)), v = 2.0 * PI * prng_db(s);
	return vec_scale(u, vec_init(cos(v), sin(v), 0.0));
}
VEC disc_normal(void) {
	return vec_init(0.0, 0.0, 1.0);
}
bool disc_isect(double r, VEC o, VEC d, VEC *tuv) {
	if (fabs(d.z) < EPS) {
		return false;
	}
	double t = -o.z / d.z;
	if (t < EPS || t > tuv->x - EPS) {
		return false;
	}
	VEC x = vec_add(o, vec_scale(t, d));
	double r2 = vec_len2(x);
	if (r2 > sqr(r)) {
		return false;
	}
	*tuv = vec_init(t, sqrt(r2), atan2(x.y, x.x) + (x.y < 0.0 ? 2.0 * PI : 0.0));
	return true;
}
void bvh_traverse(VEC *verts, FACE *faces, BVH_NODE *bvh_nodes, BVH_NODE bvh_node, int *stack, int *ptr, VEC o, VEC d, VEC *tuv, int *id) {
	if (bvh_node.l < 0) {
		for (int i = 0; i < -bvh_node.l; ++i) {
			VEC v0, v1, v2;
			face_get_verts(verts, faces[bvh_node.r + i], &v0, &v1, &v2);
			if (tri_isect(v0, v1, v2, o, d, tuv)) {
				*id = bvh_node.r + i;
			}
		}
	} else {
		double tl;
		bool hl = aabb_isect(bvh_nodes[bvh_node.l].b, o, d, &tl, tuv->x);
		double tr;
		bool hr = aabb_isect(bvh_nodes[bvh_node.r].b, o, d, &tr, tuv->x);
		if (hl && hr) {
			if (tl < tr) {
				PUSH(stack, *ptr, bvh_node.r);
				PUSH(stack, *ptr, bvh_node.l);
			} else {
				PUSH(stack, *ptr, bvh_node.l);
				PUSH(stack, *ptr, bvh_node.r);
			}
		} else if (hl) {
			PUSH(stack, *ptr, bvh_node.l);
		} else if (hr) {
			PUSH(stack, *ptr, bvh_node.r);
		}
	}
	return;
}
void bvh_traverse_wrap(VEC *verts, FACE *faces, BVH_NODE *bvh_nodes, int *stack, VEC o, VEC d, VEC *tuv, int *id) {
	BVH_NODE bvh_node = bvh_nodes[0];
	double t;
	if (aabb_isect(bvh_node.b, o, d, &t, tuv->x)) {
		int ptr = 0;
		while (true) {
			bvh_traverse(verts, faces, bvh_nodes, bvh_node, stack, &ptr, o, d, tuv, id);
			if (ptr == 0) {
				break;
			}
			bvh_node = bvh_nodes[POP(stack, ptr)];
		}
	}
}
bool traverse(double r, VEC *verts, FACE *faces, BVH_NODE *bvh_nodes, int *stack, VEC o, VEC d, VEC *tuv, int *id) {
	*id = -2;
	if (disc_isect(r, o, d, tuv)) {
		*id = -1;
	}
	bvh_traverse_wrap(verts, faces, bvh_nodes, stack, o, d, tuv, id);
	return *id > -2;
}

bool bvh_visible(VEC *verts, FACE *faces, BVH_NODE *bvh_nodes, BVH_NODE bvh_node, int *stack, int *ptr, VEC o, VEC d, VEC *tuv) {
	if (bvh_node.l < 0) {
		for (int i = 0; i < -bvh_node.l; ++i) {
			VEC v0, v1, v2;
			face_get_verts(verts, faces[bvh_node.r + i], &v0, &v1, &v2);
			if (tri_isect(v0, v1, v2, o, d, tuv)) {
				return false;
			}
		}
	} else {
		double tl;
		bool hl = aabb_isect(bvh_nodes[bvh_node.l].b, o, d, &tl, tuv->x);
		double tr;
		bool hr = aabb_isect(bvh_nodes[bvh_node.r].b, o, d, &tr, tuv->x);
		if (hl && hr) {
			if (tl < tr) {
				PUSH(stack, *ptr, bvh_node.r);
				PUSH(stack, *ptr, bvh_node.l);
			} else {
				PUSH(stack, *ptr, bvh_node.l);
				PUSH(stack, *ptr, bvh_node.r);
			}
		} else if (hl) {
			PUSH(stack, *ptr, bvh_node.l);
		} else if (hr) {
			PUSH(stack, *ptr, bvh_node.r);
		}
	}
	return true;
}
bool bvh_visible_wrap(VEC *verts, FACE *faces, BVH_NODE *bvh_nodes, int *stack, VEC o, VEC d, VEC *tuv) {
	BVH_NODE bvh_node = bvh_nodes[0];
	double t;
	if (aabb_isect(bvh_node.b, o, d, &t, tuv->x)) {
		int ptr = 0;
		while (true) {
			if (!bvh_visible(verts, faces, bvh_nodes, bvh_node, stack, &ptr, o, d, tuv)) {
				return false;
			}
			if (ptr == 0) {
				break;
			}
			bvh_node = bvh_nodes[POP(stack, ptr)];
		}
	}
	return true;
}
bool visible(double r, VEC *verts, FACE *faces, BVH_NODE *bvh_nodes, int *stack, VEC o, VEC d, VEC *tuv) {
	if (disc_isect(r, o, d, tuv)) {
		return false;
	}
	if (!bvh_visible_wrap(verts, faces, bvh_nodes, stack, o, d, tuv)) {
		return false;
	}
	return true;
}

typedef struct image IMAGE;

struct image {
	int w, h;
	VEC *data;
};

IMAGE image_init(int w, int h) {
	IMAGE img = {w, h, ALLOC(VEC, w * h)};
	return img;
}
void image_del(IMAGE img) {
	free(img.data);
	return;
}
void image_fill(IMAGE img, VEC v) {
	for (int i = 0; i < img.w * img.h; ++i) {
		img.data[i] = v;
	}
	return;
}
IMAGE image_read(char *path) {
	FILE *file = fopen(path, "rb");
	if (file == NULL) {
		puts("file error");
		exit(0);
	}
	int w, h;
	int nbuf = 256;
	char buf[nbuf];
	fgets(buf, nbuf, file);
	fgets(buf, nbuf, file);
	sscanf(buf, "%d%d", &w, &h);
	fgets(buf, nbuf, file);
	float *data = ALLOC(float, 3 * w * h);
	fread(data, sizeof(float), 3 * w * h, file);
	IMAGE img = image_init(w, h);
	for (int i = 0; i < w * h; ++i) {
		img.data[i].x = data[3 * i + 0];
		img.data[i].y = data[3 * i + 1];
		img.data[i].z = data[3 * i + 2];
	}
	free(data);
	fclose(file);
	return img;
}
void image_write(char *path, IMAGE img) {
	FILE *file = fopen(path, "wb");
	if (file == NULL) {
		puts("file error");
		exit(0);
	}
	fprintf(file, "PF\n%d %d\n-1.0\n", img.w, img.h);
	float *data = ALLOC(float, 3 * img.w * img.h);
	for (int i = 0; i < img.w * img.h; ++i) {
		data[3 * i + 0] = img.data[i].x;
		data[3 * i + 1] = img.data[i].y;
		data[3 * i + 2] = img.data[i].z;
	}
	fwrite(data, sizeof(float), 3 * img.w * img.h, file);
	free(data);
	fclose(file);
	return;
}

typedef struct conf_debug CONF_DEBUG;

struct conf_debug {
	int num_alloc_bvh_nodes;
	int sample_rate;
	int image_width;
	int image_height;
	double fx, fy;
};

void film_dim(int w, int h, double fovh, double *fx, double *fy) {
	*fx = tan(fovh * PI / 360.0);
	*fy = *fx * h / w;
	return;
}
CONF_DEBUG read_conf_debug(char *path) {
	FILE *file = fopen(path, "r");
	if (file == NULL) {
		puts("file error");
		exit(0);
	}
	CONF_DEBUG conf;
	double fovh;
	int nbuf = 256, nid = 256;
	char buf[nbuf], id[nid];
	while (true) {
		if (fgets(buf, nbuf, file) == NULL) {
			break;
		}
		sscanf(buf, "%s", id);
		if (strcmp(id, "num-alloc-bvh") == 0) {
			sscanf(buf, "%*s%d", &conf.num_alloc_bvh_nodes);
			continue;
		}
		if (strcmp(id, "image-width") == 0) {
			sscanf(buf, "%*s%d", &conf.image_width);
			continue;
		}
		if (strcmp(id, "image-height") == 0) {
			sscanf(buf, "%*s%d", &conf.image_height);
			continue;
		}
		if (strcmp(id, "sample-rate") == 0) {
			sscanf(buf, "%*s%d", &conf.sample_rate);
			continue;
		}
		if (strcmp(id, "fovh") == 0) {
			sscanf(buf, "%*s%lf", &fovh);
			continue;
		}
	}
	fclose(file);
	film_dim(conf.image_width, conf.image_height, fovh, &conf.fx, &conf.fy);
	return conf;
}

VEC trafo_img2film(int w, int h, double fx, double fy, double x, double y) {
	return vec_init(fx * (1.0 - 2.0 * x / w), fy * (1.0 - 2.0 * y / h), 1.0);
}
void trafo_film2img(int w, int h, double fx, double fy, double x, double y, int *imgx, int *imgy) {
	*imgx = iclamp(0, w - 1, w * (1.0 - x / fx) / 2.0);
	*imgy = iclamp(0, h - 1, h * (1.0 - y / fy) / 2.0);
	return;
}

MATRIX matrix_tan(VEC z) {
	VEC x, y;
	if (fabs(z.x) > fabs(z.y)) {
		x = vec_scale(1.0 / sqrt(sqr(z.x) + sqr(z.z)), vec_init(-z.z, 0.0, z.x));
	} else {
		x = vec_scale(1.0 / sqrt(sqr(z.y) + sqr(z.z)), vec_init(0.0, z.z, -z.y));
	}
	y = vec_cross(z, x);
	return matrix_initb(x, y, z);
}
VEC pol2carty1(double cost, double sint, double cosp, double sinp) {
	return vec_init(sint * sinp, cost, sint * cosp);
}
VEC pol2carty0(double t, double p) {
	return pol2carty1(cos(t), sin(t), cos(p), sin(p));
}
VEC pol2cartz1(double cost, double sint, double cosp, double sinp) {
	return vec_init(sint * cosp, sint * sinp, cost);
}
VEC pol2cartz0(double t, double p) {
	return pol2cartz1(cos(t), sin(t), cos(p), sin(p));
}
VEC refl(VEC n, VEC d, double cos) {
	return vec_sub(d, vec_scale(2.0 * cos, n));
}
VEC refr(VEC n, VEC d, double eta, double cosi, double cost) {
  return vec_sub(vec_scale(eta, d), vec_scale(eta * cosi - cost, n));
}
double rgb2lum(VEC a) {
  return vec_dot(vec_init(0.2126, 0.7152, 0.0722), a);
}

double fres_diel(double eta, double cosi, double cost) {
	double perp = (eta * cosi - cost) / (eta * cosi + cost);
	double paral = (eta * cost - cosi) / (eta * cost + cosi);
	return (sqr(perp) + sqr(paral)) / 2.0;
}
double fres_cond(double etai, double k, double cosi) {
	cosi = fabs(cosi);
	double s = 1.0 - sqr(cosi);
	double i = sqr(etai) - sqr(k) - s;
	double w = sqrt(sqr(i) + sqr(2.0 * etai * k));
	double a = sqrt(fmax(0.0, (w + i) / 2.0));
	double perp = (w + sqr(cosi) - 2.0 * a * cosi) / (w + sqr(cosi) + 2.0 * a * cosi);
	double paral = (sqr(cosi) * w + sqr(s) - 2.0 * a * cosi * s) / (sqr(cosi) * w + sqr(s) + 2.0 * a * cosi * s);
	return (perp + perp * paral) / 2.0;
}
VEC fres_cond_rgb(VEC etai, VEC k, double cosi) {
	return vec_init(fres_cond(etai.x, k.x, cosi), fres_cond(etai.y, k.y, cosi), fres_cond(etai.z, k.z, cosi));
}

VEC ggx_sample_ndf(VEC n, double r, PRNG_STATE *s) {
	double u0 = prng_db(s), u1 = prng_db(s);
	double t = atan(sqrt(sqr(r) * u0 / (1.0 - u0))), p = 2.0 * PI * u1;
	return matrix_trafo(matrix_tan(n), pol2cartz0(t, p));
}
double ggx_d(double r, double cosnm) {
	return 1.0 / (PI * sqr(r * sqr(cosnm) + fabs(1.0 - sqr(cosnm)) / r));
}
double ggx_lam(double r, double cost) {
	if (cost > 1.0 - EPS) {
		return 0.0;
	}
	if (cost < EPS - 1.0) {
		return -1.0;
	}
	return (sign(cost) * sqrt(1.0 + (1.0 - sqr(cost)) * sqr(r / cost)) - 1.0) / 2.0;
}
double ggx_g1(double r, double cost) {
	if (fabs(cost) < EPS) {
		return 0.0;
	}
	return 1.0 / (1.0 + ggx_lam(r, cost));
}
double ggx_g2(double r, double cosni, double cosno) {
	return ggx_g1(r, cosni) * ggx_g1(r, cosno);
}
bool edf_insup(EDF edf, double cosno) {
	return cosno > 0.0;
}
bool edf_con(EDF edf) {
	switch (edf.t) {
		case 0 : {
			return false;
		}
		case 1 : {
			return true;
		}
		case 2 : {
			return true;
		}
		case 3 : {
			return false;
		}
		default : {
			return false;
		}
	}
}
bool edf_eval_f(EDF edf, double cosno, VEC *fs, double *psf) {
	switch (edf.t) {
		case 1 : {
			if (cosno < 0.0) {
				return false;
			}
			*fs = vec_mul(vec_init1(1.0 / PI), edf.a);
			*psf = 1.0 / PI;
			return true;
		}
		case 2 : {
			if (cosno < 0.0) {
				return false;
			}
			double d = ggx_d(edf.r, cosno);
			*fs = vec_scale(d, edf.a);
			*psf = d;
			return true;
		}
		default : {
			return false;
		}
	}
}
bool edf_eval(EDF edf, double cosno, VEC *fs) {
	switch (edf.t) {
		case 1 : {
			if (cosno < 0.0) {
				return false;
			}
			*fs = vec_mul(vec_init1(1.0 / PI), edf.a);
			return true;
		}
		case 2 : {
			if (cosno < 0.0) {
				return false;
			}
			double d = ggx_d(edf.r, cosno);
			*fs = vec_scale(d, edf.a);
			return true;
		}
		default : {
			return false;
		}
	}
}
VEC edf_sample_f(EDF edf, VEC n, PRNG_STATE *s, VEC *tp, double *psf) {
	switch (edf.t) {
		case 1 : {
			double p = 2.0 * PI * prng_db(s);
			double cost = sqrt(prng_db(s)), sint = opst(cost), cosp = cos(p), sinp = sin(p);
			VEC wo = matrix_trafo(matrix_tan(n), pol2cartz1(cost, sint, cosp, sinp));
			*tp = vec_mul(edf.a, *tp);
			*psf = 1.0 / PI;
			return wo;
		}
		case 2 : {
			VEC wo = ggx_sample_ndf(n, edf.r, s);
			*tp = vec_mul(edf.a, *tp);
			*psf = ggx_d(edf.r, vec_dot(n, wo));
			return wo;
		}
		case 3 : {
			VEC wo = n;
			*tp = vec_mul(edf.a, *tp);
			*psf = 1.0;
			return wo;
		}
		default : {
			return vec_init1(0.0);
		}
	}
}

VEC bsdf_eval(BSDF bsdf, double ni, double nt, VEC n, VEC wi, VEC wo) {
	switch (bsdf.t) {
		case 1 : {
			double cosi = vec_dot(n, wi), coso = vec_dot(n, wo);
			return vec_scale(cosi * coso < 0.0 ? 1.0 / PI : 0.0, bsdf.a);
		}
		default : {
			return vec_init1(0.0);
		}
	}
}
double bsdf_pdf_eye(BSDF bsdf, double ni, double nt, VEC n, VEC wi, VEC wo) {
	switch (bsdf.t) {
		case 1 : {
			double cosi = vec_dot(n, wi), coso = vec_dot(n, wo);
			return cosi * coso < 0.0 ? 1.0 / PI : 0.0;
		}
		default : {
			return 0.0;
		}
	}
}
VEC ggx_sample_vndf(VEC n, double r, VEC wi, PRNG_STATE *s) {
	MATRIX m = matrix_tan(n);
	MATRIX mi = matrix_trapo(m);
	wi = matrix_trafo(mi, wi);
	VEC Vh = vec_norm(vec_init(r * wi.x, r * wi.y, wi.z));
	double lensq = sqr(Vh.x) + sqr(Vh.y);
	VEC T1 = lensq > 0.0 ? vec_scale(1.0 / sqrt(lensq), vec_init(-Vh.y, Vh.x, 0.0)) : vec_init(1.0, 0.0, 0.0);
	VEC T2 = vec_cross(Vh, T1);
	double rad = sqrt(prng_db(s));
	double phi = 2.0 * PI * prng_db(s);
	double t1 = rad * cos(phi);
	double t2 = rad * sin(phi);
	double k = 0.5 * (1.0 + Vh.z);
	t2 = (1.0 - k) * sqrt(1.0 - sqr(t1)) + k * t2;
	VEC Nh = vec_add(vec_add(vec_scale(t1, T1), vec_scale(t2, T2)), vec_scale(sqrt(fmax(0.0, 1.0 - sqr(t1) - sqr(t2))), Vh));
	VEC Ne = vec_norm(vec_init(r * Nh.x, r * Nh.y, fmax(0.0, Nh.z)));
	Ne = matrix_trafo(m, Ne);
	return Ne;
}
bool bsdf_con(BSDF bsdf) {
	switch (bsdf.t) {
		case 0 : {
			return false;
		}
		case 1 : {
			return true;
		}
		case 2 : {
			return true;
		}
		case 3 : {
			return true;
		}
		case 4 : {
			return false;
		}
		case 5 : {
			return false;
		}
		default : {
			return false;
		}
	}
}
bool bsdf_eval_fb(BSDF bsdf, double nr, double nt, VEC n, VEC wi, VEC wo, double cosni, double cosno, VEC *fs, VEC *rf, VEC *rb, double *psf, double *psb) {
	switch (bsdf.t) {
		case 1 : {
			if (cosni * cosno > 0.0) {
				return false;
			}
			*fs = vec_scale(1.0 / PI, bsdf.a);
			*rf = bsdf.a;
			*rb = bsdf.a;
			*psf = 1.0 / PI;
			*psb = 1.0 / PI;
			return true;
		}
		case 2 : {
			if (cosni * cosno > 0.0) {
				return false;
			}
			VEC m = vec_norm(vec_sub(wo, wi));
			double cosnm = vec_dot(n, m), cosmi = vec_dot(m, wi);
			VEC fr = vec_init(fres_cond(bsdf.n.x / nr, bsdf.k.x, cosmi), fres_cond(bsdf.n.y / nr, bsdf.k.y, cosmi), fres_cond(bsdf.n.z / nr, bsdf.k.z, cosmi));
			double g1i = ggx_g1(bsdf.r, cosni), g1o = ggx_g1(bsdf.r, cosno);
			double d = ggx_d(bsdf.r, cosnm);
			*fs = vec_mul(vec_scale(g1i * g1o * d / fabs(4.0 * cosni * cosno), fr), bsdf.a);
			*rf = vec_mul(vec_scale(g1o, fr), bsdf.a);
			*rb = vec_mul(vec_scale(g1i, fr), bsdf.a);
			*psf = g1i * d / fabs(4.0 * cosni * cosno);
			*psb = g1o * d / fabs(4.0 * cosni * cosno);
			return true;
		}
		default : {
			return false;
		}
	}
}
VEC bsdf_sample_fb(BSDF bsdf, double nr, double nt, VEC n, VEC wi, PRNG_STATE *s, double cosni, bool rad, VEC *tp, VEC *rf, VEC *rb, double *psf, double *psb) {
	switch (bsdf.t) {
		case 0 : {
			*psf = 1.0 / fabs(cosni);
			*psb = 1.0 / fabs(cosni);
			return wi;
		}
		case 1 : {
			*tp = bsdf.a;
			*rf = bsdf.a;
			*rb = bsdf.a;
			*psf = 1.0 / PI;
			*psb = 1.0 / PI;
			double p = 2.0 * PI * prng_db(s);
			double cost = sqrt(prng_db(s)), sint = opst(cost), cosp = cos(p), sinp = sin(p);
			return matrix_trafo(matrix_tan(vec_scale(-sign(cosni), n)), pol2cartz1(cost, sint, cosp, sinp));;
		}
		case 2 : {
			VEC m = ggx_sample_vndf(n, bsdf.r, vec_scale(sign(cosni), wi), s);
			double cosmi = vec_dot(m, wi);
			VEC wo = refl(m, wi, cosmi);
			double cosno = vec_dot(n, wo);
			if (cosni * cosno > 0.0) {
				*rf = vec_init1(0.0);
				return wo;
			}
			VEC fr = vec_init(fres_cond(bsdf.n.x / nr, bsdf.k.x, cosmi), fres_cond(bsdf.n.y / nr, bsdf.k.y, cosmi), fres_cond(bsdf.n.z / nr, bsdf.k.z, cosmi));
			double g1i = ggx_g1(bsdf.r, cosni), g1o = ggx_g1(bsdf.r, cosno);
			*tp = vec_mul(vec_scale(g1o, fr), bsdf.a);
			*rf = vec_mul(vec_scale(g1o, fr), bsdf.a);
			*rb = vec_mul(vec_scale(g1i, fr), bsdf.a);
			double cosnm = vec_dot(n, m);
			double d = ggx_d(bsdf.r, cosnm);
			*psf = g1i * d / fabs(4.0 * cosni * cosno);
			*psb = g1o * d / fabs(4.0 * cosni * cosno);
			return wo;
		}
		case 4 : {
			VEC fr = vec_init(fres_cond(bsdf.n.x / nr, bsdf.k.x, cosni), fres_cond(bsdf.n.y / nr, bsdf.k.y, cosni), fres_cond(bsdf.n.z / nr, bsdf.k.z, cosni));
			*tp = vec_mul(fr, bsdf.a);
			*rf = vec_mul(fr, bsdf.a);
			*rb = vec_mul(fr, bsdf.a);
			*psf = 1.0 / fabs(cosni);
			*psb = 1.0 / fabs(cosni);
			return refl(n, wi, cosni);
		}
		case 5 : {
			double eta = nr / nt;
			double eta2 = sqr(eta);
			double det = 1.0 - eta2 * (1.0 - sqr(cosni));
			double fr = 1.0;
			if (det > 0.0) {
				double cosnt = sign(cosni) * sqrt(det);
				fr = fres_diel(eta, cosni, cosnt);
				if (prng_db(s) > fr) {
					*tp = vec_scale(rad ? eta2 : 1.0, bsdf.a);
					*rf = bsdf.a;
					*rb = bsdf.a;
					*psf = (1.0 - fr) / fabs(cosnt);
					*psb = (1.0 - fr) / fabs(cosni);
					return refr(n, wi, eta, cosni, cosnt);
				}
			}
			*tp = bsdf.a;
			*rf = bsdf.a;
			*rb = bsdf.a;
			*psf = fr / fabs(cosni);
			*psb = fr / fabs(cosni);
			return refl(n, wi, cosni);
		}
		default : {
			return vec_init1(0.0);
		}
	}
}
bool bsdf_insup(BSDF bsdf, double ni, double nt, VEC n, VEC wi, VEC wo) {
	switch (bsdf.t) {
		case 1 : {
			double cosi = vec_dot(n, wi), coso = vec_dot(n, wo);
			return cosi * coso < 0.0;
		}
		default : {
			return false;
		}
	}
}
double p1(double h) {
	return h > -1.0 && h < 1.0 ? 0.5 : 0.0;
}
double c1(double h) {
	return clamp(0.0, 1.0, (h + 1.0) / 2.0);
}
double c1i(double h) {
	return clamp(-1.0, 1.0, 2.0 * h - 1.0);
}
double ggx_g1h(double r, double cos, double h) {
	if (cos < EPS) {
		return 0.0;
	}
	return pow(c1(h), ggx_lam(r, cos));
}
bool sample_height(double r, double cosnr, double hr, PRNG_STATE *s, double *hrd) {
	if (fabs(cosnr) < EPS) {
		*hrd = hr;
		return false;
	}
	if (cosnr > 1.0 - EPS) {
		return true;
	}
	double u = prng_db(s);
	if (cosnr < EPS - 1.0) {
		*hrd = c1i(u * c1(hr));
		return false;
	}
	double g1h = ggx_g1h(r, cosnr, hr);
	if (u < g1h) {
		return true;
	}
	*hrd = c1i(c1(hr) / pow(1.0 - u, 1.0 / ggx_lam(r, cosnr)));
	return false;
}
VEC bsdf_sample(BSDF bsdf, double nr, double nt, VEC n, VEC wi, PRNG_STATE *s, double cosni, bool fow, VEC *tp, VEC *rf) {
	switch (bsdf.t) {
		case 0 : {
			return wi;
		}
		case 1 : {
			*tp = vec_mul(bsdf.a, *tp);
			*rf = vec_mul(bsdf.a, *rf);
			double p = 2.0 * PI * prng_db(s);
			double cost = sqrt(prng_db(s)), sint = opst(cost), cosp = cos(p), sinp = sin(p);
			return matrix_trafo(matrix_tan(vec_scale(-sign(cosni), n)), pol2cartz1(cost, sint, cosp, sinp));;
		}
		case 2 : {
			VEC m = ggx_sample_vndf(n, bsdf.r, vec_scale(sign(cosni), wi), s);
			double cosmi = vec_dot(m, wi);
			VEC wo = refl(m, wi, cosmi);
			double cosno = vec_dot(n, wo);
			if (cosni * cosno > 0.0) {
				*rf = vec_init1(0.0);
				return wo;
			}
			VEC fr = fres_cond_rgb(vec_scale(1.0 / nr, bsdf.n), bsdf.k, cosmi);
			double g1o = ggx_g1(bsdf.r, cosno);
			*tp = vec_mul(vec_mul(vec_scale(g1o, fr), bsdf.a), *tp);
			*rf = vec_mul(vec_mul(vec_scale(g1o, fr), bsdf.a), *rf);
			return wo;
		}
		case 4 : {
			VEC fr = fres_cond_rgb(vec_scale(1.0 / nr, bsdf.n), bsdf.k, cosni);
			*tp = vec_mul(vec_mul(fr, bsdf.a), *tp);
			*rf = vec_mul(vec_mul(fr, bsdf.a), *rf);
			return refl(n, wi, cosni);
		}
		case 5 : {
			double eta = nr / nt;
			double eta2 = sqr(eta);
			double det = 1.0 - eta2 * (1.0 - sqr(cosni));
			double fr = 1.0;
			if (det > 0.0) {
				double cosnt = sign(cosni) * sqrt(det);
				fr = fres_diel(eta, cosni, cosnt);
				if (prng_db(s) > fr) {
					*tp = vec_mul(vec_scale(fow ? eta2 : 1.0, bsdf.a), *tp);
					*rf = vec_mul(bsdf.a, *rf);
					return refr(n, wi, eta, cosni, cosnt);
				}
			}
			*tp = vec_mul(bsdf.a, *tp);
			*rf = vec_mul(bsdf.a, *rf);
			return refl(n, wi, cosni);
		}
		case 6 : {
			VEC wr = wi;
			double hr = c1i(1.0 - EPS);
			VEC e = vec_init1(1.0);
			while (true) {
				double cosnr = vec_dot(n, wr);
				if (sample_height(bsdf.r, cosnr, hr, s, &hr)) {
					*tp = vec_mul(e, *tp);
					*rf = vec_mul(e, *rf);
					return wr;
				}
				VEC m = ggx_sample_vndf(n, bsdf.r, vec_scale(sign(cosni), wr), s);
				double cosmi = vec_dot(m, wr);
				e = vec_mul(e, fres_cond_rgb(vec_scale(1.0 / nr, bsdf.n), bsdf.k, cosmi));
				wr = refl(m, wr, cosmi);
			}
		}
		case 7 : {
		}
		default : {
			return vec_init1(0.0);
		}
	}
}
VEC bsdf_sample_eye0(BSDF bsdf, double ni, double nt, VEC n, VEC wi, PRNG_STATE *s, VEC *tp, VEC *r) {
	switch (bsdf.t) {
		case 1 : {
			double phi = 2.0 * PI * prng_db(s);
			double cost, sint, cosp, sinp;
			cost = sqrt(prng_db(s));
			sint = opst(cost);
			cosp = cos(phi);
			sinp = sin(phi);
			n = vec_scale(-sign(vec_dot(n, wi)), n);
			MATRIX b = matrix_tan(n);
			VEC wo = matrix_trafo(b, pol2cartz1(cost, sint, cosp, sinp));
			*tp = vec_mul(*tp, bsdf.a);
			*r = vec_mul(*r, bsdf.a);
			return wo;
		}
		case 2 : {
			double cosni = vec_dot(n, wi);
			VEC m = ggx_sample_vndf(n, bsdf.r, vec_scale(sign(cosni), wi), s);
			double cosmi = vec_dot(m, wi);
			VEC wo = refl(m, wi, cosmi);
			double cosno = vec_dot(n, wo);
			if (cosni * cosno > 0.0) {
				*tp = vec_init1(0.0);
				*r = vec_init1(0.0);
				return wi;
			}
			VEC fr = vec_init(fres_cond(bsdf.n.x / ni, bsdf.k.x, cosmi), fres_cond(bsdf.n.y / ni, bsdf.k.y, cosmi), fres_cond(bsdf.n.z / ni, bsdf.k.z, cosmi));
			double g1o = ggx_g1(bsdf.r, cosno);
			*tp = vec_mul(vec_mul(vec_scale(g1o, fr), bsdf.a), *tp);
			*r = vec_mul(vec_mul(vec_scale(g1o, fr), bsdf.a), *r);
			return wo;
		}
		case 3 : {
			double cosni = vec_dot(n, wi);
			VEC m = ggx_sample_vndf(n, bsdf.r, vec_scale(sign(cosni), wi), s);
			double cosmi = vec_dot(m, wi);
			double eta = ni / nt;
			double eta2 = sqr(eta);
			double det = 1.0 - eta2 * (1.0 - sqr(cosmi));
			double fr = 1.0;
			double pr = 1.0;
			if (det > 0.0) {
				double cosmt = sign(cosmi) * sqrt(det);
				fr = fres_diel(eta, cosmi, cosmt);
				double ft = eta2 * (1.0 - fr);
				pr = fr / (fr + ft);
				double pt = 1.0 - pr;
				if (prng_db(s) < pt) {
					VEC wo = refr(m, wi, eta, cosmi, cosmt);
					double cosno = vec_dot(n, wo);
					if (cosni * cosno < 0.0) {
						*tp = vec_init1(0.0);
						*r = vec_init1(0.0);
						return wi;
					}
					double g1o = ggx_g1(bsdf.r, cosno);
					*tp = vec_mul(vec_scale(ft * g1o / pt, bsdf.a), *tp);
					*r = vec_mul(bsdf.a, *r);
					return wo;
				}
			}
			VEC wo = refl(m, wi, cosmi);
			double cosno = vec_dot(n, wo);
			if (cosni * cosno > 0.0) {
				*tp = vec_init1(0.0);
				*r = vec_init1(0.0);
				return wi;
			}
			double g1o = ggx_g1(bsdf.r, cosno);
			*tp = vec_mul(vec_scale(fr * g1o / pr, bsdf.a), *tp);
			*r = vec_mul(bsdf.a, *r);
			return wo;
			
		}
		case 4 : {
			double cosi = vec_dot(n, wi);
			VEC fr = vec_init(fres_cond(bsdf.n.x / ni, bsdf.k.x, cosi), fres_cond(bsdf.n.y / ni, bsdf.k.y, cosi), fres_cond(bsdf.n.z / ni, bsdf.k.z, cosi));
			*tp = vec_mul(vec_mul(fr, bsdf.a), *tp);
			*r = vec_mul(vec_mul(fr, bsdf.a), *r);
			return refl(n, wi, cosi);
		}
		case 5 : {
			double cosi = vec_dot(n, wi);
			double eta = ni / nt;
			double eta2 = sqr(eta);
			double det = 1.0 - eta2 * (1.0 - sqr(cosi));
			double fr = 1.0;
			if (det > 0.0) {
				double cost = sign(cosi) * sqrt(det);
				fr = fres_diel(eta, cosi, cost);
				if (prng_db(s) > fr) {
					*tp = vec_mul(*tp, bsdf.a);
					*r = vec_mul(*r, bsdf.a);
					return refr(n, wi, eta, cosi, cost);
				}
			}
			*tp = vec_mul(*tp, bsdf.a);
			*r = vec_mul(*r, bsdf.a);
			return refl(n, wi, cosi);
		}
		default : {
			return wi;
		}
	}
}
VEC bsdf_sample_eye(BSDF bsdf, double ni, double nt, VEC n, VEC wom, PRNG_STATE *s, VEC *tp, double *psf, double *psb) {
	switch (bsdf.t) {
		case 1 : {
			double u0, u1;
			u0 = prng_db(s);
			u1 = prng_db(s);
			double phi = 2.0 * PI * u1;
			double cost, sint, cosp, sinp;
			cost = sqrt(u0);
			sint = opst(cost);
			cosp = cos(phi);
			sinp = sin(phi);
			n = vec_scale(-sign(vec_dot(n, wom)), n);
			MATRIX b = matrix_tan(n);
			VEC wim = matrix_trafo(b, pol2cartz1(cost, sint, cosp, sinp));
			*tp = vec_mul(*tp, bsdf.a);
			*psf = 1.0 / PI;
			*psb = 1.0 / PI;
			return wim;
		}
		case 2 : {
			double cosi = vec_dot(n, wom);
			*tp = vec_mul(vec_init(fres_cond(bsdf.n.x / ni, bsdf.k.x, cosi), fres_cond(bsdf.n.y / ni, bsdf.k.y, cosi), fres_cond(bsdf.n.z / ni, bsdf.k.z, cosi)), *tp);
			*psf = 1.0 / fabs(cosi);
			*psb = 1.0 / fabs(cosi);
			return refl(n, wom, cosi);
		}
		case 3 : {
			double cosi = vec_dot(n, wom);
			double eta = ni / nt;
			double eta2 = sqr(eta);
			double det = 1.0 - eta2 * (1.0 - sqr(cosi));
			if (det > 0.0) {
				double cost = sign(cosi) * sqrt(det);
				double fr = fres_diel(eta, cosi, cost);
				double ft = eta2 * (1.0 - fr);
				double pr = fr / (fr + ft);
				double pt = 1.0 - pr;
				if (prng_db(s) < pt) {
					*tp = vec_mul(*tp, vec_scale((ft / pt), bsdf.a));
					*psf = pt / fabs(cost);
					*psb = pt / fabs(cosi);
					return refr(n, wom, eta, cosi, cost);
				} else {
					*tp = vec_mul(*tp, vec_scale(fr / pr, bsdf.a));
					*psf = pr / fabs(cosi);
					*psb = pr / fabs(cosi);
					return refl(n, wom, cosi);
				}
			} else {
				*tp = vec_mul(*tp, bsdf.a);
				*psf = 1.0 / fabs(cosi);
				*psb = 1.0 / fabs(cosi);
				return refl(n, wom, cosi);
			}
		}
		default : {
			double cos = vec_dot(n, wom);
			*tp = vec_mul(*tp, bsdf.a);
			*psf = 1.0 / fabs(cos);
			*psb = 1.0 / fabs(cos);
			return wom;
		}
	}
}
void medr(double cos, MEDIUM mp, MEDIUM mn, MEDIUM *m) {
	*m = cos < 0.0 ? mp : mn;
	return;
}
void medt(double cos, MEDIUM mp, MEDIUM mn, MEDIUM *m) {
	*m = cos < 0.0 ? mn : mp;
	return;
}
void medrt(double cos, MEDIUM mp, MEDIUM mn, MEDIUM *mr, MEDIUM *mt) {
	medr(cos, mp, mn, mr);
	medt(cos, mp, mn, mt);
	return;
}
int discrete_sample(int n, double *cdf, PRNG_STATE *s) {
	double x = prng_db(s);
	for(int i = 0; i < n; ++i) {
		if (x < cdf[i]) {
			return i;
		}
	}
	return -1;
}
double discrete_prob(double *cdf, int i) {
	return i == 0 ? cdf[0] : cdf[i] - cdf[i - 1];
}

VEC image_interp(IMAGE img, double x, double y) {
	double x0 = floor(x - 0.5) + 0.5;
	double y0 = floor(y - 0.5) + 0.5;
	double dx = x - x0;
	double dy = y - y0;
	double six0iy0 = (1.0 - dx) * (1.0 - dy);
	double six1iy0 = dx * (1.0 - dy);
	double six0iy1 = (1.0 - dx) * dy;
	double six1iy1 = dx * dy;
	int ix0 = floor(x0);
	int ix1 = ix0 + 1;
	int iy0 = floor(y0);
	int iy1 = iy0 + 1;
	ix0 = imax(0, ix0);
	ix1 = imin(img.w - 1, ix1);
	iy0 = imax(0, iy0);
	iy1 = imin(img.h - 1, iy1);
	VEC cix0iy0 = img.data[img.w * iy0 + ix0];
	VEC cix1iy0 = img.data[img.w * iy0 + ix1];
	VEC cix0iy1 = img.data[img.w * iy1 + ix0];
	VEC cix1iy1 = img.data[img.w * iy1 + ix1];
	VEC c = vec_init1(0.0);
	c = vec_add(c, vec_scale(six0iy0, cix0iy0));
	c = vec_add(c, vec_scale(six1iy0, cix1iy0));
	c = vec_add(c, vec_scale(six0iy1, cix0iy1));
	c = vec_add(c, vec_scale(six1iy1, cix1iy1));
	return c;
}

typedef struct conf_ibl CONF_IBL;

struct conf_ibl {
	int num_alloc_bvh_nodes;
	int num_threads;
	int sample_rate;
	int image_width;
	int image_height;
	double fx, fy;
	double lens_radius;
	double lens_flen;
	double len_path;
	MATRIX rot, roti;
	IMAGE env;
	double *cdfx, *cdfy;
};

void ibl_init_cdfs(IMAGE env, double **cdfx, double **cdfy) {
	double *pdf = ALLOC(double, env.w * env.h);
	double sum = 0.0;
	for (int iy = 0; iy < env.h; ++iy) {
		double w = sqr(0.5 * PI * (iy + 1) / env.h) - sqr(0.5 * PI * iy / env.h);
		for (int ix = 0; ix < env.w; ++ix) {
			double wl = w * rgb2lum(env.data[env.w * iy + ix]);
			pdf[env.w * iy + ix] = wl;
			sum += wl;
		}
	}
	for (int i = 0; i < env.w * env.h; ++i) {
		pdf[i] /= sum;
	}
	*cdfx = ALLOC(double, env.w);
	double sumx = 0.0;
	for (int ix = 0; ix < env.w; ++ix) {
		double sumy = 0.0;
		for (int iy = 0; iy < env.h; ++iy) {
			sumy += pdf[env.w * iy + ix];
		}
		sumx += sumy;
		(*cdfx)[ix] = sumx;
	}
	*cdfy = ALLOC(double, env.h);
	double sumy = 0.0;
	for (int iy = 0; iy < env.h; ++iy) {
		double sumx = 0.0;
		for (int ix = 0; ix < env.w; ++ix) {
			sumx += pdf[env.w * iy + ix];
		}
		sumy += sumx;
		(*cdfy)[iy] = sumy;
	}
	free(pdf);
	return;
	//IMAGE img = image_init(env.w, env.h);
	//for (int iy = 0; iy < img.h; ++iy) {
	//	for (int ix = 0; ix < img.w; ++ix) {
	//		img.data[img.w * iy + ix] = vec_init1(discrete_prob(*cdfx, ix) * discrete_prob(*cdfy, iy));
	//	}
	//}
	//image_write("pdfmul.pfm", img);
	//for (int i = 0; i < img.w * img.h; ++i) {
	//	img.data[i] = vec_init1(pdf[i]);
	//}
	//image_write("pdf.pfm", img);
	//puts("A");
}
double deg2rad(double x) {
	return PI * x / 180.0;
}
double focal_length(double d, double v) {
	return d * v / (d + v);
}
double focal_distance(double d, double f) {
	return d * f / (d - f);
}
CONF_IBL read_conf_ibl(char *path) {
	FILE *file = fopen(path, "r");
	if (file == NULL) {
		puts("file error");
		exit(0);
	}
	CONF_IBL conf;
	double fovh;
	char path_env[256];
	double p, t;
	double fdist;
	int nbuf = 256, nid = 256;
	char buf[nbuf], id[nid];
	while (true) {
		if (fgets(buf, nbuf, file) == NULL) {
			break;
		}
		sscanf(buf, "%s", id);
		if (strcmp(id, "num-alloc-bvh") == 0) {
			sscanf(buf, "%*s%d", &conf.num_alloc_bvh_nodes);
			continue;
		}
		if (strcmp(id, "num-threads") == 0) {
			sscanf(buf, "%*s%d", &conf.num_threads);
			continue;
		}
		if (strcmp(id, "image-width") == 0) {
			sscanf(buf, "%*s%d", &conf.image_width);
			continue;
		}
		if (strcmp(id, "image-height") == 0) {
			sscanf(buf, "%*s%d", &conf.image_height);
			continue;
		}
		if (strcmp(id, "sample-rate") == 0) {
			sscanf(buf, "%*s%d", &conf.sample_rate);
			continue;
		}
		if (strcmp(id, "fovh") == 0) {
			sscanf(buf, "%*s%lf", &fovh);
			continue;
		}
		if (strcmp(id, "lens-radius") == 0) {
			sscanf(buf, "%*s%lf", &conf.lens_radius);
			continue;
		}
		if (strcmp(id, "lens-fdist") == 0) {
			sscanf(buf, "%*s%lf", &fdist);
			continue;
		}
		if (strcmp(id, "path-env") == 0) {
			sscanf(buf, "%*s%s", path_env);
			continue;
		}
		if (strcmp(id, "env-pos") == 0) {
			sscanf(buf, "%*s%lf%lf", &t, &p);
			continue;
		}
		if (strcmp(id, "len-path") == 0) {
			sscanf(buf, "%*s%lf", &conf.len_path);
			continue;
		}
	}
	fclose(file);
	film_dim(conf.image_width, conf.image_height, fovh, &conf.fx, &conf.fy);
	conf.env = image_read(path_env);
	conf.rot = matrix_rot_y(deg2rad(t), deg2rad(p));
	conf.roti = matrix_rot_y_inv(deg2rad(t), deg2rad(p));
	conf.lens_flen = focal_length(1.0, fdist);
	return conf;
}

void conf_ibl_del(CONF_IBL conf) {
	image_del(conf.env);
	free(conf.cdfx);
	free(conf.cdfy);
	return;
}

void dir2img(int w, int h, MATRIX rot, VEC *a, double *x, double *y) {
	*a = matrix_trafo(rot, *a);
	double t = acos(clamp(-1.0, 1.0, a->y));
	double p = atan2(a->x, a->z) + (a->x < 0.0 ? 2.0 * PI : 0.0);
	*x = w * (2.0 * PI - p) / (2.0 * PI);
	*y = h * (PI - t) / PI;
	return;
}
VEC env_eval(IMAGE env, MATRIX rot, VEC wim) {
	double x, y;
	dir2img(env.w, env.h, rot, &wim, &x, &y);
	return image_interp(env, x, y);
}
double env_pdf(int nx, double *cdfx, int ny, double *cdfy, MATRIX rot, VEC wim) {
	double x, y;
	dir2img(nx, ny, rot, &wim, &x, &y);
	int ix = iclamp(0, nx - 1, x), iy = iclamp(0, ny - 1, y);
	return nx * ny * discrete_prob(cdfx, ix) * discrete_prob(cdfy, iy) / (2.0 * sqr(PI) * fmax(opst(wim.y), EPS));
}
VEC env_sample(IMAGE env, double *cdfx, double *cdfy, MATRIX roti, PRNG_STATE *s, VEC *tp, double *psf) {
	int ix = discrete_sample(env.w, cdfx, s);
	int iy = discrete_sample(env.h, cdfy, s);
	double x = ix + prng_db(s);
	double y = iy + prng_db(s);
	VEC c = image_interp(env, x, y);
	double t = PI * (1.0 - y / env.h), p = 2.0 * PI * (1.0 - x / env.w);
	VEC wim = pol2carty0(t, p);
	*psf = env.w * env.h * discrete_prob(cdfx, ix) * discrete_prob(cdfy, iy) / (2.0 * sqr(PI) * fmax(opst(wim.y), EPS));
	wim = matrix_trafo(roti, wim);
	*tp = vec_scale(1.0 / *psf, c);
	return vec_scale(-1.0, wim);
}

VEC transmittance(VEC a, double t) {
	return vec_exp(vec_scale(-t, a));
}
bool bsdf_kernel(VEC ng, VEC ns, VEC wi, VEC wo) {
	return vec_dot(ng, wi) * vec_dot(ns, wi) > 0.0 && vec_dot(ng, wo) * vec_dot(ns, wo) > 0.0;
}
VEC scene_traverse_ibl(SCENE scene, CONF_IBL conf, BVH_NODE *bvh_nodes, int *stack, PRNG_STATE *s, VEC x, VEC w, VEC tp) {
	bool dir = false;
	double psfc = 0.0;
	VEC I = vec_init1(0.0);
	for (int i = 0; i < conf.len_path; ++i) {
		VEC tuv = vec_init1(INFINITY);
		int id;
		if (!traverse(conf.lens_radius, scene.verts, scene.faces, bvh_nodes, stack, x, w, &tuv, &id)) {
			double misw = 1.0;
			if (dir) {
				misw = 1.0 / (1.0 + sqr(env_pdf(conf.env.w, conf.cdfx, conf.env.h, conf.cdfy, conf.rot, w) / psfc));
			}
			I = vec_add(I, vec_mul(tp, vec_scale(misw, env_eval(conf.env, conf.rot, w))));
			break;
		}
		x = vec_add(x, vec_scale(tuv.x, w));
		FACE face = scene.faces[id];
		MAT mat = scene.mats[face.im];
		BSDF bsdf = scene.bsdfs[mat.ib];
		MEDIUM mp = scene.meds[mat.imp], mn = scene.meds[mat.imn];
		VEC v0, v1, v2;
		face_get_verts(scene.verts, face, &v0, &v1, &v2);
		VEC ng = tri_normg(v0, v1, v2);
		VEC n0, n1, n2;
		face_get_norms(scene.norms, face, &n0, &n1, &n2);
		VEC ns = tri_norms(n0, n1, n2, tuv.y, tuv.z);
		MEDIUM mr, mt;
		medrt(vec_dot(ng, w), mp, mn, &mr, &mt);
		VEC tr = transmittance(mr.a, tuv.x);
		tp = vec_mul(vec_scale(fabs(vec_dot(ns, w) / vec_dot(ng, w)), tr), tp);
		double pr = rgb2lum(vec_mul(tr, bsdf.a));
		dir = bsdf_con(bsdf);
		if (dir && (i + 1 < conf.len_path)) {
			VEC tpd;
			double psfd;
			VEC wi = env_sample(conf.env, conf.cdfx, conf.cdfy, conf.roti, s, &tpd, &psfd);
			VEC wim = vec_scale(-1.0, wi);
			VEC tuv = vec_init1(INFINITY);
			if (bsdf_insup(bsdf, mr.n, mt.n, ns, wi, vec_scale(-1.0, w)) && bsdf_kernel(ng, ns, w, wim) && visible(conf.lens_radius, scene.verts, scene.faces, bvh_nodes, stack, x, wim, &tuv)) {
				double psf = bsdf_pdf_eye(bsdf, mr.n, mt.n, ns, w, wim);
				VEC fsc = vec_scale(fabs(vec_dot(ns, wi)), bsdf_eval(bsdf, mr.n, mt.n, ns, vec_scale(-1.0, w), wi));
				VEC f = vec_mul(vec_mul(tp, fsc), tpd);
				double misw = 1.0 / (1.0 + sqr(fabs(vec_dot(ns, wim)) * psf * pr / psfd));
				I = vec_add(I, vec_scale(misw, f));
			}
		}
		if (prng_db(s) >= pr) {
			break;
		}
		tp = vec_scale(1.0 / pr, tp);
		double psf, psb;
		VEC wd = bsdf_sample_eye(bsdf, mr.n, mt.n, ns, w, s, &tp, &psf, &psb);
		if (!bsdf_kernel(ng, ns, w, wd)) {
			break;
		}
		w = wd;
		psfc = fabs(vec_dot(ns, wd)) * psf * pr;
	}
	return I;
}
void sample_eye(int imgw, int imgh, double fx, double fy, double lensr, double flen, PRNG_STATE *st, double imgx, double imgy, VEC *tp, VEC *pos, VEC *dir) {
	*tp = vec_init1(imgw * imgh);
	*pos = disc_sample(lensr, st);
	*dir = vec_norm(vec_sub(vec_scale(-focal_distance(1.0, flen), trafo_img2film(imgw, imgh, fx, fy, imgx, imgy)), *pos));
	return;
}

typedef struct conf_pt CONF_PT;

struct conf_pt {
	int num_alloc_bvh_nodes;
	int num_threads;
	int sample_rate;
	int image_width;
	int image_height;
	double fx, fy;
	double lens_radius;
	double lens_flen;
	double len_path;
	MEDIUM m;
};

CONF_PT read_conf_pt(char *path) {
	FILE *file = fopen(path, "r");
	if (file == NULL) {
		puts("file error");
		exit(0);
	}
	CONF_PT conf;
	double fovh;
	double fdist;
	int nbuf = 256, nid = 256;
	char buf[nbuf], id[nid];
	while (true) {
		if (fgets(buf, nbuf, file) == NULL) {
			break;
		}
		sscanf(buf, "%s", id);
		if (strcmp(id, "num-alloc-bvh") == 0) {
			sscanf(buf, "%*s%d", &conf.num_alloc_bvh_nodes);
			continue;
		}
		if (strcmp(id, "num-threads") == 0) {
			sscanf(buf, "%*s%d", &conf.num_threads);
			continue;
		}
		if (strcmp(id, "image-width") == 0) {
			sscanf(buf, "%*s%d", &conf.image_width);
			continue;
		}
		if (strcmp(id, "image-height") == 0) {
			sscanf(buf, "%*s%d", &conf.image_height);
			continue;
		}
		if (strcmp(id, "sample-rate") == 0) {
			sscanf(buf, "%*s%d", &conf.sample_rate);
			continue;
		}
		if (strcmp(id, "fovh") == 0) {
			sscanf(buf, "%*s%lf", &fovh);
			continue;
		}
		if (strcmp(id, "lens-radius") == 0) {
			sscanf(buf, "%*s%lf", &conf.lens_radius);
			continue;
		}
		if (strcmp(id, "lens-fdist") == 0) {
			sscanf(buf, "%*s%lf", &fdist);
			continue;
		}
		if (strcmp(id, "len-path") == 0) {
			sscanf(buf, "%*s%lf", &conf.len_path);
			continue;
		}
		if (strcmp(id, "medium") == 0) {
			sscanf(buf, "%*s%lf%lf%lf%lf", &conf.m.n, &conf.m.a.x, &conf.m.a.y, &conf.m.a.z);
			continue;
		}
	}
	fclose(file);
	film_dim(conf.image_width, conf.image_height, fovh, &conf.fx, &conf.fy);
	conf.lens_flen = focal_length(1.0, fdist);
	return conf;
}
bool surf_kern1(double cosg, double coss) {
	return cosg * coss > 0.0 && fabs(cosg) > EPS && fabs(coss) > EPS;
}
VEC scene_traverse_pt(SCENE scene, CONF_PT conf, BVH_NODE *bvh_nodes, int *stack, PRNG_STATE *s, MEDIUM m, VEC x, VEC w, VEC tp) {
	VEC I = vec_init1(0.0);
	for (int i = 0; i < conf.len_path; ++i) {
		VEC tuv = vec_init1(INFINITY);
		int id;
		if (!traverse(conf.lens_radius, scene.verts, scene.faces, bvh_nodes, stack, x, w, &tuv, &id)) {
			break;
		}
		x = vec_add(x, vec_scale(tuv.x, w));
		EDF edf;
		BSDF bsdf;
		VEC ng, ns;
		MEDIUM mp, mn;
		if (id == -1) {
			edf.t = 0;
			bsdf.t = 0;
			ng = disc_normal();
			ns = disc_normal();
			mp = conf.m;
			mn = conf.m;
		} else {
			FACE face = scene.faces[id];
			MAT mat = scene.mats[face.im];
			edf = scene.edfs[mat.ie];
			bsdf = scene.bsdfs[mat.ib];
			mp = scene.meds[mat.imp];
			mn = scene.meds[mat.imn];
			VEC v0, v1, v2;
			face_get_verts(scene.verts, face, &v0, &v1, &v2);
			ng = tri_normg(v0, v1, v2);
			VEC n0, n1, n2;
			face_get_norms(scene.norms, face, &n0, &n1, &n2);
			ns = tri_norms(n0, n1, n2, tuv.y, tuv.z);
		}
		double cosgi = vec_dot(ng, w), cossi = vec_dot(ns, w);
		if (!surf_kern1(cosgi, cossi)) {
			break;
		}
		VEC tr = transmittance(m.a, tuv.x);
		tp = vec_mul(vec_scale(fabs(cossi / cosgi), tr), tp);
		VEC fs;
		if (edf_con(edf) && edf_eval(edf, -cossi, &fs)) {
			I = vec_add(I, vec_mul(tp, fs));
		}
		MEDIUM mr, mt;
		medrt(cosgi, mp, mn, &mr, &mt);
		VEC wo = bsdf_sample(bsdf, mr.n, mt.n, ns, w, s, cossi, true, &tp, &tr);
		double cosgo = vec_dot(ng, wo), cosso = vec_dot(ns, wo);
		if (!surf_kern1(cosgo, cosso)) {
			break;
		}
		double pr = rgb2lum(tr);
		if (prng_db(s) >= pr) {
			break;
		}
		tp = vec_scale(1.0 / pr, tp);
		w = wo;
		medt(cosgo, mp, mn, &m);
	}
	return I;
}
void render_pt(char *path_scene, char *path_conf) {
	SCENE scene = read_scene(path_scene);
	CONF_PT conf = read_conf_pt(path_conf);
	BVH_NODE *bvh_nodes = ALLOC(BVH_NODE, conf.num_alloc_bvh_nodes);
	int nbvh_nodes = bvh_build_wrap(&scene, conf.num_alloc_bvh_nodes, bvh_nodes);
	omp_set_num_threads(conf.num_threads);
	int *stacks[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		stacks[i] = ALLOC(int, nbvh_nodes);
	}
	IMAGE imgs[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		imgs[i] = image_init(conf.image_width, conf.image_height);
	}
	double t = time_get();
	#pragma omp parallel
	{
		int *stack = stacks[omp_get_thread_num()];
		IMAGE img = imgs[omp_get_thread_num()];
		PRNG_STATE s = prng_seed(SEED);
		for (int i = 0; i < omp_get_thread_num(); ++i) {
			s = prng_jump(s);
		}
		for (int i = 0; i < img.w * img.h; ++i) {
			int ix = i % img.w, iy = i / img.w;
			VEC sum = vec_init1(0.0);
			for (int j = 0; j < conf.sample_rate * conf.sample_rate; ++j) {
				int dx = j % conf.sample_rate, dy = j / conf.sample_rate;
				double x = ix + (dx + prng_db(&s)) / conf.sample_rate;
				double y = iy + (dy + prng_db(&s)) / conf.sample_rate;
				VEC tp, pos, dir;
				MEDIUM m = conf.m;
				sample_eye(img.w, img.h, conf.fx, conf.fy, conf.lens_radius, conf.lens_flen, &s, x, y, &tp, &pos, &dir);
				sum = vec_add(sum, scene_traverse_pt(scene, conf, bvh_nodes, stack, &s, m, pos, dir, tp));
			}
			img.data[i] = vec_scale(1.0 / (img.w * img.h * sqr(conf.sample_rate)), sum);
		}
	}
	t = time_get() - t;
	printf("rendered: %f\n", t);
	IMAGE img = image_init(conf.image_width, conf.image_height);
	image_fill(img, vec_init1(0.0));
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		for (int j = 0; j < img.w * img.h; ++j) {
			img.data[j] = vec_add(img.data[j], imgs[i].data[j]);
		}
	}
	for (int i = 0; i < img.w * img.h; ++i) {
		img.data[i] = vec_scale(1.0 / conf.num_threads, img.data[i]);
	}
	image_write("pt.pfm", img);
}
void render_ibl(char *path_scene, char *path_conf) {
	SCENE scene = read_scene(path_scene);
	CONF_IBL conf = read_conf_ibl(path_conf);
	ibl_init_cdfs(conf.env, &conf.cdfx, &conf.cdfy);
	BVH_NODE *bvh_nodes = ALLOC(BVH_NODE, conf.num_alloc_bvh_nodes);
	int nbvh_nodes = bvh_build_wrap(&scene, conf.num_alloc_bvh_nodes, bvh_nodes);
	omp_set_num_threads(conf.num_threads);
	int *stacks[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		stacks[i] = ALLOC(int, nbvh_nodes);
	}
	IMAGE imgs[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		imgs[i] = image_init(conf.image_width, conf.image_height);
	}
	double t = time_get();
	#pragma omp parallel
	{
		int *stack = stacks[omp_get_thread_num()];
		IMAGE img = imgs[omp_get_thread_num()];
		PRNG_STATE s = prng_seed(SEED);
		for (int i = 0; i < omp_get_thread_num(); ++i) {
			s = prng_jump(s);
		}
		for (int i = 0; i < img.w * img.h; ++i) {
			int ix = i % img.w, iy = i / img.w;
			VEC sum = vec_init1(0.0);
			for (int j = 0; j < conf.sample_rate * conf.sample_rate; ++j) {
				int dx = j % conf.sample_rate, dy = j / conf.sample_rate;
				double x = ix + (dx + prng_db(&s)) / conf.sample_rate;
				double y = iy + (dy + prng_db(&s)) / conf.sample_rate;
				VEC tp, pos, dir;
				sample_eye(img.w, img.h, conf.fx, conf.fy, conf.lens_radius, conf.lens_flen, &s, x, y, &tp, &pos, &dir);
				sum = vec_add(sum, scene_traverse_ibl(scene, conf, bvh_nodes, stack, &s, pos, dir, tp));
			}
			img.data[i] = vec_scale(1.0 / (img.w * img.h * sqr(conf.sample_rate)), sum);
		}
	}
	t = time_get() - t;
	printf("rendered: %f\n", t);
	IMAGE img = image_init(conf.image_width, conf.image_height);
	image_fill(img, vec_init1(0.0));
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		for (int j = 0; j < img.w * img.h; ++j) {
			img.data[j] = vec_add(img.data[j], imgs[i].data[j]);
		}
	}
	for (int i = 0; i < img.w * img.h; ++i) {
		img.data[i] = vec_scale(1.0 / conf.num_threads, img.data[i]);
	}
	image_write("ibl.pfm", img);
}

void render_debug(char *path_scene, char *path_conf) {
	SCENE scene = read_scene(path_scene);
	CONF_DEBUG conf = read_conf_debug(path_conf);
	BVH_NODE *bvh_nodes = ALLOC(BVH_NODE, conf.num_alloc_bvh_nodes);
	int nbvh_nodes = bvh_build_wrap(&scene, conf.num_alloc_bvh_nodes, bvh_nodes);
	IMAGE img = image_init(conf.image_width, conf.image_height);
	{
		image_fill(img, vec_init1(0.0));
		for (int i = 0; i < nbvh_nodes; ++i) {
			AABB box = bvh_nodes[i].b;
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
				int w = img.w, h = img.h;
				double fx = conf.fx, fy = conf.fy;
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
				int ix = imin(w - 1, a.x);
				int iy = imin(h - 1, a.y);
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
					img.data[ix + w * iy] = vec_add(img.data[ix + w * iy], vec_init1((tmax - tmin) * vec_len(d)));
					if (1.0 - tmax < EPS) {
						break;
					}
					ix += sx * px;
					iy += sy * py;
				}
			}
		}
	}
	image_write("bvh.pfm", img);
		image_fill(img, vec_init1(0.0));
		int *stack = ALLOC(int, nbvh_nodes);
		PRNG_STATE s = prng_seed(SEED);
		for (int i = 0; i < img.w * img.h; ++i) {
			int ix = i % img.w, iy = i / img.w;
			VEC sum = vec_init1(0.0);
			for (int j = 0; j < conf.sample_rate * conf.sample_rate; ++j) {
				int dx = j % conf.sample_rate, dy = j / conf.sample_rate;
				double x = ix + (dx + prng_db(&s)) / conf.sample_rate;
				double y = iy + (dy + prng_db(&s)) / conf.sample_rate;
				VEC o = vec_init1(0.0);
				VEC d = vec_norm(vec_scale(-1.0, trafo_img2film(img.w, img.h, conf.fx, conf.fy, x, y)));
				VEC tuv = vec_init1(INFINITY);
				int id;
				if (traverse(0.0, scene.verts, scene.faces, bvh_nodes, stack, o, d, &tuv, &id)) {
					VEC n0, n1, n2;
					face_get_norms(scene.norms, scene.faces[id], &n0, &n1, &n2);
					VEC ns = tri_norms(n0, n1, n2, tuv.y, tuv.z);
					sum = vec_add(sum, vec_abs(ns));
				}
			}
			img.data[img.w * iy + ix] = vec_scale(1.0 / sqr(conf.sample_rate), sum);
		}
		free(stack);
	image_write("tmp.pfm", img);
}
typedef struct bpt_vertex BPT_VERTEX;

struct bpt_vertex {
	VEC x, ng, ns, wi, tp, tr;
	double t, cossi, cosgo, paf, pab;
	int id;
	bool con, con0;
};

typedef struct conf_bpt CONF_BPT;

struct conf_bpt {
	int num_alloc_bvh_nodes;
	int num_threads;
	int sample_rate;
	int image_width, image_height;
	double fx, fy;
	double lens_radius;
	double lens_flen;
	int depth;
	int nls;
	int *ils;
	double *cdf;
	MEDIUM m;
};

void conf_bpt_del(CONF_BPT conf) {
	free(conf.ils);
	free(conf.cdf);
}
CONF_BPT read_conf_bpt(char *path) {
	FILE *file = fopen(path, "r");
	if (file == NULL) {
		puts("file error");
		exit(0);
	}
	CONF_BPT conf;
	double fovh;
	double fdist;
	int nbuf = 256, nid = 256;
	char buf[nbuf], id[nid];
	while (true) {
		if (fgets(buf, nbuf, file) == NULL) {
			break;
		}
		sscanf(buf, "%s", id);
		if (strcmp(id, "num-alloc-bvh") == 0) {
			sscanf(buf, "%*s%d", &conf.num_alloc_bvh_nodes);
		} else if (strcmp(id, "num-threads") == 0) {
			sscanf(buf, "%*s%d", &conf.num_threads);
		} else if (strcmp(id, "image-dim") == 0) {
			sscanf(buf, "%*s%d%d", &conf.image_width, &conf.image_height);
		} else if (strcmp(id, "sample-rate") == 0) {
			sscanf(buf, "%*s%d", &conf.sample_rate);
		} else if (strcmp(id, "fovh") == 0) {
			sscanf(buf, "%*s%lf", &fovh);
		} else if (strcmp(id, "lens-radius") == 0) {
			sscanf(buf, "%*s%lf", &conf.lens_radius);
		} else if (strcmp(id, "lens-fdist") == 0) {
			sscanf(buf, "%*s%lf", &fdist);
		} else if (strcmp(id, "depth") == 0) {
			sscanf(buf, "%*s%d", &conf.depth);
		} else if (strcmp(id, "medium") == 0) {
			sscanf(buf, "%*s%lf%lf%lf%lf", &conf.m.n, &conf.m.a.x, &conf.m.a.y, &conf.m.a.z);
		}
	}
	fclose(file);
	film_dim(conf.image_width, conf.image_height, fovh, &conf.fx, &conf.fy);
	conf.lens_flen = focal_length(1.0, fdist);
	return conf;
}
int init_cdf(SCENE scene, int **ils, double **cdf) {
	printf("proc init-cdf\n");
	double t = time_get();
	int nls = 0;
	for (int i = 0; i < scene.nfaces; ++i) {
		if (scene.edfs[scene.mats[scene.faces[i].im].ie].t > 0) {
			++nls;
		}
	}
	*ils = ALLOC(int, nls);
	*cdf = ALLOC(double, nls);
	int cnt = 0;
	double sum = 0.0;
	for (int i = 0; i < scene.nfaces; ++i) {
		FACE face = scene.faces[i];
		EDF edf = scene.edfs[scene.mats[face.im].ie];
		if (edf.t > 0) {
			VEC v0, v1, v2;
			face_get_verts(scene.verts, face, &v0, &v1, &v2);
			sum += rgb2lum(edf.a) * tri_area(v0, v1, v2);
			(*ils)[cnt] = i;
			(*cdf)[cnt] = sum;
			++cnt;
		}
	}
	for (int i = 0; i < nls; ++i) {
		(*cdf)[i] /= sum;
	}
	t = time_get() - t;
	printf("time %f\n", t);
	printf("num-lights %d\n", nls);
	return nls;
}

int bpt_traverse(SCENE scene, int depth, double lensr, MEDIUM me, BVH_NODE *bvh_nodes, int *stack, PRNG_STATE *s, bool type, MEDIUM m, VEC x, VEC w, VEC tp, double ptmp, BPT_VERTEX *path) {
	VEC rb;
	double psb;
	for (int i = 0; i < depth; ++i) {
		VEC tuv = vec_init1(INFINITY);
		int id;
		if (!traverse(lensr, scene.verts, scene.faces, bvh_nodes, stack, x, w, &tuv, &id)) {
			return i;
		}
		x = vec_add(x, vec_scale(tuv.x, w));
		BSDF bsdf;
		VEC ng, ns;
		MEDIUM mp, mn;
		if (id == -1) {
			bsdf.t = 0;
			ng = disc_normal();
			ns = disc_normal();
			mp = me;
			mn = me;
		} else {
			FACE face = scene.faces[id];
			MAT mat = scene.mats[face.im];
			bsdf = scene.bsdfs[mat.ib];
			mp = scene.meds[mat.imp];
			mn = scene.meds[mat.imn];
			VEC v0, v1, v2;
			face_get_verts(scene.verts, face, &v0, &v1, &v2);
			ng = tri_normg(v0, v1, v2);
			VEC n0, n1, n2;
			face_get_norms(scene.norms, face, &n0, &n1, &n2);
			ns = tri_norms(n0, n1, n2, tuv.y, tuv.z);
		}
		double cosgi = vec_dot(ng, w), cossi = vec_dot(ns, w);
		if (!surf_kern1(cosgi, cossi)) {
			return i;
		}
		VEC tr = transmittance(m.a, tuv.x);
		tp = vec_mul(vec_scale(fabs(cossi / cosgi), tr), tp);
		
		if (i > 0) {
			BPT_VERTEX *vm2 = path + i - 1, *vm1 = path + i;
			vm2->pab = rgb2lum(vec_mul(tr, rb)) * psb * fabs(vm1->cossi * vm2->cosgo) / sqr(vm1->t);
		}
		
		BPT_VERTEX v;
		v.x = x;
		v.ng = ng;
		v.ns = ns;
		v.wi = w;
		v.tp = tp;
		v.tr = tr;
		v.t = tuv.x;
		v.cossi = cossi;
		v.paf = ptmp * fabs(cosgi) / sqr(tuv.x);
		v.id = id;
		v.con = bsdf_con(bsdf);
		v.con0 = type ? (id == -1 ? false : edf_con(scene.edfs[scene.mats[scene.faces[id].im].ie])) : (id == -1 ? true : false);
		
		MEDIUM mr, mt;
		medrt(cosgi, mp, mn, &mr, &mt);
		VEC tpd, rf;
		double psf;
		VEC wo = bsdf_sample_fb(bsdf, mr.n, mt.n, ns, w, s, cossi, type, &tpd, &rf, &rb, &psf, &psb);
		double cosgo = vec_dot(ng, wo), cosso = vec_dot(ns, wo);
		
		v.cosgo = cosgo;
		path[i + 1] = v;
		
		if (!surf_kern1(cosgo, cosso)) {
			return i + 1;
		}
		double pr = rgb2lum(vec_mul(tr, rf));
		if (prng_db(s) >= pr) {
			return i + 1;
		}
		
		tp = vec_scale(1.0 / pr, vec_mul(tp, tpd));
		ptmp = psf * fabs(cosso) * pr;
		w = wo;
		medt(cosgo, mp, mn, &m);
	}
	return depth;
}

int search(int nis, int *is, int x) {
	int a = 0, b = nis;
	while (a < b) {
		int c = (a + b) / 2;
		if (is[c] < x) {
			a = c + 1;
		} else if (is[c] > x) {
			b = c;
		} else {
			return c;
		}
	}
	return -1;
}
BPT_VERTEX bpt_sample_eye(double lensr, double flen, int imgw, int imgh, double fx, double fy, PRNG_STATE *s, double imgx, double imgy, VEC *tp, VEC *x, VEC *wo, double *ptmp) {
	BPT_VERTEX vert;
	vert.id = -1;
	*x = disc_sample(lensr, s);
	vert.x = *x;
	vert.ng = disc_normal();
	vert.ns = disc_normal();
	vert.paf = disc_pdf(lensr);
	vert.tp = vec_init1(1.0 / vert.paf);
	*tp = vec_init1(imgw * imgh);
	*wo = vec_norm(vec_sub(vec_scale(-focal_distance(1.0, flen), trafo_img2film(imgw, imgh, fx, fy, imgx, imgy)), *x));
	vert.cosgo = wo->z;
	*ptmp = 1.0 / (4.0 * fx * fy * cube(fabs(wo->z)));
	return vert;
}
int bpt_gen_eyepath(SCENE scene, int depth, double lensr, double flen, int imgw, int imgh, double fx, double fy, MEDIUM me, BVH_NODE *bvh_nodes, int *bvh_stack, PRNG_STATE *s, double imgx, double imgy, BPT_VERTEX *path) {
	VEC tp, x, w;
	double ptmp;
	MEDIUM m = me;
	path[0] = bpt_sample_eye(lensr, flen, imgw, imgh, fx, fy, s, imgx, imgy, &tp, &x, &w, &ptmp);
	return bpt_traverse(scene, depth, lensr, m, bvh_nodes, bvh_stack, s, true, m, x, w, tp, ptmp, path) + 1;
}
//struct bpt_vertex {
//	VEC x, ng, ns, wi, tp, tr;
//	double t, cossi, cosgo, paf, pab;
//	int id;
//};

BPT_VERTEX bpt_sample_light(SCENE scene, int nls, int *ils, double *cdf, PRNG_STATE *s, VEC *tp, VEC *x, VEC *wo, double *ptmp, bool *kern, MEDIUM *m) {
	BPT_VERTEX vert;
	vert.id = ils[discrete_sample(nls, cdf, s)];
	FACE face = scene.faces[vert.id];
	VEC v0, v1, v2;
	face_get_verts(scene.verts, face, &v0, &v1, &v2);
	vert.ng = tri_normg(v0, v1, v2);
	double u, v;
	tri_sample(s, &u, &v);
	*x = tri_point(v0, v1, v2, u, v);
	vert.x = *x;
	VEC n0, n1, n2;
	face_get_norms(scene.norms, face, &n0, &n1, &n2);
	vert.ns = tri_norms(n0, n1, n2, u, v);
	vert.paf = discrete_prob(cdf, search(nls, ils, vert.id)) * tri_pdf(v0, v1, v2);
	*tp = vec_init1(1.0 / vert.paf);
	vert.tp = *tp;
	MAT mat = scene.mats[face.im];
	EDF edf = scene.edfs[mat.ie];
	*wo = edf_sample_f(edf, vert.ns, s, tp, ptmp);
	vert.cosgo = vec_dot(vert.ng, *wo);
	vert.con = edf_con(edf);
	double cosso = vec_dot(vert.ns, *wo);
	*ptmp = fabs(cosso) * *ptmp;
	*kern = surf_kern1(vert.cosgo, cosso);
	medt(vert.cosgo, scene.meds[mat.imp], scene.meds[mat.imn], m);
	return vert;
}

int bpt_gen_lightpath(SCENE scene, int depth, double lensr, MEDIUM me, int nls, int *ils, double *cdf, BVH_NODE *bvh_nodes, int *bvh_stack, PRNG_STATE *s, BPT_VERTEX *path) {
	VEC tp, x, w;
	double ptmp;
	bool kern;
	MEDIUM m;
	path[0] = bpt_sample_light(scene, nls, ils, cdf, s, &tp, &x, &w, &ptmp, &kern, &m);
	if (!kern) {
		return 1;
	}
	return bpt_traverse(scene, depth, lensr, me, bvh_nodes, bvh_stack, s, false, m, x, w, tp, ptmp, path) + 1;
}

bool wdf_eval_b(double lensr, double flen, int imgw, int imgh, double fx, double fy, VEC x, VEC xd, VEC *fs, double *psb, int *imgx, int *imgy) {
	if (xd.z > 0.0) {
		return false;
	}
	VEC d;
	if (fabs(xd.z + flen) > EPS) {
		d = vec_sub(vec_scale(focal_distance(-xd.z, flen) / xd.z, xd), x);
		d = vec_add(x, vec_scale(1.0 / d.z, d));
	} else {
		d = vec_add(x, vec_scale(1.0 / xd.z, xd));
	}
	if (fabs(d.x) > fx || fabs(d.y) > fy) {
		return false;
	}
	VEC v = vec_sub(x, xd);
	double cos4 = quad(v.z) / sqr(vec_len2(v));
	*fs = vec_init1(imgw * imgh / (PI * sqr(lensr) * 4.0 * fx * fy * cos4));
	*psb = 1.0 / (4.0 * fx * fy * cos4);
	trafo_film2img(imgw, imgh, fx, fy, d.x, d.y, imgx, imgy);
	return true;
}
bool eye_eval_b(double lensr, double flen, int imgw, int imgh, double fx, double fy, VEC x, VEC xd, VEC *fs, double *pab, double *psb, int *imgx, int *imgy) {
	if (!wdf_eval_b(lensr, flen, imgw, imgh, fx, fy, x, xd, fs, psb, imgx, imgy)) {
		return false;
	}
	*pab = disc_pdf(lensr);
	return true;
}
bool light_eval_f(SCENE scene, int nls, int *ils, double *cdf, int id, double cosno, VEC *fs, double *paf, double *psf) {
	if (!edf_eval_f(scene.edfs[scene.mats[scene.faces[id].im].ie], cosno, fs, psf)) {
		return false;
	}
	VEC v0, v1, v2;
	face_get_verts(scene.verts, scene.faces[id], &v0, &v1, &v2);
	*paf = discrete_prob(cdf, search(nls, ils, id)) * tri_pdf(v0, v1, v2);
	return true;
}

void bpt_vc(SCENE scene, BVH_NODE *bvh_nodes, int *bvh_stack, CONF_BPT conf, IMAGE img, int ip, int N, int nl, BPT_VERTEX *pl, int ne, BPT_VERTEX *pe) {
	VEC sum0 = vec_init1(0.0);
	
	//s = 1, t = 1
	for (int i = 0; i < 1; ++i) {
		if (1 > conf.depth) {
			break;
		}
		BPT_VERTEX ve = pe[i];
		if (!ve.con) {
			break;
		}
		for (int j = 0; j < 1; ++j) {
			BPT_VERTEX vl = pl[j];
			if (!vl.con) {
				break;
			}
			
			VEC d = vec_sub(ve.x, vl.x);
			double t = vec_len(d);
			d = vec_scale(1.0 / t, d);
			
			double coslgo = vec_dot(vl.ng, d), coslso = vec_dot(vl.ns, d);
			if (!surf_kern1(coslgo, coslso)) {
				continue;
			}
			MAT ml = scene.mats[scene.faces[vl.id].im];
			MEDIUM mlp = scene.meds[ml.imp], mln = scene.meds[ml.imn];
			EDF dfl = scene.edfs[ml.ie];
			VEC fsl;
			double pseb0;
			if (!edf_eval_f(dfl, coslso, &fsl, &pseb0)) {
				continue;
			}
			
			double cosegi = vec_dot(ve.ng, d), cosesi = vec_dot(ve.ns, d);
			if (!surf_kern1(cosegi, cosesi)) {
				continue;
			}
			VEC fse;
			double pslb0;
			int ix, iy;
			if (!wdf_eval_b(conf.lens_radius, conf.lens_flen, img.w, img.h, conf.fx, conf.fy, ve.x, vl.x, &fse, &pslb0, &ix, &iy)) {
				continue;
			}
			
			VEC tuv = vec_init1(t);
			if (!visible(conf.lens_radius, scene.verts, scene.faces, bvh_nodes, bvh_stack, vl.x, d, &tuv)) {
				continue;
			}
			
			MEDIUM m;
			medt(coslso, mlp, mln, &m);
			
			VEC tr = transmittance(m.a, t);
			
			double palb0 = pslb0 * fabs(cosesi * coslgo) / sqr(t);
			double paeb0 = pseb0 * fabs(coslso * cosegi) / sqr(t);
			
			double w = 1.0;
			double r;
			
			r = 1.0;
			r *= palb0 / vl.paf;
			w += sqr(r);
			
			r = 1.0;
			r *= paeb0 / ve.paf;
			w += sqr(r);
			
			w = 1.0 / w;
			
			img.data[img.w * iy + ix] = vec_add(img.data[img.w * iy + ix], vec_scale(w * fabs(coslso * cosesi) / sqr(t) / N, vec_mul(vec_mul(vec_mul(fsl, vl.tp), vec_mul(fse, ve.tp)), tr)));
		}
	}
	
	//s >= 2, t = 1
	for (int j = 1; j < nl; ++j) {
		if (j + 1 > conf.depth) {
			break;
		}
		BPT_VERTEX ve = pe[0];
		if (!ve.con) {
			break;
		}
		BPT_VERTEX vl = pl[j];
		if (!vl.con) {
			continue;
		}
		
		VEC d = vec_sub(ve.x, vl.x);
		double t = vec_len(d);
		d = vec_scale(1.0 / t, d);
		
		double coslgo = vec_dot(vl.ng, d), coslso = vec_dot(vl.ns, d);
		if (!surf_kern1(coslgo, coslso)) {
			continue;
		}
		MAT ml = scene.mats[scene.faces[vl.id].im];
		MEDIUM mlp = scene.meds[ml.imp], mln = scene.meds[ml.imn], mlr, mlt;
		medrt(vl.cossi, mlp, mln, &mlr, &mlt);
		BSDF dfl = scene.bsdfs[ml.ib];
		VEC fsl, reb0, rlb1;
		double pseb0, pslb1;
		if (!bsdf_eval_fb(dfl, mlr.n, mlt.n, vl.ns, vl.wi, d, vl.cossi, coslso, &fsl, &reb0, &rlb1, &pseb0, &pslb1)) {
			continue;
		}
		
		double cosegi = vec_dot(ve.ng, d), cosesi = vec_dot(ve.ns, d);
		if (!surf_kern1(cosegi, cosesi)) {
			continue;
		}
		VEC fse;
		double pslb0;
		int ix, iy;
		if (!wdf_eval_b(conf.lens_radius, conf.lens_flen, img.w, img.h, conf.fx, conf.fy, ve.x, vl.x, &fse, &pslb0, &ix, &iy)) {
			continue;
		}
		
		VEC tuv = vec_init1(t);
		if (!visible(conf.lens_radius, scene.verts, scene.faces, bvh_nodes, bvh_stack, vl.x, d, &tuv)) {
			continue;
		}
		
		MEDIUM m;
		medt(coslso, mlp, mln, &m);
		
		VEC tr = transmittance(m.a, t);
		
		BPT_VERTEX vlm1 = pl[j - 1];
		double palb0 = pslb0 * fabs(cosesi * coslgo) / sqr(t);
		double palb1 = rgb2lum(vec_mul(tr, rlb1)) * pslb1 * fabs(vl.cossi * vlm1.cosgo) / sqr(vl.t);
		
		double paeb0 = rgb2lum(vec_mul(vl.tr, reb0)) * pseb0 * fabs(coslso * ve.cosgo) / sqr(t);
		
		double w = 1.0;
		double r;
		
		r = 1.0;
		r *= palb0 / vl.paf;
		if (vlm1.con) {
			w += sqr(r);
		}
		r *= palb1 / vlm1.paf;
		if (vlm1.con && (j > 1 ? pl[j - 2].con : true)) {
			w += sqr(r);
		}
		for (int k = j - 2; k >= 0; --k) {
			BPT_VERTEX vl = pl[k];
			r *= vl.pab / vl.paf;
			if (vl.con && (k > 0 ? pl[k - 1].con : true)) {
				w += sqr(r);
			}
		}
		
		r = 1.0;
		r *= paeb0 / ve.paf;
		w += sqr(r);
		
		w = 1.0 / w;
		
		img.data[img.w * iy + ix] = vec_add(img.data[img.w * iy + ix], vec_scale(w * fabs(coslso * cosesi) / sqr(t) / N, vec_mul(vec_mul(vec_mul(fsl, vl.tp), vec_mul(fse, ve.tp)), tr)));
	}
	
	//s >= 2, t = 0
	for (int j = 1; j < nl; ++j) {
		BPT_VERTEX vl = pl[j];
		if (!vl.con0) {
			continue;
		}
		BPT_VERTEX vlm1 = pl[j - 1];
		VEC fsl;
		double palb0, pslb1;
		int ix, iy;
		if (!eye_eval_b(conf.lens_radius, conf.lens_flen, img.w, img.h, conf.fx, conf.fy, vl.x, vlm1.x, &fsl, &palb0, &pslb1, &ix, &iy)) {
			continue;
		}
		double palb1 = pslb1 * fabs(vl.cossi * vlm1.cosgo) / sqr(vl.t);
		double w = 1.0;
		double r;
		r = 1.0;
		r *= palb0 / vl.paf;
		w += vlm1.con ? sqr(r) : 0.0;
		r *= palb1 / vlm1.paf;
		w += (vlm1.con && (j > 1 ? pl[j - 2].con : true)) ? sqr(r) : 0.0;
		for (int k = j - 2; k >= 0; --k) {
			BPT_VERTEX vl = pl[k];
			r *= vl.pab / vl.paf;
			w += (vl.con && (k > 0 ? pl[k - 1].con : true)) ? sqr(r) : 0.0;
		}
		w = 1.0 / w;
		img.data[img.w * iy + ix] = vec_add(img.data[img.w * iy + ix], vec_scale(w / N, vec_mul(fsl, vl.tp)));
	}
	
	//s = 0, t >= 2
	for (int i = 1; i < ne; ++i) {
		BPT_VERTEX ve = pe[i];
		if (!ve.con0) {
			continue;
		}
		VEC fse;
		double paeb0, pseb1;
		if (!light_eval_f(scene, conf.nls, conf.ils, conf.cdf, ve.id, -ve.cossi, &fse, &paeb0, &pseb1)) {
			continue;
		}
		BPT_VERTEX vem1 = pe[i - 1];
		double paeb1 = pseb1 * fabs(ve.cossi * vem1.cosgo) / sqr(ve.t);
		
		double w = 1.0;
		double r;
		
		r = 1.0;
		r *= paeb0 / ve.paf;
		if (vem1.con) {
			w += sqr(r);
		}
		r *= paeb1 / vem1.paf;
		if (vem1.con && (i > 1 ? pe[i - 2].con : true)) {
			w += sqr(r);
		}
		for (int k = i - 2; k >= 0; --k) {
			BPT_VERTEX ve = pe[k];
			r *= ve.pab / ve.paf;
			if (ve.con && (k > 0 ? pe[k - 1].con : true)) {
				w += sqr(r);
			}
		}
		
		w = 1.0 / w;
		
		sum0 = vec_add(sum0, vec_scale(w, vec_mul(fse, ve.tp)));
	}
	
	//s = 1, t >= 2
	for (int i = 1; i < ne; ++i) {
		if (i + 1 > conf.depth) {
			break;
		}
		BPT_VERTEX ve = pe[i];
		if (!ve.con) {
			continue;
		}
		BPT_VERTEX vl = pl[0];
		if (!vl.con) {
			break;
		}
			
		VEC d = vec_sub(ve.x, vl.x);
		double t = vec_len(d);
		d = vec_scale(1.0 / t, d);
		
		double coslgo = vec_dot(vl.ng, d), coslso = vec_dot(vl.ns, d);
		if (!surf_kern1(coslgo, coslso)) {
			continue;
		}
		MAT ml = scene.mats[scene.faces[vl.id].im];
		EDF dfl = scene.edfs[ml.ie];
		VEC fsl;
		double pseb0;
		if (!edf_eval_f(dfl, coslso, &fsl, &pseb0)) {
			continue;
		}
		
		double cosegi = vec_dot(ve.ng, d), cosesi = vec_dot(ve.ns, d);
		if (!surf_kern1(cosegi, cosesi)) {
			continue;
		}
		MAT me = scene.mats[scene.faces[ve.id].im];
		MEDIUM mep = scene.meds[me.imp], men = scene.meds[me.imn], mer, met;
		medrt(cosesi, mep, men, &mer, &met);
		BSDF dfe = scene.bsdfs[me.ib];
		VEC fse, rlb0, reb1;
		double pslb0, pseb1;
		if (!bsdf_eval_fb(dfe, mer.n, met.n, ve.ns, d, vec_scale(-1.0, ve.wi), cosesi, -ve.cossi, &fse, &reb1, &rlb0, &pseb1, &pslb0)) {
			continue;
		}
		VEC tuv = vec_init1(t);
		if (!visible(conf.lens_radius, scene.verts, scene.faces, bvh_nodes, bvh_stack, vl.x, d, &tuv)) {
			continue;
		}
		
		MEDIUM m;
		medr(cosesi, mep, men, &m);
		
		VEC tr = transmittance(m.a, t);
		
		double palb0 = rgb2lum(vec_mul(ve.tr, rlb0)) * pslb0 * fabs(cosesi * coslgo) / sqr(t);
		
		BPT_VERTEX vem1 = pe[i - 1];
		double paeb0 = pseb0 * fabs(coslso * cosegi) / sqr(t);
		double paeb1 = rgb2lum(vec_mul(tr, reb1)) * pseb1 * fabs(ve.cossi * vem1.cosgo) / sqr(ve.t);
		
		double w = 1.0;
		double r;
		
		r = 1.0;
		r *= palb0 / vl.paf;
		w += sqr(r);
		
		r = 1.0;
		r *= paeb0 / ve.paf;
		if (vem1.con) {
			w += sqr(r);
		}
		r *= paeb1 / vem1.paf;
		if (vem1.con && (i > 1 ? pe[i - 2].con : true)) {
			w += sqr(r);
		}
		for (int k = i - 2; k >= 0; --k) {
			BPT_VERTEX ve = pe[k];
			r *= ve.pab / ve.paf;
			if (ve.con && (k > 0 ? pe[k - 1].con : true)) {
				w += sqr(r);
			}
		}
		
		w = 1.0 / w;
		
		sum0 = vec_add(sum0, vec_scale(w * fabs(coslso * cosesi) / sqr(t), vec_mul(vec_mul(vec_mul(fsl, vl.tp), vec_mul(fse, ve.tp)), tr)));
	}
	
	//s >= 2, t >= 2
	for (int i = 1; i < ne; ++i) {
		if (i + 2 > conf.depth) {
			break;
		}
		BPT_VERTEX ve = pe[i];
		if (!ve.con) {
			continue;
		}
		for (int j = 1; j < nl; ++j) {
			if (i + j + 1 > conf.depth) {
				break;
			}
			BPT_VERTEX vl = pl[j];
			if (!vl.con) {
				continue;
			}
			
			VEC d = vec_sub(ve.x, vl.x);
			double t = vec_len(d);
			d = vec_scale(1.0 / t, d);
			
			double coslgo = vec_dot(vl.ng, d), coslso = vec_dot(vl.ns, d);
			if (!surf_kern1(coslgo, coslso)) {
				continue;
			}
			MAT ml = scene.mats[scene.faces[vl.id].im];
			MEDIUM mlp = scene.meds[ml.imp], mln = scene.meds[ml.imn], mlr, mlt;
			medrt(vl.cossi, mlp, mln, &mlr, &mlt);
			BSDF dfl = scene.bsdfs[ml.ib];
			VEC fsl, reb0, rlb1;
			double pseb0, pslb1;
			if (!bsdf_eval_fb(dfl, mlr.n, mlt.n, vl.ns, vl.wi, d, vl.cossi, coslso, &fsl, &reb0, &rlb1, &pseb0, &pslb1)) {
				continue;
			}
			
			double cosegi = vec_dot(ve.ng, d), cosesi = vec_dot(ve.ns, d);
			if (!surf_kern1(cosegi, cosesi)) {
				continue;
			}
			MAT me = scene.mats[scene.faces[ve.id].im];
			MEDIUM mep = scene.meds[me.imp], men = scene.meds[me.imn], mer, met;
			medrt(cosesi, mep, men, &mer, &met);
			BSDF dfe = scene.bsdfs[me.ib];
			VEC fse, rlb0, reb1;
			double pslb0, pseb1;
			if (!bsdf_eval_fb(dfe, mer.n, met.n, ve.ns, d, vec_scale(-1.0, ve.wi), cosesi, -ve.cossi, &fse, &reb1, &rlb0, &pseb1, &pslb0)) {
				continue;
			}
			
			VEC tuv = vec_init1(t);
			if (!visible(conf.lens_radius, scene.verts, scene.faces, bvh_nodes, bvh_stack, vl.x, d, &tuv)) {
				continue;
			}
			
			MEDIUM m;
			medt(coslso, mlp, mln, &m);
			
			VEC tr = transmittance(m.a, t);
			
			BPT_VERTEX vlm1 = pl[j - 1];
			double palb0 = rgb2lum(vec_mul(ve.tr, rlb0)) * pslb0 * fabs(cosesi * coslgo) / sqr(t);
			double palb1 = rgb2lum(vec_mul(tr, rlb1)) * pslb1 * fabs(vl.cossi * vlm1.cosgo) / sqr(vl.t);
			
			BPT_VERTEX vem1 = pe[i - 1];
			double paeb0 = rgb2lum(vec_mul(vl.tr, reb0)) * pseb0 * fabs(coslso * cosegi) / sqr(t);
			double paeb1 = rgb2lum(vec_mul(tr, reb1)) * pseb1 * fabs(ve.cossi * vem1.cosgo) / sqr(ve.t);
			
			double w = 1.0;
			double r;
			
			r = 1.0;
			r *= palb0 / vl.paf;
			w += vlm1.con ? sqr(r) : 0.0;
			r *= palb1 / vlm1.paf;
			w += (vlm1.con && (j > 1 ? pl[j - 2].con : true)) ? sqr(r) : 0.0;
			for (int k = j - 2; k >= 0; --k) {
				BPT_VERTEX vl = pl[k];
				r *= vl.pab / vl.paf;
				w += (vl.con && (k > 0 ? pl[k - 1].con : true)) ? sqr(r) : 0.0;
			}
			
			r = 1.0;
			r *= paeb0 / ve.paf;
			w += vem1.con ? sqr(r) : 0.0;
			r *= paeb1 / vem1.paf;
			w += (vem1.con && (i > 1 ? pe[i - 2].con : true)) ? sqr(r) : 0.0;
			for (int k = i - 2; k >= 0; --k) {
				BPT_VERTEX ve = pe[k];
				r *= ve.pab / ve.paf;
				w += (ve.con && (k > 0 ? pe[k - 1].con : true)) ? sqr(r) : 0.0;
			}
			w = 1.0 / w;
			sum0 = vec_add(sum0, vec_scale(w * fabs(coslso * cosesi) / sqr(t), vec_mul(vec_mul(vec_mul(fsl, vl.tp), vec_mul(fse, ve.tp)), tr)));
		}
	}
	img.data[ip] = vec_add(img.data[ip], vec_scale(1.0 / N, sum0));
}
void render_bpt(char *path_scene, char *path_conf) {
	printf("proc render-bpt\n");
	double t = time_get();
	SCENE scene = read_scene(path_scene);
	CONF_BPT conf = read_conf_bpt(path_conf);
	BVH_NODE *bvh_nodes = ALLOC(BVH_NODE, conf.num_alloc_bvh_nodes);
	int nbvh_nodes = bvh_build_wrap(&scene, conf.num_alloc_bvh_nodes, bvh_nodes);
	conf.nls = init_cdf(scene, &conf.ils, &conf.cdf);
	omp_set_num_threads(conf.num_threads);
	IMAGE imgs[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		imgs[i] = image_init(conf.image_width, conf.image_height);
	}
	#pragma omp parallel
	{
		int *stack = ALLOC(int, nbvh_nodes);
		BPT_VERTEX *pe = ALLOC(BPT_VERTEX, conf.depth + 1);
		BPT_VERTEX *pl = ALLOC(BPT_VERTEX, conf.depth + 1);
		IMAGE img = imgs[omp_get_thread_num()];
		image_fill(img, vec_init1(0.0));
		PRNG_STATE s = prng_seed(SEED);
		for (int i = 0; i < omp_get_thread_num(); ++i) {
			s = prng_jump(s);
		}
		for (int i = 0; i < img.w * img.h; ++i) {
			int ix = i % img.w, iy = i / img.w;
			for (int j = 0; j < conf.sample_rate * conf.sample_rate; ++j) {
				int dx = j % conf.sample_rate, dy = j / conf.sample_rate;
				double x = ix + (dx + prng_db(&s)) / conf.sample_rate, y = iy + (dy + prng_db(&s)) / conf.sample_rate;
				int ne = bpt_gen_eyepath(scene, conf.depth, conf.lens_radius, conf.lens_flen, conf.image_width, conf.image_height, conf.fx, conf.fy, conf.m, bvh_nodes, stack, &s, x, y, pe);
				int nl = bpt_gen_lightpath(scene, conf.depth, conf.lens_radius, conf.m,conf.nls, conf.ils, conf.cdf, bvh_nodes, stack, &s, pl);
				bpt_vc(scene, bvh_nodes, stack, conf, img, i, img.w * img.h * conf.sample_rate * conf.sample_rate, nl, pl, ne, pe);
			}
		}
		free(stack);
		free(pl);
		free(pe);
	}
	free(bvh_nodes);
	scene_del(scene);
	IMAGE img = image_init(conf.image_width, conf.image_height);
	image_fill(img, vec_init1(0.0));
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		for (int j = 0; j < img.w * img.h; ++j) {
			img.data[j] = vec_add(img.data[j], imgs[i].data[j]);
		}
	}
	for (int i = 0; i < img.w * img.h; ++i) {
		img.data[i] = vec_scale(1.0 / conf.num_threads, img.data[i]);
	}
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		image_del(imgs[i]);
	}
	image_write("bpt.pfm", img);
	t = time_get() - t;
	printf("time %f\n", t);
	return;
}
typedef struct conf_bpm CONF_BPM;

struct conf_bpm {
	int num_alloc_bvh_nodes;
	int num_threads;
	int num_samples;
	int image_width, image_height;
	double fx, fy;
	double lens_radius;
	double lens_flen;
	int depth;
	int nls;
	int *ils;
	double *cdf;
	MEDIUM m;
	double init_rad;
	double alpha;
	int num_verts;
	int num_nodes;
	int num_iter;
	int num_alloc_bvh_bpm_nodes;
};

typedef struct bpm_vertex BPM_VERTEX;

struct bpm_vertex {
	VEC tp;
	double u, v, paf, pab;
	int id;
};
typedef struct bpm_node BPM_NODE;

struct bpm_node {
	int m, n;
};

void conf_bpm_del(CONF_BPM conf) {
	free(conf.ils);
	free(conf.cdf);
}
CONF_BPM read_conf_bpm(char *path) {
	FILE *file = fopen(path, "r");
	if (file == NULL) {
		puts("file error");
		exit(0);
	}
	CONF_BPM conf;
	double fovh;
	double fdist;
	int nbuf = 256, nid = 256;
	char buf[nbuf], id[nid];
	while (true) {
		if (fgets(buf, nbuf, file) == NULL) {
			break;
		}
		sscanf(buf, "%s", id);
		if (strcmp(id, "num-alloc-bvh") == 0) {
			sscanf(buf, "%*s%d", &conf.num_alloc_bvh_nodes);
		} else if (strcmp(id, "num-threads") == 0) {
			sscanf(buf, "%*s%d", &conf.num_threads);
		} else if (strcmp(id, "image-dim") == 0) {
			sscanf(buf, "%*s%d%d", &conf.image_width, &conf.image_height);
		} else if (strcmp(id, "num-samples") == 0) {
			sscanf(buf, "%*s%d", &conf.num_samples);
		} else if (strcmp(id, "fovh") == 0) {
			sscanf(buf, "%*s%lf", &fovh);
		} else if (strcmp(id, "lens-radius") == 0) {
			sscanf(buf, "%*s%lf", &conf.lens_radius);
		} else if (strcmp(id, "lens-fdist") == 0) {
			sscanf(buf, "%*s%lf", &fdist);
		} else if (strcmp(id, "depth") == 0) {
			sscanf(buf, "%*s%d", &conf.depth);
		} else if (strcmp(id, "medium") == 0) {
			sscanf(buf, "%*s%lf%lf%lf%lf", &conf.m.n, &conf.m.a.x, &conf.m.a.y, &conf.m.a.z);
		} else if (strcmp(id, "init-rad") == 0) {
			sscanf(buf, "%*s%lf", &conf.init_rad);
		} else if (strcmp(id, "alpha") == 0) {
			sscanf(buf, "%*s%lf", &conf.alpha);
		} else if (strcmp(id, "num-verts") == 0) {
			sscanf(buf, "%*s%d", &conf.num_verts);
		} else if (strcmp(id, "num-nodes") == 0) {
			sscanf(buf, "%*s%d", &conf.num_nodes);
		} else if (strcmp(id, "num-iter") == 0) {
			sscanf(buf, "%*s%d", &conf.num_iter);
		} else if (strcmp(id, "num-alloc-bvh-bpm-nodes") == 0) {
			sscanf(buf, "%*s%d", &conf.num_alloc_bvh_bpm_nodes);
		}
	}
	fclose(file);
	film_dim(conf.image_width, conf.image_height, fovh, &conf.fx, &conf.fy);
	conf.lens_flen = focal_length(1.0, fdist);
	return conf;
}
int bpm_traverse(SCENE scene, int depth, double lensr, MEDIUM me, BVH_NODE *bvh_nodes, int *stack, PRNG_STATE *s, bool type, MEDIUM m, VEC x, VEC w, VEC tp, double ptmp, BPM_VERTEX *verts, int ptr, double cosgom1) {
	VEC rbm1;
	double psbm1;
	double cosgom2;
	double cossim1;
	double tm1;
	for (int i = 1; i <= depth; ++i) {
		VEC tuv = vec_init1(INFINITY);
		int id;
		if (!traverse(lensr, scene.verts, scene.faces, bvh_nodes, stack, x, w, &tuv, &id)) {
			return ptr + i;
		}
		x = vec_add(x, vec_scale(tuv.x, w));
		BSDF bsdf;
		VEC ng, ns;
		MEDIUM mp, mn;
		if (id == -1) {
			bsdf.t = 0;
			ng = disc_normal();
			ns = disc_normal();
			mp = me;
			mn = me;
		} else {
			FACE face = scene.faces[id];
			MAT mat = scene.mats[face.im];
			bsdf = scene.bsdfs[mat.ib];
			mp = scene.meds[mat.imp];
			mn = scene.meds[mat.imn];
			VEC v0, v1, v2;
			face_get_verts(scene.verts, face, &v0, &v1, &v2);
			ng = tri_normg(v0, v1, v2);
			VEC n0, n1, n2;
			face_get_norms(scene.norms, face, &n0, &n1, &n2);
			ns = tri_norms(n0, n1, n2, tuv.y, tuv.z);
		}
		double cosgi = vec_dot(ng, w), cossi = vec_dot(ns, w);
		if (!surf_kern1(cosgi, cossi)) {
			return ptr + i;
		}
		VEC tr = transmittance(m.a, tuv.x);
		tp = vec_mul(vec_scale(fabs(cossi / cosgi), tr), tp);
		
		if (i > 1) {
			BPM_VERTEX *vm2 = verts + ptr + i - 2;
			vm2->pab = rgb2lum(vec_mul(tr, rbm1)) * psbm1 * fabs(cossim1 * cosgom2) / sqr(tm1);
		}
		
		BPM_VERTEX v;
		v.tp = tp;
		v.u = tuv.y;
		v.v = tuv.z;
		v.paf = ptmp * fabs(cosgi) / sqr(tuv.x);
		v.id = id;
		verts[ptr + i] = v;
		
		MEDIUM mr, mt;
		medrt(cosgi, mp, mn, &mr, &mt);
		VEC tpd, rf;
		VEC wo = bsdf_sample_fb(bsdf, mr.n, mt.n, ns, w, s, cossi, type, &tpd, &rf, &rbm1, &ptmp, &psbm1);
		double cosgo = vec_dot(ng, wo), cosso = vec_dot(ns, wo);
		
		if (!surf_kern1(cosgo, cosso)) {
			return ptr + i + 1;
		}
		double pr = rgb2lum(vec_mul(tr, rf));
		if (prng_db(s) >= pr) {
			return ptr + i + 1;
		}
		
		tp = vec_scale(1.0 / pr, vec_mul(tp, tpd));
		ptmp *= fabs(cosso) * pr;
		w = wo;
		medt(cosgo, mp, mn, &m);
		cossim1 = cossi;
		tm1 = tuv.x;
		cosgom2 = cosgom1;
		cosgom1 = cosgo;
	}
	return ptr + depth + 1;
}
BPM_VERTEX bpm_sample_light(SCENE scene, int nls, int *ils, double *cdf, PRNG_STATE *s, VEC *tp, VEC *x, VEC *wo, double *ptmp, bool *kern, MEDIUM *m, double *cosgo) {
	BPM_VERTEX vert;
	vert.id = ils[discrete_sample(nls, cdf, s)];
	FACE face = scene.faces[vert.id];
	VEC v0, v1, v2;
	face_get_verts(scene.verts, face, &v0, &v1, &v2);
	VEC ng = tri_normg(v0, v1, v2);
	double u, v;
	tri_sample(s, &u, &v);
	vert.u = u;
	vert.v = v;
	*x = tri_point(v0, v1, v2, u, v);
	VEC n0, n1, n2;
	face_get_norms(scene.norms, face, &n0, &n1, &n2);
	VEC ns = tri_norms(n0, n1, n2, u, v);
	vert.paf = discrete_prob(cdf, search(nls, ils, vert.id)) * tri_pdf(v0, v1, v2);
	*tp = vec_init1(1.0 / vert.paf);
	vert.tp = *tp;
	MAT mat = scene.mats[face.im];
	EDF edf = scene.edfs[mat.ie];
	*wo = edf_sample_f(edf, ns, s, tp, ptmp);
	*cosgo = vec_dot(ng, *wo);
	double cosso = vec_dot(ns, *wo);
	*ptmp = fabs(cosso) * *ptmp;
	*kern = surf_kern1(*cosgo, cosso);
	medt(*cosgo, scene.meds[mat.imp], scene.meds[mat.imn], m);
	return vert;
}
int bpm_gen_lightpath(SCENE scene, int depth, double lensr, MEDIUM me, int nls, int *ils, double *cdf, BVH_NODE *bvh_nodes, int *bvh_stack, PRNG_STATE *s, BPM_VERTEX *verts, int ptr) {
	VEC tp, x, w;
	double ptmp;
	bool kern;
	MEDIUM m;
	double cosgo;
	verts[ptr] = bpm_sample_light(scene, nls, ils, cdf, s, &tp, &x, &w, &ptmp, &kern, &m, &cosgo);
	if (!kern) {
		return ptr + 1;
	}
	return bpm_traverse(scene, depth, lensr, me, bvh_nodes, bvh_stack, s, false, m, x, w, tp, ptmp, verts, ptr, cosgo);
}
typedef struct bvh_bpm_build_entry BVH_BPM_BUILD_ENTRY;
struct bvh_bpm_build_entry {
	AABB b;
	int i;
};

int bvh_sphere_cmp_x(void const *_a, void const *_b) {
	BVH_BPM_BUILD_ENTRY const *a = _a;
	BVH_BPM_BUILD_ENTRY const *b = _b;
	double ca = 0.5 * (a->b.min.x + a->b.max.x);
	double cb = 0.5 * (b->b.min.x + b->b.max.x);
	return ca < cb ? -1 : ca > cb ? 1 : 0;
}
int bvh_sphere_cmp_y(void const *_a, void const *_b) {
	BVH_BPM_BUILD_ENTRY const *a = _a;
	BVH_BPM_BUILD_ENTRY const *b = _b;
	double ca = 0.5 * (a->b.min.y + a->b.max.y);
	double cb = 0.5 * (b->b.min.y + b->b.max.y);
	return ca < cb ? -1 : ca > cb ? 1 : 0;
}
int bvh_sphere_cmp_z(void const *_a, void const *_b) {
	BVH_BPM_BUILD_ENTRY const *a = _a;
	BVH_BPM_BUILD_ENTRY const *b = _b;
	double ca = 0.5 * (a->b.min.z + a->b.max.z);
	double cb = 0.5 * (b->b.min.z + b->b.max.z);
	return ca < cb ? -1 : ca > cb ? 1 : 0;
}
AABB sphere_aabb(double r, VEC c) {
	AABB box = {vec_init(c.x - r, c.y - r, c.z - r), vec_init(c.x + r, c.y + r, c.z + r)};
	return box;
}
double aabb_volume(AABB box) {
	VEC d = vec_sub(box.max, box.min);
	return d.x * d.y * d.z;
}
void bvh_bpm_build(BVH_BPM_BUILD_ENTRY *bes, double *al, double *ar, BVH_NODE *bvh_nodes, BVH_BUILD_QUERY *stack, BVH_BUILD_QUERY q, int ptr_deq, int *ptr_enq) {
	AABB box = aabb_default();
	for (int i = q.s; i < q.e; ++i) {
		box = aabb_merge(box, bes[i].b);
	}
	double voli = 1.0 / aabb_volume(box);
	double cost = q.e - q.s;
	int axis = -1, piv = -1;
	for (int i = 0; i < 3; ++i) {
		switch (i) {
			case 0 : {
				qsort(bes + q.s, q.e - q.s, sizeof(BVH_BPM_BUILD_ENTRY), bvh_sphere_cmp_x);
			} break;
			case 1 : {
				qsort(bes + q.s, q.e - q.s, sizeof(BVH_BPM_BUILD_ENTRY), bvh_sphere_cmp_y);
			} break;
			case 2 : {
				qsort(bes + q.s, q.e - q.s, sizeof(BVH_BPM_BUILD_ENTRY), bvh_sphere_cmp_z);
			} break;
		}
		AABB boxc;
		boxc = aabb_default();
		al[0] = 0.0;
		for (int j = q.s; j < q.e; ++j) {
			boxc = aabb_merge(boxc, bes[j].b);
			al[j - q.s + 1] = aabb_volume(boxc);
		}
		boxc = aabb_default();
		ar[q.e - q.s] = 0.0;
		for (int j = q.e - 1; j >= q.s; --j) {
			boxc = aabb_merge(boxc, bes[j].b);
			ar[j - q.s] = aabb_volume(boxc);
		}
		for (int j = q.s; j <= q.e; ++j) {
			double cost1 = 2.0 + voli * (al[j - q.s] * (j - q.s) + ar[j - q.s] * (q.e - j));
			if (cost1 < cost) {
				cost = cost1;
				axis = i;
				piv = j;
			}
		}
	}
	if (axis == -1) {
		BVH_NODE node = {box, q.s - q.e, q.s};
		bvh_nodes[ptr_deq] = node;
		return;
	}
	switch (axis) {
		case 0 : {
			qsort(bes + q.s, q.e - q.s, sizeof(BVH_BPM_BUILD_ENTRY), bvh_sphere_cmp_x);
		} break;
		case 1 : {
			qsort(bes + q.s, q.e - q.s, sizeof(BVH_BPM_BUILD_ENTRY), bvh_sphere_cmp_y);
		} break;
		case 2 : {
			qsort(bes + q.s, q.e - q.s, sizeof(BVH_BPM_BUILD_ENTRY), bvh_sphere_cmp_z);
		} break;
	}
	BVH_NODE node = {box, *ptr_enq + 1, *ptr_enq + 2};
	bvh_nodes[ptr_deq] = node;
	BVH_BUILD_QUERY ql = {q.s, piv};
	PUSH(stack, *ptr_enq, ql);
	BVH_BUILD_QUERY qr = {piv, q.e};
	PUSH(stack, *ptr_enq, qr);
	return;
}
double vcm_radius(double a, double ri, int i) {
	return pow(i + 1, 0.5 * (a - 1.0)) * ri;
}
int bvh_bpm_build_wrap(VEC *verts, FACE *faces, double r, int num_samples, int num_bpm_nodes, BPM_NODE **bpm_nodes, BPM_NODE **bpm_nodes_cp, int num_alloc_bvh_bpm_nodes, BVH_NODE *bvh_bpm_nodes, BVH_BPM_BUILD_ENTRY *bes, BVH_BUILD_QUERY *bqs, double *al, double *ar, BPM_VERTEX **vss, int **ess) {
	printf("proc bvh-bpm-build\n");
	double t = time_get();
	for (int i = 0; i < num_bpm_nodes; ++i) {
		BPM_NODE bpm_node = (*bpm_nodes)[i];
		int iz = bpm_node.m / num_samples, iy = bpm_node.m % num_samples, ix = bpm_node.n;
		BPM_VERTEX v = vss[iz][iy == 0 ? 0 : ess[iz][iy - 1] + ix];
		VEC v0, v1, v2;
		face_get_verts(verts, faces[v.id], &v0, &v1, &v2);
		BVH_BPM_BUILD_ENTRY e = {sphere_aabb(r, tri_point(v0, v1, v2, v.u, v.v)), i};
		bes[i] = e;
	}
	int ptr_deq = 0, ptr_enq = 0;
	BVH_BUILD_QUERY q = {0, num_bpm_nodes};
	while (true) {
		bvh_bpm_build(bes, al, ar, bvh_bpm_nodes, bqs, q, ptr_deq, &ptr_enq);
		if (ptr_deq == ptr_enq) {
			break;
		}
		q = DEQ(bqs, ptr_deq);
	}
	for (int i = 0; i < num_bpm_nodes; ++i) {
		(*bpm_nodes_cp)[i] = (*bpm_nodes)[bes[i].i];
	}
	BPM_NODE *cp = *bpm_nodes;
	*bpm_nodes = *bpm_nodes_cp;
	*bpm_nodes_cp = cp;
	t = time_get() - t;
	printf("time %f\n", t);
	printf("usage %d/%d\n", ptr_deq + 1, num_alloc_bvh_bpm_nodes);
	return ptr_deq + 1;
}
void render_vcm(char *path_scene, char *path_conf) {
	printf("proc render-vcm\n");
	double t = time_get();
	SCENE scene = read_scene(path_scene);
	CONF_BPM conf = read_conf_bpm(path_conf);
	BVH_NODE *bvh_nodes = ALLOC(BVH_NODE, conf.num_alloc_bvh_nodes);
	int num_bvh_nodes = bvh_build_wrap(&scene, conf.num_alloc_bvh_nodes, bvh_nodes);
	conf.nls = init_cdf(scene, &conf.ils, &conf.cdf);
	omp_set_num_threads(conf.num_threads);
	IMAGE imgs[omp_get_max_threads()];
	int *bvh_stacks[omp_get_max_threads()];
	BPM_VERTEX *vss[omp_get_max_threads()];
	BPT_VERTEX *pes[omp_get_max_threads()];
	int *ess[omp_get_max_threads()];
	PRNG_STATE ss[omp_get_max_threads()];
	BPM_NODE *nodes = ALLOC(BPM_NODE, conf.num_nodes);
	BPM_NODE *nodes_cp = ALLOC(BPM_NODE, conf.num_nodes);
	BVH_NODE *bvh_bpm_nodes = ALLOC(BVH_NODE, conf.num_alloc_bvh_bpm_nodes);
	double *al = ALLOC(double, conf.num_nodes + 1);
	double *ar = ALLOC(double, conf.num_nodes + 1);
	BVH_BPM_BUILD_ENTRY *bes = ALLOC(BVH_BPM_BUILD_ENTRY, conf.num_nodes);
	BVH_BUILD_QUERY *bqs = ALLOC(BVH_BUILD_QUERY, conf.num_nodes);
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		imgs[i] = image_init(conf.image_width, conf.image_height);
		image_fill(imgs[i], vec_init1(0.0));
		bvh_stacks[i] = ALLOC(int, num_bvh_nodes);
		vss[i] = ALLOC(BPM_VERTEX, conf.num_verts);
		pes[i] = ALLOC(BPT_VERTEX, conf.depth + 1);
		ess[i] = ALLOC(int, conf.num_samples);
		ss[i] = prng_seed(SEED);
		for (int j = 0; j < i; ++j) {
			ss[i] = prng_jump(ss[i]);
		}
	}
	for (int i = 0; i < conf.num_iter; ++i) {
		printf("iteration %d\n", i);
		#pragma omp parallel
		{
			int *bvh_stack = bvh_stacks[omp_get_thread_num()];
			BPM_VERTEX *vs = vss[omp_get_thread_num()];
			int *es = ess[omp_get_thread_num()];
			PRNG_STATE *s = ss + omp_get_thread_num();
			int ptr = 0;
			for (int j = 0; j < conf.num_samples; ++j) {
				ptr = bpm_gen_lightpath(scene, conf.depth, conf.lens_radius, conf.m, conf.nls, conf.ils, conf.cdf, bvh_nodes, bvh_stack, s, vs, ptr);
				es[j] = ptr;
			}
		}
		printf("vertex usage\n");
		for (int j = 0; j < omp_get_max_threads(); ++j) {
			printf("thread %d %d/%d\n", j, ess[j][conf.num_samples - 1], conf.num_verts);
		}
		int num_nodes = 0;
		for (int j = 0; j < conf.num_samples * omp_get_max_threads(); ++j) {
			int iz = j / conf.num_samples, iy = j % conf.num_samples;
			int s = iy == 0 ? 0 : ess[iz][iy - 1], e = ess[iz][iy];
			for (int k = 1; k < e - s; ++k) {
				BPM_VERTEX v = vss[iz][s + k];
				if (v.id == -1) {
					continue;
				}
				if (!bsdf_con(scene.bsdfs[scene.mats[scene.faces[v.id].im].ib])) {
					continue;
				}
				BPM_NODE node = {j, k};
				PUSH(nodes, num_nodes, node);
			}
		}
		printf("node usage %d/%d\n", num_nodes, conf.num_nodes);
		double r = vcm_radius(conf.alpha, conf.init_rad, i);
		int num_bvh_bpm_nodes = bvh_bpm_build_wrap(scene.verts, scene.faces, r, conf.num_samples, num_nodes, &nodes, &nodes_cp, conf.num_alloc_bvh_bpm_nodes, bvh_bpm_nodes, bes, bqs, al, ar, vss, ess);
		#pragma omp parallel
		{
			int *bvh_stack = bvh_stacks[omp_get_thread_num()];
			PRNG_STATE *s = ss + omp_get_thread_num();
			BPT_VERTEX *pe = pes[omp_get_thread_num()];
			IMAGE img = imgs[omp_get_thread_num()];
			for (int j = 0; j < conf.num_samples; ++j) {
				double x = img.w * prng_db(s), y = img.h * prng_db(s);
				int ne = bpt_gen_eyepath(scene, conf.depth, conf.lens_radius, conf.lens_flen, conf.image_width, conf.image_height, conf.fx, conf.fy, conf.m, bvh_nodes, bvh_stack, s, x, y, pe);
				
			}
		}
	}
	IMAGE img = image_init(conf.image_width, conf.image_height);
	image_fill(img, vec_init1(0.0));
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		for (int j = 0; j < img.w * img.h; ++j) {
			img.data[j] = vec_add(img.data[j], vec_scale(1.0 / omp_get_max_threads(), imgs[i].data[j]));
		}
	}
	image_write("vcm.pfm", img);
	scene_del(scene);
	conf_bpm_del(conf);
	free(bvh_nodes);
	free(nodes);
	free(nodes_cp);
	free(bvh_bpm_nodes);
	free(al);
	free(ar);
	free(bes);
	free(bqs);
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		image_del(imgs[i]);
		free(bvh_stacks[i]);
		free(vss[i]);
		free(pes[i]);
		free(ess[i]);
	}
	t = time_get() - t;
	printf("time %f\n", t);
	return;
}
void dispatch(int argc, char **argv) {
	if (argc < 4) {
		puts("insufficient arguments");
		exit(0);
	}
	char *method = argv[1];
	char *path_scene = argv[2];
	char *path_conf = argv[3];
	if (strcmp(method, "debug") == 0) {
		render_debug(path_scene, path_conf);
	} else if (strcmp(method, "ibl") == 0) {
		render_ibl(path_scene, path_conf);
	} else if (strcmp(method, "pt") == 0) {
		render_pt(path_scene, path_conf);
	} else if (strcmp(method, "bpt") == 0) {
		render_bpt(path_scene, path_conf);
	} else if (strcmp(method, "vcm") == 0) {
		render_vcm(path_scene, path_conf);
	}
	return;
}
double D(double r, double cos) {
	return 1.0 / (PI * sqr(r * sqr(cos) + fabs(1.0 - sqr(cos)) / r));
}
VEC spl(VEC n, PRNG_STATE *s) {
	double phi = 2.0 * PI * prng_db(s);
	double cost = 2.0 * prng_db(s) - 1.0, sint = opst(cost), cosp = cos(phi), sinp = sin(phi);
	VEC wo = pol2cartz1(cost, sint, cosp, sinp);
	return matrix_trafo(matrix_tan(n), wo);
}
VEC eval_cond(BSDF bsdf, double ni, double no, VEC n, VEC wi, VEC wo) {
	double cosni = vec_dot(n, wi), cosno = vec_dot(n, wo);
	if (cosni * cosno > 0.0) {
		return vec_init1(0.0);
	}
	VEC m = vec_norm(vec_sub(wo, wi));
	double cosnm = vec_dot(n, m);
	double cosmi = vec_dot(m, wi);
	VEC fr = vec_init(fres_cond(bsdf.n.x / ni, bsdf.k.x, cosmi), fres_cond(bsdf.n.y / ni, bsdf.k.y, cosmi), fres_cond(bsdf.n.z / ni, bsdf.k.z, cosmi));
	double g = ggx_g2(bsdf.r, cosni, cosno);
	double d = D(bsdf.r, fabs(cosnm));
	return vec_scale(g * d / fabs(4.0 * cosni * cosno), fr);
}
VEC eval_diel(BSDF bsdf, double nr, double nt, VEC n, VEC wi, VEC wo) {
	double cosni = vec_dot(n, wi), cosno = vec_dot(n, wo);
	if (cosni * cosno > 0.0) {
		VEC m = vec_norm(vec_sub(vec_scale(nt, wo), vec_scale(nr, wi)));
		m = vec_scale(sign(vec_dot(n, m)), m);
		double cosnm = vec_dot(n, m);
		double cosmi = vec_dot(m, wi), cosmo = vec_dot(m, wo);
		if (cosni * cosmi < 0.0 || cosno * cosmo < 0.0) {
			return vec_init1(0.0);
		}
		double eta = nr / nt;
		double eta2 = sqr(eta);
		double det = 1.0 - eta2 * (1.0 - sqr(cosmi));
		double ft = 0.0;
		if (det > 0.0) {
			double cosmt = sign(cosmi) * sqrt(det);
			ft = 1.0 - fres_diel(eta, cosmi, cosmt);
		}
		double g = ggx_g2(bsdf.r, cosni, cosno);
		double d = D(bsdf.r, fabs(cosnm));
		return vec_init1(g * d * ft * fabs(cosmi * cosmo / (cosni * cosno)) * sqr(nr) / sqr(nt * cosmo - nr * cosmi) / eta2);
	} else {
		VEC m = vec_norm(vec_sub(wo, wi));
		double cosnm = vec_dot(n, m);
		double cosmi = vec_dot(m, wi);
		double eta = nr / nt;
		double eta2 = sqr(eta);
		double det = 1.0 - eta2 * (1.0 - sqr(cosmi));
		double fr = 1.0;
		if (det > 0.0) {
			double cosmt = sign(cosmi) * sqrt(det);
			fr = fres_diel(eta, cosmi, cosmt);
		}
		double g = ggx_g2(bsdf.r, cosni, cosno);
		double d = D(bsdf.r, fabs(cosnm));
		return vec_init1(g * d * fr / fabs(4.0 * cosni * cosno));
	}
}
void test_eval(void) {
	int N = 100000000;
	MEDIUM mp = {1.3, {0.0, 0.0, 0.0}};
	MEDIUM mn = {1.5, {0.0, 0.0, 0.0}};
	BSDF bsdf = {2, {0.9, 0.9, 0.9}, 0.1, {0.1431189557, 0.3749570432, 1.4424785571}, {3.9831604247, 2.3857207478, 1.603215289}};
	VEC ns = vec_init(0.0, 0.0, 1.0);
	VEC wi = vec_scale(-1.0, pol2cartz0(PI / 4.0, PI / 4.0));
	IMAGE img = image_init(512, 256);
	PRNG_STATE s = prng_seed(SEED);
	image_fill(img, vec_init1(0.0));
	if (true)
	for (int i = 0; i < N; ++i) {
		double cosni = vec_dot(ns, wi);
		MEDIUM mr, mt;
		medrt(cosni, mp, mn, &mr, &mt);
		VEC tp = vec_init1(1.0), rf = vec_init1(1.0), rb = vec_init1(1.0);
		double psf, psb;
		VEC wo = bsdf_sample_fb(bsdf, mr.n, mt.n, ns, wi, &s, cosni, false, &tp, &rf, &rb, &psf, &psb);
		double pr = rgb2lum(rf);
		if (prng_db(&s) >= pr) {
			continue;
		}
		tp = vec_scale(1.0 / pr, tp);
		int ix = iclamp(0, img.w - 1, img.w * (atan2(wo.y, wo.x) + (wo.y < 0.0 ? 2.0 * PI : 0.0)) / (2.0 * PI));
		int iy = iclamp(0, img.h - 1, img.h * acos(clamp(-1.0, 1.0, wo.z)) / PI);
		img.data[iy * img.w + ix] = vec_add(img.data[iy * img.w + ix], vec_scale(img.w * img.h / CAST(double, N), tp));
	}
	image_write("spl.pfm", img);
	image_fill(img, vec_init1(0.0));
	if (true)
	for (int i = 0; i < N; ++i) {
		VEC wo = spl(ns, &s);
		double cosni = -vec_dot(ns, wo), cosno = -vec_dot(ns, wi);
		MEDIUM mr, mt;
		medrt(cosni, mp, mn, &mr, &mt);
		VEC tp = vec_init1(1.0), rf = vec_init1(1.0), rb = vec_init1(1.0);
		double psf, psb;
		if (!bsdf_eval_fb(bsdf, mr.n, mt.n, ns, vec_scale(-1.0, wo), vec_scale(-1.0, wi), cosni, cosno, &tp, &rf, &rb, &psf, &psb)) {
			continue;
		}
		tp = vec_scale(4.0 * PI * fabs(cosni), tp);
		int ix = iclamp(0, img.w - 1, img.w * (atan2(wo.y, wo.x) + (wo.y < 0.0 ? 2.0 * PI : 0.0)) / (2.0 * PI));
		int iy = iclamp(0, img.h - 1, img.h * acos(clamp(-1.0, 1.0, wo.z)) / PI);
		img.data[iy * img.w + ix] = vec_add(img.data[iy * img.w + ix], vec_scale(img.w * img.h / CAST(double, N), tp));
	}
	image_write("eval.pfm", img);
}
int main(int argc, char **argv) {
	dispatch(argc, argv);
	//test_eval();
	return 0;
}