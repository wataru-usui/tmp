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
typedef struct face face;
typedef struct med med;
typedef struct sdf sdf;
typedef struct scene scene;
f8 timer(void) {
	return clock() / (f8)(CLOCKS_PER_SEC);
}
struct i4s2 {
	i4s x, y;
};
struct i4s3 {
	i4s x, y, z;
};
struct f83 {
	f8 x, y, z;
};
struct face {
	i4s ei, bi;
	i4s2 mis;
	i4s3 vis, nis;
};
struct med {
	f8 n;
	f83 a;
};
struct sdf {
	i4s t;
	f8 ps[8];
};
struct scene {
	i4s esn, bsn, msn, vsn, nsn, fsn;
	sdf *es, *bs;
	med *ms;
	f83 *vs, *ns;
	face *fs;
};
void scene_free(scene s) {
	free(s.es);
	free(s.bs);
	free(s.ms);
	free(s.vs);
	free(s.ns);
	free(s.fs);
	return;
}
void scene_print(scene s) {
	printf("%d\n", s.esn);
	for (i4s i = 0; i < s.esn; ++i) {
		sdf a = s.es[i];
		printf("%d %f %f %f %f %f %f %f %f\n", a.t, a.ps[0], a.ps[1], a.ps[2], a.ps[3], a.ps[4], a.ps[5], a.ps[6], a.ps[7]);
	}
	printf("%d\n", s.bsn);
	for (i4s i = 0; i < s.bsn; ++i) {
		sdf a = s.bs[i];
		printf("%d %f %f %f %f %f %f %f %f\n", a.t, a.ps[0], a.ps[1], a.ps[2], a.ps[3], a.ps[4], a.ps[5], a.ps[6], a.ps[7]);
	}
	printf("%d\n", s.msn);
	for (i4s i = 0; i < s.msn; ++i) {
		med a = s.ms[i];
		printf("%f %f %f %f\n", a.n, a.a.x, a.a.y, a.a.z);
	}
	printf("%d\n", s.vsn);
	for (i4s i = 0; i < s.vsn; ++i) {
		f83 a = s.vs[i];
		printf("%f %f %f\n", a.x, a.y, a.z);
	}
	printf("%d\n", s.nsn);
	for (i4s i = 0; i < s.nsn; ++i) {
		f83 a = s.ns[i];
		printf("%f %f %f\n", a.x, a.y, a.z);
	}
	printf("%d\n", s.fsn);
	for (i4s i = 0; i < s.fsn; ++i) {
		face a = s.fs[i];
		printf("%d %d %d %d %d %d %d %d %d %d \n", a.ei, a.bi, a.mis.x, a.mis.y, a.vis.x, a.vis.y, a.vis.z, a.nis.x, a.nis.y, a.nis.z);
	}
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
scene scene_read_obj(i1 *p) {
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
	s.esn = 1;
	s.bsn = 1;
	s.msn = 2;
	s.vsn = 0;
	s.nsn = 0;
	s.fsn = 0;
	while (fgets(l, ln, f)) {
		sscanf(l, "%s", k);
		if (strcmp(k, "v") == 0) {
			++s.vsn;
			continue;
		}
		if (strcmp(k, "vn") == 0) {
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
	sdf e = {0, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
	s.es[0] = e;
	sdf b = {0, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
	s.bs[0] = b;
	med m0 = {1.0, {0.0, 0.0, 0.0}};
	s.ms[0] = m0;
	med m1 = {1.0, {0.0, 0.0, 0.0}};
	s.ms[1] = m1;
	s.vsn = 0;
	s.nsn = 0;
	s.fsn = 0;
	rewind(f);
	while (fgets(l, ln, f)) {
		sscanf(l, "%s", k);
		if (strcmp(k, "v") == 0) {
			f83 a;
			sscanf(l, "%*s%lf%lf%lf", &a.x, &a.y, &a.z);
			s.vs[s.vsn++] = a;
			continue;
		}
		if (strcmp(k, "vn") == 0) {
			f83 a;
			sscanf(l, "%*s%lf%lf%lf", &a.x, &a.y, &a.z);
			s.ns[s.nsn++] = a;
			continue;
		}
		if (strcmp(k, "f") == 0) {
			face f;
			f.ei = 0;
			f.bi = 0;
			i4s2 mis = {0, 1};
			f.mis = mis;
			for (i4s i = 0; l[i] != '\0'; ++i) {
				if (l[i] == '/') {
					l[i] = ' ';
				}
			}
			sscanf(l, "%*s%d%d%d%d%d%d", &f.vis.x, &f.nis.x, &f.vis.y, &f.nis.y, &f.vis.z, &f.nis.z);
			--f.vis.x;
			--f.vis.y;
			--f.vis.z;
			--f.nis.x;
			--f.nis.y;
			--f.nis.z;
			s.fs[s.fsn++] = f;
			continue;
		}
	}
	fclose(f);
	return s;
}
scene scene_merge(i4s psn, i1 **ps) {
	scene ss[psn];
	for (i4s i = 0; i < psn; ++i) {
		ss[i] = scene_read(ps[i]);
	}
	scene s;
	s.esn = 0;
	s.bsn = 0;
	s.msn = 0;
	s.vsn = 0;
	s.nsn = 0;
	s.fsn = 0;
	for (i4s i = 0; i < psn; ++i) {
		scene a = ss[i];
		s.esn += a.esn;
		s.bsn += a.bsn;
		s.msn += a.msn;
		s.vsn += a.vsn;
		s.nsn += a.nsn;
		s.fsn += a.fsn;
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
	for (i4s i = 0; i < psn; ++i) {
		scene a = ss[i];
		for (i4s j = 0; j < a.esn; ++j) {
			s.es[s.esn + j] = a.es[j];
		}
		for (i4s j = 0; j < a.bsn; ++j) {
			s.bs[s.bsn + j] = a.bs[j];
		}
		for (i4s j = 0; j < a.msn; ++j) {
			s.ms[s.msn + j] = a.ms[j];
		}
		for (i4s j = 0; j < a.vsn; ++j) {
			s.vs[s.vsn + j] = a.vs[j];
		}
		for (i4s j = 0; j < a.nsn; ++j) {
			s.ns[s.nsn + j] = a.ns[j];
		}
		for (i4s j = 0; j < a.fsn; ++j) {
			face f = a.fs[j];
			f.ei += s.esn;
			f.bi += s.bsn;
			f.mis.x += s.msn;
			f.mis.y += s.msn;
			f.vis.x += s.vsn;
			f.vis.y += s.vsn;
			f.vis.z += s.vsn;
			f.nis.x += s.nsn;
			f.nis.y += s.nsn;
			f.nis.z += s.nsn;
			s.fs[s.fsn + j] = f;
		}
		s.esn += a.esn;
		s.bsn += a.bsn;
		s.msn += a.msn;
		s.vsn += a.vsn;
		s.nsn += a.nsn;
		s.fsn += a.fsn;
	}
	for (i4s i = 0; i < psn; ++i) {
		scene_free(ss[i]);
	};
	return s;
}
void scene_write(scene s, i1 *p) {
	FILE *f = fopen(p, "wb");
	if (f == NULL) {
		puts("i/o err");
		exit(1);
	}
	for (i4s i = 0; i < s.esn; ++i) {
		sdf a = s.es[i];
		fprintf(f, "e %d %.*f %.*f %.*f %.*f %.*f %.*f %.*f %.*f \n", a.t, DECIMAL_DIG, a.ps[0], DECIMAL_DIG, a.ps[1], DECIMAL_DIG, a.ps[2], DECIMAL_DIG, a.ps[3], DECIMAL_DIG, a.ps[4], DECIMAL_DIG, a.ps[5], DECIMAL_DIG, a.ps[6], DECIMAL_DIG, a.ps[7]);
	}
	for (i4s i = 0; i < s.bsn; ++i) {
		sdf a = s.bs[i];
		fprintf(f, "b %d %.*f %.*f %.*f %.*f %.*f %.*f %.*f %.*f \n", a.t, DECIMAL_DIG, a.ps[0], DECIMAL_DIG, a.ps[1], DECIMAL_DIG, a.ps[2], DECIMAL_DIG, a.ps[3], DECIMAL_DIG, a.ps[4], DECIMAL_DIG, a.ps[5], DECIMAL_DIG, a.ps[6], DECIMAL_DIG, a.ps[7]);
	}
	for (i4s i = 0; i < s.msn; ++i) {
		med a = s.ms[i];
		fprintf(f, "m %.*f %.*f %.*f %.*f \n", DECIMAL_DIG, a.n, DECIMAL_DIG, a.a.x, DECIMAL_DIG, a.a.y, DECIMAL_DIG, a.a.z);
	}
	for (i4s i = 0; i < s.vsn; ++i) {
		f83 a = s.vs[i];
		fprintf(f, "v %.*f %.*f %.*f \n", DECIMAL_DIG, a.x, DECIMAL_DIG, a.y, DECIMAL_DIG, a.z);
	}
	for (i4s i = 0; i < s.nsn; ++i) {
		f83 a = s.ns[i];
		fprintf(f, "n %.*f %.*f %.*f \n", DECIMAL_DIG, a.x, DECIMAL_DIG, a.y, DECIMAL_DIG, a.z);
	}
	for (i4s i = 0; i < s.fsn; ++i) {
		face a = s.fs[i];
		fprintf(f, "f %d %d %d %d %d %d %d %d %d %d \n", a.ei, a.bi, a.mis.x, a.mis.y, a.vis.x, a.vis.y, a.vis.z, a.nis.x, a.nis.y, a.nis.z);
	}
	fclose(f);
	return;
}
i4s main(i4s argc, i1 **args) {
	//scene s = scene_read_obj("./light/1.obj");
	//scene_write(s, "./light/1.txt");
	i4s nps = 3;
	i1 *ps[] = {"./light/1.txt", "./wall/wall.txt", "./bunny/bunny.txt"};
	scene s = scene_merge(nps, ps);
	scene_write(s, "scene.txt");
	scene_free(s);
	return 0;
}