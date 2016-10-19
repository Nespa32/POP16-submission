

#ifndef _HIMENO_CONTEST_
#define _HIMENO_CONTEST_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>


extern const double THRESHOLD_MS;
extern const long unsigned int MAXIMUM_ITERARIONS;
extern const double REFERENCE_MS[];

#define MR(mt, n, r, c, d)  mt->m[(n) * mt->mrows * mt->mcols * mt->mdeps + (r) * mt->mcols* mt->mdeps + (c) * mt->mdeps + (d)]

struct Mat {
    float *m;
    int mrows;
    int mcols;
    int mdeps;
};

/* prototypes */
typedef struct Mat Matrix;

float jacobi(Matrix *a, Matrix *b, Matrix *c, Matrix *p, Matrix *bnd, Matrix *wrk1, Matrix *wrk2);

int newMat(Matrix *Mat, int mnums, int mrows, int mcols, int mdeps);

void clearMat(Matrix *Mat);

void set_param(int i[], char *size);

void mat_set(Matrix *Mat, int l, float z);

void mat_set_init(Matrix *Mat);

void get_matrix_size(int argc, char **argv, int *mimax, int *mjmax, int *mkmax);

void clear_matrices(Matrix *p, Matrix *bnd, Matrix *wrk1, Matrix *wrk2, Matrix *a, Matrix *b, Matrix *c);

void init_matrices(int mimax, int mjmax, int mkmax, Matrix *a, Matrix *b, Matrix *c, Matrix *p, Matrix *bnd, Matrix *wrk1, Matrix *wrk2);

void verify_input_matrices(Matrix *a, Matrix *b, Matrix *c, Matrix *p, Matrix *bnd, Matrix *wrk1, Matrix *wrk2);

void mat_verify(Matrix *Mat, int l, float val);

void mat_verify_init(Matrix *Mat);

bool verify_outputs(float gosa, Matrix *matrix);

void write_outputs(float gosa, Matrix *wrk2);

double get_reference_time();

#endif
