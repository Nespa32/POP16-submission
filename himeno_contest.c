/********************************************************************

 This benchmark test program is measuring a cpu performance
 of floating point operation by a Poisson equation solver.

 If you have any question, please ask me via email.
 written by Ryutaro HIMENO, November 26, 2001.
 Version 3.0
 ----------------------------------------------
 Ryutaro Himeno, Dr. of Eng.
 Head of Computer Information Division,
 RIKEN (The Institute of Pysical and Chemical Research)
 Email : himeno@postman.riken.go.jp
 ---------------------------------------------------------------
 You can adjust the size of this benchmark code to fit your target
 computer. In that case, please chose following sets of
 [mimax][mjmax][mkmax]:
 small : 33,33,65
 small : 65,65,129
 midium: 129,129,257
 large : 257,257,513
 ext.large: 513,513,1025
 This program is to measure a computer performance in MFLOPS
 by using a kernel which appears in a linear solver of pressure
 Poisson eq. which appears in an incompressible Navier-Stokes solver.
 A point-Jacobi method is employed in this solver as this method can
 be easyly vectrized and be parallelized.
 ------------------
 Finite-difference method, curvilinear coodinate system
 Vectorizable and parallelizable on each grid point
 No. of grid points : imax x jmax x kmax including boundaries
 ------------------
 A,B,C:coefficient matrix, wrk1: source term of Poisson equation
 wrk2 : working area, OMEGA : relaxation parameter
 BND:control variable for boundaries and objects ( = 0 or 1)
 P: pressure
********************************************************************/


#include <stdbool.h>
#include "himeno_contest.h"

const double THRESHOLD_MS = 500.0;
const long unsigned int MAXIMUM_ITERARIONS = (long unsigned int) 1e7;

static const int STR_SIZE = 10;
static const int FILENAME_SIZE = 50;
static const double EPSILON = 0.001;

const double REFERENCE_MS[] = {7.9, 50.86, 462.02, 3816.01, 57708.81};
static int matrix_dims[3];
static size_t matrix_size;


void get_matrix_size(int argc, char **argv, int *mimax, int *mjmax, int *mkmax) {

    char size[STR_SIZE];
    int msize[3];

    if (argc == 2) {
        strncpy(size, argv[1], STR_SIZE);
    } else {

        printf("Possible Sizes:\n");
        printf("\t    XS (32x32x64)\n");
        printf("\t    S  (64x64x128)\n");
        printf("\t    M  (128x128x256)\n");
        printf("\t    L  (256x256x512)\n");
        printf("\t    XL (512x512x1024)\n\n");
        printf("Grid Size = ");
        int res = scanf("%s", size);
        printf("\n");

        if (res == EOF) {
            printf("Could not read matrix size !!\n");
            fflush(stdout);
            exit(1);
        }
    }

    set_param(msize, size);

    *mimax = msize[0];
    *mjmax = msize[1];
    *mkmax = msize[2];

    matrix_dims[0] = msize[0];
    matrix_dims[1] = msize[1];
    matrix_dims[2] = msize[2];

    matrix_size = (size_t) msize[0] * msize[1] * msize[2];
}

void set_param(int is[], char *size) {

    if (!strcmp(size, "XS") || !strcmp(size, "xs")) {
        is[0] = 32;
        is[1] = 32;
        is[2] = 64;

        return;
    }
    if (!strcmp(size, "S") || !strcmp(size, "s")) {
        is[0] = 64;
        is[1] = 64;
        is[2] = 128;


        return;
    }
    if (!strcmp(size, "M") || !strcmp(size, "m")) {
        is[0] = 128;
        is[1] = 128;
        is[2] = 256;


        return;
    }
    if (!strcmp(size, "L") || !strcmp(size, "l")) {
        is[0] = 256;
        is[1] = 256;
        is[2] = 512;


        return;
    }
    if (!strcmp(size, "XL") || !strcmp(size, "xl")) {
        is[0] = 512;
        is[1] = 512;
        is[2] = 1024;

        return;
    } else {
        printf("Invalid input character !!\n");
        fflush(stdout);
        exit(1);
    }
}

int newMat(Matrix *Mat, int mnums, int mrows, int mcols, int mdeps) {

    Mat->mrows = mrows;
    Mat->mcols = mcols;
    Mat->mdeps = mdeps;
    Mat->m = NULL;
    Mat->m = (float *)
            malloc(mnums * mrows * mcols * mdeps * sizeof(float));

    return (Mat->m != NULL) ? 1 : 0;
}

void clear_matrices(Matrix *p, Matrix *bnd, Matrix *wrk1, Matrix *wrk2, Matrix *a, Matrix *b, Matrix *c) {

    clearMat(p);
    clearMat(bnd);
    clearMat(wrk1);
    clearMat(wrk2);
    clearMat(a);
    clearMat(b);
    clearMat(c);
}

void clearMat(Matrix *Mat) {
    if (Mat->m)
        free(Mat->m);
    Mat->m = NULL;
    Mat->mcols = 0;
    Mat->mrows = 0;
    Mat->mdeps = 0;
}

void
init_matrices(int mimax, int mjmax, int mkmax, Matrix *a, Matrix *b, Matrix *c, Matrix *p, Matrix *bnd, Matrix *wrk1,
              Matrix *wrk2) {

    newMat(p, 1, mimax, mjmax, mkmax);
    newMat(bnd, 1, mimax, mjmax, mkmax);
    newMat(wrk1, 1, mimax, mjmax, mkmax);
    newMat(wrk2, 1, mimax, mjmax, mkmax);
    newMat(a, 4, mimax, mjmax, mkmax);
    newMat(b, 3, mimax, mjmax, mkmax);
    newMat(c, 3, mimax, mjmax, mkmax);

    mat_set_init(p);
    mat_set(bnd, 0, 1.0f);
    mat_set(wrk1, 0, 0.0f);
    mat_set(wrk2, 0, 0.0f);
    mat_set(a, 0, 1.0f);
    mat_set(a, 1, 1.0f);
    mat_set(a, 2, 1.0f);
    mat_set(a, 3, 1.0f / 6.0f);
    mat_set(b, 0, 0.0f);
    mat_set(b, 1, 0.0f);
    mat_set(b, 2, 0.0f);
    mat_set(c, 0, 1.0f);
    mat_set(c, 1, 1.0f);
    mat_set(c, 2, 1.0f);
}

void verify_input_matrices(Matrix *a, Matrix *b, Matrix *c, Matrix *p, Matrix *bnd, Matrix *wrk1, Matrix *wrk2) {
	mat_verify_init(p);
    mat_verify(bnd, 0, 1.0f);
    mat_verify(wrk1, 0, 0.0f);
    mat_verify(wrk2, 0, 0.0f);
    mat_verify(a, 0, 1.0f);
    mat_verify(a, 1, 1.0f);
    mat_verify(a, 2, 1.0f);
    mat_verify(a, 3, 1.0f / 6.0f);
    mat_verify(b, 0, 0.0f);
    mat_verify(b, 1, 0.0f);
    mat_verify(b, 2, 0.0f);
    mat_verify(c, 0, 1.0f);
    mat_verify(c, 1, 1.0f);
    mat_verify(c, 2, 1.0f);	
}

void mat_verify_init(Matrix *Mat) {
	
    int i, j, k;

    for (i = 0; i < Mat->mrows; i++)
        for (j = 0; j < Mat->mcols; j++)
            for (k = 0; k < Mat->mdeps; k++) {
				if(MR(Mat, 0, i, j, k) != (float) (i * i) / (float) ((Mat->mrows - 1) * (Mat->mrows - 1))) {
					printf("Verification of input matrices (stage 1) failed!\n");
					exit(-1);
				}
			}
}



void mat_verify(Matrix *Mat, int l, float val) {
    
	int i, j, k;

    for (i = 0; i < Mat->mrows; i++)
        for (j = 0; j < Mat->mcols; j++)
            for (k = 0; k < Mat->mdeps; k++) {
				if(MR(Mat, l, i, j, k) != val) {
					printf("Verification of input matrices (stage 2) failed!\n");
					exit(-1);					
				}
			}
}



void write_outputs(float gosa, Matrix *wrk2) {

    char filename[FILENAME_SIZE];

    switch (matrix_dims[2]) {
        case 64:
            strncpy(filename, "outputs/xs.out", FILENAME_SIZE);
            break;
        case 128:
            strncpy(filename, "outputs/s.out", FILENAME_SIZE);
            break;
        case 256:
            strncpy(filename, "outputs/m.out", FILENAME_SIZE);
            break;
        case 512:
            strncpy(filename, "outputs/l.out", FILENAME_SIZE);
            break;
        case 1024:
            strncpy(filename, "outputs/xl.out", FILENAME_SIZE);
            break;
        default:
            break;
    }


    FILE *output_file = fopen(filename, "wb");
    if (output_file == NULL) {

        perror("Could not open output file !!\n");
        exit(1);
    }

    fwrite(&gosa, sizeof(float), 1, output_file);
    fwrite(wrk2->m, sizeof(float), matrix_size, output_file);

    fclose(output_file);
}

bool verify_outputs(float gosa, Matrix *matrix) {

    char filename[FILENAME_SIZE];
    int i, j, k;

    switch (matrix_dims[2]) {
        case 64:
            strncpy(filename, "outputs/xs.out", FILENAME_SIZE);
            break;
        case 128:
            strncpy(filename, "outputs/s.out", FILENAME_SIZE);
            break;
        case 256:
            strncpy(filename, "outputs/m.out", FILENAME_SIZE);
            break;
        case 512:
            strncpy(filename, "outputs/l.out", FILENAME_SIZE);
            break;
        case 1024:
            strncpy(filename, "outputs/xl.out", FILENAME_SIZE);
            break;
        default:
            break;
    }

    FILE *output_file = fopen(filename, "rb");
    if (output_file == NULL) {

        char text[] = "Could not open output file:";
        const size_t message_size = strlen(filename) + strlen(text) + 3;
        char *message = (char *) malloc(sizeof(char) * message_size);
        snprintf(message, message_size, "%s %s\n", text, filename);
        perror(message);

        return false;
    }

    /* read expected gosa */
    float expected_gosa;

    size_t res = fread(&expected_gosa, sizeof(float), 1, output_file);
    if (res != 1) {

        char text[] = "Failed to read float from output file:";
        const size_t message_size = strlen(filename) + strlen(text) + 3;
        char *message = (char *) malloc(sizeof(char) * message_size);
        snprintf(message, message_size, "%s %s\n", text, filename);
        perror(message);

        fclose(output_file);

        return false;
    }

	/* test output gosa */
	if (gosa < expected_gosa) {

        printf("Verification failed! Expected gosa to be AT LEAST an expected value\n");
        printf("  gosa     : %e\n", gosa);
        printf("  expected : %e\n", expected_gosa);

        fclose(output_file);

        return false;
    }
	

    /* read expected matrix */
    Matrix *expected_matrix = (Matrix *) malloc(sizeof(Matrix));
    newMat(expected_matrix, 1, matrix->mrows, matrix->mcols, matrix->mdeps);

    res = fread(expected_matrix->m, sizeof(float), matrix_size, output_file);
    if (res != matrix_size) {

        char text[] = "Failed to read floats from output file:";
        const size_t message_size = strlen(filename) + strlen(text) + 3;
        char *message = (char *) malloc(sizeof(char) * message_size);
        snprintf(message, message_size, "%s %s\n", text, filename);
        perror(message);

        fclose(output_file);
        clearMat(expected_matrix);
        free(expected_matrix);

        return false;
    }

    /* test output matrix */
    for (i = 0; i < matrix->mrows; i++) {
        for (j = 0; j < matrix->mcols; j++) {
            for (k = 0; k < matrix->mdeps; k++) {

                const float value = MR(matrix, 0, i, j, k);
                const float expected_value = MR(expected_matrix, 0, i, j, k);

                if (fabsf(value - expected_value) > EPSILON * value) {

                    printf("Verification failed!\n");
                    printf("  wrk2[%d,%d,%d]          : %e\n", i, j, k, value);
                    printf("  expected wrk2[%d,%d,%d] : %e\n", i, j, k, expected_value);

                    fclose(output_file);
                    clearMat(expected_matrix);
                    free(expected_matrix);

                    return false;
                }
            }
        }
    }

    printf("Verification successful!\n");

    fclose(output_file);
    clearMat(expected_matrix);
    free(expected_matrix);

    return true;
}

double get_reference_time() {

    switch (matrix_dims[2]) {
        case 64:
            return REFERENCE_MS[0];
        case 128:
            return REFERENCE_MS[1];
        case 256:
            return REFERENCE_MS[2];
        case 512:
            return REFERENCE_MS[3];
        case 1024:
            return REFERENCE_MS[4];
        default:
            return 0.0;
    }
}