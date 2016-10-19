#include "himeno_contest.h"

/** SOLUTION
  
  Names of the participants:
  GCC flags:
  SCORE in the ANTAREX machine:
  SCORE in the participants PC:
  Characteristics of participants PC (microprocessor, memory, clock frequency, etc.):
  
  
  Additional comments about your solution:
  
  
**/


/**
 * Jabobi Kernel
 */ 
float jacobi(Matrix *a, Matrix *b, Matrix *c, Matrix *p, Matrix *bnd, Matrix *wrk1, Matrix *wrk2) {

    int i, j, k, imax, jmax, kmax;
    float gosa, s0, ss;
    const float omega = 0.8;

    imax = p->mrows - 1;
    jmax = p->mcols - 1;
    kmax = p->mdeps - 1;

    gosa = 0.0;

    for (i = 1; i < imax; i++) {
        for (j = 1; j < jmax; j++) {
            for (k = 1; k < kmax; k++) {
                s0 = MR(a, 0, i, j, k) * MR(p, 0, i + 1, j, k)
                     + MR(a, 1, i, j, k) * MR(p, 0, i, j + 1, k)
                     + MR(a, 2, i, j, k) * MR(p, 0, i, j, k + 1)
                     + MR(b, 0, i, j, k)
                       * (MR(p, 0, i + 1, j + 1, k) - MR(p, 0, i + 1, j - 1, k)
                          - MR(p, 0, i - 1, j + 1, k) + MR(p, 0, i - 1, j - 1, k))
                     + MR(b, 1, i, j, k)
                       * (MR(p, 0, i, j + 1, k + 1) - MR(p, 0, i, j - 1, k + 1)
                          - MR(p, 0, i, j + 1, k - 1) + MR(p, 0, i, j - 1, k - 1))
                     + MR(b, 2, i, j, k)
                       * (MR(p, 0, i + 1, j, k + 1) - MR(p, 0, i - 1, j, k + 1)
                          - MR(p, 0, i + 1, j, k - 1) + MR(p, 0, i - 1, j, k - 1))
                     + MR(c, 0, i, j, k) * MR(p, 0, i - 1, j, k)
                     + MR(c, 1, i, j, k) * MR(p, 0, i, j - 1, k)
                     + MR(c, 2, i, j, k) * MR(p, 0, i, j, k - 1)
                     + MR(wrk1, 0, i, j, k);

                ss = (s0 * MR(a, 3, i, j, k) - MR(p, 0, i, j, k)) * MR(bnd, 0, i, j, k);

                gosa += ss * ss;

                MR(wrk2, 0, i, j, k) = MR(p, 0, i, j, k) + omega * ss;
            }
        }
    }

    return (gosa);
}

/**
 * Matrix initialization - Type 1 
 */
void mat_set_init(Matrix *Mat) {

    int i, j, k;

    for (i = 0; i < Mat->mrows; i++)
        for (j = 0; j < Mat->mcols; j++)
            for (k = 0; k < Mat->mdeps; k++)
                MR(Mat, 0, i, j, k) = (float) (i * i)
                                      / (float) ((Mat->mrows - 1) * (Mat->mrows - 1));
}

/**
 * Matrix initialization - Type 2 
 */
void mat_set(Matrix *Mat, int l, float val) {
    int i, j, k;

    for (i = 0; i < Mat->mrows; i++)
        for (j = 0; j < Mat->mcols; j++)
            for (k = 0; k < Mat->mdeps; k++)
                MR(Mat, l, i, j, k) = val;
}


