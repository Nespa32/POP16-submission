#include "himeno_contest.h"
#include "pthread.h"

/** SOLUTION
  
  Names of the participants:
  GCC flags:
  SCORE in the ANTAREX machine:
  SCORE in the participants PC:
  Characteristics of participants PC (microprocessor, memory, clock frequency, etc.):
  
  
  Additional comments about your solution:
  
  
**/

#define THREAD_COUNT 4

typedef struct
{
    int startIdx;
    int endIdx;
    float gosa;
} t_args;

// globals, just for ease of use
Matrix* mat_a;
Matrix* mat_b;
Matrix* mat_c;
Matrix* mat_p;
Matrix* mat_bnd;
Matrix* mat_wrk1;
Matrix* mat_wrk2;

void *threadWorker(void* args)
{
    float const omega = 0.8;

    int const jmax = mat_p->mcols - 1;
    int const kmax = mat_p->mdeps - 1;

    t_args* t = (t_args*)args;
    int const startIdx = t->startIdx;
    int const endIdx = t->endIdx;
    float gosa = 0.0f;

    for (int i = startIdx; i < endIdx; i++)
    {
        for (int j = 1; j < jmax; j++)
        {
            for (int k = 1; k < kmax; k++)
            {
                float s0 = MR(mat_a, 0, i, j, k) * MR(mat_p, 0, i + 1, j, k)
                     + MR(mat_a, 1, i, j, k) * MR(mat_p, 0, i, j + 1, k)
                     + MR(mat_a, 2, i, j, k) * MR(mat_p, 0, i, j, k + 1)
                     + MR(mat_b, 0, i, j, k)
                       * (MR(mat_p, 0, i + 1, j + 1, k) - MR(mat_p, 0, i + 1, j - 1, k)
                          - MR(mat_p, 0, i - 1, j + 1, k) + MR(mat_p, 0, i - 1, j - 1, k))
                     + MR(mat_b, 1, i, j, k)
                       * (MR(mat_p, 0, i, j + 1, k + 1) - MR(mat_p, 0, i, j - 1, k + 1)
                          - MR(mat_p, 0, i, j + 1, k - 1) + MR(mat_p, 0, i, j - 1, k - 1))
                     + MR(mat_b, 2, i, j, k)
                       * (MR(mat_p, 0, i + 1, j, k + 1) - MR(mat_p, 0, i - 1, j, k + 1)
                          - MR(mat_p, 0, i + 1, j, k - 1) + MR(mat_p, 0, i - 1, j, k - 1))
                     + MR(mat_c, 0, i, j, k) * MR(mat_p, 0, i - 1, j, k)
                     + MR(mat_c, 1, i, j, k) * MR(mat_p, 0, i, j - 1, k)
                     + MR(mat_c, 2, i, j, k) * MR(mat_p, 0, i, j, k - 1)
                     + MR(mat_wrk1, 0, i, j, k);

                float ss = (s0 * MR(mat_a, 3, i, j, k) -
                    MR(mat_p, 0, i, j, k)) * MR(mat_bnd, 0, i, j, k);

                gosa += ss * ss;

                MR(mat_wrk2, 0, i, j, k) = MR(mat_p, 0, i, j, k) + omega * ss;
            }
        }
    }

    printf("startIdx %d, endIdx %d, gosa %f\n", startIdx, endIdx, gosa);
    t->gosa = gosa;
    return NULL;
}

/**
 * Jabobi Kernel
 */ 
float jacobi(Matrix *a, Matrix *b, Matrix *c, Matrix *p, Matrix *bnd,
    Matrix *wrk1, Matrix *wrk2) {

    mat_a = a;
    mat_b = b;
    mat_c = c;
    mat_p = p;
    mat_bnd = bnd;
    mat_wrk1 = wrk1;
    mat_wrk2 = wrk2;

    float gosa = 0.0f;

    pthread_t threads[THREAD_COUNT];
    t_args thread_args[THREAD_COUNT];

    if (p->mrows % THREAD_COUNT != 0)
    {
        fprintf(stderr, "Thread count %d not compatible with row division (rows %d)",
            THREAD_COUNT, p->mrows);
    }

    int work_per_thread = p->mrows / THREAD_COUNT;
    for (int i = 0; i < THREAD_COUNT; ++i)
    {
        thread_args[i].startIdx = i * work_per_thread;
        thread_args[i].endIdx = (i + 1) * work_per_thread;
        thread_args[i].gosa = 0.0f;
    }

    // first and last thread have a different workload
    thread_args[0].startIdx = 1;
    thread_args[THREAD_COUNT - 1].endIdx = p->mrows - 1;

    // thread work
    for (int i = 0; i < THREAD_COUNT; ++i)
        pthread_create(&threads[i], NULL, threadWorker, (void*)(&thread_args[i]));

    for (int i = 0; i < THREAD_COUNT; ++i)
    {
        pthread_join(threads[i], NULL);
        gosa += thread_args[i].gosa;
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

    size_t size = Mat->mrows * Mat->mcols * Mat->mdeps;
    float* ptr_start = (float*)(Mat->m + (size * l));
    for (size_t i = 0; i < size; ++i)
        ptr_start[i] = val;
}


