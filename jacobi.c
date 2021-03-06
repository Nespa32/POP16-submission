#include "himeno_contest.h"
#include <unistd.h>
#include "pthread.h"

/** SOLUTION
  
  Names of the participants: Ricardo Leite, Andre Cascais
  GCC flags: -O3 -pthread
  SCORE in the ANTAREX machine: 3100
  SCORE in the participants PC: 1700
  Characteristics of participants PC (microprocessor, memory, clock frequency, etc.):
  It's a VM :D, 4 'dedicated' logical cores assigned out of 4 total,
    1x Intel i7-4700MQ 2.40GHz, 2GB RAM
  
  Additional comments about your solution:
  
  
**/

#define THREAD_COUNT 32

typedef struct
{
    int startIdx;
    int endIdx;
    float gosa;
    pthread_cond_t cond_var;
    pthread_mutex_t cond_var_lock;
} t_args;

// globals, just for ease of use
Matrix* mat_a;
Matrix* mat_b;
Matrix* mat_c;
Matrix* mat_p;
Matrix* mat_bnd;
Matrix* mat_wrk1;
Matrix* mat_wrk2;

pthread_mutex_t thread_count_lock;
int done_thread_count = 0;

void *threadWorker(void* args)
{
    float const omega = 0.8;

    t_args* t = (t_args*)args;

    while (true)
    {
        int const jmax = mat_p->mcols - 1;
        int const kmax = mat_p->mdeps - 1;

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
                        + MR(mat_c, 2, i, j, k) * MR(mat_p, 0, i, j, k - 1)
                        + MR(mat_c, 0, i, j, k) * MR(mat_p, 0, i - 1, j, k)
                        + MR(mat_c, 1, i, j, k) * MR(mat_p, 0, i, j - 1, k)
                        + MR(mat_b, 0, i, j, k)
                          * (MR(mat_p, 0, i + 1, j + 1, k) - MR(mat_p, 0, i + 1, j - 1, k)
                            - MR(mat_p, 0, i - 1, j + 1, k) + MR(mat_p, 0, i - 1, j - 1, k))
                        + MR(mat_b, 1, i, j, k)
                          * (MR(mat_p, 0, i, j + 1, k + 1) - MR(mat_p, 0, i, j - 1, k + 1)
                            - MR(mat_p, 0, i, j + 1, k - 1) + MR(mat_p, 0, i, j - 1, k - 1))
                        + MR(mat_b, 2, i, j, k)
                          * (MR(mat_p, 0, i + 1, j, k + 1) - MR(mat_p, 0, i - 1, j, k + 1)
                            - MR(mat_p, 0, i + 1, j, k - 1) + MR(mat_p, 0, i - 1, j, k - 1))
                        + MR(mat_wrk1, 0, i, j, k);

                    float ss = (s0 * MR(mat_a, 3, i, j, k) -
                        MR(mat_p, 0, i, j, k)) * MR(mat_bnd, 0, i, j, k);

                    gosa += ss * ss;

                    MR(mat_wrk2, 0, i, j, k) = MR(mat_p, 0, i, j, k) + omega * ss;
                }
            }
        }

        t->gosa = gosa;

        {
            pthread_mutex_lock(&thread_count_lock);
            ++done_thread_count;
            pthread_mutex_unlock(&thread_count_lock);
        }

        pthread_cond_wait(&t->cond_var, &t->cond_var_lock);
    }

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

    static float gosa = 0.0f;
    static int is_pool_init = 0;
    static pthread_t threads[THREAD_COUNT];
    static t_args thread_args[THREAD_COUNT];

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

    // init/activate threads
    done_thread_count = 0;
    if (is_pool_init == 0)
    {
        is_pool_init = 1;
        for (int i = 0; i < THREAD_COUNT; ++i)
        {
            pthread_cond_init(&thread_args[i].cond_var, NULL);
            pthread_mutex_init(&thread_args[i].cond_var_lock, NULL);
            pthread_create(&threads[i], NULL, threadWorker, (void*)(&thread_args[i]));
        }
    }
    else
    {
        // poke threads for a new job
        for (int i = 0; i < THREAD_COUNT; ++i)
            pthread_cond_signal(&thread_args[i].cond_var);
    }

    // 'busy' wait
    while (done_thread_count < THREAD_COUNT)
        usleep(1000);

    for (int i = 0; i < THREAD_COUNT; ++i)
        gosa += thread_args[i].gosa;

    return (gosa);
}

/**
 * Matrix initialization - Type 1 
 */
void mat_set_init(Matrix *Mat) {

    int i, j, k;

    for (i = 0; i < Mat->mrows; i++)
    {
        float sqr_i = (float)(i * i);
        for (j = 0; j < Mat->mcols; j++)
            for (k = 0; k < Mat->mdeps; k++)
                MR(Mat, 0, i, j, k) = sqr_i
                                      / (float) ((Mat->mrows - 1) * (Mat->mrows - 1));
    }
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


