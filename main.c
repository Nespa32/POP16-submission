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

#include "himeno_contest.h"

#include "timer.h"

int main(int argc, char *argv[]) {

    Matrix a, b, c, p, bnd, wrk1, wrk2;
    float gosa;
    int mimax, mjmax, mkmax;

    /* Variable declaration for time measuring */
    long unsigned int iterations = 1;

    Timer *timer = timer_init();

    /* get matrix size */
    get_matrix_size(argc, argv, &mimax, &mjmax, &mkmax);

    /* alloc and init matrices */
    init_matrices(mimax, mjmax, mkmax, &a, &b, &c, &p, &bnd, &wrk1, &wrk2);

	/* verify if inputs are correct */
	verify_input_matrices(&a, &b, &c, &p, &bnd, &wrk1, &wrk2);
	
    printf("Now, start the actual measurement process\n");
    printf("Wait for a while...\n\n");

    /* single test run */
    timer_start(timer);
    gosa = jacobi(&a, &b, &c, &p, &bnd, &wrk1, &wrk2);
    timer_stop(timer);

    /* iterate if it doesn't run long enough */
    if (timer_get_ms(timer) < THRESHOLD_MS) {

        iterations = 10;

        while (iterations < MAXIMUM_ITERARIONS) {

            timer_start(timer);
            for (unsigned int iter_index = 0; iter_index < iterations; iter_index++) {

                gosa = jacobi(&a, &b, &c, &p, &bnd, &wrk1, &wrk2);
            }
            timer_stop(timer);

            if (timer_get_ms(timer) >= THRESHOLD_MS) {

                break;
            }

            iterations *= 2;
        }
    }

    /* report */
    double iteration_time = timer_get_ms(timer) / iterations;
    printf("Time  : %e ms\n", iteration_time);
    printf("Score : %u\n", (unsigned int) (get_reference_time() / iteration_time * 100.0));
    printf("Gosa  : %e \n\n", gosa);

    /* verify outputs */
    verify_outputs(gosa, &wrk2);

    /* free matrices */
    clear_matrices(&p, &bnd, &wrk1, &wrk2, &a, &b, &c);

    return (0);
}

