#include <stdio.h>
#include <cholmod.h>
#include <time.h>

#define DEBUG

int main(void) {
  FILE * f_A, * f_b;

  // defination
  cholmod_sparse * A;
  cholmod_dense * x, * b;
  cholmod_factor * L;
  cholmod_common c;

#ifdef DEBUG
  // measure the running time
  clock_t startt, endt;
  startt = clock();
#endif

  // start cholmod
  cholmod_start(&c);

  // fill A and b
  if ((f_A = fopen ("data_A.mtx", "r")) != NULL) {
    printf("reading A:\n");
    A = cholmod_read_sparse(f_A, &c);
  } else {
    cholmod_free_sparse(&A, &c);
    cholmod_finish(&c);
    printf("error while reading A.\n");
    return 0;
  }

  if ((f_b = fopen ("data_b.mtx", "r")) != NULL) {
    printf("reading b:\n");
    b = cholmod_read_dense(f_b, &c);
  } else {
    cholmod_free_sparse(&A, &c);
    cholmod_free_dense(&b, &c);
    cholmod_finish(&c);
    printf("error while reading b.\n");
    return 0;
  }

  // solve
  L = cholmod_analyze(A, &c);
  cholmod_factorize(A, L, &c);
  x = cholmod_solve(CHOLMOD_A, L, b, &c);

  // print x
  cholmod_print_dense(x, "x", &c);

  // finish cholmod
  cholmod_free_factor(&L, &c);
  cholmod_free_sparse(&A, &c);
  cholmod_free_dense(&b, &c);
  cholmod_free_dense(&x, &c);
  cholmod_finish(&c);

#ifdef DEBUG
  endt = clock();
  printf("\ntime: %fs\n", (double) (endt - startt) / CLOCKS_PER_SEC);
#endif

  return 0;
}

