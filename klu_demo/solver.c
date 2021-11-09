#include "solver.h"

int klusolver(const int n, const int * Ap, const int * Ai, const double * Ax, double * x) {
  klu_common c;
  klu_symbolic * sy;
  klu_numeric * nu;

  // set defaults
  klu_defaults(&c);

  if (!Ap || !Ai || !Ax || !x)
    goto PrintError;

  if (!(sy = klu_analyze(n, Ap, Ai, c)))
    goto PrintError;

  if (!(nu = klu_factor(Ap, Ai, Ax, sy, c))) {
    klu_free_symbolic(&sy, c) ;
    goto PrintError;
  }

  // solver Ax=b, A is n-by-n and b is size-n and stored in x
  klu_solve(sy, nu, n, 1, x, c);
  
  // free & return sucess
  klu_free_numeric(&nu, c);
  klu_free_symbolic(&sy, c);
  return 1;

PrintError:
  printf ("error occurred!\n");
  return 0;
}

// TODO
cholmod_sparse * cholmod_generate(..., cholmod_common * cc) {
  int nrows, ncols, nnz;
  cholmod_triplet * T = cholmod_allocate_triplet(nrows, ncols, nnz, 1, CHOLMOD_REAL, cc);
  // fill the T->i T->j and T->x
  T->nnz = nnz;

  // convert triplet to sparse matrix
  cholmod_sparse * A = cholmod_triplet_to_sparse(T, nnz, cc);

  // free the triplet
  cholmod_free_triplet(&T, cc);

  // print the CSC format matrix A
  cholmod_print_sparse (A, "A", cc);

  return A;
}

// TODO
double * rhs_generate(int n, ...) {
  double * x = (double *) malloc(n * sizeof(double));
  if (x) {
    for (int i = 0; i < n; i++) {
      x[i] = 0;
    }
  } else {
    printf("failed to create the rhs\n");
    return NULL;
  }

  // build b here

  // fill the b into x
  for (int i = 0; i < bn; i++) {
    x[bidx[i]] = 1;
  }

  return x;
}

// TODO
int main(void) {
  // A common struct that cholmod always needs
  cholmod_common c;

  // start cholmod
  cholmod_start(&c);

  // generate the matrix A
  cholmod_sparse * A = cholmod_generate(..., &c);

  // generate the right hand side b
  double * x = rhs_generate(A->nrow, ...);

  if (A && x) {
    // solve Ax=b using klu
    klusolver(A->nrow, A->p, A->i, A->x, x);
    cholmod_free_sparse(&A, &c);
  }

  // some operations on x

  free(x);
  cholmod_finish(&c);

  return 0;
}

