#ifndef _SOLVER_H
#define _SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "cholmod.h"
#include "klu.h"


/**
 * klusolver - solve the linear system Ax = b, the sparse matrix
 * A is squared ccs-format sparse matrix and the right-hand-side
 * b is stored in x
 * @n: the rows(or columns) of the squared matrix A
 * @Ap: the column pointers of A
 * @Ai: the row indices
 * @Ax: the elements in A
 * @x: the right-hand-side b before solving
 *
 * klusolver will write the solution back to x and it will return
 * 0 if fails to solve the equation
 */
int klusolver(int n, int *Ap, int *Ai, double *Ax, double *x);

/**
 * cholmod_generate - generate the squared sparse matrix A using
 * the cholmod package
 * @cc: the pointer to cholmod_common object
 *
 */
cholmod_sparse * cholmod_generate(..., cholmod_common * cc);

/**
 * rhs_generate - generate the right-hand-side b for solving the
 * linear system Ax=b
 * @n: the size of b
 *
 */
double * rhs_generate(int n, ...);

#ifdef __cplusplus
}
#endif

#endif

