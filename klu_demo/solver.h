#ifndef _SOLVER_H
#define _SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "cholmod.h"
#include "klu.h"

int klusolver(int n, int *Ap, int *Ai, double *Ax, double *x);
cholmod_sparse * cholmod_generate(..., cholmod_common * cc);
double * rhs_generate(int n, ...);

#ifdef __cplusplus
}
#endif

#endif

