/* ========================================================================== */
/* === KLU DEMO ============================================================= */
/* ========================================================================== */

/* Read in a Matrix Market matrix (using CHOLMOD) and solve a linear system. */

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "klu.h"

/* for handling complex matrices */
#define REAL(X,i) (X [2*(i)])
#define IMAG(X,i) (X [2*(i)+1])
#define CABS(X,i) (sqrt (REAL (X,i) * REAL (X,i) + IMAG (X,i) * IMAG (X,i)))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

int klusolver(int n, int * Ap, int * Ai, double * Ax, double * x) {
  klu_common c;
  klu_symbolic * sy;
  klu_numeric * nu;

  // set defaults
  klu_defaults(&c);

  if (!Ap || !Ai || !Ax || !x)
    goto PrintError;

  if (!(sy = klu_analyze(n, Ap, Ai, &c)))
    goto PrintError;

  if (!(nu = klu_factor(Ap, Ai, Ax, sy, &c))) {
    klu_free_symbolic(&sy, &c) ;
    goto PrintError;
  }

  // solver Ax=b, A is n-by-n and b is size-n and stored in x
  klu_solve(sy, nu, n, 1, x, &c);
  
  // free & return sucess
  klu_free_numeric(&nu, &c);
  klu_free_symbolic(&sy, &c);
  return 1;

PrintError:
  printf ("error occurred!\n");
  return 0;
}

/* ========================================================================== */
/* === klu_backslash ======================================================== */
/* ========================================================================== */

static int klu_backslash    /* return 1 if successful, 0 otherwise */
(
    /* --- input ---- */
    int n,              /* A is n-by-n */
    int *Ap,            /* size n+1, column pointers */
    int *Ai,            /* size nz = Ap [n], row indices */
    double *Ax,         /* size nz, numerical values */
    int isreal,         /* nonzero if A is real, 0 otherwise */
    double *B,          /* size n, right-hand-side */

    /* --- output ---- */
    double *X,          /* size n, solution to Ax=b */
    double *R,          /* size n, residual r = b-A*x */

    /* --- scalar output --- */
    int *lunz,          /* nnz (L+U+F) */
    double *rnorm,      /* norm (b-A*x,1) / norm (A,1) */

    /* --- workspace - */

    klu_common *Common  /* default parameters and statistics */
)
{
    double anorm = 0, asum ;
    klu_symbolic *Symbolic ;
    klu_numeric *Numeric ;
    int i, j, p ;

    if (!Ap || !Ai || !Ax || !B || !X || !B) return (0) ;

    /* ---------------------------------------------------------------------- */
    /* symbolic ordering and analysis */
    /* ---------------------------------------------------------------------- */

    Symbolic = klu_analyze (n, Ap, Ai, Common) ;
    if (!Symbolic) return (0) ;

    if (isreal)
    {

        /* ------------------------------------------------------------------ */
        /* factorization */
        /* ------------------------------------------------------------------ */

        Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common) ;
        if (!Numeric)
        {
            klu_free_symbolic (&Symbolic, Common) ;
            return (0) ;
        }

        /* ------------------------------------------------------------------ */
        /* statistics (not required to solve Ax=b) */
        /* ------------------------------------------------------------------ */

        klu_rgrowth (Ap, Ai, Ax, Symbolic, Numeric, Common) ;
        klu_condest (Ap, Ax, Symbolic, Numeric, Common) ;
        klu_rcond (Symbolic, Numeric, Common) ;
        klu_flops (Symbolic, Numeric, Common) ;
        *lunz = Numeric->lnz + Numeric->unz - n + 
            ((Numeric->Offp) ? (Numeric->Offp [n]) : 0) ;

        /* ------------------------------------------------------------------ */
        /* solve Ax=b */
        /* ------------------------------------------------------------------ */

        for (i = 0 ; i < n ; i++)
        {
            X [i] = B [i] ;
        }
        klu_solve (Symbolic, Numeric, n, 1, X, Common) ;

        /* ------------------------------------------------------------------ */
        /* compute residual, rnorm = norm(b-Ax,1) / norm(A,1) */
        /* ------------------------------------------------------------------ */

        for (i = 0 ; i < n ; i++)
        {
            R [i] = B [i] ;
        }
        for (j = 0 ; j < n ; j++)
        {
            asum = 0 ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                /* R (i) -= A (i,j) * X (j) */
                R [Ai [p]] -= Ax [p] * X [j] ;
                asum += fabs (Ax [p]) ;
            }
            anorm = MAX (anorm, asum) ;
        }
        *rnorm = 0 ;
        for (i = 0 ; i < n ; i++)
        {
            *rnorm = MAX (*rnorm, fabs (R [i])) ;
        }

        /* ------------------------------------------------------------------ */
        /* free numeric factorization */
        /* ------------------------------------------------------------------ */

        klu_free_numeric (&Numeric, Common) ;

    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* statistics (not required to solve Ax=b) */
        /* ------------------------------------------------------------------ */

        Numeric = klu_z_factor (Ap, Ai, Ax, Symbolic, Common) ;
        if (!Numeric)
        {
            klu_free_symbolic (&Symbolic, Common) ;
            return (0) ;
        }

        /* ------------------------------------------------------------------ */
        /* statistics */
        /* ------------------------------------------------------------------ */

        klu_z_rgrowth (Ap, Ai, Ax, Symbolic, Numeric, Common) ;
        klu_z_condest (Ap, Ax, Symbolic, Numeric, Common) ;
        klu_z_rcond (Symbolic, Numeric, Common) ;
        klu_z_flops (Symbolic, Numeric, Common) ;
        *lunz = Numeric->lnz + Numeric->unz - n + 
            ((Numeric->Offp) ? (Numeric->Offp [n]) : 0) ;

        /* ------------------------------------------------------------------ */
        /* solve Ax=b */
        /* ------------------------------------------------------------------ */

        for (i = 0 ; i < 2*n ; i++)
        {
            X [i] = B [i] ;
        }
        klu_z_solve (Symbolic, Numeric, n, 1, X, Common) ;

        /* ------------------------------------------------------------------ */
        /* compute residual, rnorm = norm(b-Ax,1) / norm(A,1) */
        /* ------------------------------------------------------------------ */

        for (i = 0 ; i < 2*n ; i++)
        {
            R [i] = B [i] ;
        }
        for (j = 0 ; j < n ; j++)
        {
            asum = 0 ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                /* R (i) -= A (i,j) * X (j) */
                i = Ai [p] ;
                REAL (R,i) -= REAL(Ax,p) * REAL(X,j) - IMAG(Ax,p) * IMAG(X,j) ;
                IMAG (R,i) -= IMAG(Ax,p) * REAL(X,j) + REAL(Ax,p) * IMAG(X,j) ;
                asum += CABS (Ax, p) ;
            }
            anorm = MAX (anorm, asum) ;
        }
        *rnorm = 0 ;
        for (i = 0 ; i < n ; i++)
        {
            *rnorm = MAX (*rnorm, CABS (R, i)) ;
        }

        /* ------------------------------------------------------------------ */
        /* free numeric factorization */
        /* ------------------------------------------------------------------ */

        klu_z_free_numeric (&Numeric, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* free symbolic analysis, and residual */
    /* ---------------------------------------------------------------------- */

    klu_free_symbolic (&Symbolic, Common) ;
    return (1) ;
}


/* ========================================================================== */
/* === klu_demo ============================================================= */
/* ========================================================================== */

/* Given a sparse matrix A, set up a right-hand-side and solve X = A\b */

static void klu_demo (int n, int *Ap, int *Ai, double *Ax, int isreal)
{
    FILE * f ;
    clock_t startt, endt;
    int cd;
    double rnorm ;
    klu_common Common ;
    double *B, *X, *R ;
    int i, lunz ;

    printf ("KLU: %s, version: %d.%d.%d\n", KLU_DATE, KLU_MAIN_VERSION,
        KLU_SUB_VERSION, KLU_SUBSUB_VERSION) ;

    /* ---------------------------------------------------------------------- */
    /* set defaults */
    /* ---------------------------------------------------------------------- */

    klu_defaults (&Common) ;

    /* ---------------------------------------------------------------------- */
    /* create a right-hand-side */
    /* ---------------------------------------------------------------------- */

    if (isreal)
    {
        /* B = 1 + (1:n)/n */
        B = klu_malloc (n, sizeof (double), &Common) ;
        X = klu_malloc (n, sizeof (double), &Common) ;
        R = klu_malloc (n, sizeof (double), &Common) ;
        if (B)
        {
            for (i = 0 ; i < n ; i++)
            {
                B [i] = 0;
            }
            if ((f = fopen ("../cb", "r")) != NULL) {
                while (fscanf(f, "%d", &cd) == 1) {
                    B[cd] = 1;
                }
            } else {
                printf("error while reading b.\n");
            }
            fclose(f);
        }
    }
    else
    {
        /* real (B) = 1 + (1:n)/n, imag(B) = (n:-1:1)/n */
        B = klu_malloc (n, 2 * sizeof (double), &Common) ;
        X = klu_malloc (n, 2 * sizeof (double), &Common) ;
        R = klu_malloc (n, 2 * sizeof (double), &Common) ;
        if (B)
        {
            for (i = 0 ; i < n ; i++)
            {
                REAL (B, i) = 1 + ((double) i+1) / ((double) n) ;
                IMAG (B, i) = ((double) n-i) / ((double) n) ;
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* X = A\b using KLU and print statistics */
    /* ---------------------------------------------------------------------- */

    // record start time
    startt = clock();

    if (!klu_backslash (n, Ap, Ai, Ax, isreal, B, X, R, &lunz, &rnorm, &Common))
    {
        printf ("KLU failed\n") ;
    }
    else
    {
        printf ("n %d nnz(A) %d nnz(L+U+F) %d resid %g\n"
            "recip growth %g condest %g rcond %g flops %g\n",
            n, Ap [n], lunz, rnorm, Common.rgrowth, Common.condest,
            Common.rcond, Common.flops) ;
    }

    // record end time
    endt = clock();
    printf("\nUsed time: %fs\n", (double) (endt - startt) / CLOCKS_PER_SEC);

    if ((f = fopen ("x", "w")) != NULL) {
        printf("writing x to file:\n");
    } else {
        printf("error while creating file x.\n");
    }
    // print x
    for (int i = 0; i < n; i++) {
        fprintf(f, "%g\n", *(X + i));
    }
    fclose(f);
    printf("done\n");

    /* ---------------------------------------------------------------------- */
    /* free the problem */
    /* ---------------------------------------------------------------------- */

    if (isreal)
    {
        klu_free (B, n, sizeof (double), &Common) ;
        klu_free (X, n, sizeof (double), &Common) ;
        klu_free (R, n, sizeof (double), &Common) ;
    }
    else
    {
        klu_free (B, 2*n, sizeof (double), &Common) ;
        klu_free (X, 2*n, sizeof (double), &Common) ;
        klu_free (R, 2*n, sizeof (double), &Common) ;
    }
    printf ("peak memory usage: %g bytes\n\n", (double) (Common.mempeak)) ;
}


/* ========================================================================== */
/* === main ================================================================= */
/* ========================================================================== */

/* Read in a sparse matrix in Matrix Market format using CHOLMOD, and then
 * solve Ax=b with KLU.  Note that CHOLMOD is only used to read the matrix. */

#include "cholmod.h"

int main (void)
{
    cholmod_sparse *A ;
    cholmod_common ch ;
    cholmod_start (&ch) ;
    A = cholmod_read_sparse (stdin, &ch) ;
    FILE * f ;
    clock_t startt, endt;
    if (A)
    {
        if (A->nrow != A->ncol || A->stype != 0
            || (!(A->xtype == CHOLMOD_REAL || A->xtype == CHOLMOD_COMPLEX)))
        {
            printf ("invalid matrix\n") ;
        }
        else
        {
            int n = A->nrow;
            int cd;

            double * x = malloc(n * sizeof(double)) ;
            for (int i = 0; i < n; i++) {
                *(x + i) = 0;
            }
            if ((f = fopen ("../cb", "r")) != NULL) {
                while (fscanf(f, "%d", &cd) == 1) {
                    x[cd] = 1;
                }
            } else {
                printf("error while reading b.\n");
            }
            fclose(f);

            // record start time
            startt = clock();

            klusolver(A->nrow, A->p, A->i, A->x, x) ;

            // record end time
            endt = clock();
            printf("\nUsed time: %fs\n", (double) (endt - startt) / CLOCKS_PER_SEC);

            if ((f = fopen ("x", "w")) != NULL) {
                printf("writing x to file:\n");
            } else {
                printf("error while creating file x.\n");
            }
            // print x
            for (int i = 0; i < n; i++) {
                fprintf(f, "%g\n", *(x + i));
            }
            fclose(f);
            printf("done\n");

        }
        cholmod_free_sparse (&A, &ch) ;
    }
    cholmod_finish (&ch) ;
    return (0) ;
}
