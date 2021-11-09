#include <stdio.h>
#include <cs.h>
#include <time.h>

int main(void) {
  FILE * f;
  cs * T;
  clock_t startt, endt;

  // read matrix from file
  if ((f = fopen ("A", "r")) != NULL) {
    printf("reading A:\n");
    T = cs_load(f);
  } else {
    printf("error while reading A.\n");
    return 1;
  }
  fclose(f);

  // the size of the matrix
  int size = T -> m;

  // fill the right side hand b
  double * b = (double *) malloc (sizeof(double) * size);
  for (int i = 0; i < size; i++) {
    b[i] = 0;
  }

  // coordinates of the right hand side b
  int cd;

  // read right side hand from file
  if ((f = fopen ("b", "r")) != NULL) {
    while (fscanf(f, "%d", &cd) == 1) {
      b[cd] = 1;
    }
  } else {
    printf("error while reading b.\n");
    return 1;
  }
  fclose(f);

  // convert T from triplet form into compressed-column matrix
  cs * A = cs_compress(T);
  cs_spfree(T);

  // solve x
  if ((f = fopen ("x", "w")) != NULL) {
    printf("start solving x:\n");
  } else {
    printf("error while creating file x.\n");
  }

  // record start time
  startt = clock();

  if (!cs_lusol(3, A, b, 1)) {
    printf("error when solving x\n");
  }

  // recore end time
  endt = clock();
  printf("\nUsed time: %fs\n", (double) (endt - startt) / CLOCKS_PER_SEC);

  // write x to file
  printf("\nstart writing x:\n");
  for (int i = 0; i < size; i++) {
    fprintf(f, "%g\n", *(b + i));
  }
  fclose(f);
  printf("done\n");

  // free
  cs_spfree(A);
  free(b);
  
  return 0;
}

