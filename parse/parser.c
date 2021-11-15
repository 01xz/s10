#include "parser.h"

lf * get_lf(FILE * fp) {
  char word[WORD_LEN];
  char shpn[SHAPEN_LEN] = INIT_SHAPE;

  lf * f = (lf *) malloc(sizeof(lf));
  if (!f) goto PrintError;

  // nets counter and shapes counter
  f->ns = 0, f->nn = 0;
  while (fscanf(fp, "%s", word) == 1) {
    if (!strcmp(word, shpn)) {
      f->ns++;
      if (f->ns < MAX_NS) {
        sprintf(shpn, SHAPEG_NAME "%d", f->ns);
      } else {
        printf("number of shapes are out of range\n");
        goto PrintError;
      }
    }
    if (!strcmp(word, SHAPEG_NAME)) {
      f->nn++;
    }
  }
  rewind(fp);

  // create memory for 'boundary' and 'dielectric'
  f->b = (double *) malloc(sizeof(double) * RECT_SIZE);
  f->d = (double *) malloc(sizeof(double));
  /**
   * create 'nets' for recoding coordinates and 'shps'
   * for recording different shapes
   */
  f->nets = (double *) malloc(sizeof(double) * f->nn * RECT_SIZE);
  f->shps = (int *) malloc(sizeof(int) * f->ns);

  if (!f->b || !f->d || !f->nets || !f->shps)
    goto PrintError;

  // file parsing
  double * nets = f->nets;
  while (fscanf(fp, "%s", word) == 1) {
    if (!strcmp(word, "boundary")) {
      for (int i = 0; i < RECT_SIZE; i++) {
        fscanf(fp, "%lf", f->b + i);
      }
      continue;
    }
    if (!strcmp(word, "dielectric")) {
      fscanf(fp, "%lf", f->d);
      continue;
    }
    if (!strcmp(word, SHAPEG_NAME)) {
      continue;
    }
    if (!strncmp(word, SHAPEG_NAME, SHAPEGN_LEN)) {
      (f->shps)[atoi(word + SHAPEGN_LEN)]++;
      for (int i = 0; i < RECT_SIZE; i++) {
        fscanf(fp, "%lf", nets + i);
      }
      nets += RECT_SIZE;
    }
  }

  return f;

PrintError:
  printf("error occurred in get_lf()\n");
  exit(EXIT_FAILURE);
}

void * free_lf(lf * f) {
  if (f != NULL) {
    if (f->b) {
      free(f->b);
      f->b = NULL;
    }
    if (f->d) {
      free(f->d);
      f->d = NULL;
    }
    if (f->nets) {
      free(f->nets);
      f->nets = NULL;
    }
    if (f->shps) {
      free(f->shps);
      f->shps = NULL;
    }
    free(f);
  }
  return NULL;
}

