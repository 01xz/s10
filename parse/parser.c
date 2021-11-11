#include "parser.h"

lf * get_lf(FILE * fp) {
  lf * field = (lf *) malloc(sizeof(lf));
  int ns = 0, nn = 0;
  char word[WORD_LEN];
  char shpn[SHAPEN_LEN] = INIT_SHAPE;

#ifdef DEBUG
  // check the input file
  int ch;
  printf("\ncheck the input file:\n");
  while ((ch = getc(fp)) != EOF) {
    printf("%c", ch);
  }
  rewind(fp);
#endif
  
  // nets counter
  while (fscanf(fp, "%s", word) == 1) {
    if (!strcmp(word, shpn)) {
      ns++;
      sprintf(shpn, SHAPEG_NAME "%d", ns);
    }
    if (!strcmp(word, SHAPEG_NAME)) {
      nn++;
    }
  }
  rewind(fp);
#ifdef DEBUG
  printf("%d shapes and %d nets in total.\n", ns, nn);
#endif

  // create 'boundary' and 'dielectric'
  double * b = (double *) malloc(sizeof(double) * RECT_SIZE);
  double * d = (double *) malloc(sizeof(double));
  // create 'nets' for recoding coordinates and 'shps' for recording different shapes
  double * nets = (double *) malloc(sizeof(double) * nn * RECT_SIZE);
  int * shps = (int *) malloc(sizeof(int) * ns);

  // file parsing
  double * nets_p = nets;
  while (fscanf(fp, "%s", word) == 1) {
    if (!strcmp(word, "boundary")) {
      for (int i = 0; i < RECT_SIZE; i++) {
        fscanf(fp, "%lf", b + i);
      }
      continue;
    }
    if (!strcmp(word, "dielectric")) {
      fscanf(fp, "%lf", d);
      continue;
    }
    if (!strcmp(word, SHAPEG_NAME)) {
      continue;
    }
    if (!strncmp(word, SHAPEG_NAME, SHAPEGN_LEN)) {
      shps[atoi(word + SHAPEGN_LEN)]++;
      for (int i = 0; i < RECT_SIZE; i++) {
        fscanf(fp, "%lf", nets_p + i);
      }
      nets_p += RECT_SIZE;
    }
  }

  // fill in field
  field->nn   = nn;
  field->ns   = ns;
  field->b    = b;
  field->d    = d;
  field->nets = nets;
  field->shps = shps;

  return field;
}
