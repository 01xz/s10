#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#define WORD_LEN   12
#define SHAPEN_LEN 6
#define RECT_SIZE  4
#define TPZD_SIZE  16
#define SHAPE_NAME "net"
#define INIT_SHAPE "net0"

#define DEBUG

int main(int argc, char **argv) {
  // the Input Path INDex and the Output Path INDex
  int ipind, opind;
  int ch;
  char word[WORD_LEN];
  char shpn[SHAPEN_LEN] = INIT_SHAPE;
  int shpnm = 0, netnm = 0;
  FILE * fp;

  // options parsing related
  int opt;
  int opt_idx;
  char * opstr = "i:o:";
  const struct option loptstr[] = {
    {"in",  required_argument, NULL, 'i'},
    {"out", required_argument, NULL, 'o'}
  };

  // quit if using the wrong options
  if ( argc != 5 ) {
    printf("Usage: fieldsolver2d -in input.data –out result.out\n");
    exit(EXIT_FAILURE);
  }

  // options parsing
  while ((opt = getopt_long_only(argc, argv, opstr, loptstr, &opt_idx)) != -1) {
#ifdef DEBUG
    printf("opt = %c\t\t", opt);
    printf("optarg = %s\t\t", optarg);
    printf("optind = %d\t\t", optind);
    printf("argv[optind] = %s\t\t", argv[optind]);
    printf("option_index = %d\n", opt_idx);
#endif
    if (opt == 'i') {
      ipind = optind - 1;
    } else {
      opind = optind - 1;
    }
  }
#ifdef DEBUG
  printf("the input file: %s\n", argv[ipind]);
  printf("the output file: %s\n", argv[opind]);
#endif

  // open the input file
  if ((fp = fopen(argv[ipind], "r")) == NULL) {
    printf("Failed to open the input data!\n");
    exit(EXIT_FAILURE);
  }

#ifdef DEBUG
  // check the input file
  printf("\ncheck the input file:\n");
  while ((ch = getc(fp)) != EOF) {
    printf("%c", ch);
  }
  rewind(fp);
#endif
  
  // nets counter
  while (fscanf(fp, "%s", word) == 1) {
    if (!strcmp(word, shpn)) {
      shpnm++;
      sprintf(shpn, "net%d", shpnm);
    }
    if (!strcmp(word, "net")) {
      netnm++;
    }
  }
  rewind(fp);
  sprintf(shpn, INIT_SHAPE);
#ifdef DEBUG
  printf("%d shapes and %d nets in total.\n", shpnm, netnm);
#endif

  // 'boundary' and 'dielectric'
  double * boundary = (double *) malloc(sizeof(double) * RECT_SIZE);
  double * dielectric = (double *) malloc(sizeof(double));
  // create 'nets' for recoding coords and 'shps' for recording different shapes
  double * nets = (double *) malloc(sizeof(double) * netnm * RECT_SIZE);
  int * shps = (int *) malloc(sizeof(int) * shpnm);

  // file parsing
  double * nets_p = nets;
  while (fscanf(fp, "%s", word) == 1) {
    if (!strcmp(word, "boundary")) {
      for (int i = 0; i < RECT_SIZE; i++) {
        fscanf(fp, "%lf", boundary + i);
      }
      continue;
    }
    if (!strcmp(word, "dielectric")) {
      fscanf(fp, "%lf", dielectric);
      continue;
    }
    if (!strcmp(word, "net")) {
      continue;
    }
    if (!strncmp(word, "net", 3)) {
      shps[atoi(word + 3)]++;
      for (int i = 0; i < RECT_SIZE; i++) {
        fscanf(fp, "%lf", nets_p + i);
      }
      nets_p += RECT_SIZE;
    }
  }

  fclose(fp);

#ifdef DEBUG
  // check the obtained data structure
  printf("check the obtained data structure:\n");
  printf("the boundary:\n");
  for (int i = 0; i < 4; i++) {
   printf("%f\t", *(boundary + i));
  }
  printf("\nthe dielectric:\n%f\n", *dielectric);
  printf("the shapes:\n");
  for (int i = 0; i < shpnm; i++) {
    printf("%d\t", shps[i]);
  }
  printf("\nthe nets:\n");
  for (int i = 0; i < netnm * RECT_SIZE; i++) {
    printf("%f\n", nets[i]);
  }
#endif

  return 0;
}
