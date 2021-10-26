#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#define DEBUG

#define WORDS_LEN 12
#define NETSN_LEN 6

int main(int argc, char **argv) {
  // the Input Path INDex and the Output Path INDex
  int ipind, opind;
  int ch;
  char words[WORDS_LEN];
  char netsn[NETSN_LEN] = "net0";
  int netsnm = 0;
  FILE * fp;

  // options parsing related
  int opt;
  int opt_idx;
  char * opstr = "i:o:";
  const struct option loptstr[] = {
    {"in",  required_argument, NULL, 'i'},
    {"out", required_argument, NULL, 'o'}
  };

  // data
  double boundary[4] = {0.0};
  double dielectric  = 0.0;

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
  while ((ch = getc(fp)) != EOF) {
    printf("%c", ch);
  }
  rewind(fp);
#endif
  
  // nets counter
  while (fscanf(fp, "%s", words) == 1) {
    while (strcmp(words, netsn))
  }
  rewind(fp);

  // file parsing
  while (fscanf(fp, "%s", words) == 1) {
    if (!strcmp(words, "boundary")) {
      printf("obtained: %s\n", words);
      for (int i = 0; i < 4; i++) {
        fscanf(fp, "%lf", boundary + i);
      }
      continue;
    }
    if (!strcmp(words, "dielectric")) {
      printf("obtained: %s\n", words);
      fscanf(fp, "%lf", &dielectric);
      continue;
    }
    if (!strcmp(words, "net")) {
      printf("obtained: %s\n", words);
      continue;
    }
    if (!strcmp(words, netsn)) {
      netsnm++;
      strcpy();
    }
  }

  fclose(fp);

#ifdef DEBUG
  // check the obtained data structure
  printf("check the obtained data structure:\n");
  printf("the boundary:\n");
  for (int i = 0; i < 4; i++) {
   printf("%f\t", boundary[i]);
  }
  printf("\nthe dielectric:\n%f\n", dielectric);
#endif

  return 0;
}
