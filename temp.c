#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

int main(int argc, char **argv) {
  // file operation related
  int ch;
  char * words;
  FILE * fp;
  char * infile, * outfile;

  // options parsing related
  int opt;
  int opt_idx;
  char * opstr = "i:o:";
  static struct option loptstr[] = {
    {"in",  required_argument, NULL, 'i'},
    {"out", required_argument, NULL, 'o'}
  };

  // obtained data
  double boundary[4] = {0.0};
  double dielectric  = 0.0;

  // quit if using the wrong options
  if ( argc != 5 ) {
    printf("Usage: fieldsolver2d -in input.data –out result.out\n");
    exit(EXIT_FAILURE);
  }

  // options parsing
  while ((opt = getopt_long_only(argc, argv, opstr, loptstr, &opt_idx)) != -1) {
    printf("opt = %c\t\t", opt);
    printf("optarg = %s\t\t", optarg);
    printf("optind = %d\t\t", optind);
    printf("argv[optind] = %s\t\t", argv[optind]);
    printf("option_index = %d\n", opt_idx);
    if (opt == 'i') {
      infile = optarg;
    } else {
      outfile = optarg;
    }
  }
  printf("infile: %s\n", infile);
  printf("outfile: %s\n", outfile);

  // open the input file
  if ((fp = fopen(infile, "r")) == NULL) {
    printf("Failed to open the input data!\n");
    exit(EXIT_FAILURE);
  }

  // check the input file
  while ((ch = getc(fp)) != EOF) {
    printf("%c", ch);
  }

  rewind(fp);

  // obtain data from the input file
  while (fscanf(fp, "%s", words) == 1) {
    printf("%s\n", words);
    //if (strcmp(words, "boundary") == 0) {
    //  for (int i = 0; i < 4; i++) {
    //    fscanf(fp, "%lf", boundary[i]);
    //  }
    //}
  }

  fclose(fp);

  //for (int i = 0; i < 4; i++) {
  //  printf("%f\n", boundary[i]);
  //}

  return 0;
}
