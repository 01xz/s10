#include "parser.h"

int main(int argc, char **argv) {

#ifdef DEBUG
  // measure the running time
  clock_t startt, endt;
  startt = clock();
#endif

  // the Input Path INDex and the Output Path INDex
  int ipind, opind;
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
  if (argc != 5) {
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

  lf * field = get_lf(fp);

  fclose(fp);

#ifdef DEBUG
  // check the obtained data structure
  printf("check the obtained data structure:\n");
  printf("the boundary:\n");
  for (int i = 0; i < 4; i++) {
    printf("%f\t", (field->b)[i]);
  }
  printf("\nthe dielectric:\n%f\n", *(field->d));
  printf("the shapes:\n");
  for (int i = 0; i < field->ns; i++) {
    printf("%d\t", (field->shps)[i]);
  }
  printf("\nthe nets:\n");
  for (int i = 0; i < (field->nn) * RECT_SIZE; i++) {
    printf("%f\n", (field->nets)[i]);
  }
#endif

  // TODO: insert the solver here:

  // create the output data
  if ((fp = fopen(argv[opind], "w")) == NULL) {
    printf("Failed to create the output data!\n");
    exit(EXIT_FAILURE);
  }

#ifdef DEBUG
  // output data test
  fprintf(fp, "# length = 1um\n");
  fprintf(fp, "# master = net0\n");
  for (int i = 0; i < (field->ns); i++) {
    fprintf(fp, "%*s" SHAPEG_NAME "%d", 10, "", i);
  }
  fprintf(fp, "\n " INIT_SHAPE ":");
  for (int i = 0; i < (field->ns); i++) {
    fprintf(fp, "%*s%gfF", i ? 3 : 1, "", 0.00107424);
  }
  fprintf(fp, "\n");
#endif

  fclose(fp);

#ifdef DEBUG
  endt = clock();
  printf("\ntime: %fs\n", (double) (endt - startt) / CLOCKS_PER_SEC);
#endif

  return 0;
}

