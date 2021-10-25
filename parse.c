#include <stdio.h>
#include <getopt.h>

int main(int argc, char **argv) {
  char * infile, * outfile;
  int opt;
  int opt_idx;
  char * opstr = "i:o:";
  static struct option loptstr[] = {
    {"in",  required_argument, NULL, 'i'},
    {"out", required_argument, NULL, 'o'}
  };
  
  printf("%d\n", argc);
  while ((opt = getopt_long_only(argc, argv, opstr, loptstr, &opt_idx)) != -1) {
    printf("opt = %c\t\t", opt);
    printf("optarg = %s\t\t", optarg);
    printf("optind = %d\t\t", optind);
    printf("argv[optind] = %s\t\t", argv[optind]);
    printf("option_index = %d\n", opt_idx);
    switch (opt) {
      case 'i':
        infile = optarg;
        printf("infile: %s\n", infile);
        break;
      case 'o':
        outfile = optarg;
        printf("outfile: %s\n", outfile);
        break;
      default:
        printf("Usage: fieldsolver2d -in input.data –out result.out\n");
    }
  }
  return 0;
}
