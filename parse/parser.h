#ifndef _PARSER_H
#define _PARSER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

#define WORD_LEN    12
#define SHAPEN_LEN  15
#define MAX_NS      100
#define RECT_SIZE   4
#define TPZD_SIZE   16
#define SHAPEG_NAME "net"
#define SHAPEGN_LEN strlen(SHAPEG_NAME)
#define INIT_SHAPE  SHAPEG_NAME "0"

#define DEBUG

typedef struct layout_field {
  // the number of nets
  int nn; 
  // the number of shapes
  int ns;
  // the boundary of field
  double * b;
  // the dielectric
  double * d;
  // the nets list
  double * nets;
  // the shapes indicator
  int    * shps;
} lf;

/**
 * get_lf - get the field infomation from the input.data 
 * @fp: the file pointer of the input.data
 *
 * get_lf will not close @fp, you may have to close @fp
 * manually
 */
lf * get_lf(FILE * fp);

/**
 * free_lf - free the lf pointer
 * @f: the lf pointer to be freed
 *
 * usage: f = free_lf(f)
 */
void * free_lf(lf * f);

#ifdef __cplusplus
}
#endif

#endif
