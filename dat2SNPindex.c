/*
 * dat2SNPindex.c
 *
 *
 * Converts text file containing a bitstring to indices of ones. 
 * 
 *
 * Ivo Siekmann, 
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h>

#define OUT stdout
#define ERR stderr

/* Reads a list of doubles. */
long bits2Indices(FILE *fp) {
  long k=0;
  int c;
  while ((c = fgetc(fp)) != EOF) {
    /* Ignore everything that is neither 0 nor 1 */
    if( c=='0') k++;
    else if (c=='1')  {k++; fprintf(OUT, "%li\n", k);}
  }
  return k;
}

int main(int argc, char **argv) {
  long nTrace;
  FILE *data_fp=stdin;

  fprintf(ERR, "Processing data...\n");
  nTrace=bits2Indices(data_fp);
  fprintf(ERR, "Exiting after reading %li ints...\n", nTrace);


  return 0;  
}
 
