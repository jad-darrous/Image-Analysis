/* Headers for PBM+ */
/* E.B.             */

#ifndef   	_UTIL_
#define   	_UTIL_

#include <stdio.h>
#include <stdlib.h>

typedef unsigned char bit;
typedef unsigned char gray;

char pm_getc(FILE* file);
bit pm_getbit(FILE* file);
unsigned char pm_getrawbyte(FILE* file);
int pm_getint( FILE* file);

void pm_erreur(char *);

#endif