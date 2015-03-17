#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

int * readPixmap(char* filename,
		 int * type,
		 int * nbColumns,
		 int * nbRows);

int writePixmap(int * pixmap,
		int cols, int rows,
		int type,
		char * filename);
