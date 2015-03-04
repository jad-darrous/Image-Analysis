#ifndef   	_IMG_UTIL_
#define   	_IMG_UTIL_

#include <stdio.h>
#include <stdlib.h>
#include "Util.h"

typedef struct _image_header {
	char type;
	int rows, cols;
	int maxval;
	int pgmraw;
} image_header;

typedef struct _gs_image {
	image_header* hdr;
	gray* data; 
} gs_image;

struct cl_image {
	image_header* hdr;
	int*** m; 
};


void print_image_header(image_header* hdr);

gs_image* read_image(const char* filename);

void write_image(const char* filename, gs_image* img);

#endif