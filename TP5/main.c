#include "pixmapIO.h"

#include <math.h>
#include <memory.h>

#include <assert.h>

#include "kmeans.h"


int main()
{
	int *image, dimx, dimy, type;

	image = readPixmap("znonstop01.ppm", &type, &dimx, &dimy);

	printf("type: %d\n", type);
	printf("dimx: %d\n", dimx);
	printf("dimy: %d\n", dimy);

	return 0;
}
