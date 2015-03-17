#include "pixmapIO.h"


int main()
{
	int *image, dimx = 513, dimy = 511, type = 5;
	// int i, j;

	image = readPixmap("boat.pgm", &type, &dimx, &dimy);

	writePixmap(image, dimx, dimy, type, "boatXX.pgm");

/*for(i = 0 ; i < dimx ; i++)
	  for(j = 0 ; j < dimy ; j++)
	  printf("%d ", image[i * dimy + j]);	
*/

	return 0;
}
