#include "pixmapIO.h"

#include <math.h>
#include <memory.h>

#include <assert.h>

double B22[9] = {1/16.0,2/16.0,1/16.0, 2/16.0,4/16.0,2/16.0, 1/16.0,2/16.0,1/16.0};
double sobel_h1[9] = {-1.0/4.0, 0, 1.0/4.0,   -2.0/4.0, 0, 2.0/4.0,  -1.0/4.0, 0, 1.0/4.0 };
double sobel_h2[9] = {-1.0/4.0, -2.0/4.0, -1.0/4.0,   0, 0, 0,  1.0/4.0, 2.0/4.0, 1.0/4.0 };

#define GC 159.0
double G55[25] = { 
	2/GC, 4/GC, 5/GC, 4/GC, 2/GC,  
	4/GC, 9/GC, 12/GC, 9/GC, 4/GC,
	5/GC, 12/GC, 15/GC, 12/GC, 5/GC,
	4/GC, 9/GC, 12/GC, 9/GC, 4/GC,
	2/GC, 4/GC, 5/GC, 4/GC, 2/GC};

int* double_to_int(double *arr, int size) {

	int i;
	int* data = (int *) malloc(size * sizeof(int));
	for (i=0; i < size; i++)
		data[i] = (int) arr[i];
	return data; 
}

double *int_to_double(int *arr, int size) {
	int i;
	double* data = (double *) malloc(size * sizeof(double));
	for (i=0; i < size; i++)
		data[i] = (double) arr[i];
	return data; 
}

void print_int_arr(int *img, int cols, int rows) {
	int i, j;
	for(i=0; i < rows; i++, printf("\n"))
		for(j=0; j < cols; j++)
			printf("%d ", img[i*cols + j]);
	printf("----------------------------\n");
}

double* Convolute(double* pixmap, int cols, int rows, double *filter, int dim)
{
	int i, j, k, m;
	int sz = dim>>1; 
	double* data = (double *) malloc(cols * rows * sizeof(double));
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++)
			data[i*cols + j] = 0;
	for(i=sz; i < rows-sz; i++) {
		for(j=sz; j < cols-sz; j++) {
			double sum = 0;
			for(k=0; k<dim; k++) {
				for(m=0; m<dim; m++) {
					sum += pixmap[(i+k-sz)*cols+(j+m-sz)] * filter[k*dim+m];
				}
			}
			data[i*cols + j] = sum;
		}
	}
	return data;
}

double* Convolute_n(double *img, int cols, int rows, double *filter, int dim, int ntime)
{
	int i;
	double *conv = img;
	for (i=0; i < ntime; i++) {
		conv = Convolute(conv, cols, rows, filter, dim);
	}
	return conv;
}

int* normalize(double *img, int cols, int rows) {
	int i, j;
	double mn=1e+10, mx=-1e+10;
	int* data = (int *) malloc(cols * rows * sizeof(int));
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++)
			if (img[i*cols + j] < mn)
				mn = img[i*cols + j];
			else if (img[i*cols + j] > mx)
				mx = img[i*cols + j];
	double rng = mx - mn;
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++) {
			data[i*cols + j] = (int)(255 * (img[i*cols + j] - mn) / rng);
		}
	return data;
}

void normalize_and_print(double* img, int dimx, int dimy, int type, char* name) {

	int *norm_img = normalize(img, dimx, dimy);
	writePixmap(norm_img, dimx, dimy, type, name);
	free(norm_img);
}

void compute_Ix2_Iy2_IxIy(double* Ix, double* Iy, int cols, int rows, 
	double** Ix2, double** Iy2, double** IxIy) {

	int size = cols * rows;
	*Ix2 = (double *) malloc(size * sizeof(double));
	*Iy2 = (double *) malloc(size * sizeof(double));
	*IxIy = (double *) malloc(size * sizeof(double));

	int i;
	for (i=0; i < size; i++) {
		*Ix2[i] = Ix[i] * Ix[i];
		*Iy2[i] = Iy[i] * Iy[i];
		*IxIy[i] = Ix[i] * Iy[i];
	}
}

int main()
{
	int *image, dimx = 513, dimy = 511, type = 5;

	image = readPixmap("boat.pgm", &type, &dimx, &dimy);

	double *image_dbl = int_to_double(image, dimx * dimy);
	double *flt_img = Convolute_n(image_dbl, dimx, dimy, B22, 3, 2);
	// int *flt_img = Convolute_n(image, dimx, dimy, G55, 5, 3);
	// writePixmap(flt_img, dimx, dimy, type, "boat-filtered.pgm");

	double *Ix = Convolute_n(flt_img, dimx, dimy, sobel_h1, 3, 1);
	double *Iy = Convolute_n(flt_img, dimx, dimy, sobel_h2, 3, 1);

	double *Ix2, *Iy2, *IxIy;

	compute_Ix2_Iy2_IxIy(Ix, Iy, dimx, dimy, &Ix2, &Iy2, &IxIy);


	double *Ix2_s = Convolute_n(Ix2, dimx, dimy, B22, 3, 1);
	double *Iy2_s = Convolute_n(Iy2, dimx, dimy, B22, 3, 1);
	double *IxIy_s = Convolute_n(IxIy, dimx, dimy, B22, 3, 1);

	free (image);
	free (flt_img);
	free (image_dbl);
	free (Ix);
	free (Iy);
	free (Ix2);
	free (Iy2);
	free (IxIy);
	free (Ix2_s);
	free (Iy2_s);
	free (IxIy_s);

	return 0;
}
