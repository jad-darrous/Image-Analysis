#include "pixmapIO.h"

#include <math.h>
#include <memory.h>


double B22[9] = {1/16.0,2/16.0,1/16.0, 2/16.0,4/16.0,2/16.0, 1/16.0,2/16.0,1/16.0};
double sobel_h1[9] = {-1.0/4.0, 0, 1.0/4.0,   -2.0/4.0, 0, 2.0/4.0,  -1.0/4.0, 0, 1.0/4.0 };
double sobel_h2[9] = {-1.0/4.0, -2.0/4.0, -1.0/4.0,   0, 0, 0,  1.0/4.0, 2.0/4.0, 1.0/4.0 };


int* Convolute(int* pixmap, int cols, int rows, double *filter, int dim)
{
	int i, j, k, m;
	int sz = dim>>1; 
	int* data = (int *) malloc(cols * rows * sizeof(int));
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
			data[i*cols + j] = (int)sum;
		}
	}
	return data;
}

int* Convolute_n(int *img, int cols, int rows, double *filter, int dim, int ntime)
{
	int i;
	int *conv = img;
	for (i=0; i < ntime; i++) {
		conv = Convolute(conv, cols, rows, filter, dim);
	}
	return conv;
}

int* normalize(int *img, int cols, int rows) {
	int mn=0, mx=0, rng, i, j;
	int* data = (int *) malloc(cols * rows * sizeof(int));
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++)
			if (img[i*cols + j] < mn)
				mn = img[i*cols + j];
			else if (img[i*cols + j] > mx)
				mx = img[i*cols + j];
	rng = mx - mn;
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++) {
			data[i*cols + j] = (int)(255 * (img[i*cols + j] - mn) / rng);
		}
	return data;
}

int* gradientXY(int *Ix, int *Iy, int cols, int rows) {
	int i, j;
	int* data = (int *) malloc(cols * rows * sizeof(int));
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++)
			data[i*cols + j] = (int) sqrt(pow(Ix[i*cols + j], 2) + pow(Iy[i*cols + j], 2));
	return data;
}

#define PI 3.14159265

int* gradientDirection(int *Ix, int *Iy, int cols, int rows) {
	int i, j;
	int* data = (int *) malloc(cols * rows * sizeof(int));
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++) {
			double result = atan2(Iy[i*cols + j], Ix[i*cols + j]) * 180 / PI;
			if (result < 0) result += 180;
  			data[i*cols + j] = ((int)(result) / 45) % 4;
		}
	return data;
}

void non_maximum_suppression(int* img, int* grad_dir, int cols, int rows) {
	int i, j, is_maxima;
	for(i=1; i < rows-1; i++)
		for(j=1; j < cols-1; j++) {
			int val = img[i*cols + j];
			if (grad_dir[i*cols + j] == 0) {
				is_maxima = val > img[i*cols + j-1] && val > img[i*cols + j+1];
			} else if (grad_dir[i*cols + j] == 1) {
				is_maxima = val > img[(i-1)*cols + j+1] && val > img[(i+1)*cols + j-1];
			} else if (grad_dir[i*cols + j] == 2) {
				is_maxima = val > img[(i-1)*cols + j] && val > img[(i+1)*cols + j];
			} else if (grad_dir[i*cols + j] == 3) {
				is_maxima = val > img[(i-1)*cols + j-1] && val > img[(i+1)*cols + j+1];
			}
			if (!is_maxima) img[i*cols + j] = 0;
  			// printf("%d ", data[i*cols + j]);
		}
}

void thresholding(int *data, int cols, int rows, int threshold) {
	int i, j;
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++)
			data[i*cols + j] = data[i*cols + j] > threshold ? data[i*cols + j] : 0;
}

int is_surrounded_by_black(int* img, int cols, int i, int j) {
	int k, l;
	for (k = -1; k <= 1; ++k)
		for (l = -1; l <= 1; ++l)
			if (img[(i+k)*cols + j+l] == 0) {
				return 1;
			}
	return 0;
}

void hysteresis_thresholding(int* img, int cols, int rows, int t_low, int t_high) {
	int i, j;
	thresholding(img, cols, rows, t_low);
	for(i=1; i < rows-1; i++)
		for(j=1; j < cols-1; j++) {
			if (img[i*cols + j] > t_high)
				continue;
			if (is_surrounded_by_black(img, cols, i, j))
				img[i*cols + j] = 0;
		}
}

void try_multiple_thresholding(int *grd, int cols, int rows) {
	int i;
	int* th = (int *) malloc(cols * rows * sizeof(int));
	for (i=0; i<7; i++) {
		memcpy(th, grd, cols * rows * sizeof(int));
		thresholding(th, cols, rows, 10 + i*10);
		char name[100];
		sprintf(name, "boat-threshold-%d.pgm", 10 + i*10);
		writePixmap(th, cols, rows, 5, name);
	}
}

void normalize_and_print(int* img, int dimx, int dimy, int type, char* name) {

	int *norm_img = normalize(img, dimx, dimy);
	writePixmap(norm_img, dimx, dimy, type, name);
	free(norm_img);
}

int main()
{
	int *image, dimx = 513, dimy = 511, type = 5;

	image = readPixmap("boat.pgm", &type, &dimx, &dimy);

	int *flt_img = Convolute_n(image, dimx, dimy, B22, 3, 2);
	writePixmap(flt_img, dimx, dimy, type, "boat-filtered.pgm");

	int *Ix = Convolute_n(flt_img, dimx, dimy, sobel_h1, 3, 1);
	int *Iy = Convolute_n(flt_img, dimx, dimy, sobel_h2, 3, 1);

	int* grd = gradientXY(Ix, Iy, dimx, dimy);

	normalize_and_print(grd, dimx, dimy, type, "boat-grad.pgm");

	// try_multiple_thresholding(grd, dimx, dimy);

	int* theta = gradientDirection(Ix, Iy, dimx, dimy);
	non_maximum_suppression(grd, theta, dimx, dimy);

	normalize_and_print(grd, dimx, dimy, type, "boat-grad-supp.pgm");

	hysteresis_thresholding(grd, dimx, dimy, 30, 70);

	normalize_and_print(grd, dimx, dimy, type, "boat-grad-hysteresis.pgm");


	free (image);
	free (flt_img);
	free (Ix);
	free (Iy);
	free (grd);
	free (theta);

	return 0;
}
