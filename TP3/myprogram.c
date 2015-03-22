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

double* gradient_magnitude(int *Ix, int *Iy, int cols, int rows) {
	int i, j;
	double* data = (double *) malloc(cols * rows * sizeof(double));
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++)
			data[i*cols + j] = hypot(Ix[i*cols + j], Iy[i*cols + j]);
	return data;
}

#define PI 3.14159265

int discretize_angle(double a) {

	// printf("%0.2f\n", a);
	a -= 22.5;
	if (a < 0) a += 180;
	int e = (int) (a + 0.51);
	e = (e - 1 + 180) % 180;
	return ((e / 45) + 1) % 4;
}

int* gradient_direction(int *Ix, int *Iy, int cols, int rows) {
	int i, j;
	int* data = (int *) malloc(cols * rows * sizeof(int));
	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++) {
			double a = atan2(Iy[i*cols + j], Ix[i*cols + j]) * 180 / PI;
			data[i*cols + j] = discretize_angle(a);
		}
	return data;
}

void non_maximum_suppression(double* img, int* grad_dir, int cols, int rows) {
	int i, j, is_maxima;
	for(i=1; i < rows-1; i++)
		for(j=1; j < cols-1; j++) {
			int val = img[i*cols + j];
			if (grad_dir[i*cols + j] == 0) {
				is_maxima = val > img[i*cols + j-1] && val > img[i*cols + j+1];
			} else if (grad_dir[i*cols + j] == 1) {
				is_maxima = val > img[(i-1)*cols + j-1] && val > img[(i+1)*cols + j+1];
			} else if (grad_dir[i*cols + j] == 2) {
				is_maxima = val > img[(i-1)*cols + j] && val > img[(i+1)*cols + j];
			} else if (grad_dir[i*cols + j] == 3) {
				is_maxima = val > img[(i-1)*cols + j+1] && val > img[(i+1)*cols + j-1];
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
			if ((k || l) && img[(i+k)*cols + j+l]) {
				return 0;
			}
	return 1;
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

void try_multiple_thresholding(double *m, int cols, int rows) {
	int i;
	int *grd = double_to_int(m, cols * rows);
	int* th = (int *) malloc(cols * rows * sizeof(int));
	for (i=0; i<7; i++) {
		memcpy(th, grd, cols * rows * sizeof(int));
		thresholding(th, cols, rows, 10 + i*10);
		char name[100];
		sprintf(name, "boat-threshold-%d.pgm", 10 + i*10);
		writePixmap(th, cols, rows, 5, name);
	}
	free(th); free(grd);
}

void normalize_and_print(double* img, int dimx, int dimy, int type, char* name) {

	int *norm_img = normalize(img, dimx, dimy);
	writePixmap(norm_img, dimx, dimy, type, name);
	free(norm_img);
}

void apply_hysteresis_thresholding(double* img, int dimx, int dimy, int type, int t_low, int t_high) {
	int *grd = double_to_int(img, dimx * dimy);
	hysteresis_thresholding(grd, dimx, dimy, t_low, t_high);
	double* m = int_to_double(grd, dimx * dimy);
	char name[50];
	sprintf(name, "boat-hysteresis-%d-%d.pgm", t_low, t_high);
	normalize_and_print(m, dimx, dimy, type, name);
	free(m); free(grd);
}

void test_discretize_angle()
{
	assert (discretize_angle(-.5) == 0);
	assert (discretize_angle(0) == 0);
	assert (discretize_angle(.5) == 0);
	assert (discretize_angle(20) == 0);
	assert (discretize_angle(22) == 0);
	assert (discretize_angle(30) == 1);
	assert (discretize_angle(44) == 1);
	assert (discretize_angle(60) == 1);
	assert (discretize_angle(67) == 1);
	assert (discretize_angle(89.5) == 2);
	assert (discretize_angle(90) == 2);
	assert (discretize_angle(112) == 2);
	assert (discretize_angle(140) == 3);
	assert (discretize_angle(157) == 3);
	assert (discretize_angle(158) == 0);
	assert (discretize_angle(177) == 0);
	assert (discretize_angle(180) == 0);
	printf("test_discretize_angle pass\n");
	exit(1);
}

void test_gradient()
{
	int *image, dimx = 513, dimy = 511, type = 5;
	image = readPixmap("square.pgm", &type, &dimx, &dimy);
	printf("%d %d %d \n", type, dimx, dimy);
	int *Ix = Convolute_n(image, dimx, dimy, sobel_h1, 3, 1);
	int *Iy = Convolute_n(image, dimx, dimy, sobel_h2, 3, 1);
	int* theta = gradient_direction(Ix, Iy, dimx, dimy);
	print_int_arr(theta, dimx, dimy);
	exit(1);
}

int main()
{
	// test_discretize_angle();
	// test_gradient();

	int *image, dimx = 513, dimy = 511, type = 5;

	image = readPixmap("boat.pgm", &type, &dimx, &dimy);

	int *flt_img = Convolute_n(image, dimx, dimy, B22, 3, 2);
	// int *flt_img = Convolute_n(image, dimx, dimy, G55, 5, 3);
	writePixmap(flt_img, dimx, dimy, type, "boat-filtered.pgm");

	int *Ix = Convolute_n(flt_img, dimx, dimy, sobel_h1, 3, 1);
	int *Iy = Convolute_n(flt_img, dimx, dimy, sobel_h2, 3, 1);

	double* m = gradient_magnitude(Ix, Iy, dimx, dimy);

	normalize_and_print(m, dimx, dimy, type, "boat-grad-mag.pgm");

	try_multiple_thresholding(m, dimx, dimy);

	int* theta = gradient_direction(Ix, Iy, dimx, dimy);
	non_maximum_suppression(m, theta, dimx, dimy);

	normalize_and_print(m, dimx, dimy, type, "boat-max-supp.pgm");

	apply_hysteresis_thresholding(m, dimx, dimy, type, 10, 50);
	apply_hysteresis_thresholding(m, dimx, dimy, type, 30, 70);
	apply_hysteresis_thresholding(m, dimx, dimy, type, 50, 90);


	free (image);
	free (flt_img);
	free (Ix);
	free (Iy);
	free (m);
	free (theta);

	return 0;
}
