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

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

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
		(*Ix2)[i] = Ix[i] * Ix[i];
		(*Iy2)[i] = Iy[i] * Iy[i];
		(*IxIy)[i] = Ix[i] * Iy[i];
		// printf("%d\n", i);
	}
}


double* harris(int cols, int rows, double* A, double* B, double* C, double a) {
	int size = cols * rows;
	double* img = (double *) malloc(size * sizeof(double));
	int i;
	for (i=0; i < size; i++) {
		double detM = A[i]*B[i] - C[i]*C[i];
		double traceM = A[i]+B[i]; 
		double H = detM - a * traceM * traceM;
		img[i] = H;
	}
	return img;
}

int is_local_maxima(double* img, int cols, int i, int j) {
	int k, l;
	double val = img[i * cols + j];
	for (k = -1; k <= 1; ++k)
		for (l = -1; l <= 1; ++l)
			if ((k || l) && img[(i+k)*cols + j+l] > val) {
				return 0;
			}
	return 1;
}

int cmpfunc (const void * a, const void * b)
{
   return *(double*)b - *(double*)a;
}

double getThreshold(double* img, int cols, int rows, int n) {
	double* maxi = (double *) malloc(cols * rows * sizeof(double));
	int i, j, c=0;
	for(i=1; i < rows-1; i++)
		for(j=1; j < cols-1; j++) {
			if (is_local_maxima(img, cols, i, j))
				maxi[c++] = img[i*cols+j];
	}
	qsort(maxi, c, sizeof(double), cmpfunc);
	int val = maxi[n];
	free (maxi);
	return val;
}


int* add_red(double* ff, double* img, int cols, int rows, double threshold) {
	int size = cols * rows;
	int* res = (int *) malloc(3 * size * sizeof(int));
	int i;
	for (i=0; i < size; i++) {
		if (ff[i] > threshold){
			res[3*i] = 255;
			res[3*i+1] = 0;
			res[3*i+2] = 0;
		} else {
			res[3*i] = res[3*i+1] = res[3*i+2] = (int)img[i];
		}
	}
	return res;
}

int* norm_harris_x(double* img, double* filter, int cols, int rows, int n) {
	double d = getThreshold(img, cols, rows, n);
	return add_red(img, filter, cols, rows, d);
}
int* norm_harris(double* img, int cols, int rows) {
/*	double d = getThreshold(img, cols, rows, 20);
	return add_red(img, cols, rows, d);*/
	int size = cols * rows;
	int* res = (int *) malloc(size * sizeof(int));
	double* maxi = (double *) malloc(size * sizeof(double));
	int i, j, c=0;
	for (i=0; i < size; i++) {
		if (img[i] > 0.3)
			res[i] = 255;
		else if (img[i] < -0.3)
			res[i] = 125;
		else
			res[i] = 0;
		// res[i] = (int)img[i];
	}
	for (i=0; i < size; i++) res[i] = 0;
	for(i=1; i < rows-1; i++)
		for(j=1; j < cols-1; j++) {
			// img[i*cols+j] = eign[i*cols+j] > 50 ? 255 : 0;
			// res[i*cols+j] = is_local_maxima(img, cols, i, j) ? 255 : 0;
			if (is_local_maxima(img, cols, i, j))
				maxi[c++] = img[i*cols+j];
			// printf("%0.2f\t", img[i]);
	}
	int n = 20;
	qsort(maxi, c, sizeof(double), cmpfunc);
	for (i=0; i < size; i++) 
		res[i] = img[i] > maxi[n] ? 255 : 0;
	return res;
}

#define tolerance 0.1e-20

void eigenvalues(double A, double B, double C, double D, double* lambda1, double* lambda2) {
	if(B*C <= tolerance)  {
		*lambda1 = A; // *v1x  =  1;  *v1y  =  0;
		*lambda2 = D; // *v2x  =  0;  *v2y  =  1;
		return;
	}
	double tr =  A  +  D;
	double det =  A  *  D  -  B  *  C;
	double S = sqrt(pow(tr/2, 2) - det);
	*lambda1 = tr/2 + S;
	*lambda2 = tr/2 - S;
}

int* ShiTomasi(int cols, int rows, double* A, double* B, double* C) {
	int size = cols * rows;
	double lambda1, lambda2;
	double* eign = (double *) malloc(size * sizeof(double));
	int i, j;
	for (i=0; i < size; i++) {
		eigenvalues(A[i], C[i], C[i], B[i], &lambda1, &lambda2);
		eign[i] = MIN(lambda1, lambda2);
	}
	int* img = (int *) malloc(size * sizeof(int));
	for (i=0; i < size; i++) img[i] = 0;
	for(i=1; i < rows-1; i++)
		for(j=1; j < cols-1; j++) {
			// img[i*cols+j] = eign[i*cols+j] > 50 ? 255 : 0;
			img[i*cols+j] = is_local_maxima(eign, cols, i, j) ? 255 : 0;
			// printf("%0.2f\t", img[i]);
	}
	return img;
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

	normalize_and_print(Ix, dimx, dimy, type, "boat-grad-Ix.pgm");
	normalize_and_print(Iy, dimx, dimy, type, "boat-grad-Iy.pgm");

	double *Ix2, *Iy2, *IxIy;

	compute_Ix2_Iy2_IxIy(Ix, Iy, dimx, dimy, &Ix2, &Iy2, &IxIy);

	double *Ix2_s = Convolute_n(Ix2, dimx, dimy, B22, 3, 2);
	double *Iy2_s = Convolute_n(Iy2, dimx, dimy, B22, 3, 2);
	double *IxIy_s = Convolute_n(IxIy, dimx, dimy, B22, 3, 2);

	double a;
	for (a = 0.01; a < 0.2; a += 0.02) {
		double* img = harris(dimx, dimy, Ix2_s, Iy2_s, IxIy_s, a);
		char name[32];
		sprintf(name, "boat-harris-%.2f.pgm", a);
		// normalize_and_print(img, dimx, dimy, type, name);

		sprintf(name, "boat-harris-norm-%.2f.pgm", a);
		// int* norm_img = norm_harris(img, dimx, dimy);
		int* norm_img = norm_harris_x(img, flt_img, dimx, dimy, 100);
		// writePixmap(norm_img, dimx, dimy, type, name);
		writePixmap(norm_img, dimx, dimy, 6, name);
		
		free (img);
		free (norm_img);
	}

	int* shi2 = ShiTomasi(dimx, dimy, Ix2_s, Iy2_s, IxIy_s);
	writePixmap(shi2, dimx, dimy, type, "boat-shi-tomasi-2.pgm");

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
