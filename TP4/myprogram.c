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

const int nTopCorner = 100;

typedef struct _pixel {
	int i, j;
	double val;
} pixel;

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

double* convolute(double* pixmap, int cols, int rows, double *filter, int dim)
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

double* convolute_n(double *img, int cols, int rows, double *filter, int dim, int ntime)
{
	int i;
	double *conv = img;
	for (i=0; i < ntime; i++) {
		conv = convolute(conv, cols, rows, filter, dim);
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
	}
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

int cmpfunc_pixel (const void * a, const void * b)
{
	return ((pixel*)b)->val - ((pixel*)a)->val;
}

pixel* get_top_n_local_maxima(double* img, int cols, int rows, int n) {
	pixel* maxi = (pixel *) malloc(cols * rows * sizeof(pixel));
	int i, j, c=0;
	for(i=1; i < rows-1; i++)
		for(j=1; j < cols-1; j++) {
			if (is_local_maxima(img, cols, i, j)) {
				pixel p = {i, j, img[i*cols+j]};
				maxi[c++] = p;
			}
	}
	qsort(maxi, c, sizeof(pixel), cmpfunc_pixel);
	printf("Detected corners %d\n", c);
	pixel* px = (pixel *) malloc(n * sizeof(pixel));
	for (i=0; i < n; i++)
		px[i] = maxi[i];
	free (maxi);
	return px;
}

int* add_red_pixel(double* hrs_img, double* flt_img, int cols, int rows, pixel* px, int n) {
	int size = cols * rows;
	int* res = (int *) malloc(3 * size * sizeof(int));
	int i;
	for (i=0; i < size; i++) {
		res[3*i] = res[3*i+1] = res[3*i+2] = (int)flt_img[i];
	}
	for (i=0; i < n; i++) {
		int ind = 3 * (px[i].i * cols + px[i].j);
		res[ind] = 255;
		res[ind+1] = 0;
		res[ind+2] = 0;
	}
	return res;
}

#define tolerance 0.1e-20

void eigenvalues(double A, double B, double C, double D, double* lambda1, double* lambda2) {
	if (B*C <= tolerance)  {
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


/*
 * http://en.wikipedia.org/wiki/Corner_detection#The_Harris_.26_Stephens_.2F_Plessey_.2F_Shi.E2.80.93Tomasi_corner_detection_algorithm
 */
double* Harris(int cols, int rows, double* A, double* B, double* C, double a) {
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

/*
 * http://www.aishack.in/tutorials/the-shitomasi-corner-detector/
 */
double* ShiTomasi(int cols, int rows, double* A, double* B, double* C) {
	int size = cols * rows, i;
	double lambda1, lambda2;
	double* eign = (double *) malloc(size * sizeof(double));
	for (i=0; i < size; i++) {
		eigenvalues(A[i], C[i], C[i], B[i], &lambda1, &lambda2);
		eign[i] = MIN(lambda1, lambda2);
	}
	return eign;
}


int main()
{
	int *image, dimx = 513, dimy = 511, type = 5;

	image = readPixmap("boat.pgm", &type, &dimx, &dimy);

	double *image_dbl = int_to_double(image, dimx * dimy);
	double *flt_img = convolute_n(image_dbl, dimx, dimy, B22, 3, 2);
	// int *flt_img = convolute_n(image_dbl, dimx, dimy, G55, 5, 3);
	// writePixmap(flt_img, dimx, dimy, type, "boat-filtered.pgm");

	double *Ix = convolute_n(flt_img, dimx, dimy, sobel_h1, 3, 1);
	double *Iy = convolute_n(flt_img, dimx, dimy, sobel_h2, 3, 1);

	normalize_and_print(Ix, dimx, dimy, type, "boat-grad-Ix.pgm");
	normalize_and_print(Iy, dimx, dimy, type, "boat-grad-Iy.pgm");

	double *Ix2, *Iy2, *IxIy;

	compute_Ix2_Iy2_IxIy(Ix, Iy, dimx, dimy, &Ix2, &Iy2, &IxIy);

	double *Ix2_s = convolute_n(Ix2, dimx, dimy, B22, 3, 2);
	double *Iy2_s = convolute_n(Iy2, dimx, dimy, B22, 3, 2);
	double *IxIy_s = convolute_n(IxIy, dimx, dimy, B22, 3, 2);

	normalize_and_print(Ix2_s, dimx, dimy, type, "boat-grad-Ix2_s.pgm");
	normalize_and_print(Iy2_s, dimx, dimy, type, "boat-grad-Iy2_s.pgm");
	normalize_and_print(IxIy_s, dimx, dimy, type, "boat-grad-IxIy_s.pgm");

	double a;
	char name[32];
	for (a = 0.00; a <= 0.2; a += 0.02) {
		double* hrs_img = Harris(dimx, dimy, Ix2_s, Iy2_s, IxIy_s, a);

		sprintf(name, "boat-harris-%.2f.pgm", a);
		// normalize_and_print(img, dimx, dimy, type, name);

		pixel* px = get_top_n_local_maxima(hrs_img, dimx, dimy, nTopCorner);
		int* norm_img = add_red_pixel(hrs_img, flt_img, dimx, dimy, px, nTopCorner);

		sprintf(name, "boat-harris-norm-%.2f.pgm", a);
		writePixmap(norm_img, dimx, dimy, 6, name);
		
		free (px);
		free (hrs_img);
		free (norm_img);
	}

	double* shi_tomasi = ShiTomasi(dimx, dimy, Ix2_s, Iy2_s, IxIy_s);
	pixel* px = get_top_n_local_maxima(shi_tomasi, dimx, dimy, nTopCorner);
	int* shi_tomasi_org = add_red_pixel(shi_tomasi, flt_img, dimx, dimy, px, nTopCorner);
	writePixmap(shi_tomasi_org, dimx, dimy, 6, "boat-shi-tomasi-100.pgm");

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
	free (shi_tomasi);
	free (shi_tomasi_org);

	return 0;
}
