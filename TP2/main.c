#include <stdio.h>
#include <stdlib.h>
#include "ImgUtils.h"


typedef struct _kernel {
	int dim;
	double* data;
} kernel;


kernel* create_kernel(int n) {
	kernel* k = malloc(sizeof(kernel));
	k->dim = n;
	k->data = malloc(n*n*sizeof(double));
	return k;
}

int cmpfunc_gray(const void* a, const void* b)
{
	return *(gray*)a - *(gray*)b;
}

kernel* read_binominal_filter(char* fname)
{
	int i, j, dim, val, sum = 0;
	FILE* ifp = fopen(fname,"r");
	if (ifp == NULL) {
		printf("error in opening file %s\n", fname);
		exit(1);
	}
	fscanf(ifp, "%d", &dim);
	kernel* k = create_kernel(dim);
	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			fscanf(ifp, "%d", &val);
			k->data[i*dim+j] = val;
			sum += val;
		}
	}
	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			k->data[i*dim+j] /= sum;
		}
	}
	return k;
}


gs_image* Convolute(gs_image *img, kernel *filter)
{
	int i, j, k, m;
	int rows = img->hdr->rows, cols = img->hdr->cols;
	int k_sz = filter->dim;
	int _sz = k_sz>>1; 

	gs_image *conv = malloc(sizeof(gs_image));
	conv->hdr = img->hdr;
	conv->data = (gray *) malloc(cols * rows * sizeof(gray));

	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++)
			conv->data[i*cols + j] = 0;

	for(i=_sz; i < rows-_sz; i++) {
		for(j=_sz; j < cols-_sz; j++) {
			
			double sum =0;
			for(k=0; k<k_sz; k++){
				for(m=0; m<k_sz; m++){
					sum += img->data[(i+k-_sz)*cols+(j+m-_sz)] * filter->data[k*k_sz+m];
				}
			}
			conv->data[i*cols + j] = (gray)sum;
			// printf("%d ", conv->data[i*cols + j]);
		}
	}
	return conv;
}

gs_image* Convolute_n(gs_image *img, kernel *filter, int ntime)
{
	int i;
	gs_image *conv = img;
	for (i=1; i < ntime; i++) {
		conv = Convolute(conv, filter);
	}
	return conv;
}

gs_image* Convolute_median(gs_image *img, int filter_size)
{
	int i, j, k, m;
	int rows = img->hdr->rows, cols = img->hdr->cols;
	int _sz = filter_size>>1; 
	int sz = filter_size;
	int sz2 = sz * sz;

	gray* buff = malloc(sz2 * sizeof(gray));

	gs_image *conv = malloc(sizeof(gs_image));
	conv->hdr = img->hdr;
	conv->data = (gray *) malloc(cols * rows * sizeof(gray));

	for(i=0; i < rows; i++)
		for(j=0; j < cols; j++)
			conv->data[i*cols + j] = 0;

	for(i=_sz; i < rows-_sz; i++) {
		for(j=_sz; j < cols-_sz; j++) {
			
			int t = 0;
			for(k=0; k<sz; k++){
				for(m=0; m<sz; m++){
					buff[t++] = img->data[(i+k-_sz) * cols+(j+m-_sz)];
				}
			}
   			qsort(buff, sz2, sizeof(gray), cmpfunc_gray);
			conv->data[i*cols + j] = buff[sz2>>1]; 
			// printf("%d ", conv->data[i*cols + j]);
		}
	}
	free(buff);
	return conv;
}

void Filtering()
{
	printf("Read the image\n");
	gs_image* img = read_image("grenoble_noise.pgm");
	print_image_header(img->hdr);

	printf("Read the filters\n");
	kernel* k3 = read_binominal_filter("binominal_3x3.dat");
	kernel* k5 = read_binominal_filter("binominal_5x5.dat");

	printf("Convolute - Binomial filters\n");
	gs_image *conv_bin3 = Convolute(img, k3);
	gs_image *conv_bin5 = Convolute(img, k5);
	gs_image *conv_bin5_10 = Convolute_n(img, k5, 10);

	printf("Write new image\n");
	write_image("grenoble-bin3.pgm", conv_bin3);
	write_image("grenoble-bin5.pgm", conv_bin5);
	write_image("grenoble-bin5_10.pgm", conv_bin5_10);

	free(conv_bin3);
	free(conv_bin5);
	free(conv_bin5_10);

	printf("Convolute - Median filter\n");
	gs_image *conv_median = Convolute_median(img, 3);
	printf("Write new image\n");
	write_image("grenoble-median-3.pgm", conv_median);
	free(conv_median);

	free(img);
	free(k3);
	free(k5);
}

void Histograms()
{

}

int main()
{
	Filtering();
	Histograms();
	return 0;
}
