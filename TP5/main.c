#include "pixmapIO.h"
#include "kmeans.h"


Sample* Convert_image_to_color_feature(int* img, int size) 
{
	int i, j;
	Sample *list = malloc(sizeof(Sample) * size);
	for (i=0; i < size; i++) {
		Sample *smpl = &list[i];
		smpl->f = malloc(sizeof(double) * 3);
		for (j=0; j<3; j++)
			smpl->f[j] = img[3*i+j] / 255.0;
	}
	return list;
}

Sample* Convert_image_to_color_and_xy_feature(int* img, int dimx, int dimy) 
{
	int i, j;
	int size = dimx * dimy;
	Sample *list = malloc(sizeof(Sample) * size);
	for (i=0; i < size; i++) {
		Sample *smpl = &list[i];
		smpl->f = malloc(sizeof(double) * 5);
		for (j=0; j<3; j++)
			smpl->f[j] = img[3*i+j] / 255.0;
		smpl->f[3] = 1.0*(i/dimx)/dimy;
		smpl->f[4] = 1.0*(i%dimx)/dimx;
	}
	return list;
}

int* build_image(const kmeans_result* res, int size)
{
	int i, j;
	int* out = (int *) malloc(3 * size * sizeof(int));
	for (i=0; i<size; i++) {
		Sample* smpl = &(res->C[res->k[i]]);
		for (j=0; j<3; j++)
			out[3*i+j] = (int)smpl->f[j];
	}
	return out;
}


void normalize_centers_colors(const kmeans_result* res, int K)
{
	int i, j;
	for (i=0; i<K; i++) {
		for (j=0; j<3; j++) {
			res->C[i].f[j] *= 255;
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc < 3) {
		printf("usage: ./main <image-path> <#clusters>\n");
		exit(1);
	}

	int *image, dimx, dimy, type;
	char name[64];

	char *img_path = argv[1];
	int K = atoi(argv[2]);

	// Read the original image
	sprintf(name, "%s.ppm", img_path);
	printf("[Read the original image %s]\n", name);
	image = readPixmap(name, &type, &dimx, &dimy);

	printf("type: %d\n", type);
	printf("dimx: %d\n", dimx);
	printf("dimy: %d\n", dimy);


	// Cluster based on colors only
	printf("[Convert the image to list of samples]\n");
	Sample* list_clr =  Convert_image_to_color_feature(image, dimx * dimy);

	kmeans_result* res_cl = cluster(list_clr, dimx*dimy, 3, K);

	printf("[Normalize clusters centers]\n");
	normalize_centers_colors(res_cl, K);

	printf("[Build the new image]\n");
	int* c_img = build_image(res_cl, dimx * dimy);

	printf("[Write the new image]\n");
	sprintf(name, "%s-clustered-%d-color.ppm", img_path, K);
	writePixmap(c_img, dimx, dimy, type, name);


	// Cluster based on colors and XY
	printf("[Convert the image to list of samples]\n");
	Sample* list_clr_xy =  Convert_image_to_color_and_xy_feature(image, dimx, dimy);

	kmeans_result* res_cl_xy = cluster(list_clr_xy, dimx*dimy, 5, K);

	printf("[Normalize clusters centers]\n");
	normalize_centers_colors(res_cl_xy, K);

	printf("[Build the new image]\n");
	int* c_img_xy = build_image(res_cl_xy, dimx * dimy);

	printf("[Write the new image]\n");
	sprintf(name, "%s-clustered-%d-color-xy.ppm", img_path, K);
	writePixmap(c_img_xy, dimx, dimy, type, name);

	return 0;
}
