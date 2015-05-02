#include "kmeans.h"
#include "time.h"
#include "stdlib.h"
#include "stdio.h"
#include <memory.h>


void print_samples(Sample* F, int N, int fnb, const char* msg) {

	int i, j;
	printf("----- %s -----\n", msg);
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < fnb; ++j)
		{
			printf("%0.2f ", F[i].f[j]);
		}
		printf("\n");
	}
	printf("\n");
}


int find_min_index(double *arr, int K)
{
	double min = arr[0];
	int i, ind = 0;
	for (i = 1; i < K; ++i)
	{
		if (arr[i] < min) {
			min = arr[i];
			ind = i;
		}
	}
	return ind;
}

double distance(Sample p1, Sample p2, int fnb) {
	int i;
	double sum = 0;
	for (i=0; i<fnb; i++)	
	{
		sum += (p1.f[i]-p2.f[i])*(p1.f[i]-p2.f[i]);
	}
	return sum;
}

kmeans_result* cluster(Sample* samples, int N, int fnb, int K)
{
	int i, j, k, itr;
	srand (time(NULL));

	printf("[Clustering of %d samples with %d features into %d cluster]\n", N, fnb, K);

	Sample *clusters = malloc(sizeof(Sample) * K);
	
	// Init clusters centers
	for (j=0; j<K; j++)
	{
		Sample *smpl = &clusters[j];
		smpl->f = malloc(sizeof(double) * fnb);
		for (i = 0; i < fnb; ++i)
		{
			smpl->f[i] = rand()*1.0/RAND_MAX;
		}
	}

	double diff_feat[K];
	int* cp = malloc(sizeof(int) * N);
	int* ct = malloc(sizeof(int) * N);
	
	int* count = malloc(sizeof(int) * K);

	for (j=0; j<N; j++) cp[j] = -1;

	for (itr = 0; itr < MAX_ITR; ++itr)
	{
		// Assign samples to clusters
		for (j=0; j<N; j++)
		{
			for(i=0; i<K; i++)
			{
				diff_feat[i] = distance(samples[j], clusters[i], fnb);
			}
			ct[j] = find_min_index(diff_feat, K);	
		}

		// Update clusters centers
#if 0
		for(i=0; i<K; i++)
		{
			for (k=0; k<fnb; k++) {
				clusters[i].f[k] = 0;
			}
			int count = 0;
			for (j=0; j<N; j++)
			{
				if (ct[j] == i)
				{
					for (k=0; k<fnb; k++) {
						clusters[i].f[k] += samples[j].f[k];	
					}
					count++;
				}
			}
			if (!count)
				continue;
			for (k=0; k<fnb; k++) {
				clusters[i].f[k] /= count;
			}
		}
#else
		memset(count, 0, sizeof(int) * K);
		for(i=0; i<K; i++)
		{
			for (k=0; k<fnb; k++) {
				clusters[i].f[k] = 0;
			}
		}
		for (j=0; j<N; j++)
		{
			for (k=0; k<fnb; k++) {
				clusters[ct[j]].f[k] += samples[j].f[k];	
			}
			count[ct[j]]++;
		}
		for(i=0; i<K; i++)
		{
			if (!count[i])
				continue;
			for (k=0; k<fnb; k++) {
				clusters[i].f[k] /= count[i];
			}
		}
#endif
		// Check for convergence
		if (!memcmp(cp, ct, sizeof(int) * N))
			break;
		memcpy(cp, ct, sizeof(int) * N);
	}

	if (itr == MAX_ITR)
		printf("[Stopped after reaching the MAX_ITR]\n");
	else
		printf("[Clusters center are converged after %d iteration(s)]\n", itr);

	print_samples(clusters, K, fnb, "Clusters - After");

	free(ct);
	free(count);

	kmeans_result* res = malloc(sizeof(kmeans_result));
	res->k = cp;
	res->C = clusters;
	return res;
}