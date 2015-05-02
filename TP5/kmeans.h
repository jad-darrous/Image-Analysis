#ifndef _KMEANS_
#define _KMEANS_


#define MAX_ITR 100


typedef struct _sample {
	double *f; // list of features
} Sample;


typedef struct _kmeans_result {
	Sample* C; // Clusters centers
	int* k; // Cluster of sample i
} kmeans_result;

/*
 * samples: the samples to cluster
 * N: number of samples
 * fnb: number of features in each sample
 * K: number of clusters
 */
kmeans_result* cluster(Sample* samples, int N, int fnb, int K);


#endif