

typedef struct _feature {
	double *x;
} Feature;

void cluster(Feature* samples, int N, int fnb, int K);

double** get_clusters();

int* get_clustered_samples();

