#include "poc.h"

int predivide(double **a, int *compnt, int n);
void search(double **a, int *flag, int m1, int m2, int n);

int main(int argc, char *argv[]) {
	FILE *fp;
	char fn[200];

	int n = atoi(argv[1]);

	double **A;
	A = (double **)malloc(n * sizeof(double *));
	A[0] = (double *)malloc(n * n * sizeof(double));
	for (int i = 1; i < n; i++) {
		A[i] = A[i-1] + n;
	}

	double max = 0;
	fp = fopen(argv[3],"r");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fscanf(fp, "%lf", &A[i][j]);
			if (A[i][j] > max) max = A[i][j];
		}
		A[i][i] = 0;
	}
	fclose(fp);

	printf("max weight is %g\n", max);

	double thres = atof(argv[2]);
	double thrV = max * thres;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (A[i][j] >= thrV) continue;
			A[i][j] = 0;
			A[j][i] = 0;
		}
	}

	int *partition;
	partition = malloc(n * sizeof(int));
	int *compnt;
	compnt = malloc(n * sizeof(int));
	int compnum = predivide(A, compnt, n);
	printf("there are %d components\n", compnum);

	sprintf(fn, "Qx_t=%g_partition_%s", thres, argv[3]);

	if (compnum == 1) {
		com_det(n, A, partition);
		fp = fopen(fn, "w");
		for (int i = 0; i < n; i++) fprintf(fp, "%d ",partition[i]);
		fclose(fp);
		
		free(A[0]);
		free(A);
		return 0;
	}

	int *compsize, *temp;
	compsize = (int *)malloc(sizeof(int) * compnum);
	temp = (int *)malloc(sizeof(int) * compnum);

    for (int i = 0; i < compnum; i++) {
        compsize[i] = 0;
        temp[i] = 0;
    }

    for (int i = 0; i < n; i++) compsize[compnt[i]-1] += 1;

    int **comp;
    comp = (int **)malloc(sizeof(int *) * compnum);
    for (int i = 0; i < compnum; i++) comp[i] = malloc(sizeof(int) * compsize[i]);
    
    for (int i = 0; i < n; i++) {
        //printf("%d ", compnt[i]);
        comp[compnt[i]-1][temp[compnt[i]-1]] = i;
        temp[compnt[i]-1] += 1;
    }
    free(temp);

    int ncomm = 0;
    for (int i = 0; i < compnum; i++) {
    	int size = compsize[i];
    	printf("this component's size is %d\n", size);
    	if (size < 3) {
    		for (int j = 0; j < size; j++) partition[comp[i][j]] = ncomm + 1;
    		ncomm++;
    		continue;
    	}

    	double **a;
    	a = (double **)malloc(size * sizeof(double *));
    	a[0] = (double *)malloc(size * size * sizeof(double));
    	for (int j = 1; j < size; j++) a[j] = a[j - 1] + size;
    	for (int j = 0; j < size; j++) {
    		a[j][j] = 0;
    		for (int k = j + 1; k < size; k++)
    			a[j][k] = a[k][j] = A[comp[i][j]][comp[i][k]];
    	}

    	int *par;
    	par = (int *)malloc(sizeof(int) * size);
    	com_det(size, a, par);
    	
    	int temp_nc = 1;
    	for (int j = 0; j < size; j++) {
    		if (par[j] > temp_nc) temp_nc = par[j];
    		partition[comp[i][j]] = ncomm + par[j];
    	}
    	ncomm += temp_nc;

    	free(par);
    	free(a[0]);
    	free(a);
    }

	fp = fopen(fn, "w");
	for (int i = 0; i < n; i++) 
		fprintf(fp, "%d ",partition[i]);
	fclose(fp);

	free(compsize);
	for (int i = 0; i < compnum; i++) free(comp[i]);
	free(comp);
		
	free(A[0]);
	free(A);
	return 0;
}

int predivide(double **a, int *compnt, int n){
	int k = 1;
	
	for (int i = 0; i < n; i++) compnt[i] = 0;
	
	for (int i = 0; i < n; i++) {
		if (compnt[i] == 0) {
			//printf("%d\n",i+1);
			compnt[i] = k;
			search(a, compnt, i, k, n);
			k++;
			//printf("\n\n\n");
		}
	}
	return k-1;	
}

void search(double **a, int *flag, int m1, int m2, int n){
	for (int i = 0; i < n; i++) {
		if (flag[i] == 0 && a[m1][i] > 0){
			//printf("%d ",i+1);
			flag[i] = m2;
			search(a, flag, i, m2, n);
		}
	}
	return;
}
