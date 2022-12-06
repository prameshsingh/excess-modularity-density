#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

typedef struct communities {
	int size;
	int *members;
} comm;


double com_det(int n, double **A, int *par);
double smprng(void);
