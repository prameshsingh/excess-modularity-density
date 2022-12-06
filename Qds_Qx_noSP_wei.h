#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define eps 1E-6
#define alpha 1E-4
#define beta 1E-3
#define toler 1E-10

typedef struct communities {
	int size;
	int *members;
} comm;


double com_det(int n, double **A, int *par);
void init(int n, double **A);
double bisection(int *par, int n, int *ng, double **A);
int bisec_rand(int nc, int *s);
int bisec_eigen(int nc, comm g, int *s, double **A);
double bisec_powmtd(double **mat, double *v, int n);
double bisec_kl_tuning(int this, double **A, comm *G, int ng, int *s);
double final_tuning(int *par, int n, int *ng, double **A);
double agglomerative_tuning(int *par, int n, int *ng, double **A);
int store_state_ft(int *best1, double *best2, int ng, int n, int *Ec_in, double *weiEc_in, double *Ec, int *nc, int **Ex_c, double **weiEx_c, int *par);
double dQds_ft_moveone(int tryx, int xg1, int xg2, int ng, int n, int *Ec_in, double *weiEc_in, int *nc, double *Ec, int *Exc, double *weiExc);
double Mod_dens(double **A, comm *G, int ng);
int store_state(int *state1, double *state, int Ec1_in, int Ec2_in, double weiEc1_in, double weiEc2_in, int nc1, int nc2, double Ec1, double Ec2, int *Ex_c1, int *Ex_c2, double *weiEx_c1, double *weiEx_c2, int *s);
double part_Q_ds_kl_moveone(int Ec1_in, int Ec2_in, double weiEc1_in, double weiEc2_in, int nc1, int nc2, double Ec1, double Ec2, int tryx, int xg, int x1, int x2, double weix1, double weix2);
double term1and2(int Ec_in, double weiEc_in, int nc, double Ec);
int par_to_G(int n, int *par, comm *G, int ng);
int G_to_par(int n, int *par, comm *G, int *ng);
double abmax(double *v,int n);
double smprng(void);

static double m, *deg, **B, ds_one=0, ds_aver;