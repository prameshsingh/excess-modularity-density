/* Two abstract data structures are normally used to represent community structure:
one way is to use a 1-d array, array index represents the node, element value represents 
community #; second way is a 2d array(or linked-list), ng groups would have ng arrays, each array one 
community, and the element value represents the nodes.*/

/* The first data structure is convenient when you are keeping moving nodes from one 
community to another, in other words, it's dynamic-friendly; the second is convinient
when you want to get some statistics or do some calculations based on a particular 
partition, that being said, it's static-friendly. */

#include "Qds_Qx_noSP_wei.h"

double com_det(int n, double **A, int *par) {
	int ng;
	double Q_ds = 0;
	comm oricluster;
	comm *G;

	init(n, A);
	srand((long int) time(NULL));
	ds_aver = 2 * m / ((double)n * (n - 1));

	G = (comm *)malloc(sizeof(comm) * n);
    for (int i = 1; i < n; i++) G[i].size = 0;   
	
	for (int i =0; i < n; i++) par[i] = 1;
    oricluster.size = n;
    oricluster.members = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) oricluster.members[i] = i;
    G[0] = oricluster;
	ng = 1;
    /*
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			printf("%f\t", A[i][j]);
		printf("\n");
	}
    */
	Q_ds = Mod_dens(A, G, ng);
	printf("size is %d, before partition, Q_ds is %f\n", n, Q_ds);

	double dQ_ds;
	do {
		dQ_ds = 0;
	    double dQ_ds_bi = bisection(par, n, &ng, A);
	    dQ_ds += dQ_ds_bi;
	    Q_ds += dQ_ds_bi;
	    printf("dQ_ds_bi is %g\n", dQ_ds_bi);
	    printf("Based on sum, Q_ds is %f\n", Q_ds);
	    par_to_G(n, par, G, ng);
	    printf("According to definition, after bisection, Q_ds is now %f\n", Mod_dens(A, G, ng));
		printf("Community sizes are:\n");
		for (int i = 0; i < ng; i++)
			printf("%d ", G[i].size);
		printf("\n\n");

		G_to_par(n, par, G, &ng);	 
	    double dQ_ds_ft = final_tuning(par, n, &ng, A);
	    dQ_ds += dQ_ds_ft;
	    Q_ds += dQ_ds_ft;
	    printf("dQ_ds_ft is %g\n", dQ_ds_ft);
	    printf("Based on sum, Q_ds is %f\n", Q_ds);
	    par_to_G(n, par, G, ng);
		printf("According to definition, after Final tuning, Q_ds is now %f\n", Mod_dens(A, G, ng));
		printf("Community sizes are:\n");
		for (int i = 0; i < ng; i++)
			printf("%d ", G[i].size);
		printf("\n\n");


		G_to_par(n, par, G, &ng);
		double dQ_ds_ag = agglomerative_tuning(par, n, &ng, A);
		dQ_ds += dQ_ds_ag;
		Q_ds += dQ_ds_ag;
		printf("dQ_ds_ag is %g\n", dQ_ds_ag);
	    printf("Based on sum, Q_ds is %f\n", Q_ds);
	    par_to_G(n, par, G, ng);
	    printf("According to definition, after agglomerative tuning, Q_ds is now %f\n", Mod_dens(A, G, ng));
		printf("Community sizes are:\n");
		for (int i = 0; i < ng; i++)
			printf("%d ", G[i].size);
		printf("\n\n");
		
	} while (dQ_ds > toler);

	for (int i = 0; i < n; i++) {
		if (deg[i] < 0.001) par[i] = 0;
	}

	
	int csize;
	for (int i = 0; i < ng; i++) {
		csize = G[i].size;
		for (int j = 0; j < csize; j++)
			printf("%d ", G[i].members[j]+1);
		printf("\n");
	}
	/*
	for (int i = 0; i < n; i++)
		printf("%d ", par[i]);
	printf("\n");
	*/

    return 0;	
}

void init(int n, double **A) {
	deg = (double *)malloc(sizeof(double) * n);
	B = (double **)malloc(sizeof(double *) * n);

	for (int i = 0; i < n; i++) {
		B[i] = (double *)malloc(sizeof(double) * n);
		A[i][i] = 0;
	}

	m = 0;
	for (int i = 0; i < n; i++){
		deg[i] = 0;
		for (int j = 0; j < n; j++)
			deg[i] += A[i][j];
		m += deg[i];
	}
	m = m / 2;

	if (m == 0) {
		printf("A is 0 Matrix!!");
	}

	int flag = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			B[i][j] = A[i][j] - deg[i] * deg[j] / (2 * m);
			//printf("%.2f ", B[i][j]);
			if (B[i][j] != 0) flag = 1;
		}
		//printf("\n");
	}
	if (flag == 0) {
		printf("B is 0 Matrix!");
	}
}

double bisection(int *par, int n, int *ng, double **A) {
	comm *G;
	G = (comm *)malloc((*ng*2) * sizeof(comm));
	for (int i = 0; i < *ng*2; i++) G[i].size = 0;
	par_to_G(n, par, G, *ng);
	
	/*
	printf("now there are %d groups, their sizes are:\n", *ng);
	for (int i = 0; i < *ng; i++)
		printf("%d ", G[i].size);
	printf("\n\n");
	*/
	double dQ_ds_bi = 0;
	int ng_temp = *ng;

	for (int this = 0; this < ng_temp; this++) {
		int nc = G[this].size;
		//printf("this community size is %d\n", nc);
		if (nc < 2) continue;
		int *s;
		s = (int *)malloc(nc * sizeof(int));
		for (int i = 0; i < nc; i++) 
			s[i] = 1;
		//printf("group to be bisected has size %d\n", nc);
		bisec_rand(nc, s);
		//bisec_eigen(nc, G[this], s, A);

		double ddQ_ds_bi = bisec_kl_tuning(this, A, G, *ng, s);
		//printf("ddQ_ds_bi is %f\n", ddQ_ds_bi);

		if (ddQ_ds_bi > 0) {
			int np = 0, nn;
			for (int i = 0; i < nc; i++) {
				if (s[i] > 0)
					np++;
			}
			nn = nc - np;
			comm c1, c2;
			c1.members = (int *)malloc(np * sizeof(int));
			c2.members = (int *)malloc(nn * sizeof(int));

			np = nn = 0;
			for (int i = 0; i < nc; i++) {
				if (s[i] == 1) {
					c1.members[np] = G[this].members[i];
					np++;
				}
				else {
					c2.members[nn] = G[this].members[i];
					nn++;
				}
			}
			
			/*
			printf("c1 and c2 are:\n");
			for (int j = 0; j < np; j++)
				printf("%d ", c1.members[j]);
			printf("\n");
			for (int j = 0; j < nn; j++)
				printf("%d ", c2.members[j]);
			printf("\n");*/

			free(G[this].members);

			G[this].members = c1.members; //later can try another way
			G[this].size = np;
			G[*ng].members = c2.members;
			G[*ng].size = nn;
			*ng = *ng + 1;
			dQ_ds_bi += ddQ_ds_bi;
		}
		//printf("ddQ_ds_bi of bisection of community %d is %f\n", this, ddQ_ds_bi);
		free(s);
	}

	G_to_par(n, par, G, ng);

	for (int i = 0; i < *ng; i++)
		free(G[i].members);
	free(G);
	return dQ_ds_bi;
}

int bisec_rand(int nc, int *s) {
	double pro = smprng();
	for (int i = 0; i < nc; i++) {
		if (smprng() >= pro) s[i] = 1;
		else s[i] = -1;
	}
	return 0;
}

int bisec_eigen(int nc, comm g, int *s, double **A) {
	double **c;
	c = (double **)malloc(sizeof(double *) * nc);
	for (int i=0; i<nc; i++) 
		c[i] = (double *)malloc(sizeof(double) * nc);
	
	int flag = 0;
	for (int i = 0; i < nc; i++) 
		for (int j = i+1; j < nc; j++) {
			c[i][j] = c[j][i] = B[g.members[i]][g.members[j]];
			if (c[i][j] != 0) flag = 1;
		}

	for (int i = 0; i < nc; i++) {
		double bik = 0;
		for (int j = 0; j < nc; j++) 
			bik = bik+B[g.members[i]][g.members[j]];
		c[i][i] = B[g.members[i]][g.members[i]]-bik;
		if (c[i][i] != 0) flag = 1;
	}

	if (flag == 0) {
		for (int i=0; i<nc; i++) 
			free(c[i]);
		free(c);
		return 0;
	}

	double x0, x, *V;
		V = (double *)malloc(sizeof(double)*nc);
		x0 = bisec_powmtd(c, V, nc);
		if (x0 < 0) {
			for (int i = 0; i < nc; i++) c[i][i]=c[i][i]-x0;
			x = bisec_powmtd(c, V, nc)+x0;
			for (int i = 0; i < nc; i++) c[i][i]=c[i][i]+x0;
		}
		else x = x0; //leading eigenvalue

		for (int i = 0; i < nc; i++) 
			free(c[i]);
		free(c);
		
		int np = nc, nn = 0;	
		for (int i = 0; i < nc; i++) {
			if (V[i] < 0 || (V[i] < alpha && smprng() > 0.5)) {
				nn++;
				s[i] = -1;
			}
		}
		np -= nn;
		free(V);
		
		return 0;
}

double bisec_powmtd(double **mat, double *v, int n) {
	int i,j,N=0;
	double w1=0,w2=10,sum;
	double *u,*ad,*temp;
	
	//printf("you called pwmtd！\n");
	u=(double *)malloc(sizeof(double)*n);
	ad=(double *)malloc(sizeof(double)*n);
	temp=(double *)malloc(sizeof(double)*n);
	
	for (i=0; i<n; i++) 
		ad[i]=mat[i][i];//to preserve the diagonal elements of **a
	
	for (i=0; i<(int)(n/2); i++) {//the initial u and v are random
		u[i]=1;
		v[i]=1;
	}
	for (i=(int)(n/2); i<n; i++) {
		u[i]=2;
		v[i]=2;
	}
	for (i=0; i<n; i++) temp[i]=v[i]-u[i];
	w1=abmax(u,n);
	
	while (fabs(w1-w2)>eps|| fabs(abmax(temp,n))>eps) {//|| fabs(abmax(temp,n))>eps 
		N++;
		if (N>100) {
			for (i=0; i<n; i++) 
				mat[i][i]=mat[i][i]+alpha;
			N=0;
		}
		for (i=0; i<n; i++) {
			sum=0;
			for (j=0; j<n; j++) 
				sum+=mat[i][j]*v[j];
			u[i]=sum;
		}
		w1=w2;
		w2=abmax(u, n);
		while (w2==0) {
			for (i=0; i<n; i++)
				v[i]=i+0.5;//change 
			for (i=0; i<n; i++) {
				sum=0;
				for (j=0; j<n; j++) 
					sum+=mat[i][j]*v[j];
				u[i]=sum;
			}
			w2=abmax(u, n);
		}
		
		//printf("N is %d, w1 is %.2f,w2 is %.2f\n\n",N,w1,w2);
		//printf("v is :");
		//for (i=0; i<n; i++) printf("%.2g ",v[i]);
		//printf("\nu is :");
		
		for (i=0; i<n; i++) {
			u[i]=u[i]/w2;
			//printf("%.2g ",u[i]);
			temp[i]=u[i]-v[i];
			v[i]=u[i];
		}
		
		//printf("\ntemp is:");
		//for (i=0; i<n; i++)  printf("%.2g ",temp[i]);
		//printf("\n\n");
	}
	//printf("N is %d \n\n",N);
	
	for (i=0; i<n; i++) 
		mat[i][i]=ad[i];
	
	//printf("pwmtd is done!!!!\n\n");
	free(u);
	free(ad);
	free(temp);
	return w2;	
}


/*if we directly operate on the two 'members' arrays, it would be quite compicated; any 
operation involving instantly changing the groups should use the...  */

double bisec_kl_tuning(int this, double **A, comm *G, int ng, int *s) {
	double dQ_ds = 0;

	int nc = G[this].size;

	int nc1 = 0;
	int nc2 = 0;
	for (int i = 0; i < nc; i++) {
		if (s[i] > 0)
			nc1++;
	}
	nc2 = nc - nc1;
	comm c1, c2;
	
	c1.size = nc1;
	c2.size = nc2;
	c1.members = (int *)malloc(nc1 * sizeof(int));
	c2.members = (int *)malloc(nc2 * sizeof(int));
	nc1 = nc2 = 0;
	for (int i = 0; i < nc; i++) {
		if(s[i] > 0)
			c1.members[nc1++] = G[this].members[i];
		else
			c2.members[nc2++] = G[this].members[i];
	}

	/*number of links between node x and group c1, c2*/
	int *Ex_c1, *Ex_c2;
	double *weiEx_c1, *weiEx_c2;
	int Ec1_in = 0, Ec2_in = 0, Ec1c2 = 0;
	double weiEc1_in = 0, weiEc2_in = 0, weiEc1c2 = 0;

	Ex_c1 = (int *)malloc(nc * sizeof(int));
	Ex_c2 = (int *)malloc(nc * sizeof(int));
	weiEx_c1 = (double *)malloc(nc * sizeof(double));
	weiEx_c2 = (double *)malloc(nc * sizeof(double));

	for (int i = 0; i < nc; i++) {
		int node = G[this].members[i];
		Ex_c1[i] = 0;
		weiEx_c1[i] = 0;
		for (int j = 0; j < nc1; j++) {
			if (A[node][c1.members[j]] > 0) {
				Ex_c1[i] += 1;
				weiEx_c1[i] += A[node][c1.members[j]];
			}
		}
		//printf("node %d is fine\n", node);
		Ex_c2[i] = 0;
		weiEx_c2[i] = 0;
		for (int j = 0; j < nc2; j++) {
			if (A[node][c2.members[j]] > 0) {
				Ex_c2[i] += 1;
				weiEx_c2[i] += A[node][c2.members[j]];
			}
		}
		if (s[i] > 0) {
			Ec1_in += Ex_c1[i];
			weiEc1_in += weiEx_c1[i];
			Ec1c2 += Ex_c2[i];
			weiEc1c2 += weiEx_c2[i];
		} else {
			Ec2_in += Ex_c2[i];
			weiEc2_in += weiEx_c2[i];
		}
	}
	Ec1_in /= 2;
	Ec2_in /= 2;
	weiEc1_in /= 2.0;
	weiEc2_in /= 2.0;
	int E_this_in = Ec1_in + Ec2_in + Ec1c2;
	double weiE_this_in = weiEc1_in + weiEc2_in + weiEc1c2;

	/*inner links*/
	/*
	edge_count(A, c1, c1, &weiEc1_in, Ec1_in) / 2;
	edge_count(A, c2, c2, &weiEc2_in, Ec2_in) / 2;
	edge_count(A, c1, c2, &weiEc1c2, Ec1c2);
	E_this_in = Ec1_in + Ec2_in + Ec1c2;
	weiE_this_in = weiEc1_in + weiEc2_in + weiEc1c2;*/

	/*sum of degree (2mc + ec)*/
	double Ec1 = 0;
	for (int i = 0; i < nc1; i++)
		Ec1 += deg[c1.members[i]];//sum of degree 
	double Ec2 = 0;
	for (int i = 0; i < nc2; i++)
		Ec2 += deg[c2.members[i]];
	double E_this = Ec1 + Ec2;

	double ori_state, start_state = 0;
	ori_state = term1and2(E_this_in, weiE_this_in, nc, E_this); //part Q_ds for the state before bisection
	if (nc1 != 0) start_state += term1and2(Ec1_in, weiEc1_in, nc1, Ec1); 
	if (nc2 != 0) start_state += term1and2(Ec2_in, weiEc2_in, nc2, Ec2); 
	dQ_ds = start_state - ori_state;

	int *best_states1;
	int state_len1 = 4 + 3 * nc;
	double *best_states;
	int state_len = 4 + 2 * nc;
	best_states1 = (int *)malloc(state_len1 * sizeof(int));
	best_states = (double *)malloc(state_len * sizeof(double));
	
	int *moved; 
	moved = (int *)malloc(nc * sizeof(int));

	double maxstate = start_state;
	do {
		//printf("start state is: nc1 %d, nc2 %d\n\n", nc1, nc2);
		start_state = maxstate;
		for (int i = 0; i < nc; i++) moved[i] = 0;
		store_state(best_states1, best_states, Ec1_in, Ec2_in, weiEc1_in, weiEc2_in, nc1, nc2, Ec1, Ec2, Ex_c1, Ex_c2, weiEx_c1, weiEx_c2, s); //store the start as best state
		int tie_state_num = 1;

		for (int i = 0; i < nc; i++) {
			int tomove = -1;
			double maxmove = 0;
			double part_Q_ds_try;
			int tie_num = 0;
			for (int j = 0; j < nc; j++) {
				if(moved[j] == 0) {
					part_Q_ds_try = part_Q_ds_kl_moveone(Ec1_in, Ec2_in, weiEc1_in, weiEc2_in, nc1, nc2, Ec1, Ec2, G[this].members[j], s[j], Ex_c1[j], Ex_c2[j], weiEx_c1[j], weiEx_c2[j]);
					if (tomove == -1 || part_Q_ds_try > maxmove) {
						maxmove = part_Q_ds_try;
						tie_num = 1;
						tomove = j;
					} else if (part_Q_ds_try == maxmove) {
						double prob_replace = 1.0 / (++tie_num);
						if (smprng() < prob_replace) tomove = j;
					}
					//printf("try is %lf\n", part_Q_ds_try);
				}
			}
			//printf("tie_num is %d, tomove is %d\n", tie_num, tomove);
			//printf("nc1 is %d, nc2 is %d, tomove is %d, maxmove is %g, maxstate is %g\n", nc1, nc2, tomove, maxmove, maxstate);
		
			/*update parameters for calculation after each move*/
			//printf("tomove is %d, Ec1_in %d, Ec2_in %d, nc1 %d, nc2 %d, Ec1 %d, Ec2 %d, Ec1c2 %d, Ex_c1 %d, Ex_c2 %d\n", G[this].members[tomove], Ec1_in, Ec2_in, nc1, nc2, Ec1, Ec2, Ec1c2, Ex_c1[tomove], Ex_c2[tomove]);
			moved[tomove] = 1;

			
			Ec1_in -= s[tomove] * Ex_c1[tomove];
			Ec2_in += s[tomove] * Ex_c2[tomove];
			weiEc1_in -= s[tomove] * weiEx_c1[tomove];
			weiEc2_in += s[tomove] * weiEx_c2[tomove];
			nc1 -= s[tomove]; nc2 += s[tomove];
			Ec1 -= s[tomove] * deg[G[this].members[tomove]];
			Ec2 += s[tomove] * deg[G[this].members[tomove]];
			for (int j = 0; j < nc; j++) {
				double wei = A[G[this].members[tomove]][G[this].members[j]];
				if (wei < 0.001) continue;
				Ex_c1[j] -= s[tomove];
				Ex_c2[j] += s[tomove];
				weiEx_c1[j] -= s[tomove] * wei;					
				weiEx_c2[j] += s[tomove] * wei;
			}
			s[tomove] = -s[tomove];

			/*update the best state*/
			if (maxmove > maxstate) {
				store_state(best_states1, best_states, Ec1_in, Ec2_in, weiEc1_in, weiEc2_in, nc1, nc2, Ec1, Ec2, Ex_c1, Ex_c2, weiEx_c1, weiEx_c2, s);
				tie_state_num = 1;
				maxstate = maxmove;
			} else if (maxmove == maxstate  &&  i != nc - 1) {
				double prob_replace = 1.0 / (++tie_state_num);
				if (smprng() < prob_replace)
					store_state(best_states1, best_states, Ec1_in, Ec2_in, weiEc1_in, weiEc2_in, nc1, nc2, Ec1, Ec2, Ex_c1, Ex_c2, weiEx_c1, weiEx_c2, s);
			}
		}
		if (maxstate - start_state > toler) {
			Ec1_in = best_states1[0];
			Ec2_in = best_states1[1];
			nc1 = best_states1[2];
			nc2 = best_states1[3];

			weiEc1_in = best_states[0];
			weiEc2_in = best_states[1];
			Ec1 = best_states[2];
			Ec2 = best_states[3];
			
			for (int i = 0; i < nc; i++) {
				Ex_c1[i] = best_states1[i+4];
				Ex_c2[i] = best_states1[i+4+nc];
				s[i] = best_states1[i+4+nc+nc];
				weiEx_c1[i] = best_states[i + 4];
				weiEx_c2[i] = best_states[i + nc + 4];
			}
			dQ_ds += maxstate - start_state;
		}
		//printf("tie_state_num is %d\n", tie_state_num);
		//printf("Ec1_in=%d, Ec2_in=%d, nc1=%d, nc2=%d, Ec1=%d, Ec2=%d, Ec1c2=%d\n", Ec1_in, Ec2_in, nc1, nc2, Ec1, Ec2, Ec1c2);
		//printf("total number of nodes to be moved is %d, start with nc1=%d and nc2=%d, incresement of Q_ds is %g\n", nc, nc1, nc2, maxstate - start_state);
	} while (maxstate - start_state > toler);

	free(Ex_c1);
	free(Ex_c2);
	free(weiEx_c1);
	free(weiEx_c2);
	free(moved);
	free(best_states1);
	free(best_states);
	if (c1.size != nc) {
		free(c1.members);
		free(c2.members);
	}

	return dQ_ds;
}

/*not moving to a new community*/

double final_tuning(int *par, int n, int *ng, double **A) {
	if (*ng == 1) return 0;
	comm *G;
	G = (comm *)malloc(*ng * sizeof(comm));
	for (int i = 0; i < *ng; i++) G[i].size = 0;
	par_to_G(n, par, G, *ng);

	int *nc, *Ec_in, **Ex_c;
	double *Ec, *weiEc_in, **weiEx_c;
	nc = (int *)malloc(*ng * sizeof(int));          //size of each group
	Ec_in = (int *)malloc(*ng * sizeof(int));       //# of links inside each group
	weiEc_in = (double *)malloc(*ng * sizeof(double));  //weighted version of internal links
	Ec = (double *)malloc(*ng * sizeof(double));           //sum of degree within each group
	Ex_c = (int **)malloc(n * sizeof(int *));        //links between each node and each group
	weiEx_c = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++) {
		Ex_c[i] = (int *)malloc(*ng * sizeof(int));
		weiEx_c[i] = (double *)malloc(*ng * sizeof(double));
	}

	for (int i = 0; i < *ng; i++) {
		Ec_in[i] = 0;
		weiEc_in[i] = 0;
		Ec[i] = 0;
		int *c = G[i].members;
		for (int j = 0; j < n; j++) {
			Ex_c[j][i] = 0;
			weiEx_c[j][i] = 0;
			for (int k = 0; k < G[i].size; k++) {
				double wei = A[j][c[k]];
				if (wei < 0.001) continue;
				Ex_c[j][i] += 1;
				weiEx_c[j][i] += wei;
			}
			if (par[j] - 1 == i) { //node j belongs to group i
				Ec_in[i] += Ex_c[j][i];
				weiEc_in[i] += weiEx_c[j][i];
			}
		}
		Ec_in[i] /= 2;
		weiEc_in[i] /= 2.0; 
	}


	for (int i = 0; i < *ng; i++) {
		nc[i] = G[i].size;
		Ec[i] = 0;
		int *c = G[i].members;
		for (int j = 0; j < G[i].size; j++) {
			Ec[i] += deg[c[j]];
		}
	}

	for (int i = 0; i < *ng; i++) free(G[i].members);
	free(G);


	double dQds_ft = 0;

	int *moved;
	moved = (int *)malloc(n * sizeof(int));

	int *best_states1;
	double *best_states2;
	int len1 = *ng * (n + 2) + n;
	int len2 = *ng * (n + 2);
	best_states1 = (int *)malloc(len1 * sizeof(int));
	best_states2 = (double *)malloc(len2 * sizeof(double));
	
	double maxstate, sum_move;
	do {
		maxstate = 0;
		sum_move = 0;
		store_state_ft(best_states1, best_states2, *ng, n, Ec_in, weiEc_in, Ec, nc, Ex_c, weiEx_c, par);

        
		for (int i = 0; i < n; i++) moved[i] = 0;
		int tie_state_num = 1;

		for (int i = 0; i < n; i++) {
			int tomove = -1;
			int movein;
			double maxmove;
			double dQds_ft_try;
			int tie_move_num = 0;
			for (int j = 0; j < n; j++) {    //try to move node j to group k
				if (moved[j] == 0) {
					for (int k = 1; k < *ng + 1; k++) {
						if (k == par[j] || nc[k-1] == 0) continue;

						dQds_ft_try = dQds_ft_moveone(j, par[j]-1, k-1, *ng, n, Ec_in, weiEc_in, nc, Ec, Ex_c[j], weiEx_c[j]);  //
						//printf("dQds try is %g\n", dQds_ft_try);
						if (tomove == -1 || dQds_ft_try > maxmove) {
							maxmove = dQds_ft_try;
							tie_move_num = 1;
							tomove = j;
							movein = k - 1;
						} else if (dQds_ft_try == maxmove) {
							double prob_replace = 1.0 / (tie_move_num + 1);
							if (smprng() < prob_replace) {
								tomove = j;
								movein = k - 1;
							}
							tie_move_num++;
						}
					}
				}
			}
			//printf("tomove is %d\n", tomove); 

			moved[tomove] = 1;

			int moveout = par[tomove] - 1;
			nc[moveout] -= 1;
			nc[movein] += 1;
			Ec_in[moveout] -= Ex_c[tomove][moveout];
			Ec_in[movein] += Ex_c[tomove][movein];
			weiEc_in[moveout] -= weiEx_c[tomove][moveout];
			weiEc_in[movein] += weiEx_c[tomove][movein];
			Ec[moveout] -= deg[tomove];
			Ec[movein] += deg[tomove];
			
			for (int j = 0; j < n; j++) {
				double wei = A[j][tomove];
				if (wei < 0.001) continue; 
				Ex_c[j][moveout] -= 1;
				Ex_c[j][movein] += 1;
				weiEx_c[j][moveout] -= wei;
				weiEx_c[j][movein] += wei;
			}

			par[tomove] = movein + 1;
			sum_move += maxmove;
			//printf("maxmove is %g, sum_move is %f\n", maxmove, sum_move);
			
			/*
			comm *G = (comm *)malloc(sizeof(comm) * n);
    		for (int j = 0; j < n; j++) G[j].size = 0; 
    		int ngtemp = 0;
			for (int j = 0; j < n; j++)
				if (par[j] > ngtemp)
					ngtemp = par[j];
			par_to_G(n, par, G, ngtemp);
			double tempQ_ds = Mod_dens(A, G, ngtemp);

			printf("now Q_ds is %f\n", tempQ_ds);*/

			if (sum_move > maxstate) {
				maxstate = sum_move;
				store_state_ft(best_states1, best_states2, *ng, n, Ec_in, weiEc_in, Ec, nc, Ex_c, weiEx_c, par);
				tie_state_num = 1;
			} else if (sum_move == maxstate) {
				double prob_replace = 1.0 / (tie_state_num + 1);
				if (smprng() < prob_replace) store_state_ft(best_states1, best_states2, *ng, n, Ec_in, weiEc_in, Ec, nc, Ex_c, weiEx_c, par);
				tie_state_num++;
			}
			//for (int j = 0; j < *ng; j++) printf("%d ", nc[j]);
			//printf("maxmove is %f, sum_move is %f\n", maxmove, sum_move);
		}
		//printf("maxstate is %f\n", maxstate);

		if (maxstate > toler) {
			for (int i = 0; i < *ng; i++) {
				Ec_in[i] = best_states1[i];
				weiEc_in[i] = best_states2[i];
				nc[i] = best_states1[i + (*ng)];
				Ec[i] = best_states2[i + (*ng)];
			}
			for (int i = 0; i < n; i++) {
				par[i] = best_states1[i + (*ng) * 2];
				for (int j = 0; j < *ng; j++) {
					Ex_c[i][j] = best_states1[(*ng) * 2 + n + (*ng) * i + j];
					weiEx_c[i][j] = best_states2[(*ng) * 2 + j + (*ng) * i];
				}
			}
			dQds_ft += maxstate;
		}
	} while(maxstate > toler);


	for (int i = 0; i < n; i++) par[i] = best_states1[(*ng) * 2 + i];
    /*
	for (int i = 0; i < n; i++)
	    	printf("%d ", par[i]);
	    printf("\n");
	    */

	free(nc);
	free(Ec_in);
	free(weiEc_in);
	free(Ec);
	for (int i = 0; i < n; i++) {
		free(Ex_c[i]);
		free(weiEx_c[i]);
	}
	free(Ex_c);
	free(weiEx_c);
	free(moved);
	free(best_states1);
	free(best_states2);

	return dQds_ft;
}

double agglomerative_tuning(int *par, int n, int *ng, double **A) {
	/*
	comm *G;
	G = (comm *)malloc(*ng * sizeof(comm));
	for (int i = 0; i < *ng; i++) G[i].size = 0;
	par_to_G(n, par, G, *ng);*/

	int *nc, **Ecc;
	double *Ec, **weiEcc;
	nc = (int *)malloc(*ng * sizeof(int));          //size of each group
	Ec = (double *)malloc(*ng * sizeof(double));           //sum of degree within each group
	Ecc = (int **)malloc(*ng * sizeof(int *));
	weiEcc = (double **)malloc(*ng * sizeof(double *));       //# of links between each pair of group
	for (int i = 0; i < *ng; i++) {
		Ecc[i] = (int *)malloc(*ng * sizeof(int));
		weiEcc[i] = (double *)malloc(*ng * sizeof(double));
	}

	for (int i = 0; i < *ng; i++) {
		nc[i] = 0;
		Ec[i] = 0;
		for (int j = 0; j < *ng; j++){
			Ecc[i][j] = 0;
			weiEcc[i][j] = 0;
		}
	}

	for (int i = 0; i < n; i++) {
		int gi = par[i] - 1;
		nc[gi] += 1;
		Ec[gi] += deg[i];
		for (int j = i + 1; j < n; j++) {
			if (A[i][j] > 0.001) {
				int gj = par[j] - 1;
				Ecc[gi][gj] += 1;
				weiEcc[gi][gj] += A[i][j];
			}
		}
	}
	//printf("hellolllllloooldllfslfsld\n");

	for (int i = 0; i < *ng; i++) {
		for (int j = i + 1; j < *ng; j++) {
			Ecc[i][j] += Ecc[j][i];
			Ecc[j][i] = Ecc[i][j];
			weiEcc[i][j] += weiEcc[j][i];
			weiEcc[j][i] = weiEcc[i][j];
		}
	}

	int *current_state, *best_state;
	current_state = (int *)malloc(*ng * sizeof(int));
	best_state = (int *)malloc(*ng * sizeof(int));
	for (int i = 0; i < *ng; i++) {
		current_state[i] = i;
		best_state[i] = i;
	}

	/*apply an array to maintain all the remaining groups' ID, each time a group C is merged,
	move the last element to C's position, in this way, all remaining IDs are in the left part
	of the array.And you decrease the time complexity without having to change the size of 
	Ec_in, Ec, and Ecc*/ 
	int *c_remain;
	c_remain = (int *)malloc(*ng * sizeof(int));
	for (int i = 0; i < *ng; i++) c_remain[i] = i;   
	double sum_dQds = 0;
	double best_dQds = 0;
	for (int ngnow = *ng; ngnow > 1; ngnow--) {
		double max_dQ_ds = -10000;
		int tie_num = 1;
		int try1 = 0, try2 = 0;
		for (int i = 0; i < ngnow; i++) {
			for (int j = i+1; j < ngnow; j++) {
				int ii = c_remain[i];
				int jj = c_remain[j];
				double part_Qds_m = term1and2(Ecc[ii][ii]+Ecc[jj][jj]+Ecc[ii][jj], weiEcc[ii][ii] + weiEcc[jj][jj] + weiEcc[ii][jj], nc[ii]+nc[jj], Ec[ii]+Ec[jj]);
				double part_Qds_s = term1and2(Ecc[ii][ii], weiEcc[ii][ii], nc[ii], Ec[ii]) + term1and2(Ecc[jj][jj], weiEcc[jj][jj], nc[jj], Ec[jj]);

				double dQds_merge = part_Qds_m - part_Qds_s;
				//printf("dQds_merge is %f\n", dQds_merge);
				if (dQds_merge > max_dQ_ds) {
					max_dQ_ds = dQds_merge;
					try1 = i;
					try2 = j;
					tie_num = 1;
				} else if (dQds_merge == max_dQ_ds) {
					double prob_replace = 1.0 / (++tie_num);
					if (smprng() < prob_replace) {
						try1 = i;
						try2 = j;
					}
				}
			}
		}

		int to_merge1 = c_remain[try1];
		int to_merge2 = c_remain[try2];

		Ecc[to_merge1][to_merge1] += Ecc[to_merge2][to_merge2] + Ecc[to_merge1][to_merge2];
		weiEcc[to_merge1][to_merge1] += weiEcc[to_merge2][to_merge2] + weiEcc[to_merge1][to_merge2];
		Ecc[to_merge2][to_merge2] = 0;
		Ec[to_merge1] += Ec[to_merge2];
		Ec[to_merge2] = 0;
		nc[to_merge1] += nc[to_merge2];
		nc[to_merge2] = 0;
		
		for (int i = 0; i < ngnow; i++) {
			int cur = c_remain[i];
			if (cur == to_merge1 || cur == to_merge2) continue;
			Ecc[to_merge1][cur] += Ecc[to_merge2][cur];
			Ecc[cur][to_merge1] = Ecc[to_merge1][cur];
			weiEcc[to_merge1][cur] += weiEcc[to_merge2][cur];
			weiEcc[cur][to_merge1] = weiEcc[to_merge1][cur];
		}

		if (try2 != ngnow-1) c_remain[try2] = c_remain[ngnow-1];
	
		sum_dQds += max_dQ_ds;
		if (sum_dQds >= best_dQds) {
			best_dQds = sum_dQds;
			for (int i = 0; i < *ng; i++) {
				if (current_state[i] == to_merge2)
					current_state[i] = to_merge1;
				best_state[i] = current_state[i];
			}
		}
		//printf("max_dQ_ds is %g, merge1 is %d, merge2 is %d\n", max_dQ_ds, to_merge1, to_merge2);
		//printf("sum_dQds is %f\n", sum_dQds);
	}
	for (int i = 0; i < n; i++) 
		par[i] = best_state[par[i]-1] + 1;

	free(best_state);
	free(current_state);
	for (int i = 0; i < *ng; i++) {
		free(Ecc[i]);
		free(weiEcc[i]);
	}
	free(Ecc);
	free(weiEcc);
	free(nc);
	free(Ec);

	return best_dQds;
}


int store_state_ft(int *best1, double *best2, int ng, int n, int *Ec_in, double *weiEc_in, double *Ec, int *nc, int **Ex_c, double **weiEx_c, int *par) {
	for (int i = 0; i < ng; i++) {
		best1[i] = Ec_in[i];
		best1[i + ng] = nc[i];
		best2[i] = weiEc_in[i];
		best2[i + ng] = Ec[i];
	}

	for (int i = 0; i < n; i++) {
		best1[i + ng * 2] = par[i];
		for (int j = 0; j < ng; j++) {
			best1[i * ng + ng * 2  + n + j] = Ex_c[i][j];
			best2[i * ng + ng * 2 + j] = weiEx_c[i][j];
		}
	}
	return 0;
}

double dQds_ft_moveone(int tryx, int xg1, int xg2, int ng, int n, int *Ec_in, double *weiEc_in, int *nc, double *Ec, int *Exc, double *weiExc) {
	double pre_part_dQds = 0;
	double post_part_dQds = 0;
	
	pre_part_dQds = term1and2(Ec_in[xg1], weiEc_in[xg1], nc[xg1], Ec[xg1]) + term1and2(Ec_in[xg2], weiEc_in[xg2], nc[xg2], Ec[xg2]);
	post_part_dQds = term1and2(Ec_in[xg1]-Exc[xg1], weiEc_in[xg1] - weiExc[xg1], nc[xg1]-1, Ec[xg1]-deg[tryx]) + term1and2(Ec_in[xg2]+Exc[xg2], weiEc_in[xg2] + weiExc[xg2], nc[xg2]+1, Ec[xg2]+deg[tryx]);
	
	return (post_part_dQds - pre_part_dQds);
}

void edge_count(double **A, comm c1, comm c2, double *weiSum, int *n_edge) {
        *weiSum = 0;
        *n_edge = 0;
        if(c2.size == 0) return;
        for (int i = 0; i < c1.size; i++)
                for (int j = 0; j < c2.size; j++) {
                        double wei = A[c1.members[i]][c2.members[j]];
                        if (wei > 0) {
                                *weiSum += wei;
                                *n_edge += 1;
                        }
                }
}

double Mod_dens(double **A, comm *G, int ng) {
	double Q_ds = 0;	
	
	for (int i = 0; i < ng; i++) {
		int nc = G[i].size;
		if (nc == 0) continue;
		int Ec_in;
		double weiEc_in;
		edge_count(A, G[i], G[i], &weiEc_in, &Ec_in);
		Ec_in /= 2;
		weiEc_in /= 2;
		double Ec = 0;
		for (int j = 0; j < nc; j++) 
			Ec += deg[G[i].members[j]];

		double dens;
		if (nc == 1) dens = ds_one;
		else dens = (double) Ec_in  * 2 / (nc * (nc - 1));

		/*Qx needs to substract average density*/
		dens -= ds_aver;

		Q_ds += weiEc_in / m * dens;
		dens *= Ec / (2 * m);
		Q_ds -= dens * dens;

	}

	return Q_ds;
}

int store_state(int *state1, double *state, int Ec1_in, int Ec2_in, double weiEc1_in, double weiEc2_in, int nc1, int nc2, double Ec1, double Ec2, int *Ex_c1, int *Ex_c2, double *weiEx_c1, double *weiEx_c2, int *s) {
	int nc = nc1 + nc2;
	state1[0] = Ec1_in; state1[1] = Ec2_in;
	state1[2] = nc1; state1[3] = nc2;

	state[0] = weiEc1_in; state[1] = weiEc2_in;
	state[2] = Ec1; state[3] = Ec2;

	for (int i = 0; i < nc; i++) {
		state1[i + 4] = Ex_c1[i];
		state1[i + nc + 4] = Ex_c2[i];
		state1[i+ 2 * nc + 4] = s[i];
		state[i + 4] = weiEx_c1[i];
		state[i + nc + 4] = weiEx_c2[i];
	}
	return 0;
}

double part_Q_ds_kl_moveone(int Ec1_in, int Ec2_in, double weiEc1_in, double weiEc2_in, int nc1, int nc2, double Ec1, double Ec2, int tryx, int xg, int x1, int x2, double weix1, double weix2) {
	//int nc = nc1 + nc2;

	Ec1_in -= xg * x1; 
	Ec2_in += xg * x2;
	weiEc1_in -= xg * weix1;
	weiEc2_in += xg * weix2;
	nc1 -= xg; 
	nc2 += xg;
	Ec1 -= xg * deg[tryx];
	Ec2 += xg * deg[tryx];

	double moved = 0;
	if (nc1 != 0) moved += term1and2(Ec1_in, weiEc1_in, nc1, Ec1);
	if (nc2 != 0) moved += term1and2(Ec2_in, weiEc2_in, nc2, Ec2);
	//if (moved < -5) printf("xg is %d x1 is %d x2 is %d\n\n\n", xg, x1, x2);

	return moved;
}

double term1and2(int Ec_in, double weiEc_in, int nc, double Ec) {
	if (nc == 0) return 0;
	double dens;
	if (nc == 1)
		dens = ds_one;
	else
		dens = (double)Ec_in * 2 / (nc * (nc -1));

	/*Qx needs to substract average density*/
	dens -= ds_aver;

	double term1 = dens * weiEc_in / m;

	dens *= Ec / (2 * m);
	double term2 = - dens * dens;

	return (term1 + term2);
}

//time complexity O(n)
int par_to_G(int n, int *par, comm *G, int ng) {
	int cn; //cluster name
	int *sizecount;

	sizecount = (int *)malloc(ng * sizeof(int));
	for (int i = 0; i < ng; i++) {
		G[i].size = 0;
		sizecount[i] = 0;
	}
	for (int i = 0; i < n; i++) 
		G[par[i] - 1].size++;

	for (int i = 0; i < ng; i++) 
		G[i].members = (int *)malloc(G[i].size * sizeof(int));
	
	for (int i = 0; i < n; i++) {
		cn = par[i] - 1;
		G[cn].members[sizecount[cn]] = i;
		sizecount[cn]++;
	}

	free(sizecount);

	return 0;
}

//time complexity O(n), would ignore void community and change ng
int G_to_par(int n, int *par, comm *G, int *ng) {
	//int temp;
	int *t;
	int nc = 0;

	for (int i = 0; i < *ng; i++) {
		if (G[i].size != 0) {
			nc++;
			t = G[i].members;
		}
		for (int j = 0; j < G[i].size; j++)
			par[t[j]] = nc;
	}

	*ng = nc;

	return 0;
}

double abmax(double *v,int n) {
	int i;
	double temp=v[0];
	for(i=1;i<n;i++)
		if(fabs(v[i])>fabs(temp) || (fabs(v[i])==fabs(temp) && v[i]>temp))
			temp=v[i];
	return temp;
}

double smprng(void) {

	return ((double) rand())/((double) RAND_MAX);	
}


