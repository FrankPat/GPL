/*
              = = = = = = = = = = = = = = = = = = = = = = = = = =
              =             _______     ______                  =
              =            /           /     /    /             =
              =           /   ___     /_____/    /              =
              =          /      /    /          /               =
              =         /______/    /          /______          =
              =                                                 =
              =       G e n e r a l   P o w e r   L a w         =
              =                                                 =
              = = = = = = = = = = = = = = = = = = = = = = = = = =

	Copyright (c) 1998, Frank Pattyn and Wim Van Huele (fpattyn@vub.ac.be)

	You may refer to this software as:

	Pattyn, F. and W. Van Huele (1998) Power Law or Power Flaw?, Earth Surface
	Processes and Landforms 23: 761-767.

	This software is written by:

		Frank PATTYN and Wim VAN HUELE
		Department of Geography
		Vrije Universiteit Brussel
		Pleinlaan 2, B-1050 Brussel, BELGIUM
		e-mail: fpattyn@vub.ac.be

	If you use this software for your scientific work or your publications,
	please don't forget to acknowledge explicitly the use of it.

	Send me an e-mail if:

		-	you want to be informed (by e-mail) about eventual new versions of the
			software.
		-	you have found a bug AND you have CLEARLY idendified the circumstances
			for it to occur.

	Legal matters, disclaimer of warranty, merchantability or fitness for any
	purpose: Permission to use, copy and distribute this software for any purpose
	without fee is hereby granted, provided that both the copyright notice and
	this permission notice appear in all copies, and that the name of GPL is not
	used in advertising or publicity pertaining to distribution of the software
	without specific, written prior permission. The Department of Geography
	(DGGF) and the Vrije Universiteit Brussel (VUB) make no representation about
	the suitability of this software for any purpose. It is provided "as is"
	without expressed or implied warranty. It is provided with no support and
	without obligation on the part of DGGF/VUB to assist in its use, correction,
	modification or enhancement.

	Frank Pattyn and Wim Van Huele
	September 18, 1998

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SWAP(a, b) {temp = (a); (a) = (b); (b) = temp;}
#define SMALL 1.0e-10

void readfiles();
float *mat1d();
int *int1d();
void int1d_init();
float **mat2d();
void mat1d_init();
void mat2d_init();
void mat1d_free();
void int1d_free();
void mat2d_free();
char *char1d();
void power_estimate();
float estimate();
float correlation();
float rmse();
void solve1();
void solve2();
float str_to_f();

float accuracy, orig_x, orig_y;
int variance, origin_calc;

int main(argc, argv)

	int argc;
	char **argv;
{
	float *oldx, *oldy, xmin, xmax, ymin, ymax, width, errorsum;
	float macht, *newx, *newy, *newvar, *covar;
	int max, i, k, l;
	char *inputname, *outputname, *att;
	FILE *infile, *outfile;

	origin_calc = 0;
	variance = 0;
	accuracy = 0.01;
	att = char1d(0, 40);
	printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	printf("*              G P L : G e n e r a l   P o w e r   L a w                  *\n");
	printf("*                                                                         *\n");
	printf("* Copyright (c) 1998, Frank Pattyn and Wim Van Huele (fpattyn@vub.ac.be)  *\n");
	printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	if (argc == 1){
		inputname = char1d(0, 40);
		outputname = char1d(0, 40);
		readfiles(inputname, outputname);
		infile = fopen(inputname, "r");
		outfile = fopen(outputname, "w");
	}
	else{
		if (argc >= 3){
			infile = fopen(argv[1], "r");
			outfile = fopen(argv[2], "w");
		}
		else{
			printf("\n Wrong number of arguments !!\n\n");
			printf("    USAGE:   gpl\n");
			printf("             gpl inputfile outputfile [-Aaccuracy] [-Ox0/y0] [-C]\n\n");
			exit(1);
		}
	}
	for (i = 3; i <= argc - 1; i++){
		if (argv[i][0] == '-'){
			if (argv[i][1] == 'A'){
				k = 2;
				l = 0;
				while (argv[i][k] != '\0'){
					att[l++] = argv[i][k];
					++k;
				}
				att[l] = '\0';
				accuracy = str_to_f(att);
			}
			if (argv[i][1] == 'O'){
				origin_calc = 1;
				k = 2;
				l = 0;
				while (argv[i][k] != '/'){
					att[l++] = argv[i][k];
					++k;
				}
				att[l] = '\0';
				orig_x = str_to_f(att);
				++k;
				l = 0;
				while (argv[i][k] != '\0'){
					att[l++] = argv[i][k];
					++k;
				}
				att[l] = '\0';
				orig_y = str_to_f(att);
			}
			if (argv[i][1] == 'C'){
				variance = 1;
			}
		}
	}
	fscanf(infile, "%d\n", &max);
	oldx = mat1d(1, max);
	oldy = mat1d(1, max);
	covar = mat1d(1, max);
	mat1d_init(oldx, 1, max);
	mat1d_init(oldy, 1, max);

	for (i = 1; i <= max; i++){
		if (variance == 1){
			fscanf(infile, "%f\t%f\t%f\n", &oldx[i], &oldy[i], &covar[i]);
		}
		else{
			fscanf(infile, "%f\t%f\n", &oldx[i], &oldy[i]);
			covar[i] = 1.0;
		}
		if (i == 1){
			xmin = oldx[1];
			xmax = oldx[1];
			ymin = oldy[1];
			ymax = oldy[1];
		}
		else{
			xmin = (oldx[i] < xmin ? oldx[i] : xmin);
			xmax = (oldx[i] > xmax ? oldx[i] : xmax);
			ymin = (oldy[i] < ymin ? oldy[i] : ymin);
			ymax = (oldy[i] > ymax ? oldy[i] : ymax);
		}
	}
	width = (xmax - xmin) / 2.0;
	power_estimate(outfile, oldx, oldy, max, xmin, ymin, width, covar);
	mat1d_free(oldx, 1, max);
	mat1d_free(oldy, 1, max);
	mat1d_free(covar, 1, max);
	fclose(infile);
	fclose(outfile);
	printf("\n\n PRESS ENTER TO FINISH");
	i = getchar();
	i = getchar();
	return 0;
}


/*
===================================================================
Power estimate					y - y0 = a * |x - x0|^b		of the data
===================================================================
*/

void power_estimate(outfile, oldx, oldy, max, xmin, ymin, width, covar)

	float *oldx, *oldy, xmin, ymin, width, *covar;
	int max;
	FILE *outfile;
{
	float rib, ric, rmb, rmc, deltac, flag, flagg, corcoef, rmserr, macht = 2.0;
	float *newy, *coeff, *error;
	float step = 0.1;
	int i, n, count = 0;

	n = (origin_calc == 1 ? 1 : 3);
	coeff = mat1d(1, n);
	error = mat1d (1, max);
	newy = mat1d (1, max);
	mat1d_init(coeff, 1, n);
	mat1d_init(error, 1, max);
	mat1d_init(newy, 1, max);
	flagg = 1;
	coeff[1] = 1.0;
	if (origin_calc != 1){
		coeff[2] = xmin + width;
		coeff[3] = ymin;
	}
	fprintf(outfile, "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
	fprintf(outfile, "*              G P L : G e n e r a l   P o w e r   L a w                  *\n");
	fprintf(outfile, "*                                                                         *\n");
	fprintf(outfile, "* Copyright (c) 1998, Frank Pattyn and Wim Van Huele (fpattyn@vub.ac.be)  *\n");
	fprintf(outfile, "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n\n");
	printf("       Iteration      Power estimate           RMSE\n");
	printf("-----------------------------------------------------------\n");
	do {
		++count;
		estimate(oldx, oldy, max, n, newy, error, coeff, macht, covar);
		rmserr = rmse(oldy, newy, max);
		printf("       %5d          %10.5f      %15.5f\n", count, macht, rmserr);
		++count;
		rib = macht;
		ric = rmse(oldy, newy, max);
		macht += step;
		estimate(oldx, oldy, max, n, newy, error, coeff, macht, covar);
		rmserr = rmse(oldy, newy, max);
		printf("       %5d          %10.5f      %15.5f\n", count, macht, rmserr);
		deltac = rmserr - ric;
		if (deltac > 0){
			step = -step;
		}
		else{
			rmc = ric;
			rmb = rib;
			ric = rmserr;
			rib = macht;
		}
		flag = 1;
		do {
			++count;
			macht += step;
			estimate(oldx, oldy, max, n, newy, error, coeff, macht, covar);
			rmserr = rmse(oldy, newy, max);
			printf("       %5d          %10.5f      %15.5f\n", count, macht, rmserr);
			deltac = rmserr - ric;
			if (deltac > 0){
				flag = 0;
			}
			else{
				rmc = ric;
				rmb = rib;
				ric = rmserr;
				rib = macht;
			}
			if (count > 200) flag = 0;
		} while (flag == 1);
		step = step / 2.0;
		if ((fabs(step) < accuracy) || (count > 200)) flagg = 0;
	} while (flagg == 1);

	macht = rib;
	estimate(oldx, oldy, max, n, newy, error, coeff, macht, covar);
	rmserr = rmse(oldy, newy, max);
	corcoef = correlation(oldy, newy, max);
	++count;
	printf("       %5d          %10.5f      %15.5f\n", count, macht, rmserr);
	if (origin_calc == 1){
		fprintf(outfile, "Power Estimate : y - y0 = a * |x - x0|^b with fixed origin\n");
		fprintf(outfile, "----------------------------------------------------------\n\n");
	}
	else{
		fprintf(outfile, "Power Estimate : y - y0 = a * |x - x0|^b\n");
		fprintf(outfile, "----------------------------------------\n\n");
	}
	fprintf(outfile, "  correlation coefficient  :  %f\n", corcoef);
	fprintf(outfile, "  RMSE                     :  %f\n", rmserr);
	fprintf(outfile, "  coefficient a            :  %e\n", coeff[1]);
	fprintf(outfile, "  coefficient b            :  %f\n", macht);
	if (origin_calc == 1){
		fprintf(outfile, "  coefficient x0            :  %e\n", coeff[2]);
		fprintf(outfile, "  coefficient y0            :  %e\n\n", coeff[3]);
	}
	else{
		fprintf(outfile, "  coefficient x0            :  %e\n", coeff[2]);
		fprintf(outfile, "  coefficient y0            :  %e\n\n", coeff[3]);
	}
	fprintf(outfile, " orig-x\t orig-y\t new-y\t error-y\t log|x-x0|\t log(y-y0)\t log(newy-y0)\n");
	for (i = 1; i <= max; i++){
		if (origin_calc == 1){
			fprintf(outfile,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n", oldx[i], oldy[i], newy[i],
				error[i], log(fabs(oldx[i] - orig_x)), log(fabs(oldy[i] - orig_y)),
				log(fabs(newy[i] - orig_y)));
		}
		else{
			fprintf(outfile,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n", oldx[i], oldy[i], newy[i],
				error[i], log(fabs(oldx[i] - coeff[2])), log(fabs(oldy[i] - coeff[3])),
				log(fabs(newy[i] - coeff[3])));
		}
	}
	fprintf(outfile,"\n");
	mat1d_free(coeff, 1, n);
	mat1d_free(error, 1, max);
	mat1d_free(newy, 1, max);
}


/*
===================================================================
Mathematical solution : y - y0 = a * |x - x0|^b
===================================================================
*/

float estimate(x, y, max, n, newy, error, coeff, macht, covar)

	float *x, *y, *newy, *error, *coeff, macht, *covar;
	int max, n;
{
	float *vect, **em, **a, **aq, **b, *f, **ainv, **q;
	float *new_coef, dx, dummy, covsum;
	float **aat, *ainvf, **ainvb, *col, d, errorsum;
	int i, j, k, *index, teller;

	vect = mat1d(1, n);
	f = mat1d(1, max);
	ainvf = mat1d(1, max);
	col = mat1d(1, max);
	index = int1d(1, max);
	em = mat2d(1, n, 1, n);
	a = mat2d(1, max, 1, (2 * max));
	aq = mat2d(1, max, 1, (2 * max));
	q = mat2d(1, (2 * max), 1, (2 * max));
	b = mat2d(1, max, 1, n);
	ainv = mat2d(1, max, 1, max);
	ainvb = mat2d(1, max, 1, n);
	aat = mat2d(1, max, 1, max);

	mat1d_init(vect, 1, n);
	mat1d_init(f, 1, max);
	mat1d_init(ainvf, 1, max);
	mat1d_init(col, 1, max);
	int1d_init(index, 1, max);
	mat2d_init(em, 1, n, 1, n);
	mat2d_init(a, 1, max, 1, (2 * max));
	mat2d_init(aq, 1, max, 1, (2 * max));
	mat2d_init(q, 1, (2 * max), 1, (2 * max));
	mat2d_init(b, 1, max, 1, n);
	mat2d_init(ainv, 1, max, 1, max);
	mat2d_init(ainvb, 1, max, 1, n);
	mat2d_init(aat, 1, max, 1, max);

	for (i = 1, covsum = 0.0; i <= max; i++){
		covsum += covar[i];
	}
	covsum /= (float) max;
	for (i = 1; i <= max; i++){
		q[2 * i - 1][2 * i - 1] = covar[i] / covsum;
		q[2 * i][2 * i] = covar[i] / covsum;
	}
	for (teller = 1; teller <= 10; teller++){
		for (i = 1; i <= max; i++){
			dx = (origin_calc == 1 ? (x[i] - orig_x) : (x[i] - coeff[2]));
			dx = (fabs(dx) < SMALL ? SMALL : dx);
			if (dx >= 0){
				dummy = exp(macht * log(dx));
			}
			else{
				dummy = exp(macht * log(-dx));
			}
			for (j = 1; j <= (2 * max); j++){
				a[i][j] = 0.0;
			}
			a[i][2 * i - 1] = (-coeff[1]) * dummy * macht / dx;
			a[i][2 * i] = 1.0;
			b[i][1] = (-dummy);
			if (origin_calc == 1){
				f[i] = -y[i] + orig_y + coeff[1] * dummy;
			}
			else{
				b[i][2] = (coeff[1]) * dummy * macht / dx;
				b[i][3] = -1.0;
				f[i] = -y[i] + coeff[3] + coeff[1] * dummy;
			}
		}
		for (i = 1; i <= max; i++){
			for (j = 1; j <= 2 * max; j++){
				aq[i][j] = 0.0;
				for (k = 1; k <= (2 * max); k++){
					aq[i][j] += a[i][k] * q[k][j];
				}
			}
		}
		for (i = 1; i <= max; i++){
			for (j = 1; j <= max; j++){
				aat[i][j] = 0.0;
				for (k = 1; k <= (2 * max); k++){
					aat[i][j] += aq[i][k] * a[j][k];
				}
			}
		}
		solve1(aat, max, index, &d);
		for (j = 1; j <= max; j++){
			for (i = 1; i <= max; i++) col[i] = 0.0;
			col[j] = 1.0;
			solve2(aat, max, index, col);
			for (i = 1; i <= max; i++) ainv[i][j] = col[i];
		}
		for (i = 1; i <= max; i++){
			for (j = 1; j <= n; j++){
				ainvb[i][j] = 0.0;
				for (k = 1; k <= max; k++){
					ainvb[i][j] += ainv[i][k] * b[k][j];
				}
			}
		}
		for (i = 1; i <= max; i++){
			ainvf[i] = 0.0;
			for (j = 1; j <= max; j++){
				ainvf[i] += ainv[i][j] * f[j];
			}
		}
		for (i = 1; i <= n; i++){
			for (j = 1; j <= n; j++){
				em[i][j] = 0.0;
				for (k = 1; k <= max; k++){
					em[i][j] += b[k][i] * ainvb[k][j];
				}
			}
		}
		for (i = 1; i <= n; i++){
			vect[i] = 0.0;
			for (j = 1; j <= max; j++){
				vect[i] += b[j][i] * ainvf[j];
			}
		}
		solve1(em, n, index, &d);
		solve2(em, n, index, vect);
		for (i = 1; i <= n; i++){
			coeff[i] += vect[i];
		}
		errorsum = 0;
		for (i = 1; i <= max; i++){
			dx = (origin_calc == 1 ? (x[i] - orig_x) : (x[i] - coeff[2]));
			dx = (fabs(dx) < SMALL ? SMALL : dx);
			if (dx >= 0){
				dummy = exp(macht * log(dx));
			}
			else{
				dummy = exp(macht * log(-dx));
			}
			if (origin_calc == 1){
				newy[i] = orig_y + coeff[1] * dummy;
			}
			else{
				newy[i] = coeff[3] + coeff[1] * dummy;
			}
			error[i] = (newy[i] - y[i]);
			errorsum += error[i] * error[i];
		}
	}
	mat1d_free(vect, 1, n);
	mat1d_free(f, 1, max);
	mat1d_free(ainvf, 1, max);
	mat1d_free(col, 1, max);
	int1d_free(index, 1, max);
	mat2d_free(em, 1, n, 1, n);
	mat2d_free(a, 1, max, 1, (2 * max));
	mat2d_free(aq, 1, max, 1, (2 * max));
	mat2d_free(q, 1, (2 * max), 1, (2 * max));
	mat2d_free(b, 1, max, 1, n);
	mat2d_free(ainv, 1, max, 1, max);
	mat2d_free(ainvb, 1, max, 1, n);
	mat2d_free(aat, 1, max, 1, max);
	return sqrt(errorsum / (float) max);
}


/**********************************************************************/

float correlation(x, y, max)

	float *x, *y;
	int max;
{
	float flag, xmean, ymean, varx, vary, varxy;
	int i;

	xmean = ymean = varx = vary = varxy = 0;
	for (i = 1; i <= max; i++){
		xmean += x[i];
		ymean += y[i];
	}
	xmean /= (float) max;
	ymean /= (float) max;
	for (i = 1; i <= max; i++){
		varx += (x[i] - xmean) * (x[i] - xmean);
		vary += (y[i] - ymean) * (y[i] - ymean);
		varxy += (x[i] - xmean) * (y[i] - ymean);
	}
	flag = varxy / (sqrt(varx) * sqrt(vary));
	return (flag * flag);
}


/**********************************************************************/

float rmse(x, y, max)

	float *x, *y;
	int max;
{
	float flag, varx;
	int i;

	varx = 0;
	for (i = 1; i <= max; i++){
		varx += (x[i] - y[i]) * (x[i] - y[i]);
	}
	flag = sqrt(varx / max);
	return (flag);
}


/**********************************************************************/

void solve1(amat, n, index, b)

	float **amat, *b;
	int n, *index;
{
	int i, imax, j, k;
	float large, dum, som, temp;
	float *vv;

	vv = mat1d(1, n);
	*b = 1.0;
	for (i = 1; i <= n; i++){
		large = 0.0;
		for (j = 1; j <= n; j++)
			if ((temp = fabs(amat[i][j])) > large) large = temp;
		vv[i] = 1.0 / large;
	}
	for (j = 1; j <= n; j++){
		for (i = 1; i < j; i++){
			som = amat[i][j];
			for (k = 1; k < i; k++) som -= amat[i][k] * amat[k][j];
			amat[i][j] = som;
		}
		large = 0.0;
		for (i = j; i <= n; i++){
			som = amat[i][j];
			for (k = 1; k < j; k++) som -= amat[i][k] * amat[k][j];
			amat[i][j] = som;
			if ((dum = vv[i] * fabs(som)) >= large){
				large = dum;
				imax = i;
			}
		}
		if (j != imax){
			for (k = 1; k <= n; k++){
				dum = amat[imax][k];
				amat[imax][k] = amat[j][k];
				amat[j][k] = dum;
			}
			*b = -(*b);
			vv[imax] = vv[j];
		}
		index[j] = imax;
		if (amat[j][j] == 0.0) amat[j][j] = SMALL;
		if (j != n){
			dum = 1.0 / (amat[j][j]);
			for (i = j + 1; i <= n; i++) amat[i][j] *= dum;
		}
	}
	mat1d_free(vv, 1, n);
}


/**********************************************************************/

void solve2(amat, n, index, b)

	float **amat, *b;
	int n, *index;
{
	int i, ii = 0, ip, j;
	float som;

	for (i = 1; i <= n; i++){
		ip = index[i];
		som = b[ip];
		b[ip] = b[i];
		if (ii)
			for (j = ii; j <= i - 1; j++) som -= amat[i][j] * b[j];
		else if (som) ii = i;
		b[i] = som;
	}
	for (i = n; i >= 1; i--){
		som = b[i];
		for (j = i + 1; j <= n; j++) som -= amat[i][j] * b[j];
		b[i] = som / amat[i][i];
	}
}


/**********************************************************************/

void readfiles(sname, tname)

	char *sname, *tname;
{
	char c;
	int i;

	printf("\n\n ->  Input file name  :  ");
	for (i = 0; (c = getchar()) != '\n'; ++i){
		sname[i] = c;
	}
	sname[i] = '\0';
	printf(" ->  Output file name  :  ");
	for (i = 0; (c = getchar()) != '\n'; ++i){
		tname[i] = c;
	}
	tname[i] = '\0';
	printf(" ->  Accuracy level power estimate (default = 0.01) : ");
	scanf("%f", &accuracy);
	printf(" ->  Fix origin x0 / y0 ? (1 = yes; 0 = no [default]) : ");
	scanf("%d", &origin_calc);
	if (origin_calc == 1){
		printf("    ->  x-coordinate origin (x0) : ");
		scanf("%f", &orig_x);
		printf("    ->  y-coordinate origin (y0) : ");
		scanf("%f", &orig_y);
	}
	printf(" ->  Use covariance matrix ? (1 = yes; 0 = no [default]) : ");
	scanf("%d", &variance);
	printf("\n");
}

/**********************************************************************/

float *mat1d(p,q)

	int p,q;
{
	float *vec;

	vec=(float *)malloc((unsigned) (q-p+1)*sizeof(float));
	return vec-p;
}


/**********************************************************************/

int *int1d(p, q)

	int p, q;
{
	int *vec;

	vec=(int *)malloc((unsigned) (q-p+1)*sizeof(int));
	return (vec-p);
}


/**********************************************************************/

void int1d_init(mat, p1, p2)

	int *mat;
	int p1, p2;
{
	int i;

	for (i=p1 ; i<=p2 ;i++){
		mat[i] = 0;
	}
}


/**********************************************************************/

float **mat2d(p1,p2,q1,q2)

	int p1,p2,q1,q2;
{
	int i;
	float **m;

	m=(float **) malloc((unsigned) (p2-p1+1)*sizeof(float*));
	m -= p1;
	for(i=p1;i<=p2;i++){
		m[i]=(float *) malloc((unsigned) (q2-q1+1)*sizeof(float));
		m[i] -= q1;
	}
	return m;
}


/**********************************************************************/

void mat1d_init(mat, p1, p2)

	float *mat;
	int p1, p2;
{
	int i;
	for (i=p1 ; i<=p2 ;i++){
		mat[i] = 0.;
	}
}


/**********************************************************************/

char *char1d(p,q)

	int p,q;
{
	char *vec;

	vec=(char *)malloc((unsigned) (q-p+1)*sizeof(float));
	return vec-p;
}


/**********************************************************************/

void mat2d_free(m, p1, p2, q1, q2)

	float **m;
	int p1, p2, q1, q2;
{
	int i;

	for (i = p2; i >= p1; i--)
		free((char*) (m[i] + q1));
	free((char*) (m + p1));
}


/**********************************************************************/

void mat1d_free(vec, p, q)

	float *vec;
	int p, q;
{
	free((char*) (vec + p));
}


/**********************************************************************/

void int1d_free(vec, p, q)

	int *vec, p, q;
{
	free((char*) (vec + p));
}


/**********************************************************************/

void mat2d_init(mat, p1, p2, q1, q2)

	float **mat;
	int p1, p2, q1, q2;
{
	int i,j;

	for (i=p1 ; i<=p2 ; i++){
		for (j=q1 ; j<=q2 ; j++){
			mat[i][j] =0.;
		}
	}
}


/**********************************************************************/

float str_to_f(char *s)
{
	float val, power;
	int i, sign;

	for (i = 0; s[i] == ' ' || s[i] == '\n' || s[i] == '\t'; i++)
		;
	sign = 1;
	if (s[i] == '+' || s[i] == '-')
		sign = (s[i++] == '+') ? 1 : -1;
	for (val = 0; s[i] >= '0' && s[i] <= '9'; i++)
		val = 10 * val + s[i] - '0';
	if (s[i] == '.')
		i++;
	for (power = 1; s[i] >= '0' && s[i] <= '9'; i++){
		val = 10 * val + s[i] - '0';
		power *= 10;
	}
	return(sign * val / power);
}


/**********************************************************************/



