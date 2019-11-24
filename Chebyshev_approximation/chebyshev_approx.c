#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define N 10

#ifndef PI
#define PI 3.14159265358979323846
#endif

double v1, v2, v3, v4, v5;

double func(double x)
/*
single assignment code for computing function output:
y = F(x) = x^2 + x * sin(x)
*/
{
	double out;
	v1 = x;
	v2 = pow(v1, 2);
	v3 = sin(v1);
	v4 = v3 * v1;
	v5 = v4 + v2;
	out = v5;
	return out;	
}

double derivative(double x)
{
	double dv1, dv2, dv3, dv4, dv5;
	func(x);
	dv1 = 1;
	dv2 = 2 * dv1;
	dv3 = cos(v1) * dv1;
	dv4 = v3 * dv1 + v1 * dv3;
	dv5 = dv4 + dv2;
	return dv5;
}

double polynomial(double x, int iter)
/*
computes Chebyshev polynomial for the value of x and iter.
x : the value of x passed to the polynomial
iter : number of iterations to be performed 
Computes:
T_0 = 1
T_1 = x
.
.
T_n(x) = 2 * x * T_n-1(x) - T_n-2(x) , for n > 1
*/
{
	if(iter == 0) return 1;
	else if(iter == 1) return x;
	else return (2 * x * polynomial(x, iter-1)) - polynomial(x, iter - 2);
}

void coefficient(double (*func)(double), double x, int n, double *out)
/*
computes the coefficient c_j for j = 0 .. N-1
Formula :

c_j = 2/N *SUM_(k=1 to N) f(cos(pi*(k-0.5)/N)) * cos(pi * j * (k-0.5)/N))
*/
{
	int i, j;
	double f_out;
	// set the array values to 0 each time before calculation
	memset(out, 0, sizeof(double) * n);
	for(i = 1; i <= n; i++){
		f_out = func(cos(PI * (i - 0.5)/n))*2/n;
		for(j = 0; j < n; j++){
			*(out + j) += f_out * cos(PI * j * (i - 0.5)/n);
		}
	}
}

double approximation(double x, int n, double *c)
/*
computes the approximation function : 

f(x) = (SUM_(k=0 to n-1) c_k * T_k(x)) - 0.5 * c_0
*/
{
	int k;
	double approx = c[0] + c[1] * polynomial(x, 1);
	for(k=2; k<n; k++){
		approx += c[k] * polynomial(x, k);
	}
	approx -= 0.5 * c[0];
	return approx;
}

double limit(int a, int b, double x)
/*
 gerneralization of the equation for intervals [a, b]
*/
{
	return (x - 0.5 * (b - a))/(0.5 * (b + a));
}


int main()
{
	double res, x, x_dash, coef[N], dvcoef[N];
	int iter;
	FILE *f=NULL;

	// write the output into a file to be used for gnuplot
	f = fopen("plot.txt", "w");

	printf("\n######## Chebyshev Approximations ######### \n");
	printf("\nx \t\t x_dash\t\t F(x) \t\t Approx F(x) \t F`(x) \t\t Approx F`(x)\n");
	for(x = 1.0; x < 2.0; x += 0.1){
		x_dash = limit(1, 2, x);
		coefficient(func, x_dash, N, coef);
		coefficient(derivative, x_dash, N, dvcoef);
		printf("%lf \t %lf \t %lf \t %lf \t %lf \t %lf\n", x, x_dash, func(x_dash), approximation(x_dash, N, coef), derivative(x_dash), approximation(x_dash, N, dvcoef));
		fprintf(f, "%lf \t %lf \t %lf \t %lf\n", func(x_dash), approximation(x_dash, N, coef), derivative(x_dash), approximation(x_dash, N, dvcoef));
	}
	// close the stream
	fclose(f);
	printf("\n");
	return 0;
}
