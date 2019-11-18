#include<stdio.h>
#include<math.h>
double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13;

void sac(double x1, double x2, double x3, double *y){
	v1 = x1;
	v2 = x2;
	v3 = x3;
	v4 = v1 * v2;
	v5 = sin(v4);
	v6 = pow(v1, 2);
	v7 = v6 * v3;
	v8 = pow(v1, 3);
	v9 = v2 - v3;
	v10 = exp(v9);
	v11 = v10 * v8;
	v12 = v7 - v11;
	v13 = v5 + v12;
	*y = v13;
}

void forward_mode(double dx1, double dx2, double dx3, double *y){
	double dv1,dv2,dv3,dv4,dv5,dv6,dv7,dv8,dv9,dv10,dv11,dv12,dv13;
	dv1 = dx1;
    	dv2 = dx2;
    	dv3 = dx3;
	dv4 = dv1 * v2 + v1 * dv2;
	dv5 = dv4 * cos(v4);
	dv6 = 2 * v1 * dv1;
	dv7 = dv6 * v3 + v6 * dv3;
	dv8 = 3 * pow(v1, 2) * dv1;
	dv9 = dv2 - dv3;
	dv10 = exp(v9) * dv9;
	dv11 = dv8 * v10 + v8 * dv10;
	dv12 = dv7 - dv11;
	dv13 = dv5 + dv12;
	*y = dv13;
}

void reverse_mode(double x1, double x2, double x3){
	double v1, v2, v3, dv1,dv2,dv3,dv4,dv5,dv6,dv7,dv8,dv9,dv10,dv11,dv12,dv13;
	v1 = x1;
	v2 = x2;
	v3 = x3;

	dv13 = 1;
	dv12 = dv13;
	dv11 = -dv12;
	dv10 = dv11 * v8;
	dv9 = dv10 * exp(v9);
	dv8 = dv11 * v10;
	dv7 = dv12;
	dv6 = dv7 * v3;
	dv5 = dv13;
	dv4 = dv5 * cos(v4);
	dv3 = dv7 * v6 - dv9;
	dv2 = dv4 * v1 + dv9;
	dv1 = dv8 * 3 * pow(v1, 2) + dv4 * v2 + dv6 * 2 * v1;
	printf("\ndy/dx1 = %lf\n", dv1);
	printf("\ndy/dx2 = %lf\n", dv2);
	printf("\ndy/dx3 = %lf\n", dv3);
}

int main(){
double x1, x2, x3;
double y1, y2, y3, y4, y5, y6, y7;
printf("\nf(x) = sin(x1 * x2) + x1^2 * x3 - x1^3 * exp(x2 - x1");
printf("\nEnter x1 : ");
scanf("%lf", &x1);
printf("\nEnter x2 : ");
scanf("%lf", &x2);
printf("\nEnter x3 : ");
scanf("%lf", &x3);
sac(x1, x2, x3, &y1);
printf("\n---------SAC-----------");
printf("\nOutput : f(x) = %lf\n", y1);

printf("\n---------Forward Mode-----------");
forward_mode(1.0, 0.0, 0.0, &y2);
printf("\ndy/dx1 = %lf\n", y2);
forward_mode(0.0, 1.0, 0.0, &y3);
printf("\ndy/dx2 = %lf\n", y3);
forward_mode(0.0, 0.0, 1.0, &y4);
printf("\ndy/dx3 = %lf\n", y4);

printf("\n---------Reverse Mode-----------");
reverse_mode(x1, x2, x3);
return 0;
}
