#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define ROWS 3
#define COLS 3

double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18;

void sac(double x1, double x2, double x3, double *y1, double *y2, double *y3){
	v1 = x1;
	v2 = x2;
	v3 = x3;
	v4 = pow(v1, 2);
	v5 = v4 * v2;
	v6 = v1 * v2;
	v7 = v6 * v3;
	v8 = sin(v7);
	v9 = v5 + v8;
	v10 = v9 + 0.415145069;
	v11 = fabs(v7);
	v12 = exp(v11);
	v13 = v12 + v7;
	v14 = v13 - 1.066476035;
	v15 = v1 * v2;
	v16 = pow(v7,2);
	v17 = v16 + v15;
	v18 = v17 + 0.194335938;
	*y1 = v10;
	*y2 = v14;
	*y3 = v18;
}

void derivatives(double dx1, double dx2, double dx3, double *y1, double *y2, double *y3){
	double dv1,dv2,dv3,dv4,dv5,dv6,dv7,dv8,dv9,dv10,dv11,dv12,dv13, dv14, dv15, dv16, dv17, dv18;
	dv1 = dx1;
    	dv2 = dx2;
    	dv3 = dx3;
	dv4 = 2 * v1 * dv1;
	dv5 = v4 * dv2 + v2 * dv4;
	dv6 = v1 * dv2 + v2 * dv1;
	dv7 = v6 * dv3 + v3 * dv6;
	dv8 = cos(v7) * dv7;
	dv9 = dv5 + dv8;
	dv10 = dv9;
	dv11 = (v7/fabs(v7)) * dv7;
	dv12 = exp(v11) * dv11;
	dv13 = dv12 + dv7;
	dv14 = dv13;
	dv15 = v1 * dv2 + v2 * dv1;
	dv16 = 2 * v7 * dv7;
	dv17 = dv16 + dv15;
	dv18 = dv17;
	*y1 = dv10;
	*y2 = dv14;
	*y3 = dv18;
}

void jacobian_matrix(double matrix[ROWS][COLS]){
	//First order derivatives of f1 wrt x1, x2 and x3
	derivatives(1.0, 0.0, 0.0, (*(matrix + 0) + 0), (*(matrix + 1) + 0), (*(matrix + 2) + 0));

	//First order derivatives of f2
	derivatives(0.0, 1.0, 0.0, (*(matrix + 0) + 1), (*(matrix + 1) + 1), (*(matrix + 2) + 1));

	//First order derivatives of f3
	derivatives(0.0, 0.0, 1.0, (*(matrix + 0) + 2), (*(matrix + 1) + 2), (*(matrix + 2) + 2));
}

void print_matrix(double (*matrix)[COLS], int rows, int cols){
	int i, j;
	for(i = 0; i<rows; i++){
                for(j=0; j<cols; j++){
                        printf("%lf\t", matrix[i][j]);
                }
                printf("\n");
        }
}

void minor(double b[ROWS][COLS],double a[ROWS][COLS],int i,int n){
	int j,l,h=0,k=0;
	for(l=1;l<n;l++)
		for( j=0;j<n;j++){
			if(j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if(k == (n-1)){
				h++;
				k=0;
			}
		}
}

void cofactor(double mat[ROWS][COLS], double temp[ROWS][COLS], int p, int q, int n) 
{ 
    	int i = 0, j = 0, row, col; 
 
    	for (row = 0; row < n; row++){ 
        	for (col = 0; col < n; col++){ 
            		if (row != p && col != q){ 
                		temp[i][j++] = mat[row][col]; 

                		if (j == n - 1) { 
                    			j = 0; 
                    			i++; 
                		} 
            		}
        	}
    	}
}

double determinant(double mat[ROWS][COLS], int n)
{
        double D = 0;
        double temp[ROWS][COLS];
        double sign = 1;

        if (n == 1){
                return mat[0][0];
        }
        for (int f = 0; f < n; f++){
                cofactor(mat, temp, 0, f, n);
                D += sign * mat[0][f] * determinant(temp, n - 1);
                sign = -sign;
        }
        return D;
}

void adjoint(double A[ROWS][COLS], double adj[ROWS][COLS]) 
{ 
	int sign = 1;
	double temp[ROWS][COLS];

	if (ROWS == 1 && COLS == 1){ 
        	adj[0][0] = 1; 
        	return; 
    	} 

    	for (int i=0; i<ROWS; i++) {
        	for (int j=0; j<COLS; j++){
            		cofactor(A, temp, i, j, ROWS);
            		sign = ((i+j)%2==0)? 1: -1;
            		adj[j][i] = (sign)*(determinant(temp, ROWS - 1)); 
        	}
    	}
}
  
int inverse(double A[ROWS][COLS], double inverse[ROWS][COLS]) 
{
	double adj[ROWS][COLS];
    	double det = determinant(A, ROWS);

	if (det == 0)
    	{ 
        	printf("Can't find the inverse of singular matrix!"); 
        	return 1;
    	}

    	adjoint(A, adj);
    	for (int i=0; i<ROWS; i++) 
        	for (int j=0; j<COLS; j++) 
            		inverse[i][j] = adj[i][j]/det; 
    	return 0;
}

void newton_dampning(double x1, double x2, double x3, double alpha, int iter)
{
	double temp1[3], temp2[3];
	double matrix[ROWS][COLS], invmat[ROWS][COLS];
	int i, j, k;

	for(k=1; k<iter; k++){
		sac(x1, x2, x3, (temp1+0), (temp1+1), (temp1+2));
		jacobian_matrix(matrix);
		inverse(matrix, invmat);
		for(i=0; i<ROWS; i++){
			for(j=0; j<COLS; j++){
				temp2[i] += invmat[i][j] * temp1[j];
			}
		}
		for(i=0; i<3; i++){
			temp1[i] -= alpha * temp2[i];
		}
		printf("\nThe value after iteration %d is : %lf\t%lf\t%lf", k, temp1[0], temp1[1], temp1[2]);
	}
}

int main(){

	int inv, i, j;

	double x1, x2, x3, x1_0, x2_0, x3_0;
	double y1, y2, y3, alpha;
	double determin;

	double matrix[ROWS][COLS], invmat[ROWS][COLS];
	double sacmat[ROWS];
	
	printf("\n******* Computing the Jacobian ********\n");
	printf("\nEnter x1 : ");
	scanf("%lf", &x1);
	printf("\nEnter x2 : ");
	scanf("%lf", &x2);
	printf("\nEnter x3 : ");
	scanf("%lf", &x3);

	// call single assignment code
        sac(x1, x2, x3, (sacmat + 0), (sacmat + 1), (sacmat + 2));

        printf("\nFunction outputs:\n");
        for(i=0; i<COLS; i++)
                printf("%lf\t", *(sacmat + i));

        // compute the jacobian and inverse jacobian matrices
        jacobian_matrix(matrix);
        printf("\n\nJacboian Matrix:\n");
	print_matrix(matrix, ROWS, COLS);
	
	// inverse matrix
	inv = inverse(matrix, invmat);
	printf("\nInverse matrix:\n");
	print_matrix(invmat, ROWS, COLS);

	printf("\n********** Dampned Newton Method **********\n");
	//damped newton method
	printf("\nEnter the initial values for the vector:\n");
	printf("\nEnter x1_0 : ");
        scanf("%lf", &x1_0);
        printf("\nEnter x2_0 : ");
        scanf("%lf", &x2_0);
        printf("\nEnter x3_0 : ");
        scanf("%lf", &x3_0);
	printf("\nEnter alpha value bw the interval (0, 1] : ");
	scanf("%lf", &alpha);
	newton_dampning(x1_0, x2_0, x3_0, alpha, 10);
	printf("\n");
	return 0;
}
