#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double u(double x, double y)
{
	return 2.0 - (1.0 / (1.0 + exp(-x))) - (1.0 / (1.0 + exp(-y)));
}


int main(int argc, char** argv)
{
	int k = 0;
	int l = 1;
	double h1 = 0.3;
	int h2 = 0.25;
	double a = 0.5;
	double e = 0.1;
	int Tmax = 2;

	int m = (int)((double(l - k) / double(h1)));
	double n = pow(m - 1, 2.0);
	double r = h2 / pow(h1, 2);

	// Identity matrix
	double** I = (double**)malloc(n * sizeof(double *));
	for (int i = 0; i < n; ++i)
	{
		double* I_row = (double*)malloc(n * sizeof(double));
		I[i] = I_row;
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (i == j)
				I[i][j] = 1.0;
			else
				I[i][j] = 0.0;
		}
	}

	// A matrix
	double** A = (double**)malloc(n * sizeof(double *));
	for (int i = 0; i < n; ++i)
	{
		double* A_row = (double*)malloc(n * sizeof(double));
		A[i] = A_row;
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			A[i][j] = (-4.0) * I[i][j];
		}
	}

	// Set the values of A matrix
	for (int i = 0; i < (n - 1); ++i)
	{
		if ((i % (m - 1)) == 0)
		{
			A[i][i + 1] = 0;
			A[i + 1][i] = 0;
		}
		else
		{
			A[i][i + 1] = 1;
			A[i + 1][i] = 1;
		}
	}

	for (int i = 0; i <= (n - m); ++i)
	{
		A[i][m - 1 + i] = 1;
		A[m - 1 + i][i] = 1;
	}
	
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			A[i][j] = I[i][j] - (e * r * A[i][j]);
		}
	}

	// Create and the set the values of vector b1
	double* b1 = (double *)malloc(n * sizeof(double));
	for (int i = 0; i < n; ++i)
		b1[i] = 0;

	int index = 0;
	for (int i = 0; i <= m; i += m)
	{
		if (i == m)
		{
			k = n - m - 2;
		}
		for (int j = 0; j < (m - 1); ++j)
		{
			double x = i * h1;
			double y = j * h1;
			double o = u(x, y);
			b1[k] = o;
			k += 1;
		}
	}
	
	for (int i = 0; i < n; ++i)
		b1[i] = e * r * b1[i];

	// Create and the set the values of vector b2
	int step = m - 2;
	double* b2 = (double *)malloc(n * sizeof(double));
	for (int i = 0; i < n; ++i)
		b2[i] = 0;

	index = 0;
	for (int i = 0; i < (m - 1); ++i)
	{
		for (int j = 0; j <= m; j+=m)
		{
			double x = i * h1;
			double y = j * h1;
			double o = u(x, y);
			b2[k] = o;
			if (j == 0)
				k += step;
			else
				k += 1;
		}
	}

	for (int i = 0; i < n; ++i)
		b2[i] = e * r * b2[i];


	// Create and the set the values of vector b3
	double* b3 = (double *)malloc(n * sizeof(double));
	for (int i = 0; i < n; ++i)
		b3[i] = 0;

	index = 0;
	for (int i = 0; i < (m - 1); ++i)
	{
		for (int j = 0; j < (m - 1); ++j)
		{
			double x = i * h1;
			double y = j * h1;
			double o = u(x, y);
			b3[k] = o * (1 - o) * (o - a);
			k += 1;
		}
	}

	for (int i = 0; i < n; ++i)
		b3[i] = h2 * b3[i];

	// Create and the set the values of vector b
	double* b = (double *)malloc(n * sizeof(double));
	for (int i = 0; i < n; ++i)
		b[i] = b1[i] + b2[i] + b3[i];

	// K1 matrix of preconditioner
	double** K1 = (double**)malloc(n * sizeof(double*));
	for (int i = 0; i < n; ++i)
		double* K1_row = (double *)malloc(n * sizeof(double));

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (i == j)
				K1[i][j] = A[i][j];
			else
				K1[i][j] = 0;
		}
	}

	// TODO implement inverse of matrix analytically
	// K2 = inv(K1);

	// TODO implement the BiCGSTAB method
	// TODO complete the final solution of the Nagumo pde

	// Print the A matrix and the b vector
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			printf("%f ", A[i][j]);
		}
		printf("\n");
	}

	printf("\n");

	for (int i = 0; i < n; ++i)
	{
		printf("%f\n", b[i]);
	}

	printf("\n");

	// clear memory
	for (int i = 0; i < n; ++i)
	{
		free(I[i]);
		free(A[i]);
	}
	free(I);
	free(A);
	free(b1);
	free(b2);
	free(b3);
	free(b);

	return 0;
}