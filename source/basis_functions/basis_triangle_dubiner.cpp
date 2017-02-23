#include <iostream>
#include <cmath>
using namespace std;

#include "basis_functions.h"

void dubiner_phi(int p, int q, int number_gp, double n1[], double n2[], double phi_pq[]){
	double* psi_p = new double [number_gp];
	double* psi_pq = new double [number_gp];

	jacobi_polynomial(p,0,0,number_gp,n1,psi_p);
	jacobi_polynomial(q,2*p+1,0,number_gp,n2,psi_pq);

	for (int i=0; i<number_gp;i++){
		phi_pq[i] = psi_p[i]*pow((1-n2[i])/2,p)*psi_pq[i];
	}

	delete [] psi_p;
	delete [] psi_pq;
}

void dubiner_test(int p, int number_gp, double** phi_area, double* w){
	int number_bf = (p+1)*(p+2)/2;

	double** M = new double*[number_bf];
	
	for (int i=0; i<number_bf; i++){
		M[i] = new double[number_bf];
		for (int j=0; j<number_bf; j++){
			M[i][j] = 0;
			for (int k=0; k<number_gp; k++){
				M[i][j] = M[i][j] + w[k]*phi_area[i][k]*phi_area[j][k];
			}
		}
	}

	int m = 0;
	double M_exact;
	for (int i=0; i<=p; i++){
		for (int j=0; j<=p-i;j++){
			M_exact = 2.0/((2*i+1)*(i+j+1));
			M[m][m] = abs((M[m][m]-M_exact)/M_exact);
			m = m + 1;
		}
	}

	for (int i=0; i<number_bf; i++){
		for (int j=0; j<number_bf; j++){
			if ( abs(M[i][j]) > pow(10.0,-10.0))
			{
				cout << "\n";
				cout << "DUBINER_BASIS - Test fail!\n";
				cout << "  (i,j) = " << i << " , " << j << "\n";
				exit ( 1 );
			}
		}
	}

	cout << "\n";
	cout << "DUBINER_BASIS - Test success!\n";

	for (int i=0; i<number_bf; i++){
		delete [] M[i];
	}
	delete [] M;
}
