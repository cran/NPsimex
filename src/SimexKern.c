/*#include <stdio.h>
#include <math.h>  */
#include <R.h>
#include <Rmath.h>

void SimexKern(double *x, int *n, double *tt, int *ntt, double *ll, int *nll, double *msigma, double *result){
	int i, j, k; 
	double d, ksum;
	for(k=0; k < *nll; k++){
		for(i=0; i < *ntt; i++){
			ksum = 0;
			for(j=0; j < *n; j++){
				d = tt[i] - x[j];
				ksum += dnorm(d / ((*msigma) * ll[k]), 0, 1, 0);
			}
			result[k*(*ntt)+i] = ksum / ((*n) * ((*msigma) * ll[k]));
		}		
	}
}

void SimexKernH(double *x, int *n, double *tt, int *ntt, double *ll, int *nll, double *msigma, double *result){
	int i, j, k; 
	double d, ksum;
	for(k=0; k < *nll; k++){
		for(i=0; i < *ntt; i++){
			ksum = 0;
			for(j=0; j < *n; j++){
				d = tt[i] - x[j];
				ksum += dnorm(d / (msigma[j] * ll[k]), 0, 1, 0) / (msigma[j] * ll[k]);
			}
			result[k*(*ntt)+i] = ksum / (*n);
		}		
	}
}
