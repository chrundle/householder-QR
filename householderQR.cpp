#include <cmath>
#include <stdlib.h>
#include <stdio.h>
//#include <iostream>
#include "../linearalgebra.h"


void householder (double ** a, double ** v, int m, int n) {
    int i, j;
    double vnorm, vTa;

    for(i = 0; i < n; i++) {
        /* set v[i] equal to subvector a[i][i : m] */
        subvec_copy(a[i], v[i], m - i, i);

        /* set v[i][0] = v[i][0] + sign(v[i][0]) * ||v[i]|| */
        if(v[i][0] < 0) {
            v[i][0] -= norm(v[i], m - i);
        }
        else {
            v[i][0] += norm(v[i], m - i);
        }

        /* normalize v[i] */
        vnorm = norm(v[i], m - i);
        scalar_mult(v[i], vnorm, m - i, v[i]);
    
        for(j = i; j < n; j++) {
            /* set a[j][i:m] = a[j][i:m] - 2 * (v[i]^T a[j][i:m]) * v[i] */
            vTa = subdot_product(a[j], v[i], m - i, i);
            vTa *= 2;
            subscalar_sub(v[i], vTa, m - i, i, a[j]);
        }
    }
}


int main () {
    int i, j, m, n;
    double x;

    /* let user set the dimension of matrix A */
    printf("Enter the dimension m (where A is a m by n matrix): ");
    scanf("%i", &m);
    printf("Enter the dimension n (where A is a m by n matrix): ");
    scanf("%i", &n);

    /* check if m < n */
    if(m < n) {
        printf("For a successful factorization, we require n <= m. "
               "Terminating program.\n");
        return 0;
    }

    /* allocate memory for A and vectors v */
    double ** a = new double*[n];
    double ** v = new double*[n];
    for(i = 0; i < n; i++) {
        a[i] = new double[m];
        v[i] = new double[m - i];
    }

    /* initialize the values in matrix A */
    for(i = 0; i < n; i++) {
        for(j = 0; j < m; j++) {
            if(j < i) {
                a[i][j] = 0;
            }
            else {
                a[i][j] = j - i + 1; // this choice of values was arbitrary
            }
        }
    }

    /* print the matrix A before calling houheholder */
    printf("A = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {

            printf("%9.6g ", a[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    householder(a, v, m, n);

    /* print the matrix R (stored in A) after calling houheholder */
    printf("R = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < n; j++) {

            printf("%9.6g ", a[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    /* print the vectors v after calling householder */
    for(i = 0; i < n; i++) {
        printf("v[%i] = ", i);
        for(j = 0; j < m - i; j++) {
            printf("%9.6g ", v[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    /* print numerical evidence that v's are normalized */
    printf("Numerical verification that v_1, ..., v_%i are "
           "normalized:\n", n);
    for(i = 0; i < n; i++) {
        x = dot_product(v[i], v[i], m - i);
        printf("v_%i * v_%i = %lg\n", i + 1, i + 1, x);
    }

    /* free memory */
    for(i = 0; i < n; i++) {
        delete[] a[i];
        delete[] v[i];
    }
    delete[] a;
    delete[] v;
    return 0;
}
