#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
// LINALG
typedef struct matrix_type {
  int m, n;
  double** v;
} mat_t, *mat;

mat matrix_new(int m, int n);

void matrix_delete(mat m);

void matrix_transpose(mat m);

mat matrix_copy(int n, double** a, int m);

mat matrix_mul(mat x, mat y);

mat matrix_minor(mat x, int d);

void matrix_show(mat m);

double* vmadd(double v[], double b[], double s, double c[], int n);
// c = a + b*s

double* mcol(mat m, double* v, int c);
// take the c-th column from the matrix m and put it into vector v
// vector v which is the column c in matrix m

mat vmul(double v[], int n);
// create matrix m = I - v * v^T

double vnorm(double x[], int n);

double* vdiv(double x[], double d, double y[], int n);

//
// MATRIX FUNCTIONS
//

mat matrix_new(int m, int n) {
  int i;
  mat x = (mat)malloc(sizeof(mat_t));
  x->v = (double**)malloc(sizeof(double) * m);
  x->v[0] = (double*)calloc(sizeof(double), m * n);
  for (i = 0; i < m; i++) x->v[i] = x->v[0] + n * i;
  x->m = m;
  x->n = n;
  return x;
}

void matrix_delete(mat m) {
  free(m->v[0]);
  free(m->v);
  free(m);
}

void matrix_transpose(mat m) {
  int i, j;
  for (i = 0; i < m->m; i++) {
    for (j = 0; j < i; j++) {
      double t = m->v[i][j];
      m->v[i][j] = m->v[j][i];
      m->v[j][i] = t;
    }
  }
}

mat matrix_copy(int n, double** a, int m) {
  int i, j;
  mat x = matrix_new(m, n);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) x->v[i][j] = a[i][j];
  return x;
}

mat matrix_mul(mat x, mat y) {
  int i, j, k;
  mat r;
  if (x->n != y->m) return 0;
  r = matrix_new(x->m, y->n);
  for (i = 0; i < x->m; i++)
    for (j = 0; j < y->n; j++)
      for (k = 0; k < x->n; k++) r->v[i][j] += x->v[i][k] * y->v[k][j];
  return r;
}

mat matrix_minor(mat x, int d) {
  int i, j;
  mat m = matrix_new(x->m, x->n);
  for (i = 0; i < d; i++) m->v[i][i] = 1;
  for (i = d; i < x->m; i++)
    for (j = d; j < x->n; j++) m->v[i][j] = x->v[i][j];
  return m;
}

//
// VECTOR FUNCTIONS
//

/* c = a + b * s */
double* vmadd(double a[], double b[], double s, double c[], int n) {
  int i;
  for (i = 0; i < n; i++) c[i] = a[i] + s * b[i];
  return c;
}

/* m = I - v v^T */
mat vmul(double v[], int n) {
  int i, j;
  mat x = matrix_new(n, n);
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) x->v[i][j] = -2 * v[i] * v[j];
  for (i = 0; i < n; i++) x->v[i][i] += 1;

  return x;
}

/* ||x|| */
double vnorm(double x[], int n) {
  int i;
  double sum = 0;
  for (i = 0; i < n; i++) sum += x[i] * x[i];
  return sqrt(sum);
}

/* y = x / d */
double* vdiv(double x[], double d, double y[], int n) {
  int i;
  for (i = 0; i < n; i++) y[i] = x[i] / d;
  return y;
}

/* take c-th column of m, put in v */
double* mcol(mat m, double* v, int c) {
  int i;
  for (i = 0; i < m->m; i++) v[i] = m->v[i][c];
  return v;
}

void matrix_show(mat m) {
  int i, j;
  for (i = 0; i < m->m; i++) {
    for (j = 0; j < m->n; j++) {
      printf(" %8.3f", m->v[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}
// END LINALG

/* ----------------------- norm ----------------------- */
/*  Given an array and its length, this function
    computes the 2-norm of the array.

    Input variables:
        x     : pointer to array for which the 2-norm should
                 be computed.
        length: number of entries in x.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.      */

double norm(double* x, int length) {
  int i, length5;
  double a, sum = 0;

  length5 = length % 5;

  for (i = 0; i < length5; i++) {
    sum += x[i] * x[i];
  }
  for (; i < length; i += 5) {
    sum += x[i] * x[i] + x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2] +
           x[i + 3] * x[i + 3] + x[i + 4] * x[i + 4];
  }

  return sqrt(sum);
}

/* ----------------------- vec_copy ----------------------- */
/*  Given two arrays of the same length and their length,
    this function stores the values from the first array
    in the second array.

    Input variables:
        x     : pointer to array whose entries are to be
                 copied.
        y     : pointer to array in which the components
                 of x are to be stored.
        length: number of entries in x and in y.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.          */

void vec_copy(double* x, double* y, int length) {
  int i, length5;

  length5 = length % 5;

  for (i = 0; i < length5; i++) {
    y[i] = x[i];
  }
  for (; i < length; i += 5) {
    y[i] = x[i];
    y[i + 1] = x[i + 1];
    y[i + 2] = x[i + 2];
    y[i + 3] = x[i + 3];
    y[i + 4] = x[i + 4];
  }
}

/* ------------------- partialvec_copy ------------------- */
/*  Given two arrays, the length of the second array, and
    an index this function stores the values from the
    subarray x[index : index + length] in the array
    y[0 : length].

    Input variables:
        x     : pointer to array whose entries are to be
                 copied.
        y     : pointer to array in which the components
                 of x are to be stored.
        length: number of entries in y.
        index : starting index of subarray of x to be
                copied to y.

    Example: Suppose x is a pointer to the array
    {1, 2, 3, 4, 5}, y is a pointer to the array {0, 0, 0},
    length = 3, and index = 2. Then after executing
    partialvec_copy(x, y, 3, 2), the array pointed to by
    y is now {3, 4, 5}.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.         */

void partialvec_copy(double* x, double* y, int length, int index) {
  int i, length5;

  length5 = length % 5;

  for (i = 0; i < length5; i++) {
    y[i] = x[i + index];
  }
  for (; i < length; i += 5) {
    y[i] = x[i + index];
    y[i + 1] = x[i + index + 1];
    y[i + 2] = x[i + index + 2];
    y[i + 3] = x[i + index + 3];
    y[i + 4] = x[i + index + 4];
  }
}

/* ----------------------- scalar_div ----------------------- */
/*  Given two arrays of the same length, their length, and a
    scalar value this function divides the values from the
    first array by the scalar value and stores the computed
    number in the second array.

    Input variables:
        x     : pointer to array whose components are to be
                 divided by r and stored in second array, y.
        r     : scalar used in division.
        length: number of entries in x and in y.
        y     : pointer to array in which the components
                 of x are to be stored.


    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.            */

void scalar_div(double* x, double r, int length, double* y) {
  int i, length5;

  length5 = length % 5;

  for (i = 0; i < length5; i++) {
    y[i] = x[i] / r;
  }
  for (; i < length; i += 5) {
    y[i] = x[i] / r;
    y[i + 1] = x[i + 1] / r;
    y[i + 2] = x[i + 2] / r;
    y[i + 3] = x[i + 3] / r;
    y[i + 4] = x[i + 4] / r;
  }
}

/* ----------------------- scalar_sub ----------------------- */
/*  Given two arrays of the same length, their length, and a
    scalar value this function multiplies the values from the
    first array by the scalar value and then subtracts the
    computed components from the components the second array.

    Input variables:
        x     : pointer to array whose components are to be
                 multiplied by r then subtracted from the
                 components of the second array, y.
        r     : scalar used in multiplication.
        length: number of entries in x and in y.
        y     : pointer to array in which the components
                 of x are to be stored.


    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.            */

void scalar_sub(double* x, double r, int length, double* y) {
  int i, length5;

  length5 = length % 5;

  for (i = 0; i < length5; i++) {
    y[i] -= r * x[i];
  }
  for (; i < length; i += 5) {
    y[i] -= r * x[i];
    y[i + 1] -= r * x[i + 1];
    y[i + 2] -= r * x[i + 2];
    y[i + 3] -= r * x[i + 3];
    y[i + 4] -= r * x[i + 4];
  }
}

/* --------------------- partialscalar_sub --------------------- */
/*  Given two arrays, the length of the second array, a scalar
    value, and an index, this function multiplies the values
    starting at the given index from the first array by the
    scalar value and then subtracts the computed components from
    the components the second array.

    Input variables:
        x     : pointer to array whose components are to be
                 multiplied by r then subtracted from the
                 components of the second array, y.
        r     : scalar used in multiplication.
        length: number of entries in y.
        index :
        y     : pointer to array in which the components
                 of x are to be stored.

    Example: Suppose x is a pointer to the array
    {1, 2, 3, 4, 5}, y is a pointer to the array {0, 0, 0},
    length = 3, r = -1, and index = 2. Then after executing
    partialscalar_sub(x, -1, 3, 2, y), the array pointed to
    by y is now {-3, -4, -5}.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.               */

void partialscalar_sub(double* x, double r, int length, int index, double* y) {
  int i, length5;

  length5 = length % 5;

  for (i = 0; i < length5; i++) {
    y[i + index] -= r * x[i];
  }
  for (; i < length; i += 5) {
    y[i + index] -= r * x[i];
    y[i + index + 1] -= r * x[i + 1];
    y[i + index + 2] -= r * x[i + 2];
    y[i + index + 3] -= r * x[i + 3];
    y[i + index + 4] -= r * x[i + 4];
  }
}

/* --------------------- dot_product --------------------- */
/*  Given two arrays of the same length and their length,
    this function returns the dot product of the two
    arrays.

    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in x and in y.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.         */

double dot_product(double* x, double* y, int length) {
  int i, length5;
  double sum = 0;

  length5 = length % 5;

  for (i = 0; i < length5; i++) {
    sum += x[i] * y[i];
  }
  for (; i < length; i += 5) {
    sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] +
           x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
  }

  return sum;
}

/* ------------------ partialdot_product ------------------ */
/*  Given two arrays of the same length, their length, and
    an index this function returns the dot product of the
    two subarrays x[index : length] and y[index : length].

    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in x and in y.
        index : starting index for subarrays.

    Example: Suppose x is a pointer to the array
    {1, 2, 3, 4}, y is a pointer to the array {5, 6, 7, 8},
    length = 4, and index = 2. Then the value returned by
    executing partialdot_product(x, y, 4, 2) is 53, which
    is computed by
        x[2] * y[2] + x[3] * y[3] = 3 * 7 + 4 * 8
                                  = 21 + 32
                                  = 53.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.          */

double partialdot_product(double* x, double* y, int length, int index) {
  int i, length5;
  double sum = 0;

  length5 = length % 5;

  for (i = index; i < length5; i++) {
    sum += x[i] * y[i];
  }
  for (; i < length; i += 5) {
    sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] +
           x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
  }

  return sum;
}

/* -------------------- subdot_product -------------------- */
/*  Given two arrays, the length of the second array, and
    an index this function returns the dot product of the
    two subarrays x[index : index + length] and
    y[0 : length]. It is necessary that index + length is
    at most the length of the first array.

    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in y.
        index : starting index for subarray of x.

    Example: Suppose x is a pointer to the array
    {1, 2, 3, 4, 5}, y is a pointer to the array
    {-1, -2, -3}, length = 3, and index = 2. Then the value
    returned by executing subdot_product(x, y, 3, 2) is 53,
    which is computed by
            x[2] * y[0] + x[3] * y[1] + x[4] * y[2]

          =  3   *  -1  +  4   *  -2  +  5   *  -3

          = -    3      -      8      -      15

          = -26.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.          */

double subdot_product(double* x, double* y, int length, int index) {
  int i, length5;
  double sum = 0;

  length5 = length % 5;

  for (i = 0; i < length5; i++) {
    sum += x[i + index] * y[i];
  }
  for (; i < length; i += 5) {
    sum += x[i + index] * y[i] + x[i + index + 1] * y[i + 1] +
           x[i + index + 2] * y[i + 2] + x[i + index + 3] * y[i + 3] +
           x[i + index + 4] * y[i + 4];
  }

  return sum;
}

/* ----------------------- pth_householder ----------------------- */
/*  Given a matrix A of dimension m by n (with n <= m) and
    arrays v_i of dimension m-i, for i = 1, ..., m - 1,
    respectively, this algorithm computes n reflection vectors
    and the factor R of a full QR decomposition of A, where R
    is a m by n upper triangular matrix. The n reflection
    vectors are stored in the arrays v_1, ..., v_n and the
    columns of A are overwritten by the columns of R.

    Input variables:
        a: pointer to array of arrays, the ith array of
            which should correspond to the ith column of the
            matrix A. During the algorithm, the columns of R
            will overwrite the columns of A.
        v: pointer to array of arrays in which the ith
            reflection vector of dimension m - i will be
            stored.
        m: number of rows of A.
        n: number of columns of A.

    Features: The number of flops for this implementation is
    ~ 2 * m * n^2 - (2/3) * n^3 and requires O(1) additional
    memory.

                                                   */
///
/// PARALLEL SECTION
///
struct args_t {
  int m;
  int n;
  mat a;
  mat v;
  int id;
};
void* worker(void* data) {
  args_t* d = (args_t*)data;
  mat a = d->a;
  mat v = d->v;
  int m = d->m;
  int n = d->n;
  int i = d->id;

  double vnorm, vTa, vpartdot;
  partialvec_copy(a->v[i], v->v[i], m - i, i);

  /* vpartdot = ||v->v[i]||^2 - v->v[i][0] * v->v[i][0]; since vpartdot
         is unaffected by the change in v->v[i][0], storing this value
         prevents the need to recalculate the entire norm of v->v[i]
         after updating v->v[i][0] in the following step              */
  vpartdot = partialdot_product(v->v[i], v->v[i], m - i, 1);

  /* set v->v[i][0] = v->v[i][0] + sign(v->v[i][0]) * ||v->v[i]|| */
  if (v->v[i][0] < 0) {
    v->v[i][0] -= sqrt(v->v[i][0] * v->v[i][0] + vpartdot);
  } else {
    v->v[i][0] += sqrt(v->v[i][0] * v->v[i][0] + vpartdot);
  }

  /* normalize v->v[i] */
  vnorm = sqrt(v->v[i][0] * v->v[i][0] + vpartdot);
  scalar_div(v->v[i], vnorm, m - i, v->v[i]);

  for (int j = i; j < n; j++) {
    /* set a->v[j][i:m] = a->v[j][i:m] - 2 * (v->v[i]^T a->v[j][i:m]) * v->v[i]
     */
    vTa = subdot_product(a->v[j], v->v[i], m - i, i);
    vTa *= 2;
    partialscalar_sub(v->v[i], vTa, m - i, i, a->v[j]);
  }
  pthread_exit(NULL);
}

void pth_householder(mat a, mat v, int m, int n) {
  int i, j;
  args_t* d = (args_t*)calloc(n, sizeof(args_t));
  pthread_t* threads = (pthread_t*)calloc(n, sizeof(pthread_t));

  for (i = 0; i < n; i++) {
    /* set v->v[i] equal to subvector a->v[i][i : m] */
    d[i].a = a;
    d[i].id = i;
    d[i].v = v;
    d[i].m = m;
    d[i].n = n;
  }

  for (i = 0; i < n; i++) {
    pthread_create(&threads[i], NULL, worker, (void*)&d[i]);
  }

  for (i = 0; i < n; i++) {
    pthread_join(threads[i], NULL);
  }
  free(d);
  free(threads);
}
/// END PARALLEL SECTION

/// SEQ pth_householder
void householder(mat a, mat v, int m, int n) {
  int i, j;
  double vnorm, vTa, vpartdot;

  for (i = 0; i < n; i++) {
    double vnorm, vTa, vpartdot;
    partialvec_copy(a->v[i], v->v[i], m - i, i);

    /* vpartdot = ||v->v[i]||^2 - v->v[i][0] * v->v[i][0]; since vpartdot
       is unaffected by the change in v->v[i][0], storing this value
       prevents the need to recalculate the entire norm of v->v[i]
       after updating v->v[i][0] in the following step              */
    vpartdot = partialdot_product(v->v[i], v->v[i], m - i, 1);

    /* set v->v[i][0] = v->v[i][0] + sign(v->v[i][0]) * ||v->v[i]|| */
    if (v->v[i][0] < 0) {
      v->v[i][0] -= sqrt(v->v[i][0] * v->v[i][0] + vpartdot);
    } else {
      v->v[i][0] += sqrt(v->v[i][0] * v->v[i][0] + vpartdot);
    }

    /* normalize v->v[i] */
    vnorm = sqrt(v->v[i][0] * v->v[i][0] + vpartdot);
    scalar_div(v->v[i], vnorm, m - i, v->v[i]);

    for (int j = i; j < n; j++) {
      /* set a->v[j][i:m] = a->v[j][i:m] - 2 * (v->v[i]^T a->v[j][i:m]) *
       * v->v[i] */
      vTa = subdot_product(a->v[j], v->v[i], m - i, i);
      vTa *= 2;
      partialscalar_sub(v->v[i], vTa, m - i, i, a->v[j]);
    }
  }
}
/// END SEQ pth_householder

#include <cstring>
int main(int argc, char** argv) {
  int VER = 0;
  int i, j, m, n;
  double x;
  if (argc > 1 && !strcmp(argv[1], "-v")) {
    printf("Running Verbousely!!\n");
    VER = 1;
  }

  /* let user set the dimension of matrix A */
  printf("Enter the dimension m (where A is a m by n matrix): ");
  scanf("%i", &m);
  printf("Enter the dimension n (where A is a m by n matrix): ");
  scanf("%i", &n);

  /* check if m < n */
  if (m < n) {
    printf(
        "For a successful factorization, this implementation "
        "requires n <= m.\nTerminating program.\n");
    return 0;
  }

  /* allocate memory for A and vectors v */
  mat a = matrix_new(m, n);
  
  mat v = matrix_new(m, n);
  // double ** v = new double*[n];
  // for(i = 0; i < n; i++) {
  //     a->v[i] = new double[m];
  //     v->v[i] = new double[m - i];
  // }

  /* initialize the values in matrix A */
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      if (j < i) {
        a->v[i][j] = 0;
      } else {
        a->v[i][j] = 1;  // j - i + 1; // this choice of values was arbitrary
      }
    }
  }

  mat old_a = matrix_copy(a->m, a->v, a->n);

  /* print the matrix A before calling houheholder */
  if (VER) {
    printf("A = \n");
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        printf("%9.6g ", a->v[j][i]);
      }
      printf("\n");
    }
    printf("\n");
  }

  clock_t start, end;
  start = clock();
  pth_householder(a, v, m, n);
  end = clock();

  /* print the matrix R (stored in A) after calling houheholder */
  if (VER) {
    printf("R = \n");
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        printf("%9.6g ", a->v[j][i]);
      }
      printf("\n");
    }
    printf("\n");

    /* print the vectors v after calling pth_householder */
    for (i = 0; i < n; i++) {
      printf("v[%i] = ", i);
      for (j = 0; j < m - i; j++) {
        printf("%9.6g ", v->v[i][j]);
      }
      printf("\n");
    }
    printf("\n");

    /* print numerical evidence that v's are normalized */
    printf(
        "Numerical verification that v_1, ..., v_%i are "
        "normalized:\n",
        n);
    for (i = 1; i < n; i++) {
      x = dot_product(v->v[i - 1], v->v[i - 1], m - i + 1);
      printf("||v[%i]|| = %lg, ", i, x);
      if (i % 5 == 0) {
        printf("\n");
      }
    }
    x = dot_product(v->v[n - 1], v->v[n - 1], m - n + 1);
    printf("||v[%i]|| = %lg.", n, x);

    if (n % 5 != 0) printf("\n");
  }

  printf("\n=====================\n");

  printf("It took %f seconds for a parallel algorithm to run!\n",
         double(end - start) / CLOCKS_PER_SEC);
  
  start = clock();
  householder(old_a, v, m, n);
  end = clock();
  printf("It took %f seconds for a sequential algorithm to run!\n",
         double(end - start) / CLOCKS_PER_SEC);
  matrix_delete(a);
  matrix_delete(v);
  matrix_delete(old_a);
  // /* free memory */
  // for(i = 0; i < n; i++) {
  //     delete[] a->v[i];
  //     delete[] v[i];
  // }
  // delete[] a;
  // delete[] v;
  return 0;
}
