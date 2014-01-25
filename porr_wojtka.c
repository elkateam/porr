#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
int size = 0;
double max = -DBL_MAX;
double **input_matrix;


void read_matrix();
void display_matrix(double **matrix);
void die(char *tekst);

void read_matrix() {
    int ok = 0; // zeby sprawdzic czy dane zostaly wczytane ok

    ok = scanf("%d", &size);
    if (!ok || size < 1)
        die("\nBlad wejscia size");

    input_matrix = (double **) malloc(size * sizeof (double *));
    int i = 0;
    int j = 0;
    for (i = 0; i < size; i++) {
        input_matrix[i] = (double *) malloc((size + 1) * sizeof (double));
    }
    for (i = 0; i < size; i++) {
        for (j = 0; j < size + 1; j++) {
            ok = scanf("%lf", &input_matrix[i][j]);
            if (!ok || size < 1)
                die("\nBlad wejscia - wczytywanie tablicy");
            if (max < input_matrix[i][j] && j < size)
                max = input_matrix[i][j];
        }
    }
}

void display_matrix(double **matrix) {
    int i, j;
    // printf("Wypisuje macierz:\n");
    for (i = 0; i < size; i++) {
        for (j = 0; j < size + 1; j++) {
            printf("%lf ", matrix[i][j]);
            // printf("[%d, %d] %lf ",i , j, matrix[i][j]);
        }
        printf("\n");
    }

}

void die(char *tekst) {
    printf("%s\a", tekst); // komunikat o bledzie + sygnal dzwiekowy beep
    exit(EXIT_FAILURE); // konczymy z bledem
}

long int cheb_parallel(double** A, int size, double alfa, double beta, int s, int max_iter, double error, int parallel) {
    int i, j, k, l, current, prev, next, idx, jdx;
    double w0, wi[3], L, c, local_error;
    double** xi;
    double *axi, *product;


    xi = (double**) malloc(3 * sizeof (double));
    if (xi == NULL)
        die("\nBlad alokacji cheb");

    product = (double*) malloc(size * sizeof (double));
    if (product == NULL)
        die("\nBlad alokacji cheb");

    for (j = 0; j < 3; j++) {
        xi[j] = (double*) malloc(size * sizeof (double));
        if (xi[j] == NULL)
            die("\nBlad alokacji cheb");
    }
    axi = (double*) malloc(size * sizeof (double));
    if (axi == NULL)
        die("\nBlad alokacji cheb");

    i = 1;
    k = 0;
    w0 = (beta - alfa) / (beta + alfa);
    c = 2.0 / (alfa + beta);
    L = (2.0 * (beta + alfa)) / (beta - alfa);

    for (j = 0; j < size; j++) {
        xi[0][j] = 0;
        xi[1][j] = 0;
    }

    wi[0] = 0;
    wi[1] = w0;
double sum = 0;


#pragma acc data copy(product[0:size],axi[0:size],xi[0:3][0:size]) copyin( A[0:size][size])
    {

        while (1) {

            current = i % 3;
            prev = (i - 1) % 3;
            next = (i + 1) % 3;

#pragma acc region present(product[0:size],axi[0:size],xi[0:3][0:size],A[0:size][size]) copy(wi[0:3]) copyout(local_error)
            {


#pragma acc loop independent  vector(16)
                for (j = 0; j < size; j++) {
                     sum = 0;
                    for (l = 0; l < size; l++) {
                        sum += A[l][j] * xi[current][l];
                    }
                    axi[j] = sum;
                }
#pragma acc loop independent  vector(16)
                for (j = 0; j < size; j++) {
                    xi[next][j] = xi[current][j] + wi[current] * wi[prev] * (xi[current][j] - xi[prev][j]) - c * (1 + wi[current] * wi[prev]) * (axi[j] - A[j][size]);
                    wi[next] = 1 / (L - wi[current]);
                }

                //matrix multiply
                //matrix_multiply(A, xi[next], product, size);

#pragma acc loop independent  vector(16)
                for (idx = 0; idx < size; idx++) {
                    sum = 0;
                    for (jdx = 0; jdx < size; jdx++)
                         sum += A[idx][jdx] * xi[next][jdx];
                    product[idx] =sum;
                }

                //distance
                //local_error = cartesian_distance(product, A, size);

#pragma acc loop independent  vector(16)
                for (idx = 0; idx < size; idx++) {
                    local_error += (product[idx] - A[idx][size]) * (product[idx] - A[idx][size]);
                }
                local_error = sqrt(local_error);
            }

            k++;
            i++;

            if (i >= max_iter || error >= local_error) {
                break;
            }
            if (k < s)
                continue;
            else {
                k = 0;
                wi[current] = 0;
                wi[next] = w0;
            }
        }
    }


    long int elapsed = 1;
    if (parallel)
        printf("czas wiele watkow:\t%ld μs blad:\t%lf\n", elapsed, local_error);
    else
        printf("czas jeden watek:\t%ld μs blad:\t%lf\n", elapsed, local_error);

    for (j = 0; j < 2; j++) {
        free(xi[j]);
    }
    free(xi);
    free(axi);
    free(product);
    return elapsed;
}

int main() {
    read_matrix();
    int iter = 0;
    long max_iterations = 1;
    printf("-------------------------------------------\n");
    printf("------------------Cheb---------------------\n");
    printf("-------------------------------------------\n");

    for (iter = 0; iter < max_iterations; iter++) {
       cheb_parallel(input_matrix, size, 100, 2 * max, 30, 1000, 0.01, 0);
    }
    return 0;
}
