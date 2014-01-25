#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct{
    int rows;
    int cols;
    double **v_arr;
}matrix;

double c [][30] = {
    {3,1,2,1,4,3,5,1,2,2,0,0.01,1,2,1,1,6,7,3,5,1,2,3,3,2,4,4,3,2},
{4,2,3,1,6,7,3,5,3,1,2,2,8,5,3,4,1,6,7,3,5,6,8,1,2,4,2,1,5,5},
{0,0.01,1,2,0,0.01,1,6,7,3,5,3,2,4,4,8,1,2,1,4,4,3,2,1,4,4,3,2},
{4,2,3,1,2,1,4,1,3,1,2,1,1,6,7,3,5,1,3,5,2,1,4,3,5,1,2,2,2,5},
{1,6,7,3,5,3,1,6,7,3,5,6,8,1,2,4,2,4,4,8,1,6,7,3,5,1,2,3,3,2},
{29,9,15,4,8,12,8,22,8,26,26,2,22,25,25,17,28,11,3,13,20,14,21,30,8,10,2,7,2,9},
{6,5,24,4,26,21,5,18,16,21,28,11,6,24,1,3,20,12,22,11,23,18,11,27,6,30,2,19,25,23},
{1,3,1,21,18,9,11,27,28,1,4,12,8,27,25,16,25,19,7,8,24,30,16,29,23,12,13,24,29,26},
{9,7,23,22,4,23,10,27,27,21,8,24,26,23,15,7,3,4,28,14,17,13,1,15,3,20,7,27,2,9},
{10,12,2,18,17,23,3,17,29,23,14,4,14,16,23,3,16,4,3,11,13,23,2,6,25,14,18,29,16,2},
{10,28,4,10,22,4,28,6,7,20,1,18,22,26,1,1,14,18,13,4,22,13,20,1,21,14,19,9,10,20},
{16,22,3,23,12,11,6,4,21,25,6,11,3,8,22,8,22,1,9,29,1,9,11,15,18,5,2,9,2,16},
{18,12,7,30,16,28,10,27,21,29,21,13,30,9,10,14,11,17,7,9,11,5,1,23,18,3,14,22,25,17},
{26,1,10,19,23,23,1,11,5,30,9,25,5,18,24,24,26,2,13,19,1,28,9,27,4,1,8,20,4,10},
{7,9,30,30,4,3,12,8,18,5,11,28,15,20,27,3,28,3,5,9,14,26,4,23,26,23,24,6,20,25},
{19,15,28,24,12,14,29,6,11,8,3,21,27,29,2,13,8,8,29,26,17,7,22,1,26,11,14,12,26,7},
{18,23,15,11,12,6,29,2,13,6,10,1,9,1,18,18,27,7,7,12,19,6,25,18,16,7,17,24,16,7},
{3,30,3,16,3,7,29,28,23,4,1,7,7,28,10,28,22,4,7,19,27,3,28,14,19,17,18,26,16,17},
{19,8,7,30,28,9,13,2,24,18,15,7,29,18,9,29,28,28,8,25,25,3,14,25,28,8,5,17,20,2},
{9,28,23,22,6,19,11,9,13,4,6,29,25,16,11,27,1,5,7,18,28,5,13,27,1,9,14,4,28,14},
{25,2,5,27,20,27,13,8,21,8,18,28,2,16,30,8,8,28,9,30,21,10,17,29,10,22,1,7,6,2},
{23,18,25,16,16,18,5,2,24,10,4,17,1,9,10,23,18,7,22,18,24,13,22,3,27,24,3,9,13,2},
{30,8,4,3,16,30,21,27,22,22,20,20,14,30,5,6,11,12,8,11,5,27,1,24,21,5,21,21,5,27},
{7,22,11,30,2,16,23,22,6,19,7,14,30,29,6,30,22,7,7,26,16,19,10,8,12,15,8,15,11,7},
{23,18,1,13,23,21,15,26,13,20,26,13,24,30,20,21,3,11,19,22,13,25,17,5,23,17,3,3,23,1},
{25,14,15,23,1,27,11,17,29,15,4,3,3,1,5,24,25,30,24,13,30,25,13,12,4,3,11,23,29,30},
{29,12,28,11,1,26,13,13,16,21,14,11,21,2,27,8,26,15,19,4,27,4,9,5,23,2,1,10,14,27},
{5,2,28,10,6,14,28,6,13,10,3,21,26,17,2,27,29,12,10,16,14,28,12,23,8,26,7,8,12,24},
{26,4,26,3,13,30,12,21,15,24,12,20,8,11,7,28,29,14,27,5,24,12,22,10,18,22,15,5,23,1},
{15,2,14,23,9,8,5,7,1,20,23,27,14,21,3,26,22,19,7,13,29,15,4,16,3,15,17,15,9,29},
};

double zeros [][30] = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};


// samo zaalokowanie pamieci
matrix* init_matrix(int num_rows, int num_cols){
    matrix* m = malloc(sizeof(matrix));
    m->rows = num_rows;
    m->cols = num_cols;
    m->v_arr = malloc(sizeof(double)*num_rows);
    int i,j;
    for (i=0; i<num_rows; i++){
        m->v_arr[i]=malloc(sizeof(double)*num_cols);
    }
    return m;
}

void free_matrix_memory(matrix* m){
    int i;

    for (i=0; i<m->rows; i++){
        free(m->v_arr[i]);
    }
    free(m->v_arr);
    free(m);
}

void print_matrix(matrix* m){
    int i,j;
    for (i=0; i<m->rows; i++){
        for (j=0; j<m->cols; j++){
            printf(" %8.3f", m->v_arr[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// wypelnienie wartosciami - trzeba zawsze nowa macierz wypelnic zerami, inaczej czesto zostaja smieci w pamieci i zle sie liczy
matrix* get_matrix_data(int cols; double a[][cols], int rows, int cols){
    int i,j;
    matrix* m = init_matrix(rows, cols);
    for (i=0; i<rows; i++){
        for (j=0; j<cols; j++){
            m->v_arr[i][j] = a[i][j];
        }
    }
    return m;
}

double get_euclides_norm(int size_x; double x[size_x], int size_x){
    int i;
    double norm=0.0;
    for (i=0; i<size_x; i++){
        norm+=(x[i]*x[i]);
    }
    return sqrt(norm);
}

void get_e_vector(int size_e; double *e, int pos_of_one, int size_e){
    int i;
    for (i=0; i<size_e; i++){
        e[i] = (i==pos_of_one) ? 1.0 : 0.0;
    }
}

void get_matrix_transposition(matrix* m){
    int i,j;
    double tmp_elem;

    for (i=0; i<m->rows; i++){
        for (j=0; j<i; j++){
            tmp_elem = m->v_arr[i][j];
            m->v_arr[i][j]= m->v_arr[j][i];
            m->v_arr[j][i] = tmp_elem;
        }
    }
}

void multiply_matricies(matrix *m1, matrix *m2, matrix *result){
    //result = init_matrix(m1->rows, m2->cols);
    int i,j,k;
    #pragma omp parallel for private(i, j, k) shared (result, m1,m2)
    for (i=0; i<m1->rows; i++){
        for (j=0; j<m2->cols; j++){
            for (k=0; k<m1->cols; k++){
                result->v_arr[i][j] += m1->v_arr[i][k] * m2->v_arr[k][j];
            }
        }
    }
}

void multiply_column_by_scalar(int size_n; double *n, int size_n, double scalar){
    int i;
    for (i=0; i<size_n; i++){
        n[i]=n[i]*scalar;
    }
}

void divide_column_by_scalar(int size_n; double n[size_n], int size_n, double scalar){
    int i;
    for (i=0; i<size_n; i++){
        n[i]=n[i]/scalar;
    }
}

void get_v(int size_u; double u[size_u], int size_u){
    double norm_u = get_euclides_norm(u, size_u);
    divide_column_by_scalar(u, size_u, norm_u);
}

void subtract_matricies(matrix *m1, matrix *m2, matrix *result){
    int i,j;
    for (i=0; i<m1->rows; i++){
        for (j=0; j<m1->cols; j++){
            result->v_arr[i][j]=m1->v_arr[i][j]-m2->v_arr[i][j];
        }
    }
}

matrix* init_ones_matrix(int rows, int cols){
    if (rows != cols)
        return 0;
    matrix* ones = init_matrix(rows, cols);
    int i,j;
    for (i=0; i<ones->rows; i++){
        for (j=0; j<ones->cols; j++){
            ones->v_arr[i][j] = (i==j) ? 1.0 : 0.0;
        }
    }
    return ones;
}

void get_minor_of_matrix(matrix *m, matrix *x, int d){
    int i,j;
    for (i=0; i<d; i++){
        x->v_arr[i][i] = 1.0;
    }
    for (i=d; i<m->rows; i++){
        for (j=d; j<m->cols; j++){
            x->v_arr[i][j] = m->v_arr[i][j];
        }
    }
}

void get_nth_column_of_matrix(matrix* m, int n, double *v){
    int i;
    for (i=0; i<m->rows; i++){
        v[i] = m->v_arr[i][n];
    }
}

void add_two_columns(double *col1, double *col2, int size_cols){
    int i;

    for (i=0; i<size_cols; i++){
        col1[i]+=col2[i];
    }
}

matrix *get_multiply_v(double *v, int size_v){
    matrix *tmp = init_matrix(size_v, size_v);
    int i,j;

    for (i=0; i<tmp->rows; i++){
        for (j=0; j<tmp->cols; j++){
            tmp->v_arr[i][j] = 2*v[i]*v[j];
        }
    }
    return tmp;
}

void make_copy_of_matrix(matrix *from, matrix *to){
    int i,j;
    for (i=0; i<from->rows; i++){
        for (j=0; j<from->cols; j++){
            to->v_arr[i][j] = from->v_arr[i][j];
        }
    }
}

void householder_decomposition(matrix *A, matrix *Q, int n){
    matrix *minorA = get_matrix_data(zeros, A->rows, A->cols);
    get_minor_of_matrix(A, minorA, n);
    double *x = malloc(sizeof(double)*minorA->rows);
    get_nth_column_of_matrix(minorA, n, x);
    double alfa = get_euclides_norm(x, minorA->rows);
    alfa = (x[n] > 0) ? -alfa : alfa;
    double *e = malloc(sizeof(double)*minorA->rows);
    get_e_vector(e,n,minorA->cols);
    multiply_column_by_scalar(e, minorA->rows, alfa);
    add_two_columns(x, e, minorA->rows);
    get_v(x, minorA->rows);
    matrix *tmp = get_multiply_v(x, minorA->rows);
    matrix *ones = init_ones_matrix(minorA->rows, minorA->cols);
    subtract_matricies(ones, tmp, Q);

    free_matrix_memory(ones);
    free_matrix_memory(tmp);
    free(e);
    free(x);
    free_matrix_memory(minorA);

    get_matrix_transposition(Q);
}

void make_QR_decomposition(matrix *A, matrix *Q, matrix *R){
    int i;
    matrix *qtab[A->cols-1];
    matrix *multiply_qtab[A->cols-1];
    multiply_qtab[0] = A;
    for (i=0; i<A->cols-1; i++){
        if (i>0){
            multiply_qtab[i] = get_matrix_data(zeros, A->rows, A->cols);
            multiply_matricies(qtab[i-1],A,multiply_qtab[i]);
        }
        qtab[i] = get_matrix_data(zeros, A->rows, A->cols);
        householder_decomposition(multiply_qtab[i], qtab[i], i);
    }

    matrix *Q1 = get_matrix_data(zeros, A->rows, A->cols);
    matrix *Q2 = get_matrix_data(zeros, A->rows, A->cols);
    Q1 = qtab[0];
    Q2 = qtab[1];
    //inicjuje, zeby nie sypnelo segfaultem
    matrix *tmp_multiply_qtab[A->cols-2];
    for (i=0; i<A->cols-2; i++){
        tmp_multiply_qtab[i] = get_matrix_data(zeros, A->rows, A->cols);
    }
    //glowna petla gdzie mnoze elementy
    for (i=0; i<A->cols-2; i++){
        multiply_matricies(Q1, Q2, tmp_multiply_qtab[i]);
        if (i+2 < A->cols-2){
            make_copy_of_matrix(tmp_multiply_qtab[i], Q1);
            make_copy_of_matrix(qtab[i+2], Q2);
        }
    }
    //i kopiuje ostatni elem. z tmp tablicy do Q
    make_copy_of_matrix(tmp_multiply_qtab[A->cols-3],Q);
    free_matrix_memory(Q2);
    free_matrix_memory(Q1);
    for (i=0; i<A->cols-2; i++)
        free_matrix_memory(tmp_multiply_qtab[i]);
    //a tu znowu odwaracam to Q
    //bo jest mi potrzebne do R (R = Q'*A)
    matrix *Qtransposed = get_matrix_data(zeros, A->rows, A->cols);
    make_copy_of_matrix(Q, Qtransposed);
    get_matrix_transposition(Qtransposed);
    multiply_matricies(Qtransposed, A, R);
    free_matrix_memory(Qtransposed);
}



int main(){
    double czas = 0.0;
    int i=0;
    for (; i<100; ++i) {
    int j,k;
    clock_t tic = clock();

    matrix *A = get_matrix_data(c, 30, 30);
    matrix *Q = get_matrix_data(zeros, A->rows, A->cols);
    matrix *R = get_matrix_data(zeros, A->rows, A->cols);
    make_QR_decomposition(A, Q, R);

    matrix *QR = get_matrix_data(zeros, A->rows, A->cols);
    multiply_matricies(Q,R,QR);
    clock_t toc = clock();

    FILE *fp;
    if ((fp=fopen("out.txt", "w"))==NULL) {
        printf ("Nie mogê otworzyæ pliku test.txt do zapisu!\n");
        exit(1);
    }
    for (j = 0; j <QR->rows;++j)
    {
        for (k=0; k<QR->rows;++k) {
             fprintf (fp, "%f ", QR->v_arr[j][k]);
        }
        fprintf (fp, "\n");
    }

    fclose (fp); /* zamknij plik */

    free_matrix_memory(QR);
    free_matrix_memory(R);
    free_matrix_memory(Q);
    free_matrix_memory(A);
    czas+=(double)(toc-tic)*1000/CLOCKS_PER_SEC;
    }
    printf("Czas wykonania: %5.1f [ms]", czas/100);

    return 0;
}
