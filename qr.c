#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct{
    int rows;
    int cols;
    double **v_arr;
}matrix;

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
matrix* get_random_matrix(int num_rows, int num_cols, int zakres) {
    matrix *m = init_matrix(num_rows, num_cols);
    srand(time(NULL));
    int i,j;
    for (i=0; i<m->rows; i++){
        for (j=0; j<m->cols; j++){
            if (i==j){
                m->v_arr[i][j]=rand()%zakres+(num_cols-1)*zakres;
            } else {
                m->v_arr[i][j]=rand()%zakres;
            }
        }
    }
    return m;
}

matrix* make_zeroes_matrix(int num_rows, int num_cols){
    int i,j;
    matrix* m = init_matrix(num_rows, num_cols);
    for (i=0; i<num_rows; i++){
        for (j=0; j<num_cols; j++){
            m->v_arr[i][j] = 0.0;
        }
    }
    return m;
}

double get_frobenius_norm(matrix *m){
    double result = 0.0;
    int i, j;
    for (i=0; i<m->rows; i++){
        for (j=0; j<m->cols; j++){
            result = result + (m->v_arr[i][j]*m->v_arr[i][j]);
        }
    }
    return sqrt(result);
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
    matrix *minorA = make_zeroes_matrix(A->rows, A->cols);
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
            multiply_qtab[i] = make_zeroes_matrix(A->rows, A->cols);
            multiply_matricies(qtab[i-1],A,multiply_qtab[i]);
        }
        qtab[i] = make_zeroes_matrix(A->rows, A->cols);
        householder_decomposition(multiply_qtab[i], qtab[i], i);
    }

    matrix *Q1 = make_zeroes_matrix(A->rows, A->cols);
    matrix *Q2 = make_zeroes_matrix(A->rows, A->cols);
    Q1 = qtab[0];
    Q2 = qtab[1];
    //inicjuje, zeby nie sypnelo segfaultem
    matrix *tmp_multiply_qtab[A->cols-2];
    for (i=0; i<A->cols-2; i++){
        tmp_multiply_qtab[i] = make_zeroes_matrix(A->rows, A->cols);
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
    matrix *Qtransposed = make_zeroes_matrix(A->rows, A->cols);

    make_copy_of_matrix(Q, Qtransposed);
    get_matrix_transposition(Qtransposed);
    multiply_matricies(Qtransposed, A, R);
    free_matrix_memory(Qtransposed);
}



int main(){
    //Losowa macierz:
    matrix *a = get_random_matrix(10,10,100);
    print_matrix(a);
    double frobenius = get_frobenius_norm(a);
    printf("Frobenius: %f\n", frobenius);
    free_matrix_memory(a);
    double czas = 0.0;
    int i=0;
    //for (; i<100; ++i) {
        int j,k;
        clock_t tic = clock();

        matrix *A = get_random_matrix(300,300,100);
        matrix *Q = make_zeroes_matrix(A->rows, A->cols);
        matrix *R = make_zeroes_matrix(A->rows, A->cols);
        make_QR_decomposition(A, Q, R);

        matrix *QR = make_zeroes_matrix(A->rows, A->cols);
        multiply_matricies(Q,R,QR);
        clock_t toc = clock();

        FILE *fp;
        if ((fp=fopen("out.txt", "w"))==NULL) {
            printf ("Nie moge otworzyc pliku out.txt do zapisu!\n");
            exit(1);
        }
        for (j = 0; j <QR->rows;++j) {
            for (k=0; k<QR->rows;++k) {
                 fprintf (fp, "%f ", QR->v_arr[j][k]);
            }
            fprintf (fp, "\n");
        }

        fclose (fp); // zamknij plik

        free_matrix_memory(QR);
        free_matrix_memory(R);
        free_matrix_memory(Q);
        free_matrix_memory(A);
        czas+=(double)(toc-tic)*1000/CLOCKS_PER_SEC;
    //}
    printf("Czas wykonania: %5.1f [s]", czas/1000);

    return 0;
}
