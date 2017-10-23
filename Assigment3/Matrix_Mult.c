//
// Created by Fernando Montes on 10/23/17.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void Get_dims(int* m_p, int* local_m_p, int* n_p, int* local_n_p, int my_rank, int comm_sz, MPI_Comm comm);
void Allocate_arrays(double** local_A_pp, double** local_x_pp, double** local_y_pp, int local_m, int n, int local_n);
void Read_matrix(char prompt[], double local_A[], int m, int local_m, int n, int my_rank, MPI_Comm comm);
void Read_vector(char prompt[], double local_vec[], int n, int local_n, int my_rank, MPI_Comm comm);
void Print_matrix(char title[], double local_A[], int m, int local_m, int n, int my_rank, MPI_Comm comm);
void Print_vector(char title[], double local_vec[], int n, int local_n, int my_rank, MPI_Comm comm);
void Mat_vect_mult(double local_A[], double local_x[], double local_y[], int local_m, int n, int local_n, MPI_Comm comm);

int main(void) {
    double* local_A;
    double* local_x;
    double* local_y;
    int m, local_m, n, local_n;
    int my_rank, comm_sz;
    MPI_Comm comm;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    Get_dims(&m, &local_m, &n, &local_n, my_rank, comm_sz, comm);
    Allocate_arrays(&local_A, &local_x, &local_y, local_m, n, local_n);
    Read_matrix("A", local_A, m, local_m, n, my_rank, comm);
    Read_vector("y", local_x, n, local_n, my_rank, comm);
    Mat_vect_mult(local_A, local_x, local_y, local_m, n, local_n, comm);
    Print_vector("result", local_y, m, local_m, my_rank, comm);

    free(local_A);
    free(local_x);
    free(local_y);
    MPI_Finalize();
    return 0;
}  /* main */

void Get_dims(int* m_p, int* local_m_p, int* n_p, int* local_n_p, int my_rank, int comm_sz, MPI_Comm  comm) {

    if (my_rank == 0) {
        printf("Enter n (for n * n matrix): \n");
        scanf("%d", m_p);
        *n_p = *m_p;
    }
    MPI_Bcast(m_p, 1, MPI_INT, 0, comm);
    MPI_Bcast(n_p, 1, MPI_INT, 0, comm);

    *local_m_p = *m_p/comm_sz;
    *local_n_p = *n_p/comm_sz;
}  /* Get_dims */

void Allocate_arrays(double**  local_A_pp, double**  local_x_pp, double**  local_y_pp, int local_m, int n, int local_n) {

    *local_A_pp = malloc(local_m*n*sizeof(double));
    *local_x_pp = malloc(local_n*sizeof(double));
    *local_y_pp = malloc(local_m*sizeof(double));

}  /* Allocate_arrays */

void Read_matrix(char prompt[], double local_A[], int m, int local_m , int n, int my_rank, MPI_Comm  comm) {
    double* A = malloc(m*n*sizeof(double));;

    if (my_rank == 0) {
        printf("Enter the matrix %s\n", prompt);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                scanf("%lf", &A[i*n+j]);
        MPI_Scatter(A, local_m*n, MPI_DOUBLE, local_A, local_m*n, MPI_DOUBLE, 0, comm);
    } else {
        MPI_Scatter(A, local_m*n, MPI_DOUBLE, local_A, local_m*n, MPI_DOUBLE, 0, comm);
    }
    free(A);
}  /* Read_matrix */

void Read_vector(char prompt[], double local_vec[], int n, int local_n, int my_rank, MPI_Comm  comm) {
    double* vec = NULL;
    int i;

    if (my_rank == 0) {
        vec = malloc(n*sizeof(double));
        printf("Enter the vector %s\n", prompt);
        for (i = 0; i < n; i++)
            scanf("%lf", &vec[i]);
        MPI_Scatter(vec, local_n, MPI_DOUBLE,
                    local_vec, local_n, MPI_DOUBLE, 0, comm);
        free(vec);
    } else {
        MPI_Scatter(vec, local_n, MPI_DOUBLE,
                    local_vec, local_n, MPI_DOUBLE, 0, comm);
    }
}  /* Read_vector */

void Print_matrix(char title[], double local_A[], int m, int local_m, int n, int my_rank, MPI_Comm comm) {
    double* A = malloc(m*n*sizeof(double));

    if (my_rank == 0) {

        MPI_Gather(local_A, local_m*n, MPI_DOUBLE, A, local_m*n, MPI_DOUBLE, 0, comm);
        printf("\nThe matrix %s\n", title);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                printf("%f ", A[i*n+j]);
            printf("\n");
        }
        printf("\n");

    } else {
        MPI_Gather(local_A, local_m*n, MPI_DOUBLE, A, local_m*n, MPI_DOUBLE, 0, comm);
    }

    free(A);
}  /* Print_matrix */

void Print_vector(char title[], double local_vec[], int n, int local_n, int my_rank, MPI_Comm  comm) {
    double* vec  = malloc(n*sizeof(double));

    if (my_rank == 0) {
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vec, local_n, MPI_DOUBLE, 0, comm);
        printf("\nThe vector %s\n", title);
        for (int i = 0; i < n; i++) {
            printf("%f ", vec[i]);
        }
        printf("\n");

    }  else {
        MPI_Gather(local_vec, local_n, MPI_DOUBLE, vec, local_n, MPI_DOUBLE, 0, comm);
    }

    free(vec);
}  /* Print_vector */

void Mat_vect_mult(double local_A[], double local_x[], double local_y[], int local_m, int n, int local_n, MPI_Comm  comm) {
    double* x = malloc(n*sizeof(double));

    MPI_Allgather(local_x, local_n, MPI_DOUBLE, x, local_n, MPI_DOUBLE, comm);

    for (int local_i = 0; local_i < local_m; local_i++) {
        local_y[local_i] = 0.0;
        for (int j = 0; j < n; j++)
            local_y[local_i] += local_A[local_i*n+j]*x[j];
    }
    free(x);
}  /* Mat_vect_mult */