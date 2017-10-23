//
// Created by Fernando Montes on 10/14/17.
//
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void Read_n(int* n_p, int* local_n_p, int my_rank, int comm_sz, MPI_Comm comm);
void Read_vector(double local_a[], int local_n, int n, char vec_name[], int my_rank, MPI_Comm comm);
void Print_vector(double local_b[], int local_n, int n, char title[], int my_rank, MPI_Comm comm);


int main() {
    int n, local_n;
    int comm_sz, my_rank;
    double *local_x, *local_y, *local_z;
    MPI_Comm comm;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    Read_n(&n, &local_n, my_rank, comm_sz, comm);

    local_x = malloc(local_n*sizeof(double));
    local_y = malloc(local_n*sizeof(double));
    local_z = malloc(local_n*sizeof(double));

    //Reads in vector x
    Read_vector(local_x, local_n, n, "x", my_rank, comm);

    //Reads in vector y
    Read_vector(local_y, local_n, n, "y", my_rank, comm);

    for (int i = 0; i < local_n; i++) {
        local_z[i] = local_x[i] + local_y[i];
    }

    Print_vector(local_x, local_n, n, "x is", my_rank, comm);
    Print_vector(local_y, local_n, n, "y is", my_rank, comm);
    Print_vector(local_z, local_n, n, "The sum is", my_rank, comm);

    free(local_x);
    free(local_y);
    free(local_z);

    MPI_Finalize();

    return 0;
}

void Read_n(int* n_p, int* local_n_p, int my_rank, int comm_sz, MPI_Comm comm) {

    if (my_rank == 0) {
        printf("What's the order of the vectors?\n");
        scanf("%d", n_p);
    }

    MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
    *local_n_p = *n_p/comm_sz;
}

void Read_vector(double local_a[], int local_n , int n, char vec_name[], int my_rank, MPI_Comm comm) {
    double* a = NULL;
    int i;

    MPI_Datatype input;
    MPI_Type_contiguous(local_n, MPI_DOUBLE,&input);
    MPI_Type_commit(&input);

    if (my_rank == 0) {
        a = malloc(n*sizeof(double));
        printf("Enter the vector %s\n", vec_name);
        for (i = 0; i < n; i++) {
            scanf("%lf", &a[i]);
        }
    }

    MPI_Scatter(a, 1, input, local_a, 1, input, 0, comm);
    MPI_Type_free(&input);

    free(a);
}

void Print_vector(double local_b[], int local_n, int n, char title[], int my_rank, MPI_Comm comm) {
    double* b = NULL;
    int i;

    MPI_Datatype output;
    MPI_Type_contiguous(local_n, MPI_DOUBLE,&output);
    MPI_Type_commit(&output);

    b = malloc(n*sizeof(double));

    MPI_Gather(local_b, 1, output, b, 1, output, 0, comm);

    if (my_rank == 0) {
        printf("%s\n", title);
        for (i = 0; i < n; i++) {
            printf("%f ", b[i]);
        }
        printf("\n");
    }
    free(b);
    MPI_Type_free(&output);
}