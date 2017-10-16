//
// Created by Fernando Montes on 10/16/17.
//

#include <stdio.h>

#include <mpi.h>

#define NUM_PARAMS 12


void get_model_params(double params[NUM_PARAMS])

{

    int i;


    /* Read parameters from file.  Here we just make them up. */

    for (i = 0; i < NUM_PARAMS; i++)

        params[i] = (double)(i * i);


}


int main(int argc, char** argv)

{

    int i, rank, size;

    double params[NUM_PARAMS];

    MPI_Datatype ParameterType;


    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);


    if (!rank) /* Only the first process reads parameters */ get_model_params(params);


    MPI_Type_contiguous(NUM_PARAMS, MPI_DOUBLE, &ParameterType);

    MPI_Type_commit(&ParameterType);

    /* First process broadcasts parameters to all processes */

    MPI_Bcast(params, 1, ParameterType, 0, MPI_COMM_WORLD);

    MPI_Type_free(&ParameterType);


    printf("%d: Got parameters: ", rank);

    for (i = 0; i < NUM_PARAMS; i++)

        printf(" %lf", params[i]);

    printf("\n");


    MPI_Finalize();

}