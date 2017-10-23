//
// Created by Fernando Montes on 10/14/17.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void arrayMaker(int array[], int size);

int main(){

    int local_n;
    int comm_sz, my_rank;
    int size = 6;
    int* array;
    int* stf;
    array = malloc(size*sizeof(int));
    int* local_x = malloc(size*sizeof(int));
    stf = malloc(size*sizeof(int));

    MPI_Comm comm;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    if(my_rank==0) {
        array[0] = 1;
        array[1] = 3;
        array[2] = 4;
        array[3] = 7;
    }

    local_n = 4 / comm_sz;
    MPI_Scatter(array, local_n, MPI_INT, local_x, local_n, MPI_INT, 0, comm);


MPI_Scan(&local_x, &stf, local_n, MPI_INT, MPI_SUM, comm);
    for(int i = 0; i<=local_n; i++) {

        //printf("processor %i's value is  %i and the prefix sum currently is: " ,my_rank, local_x[i]);
        //MPI_Scan(&local_x, &stf, local_n, MPI_INT, MPI_SUM, comm);

        printf("%i \n", stf[i]);
    }


    free(array);
    free(local_x);

    MPI_Finalize();
    return 0;

}

void arrayMaker(int array[], int size){

    int* temp = malloc(size* sizeof(int));



    *array = *temp;


}
