//
// Created by Fernando Montes on 10/25/17.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

int main(int argc, char* argv[]){
    int my_rank, comm_sz,localTosses, number_in_circle= 0, global_num_in_circle =0;
    double distance_squared, x, y, estimate, number_of_tosses =0.0;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);



    //sscanf(argv[1], "%lf", &number_of_tosses);
    if(my_rank==0) {
        printf("Enter number of points to use: \n");
        scanf("%lf", &number_of_tosses);
    }
    MPI_Bcast(&number_of_tosses,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    localTosses = lround(number_of_tosses)/comm_sz;

    srand((unsigned)(time(0)));


    for (int toss = 0; toss< localTosses; toss++){
        x = (((double)rand())/((double)RAND_MAX))*2-1;
        y = (((double)rand())/((double)RAND_MAX))*2-1;
        distance_squared =( (x*x) + (y*y) );
        if (distance_squared <= 1 ){
            number_in_circle++;
        }

    }


    MPI_Allreduce(&number_in_circle,&global_num_in_circle,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    estimate = (global_num_in_circle*4)/(number_of_tosses);
    MPI_Barrier(MPI_COMM_WORLD);


    if (my_rank == 0) {
        printf("Pi Estimation:  %lf\n",estimate);
    }

    MPI_Finalize();
    return 0;
}
