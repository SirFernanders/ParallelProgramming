//
// Created by Fernando Montes on 10/3/17.
//

/* File:      histogram.c
 * Purpose:   Build a histogram from some random data
 *
 * Compile:   mpicc -g -Wall -o histogram histogram.c
 * Run:       ./histogram <bin_count> <min_meas> <max_meas> <data_count>
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void Usage(char prog_name[]);

void Get_args(
        char*    argv[]        /* in  */,
        int*     bin_count_p   /* out */,
        float*   min_meas_p    /* out */,
        float*   max_meas_p    /* out */,
        int*     data_count_p  /* out */,
        int*     local_data_count_p,
        int      my_rank,
        int      comm_sz,
        MPI_Comm comm
);

void Gen_data(
        float   min_meas    /* in  */,
        float   max_meas    /* in  */,
        float   data[]      /* out */,
        int     data_count  /* in  */,
        float   local_data[],
        int     local_data_count,
        int my_rank,
        MPI_Comm comm
);


void Gen_bins(
        float min_meas      /* in  */,
        float max_meas      /* in  */,
        float bin_maxes[]   /* out */,
        int   local_bin_count[]  /* out */,
        int   bin_count     /* in  */);

int Which_bin(
        float    data         /* in */,
        float    bin_maxes[]  /* in */,
        int      bin_count    /* in */,
        float    min_meas     /* in */);

void Print_histo(
        float    bin_maxes[]   /* in */,
        int      bin_counts[]  /* in */,
        int      bin_count     /* in */,
        float    min_meas      /* in */);

int main(int argc, char* argv[]) {
    int bin_count, i, bin;
    float min_meas, max_meas;
    float* bin_maxes;
    int* bin_counts;
    int* local_bin_count;
    int data_count;
    int local_d_count;
    float* data;
    float* local_d;
    int my_rank;
    int comm_sz;
    MPI_Comm comm;

    MPI_Init(NULL,NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    /* Check and get command line args */
    if (argc != 5) Usage(argv[0]);
    Get_args(argv, &bin_count, &min_meas, &max_meas, &data_count, &local_d_count, my_rank, comm_sz, comm);

    /* Allocate arrays needed */
    bin_maxes = malloc(bin_count*sizeof(float));
    bin_counts = malloc(bin_count*sizeof(int));
    local_bin_count = malloc(bin_count*sizeof(int));
    local_d = malloc(bin_count*sizeof(float));
    data = malloc(data_count*sizeof(float));

    /* Generate the data */
    Gen_data(min_meas, max_meas, data, data_count,local_d, local_d_count, my_rank, comm);

    /* Create bins for storing counts */
    Gen_bins(min_meas, max_meas, bin_maxes, local_bin_count, bin_count);

    /* Count number of values in each bin */
    for (i = 0; i < local_d_count; i++) {
        bin = Which_bin(local_d[i], bin_maxes, bin_count, min_meas);
        local_bin_count[bin]++;
    }

#  ifdef DEBUG
    printf("bin_counts = ");
   for (i = 0; i < bin_count; i++)
      printf("%d ", bin_counts[i]);
   printf("\n");
#  endif

    MPI_Reduce(local_bin_count,bin_counts,bin_count,MPI_INT,MPI_SUM,0,comm);

    /* Print the histogram */
    if(my_rank==0)
        Print_histo(bin_maxes, bin_counts, bin_count, min_meas);

    free(local_bin_count);
    free(local_d);
    free(data);
    free(bin_maxes);
    free(bin_counts);
    MPI_Finalize();
    return 0;

}  /* main */



void Usage(char prog_name[] /* in */) {
    fprintf(stderr, "usage: %s ", prog_name);
    fprintf(stderr, "<bin_count> <min_meas> <max_meas> <data_count>\n");
    exit(0);
}  /* Usage */



void Get_args(
        char*    argv[]        /* in  */,
        int*     bin_count_p   /* out */,
        float*   min_meas_p    /* out */,
        float*   max_meas_p    /* out */,
        int*     data_count_p  /* out */,
        int*     local_data_count_p,
        int      my_rank,
        int      comm_sz,
        MPI_Comm comm
) {

    if (my_rank == 0) {
        *bin_count_p = strtol(argv[1], NULL, 10);
        *min_meas_p = strtof(argv[2], NULL);
        *max_meas_p = strtof(argv[3], NULL);
        *data_count_p = strtol(argv[4], NULL, 10);


#  ifdef DEBUG
        printf("bin_count = %d\n", *bin_count_p);
               printf("min_meas = %f, max_meas = %f\n", *min_meas_p, *max_meas_p);
               printf("data_count = %d\n", *data_count_p);
#  endif
    }
    *local_data_count_p = *data_count_p / comm_sz;
    *data_count_p = *local_data_count_p * comm_sz;

    MPI_Bcast(bin_count_p,1,MPI_INT,0,comm);
    MPI_Bcast(min_meas_p,1,MPI_FLOAT,0,comm);
    MPI_Bcast(max_meas_p,1,MPI_FLOAT,0,comm);
    MPI_Bcast(data_count_p,1,MPI_INT,0,comm);
    MPI_Bcast(local_data_count_p,1,MPI_INT,0,comm);
}  /* Get_args */



void Gen_data(
        float   min_meas    /* in  */,
        float   max_meas    /* in  */,
        float   data[]      /* out */,
        int     data_count  /* in  */,
        float local_d[], int local_d_count,
        int my_rank,
        MPI_Comm comm) {
    int i;

    if (my_rank==0) {
        srandom(0);
        for (i = 0; i < data_count; i++)
            data[i] = min_meas + (max_meas - min_meas) * random() / ((double) RAND_MAX);

#  ifdef DEBUG
        printf("data = ");
               for (i = 0; i < data_count; i++)
                  printf("%4.3f ", data[i]);
               printf("\n");
#  endif
    }

    MPI_Scatter(data, local_d_count,MPI_FLOAT,local_d, local_d_count,MPI_FLOAT,0,comm);
    if(my_rank==0)
        free(data);
}  /* Gen_data */



void Gen_bins(
        float min_meas      /* in  */,
        float max_meas      /* in  */,
        float bin_maxes[]   /* out */,
        int   bin_counts[]  /* out */,
        int   bin_count     /* in  */) {
    float bin_width;
    int   i;

    bin_width = (max_meas - min_meas)/bin_count;

    for (i = 0; i < bin_count; i++) {
        bin_maxes[i] = min_meas + (i+1)*bin_width;
        bin_counts[i] = 0;
    }

#  ifdef DEBUG
    printf("bin_maxes = ");
   for (i = 0; i < bin_count; i++)
      printf("%4.3f ", bin_maxes[i]);
   printf("\n");
#  endif
}  /* Gen_bins */



int Which_bin(
        float   data          /* in */,
        float   bin_maxes[]   /* in */,
        int     bin_count     /* in */,
        float   min_meas      /* in */) {
    int bottom = 0, top =  bin_count-1;
    int mid;
    float bin_max, bin_min;

    while (bottom <= top) {
        mid = (bottom + top)/2;
        bin_max = bin_maxes[mid];
        bin_min = (mid == 0) ? min_meas: bin_maxes[mid-1];
        if (data >= bin_max)
            bottom = mid+1;
        else if (data < bin_min)
            top = mid-1;
        else
            return mid;
    }

    /* Whoops! */
    fprintf(stderr, "Data = %f doesn't belong to a bin!\n", data);
    fprintf(stderr, "Quitting\n");
    exit(-1);
}  /* Which_bin */


void Print_histo(
        float  bin_maxes[]   /* in */,
        int    bin_counts[]  /* in */,
        int    bin_count     /* in */,
        float  min_meas      /* in */) {
    int i, j;
    float bin_max, bin_min;

    for (i = 0; i < bin_count; i++) {
        bin_max = bin_maxes[i];
        bin_min = (i == 0) ? min_meas: bin_maxes[i-1];
        printf("%.3f-%.3f:\t", bin_min, bin_max);
        for (j = 0; j < bin_counts[i]; j++)
            printf("X");
        printf("\n");
    }
}  /* Print_histo */

